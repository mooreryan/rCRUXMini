#include <Rcpp.h>
#include <boost/container_hash/hash.hpp>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct StringLongPair {
    std::string string_value;
    long long_value{};

    StringLongPair() = default;

    StringLongPair(const std::string &str, long l)
        : string_value(str), long_value(l) {}

    bool operator==(const StringLongPair &other) const {
        return string_value == other.string_value &&
               long_value == other.long_value;
    }

    // This is a "friend" to work with boost
    friend std::size_t hash_value(const StringLongPair &key) {
        std::size_t seed = 0;
        boost::hash_combine(seed, key.string_value);
        boost::hash_combine(seed, key.long_value);
        return seed;
    }
};

struct AmpliconRegion {
    std::string saccver;
    std::string sgi;
    std::string staxids;
    long forward_start{};
    long forward_stop{};
    long forward_mismatch{};
    long reverse_start{};
    long reverse_stop{};
    long reverse_mismatch{};
    long product_length{};

    std::string toTsvString() const {
        std::ostringstream oss;
        oss << saccver << "\t" << sgi << "\t" << staxids << "\t"
            << forward_start << "\t" << forward_stop << "\t" << forward_mismatch
            << "\t" << reverse_start << "\t" << reverse_stop << "\t"
            << reverse_mismatch << "\t" << product_length;
        return oss.str();
    }

    static std::string tsvHeader() {
        return "saccver\tsgi\tstaxids\tforward_start\tforward_stop\t"
               "forward_mismatch\treverse_start\treverse_stop\t"
               "reverse_mismatch\tproduct_length";
    }
};

struct PrimerBlastHit {
    std::string qseqid;
    std::string sgi;
    std::string saccver;
    long mismatch{};
    long sstart{};
    long send{};
    std::string staxids;
};

// Replace throwing this with `Rcpp::stop` somehow...or do both?
class PrimerBlastParseException final : public std::runtime_error {
  public:
    explicit PrimerBlastParseException(const std::string &msg)
        : std::runtime_error(msg) {}
};

struct PartitionedHits {
    std::vector<PrimerBlastHit> forward_hits;
    std::vector<PrimerBlastHit> reverse_hits;
};

namespace {
void validateToken(const std::string &token, const std::string &field_name) {
    if (token.empty() ||
        token.find_first_not_of(" \t\n\r\f\v") == std::string::npos) {
        throw PrimerBlastParseException(
            "Empty or whitespace-only value for field: " + field_name);
    }
}

long safeStol(const std::string &str, const std::string &field_name) {
    if (str.empty()) {
        throw PrimerBlastParseException("Empty value for field: " + field_name);
    }

    try {
        size_t pos = 0;
        const long VALUE = std::stol(str, &pos);

        // Check if the entire string was consumed
        if (pos != str.length()) {
            throw PrimerBlastParseException("Invalid numeric value for field " +
                                            field_name + ": " + str);
        }
        return VALUE;
    } catch (const std::invalid_argument &e) {
        throw PrimerBlastParseException("Invalid numeric value for field " +
                                        field_name + ": " + str);
    } catch (const std::out_of_range &e) {
        throw PrimerBlastParseException(
            "Numeric value out of range for field " + field_name + ": " + str);
    }
}

PrimerBlastHit parsePrimerBlastHitLine(const std::string &line) {
    std::stringstream line_stream(line);
    std::vector<std::string> tokens;
    std::string token;

    // Read all fields
    while (std::getline(line_stream, token, '\t')) {
        tokens.push_back(token);
    }

    // Validate field count
    constexpr size_t EXPECTED_FIELDS = 7;
    if (tokens.size() != EXPECTED_FIELDS) {
        throw PrimerBlastParseException(
            "Expected " + std::to_string(EXPECTED_FIELDS) +
            " fields, but got " + std::to_string(tokens.size()));
    }

    PrimerBlastHit row;

    validateToken(tokens[0], "qseqid");
    row.qseqid = tokens[0];

    validateToken(tokens[1], "sgi");
    row.sgi = tokens[1];

    validateToken(tokens[2], "saccver");
    row.saccver = tokens[2];

    row.mismatch = safeStol(tokens[3], "mismatch");
    row.sstart = safeStol(tokens[4], "sstart");
    row.send = safeStol(tokens[5], "send");

    validateToken(tokens[6], "staxids");
    row.staxids = tokens[6];

    return row;
}

bool isForwardPrimer(const std::string &qseqid) {
    return qseqid.rfind("forward", 0) == 0;
}

bool isReversePrimer(const std::string &qseqid) {
    return qseqid.rfind("reverse", 0) == 0;
}

bool isHeaderLine(const std::string &line) {
    return line.rfind("qseqid", 0) == 0;
}

void removeAccessionsWithoutForwardAndReverseHits(
    std::unordered_map<std::string, PartitionedHits> &accession_hits) {
    size_t index = 0;

    auto it = accession_hits.begin();
    while (it != accession_hits.end()) {
        index++;

        if (index % 10000 == 0) {
            Rcpp::checkUserInterrupt();
        }

        auto partitioned_hits = it->second;
        auto forward_hit_count = partitioned_hits.forward_hits.size();
        auto reverse_hit_count = partitioned_hits.reverse_hits.size();
        auto has_forward_hits = forward_hit_count > 0;
        auto has_reverse_hits = reverse_hit_count > 0;

        if (has_forward_hits && has_reverse_hits) {
            // This accession has both, so it can stay.
            ++it;
        } else {
            // Erasing invalidates the iterator, so we need to
            // return a fresh one pointing to the correct next
            // element.
            //
            it = accession_hits.erase(it);
        }
    }
}

long calculate_product_length(const long forward_start,
                              const long forward_stop,
                              const long reverse_start,
                              const long reverse_stop) {
    // -1 represents invalid orientations
    long result = -1;

    //  F ---------->                        <------------ R  Primers
    // ====================================================== DNA Target
    //    ---------->************************<------------    Product
    //    ^ Forward start                  Reverse start ^

    if (forward_start < reverse_start && forward_start < forward_stop &&
        reverse_stop < reverse_start) {
        result = reverse_start - forward_start;
        assert(result >= 0);
        return result;
    }

    //  R ---------->                        <------------ F  Primers
    // ====================================================== DNA Target
    //    ---------->************************<------------    Product
    //    ^ Reverse start                  Forward start ^
    if (forward_start > reverse_start && forward_start > forward_stop &&
        reverse_stop > reverse_start) {
        result = forward_start - reverse_start;
        assert(result >= 0);
        return result;
    }

    return result;
}

// Get all the amplicon regions for one accession's partitioned hits.
std::vector<AmpliconRegion>
allAmpliconRegions(const PartitionedHits &partitioned_hits) {
    assert(!partitioned_hits.forward_hits.empty());
    assert(!partitioned_hits.reverse_hits.empty());

    std::vector<AmpliconRegion> amplicon_regions;

    for (auto &forward_hit : partitioned_hits.forward_hits) {
        for (auto &reverse_hit : partitioned_hits.reverse_hits) {
            auto product_length =
                calculate_product_length(forward_hit.sstart, forward_hit.send,
                                         reverse_hit.sstart, reverse_hit.send);

            AmpliconRegion amplicon_region{};
            // These _should_ be properly set up, so the forward and
            // reverse hit should have the same saccver.
            assert(forward_hit.saccver == reverse_hit.saccver);
            amplicon_region.saccver = forward_hit.saccver;
            amplicon_region.sgi = forward_hit.sgi;
            amplicon_region.staxids = forward_hit.staxids;

            amplicon_region.forward_start = forward_hit.sstart;
            amplicon_region.forward_stop = forward_hit.send;
            amplicon_region.forward_mismatch = forward_hit.mismatch;

            amplicon_region.reverse_start = reverse_hit.sstart;
            amplicon_region.reverse_stop = reverse_hit.send;
            amplicon_region.reverse_mismatch = reverse_hit.mismatch;

            amplicon_region.product_length = product_length;

            amplicon_regions.push_back(amplicon_region);
        }
    }

    return amplicon_regions;
}
} // namespace

// [[Rcpp::export(name = "parse_primer_blast")]]
void parsePrimerBlast(const std::string &primer_blast_tsv,
                      const std::string &output_tsv,
                      const long maximum_mismatches,
                      const long minimum_length,
                      const long maximum_length,
                      const int num_threads) {
    omp_set_num_threads(num_threads);

    std::ifstream primer_blast_stream(primer_blast_tsv);

    std::unordered_map<std::string, PartitionedHits> accession_hits;

    if (!primer_blast_stream.is_open()) {
        Rcpp::stop("Failed to open primer blast file");
    }

    std::string line;

    size_t line_index = 0;

    std::unordered_map<StringLongPair, std::vector<PrimerBlastHit>,
                       boost::hash<StringLongPair>>
        saccver_sstart_grouped_hits;

    while (std::getline(primer_blast_stream, line)) {
        line_index++;

        if (line_index % 10000 == 0) {
            Rcpp::checkUserInterrupt();
        }

        if (isHeaderLine(line)) {
            continue;
        }

        try {
            const PrimerBlastHit hit = parsePrimerBlastHitLine(line);

            // TODO: this should be >, but original code has it
            // wrong
            if (hit.mismatch >= maximum_mismatches) {
                continue;
            }

            auto key = StringLongPair{};
            key.string_value = hit.saccver;
            key.long_value = hit.sstart;
            saccver_sstart_grouped_hits[key].push_back(hit);
        } catch (const PrimerBlastParseException &e) {
            Rcpp::stop("Error parsing line: %s", line);
        }
    }

    std::vector<PrimerBlastHit> filtered_hits;
    for (auto pair : saccver_sstart_grouped_hits) {
        auto key = pair.first;
        auto hits = pair.second;
        assert(hits.size() > 0);

        // Find the minimum number of mismatches for this saccver+sstart group.
        long minimum_mismatches = hits[0].mismatch;
        for (auto hit : hits) {
            if (hit.mismatch < minimum_mismatches) {
                minimum_mismatches = hit.mismatch;
            }
        }

        // TODO: we should probably keep the longest hit with the fewest
        // mismatches rather than the first. Keep the first hit with the minimum
        // number of mismatches
        PrimerBlastHit first_hit_with_min_mismatches = hits[0];
        for (auto hit : hits) {
            if (hit.mismatch == minimum_mismatches) {
                first_hit_with_min_mismatches = hit;
                break;
            }
        }

        filtered_hits.push_back(first_hit_with_min_mismatches);
    }

    // TODO: original code has a group by saccver+send and keep distinct...but
    // idk if you can ever even get that since we are only ever keeping one of
    // the filtered hits anyway?
    std::unordered_map<StringLongPair, PrimerBlastHit,
                       boost::hash<StringLongPair>>
        unique_saccver_send_hits;
    for (const auto &hit : filtered_hits) {
        // If this key is not present in the map, add the hit
        // If the key is already present in the map, don't override it
        auto key = StringLongPair(hit.saccver, hit.send);
        unique_saccver_send_hits.emplace(key, hit);
    }

    // Extract unique hits back to original vector
    filtered_hits.clear();
    filtered_hits.reserve(unique_saccver_send_hits.size());
    for (const auto &pair : unique_saccver_send_hits) {
        filtered_hits.push_back(pair.second);
    }

    for (auto hit : filtered_hits) {
        // Add the hit if it is valid
        if (isForwardPrimer(hit.qseqid)) {
            accession_hits[hit.saccver].forward_hits.push_back(hit);
        } else if (isReversePrimer(hit.qseqid)) {
            accession_hits[hit.saccver].reverse_hits.push_back(hit);
        } else {
            // If neither a forward nor reverse hit, then
            // skip.
            ;
        }
    }

    // TODO: should this be guarded by the try/catch
    primer_blast_stream.close();

    removeAccessionsWithoutForwardAndReverseHits(accession_hits);

    // Open output file
    std::ofstream output_stream(output_tsv);
    if (!output_stream.is_open()) {
        Rcpp::stop("Error opening output file: %s", output_tsv);
    }

    output_stream << AmpliconRegion::tsvHeader() << std::endl;

    // Convert to vector for indexed access
    std::vector<std::pair<std::string, PartitionedHits>> accession_hits_vec(
        accession_hits.begin(), accession_hits.end());

    // The interrupt rate is connected to this....if it is too high,
    // interrupts won't work very well, if it is too low, then it won't use
    // cores effectively.
    const int CHUNK_SIZE = 50000;
    const int N = accession_hits_vec.size();

    for (int start = 0; start < N; start += CHUNK_SIZE) {
        Rcpp::checkUserInterrupt();

        int end = std::min(start + CHUNK_SIZE, N);

        // We really don't want to save ALL the region strings in
        // memory. So we give a new storage for each chunk we process.
        //
        // Each thread gets its own output buffer Use max_threds here in
        // case a particular parallel region overrides the current
        // number of threads.
        std::vector<std::vector<std::string>> filtered_regions(
            omp_get_max_threads());

#pragma omp parallel for schedule(static)
        for (int i = start; i < end; ++i) {
            int thread_id = omp_get_thread_num();

            const auto &pair = accession_hits_vec[i];
            auto accession = pair.first;
            auto partitioned_hits = pair.second;
            auto amplicon_regions = allAmpliconRegions(partitioned_hits);

            // Filter out any putative amplicons that aren't in the length
            // window
            std::vector<AmpliconRegion> filtered_amplicon_regions;
            for (const auto &region : amplicon_regions) {
                if (minimum_length <= region.product_length &&
                    region.product_length <= maximum_length) {
                    filtered_amplicon_regions.push_back(region);
                }
            }
            // Take the first, if there even is one...
            if (!filtered_amplicon_regions.empty()) {
                auto region = filtered_amplicon_regions.front();
                filtered_regions[thread_id].push_back(region.toTsvString());
            }

            // TODO: really should pick the best region? Or maybe do them all?
        }

        Rcpp::checkUserInterrupt();

        for (auto &chunk : filtered_regions) {
            for (auto &region : chunk) {
                output_stream << region << std::endl;
            }
        }
    }

    output_stream.close();
}
