use ahash::{HashMap, HashMapExt};
use clap::Parser;
use itertools::Itertools;
use jiff::Zoned;
use rayon::prelude::*;
use std::fs::File;
use std::io;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::ops::Deref;
use std::path::{Path, PathBuf};

fn log(message: &str) {
    let now = Zoned::now().strftime("%Y-%M-%d %H:%M:%S");

    writeln!(io::stderr(), "{} -- {}", now, message).expect("failed to log message");
}

fn main() {
    let cli = Cli::parse();
    assert!(
        cli.minimum_length <= cli.maximum_length,
        "min must be <= max"
    );
    assert!(!cli.blast_files.is_empty());

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    let out = File::create(cli.out_file).expect("could not create outfile");
    let mut writer = BufWriter::new(out);

    let mut forward_and_reverse_hits_grouped_by_saccver: HashMap<String, ForwardAndReverseHits> =
        HashMap::new();

    // NOTE: original code has a group by saccver+send and keep distinct...but
    // IDK if you can ever even get that since we are only ever keeping one of
    // the filtered hits anyway?

    log("processing files");
    let hits_grouped_by_saccver_sstart = cli
        .blast_files
        .par_iter()
        .map(|path| group_and_filter_hits(path, cli.maximum_mismatches))
        .reduce(HashMap::new, merge_grouped_hit_maps);

    // Within the saccver-sstart grouped hits, we need to keep the one with the lowest number of
    // mismatches. If there is more than one of those, the original code simply picks the first one.
    log("picking best hits");
    hits_grouped_by_saccver_sstart
        .values()
        .map(|hits| {
            hits.iter()
                .min_by(|a, b| a.mismatch.cmp(&b.mismatch))
                .expect("expected at least one primer blast hit")
        })
        .unique_by(|hit| (hit.saccver.clone(), hit.send))
        .for_each(|hit| {
            forward_and_reverse_hits_grouped_by_saccver
                .entry(hit.saccver.clone())
                .or_insert(ForwardAndReverseHits::empty())
                .add_primer_hit(hit.clone());
        });

    log("finding amplicon regions");
    writeln!(writer, "{}", AmpliconRegion::tsv_header()).expect("could not write line");
    let amplicon_regions: Vec<_> = forward_and_reverse_hits_grouped_by_saccver
        .par_iter()
        .filter(|(_, forward_and_reverse_hits)| {
            forward_and_reverse_hits.has_both_forward_and_reverse_hits()
        })
        .flat_map(|(_saccver, forward_and_reverse_hits)| {
            // For each subject accession, go through all of its putative amplicon regions and pick
            // the first one that satisfies the length requirement.
            //
            // TODO: This is the original behavior, but shouldn't we try and pick the "best"?
            forward_and_reverse_hits
                .iter_amplicon_regions()
                // Find takes the first one that satisfies the predicate
                .find(|amplicon_region| match amplicon_region.product_length() {
                    None => false,
                    Some(product_length) => {
                        cli.minimum_length <= product_length && product_length <= cli.maximum_length
                    }
                })
        })
        .collect();

    log("writing regions");
    for amplicon_region in amplicon_regions {
        writeln!(writer, "{}", amplicon_region.to_tsv_string()).expect("could not write line")
    }

    log("done!");
}

/// Process a given primer BLAST output file.
///
/// Does some filtering and then groups on subject accession version and subject start position.
///
/// We need to group the hits based on their accession and start location.
/// Additionally, we need to ignore any hits that have too many mismatches.
/// This needs to be done across all input files.
///
/// Crashes on io issues
///
fn group_and_filter_hits(
    path: &Path,
    maximum_mismatches: u64,
) -> HashMap<SaccverSstartPair, Vec<PrimerHit>> {
    let file = File::open(path).expect("failed to open file");
    let reader = BufReader::new(file);

    let mut saccver_sstart_grouped_hits: HashMap<SaccverSstartPair, Vec<PrimerHit>> =
        HashMap::new();

    for line in reader.lines() {
        let line = line.expect("failed to read line");

        if BlastHit::is_header_line(&line) {
            continue;
        }

        let hit = PrimerHit::from_line(&line);

        // NOTE: The original code had >= here, but that is incorrect.
        // TODO: change to >
        if hit.mismatch >= maximum_mismatches {
            continue;
        }

        let key = SaccverSstartPair::from_primer_blast_hit(&hit);
        saccver_sstart_grouped_hits
            .entry(key)
            .or_default()
            .push(hit);
    }

    saccver_sstart_grouped_hits
}

fn merge_grouped_hit_maps(
    mut a: HashMap<SaccverSstartPair, Vec<PrimerHit>>,
    mut b: HashMap<SaccverSstartPair, Vec<PrimerHit>>,
) -> HashMap<SaccverSstartPair, Vec<PrimerHit>> {
    // We always want to insert the smaller map into the larger one
    if b.len() > a.len() {
        std::mem::swap(&mut a, &mut b);
    }

    for (key, mut values) in b {
        a.entry(key)
            .and_modify(|e| e.append(&mut values))
            .or_insert(values);
    }

    a
}

/// Parse rCRUXMini primer BLAST results
///
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Write output here
    ///
    #[arg(long)]
    out_file: PathBuf,

    /// Number of threads to use
    ///
    #[arg(long, default_value_t = 1)]
    threads: usize,

    #[arg(long)]
    maximum_mismatches: u64,

    #[arg(long)]
    minimum_length: u64,

    #[arg(long)]
    maximum_length: u64,

    /// Output files from rCRUXMini primer BLAST
    ///
    blast_files: Vec<PathBuf>,
}

#[derive(Debug, Clone, Default)]
struct AmpliconRegion {
    saccver: String,
    sgi: String,
    staxids: String,
    forward_start: u64,
    forward_stop: u64,
    forward_mismatch: u64,
    reverse_start: u64,
    reverse_stop: u64,
    reverse_mismatch: u64,
}

impl AmpliconRegion {
    fn new(forward_hit: &ForwardPrimerHit, reverse_hit: &ReversePrimerHit) -> Self {
        let forward_hit = &forward_hit.0;
        let reverse_hit = &reverse_hit.0;

        assert_eq!(forward_hit.saccver, reverse_hit.saccver);
        assert_eq!(forward_hit.sgi, reverse_hit.sgi);
        assert_eq!(forward_hit.staxids, reverse_hit.staxids);

        AmpliconRegion {
            saccver: forward_hit.saccver.clone(),
            sgi: forward_hit.sgi.clone(),
            staxids: forward_hit.staxids.clone(),
            forward_start: forward_hit.sstart,
            forward_stop: forward_hit.send,
            forward_mismatch: forward_hit.mismatch,
            reverse_start: reverse_hit.sstart,
            reverse_stop: reverse_hit.send,
            reverse_mismatch: reverse_hit.mismatch,
        }
    }

    /// Calculate the product length of the amplicon region.
    ///
    /// Not all amplicon regions have valid product lengths. For example, the primers may be in the
    /// incorrect orientation.
    ///
    fn product_length(&self) -> Option<u64> {
        //  F ---------->                        <------------ R  Primers
        // ====================================================== DNA Target
        //    ---------->************************<------------    Product
        //    ^ Forward start                  Reverse start ^
        if self.forward_start < self.reverse_start
            && self.forward_start < self.forward_stop
            && self.reverse_stop < self.reverse_start
        {
            // Assuming there isn't a bug in the conditional logic above, this shouldn't ever
            // saturate.
            let product_length = self.reverse_start.saturating_sub(self.forward_start);
            return Some(product_length);
        }

        //  R ---------->                        <------------ F  Primers
        // ====================================================== DNA Target
        //    ---------->************************<------------    Product
        //    ^ Reverse start                  Forward start ^
        if self.forward_start > self.reverse_start
            && self.forward_start > self.forward_stop
            && self.reverse_stop > self.reverse_start
        {
            let product_length = self.forward_start.saturating_sub(self.reverse_start);
            return Some(product_length);
        }

        None
    }

    fn to_tsv_string(&self) -> String {
        let product_length = match self.product_length() {
            None => "NA".to_string(),
            Some(n) => n.to_string(),
        };

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.saccver,
            self.sgi,
            self.staxids,
            self.forward_start,
            self.forward_stop,
            self.forward_mismatch,
            self.reverse_start,
            self.reverse_stop,
            self.reverse_mismatch,
            product_length
        )
    }

    fn tsv_header() -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            "saccver",
            "sgi",
            "staxids",
            "forward_start",
            "forward_stop",
            "forward_mismatch",
            "reverse_start",
            "reverse_stop",
            "reverse_mismatch",
            "product_length"
        )
    }
}

#[derive(Eq, PartialEq, Clone, Hash)]
struct BlastHit {
    qseqid: String,
    sgi: String,
    saccver: String,
    mismatch: u64,
    sstart: u64,
    send: u64,
    staxids: String,
}

impl BlastHit {
    fn is_header_line(line: &str) -> bool {
        line.starts_with("qseqid")
    }
}

#[derive(Clone)]
struct ForwardPrimerHit(BlastHit);

impl Deref for ForwardPrimerHit {
    type Target = BlastHit;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Clone)]
struct ReversePrimerHit(BlastHit);

impl Deref for ReversePrimerHit {
    type Target = BlastHit;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Clone)]
enum PrimerHit {
    Forward(ForwardPrimerHit),
    Reverse(ReversePrimerHit),
}

impl Deref for PrimerHit {
    type Target = BlastHit;
    fn deref(&self) -> &Self::Target {
        match self {
            PrimerHit::Forward(hit) => &hit.0,
            PrimerHit::Reverse(hit) => &hit.0,
        }
    }
}

impl PrimerHit {
    fn from_line(line: &str) -> Self {
        let fields: Vec<_> = line.split('\t').collect();
        Self::from_fields(&fields)
    }

    fn from_fields(fields: &[&str]) -> Self {
        assert_eq!(
            fields.len(),
            7,
            "Expected 7 fields in BLAST output, got {}",
            fields.len()
        );

        // TODO: the original code checks for whitespace only fields.

        let blast_hit = BlastHit {
            qseqid: fields[0].to_string(),
            sgi: fields[1].to_string(),
            saccver: fields[2].to_string(),
            mismatch: fields[3].parse().expect("failed to parse mismatch"),
            sstart: fields[4].parse().expect("failed to parse sstart"),
            send: fields[5].parse().expect("failed to parse send"),
            staxids: fields[6].to_string(),
        };

        if blast_hit.qseqid.starts_with("forward") {
            PrimerHit::Forward(ForwardPrimerHit(blast_hit))
        } else if blast_hit.qseqid.starts_with("reverse") {
            PrimerHit::Reverse(ReversePrimerHit(blast_hit))
        } else {
            // If we hit this, then there is a bug in the rCRUXMini R code.
            panic!("expected either forward or reverse primer")
        }
    }
}

struct ForwardAndReverseHits {
    forward_hits: Vec<ForwardPrimerHit>,
    reverse_hits: Vec<ReversePrimerHit>,
}

impl ForwardAndReverseHits {
    fn has_both_forward_and_reverse_hits(&self) -> bool {
        !self.forward_hits.is_empty() && !self.reverse_hits.is_empty()
    }

    fn add_primer_hit(&mut self, primer_hit: PrimerHit) {
        match primer_hit {
            PrimerHit::Forward(hit) => {
                self.forward_hits.push(hit);
            }
            PrimerHit::Reverse(hit) => {
                self.reverse_hits.push(hit);
            }
        }
    }

    fn empty() -> Self {
        ForwardAndReverseHits {
            forward_hits: Vec::new(),
            reverse_hits: Vec::new(),
        }
    }

    fn iter_amplicon_regions(&self) -> impl Iterator<Item = AmpliconRegion> {
        self.forward_hits.iter().flat_map(|forward_hit| {
            self.reverse_hits
                .iter()
                .map(|reverse_hit| AmpliconRegion::new(forward_hit, reverse_hit))
        })
    }
}

#[derive(Eq, PartialEq, Hash)]
struct SaccverSstartPair {
    saccver: String,
    sstart: u64,
}

impl SaccverSstartPair {
    fn from_primer_blast_hit(primer_hit: &PrimerHit) -> Self {
        SaccverSstartPair {
            saccver: primer_hit.saccver.clone(),
            sstart: primer_hit.sstart,
        }
    }
}
