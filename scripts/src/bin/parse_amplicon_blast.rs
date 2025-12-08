use ahash::{HashMap, HashMapExt};
use clap::Parser;
use rayon::prelude::*;
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    run(cli)
}

/// Parse rCRUXMini amplicon BLAST results
///
/// Handles the first couple of steps at least...
///
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Write single taxonomy hits here
    ///
    single_taxonomy_hits_out: PathBuf,

    /// Write multi taxonomy hits here
    ///
    multi_taxonomy_hits_out: PathBuf,

    /// Output files from rCRUXMini amplicon BLAST
    ///
    blast_files: Vec<PathBuf>,

    /// Number of threads to use
    ///
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
}

#[derive(Clone)]
struct Record {
    query_accession: String,
    subject_accession_version: String,
    percent_identical_matches: f32,
    alignment_length: u64,
    expect_value: f32,
    subject_sequence_length: u64,
    subject_alignment_start: u64,
    subject_alignment_end: u64,
    subject_aligned_sequence: String,
    unique_subject_taxonomy_ids: String,
    // These are derived fields
    degapped_subject_aligned_sequence: String,
    degapped_alignment_length: usize,
}

const SUBJECT_ACCESSION_VERSION_INDEX: usize = 1;
const ALIGNMENT_LENGTH_INDEX: usize = 3;
const SUBJECT_ALIGNED_SEQUENCE_INDEX: usize = 8;

impl Record {
    fn has_multiple_taxonomy_ids(&self) -> bool {
        self.unique_subject_taxonomy_ids.contains(';')
    }

    fn from_fields(fields: &[&str]) -> Self {
        let subject_aligned_sequence = fields[SUBJECT_ALIGNED_SEQUENCE_INDEX].to_string();
        let degapped_subject_aligned_sequence = subject_aligned_sequence.replace('-', "");
        let degapped_alignment_length = degapped_subject_aligned_sequence.len();

        Record {
            query_accession: fields[0].to_string(),
            subject_accession_version: fields[SUBJECT_ACCESSION_VERSION_INDEX].to_string(),
            percent_identical_matches: fields[2]
                .parse()
                .expect("failed to parse percent_identical_matches"),
            alignment_length: fields[ALIGNMENT_LENGTH_INDEX]
                .parse()
                .expect("failed to parse alignment_length"),
            expect_value: fields[4].parse().expect("failed to parse expect_value"),
            subject_sequence_length: fields[5]
                .parse()
                .expect("failed to parse subject_sequence_length"),
            subject_alignment_start: fields[6]
                .parse()
                .expect("failed to parse subject_alignment_start"),
            subject_alignment_end: fields[7]
                .parse()
                .expect("failed to parse subject_alignment_end"),
            subject_aligned_sequence,
            unique_subject_taxonomy_ids: fields[9].to_string(),
            degapped_subject_aligned_sequence,
            degapped_alignment_length,
        }
    }

    fn to_tsv_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_accession,
            self.subject_accession_version,
            self.percent_identical_matches,
            self.alignment_length,
            self.expect_value,
            self.subject_sequence_length,
            self.subject_alignment_start,
            self.subject_alignment_end,
            self.subject_aligned_sequence,
            self.unique_subject_taxonomy_ids,
            self.degapped_subject_aligned_sequence,
            self.degapped_alignment_length,
        )
    }

    fn header_tsv_string() -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            "query_accession",
            "subject_accession_version",
            "percent_identical_matches",
            "alignment_length",
            "expect_value",
            "subject_sequence_length",
            "subject_alignment_start",
            "subject_alignment_end",
            "subject_aligned_sequence",
            "unique_subject_taxonomy_ids",
            "degapped_subject_aligned_sequence",
            "degapped_alignment_length",
        )
    }
}

fn run(cli: Cli) -> io::Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    let single_taxonomy_hits_file = File::create(cli.single_taxonomy_hits_out)?;
    let mut single_taxonomy_hits_writer = BufWriter::new(single_taxonomy_hits_file);

    let multi_taxonomy_hits_file = File::create(cli.multi_taxonomy_hits_out)?;
    let mut multi_taxonomy_hits_writer = BufWriter::new(multi_taxonomy_hits_file);

    // Go through all the blast files, filtering hits only keeping the hit with the longest alignment
    // length per subject accession.

    let best_hits = cli
        .blast_files
        .par_iter()
        // Only the Ok results will be processed by the reduce function
        .flat_map(|path| process_file(path))
        .reduce(HashMap::new, merge_best_hit_maps);

    // Write lines. Even empty files should have the expected header.
    writeln!(
        multi_taxonomy_hits_writer,
        "{}",
        Record::header_tsv_string()
    )?;
    writeln!(
        single_taxonomy_hits_writer,
        "{}",
        Record::header_tsv_string()
    )?;

    for (_subject_accession, record) in best_hits {
        if record.has_multiple_taxonomy_ids() {
            writeln!(multi_taxonomy_hits_writer, "{}", record.to_tsv_string())?;
        } else {
            writeln!(single_taxonomy_hits_writer, "{}", record.to_tsv_string())?;
        }
    }

    Ok(())
}

fn process_file(path: &Path) -> io::Result<HashMap<String, Record>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut best_hits = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if !line.starts_with("query_accession") {
            let fields: Vec<_> = line.split('\t').collect();

            let subject_accession_version = fields[SUBJECT_ACCESSION_VERSION_INDEX];
            let alignment_length = fields[ALIGNMENT_LENGTH_INDEX]
                .parse::<u64>()
                .expect("failed to parse alignment_length");

            match best_hits.entry(subject_accession_version.to_string()) {
                Entry::Occupied(mut entry) => {
                    let best_record: &Record = entry.get();

                    if alignment_length > best_record.alignment_length {
                        entry.insert(Record::from_fields(&fields));
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(Record::from_fields(&fields));
                }
            }
        }
    }

    Ok(best_hits)
}

fn merge_best_hit_maps(
    mut a: HashMap<String, Record>,
    b: HashMap<String, Record>,
) -> HashMap<String, Record> {
    // We always want to insert the smaller map into the larger one
    if b.len() > a.len() {
        return merge_best_hit_maps(b, a);
    }

    for (subject_accession, record_b) in b {
        match a.entry(subject_accession) {
            Entry::Occupied(mut entry_a) => {
                let record_a = entry_a.get();

                // Replace the record in A with the one in B
                // if the one in B has a longer alignment length
                if record_b.alignment_length > record_a.alignment_length {
                    entry_a.insert(record_b);
                }
            }
            Entry::Vacant(entry) => {
                entry.insert(record_b);
            }
        }
    }

    a
}
