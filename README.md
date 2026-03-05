<!-- TOC --><a name="rcruxmini"></a>

# rCRUXMini

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Setup](#setup)
  - [Get the code](#get-the-code)
  - [Install Rust Scripts](#install-rust-scripts)
  - [Create Taxonomy DB](#create-taxonomy-db)
  - [Config Files](#config-files)
    - [Required Top-Level Options](#required-top-level-options)
    - [Optional Top-Level Options](#optional-top-level-options)
    - [Primer & Amplicon BLAST Options](#primer-amplicon-blast-options)
      - [Primer BLAST Defaults](#primer-blast-defaults)
      - [Amplicon BLAST Defaults](#amplicon-blast-defaults)
    - [Plausible Amplicons Options](#plausible-amplicons-options)
      - [Defaults](#defaults)
    - [Example](#example)
  - [Run the rCRUXMini CLI](#run-the-rcruxmini-cli)
- [rCRUXMini Output Files](#rcruxmini-output-files)

<!-- TOC end -->

<!-- TOC --><a name="setup"></a>

## Setup

The basic flow is this:

- Get the code from GitHub: https://github.com/mooreryan/rCRUXMini
- Install Rust scripts
- Create a config file
- Create taxonomy DB
- Run CLI program

<!-- TOC --><a name="get-the-code"></a>

### Get the code

For now, you must track the main branch on GitHub. The main branch should always be in a working state, minus any possible bugs of course. That is, I don't push code that is known to be broken or unfinished to the main branch. This may not be the case for feature branches, however!

Using git:

```
git clone https://github.com/mooreryan/rCRUXMini.git
```

Or manually download [from this link](https://github.com/mooreryan/rCRUXMini/archive/refs/heads/main.zip)

_Note: If you would like to include a version number for rCRUXMini, say for like a manuscript or something, put the commit hash._

<!-- TOC --><a name="install-rust-scripts"></a>

### Install Rust Scripts

rCRUXMini has some parts that require compilation. For this, you will need a Rust toolchain installed on your computer/server. Check out the [Install Rust](https://rust-lang.org/tools/install/) page for how to set up Rust for your particular OS.

Once you have done that, you should be able to run `cargo --version` without error. If that is going well, you can install the scripts:

```
cd <path/to/rCRUXMini/repository>/scripts
cargo install --path . --root /usr/local --force
```

- You should replace `<path/to/rCRUXMini/repository>` with whatever you need to get into the rCRUXMini code repository that you downloaded earlier.
- You should replace `/usr/local` with wherever you want the execuatble scripts to live.

Alternatively, you can run `cargo build --release`, which will install the scripts to `<path/to/rCRUXMini/repository>/scripts/target/build/release`.

<!-- TOC --><a name="create-taxonomy-db"></a>

### Create Taxonomy DB

First, you will need to download the NCBI taxonomy taxdump and accession2taxid files. You can find them here:

- https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
  - https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
  - https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
- https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
  - https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
  - https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5

That FTP site has nice README files as well that can help you out.

Once you grab the required files, and put them somewhere convenient, then you will need to run a bit of R code to set up the DB in such a way that the taxonomizr package can read it. Something like this should do it:

```R
#!/usr/bin/env Rscript
library(taxonomizr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if all required arguments are provided
if (length(args) != 4) {
  stop(
    "Usage: make_taxonomizr_db.R <names.dmp> <nodes.dmp> <accession2taxid> <output.db>",
    call. = FALSE
  )
}

names_file <- args[1]
nodes_file <- args[2]
accession_file <- args[3]
sql_file <- args[4]

read.names.sql(nameFile = names_file, sqlFile = sql_file)
read.nodes.sql(nodeFile = nodes_file, sqlFile = sql_file)
read.accession2taxid(taxaFiles = accession_file, sqlFile = sql_file)
```

_Note: The so-called `output.db` from this little R script will be what you need to provide to the `taxonomy_database` option in the config file._

<!-- TOC --><a name="config-files"></a>

### Config Files

<!-- TOC --><a name="required-top-level-options"></a>

#### Required Top-Level Options

- `forward_primers`: one or more forward primers
- `reverse_primers`: one or more reverse primers
- `output_directory`: where rCRUXMini will put results
- `taxonomy_database`: path to the taxonomy DB
- `blast_databases`: one or more paths to (nucleotide) BLAST DBs to search

<!-- TOC --><a name="optional-top-level-options"></a>

#### Optional Top-Level Options

- `workers`
  - How many jobs to run in parallel?
  - More jobs (up to a point) mean less wall time
  - More jobs mean more memory (because multiple things will run at once)
  - _Default: `1`_
- `query_chunk_count`
  - Control the amount of "splits" queries will be broken up into
  - `5` would mean split the queries into five individual files to search
  - Increasing this value will keep memory down, since BLAST memory usage is also influenced by the size of the query set.
  - _Default: `1`_
- `ncbi_bin_directory`
  - Path to location of blastn, blastdbcmd, etc.
  - If you don't provide this, we assume they are on your PATH.
- `scripts_bin_directory`
  - Path to the location of rCRUXMini helper scripts.
  - If you don't provide this, we assume they are on your PATH.
- `primer_blast`: options passed to blastn for the primer blast
- `amplicon_blast`: options passed to blastn for the amplicon blast
- `plausible_amplicons`: options controlling filtering of plausible amplicons

<!-- TOC --><a name="primer-amplicon-blast-options"></a>

#### Primer & Amplicon BLAST Options

The `primer_blast` and `amplicon_blast` sections of the config file take the same options.

- `evalue`
- `num_alignments`
- `num_threads`
- `perc_identity`
- `qcov_hsp_perc`
- `reward`
- `task`
- `word_size`

These are all optional. That is, you can provide any or all of them.

<!-- TOC --><a name="primer-blast-defaults"></a>

##### Primer BLAST Defaults

- `evalue`: 3e7,
- `num_alignments`: 10000000,
- `num_threads`: 1,
- `perc_identity`: 50,
- `qcov_hsp_perc`: 90,
- `reward`: 2,
- `task`: "blastn-short",
- `word_size`: 7

_Note: You should probably set these yourself rather than relying on the defaults!_

<!-- TOC --><a name="amplicon-blast-defaults"></a>

##### Amplicon BLAST Defaults

- `evalue`: 3e7,
- `num_alignments`: 10000000,
- `num_threads`: 1,
- `task`: "megablast"

_Note: You should probably set these yourself rather than relying on the defaults!_

<!-- TOC --><a name="plausible-amplicons-options"></a>

#### Plausible Amplicons Options

- `minimum_length`: valid amplicons must be at least this long
- `maximum_length`: valid amplicons must be no more than this long
- `maximum_mismatches`: valid amplicons must have no more than this many mismatches the BLAST search
- `ambiguous_run_limit`: valid amplicons must no more than than this many ambiguous bases in a row

<!-- TOC --><a name="defaults"></a>

##### Defaults

- `minimum_length`: 150
- `maximum_length`: 650
- `maximum_mismatches`: 4
- `ambiguous_run_limit`: 5

_Again, you should set these to whatever makes sense for your amplicon!_

<!-- TOC --><a name="example"></a>

#### Example

Here is an example config file.

```yaml
# Each of these can be a list if you have more than one.
forward_primers: "ANTG"
reverse_primers:
  - "AACCTTGG"
  - "CCTTAAGG"

# Where rCRUXMini will put the output files
output_directory: "/home/ryan/rcrux_output"

# Path to your taxonomy DB
taxonomy_database: "/home/ryan/taxonomy.db"

# Path to location of blastn, blastdbcmd, etc. If you don't provide this, we
# assume they are on your PATH.
ncbi_bin_directory: "/opt/ncbi/bin"

# Path to the location of rCRUXMini helper scripts. If you don't provide this,
# we assume they are on your PATH.
scripts_bin_directory: "/opt/rCRUXMini/scripts/target/release"

# Path or paths of blast DBs to search. You can split up a BLAST DB to keep the
# memory down. It will increase the runtime (likely) and affect evalues.
blast_databases:
  - "/opt/ncbi/blastdb/nt.001"
  - "/opt/ncbi/blastdb/nt.002"
  # As many as you want!

# If you just have one blast DB, you can do this:
# blast_databases: "/path/to/my/db"

# We can split queries to keep memory down. BLAST memory usage also scales with
# the queries in addition to the subjects.
query_chunk_count: 2

# How many jobs to run at once. Higher numbers will use less wall time, but will
# also use more memory, since more than one job may be running at one time.
workers: 4

# These options are passed to the blastn command for the primer blast step.
primer_blast:
  evalue: 3e7
  num_alignments: 10000000
  num_threads: 1
  perc_identity: 50
  qcov_hsp_perc: 90
  reward: 2
  task: "blastn-short"
  word_size: 7

# These options are passed to the blastn command for the primer blast step. You
# can use the same options here as for the primer_blast.
amplicon_blast:
  evalue: 3e7
  num_alignments: 10000000
  num_threads: 1
  perc_identity: 50
  qcov_hsp_perc: 90
  task: "megablast"

# What do you consider plausible amplicons? You can filter by length, max
# mismatches in the blast, max number of ambiguious bases (eg, NNNNN).
#
# You should base these on things you know about your specific amplicon.
plausible_amplicons:
  minimum_length: 150
  maximum_length: 650
  maximum_mismatches: 4
  ambiguous_run_limit: 5
```

<!-- TOC --><a name="run-the-rcruxmini-cli"></a>

### Run the rCRUXMini CLI

Now, it's finally time to run rCRUXMini!

```
time Rscript --vanilla \
<path to rCRUXMini repo>/inst/cli/rCRUXMini.R \
<path to config file>/config.yml \
<path to rCRUXMini repo>
```

- Be sure to replace `<path to rCRUXMini repo>` with the actual location where you downloaded rCRUXMini repository.
- Be sure to replace `<path to config file>` with the actual location of your config file.

So, you might have something that looks like this:

```
time Rscript --vanilla \
/opt/rCRUXMini/inst/cli/rCRUXMini.R \
/home/ryan/cool_project/rcrux_config.yml \
/opt/rCRUXMini
```

<!-- TOC --><a name="rcruxmini-output-files"></a>

## rCRUXMini Output Files

Guide coming soon!
