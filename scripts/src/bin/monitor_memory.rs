use anyhow::{anyhow, Context, Result};
use clap::Parser;
use itertools::Itertools;
use jiff::civil::DateTime;
use jiff::tz::TimeZone;
use jiff::{Timestamp, Unit};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufWriter};
use std::path::PathBuf;
use std::process::{Command, Output};
use std::time;

const MISSING_TARGET_PID_THRESHOLD: u8 = 10;

fn main() {
    let cli = Cli::parse();

    match run(cli) {
        Ok(()) => (),
        Err(error) => {
            eprintln!("ERROR: {error}");
            std::process::exit(1);
        }
    }
}

fn run(cli: Cli) -> Result<()> {
    // TODO: would be nice to write to a gzip file for long running monitors, or else do some log rotation.
    let file = File::create(&cli.outfile)?;
    let mut writer = BufWriter::new(file);

    let mut missing_target_pid_count = 0;

    loop {
        std::thread::sleep(time::Duration::from_secs(1));
        let found_target_pid = process_ps(&mut writer, &cli)?;

        // Anytime we find the target PID, we reset the count
        if found_target_pid {
            missing_target_pid_count = 0;
        } else {
            missing_target_pid_count += 1;
        }

        // If we haven't found the target PID in the last 10 tries (about ten seconds),
        // then that process is probably gone, so stop.
        if missing_target_pid_count > MISSING_TARGET_PID_THRESHOLD {
            break;
        }
    }

    Ok(())
}

struct ProcessRecord {
    timestamp: Timestamp,
    monitored_pid: u32,
    uid: u32,
    pid: u32,
    lstart: DateTime,
    rss: u64,
    vsz: u64,
    comm: String,
    args: String,
}

impl ProcessRecord {
    fn to_tsv_string(&self) -> Result<String> {
        let format = "%Y-%M-%dT%H:%M:%S";
        let timestamp = self.timestamp.to_zoned(TimeZone::system()).strftime(format);
        let lstart = self.lstart.to_zoned(TimeZone::system())?.strftime(format);
        let formatted = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            timestamp,
            self.monitored_pid,
            self.uid,
            self.pid,
            lstart,
            self.rss,
            self.vsz,
            self.comm,
            self.args
        );

        Ok(formatted)
    }

    /// The timestamp represents the time when the ps command was run. You should get that before
    /// running `ps`, then pass it in here.
    ///
    fn parse(line: &str, timestamp: &Timestamp, monitored_pid: u32) -> Result<ProcessRecord> {
        let mut tokens = line.split_whitespace();

        // 101 12345 Thu Dec 11 10:30:14 2025 18720 426993376 blah blah -i something -o something_else

        let uid = tokens.next().context("no uid")?.parse::<u32>()?;
        let pid = tokens.next().context("no pid")?.parse::<u32>()?;
        let _weekday = tokens.next().context("no weekday")?;
        let month = tokens.next().context("no month")?;
        let day = tokens.next().context("no day")?.parse::<i8>()?;
        let hms = tokens.next().context("no time")?;
        let year = tokens.next().context("no year")?.parse::<i16>()?;
        let rss = tokens.next().context("no rss")?.parse::<u64>()?;
        let vsz = tokens.next().context("no vsz")?.parse::<u64>()?;
        let comm = tokens.next().context("no comm")?.to_string();
        let args = tokens.join(" ");

        let month = match month.to_lowercase().as_ref() {
            "jan" => 1,
            "feb" => 2,
            "mar" => 3,
            "apr" => 4,
            "may" => 5,
            "jun" => 6,
            "jul" => 7,
            "aug" => 8,
            "sep" => 9,
            "oct" => 10,
            "nov" => 11,
            "dec" => 12,
            bad_month => return Err(anyhow!("invalid month, got {bad_month}")),
        };

        let mut hms_tokens = hms.split(":");
        let hour = hms_tokens.next().context("no hours")?.parse::<i8>()?;
        let minute = hms_tokens.next().context("no minutes")?.parse::<i8>()?;
        let second = hms_tokens.next().context("no seconds")?.parse::<i8>()?;

        let lstart = jiff::civil::date(year, month, day).at(hour, minute, second, 0);

        let record = ProcessRecord {
            timestamp: *timestamp,
            monitored_pid,
            uid,
            pid,
            lstart,
            rss,
            vsz,
            comm,
            args,
        };

        Ok(record)
    }
}

fn run_ps() -> Result<Output> {
    let output = Command::new("ps")
        // ps ax -w -w -o uid= -o pid= -o lstart= -o rss= -o vsz= -o comm= -o args=
        .args([
            "ax", "-w", "-w", "-o", "uid=", "-o", "pid=", "-o", "lstart=", "-o", "rss=", "-o",
            "vsz=", "-o", "comm=", "-o", "args=",
        ])
        .output()?;

    Ok(output)
}

fn process_ps(writer: &mut BufWriter<File>, cli: &Cli) -> Result<bool> {
    let now = Timestamp::now().round(Unit::Second)?;
    let ps_output = run_ps()?;

    if !ps_output.status.success() {
        return Err(anyhow!(
            "bad exit status for ps command: {}\nstderr:\n{}\n",
            ps_output.status,
            String::from_utf8(ps_output.stderr)?
        ));
    }

    let mut found_target_pid = false;

    for line in ps_output.stdout.lines() {
        let line = line?;
        let record = ProcessRecord::parse(&line, &now, cli.pid)?;

        if record.pid == cli.pid {
            found_target_pid = true;
        }

        if record.uid == cli.uid {
            let output_line = record.to_tsv_string()?;
            writeln!(writer, "{output_line}")?;
        }
    }

    Ok(found_target_pid)
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
/// Track some usage stats for the current user's processes
///
/// When the process with the given PID is no longer running,
/// this program will terminate.
///
struct Cli {
    /// We will only keep processes owned by this user ID
    ///
    #[arg(short, long)]
    uid: u32,

    /// We will shut down when this process ID is gone
    ///
    #[arg(short, long)]
    pid: u32,

    /// Where to write the output?
    ///
    #[arg(short, long)]
    outfile: PathBuf,
}
