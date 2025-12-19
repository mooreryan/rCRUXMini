use clap::Parser;
use scripts::monitor_memory::{run, Cli};

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
