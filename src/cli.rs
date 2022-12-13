use clap::Parser;


#[derive(Parser)]
#[command(name = "QuantumSimulator")]
#[command(author = "Michael Tauraso <mtauraso@uw.edu>")]
#[command(version = "0.1.0")]
#[command(about = "Simulates a Quantum Circuit defined in an OpenQASM 2.0 file", long_about = None)]
#[command(next_line_help = false)]
pub struct Cli {
    /// The path to the OpenQASM 2.0 file specifying the circuit.
    #[arg(default_value_t = String::from("test.qasm"))]
    pub file: String,

    /// Print out amplitudes at each step of execution
    #[arg(long, short)]
    pub trace: bool,

    /// Print amplitude states in debug mode. Bits are printed in the order held internally ignoring registers
    #[arg(long, short)]
    pub debug: bool,
}
