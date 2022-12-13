Quantum Simulator is a simple command line application which can execute a program written in OpenQASM 2.0 and output the probabilities associated with measurements made as part of the program. The underlying quantum simulation is state vector based and does not incorporate sources of error arising from a physical system; however it does support all operations defined in the [OpenQASM 2.0 spec](https://arxiv.org/pdf/1707.03429.pdf)

Detail on the implementation choices and mathematical formalism can be found in [the paper](https://mtauraso.github.io/QuantumSimulator/LaTeX/QuantumSimulator.pdf)

The source code can be found [here](https://github.com/mtauraso/QuantumSimulator)

# Getting Started

## Installation
### Windows
The fastest way to get started is to download `quantum_simulator.exe` as well as the files in the `sample` directory. For most programs you will need the program file and also the file `qelib1.inc` for gate definitions. You can then open windows terminal, navigate to your download location and run `quantum_simulator.exe --help` You should see a help message.

### Other systems (OSX, Linux)
`git clone` the [source code](https://github.com/mtauraso/QuantumSimulator) and run `cargo build`. This should generate an executable which you can run directly from `/target/debug/quantum_simulator` or you can use `cargo run`. If you are using `cargo run` and want to pass command line arguments you can use `cargo run --` followed by the command line arguments for `quantum_simulator` e.g. `cargo run -- --help` will print the help messages

## Running programs
From your terminal type `quantum_simulator.exe <filename>` to run a qasm file. You should get output showing the probabilities associated with any measurements done in the program. The image below shows how this looks

![Gif showing quantum_simulator being run on a simple test program](https://mtauraso.github.io/QuantumSimulator/images/simple.gif)

You can also run quantum simulator with the `--trace` option, which will enumerate the steps in the circuit, printing out the amplitudes after each step. Here's another video showing the process:

![Gif showing quantum_simulator being run on a simple test program in trace mode](https://mtauraso.github.io/QuantumSimulator/images/trace.gif)

Note that the amplitudes are printed with the same names defined in the program. The simulator also keeps track of hidden bits, which are detailed in [the paper](https://mtauraso.github.io/QuantumSimulator/LaTeX/QuantumSimulator.pdf)

# Sample circuits
Sample circuits are in the `sample` directory. A few basic circuits are there for ease of use, many of them are originally from the OpenQASM spec with only slight modifications.


Note:Deliverables would be program, source code, readme & getting started level documentation, screen-capture demo videos, and runnable sample circuits that the simulator can execute put up on a website. 