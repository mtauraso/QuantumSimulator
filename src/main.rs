//extern crate openqasm as oq;
//extern crate ndarray as nd;
//extern crate num;

mod register;
mod check;

use openqasm as oq;

use oq::{GenericError,ProgramVisitor, translate::Linearize};
use register::{QuantumRegister, QuantumRegisterGateWriter};

use check::CheckGateWriter;

#[cfg(test)]
use std::path::Path;


fn main() {
    let program = openqasm_parse_file("test.qasm");
    let register = openqasm_run_program(program);

    
    // Print the final state
    println!("Final Amplitudes: {}", register.to_string());

    println!("Probability second qubit is true: {:.3}",register.probability(vec![None, Some(true), None]));

    println!("Should print 1 if we're still conserving probability {:.3}", register.probability(vec![None, None, None]));


    println!("Probability second qubit is true given the first qubit is also true: {:.3}", 
        register.probability(vec![Some(true), Some(true), None]) / 
        register.probability(vec![Some(true), None, None])
    );

    println!("Probability all bits are zero: {:.3}", register.probability(vec![Some(false), Some(false), Some(false)]));

}

pub fn openqasm_check_program(parser: oq::Parser) -> Result<oq::Program,oq::Errors> {
    let program = match parser.done().to_errors() {
        Ok(program) => program,
        Err(errors) => return Err(errors),
    };

    if let Err(errors) = program.type_check().to_errors() {
        return Err(errors);
    }

    if let Err(error) = Linearize::new(CheckGateWriter{}).visit_program(&program) {
        return Err(
            oq::Errors{ 
                errors: vec![error.into()]
            });
    }
    Ok(program)
}

pub fn openqasm_parse_file(path: &str) -> oq::Program{
    let mut cache = oq::SourceCache::new();
    let mut parser = oq::Parser::new(&mut cache);

    parser.parse_file(path);

    match openqasm_check_program(parser) {
        Ok(program) => program,
        Err(errors) => {
            errors.print(&mut cache).unwrap();
            panic!("Cannot continue due to errors above");
        }
    }
}

#[cfg(test)]
fn openqasm_parse_source<P: AsRef<Path>>(source: String, path: Option<P>) -> oq::Program{
    let mut cache = oq::SourceCache::new();
    let mut parser = oq::Parser::new(&mut cache);
    
    parser.parse_source(source, path);

    match openqasm_check_program(parser) {
        Ok(program) => program,
        Err(errors) => {
            errors.print(&mut cache).unwrap();
            panic!("Cannot continue due to errors above");
        }
    }
}


fn openqasm_run_program (program: oq::Program) -> QuantumRegister {
    // Run the program
    let mut register = QuantumRegister::new();
    let mut linearizer = Linearize::new(
        QuantumRegisterGateWriter::new(&mut register));
    linearizer.visit_program(&program).to_errors().unwrap();

    register
}

/// Print any sort of errors (ugly but functional if were using openqasm without ariadne)
fn _openqasm_print_errors(errors: oq::Errors) {
    for error in errors.errors {
        println!("{:?}", error);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use approx::*;
    use ndarray::prelude::*;

    #[test]
    fn state() {
        let source = r#"
        OPENQASM 2.0;

        // Standard header for usual gate set
        include "qelib1.inc";

        qreg q[1];
        creg c[1];

        // Need to have at least one statement otherwise OpenQASM translation layer never initializes state
        barrier;
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);

        // Newly created cbits and qubits should all be zero
        assert_relative_eq!(register.probability(vec![Some(false), Some(false)]), 1.0);
        assert_relative_eq!(register.probability(vec![Some(false), Some(true)]), 0.0);
        assert_relative_eq!(register.probability(vec![Some(true), Some(false)]), 0.0);
        assert_relative_eq!(register.probability(vec![Some(true), Some(true)]), 0.0);

        // All probabilities should add to 1.0
        assert_relative_eq!(register.probability(vec![None, None]), 1.0);
    }

    #[test]
    fn hadamard() {
        let source = r#"
        OPENQASM 2.0;

        // Standard header for usual gate set
        include "qelib1.inc";

        qreg q[1];
        h q;
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);
        
        // Each state of the qubit should have equal probability
        assert_relative_eq!(register.probability(vec![Some(true)]), 0.5);
        assert_relative_eq!(register.probability(vec![Some(false)]), 0.5);

        // Probabilities should sum to one
        assert_relative_eq!(register.probability(vec![None]), 1.0);
    }

    #[test]
    fn cnot() {
        let source = r#"
        OPENQASM 2.0;

        // Standard header for usual gate set
        include "qelib1.inc";

        qreg q[2];
        h q[0];
        cx q[0], q[1];
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);
        assert_relative_eq!(register.probability(vec![Some(true), Some(false)]), 0.0);
        assert_relative_eq!(register.probability(vec![Some(true), Some(true)]), 0.5);
        assert_relative_eq!(register.probability(vec![Some(false), Some(false)]), 0.5);
        assert_relative_eq!(register.probability(vec![None]), 1.0);
    }

    #[test]
    fn measure() {
        let source = r#"
        OPENQASM 2.0;

        // Standard header for usual gate set
        include "qelib1.inc";

        qreg q[2];
        creg c[1];
        h q;
        cx q[0], q[1];
        
        measure q[0] -> c[0];
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);
        

        // The cbit should have a 50/50 chance of being set
        assert_relative_eq!(register.probability(vec![None, None, Some(true)]), 0.5);
        assert_relative_eq!(register.probability(vec![None, None, Some(false)]), 0.5);

        //Each bit pattern 00, 01, 10, 11 in the qubits should have 25% prob
        assert_relative_eq!(register.probability(vec![Some(false), Some(false), None]), 0.25);
        assert_relative_eq!(register.probability(vec![Some(false), Some(true), None]), 0.25);
        assert_relative_eq!(register.probability(vec![Some(true), Some(false), None]), 0.25);
        assert_relative_eq!(register.probability(vec![Some(true), Some(true), None]), 0.25);

        // Probabilities should add to 1.0
        assert_relative_eq!(register.probability(vec![None, None, None]), 1.0);
    }

    #[test]
    fn measure_amplitudes() {
        let source = r#"
        OPENQASM 2.0;

        // Standard header for usual gate set
        include "qelib1.inc";

        qreg q[3];
        creg c[1];
        h q[0];
        h q[1];
        cx q[1], q[2];
        z q[2];
        cx q[2], q[0];
        cx q[0], q[2];
        measure q[1] -> c[0];
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);

        println!("{:?}",register.amplitudes());

        let amps = register.amplitudes().to_owned().into_dimensionality::<Ix4>().unwrap();

        // bit patterns 110 and 011 should have negative imaginary amplitudes of 0.5
        // These are the situations where q[1] was measured to be 1 so c[0] is also 1
        assert_abs_diff_eq!(amps[[1,1,0,1]].im, -0.5, epsilon = 0.0001);
        assert_abs_diff_eq!(amps[[0,1,1,1]].im, -0.5, epsilon = 0.0001);

        // bit patterns 000 and 101 should have positive imaginary amplitudes of 0.5
        // These are the situations where q[2] was 0 before reset
        assert_abs_diff_eq!(amps[[1,0,1,0]].im, 0.5, epsilon = 0.0001);
        assert_abs_diff_eq!(amps[[0,0,0,0]].im, 0.5, epsilon = 0.0001);

        //Each bit pattern 1101, 0111, 1010, 0000 in all the bits should have 25% prob
        assert_abs_diff_eq!(register.probability(vec![Some(true), Some(true), Some(false), Some(true)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(false), Some(true), Some(true), Some(true)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(true), Some(false), Some(true), Some(false)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(false), Some(false), Some(false), Some(false)]), 0.25, epsilon = 0.0001);

        // Probabilities should add to 1.0
        assert_relative_eq!(register.probability(vec![None, None, None]), 1.0, epsilon = 0.0001);
    }

    #[test]
    fn measure_overwrite() {
        let source = r#"
        OPENQASM 2.0;
        include "qelib1.inc";

        qreg q[3];
        creg c[1];
        h q[0];
        h q[1];
        measure q[1] -> c[0];
        ccx q[0], q[1], q[2];
        h q[2];
        measure q[2] -> c[0];
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);

        println!("{:?}",register.amplitudes());

        //let amps = register.amplitudes().to_owned().into_dimensionality::<Ix3>().unwrap();

        // // bit patterns 110 and 011 should have negative imaginary amplitudes of 0.5
        // // These are the situations where q[1] was measured to be 1 so c[0] is also 1
        // assert_abs_diff_eq!(amps[[1,1,0,1]].im, -0.5, epsilon = 0.0001);
        // assert_abs_diff_eq!(amps[[0,1,1,1]].im, -0.5, epsilon = 0.0001);

        // // bit patterns 000 and 101 should have positive imaginary amplitudes of 0.5
        // // These are the situations where q[2] was 0 before reset
        // assert_abs_diff_eq!(amps[[1,0,1,0]].im, 0.5, epsilon = 0.0001);
        // assert_abs_diff_eq!(amps[[0,0,0,0]].im, 0.5, epsilon = 0.0001);


        // assert_abs_diff_eq!(register.probability(vec![Some(true), Some(true), Some(true)]), 0.25, epsilon = 0.0001);
        // assert_abs_diff_eq!(register.probability(vec![Some(false), Some(false), Some(true)]), 0.25, epsilon = 0.0001);
        // assert_abs_diff_eq!(register.probability(vec![Some(true), Some(true), Some(false)]), 0.25, epsilon = 0.0001);
        // assert_abs_diff_eq!(register.probability(vec![Some(false), Some(false), Some(false)]), 0.25, epsilon = 0.0001);

        // Probabilities should add to 1.0
        assert_relative_eq!(register.probability(vec![None, None, None]), 1.0, epsilon = 0.0001);
    }


    #[test]
    fn reset() {
        let source = r#"
        OPENQASM 2.0;

        // Standard header for usual gate set
        include "qelib1.inc";

        qreg q[3];
        h q[0];
        h q[1];
        ccx q[0], q[1], q[2];
        reset q[2];
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);

        // This should have resulted in a single hidden qubit
        assert_eq!(register.bit_count() - register.visible_bit_count(), 1);

        //Each bit pattern 00, 01, 10, 11 in the qubits should have 25% prob
        assert_abs_diff_eq!(register.probability(vec![Some(true), Some(true), Some(false)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(false), Some(false), Some(false)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(true), Some(false), Some(false)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(false), Some(true), Some(false)]), 0.25, epsilon = 0.0001);

        // Probabilities should add to 1.0
        assert_relative_eq!(register.probability(vec![None, None, None]), 1.0, epsilon = 0.0001);
    }

    #[test]
    fn reset_amplitude_signs() {
        let source = r#"
        OPENQASM 2.0;

        // Standard header for usual gate set
        include "qelib1.inc";

        qreg q[3];
        h q[0];
        h q[1];
        cx q[1], q[2];
        z q[2];
        reset q[2];
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);

        println!("{:?}",register.amplitudes());

        // This should have resulted in one hidden bit
        assert_eq!(register.bit_count() - register.visible_bit_count(), 1);

        //q[2] should have zero probability to be set
        assert_abs_diff_eq!(register.probability(vec![None, None, Some(true)]), 0.0, epsilon = 0.0001);

        let amps = register.amplitudes().to_owned().into_dimensionality::<Ix4>().unwrap();

        // bit patterns 010 and 110 should have negative imaginary amplitudes of 0.5
        // These are the situations where q[2] was 1 before reset
        assert_abs_diff_eq!(amps[[0,1,0,1]].im, -0.5, epsilon = 0.0001);
        assert_abs_diff_eq!(amps[[1,1,0,1]].im, -0.5, epsilon = 0.0001);

        // bit patterns 000 and 100 should have positive imaginary amplitudes of 0.5
        // These are the situations where q[2] was 0 before reset
        assert_abs_diff_eq!(amps[[1,0,0,0]].im, 0.5, epsilon = 0.0001);
        assert_abs_diff_eq!(amps[[0,0,0,0]].im, 0.5, epsilon = 0.0001);

        //Each bit pattern 010, 110, 100, 000 in the qubits should have 25% prob
        assert_abs_diff_eq!(register.probability(vec![Some(false), Some(true), Some(false)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(true), Some(true), Some(false)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(true), Some(false), Some(false)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(false), Some(false), Some(false)]), 0.25, epsilon = 0.0001);

        // Probabilities should add to 1.0
        assert_relative_eq!(register.probability(vec![None, None, None]), 1.0, epsilon = 0.0001);
    }

    #[test]
    fn condition_measure() {
        let source = r#"
        OPENQASM 2.0;
        include "qelib1.inc";
        
        qreg q[2];
        creg c[1];
        h q[0];
        h q[1];
        measure q[1] -> c[0];
        if (c==0) z q[0];
        measure q[0] -> c[0];
        "#;

        let program = openqasm_parse_source(source.to_string(), Some(Path::new("./")));
        let register = openqasm_run_program(program);

        println!("{:?}",register.amplitudes());

        // This should have resulted in one hidden bit
        assert_eq!(register.bit_count() - register.visible_bit_count(), 1);

        // These are the amplitudes we should have
        let amps = register.amplitudes().to_owned().into_dimensionality::<Ix4>().unwrap();
        assert_abs_diff_eq!(amps[[0,0,0,0]].im,  0.5, epsilon = 0.0001);
        assert_abs_diff_eq!(amps[[0,1,0,1]].re, -0.5, epsilon = 0.0001);
        assert_abs_diff_eq!(amps[[1,0,1,0]].im, -0.5, epsilon = 0.0001);
        assert_abs_diff_eq!(amps[[1,1,1,1]].re, -0.5, epsilon = 0.0001);

        // each two bit pattern should have 25% probability
        assert_abs_diff_eq!(register.probability(vec![Some(false), Some(true)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(true), Some(true)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(true), Some(false)]), 0.25, epsilon = 0.0001);
        assert_abs_diff_eq!(register.probability(vec![Some(false), Some(false)]), 0.25, epsilon = 0.0001);

        // Probabilities should add to 1.0
        assert_relative_eq!(register.probability(vec![None, None, None]), 1.0, epsilon = 0.0001);
    }
}