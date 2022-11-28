extern crate openqasm as oq;
extern crate ndarray as nd;
extern crate num;

mod register;
mod eval;

use std::f64::consts::FRAC_PI_2;

use num::complex::Complex64;
use oq::GenericError;
use oq::ProgramVisitor;
use nd::Array;

use register::QuantumRegister;


fn main() {
    println!("Hello, world!");
    ndarray_test();
    openqasm_parse();
}

fn ndarray_test() {
    let dim = 12;
    let a = Complex64::from_polar( 1.0, FRAC_PI_2);
    let b = Complex64::from_polar( 1.0, FRAC_PI_2);
    let aa = Array::from_elem( (dim, dim), a);
    let bb = Array::from_elem( (dim, dim), b);

    let cc = aa*bb;

    println!("{:?}", cc);
}

fn openqasm_parse() {

    println!("");

    let mut cache = oq::SourceCache::new();
    let mut parser = oq::Parser::new(&mut cache);
    parser.parse_file("test.qasm");


    let program = match parser.done().to_errors() {
        Ok(program) => program,
        Err(errors) => {
            errors.print(&mut cache).unwrap();
            //openqasm_print_errors(errors);
            panic!("Cannot continue due to Parse errors above.");
        }
    };

    if let Err(errors) = program.type_check().to_errors() {
        errors.print(&mut cache).unwrap();
        //openqasm_print_errors(errors);
        panic!("Cannot continue due to type errors above.");
    }

    let register = QuantumRegister::from_program(&program);
    println!("{}",register.to_string());
    //register.foo();


    let mut l = oq::translate::Linearize::new(eval::GatePrinter,100);
    l.visit_program(&program).to_errors().unwrap();

}

/// Print any sort of errors (ugly but functional if were using openqasm without ariadne)
fn _openqasm_print_errors(errors: oq::Errors) {
    for error in errors.errors {
        println!("{:?}", error);
    }
}
