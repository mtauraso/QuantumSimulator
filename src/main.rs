
extern crate openqasm;
extern crate approx;
extern crate nalgebra as na;

use std::env;
use std::fs::File;
use std::io::prelude::*;

use openqasm as oq;
use oq::GenericError;
use oq::ProgramVisitor;


fn main() {
    println!("Hello, world!");
    openqasm_parse();
}

fn openqasm_parse() {
    let _cwd = env::current_dir().unwrap();
    let mut source = String::new();

    let mut f = File::open("test.qasm").expect("cannot find source file 'test.qasm'");
    f.read_to_string(&mut source).expect("couldn't read file 'test.qasm'");

    println!("");
    //println!("{:?}", source);

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

    RegFinder.visit_program(&program).unwrap();


    let mut l = oq::translate::Linearize::new(GatePrinter,100);
    l.visit_program(&program).to_errors().unwrap();

}

/// Print any sort of errors (ugly but functional if were using openqasm without ariadne)
fn _openqasm_print_errors(errors: oq::Errors) {
    for error in errors.errors {
        println!("{:?}", error);
    }
}

use oq::{
    ast::Symbol,
    translate::{GateWriter, Value},
};

struct GatePrinter;

impl GateWriter for GatePrinter {

    // We return results, but only an error type that can never exist
    // This is essentially a compile time proof that errors are unreachable
    type Error = std::convert::Infallible;

    fn initialize(&mut self, _: &[Symbol], _: &[Symbol]) -> Result<(), Self::Error> {
        Ok(())
    }

    fn write_cx(&mut self, copy: usize, xor: usize) -> Result<(), Self::Error> {
        println!("cx {copy} {xor}");
        Ok(())
    }

    fn write_u(&mut self, theta: Value, phi: Value, lambda: Value, reg: usize) -> Result<(), Self::Error> {
        println!("u({theta}, {phi}, {lambda}) {reg}");
        Ok(())
    }

    fn write_opaque(&mut self, name: &Symbol, _: &[Value], _: &[usize]) -> Result<(), Self::Error> {
        println!("opaque gate {}", name);
        Ok(())
    }

    fn write_barrier(&mut self, _: &[usize]) -> Result<(), Self::Error> {
        Ok(())
    }

    fn write_measure(&mut self, from: usize, to: usize) -> Result<(), Self::Error> {
        println!("measure {} -> {}", from, to);
        Ok(())
    }

    fn write_reset(&mut self, reg: usize) -> Result<(), Self::Error> {
        println!("reset {reg}");
        Ok(())
    }

    fn start_conditional(&mut self, reg: usize, count: usize, value: u64) -> Result<(), Self::Error> {
        println!("if ({reg}:{count} == {value}) {{");
        Ok(())
    }

    fn end_conditional(&mut self) -> Result<(), Self::Error> {
        println!("}}");
        Ok(())
    }
}

struct RegFinder;

impl ProgramVisitor for RegFinder {

    type Error = std::convert::Infallible;

    fn visit_qreg(&mut self, reg: &oq::Span<oq::Reg>) -> Result<(), Self::Error> {
        println!("Qreg: {} [{}]", reg.inner.name.as_str(), reg.inner.index.unwrap());
        Ok(())
    }
    fn visit_creg(&mut self, reg: &oq::Span<oq::Reg>) -> Result<(), Self::Error> {
        println!("Qreg: {} [{}]", reg.inner.name.as_str(), reg.inner.index.unwrap());
        Ok(())
    }
}


// use approx::{relative_eq};
// use na::{Vector3, Rotation3};
// fn nalgebra_test() {

// // TODO: nalgebra might be overkill

//     let axis  = Vector3::x_axis();
//     let angle = 1.57;
//     let b     = Rotation3::from_axis_angle(&axis, angle);

//     relative_eq!(b.axis().unwrap(), axis);
//     relative_eq!(b.angle(), angle);

//     println!("");
//     println!("{:?}", b);
// }



