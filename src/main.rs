
extern crate openqasm;
extern crate approx;
extern crate nalgebra as na;
mod register;
mod eval;

use openqasm as oq;
use oq::GenericError;
use oq::ProgramVisitor;


fn main() {
    println!("Hello, world!");
    openqasm_parse();
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

    register::RegFinder.visit_program(&program).unwrap();


    let mut l = oq::translate::Linearize::new(eval::GatePrinter,100);
    l.visit_program(&program).to_errors().unwrap();

}

/// Print any sort of errors (ugly but functional if were using openqasm without ariadne)
fn _openqasm_print_errors(errors: oq::Errors) {
    for error in errors.errors {
        println!("{:?}", error);
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



