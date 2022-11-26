use openqasm as oq;

use oq::{
    ast::Symbol,
    translate::{GateWriter, Value},
};

pub struct GatePrinter;

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