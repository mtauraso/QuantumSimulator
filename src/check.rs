use simple_error::SimpleError;
use openqasm as oq;

use oq::{
    translate::{GateWriter, Value}, 
    ast::Symbol
};

pub struct CheckGateWriter;
impl GateWriter for CheckGateWriter {
    
    type Error = SimpleError;

    fn initialize(&mut self, _: &[Symbol], _: &[Symbol]) -> Result<(), Self::Error> { Ok(()) }
    fn write_cx(&mut self, _: usize, _xor: usize) -> Result<(), Self::Error> { Ok(()) }
    fn write_u(&mut self, _: Value, _: Value, _: Value, _: usize) -> Result<(), Self::Error> { Ok(()) }
    fn write_barrier(&mut self, _: &[usize]) -> Result<(), Self::Error> { Ok(()) }
    fn write_measure(&mut self, _: usize, _: usize) -> Result<(), Self::Error> { Ok(()) }
    fn write_reset(&mut self, _: usize) -> Result<(), Self::Error> { Ok(()) }
    fn start_conditional(&mut self, _: usize, _: usize, _: u64) -> Result<(), Self::Error> { Ok(()) }
    fn end_conditional(&mut self) -> Result<(), Self::Error> { Ok(()) }

    fn write_opaque(&mut self, _: &Symbol, _: &[Value], _: &[usize]) -> Result<(), Self::Error> {
        Err(Self::Error::new("Opaque gates Not Implemented"))
    }
}