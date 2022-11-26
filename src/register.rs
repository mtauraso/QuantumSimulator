
use openqasm as oq;
use oq::ProgramVisitor;

pub struct RegFinder;

impl ProgramVisitor for RegFinder {

    type Error = std::convert::Infallible;

    fn visit_qreg(&mut self, reg: &oq::Span<oq::Reg>) -> Result<(), Self::Error> {
        println!("Qreg: {} [{}]", reg.inner.name.as_str(), reg.inner.index.unwrap());
        Ok(())
    }
    fn visit_creg(&mut self, reg: &oq::Span<oq::Reg>) -> Result<(), Self::Error> {
        println!("Creg: {} [{}]", reg.inner.name.as_str(), reg.inner.index.unwrap());
        Ok(())
    }
}
