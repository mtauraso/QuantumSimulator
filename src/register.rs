use std::{collections::HashMap};
use num::{complex::Complex32, ToPrimitive, Zero};
use openqasm as oq;
use oq::{ProgramVisitor, Program};
use nd::{ArrayD,IxDyn, Dimension};


// TODO:
// RegFinder should collect the information, what I'm going to need:
// 1) A count of total qubits
// 2) enough information to map a q register name/index to an row/col in my 2^n dimension monster array
//
// Then something needs to process the register spec and actually allocate the giant array
// Then we can do something with it.



// This is actually all the registers
// First round implementation is that we get a multidimensional 2x2x2x... tensor from ndarray of complex amplitudes
// Then gate ops are all contractions on the tensor
// See here: https://www.kattemolle.com/other/QCinPY.html
#[derive(Debug)]
pub struct QuantumRegister {
    amplitudes: Box<ArrayD<Complex32>>,
    qubit_count: u8,
    qreg_spec: HashMap<String, LogicalQuantumRegister>,
}

/// Stores a mapping between the indexes of the quantum registers in the quantum program 
/// language and the master quantum register for the entire simulator
#[derive(Debug)]
struct LogicalQuantumRegister {
    offset: u8,
    size: u8,
}

impl QuantumRegister {

    /// New quantum register of the required size
    /// Automatically initializes to the state |000...0>
    /// Warning: Allocates 2^dim complex32 values. I think on the stack since there is no box.
    pub fn from_program(program: &Program) -> QuantumRegister {

        let mut qreg = QuantumRegister { 
            // This is technically an invalid state, but we need to initialize amplitudes 
            // to something with the correct type to move on.
            amplitudes: Box::new(ArrayD::zeros(IxDyn(&[2]))),
            qreg_spec: HashMap::new(),
            qubit_count: 0,
        };

        // Visit the passed in program in order to know how many quantum registers we need
        qreg.visit_program(program).unwrap();

        // Create the tensor to hold amplitudes
        let dim = qreg.qubit_count.to_usize().unwrap();
        let tensor_shape = vec![2; dim];
        qreg.amplitudes = Box::new(ArrayD::zeros(tensor_shape));

        // Set the |0000...0> state only
        let zero_state = vec![0;dim];
        qreg.amplitudes[zero_state.as_slice()] = Complex32::new(1.0,0.0);

        qreg
    }

    pub fn to_string(& self) -> String{
        let mut out = String::new();

        for (indicies,amplitude) in self.amplitudes.indexed_iter() {
            if ! Complex32::eq(amplitude, &Complex32::zero()) {

                if out.len() > 0 {
                    out.push_str(" + ")
                }

                if amplitude.im.is_zero() {
                    out.push_str(format!("{} |", amplitude.re).as_str());
                } else if amplitude.re.is_zero() {
                    out.push_str(format!("{} i |", amplitude.im).as_str());
                } else {
                    out.push_str(format!("({} + {} i)|", amplitude.re, amplitude.im).as_str());
                }

                for bit in indicies.as_array_view() {
                    out.push_str(format!("{}", bit).as_str());
                }
                out.push_str(">");
                //println!("{:?} {:?}", indicies.as_array_view(), amplitude);
            }
        }
        out
    }
}

impl ProgramVisitor for QuantumRegister {

    type Error = std::convert::Infallible;

    fn visit_qreg(&mut self, reg: &oq::Span<oq::Reg>) -> Result<(), Self::Error> {
        let register_qubits = reg.inner.index.unwrap();

        self.qreg_spec.insert(
            reg.inner.name.to_string(), 
            LogicalQuantumRegister { 
                offset: self.qubit_count, 
                size: register_qubits.to_u8().unwrap() });

        self.qubit_count += register_qubits.to_u8().unwrap();

        println!("Qreg: {} [{}]", reg.inner.name.as_str(), reg.inner.index.unwrap());
        Ok(())
    }
    fn visit_creg(&mut self, reg: &oq::Span<oq::Reg>) -> Result<(), Self::Error> {
        println!("Creg: {} [{}]", reg.inner.name.as_str(), reg.inner.index.unwrap());
        Ok(())
    }
}