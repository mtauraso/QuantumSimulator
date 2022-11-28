use num::{complex::Complex32, Zero};
use openqasm as oq;
use oq::{
    translate::{GateWriter, Value}, 
    ast::Symbol};
use nd::{ArrayD,IxDyn, Dimension};




// This is actually all the registers
// First round implementation is that we get a multidimensional 2x2x2x... tensor from ndarray of complex amplitudes
// Then gate ops are all contractions on the tensor
// See here: https://www.kattemolle.com/other/QCinPY.html



#[derive(Debug)]
pub struct QuantumRegister {
    amplitudes: Box<ArrayD<Complex32>>,
    qubit_count: usize,
    qubit_names: Vec<Symbol>,
    cbit_names: Vec<Symbol>,
}

impl QuantumRegister {
    pub fn new() -> QuantumRegister {
        // Actual initialization occurs in gatewriter initialize 
        QuantumRegister { 
            amplitudes: Box::new(ArrayD::zeros(IxDyn(&[2]))),
            qubit_count: 0, 
            qubit_names: Vec::new(), 
            cbit_names: Vec::new(), 
        }
    }


    // Called from gatewriter when the number of bits is available
    fn initialize(&mut self, qubit_names: Vec<Symbol>, cbit_names: Vec<Symbol>) {
        self.qubit_count = qubit_names.len();
        self.qubit_names = qubit_names;
        self.cbit_names = cbit_names;

        let tensor_shape = vec![2; self.qubit_count];
        self.amplitudes = Box::new(ArrayD::zeros(tensor_shape));

        // Set the |0000...0> state only
        let zero_state = vec![0; self.qubit_count];
        self.amplitudes[zero_state.as_slice()] = Complex32::new(1.0,0.0);
    }

    // TODO: Implement gates
    // Apply a cx gate to the numbered qubits
    fn apply_cx(&mut self, _copy:usize, _xor:usize) {
    }

    // TODO: Implement gates
    // Apply a U gate with given parameters to the given qubit
    fn apply_u(&mut self, _theta:Value, _phi:Value, _lambda:Value, _bit: usize) {
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


#[derive(Debug)]
pub struct QuantumRegisterGateWriter<'a> {
    data: &'a mut QuantumRegister
}

impl QuantumRegisterGateWriter <'_> {
    pub fn new(data: &mut QuantumRegister) -> QuantumRegisterGateWriter {
        QuantumRegisterGateWriter { data }
    }
}

impl GateWriter for QuantumRegisterGateWriter <'_> {

    // We return results, but only an error type that can never exist
    // This is essentially a compile time proof that errors are unreachable
    type Error = std::convert::Infallible;

    fn initialize(&mut self, qubits: &[Symbol], cbits: &[Symbol]) -> Result<(), Self::Error> {
        self.data.initialize(qubits.to_vec(), cbits.to_vec());
        Ok(())
    }

    fn write_cx(&mut self, copy: usize, xor: usize) -> Result<(), Self::Error> {
        self.data.apply_cx(copy, xor);
        println!("cx {copy} {xor}");
        Ok(())
    }

    fn write_u(&mut self, theta: Value, phi: Value, lambda: Value, reg: usize) -> Result<(), Self::Error> {
        self.data.apply_u(theta, phi, lambda, reg);
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