use num::{
    complex::Complex32, 
    complex::Complex,
    Zero, 
    Rational64,
    traits::{CheckedAdd, CheckedSub,CheckedDiv}, Num
};
use openqasm as oq;
use oq::{
    translate::{GateWriter, Value}, 
    ast::Symbol};
// use nd::{ArrayD,IxDyn, Dimension};
use nd::prelude::*;
use ndarray_einsum_beta::*;


// This is actually all the registers
// First round implementation is that we get a multidimensional 2x2x2x... tensor from ndarray of complex amplitudes
// Then gate ops are all contractions on the tensor
// See here: https://www.kattemolle.com/other/QCinPY.html
// I have extended this method a bit to deal with cbits
// The layout of tensor indicies for N qubits and M cbits is: [qubit0, ..., qubitN, cbit0, ..., cbitM]
// All tensor contraction operations put the tensor indicies back in this state when they are done.
#[derive(Debug)]
pub struct QuantumRegister {
    amplitudes: Box<ArrayD<Complex32>>,
    bit_count: usize,
    qubit_count: usize,
    cbit_count: usize,
    qubit_names: Vec<Symbol>,
    cbit_names: Vec<Symbol>,
}

impl QuantumRegister {
    pub fn new() -> QuantumRegister {
        // Actual initialization occurs in gatewriter initialize 
        QuantumRegister { 
            amplitudes: Box::new(ArrayD::zeros(IxDyn(&[2]))),
            bit_count: 0,
            qubit_count: 0, 
            cbit_count: 0,
            qubit_names: Vec::new(), 
            cbit_names: Vec::new(), 
        }
    }


    // Called from gatewriter when the number of bits is available
    fn initialize(&mut self, qubit_names: Vec<Symbol>, cbit_names: Vec<Symbol>) {
        self.qubit_count = qubit_names.len();
        self.cbit_count = cbit_names.len();
        self.qubit_names = qubit_names;
        self.cbit_names = cbit_names;
        self.bit_count = self.qubit_count + self.cbit_count;

        let tensor_shape = vec![2; self.bit_count];
        self.amplitudes = Box::new(ArrayD::zeros(tensor_shape));

        // Set the |0000...0> state only
        // Also set all classical bits to zero
        let zero_state = vec![0; self.bit_count];
        self.amplitudes[zero_state.as_slice()] = Complex32::new(1.0,0.0);
        println!("{}", self.to_string());
    }

    // Apply a cx gate to the numbered qubits
    fn apply_cx(&mut self, control_qubit:usize, target_qubit:usize) {

        // TODO: Move this definition elsewhere so we can use it as a const
        // Cnot matrix from OpenQASM spec: "Open Quantum Assembly Language" Cross, Bishop, Smolin Gambetta (2017)
        let c_1 = Complex32::new(1.0, 0.0);
        let c_0 = Complex32::new(0.0, 0.0);
        let cnot_matrix = array![
            [c_1,c_0,c_0,c_0],
            [c_0,c_1,c_0,c_0],
            [c_0,c_0,c_0,c_1],
            [c_0,c_0,c_1,c_0]];

        //println!("{}", cnot_matrix);

        let cnot_tensor = cnot_matrix.into_shape([2,2,2,2]).unwrap();
        //println!("{}", cnot_tensor);

        // Tensordot does the first two
        // Create a mutable view with the axes of amplitude permuted so the two active bits are at the end
        //  e.g index order is Amplitudes_a,b,....,c,t
        // Contract this view with CNOT_nc,nt,c,t
        // Result should have form Amplitudes_a,b,...,nc,nt

        // Is the result a new array..?
        //   If so we need the stored array to have the right order. 
        //   We could store a view into the new array with the correct permutation. 
        //   How well this works will depend on the lifetime rules, and the memory is going to be views-on-views for multiple ops
        //
        //   Or we could copy all the amplitudes (probably bad for perf)


        // Construct the permutation of axes we'll need to un-permute the result of the
        // tensor dot product below
        let mut axis_permutation = Vec::with_capacity(self.bit_count);
        let mut offset = 2;
        for axis in num::range(0,self.bit_count){
            axis_permutation.push(
                if axis == control_qubit { offset -= 1; 0} 
                else if axis == target_qubit { offset -= 1; 1}
                else {axis + offset }
            );
        }

        // Calculate new amplitudes as a tensor dot product which will compute
        // CNOT_nc,nt,c,t Amplitude_x,y,...,c,...,t,...,z -> NewAmplitude_nc,nt,x,y,...,z
        // permute axes returns a,b to the old locations of c, t so we get
        // NewAmplitude_x,y,...,nc,...,nt,...,z
        // Note, this effectively copies all the amplitudes, what we get back from tensordot is owned
        let amplitude_result = 
        tensordot(
            &cnot_tensor, 
            &self.amplitudes, 
            &[Axis(2),Axis(3)],
            &[Axis(control_qubit), Axis(target_qubit)]
        ).permuted_axes(IxDyn(axis_permutation.as_slice()));

        // Replace the amplitudes with the new ones.
        // TODO TRACE: this is where we would need to keep history
        self.amplitudes = Box::new(amplitude_result);

        println!("{}", self.to_string());

    }

    // TODO: Implement gates
    // Apply a U gate with given parameters to the given qubit
    fn apply_u(&mut self, theta:Value, phi:Value, lambda:Value, target_qubit: usize) {
        // Todo: These can probably be computed at program parse time and then used during execution
        // U matrix from OpenQASM spec: "Open Quantum Assembly Language" Cross, Bishop, Smolin Gambetta (2017)
        let philambdaplus = 
            checked_div(
                checked_add(phi, lambda).unwrap(),
                Rational64::from(2)).unwrap().into_float();
        let philambdaminus = 
            checked_div(
                checked_sub(phi, lambda).unwrap(),
                Rational64::from(2)).unwrap().into_float();
        let halftheta = checked_div(theta,Rational64::from(2)).unwrap();
        let sinhalftheta = halftheta.checked_sin().unwrap().into_float();
        let coshalftheta = halftheta.checked_cos().unwrap().into_float();

        let term00 = Complex32::from_polar(coshalftheta, -philambdaplus);
        let term01 = Complex32::from_polar(-sinhalftheta, -philambdaminus);
        let term10 = Complex32::from_polar(sinhalftheta, philambdaminus);
        let term11 = Complex32::from_polar(coshalftheta, philambdaplus);

        let u_matrix = array![
            [term00,term01],
            [term10,term11],
        ];

        // Construct the permutation of axes we'll need to un-permute the result of the
        // tensor dot product below
        let mut axis_permutation = Vec::with_capacity(self.bit_count);
        let mut offset = 1;
        for axis in num::range(0,self.bit_count){
            axis_permutation.push(
                if axis == target_qubit { offset -= 1; 0}
                else {axis + offset }
            );
        }

        // Compute the new amplitudes by doing a tensor contraction on the existing axes
        let amplitude_result =
        tensordot(
            &u_matrix, 
            &self.amplitudes, 
            &[Axis(1)],
            &[Axis(target_qubit)]
        ).permuted_axes(IxDyn(axis_permutation.as_slice()));

        self.amplitudes = Box::new(amplitude_result);

        println!("{}", self.to_string());

    }


    // The simple way to keep this flexible is to treat the cbits just like they are qubits
    // Lets say you have a = cbit and b,c = qubits amplitudes are A_abc
    // Lets say you entangle b and c so your whole amplitude circus is in A_0bc because a is set to zero
    // When you measure qubit b into cbit a A_1bc becomes the projection of b/c state onto b=1 b is in |1> eignenstate
    // A_0bc becomes the projection of b/c state into b=0 and b is in |0> eigenstate
    // You can still hold everything in A_abc
    // Quantum gates work as before because you can still contract the amplitudes in b/c no matter that you have multiple copies
    // depending on what a is. The summation in the contraction still sums over those.
    // Classical "if" gates work by only manipulating a slice where the cbits have the desired values
    fn measure(&mut self, from_qubit:usize, to_cbit:usize) {

        // Construct the index of the classical bit
        let cbit_index = to_cbit + self.qubit_count;

        // Construct projection tensors for zero and one values.
        let c_1 = Complex32::new(1.0, 0.0);
        let c_0 = Complex32::new(0.0, 0.0);
        let zero_projector = array![
            [c_1,c_0],
            [c_0,c_0],
        ];
        let one_projector = array![
            [c_0,c_0],
            [c_0,c_1],
        ];

        // Construct the permutation of axes we'll need to un-permute the result of the
        // tensor dot products below. This unpermuting will help keep the tensor indexes
        // that represent our qubits and cbits in a consistent order so we can address them.
        let mut axis_permutation = Vec::with_capacity(self.bit_count);
        let mut offset = 1;
        for axis in num::range(0,self.bit_count){
            axis_permutation.push(
                if axis == from_qubit { offset -= 1; 0}
                else {axis + offset }
            );
        }

        // Compute the new amplitudes by doing a tensor contraction that projects the state
        // into the world where the qubit is measured to be zero
        let amplitude_result_zero=
        tensordot(
            &zero_projector, 
            &self.amplitudes, 
            &[Axis(1)],
            &[Axis(from_qubit)])
        .permuted_axes(IxDyn(axis_permutation.as_slice()))
        // Then sum over the possible values of the cbit before we overwrite it
        .sum_axis(Axis(cbit_index));

        // Compute the new amplitudes by doing a tensor contraction that projects the state
        // into the world where the qubit is measured to be one.
        let amplitude_result_one =
        tensordot(
            &one_projector, 
            &self.amplitudes, 
            &[Axis(1)],
            &[Axis(from_qubit)]) 
        .permuted_axes(IxDyn(axis_permutation.as_slice()))
        // Then sum over the possible values of the cbit before we overwrite it.
        .sum_axis(Axis(cbit_index));

        // Some Debug information .
        // println!("Amplitudes when the qubit is measured as zero");
        // println!("{}",Self::to_string_internal(& amplitude_result_zero.view()));
        // println!("{:?}",amplitude_result_zero);
        // println!("{:.3}", norm_sqr(&amplitude_result_zero));
        // println!("Amplitudes when the qubit is measured as one");
        // println!("{}",Self::to_string_internal(& amplitude_result_one.view()));
        // println!("{:?}",amplitude_result_one);
        // println!("{:.3}", norm_sqr(&amplitude_result_one));


        // Set the amplitude tensor to the projections we got. This effectively 
        // 1) "overwrites" whatever value the cbit had
        // 2) Collapses the state of the qubit to the basis state it was measured in
        //
        // When the cbit is zero, assign the amplitudes from the zero projection.
        self.amplitudes.index_axis_mut(Axis(cbit_index),0).assign(&amplitude_result_zero);
        // When the cbit is one, assign the amplitudes from the one projection.
        self.amplitudes.index_axis_mut(Axis(cbit_index),1).assign(&amplitude_result_one);


        println!("{}", self.to_string());
    }

    pub fn to_string(& self) -> String{
        Self::to_string_internal(&(*self.amplitudes).view())
    }

    fn to_string_internal(amplitudes:& ArrayViewD<Complex32>) -> String{
        let mut out = String::new();

        for (indicies,amplitude) in amplitudes.indexed_iter() {
            if ! Complex32::eq(amplitude, &Complex32::zero()) {

                if out.len() > 0 {
                    out.push_str(" + ")
                }

                // TODO:
                // Use approx for all zero comparisons with some accuracy level
                // Put an accuracy level in for rounding
                // Add a case for the whole number being zero and just printing zero
                if amplitude.im.is_zero() {
                    out.push_str(format!("{:.3} |", amplitude.re).as_str());
                } else if amplitude.re.is_zero() {
                    out.push_str(format!("{:.3} i |", amplitude.im).as_str());
                } else if amplitude.im.is_sign_negative(){
                    out.push_str(format!("({:.3} - {:.3} i)|", amplitude.re, -amplitude.im).as_str());
                } else {
                    out.push_str(format!("({:.3} + {:.3} i)|", amplitude.re, amplitude.im).as_str());
                }

                for bit in indicies.as_array_view() {
                    out.push_str(format!("{:.3}", bit).as_str());
                }
                out.push_str(">");
                //println!("{:?} {:?}", indicies.as_array_view(), amplitude);
            }
        }

        if out.len() == 0 {
            out.push_str("0");
        }

        out
    }
}


fn _norm_sqr<T: Num + Clone> (amplitudes:&ArrayD<Complex<T>>) -> T {
    amplitudes.fold(T::zero(), |accum:T, elem:&Complex<T>| -> T{
        accum + elem.norm_sqr()
    })
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
        println!("cx control qubit {copy} target qubit {xor}");
        self.data.apply_cx(copy, xor);
        Ok(())
    }

    fn write_u(&mut self, theta: Value, phi: Value, lambda: Value, reg: usize) -> Result<(), Self::Error> {
        println!("u({theta}, {phi}, {lambda}) on qubit {reg}");
        self.data.apply_u(theta, phi, lambda, reg);
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
        println!("measure qubit {} -> cbit {}", from, to);
        self.data.measure(from, to);
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


// Looks like OpenQasm thought math for oq::Value was carried over from num::Rational64, 
// but the crate doesn't implement these, so we have just enough here to muddle through
fn checked_add(lhs: Value, rhs: Value) -> Option<Value> {
    Some(Value {
        a: lhs.a.checked_add(&rhs.a)?,
        b: lhs.b.checked_add(&rhs.b)?,
    })
}

fn checked_sub(lhs: Value, rhs: Value) -> Option<Value> {
    Some(Value {
        a: lhs.a.checked_sub(&rhs.a)?,
        b: lhs.b.checked_sub(&rhs.b)?,
    })
}

fn checked_div(lhs: Value, rhs: Rational64) -> Option<Value>{
    Some(Value {
        a: lhs.a.checked_div(&rhs)?,
        b: lhs.b.checked_div(&rhs)?,
    })
}