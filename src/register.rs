use ndarray as nd;
use num::{
    complex::Complex32, 
    complex::Complex,
    Zero, One,
    Rational64,
    traits::{CheckedAdd,CheckedDiv}, Num
};
use openqasm as oq;
use oq::{
    translate::{GateWriter, Value}, 
    ast::Symbol};
use nd::prelude::*;
use ndarray_einsum_beta::*;
use approx::*;
use std::collections::HashMap;

use crate::cli::Cli;
use clap::Parser;


// This is actually all the registers
// First round implementation is that we get a multidimensional 2x2x2x... tensor from ndarray of complex amplitudes
// Then gate ops are all contractions on the tensor
// See here: https://www.kattemolle.com/other/QCinPY.html
// I have extended this method a bit to deal with cbits, but it is the same basic idea.
//
// The layout of tensor indicies for N qubits and M cbits is: [qubit0, ..., qubitN, cbit0, ..., cbitM]
// All tensor contraction operations put the tensor indicies back in this state when they are done.
#[derive(Debug)]
pub struct QuantumRegister {
    amplitudes: Box<ArrayD<Complex32>>,

    // used for counting indicies
    qubit_count: usize,
    cbit_count: usize,
    hidden_qubit_count: usize,

    // Used for tracking cbit state
    cbit_values: Vec<Option<bool>>,

    // Currently unused
    qubit_names: Vec<Symbol>,
    cbit_names: Vec<Symbol>,

    // Command line options that affect output
    trace: bool,
    debug: bool,

    //Start index and extent for each cregister/qregister
    // Start indicies are relative to the start of q/c registers
    // in the amplitude indexes
    cregister_extents: HashMap<String, (usize, usize)>,
    qregister_extents: HashMap<String, (usize, usize)>,

    // Used for conditionals
    amplitudes_save: Option<Box<ArrayD<Complex32>>>,
    pinned_values: Option<Vec<Option<bool>>>,
}

impl QuantumRegister {
    pub fn new() -> QuantumRegister {
        // Actual initialization occurs in gatewriter initialize 
        QuantumRegister { 
            amplitudes: Box::new(ArrayD::zeros(IxDyn(&[2]))),
            cbit_values: Vec::new(),
            amplitudes_save: None,
            pinned_values: None,
            qubit_count: 0, 
            cbit_count: 0,
            hidden_qubit_count: 0,
            qubit_names: Vec::new(), 
            cbit_names: Vec::new(),
            trace: false,
            debug: false,
            cregister_extents: HashMap::new(),
            qregister_extents: HashMap::new(),
        }
    }

    #[cfg(test)]
    pub fn amplitudes(&self) -> ArrayViewD<Complex32> { self.amplitudes.view() }
    #[cfg(test)]
    pub fn visible_bit_count(&self) -> usize{ self.qubit_count + self.cbit_count }

    pub fn bit_count(&self) -> usize { self.qubit_count + self.cbit_count + self.hidden_qubit_count }


    pub fn print_probabilities(&self) {
        // Ignores hidden qubits, simply tracing over them
        // Prints out a summary at point-in-time by classical register for all possible bit patterns

        let max_cbit_value:usize = 1 << self.cbit_count;

        for cbit_value in 0..max_cbit_value {
            let pinned_bit_pattern = self.cbit_pattern_to_pinned_bits(cbit_value, 0, self.cbit_count);
            let probability = self.probability(&pinned_bit_pattern);

            if abs_diff_eq!(probability,  0.0, epsilon = 0.0001) {
                continue;
            } else {
                let register_string = self.bit_pattern_as_register_names(&pinned_bit_pattern);
                println!("Probability for {}: {:.4}", register_string, probability);
            }
        }
    }

    fn cbit_pattern_to_pinned_bits(&self, creg_value:usize, creg_offset:usize, creg_size:usize) -> Vec<Option<bool>> {
         // Make a list of which indicies are pinned and to what values
         let cbit_index = self.qubit_count + creg_offset;
         let mut pinned_values:Vec<Option<bool>> = vec![None; cbit_index];
         for bit_index in 0..creg_size {
             match (creg_value >> bit_index) & 1 {
                 0 => pinned_values.push(Some(false)),
                 1 => pinned_values.push(Some(true)),
                 _ => panic!("bitwise and of a number with 1 didn't give 1 or zero."),
             }
         }
         pinned_values
    }

    // Given a bit pattern (vec of option bools) for values of the cbits
    // Look at the register and and me back a string describing all the classical registers
    // form should be reg = <binary val>, reg = <binary val>
    pub fn bit_pattern_as_register_names (&self, pinned_bits: &Vec<Option<bool>>) -> String {

        let qreg_string = Self::register_printable_helper(pinned_bits, 0, &self.qregister_extents);
        let creg_string = Self::register_printable_helper(pinned_bits, self.qubit_count, &self.cregister_extents);

        qreg_string + creg_string.as_str()
    }

    fn register_printable_helper(pinned_bits: &Vec<Option<bool>>, reg_extents_start:usize, register_extents: &HashMap<String, (usize, usize)> ) -> String {
        let mut returnme:Vec<String> = Vec::new();

        for (key, (offset, size)) in register_extents.iter() {
            let mut all_unset = true;

            let first_index:usize = reg_extents_start + offset;
            let last_index:usize = first_index + size;
            let bits = pinned_bits[first_index..last_index].to_vec();

            // Note we iterate in reverse order here, so the most significant bit with the largest index
            // is printed out first.
            let bit_pattern: Vec<&str> = bits.iter().rev().map(|bit_option| {
                match bit_option {
                    Some(true) => {all_unset = false; "1"},
                    Some(false) => {all_unset = false; "0"},
                    None => "_",
                }
            }).collect();

            if all_unset == false {
                returnme.push(format!("{}={}", key.as_str(), bit_pattern.join("").as_str()));
            }
        }
        returnme.join(",")
    }

    // Called from gatewriter when the number of bits is available
    fn initialize(&mut self, qubit_names: Vec<Symbol>, cbit_names: Vec<Symbol>) {

        let cli = Cli::parse();

        self.trace = cli.trace;
        self.debug = cli.debug;

        self.qubit_count = qubit_names.len();
        self.cbit_count = cbit_names.len();
        self.hidden_qubit_count = 0;

        let tensor_shape = vec![2; self.bit_count()];
        self.amplitudes = Box::new(ArrayD::zeros(tensor_shape));

        // Set the |0000...0> state only
        // Also set all classical bits to zero
        let zero_state = vec![0; self.bit_count()];
        self.amplitudes[zero_state.as_slice()] = Complex32::new(1.0,0.0);

        self.cbit_values = vec![Some(false); self.cbit_count];


        // Set up cregister extents, and save the bit names passed in by caller
        Self::populate_register_hash(&cbit_names, &mut self.cregister_extents);
        Self::populate_register_hash(&qubit_names, &mut self.qregister_extents);
        self.qubit_names = qubit_names;
        self.cbit_names = cbit_names;

        if self.trace { println!("{}", self.to_string()) };
    }


    fn populate_register_hash ( bit_names: & Vec<Symbol>, register_extents: & mut HashMap<String, (usize, usize)>) {

        let mut register_name = String::new();
        let mut start_index:usize = 0;
        for (cbit_index, cbit_symbol) in bit_names.iter().enumerate() {
            let cbit_name = cbit_symbol.to_string();
            let register_parse: Vec<&str> = cbit_name.split(&['[',']'][..]).collect();

            let current_bit:usize = register_parse[1].to_string().parse::<usize>().unwrap();

            if cbit_index == 0 {
                start_index = cbit_index;
                register_name = register_parse[0].to_string();
                continue;
            }

            match current_bit {
                 0 => {
                    // Insert the last register we saw, because this is the beginning of a new one
                    register_extents.insert(register_name, (start_index, cbit_index - start_index));

                    // Since this is a new register, save our start index and the new name
                    start_index = cbit_index;
                    register_name = register_parse[0].to_string();
                 },
                 _ => {}
            }
        }
        register_extents.insert(register_name, (start_index, bit_names.len() - start_index));
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
        let mut axis_permutation = Vec::with_capacity(self.bit_count());
        let mut offset = 2;
        for axis in num::range(0,self.bit_count()){
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

        if self.trace { println!("{}", self.to_string()) };
    }

    // Apply a U gate with given parameters to the given qubit
    fn apply_u(&mut self, theta:Value, phi:Value, lambda:Value, target_qubit: usize) {
        // Todo: These can probably be computed at program parse time and then used during execution
        //
        // U matrix from OpenQasm 3.0 standard at https://openqasm.com/language/gates.html#built-in-gates
        // I am not using the matrix from the OpenQASM 2.0 spec because it is annoyingly different
        // from the standard set of matricies, especially for 
        // hadamard which is no longer [[1,1][1,-1]] in strict OpenQASM 2.0.
        let philambdaplus = checked_add(phi, lambda).unwrap().into_float();

        let halftheta = checked_div(theta,Rational64::from(2)).unwrap();
        let sinhalftheta = halftheta.checked_sin().unwrap().into_float();
        let coshalftheta = halftheta.checked_cos().unwrap().into_float();

        let term00 = Complex32::from_polar(coshalftheta, 0.0);
        let term01 = Complex32::from_polar(-sinhalftheta, lambda.into_float());
        let term10 = Complex32::from_polar(sinhalftheta, phi.into_float());
        let term11 = Complex32::from_polar(coshalftheta, philambdaplus);

        let u_matrix = array![
            [term00,term01],
            [term10,term11],
        ];

        // Construct the permutation of axes we'll need to un-permute the result of the
        // tensor dot product below
        let mut axis_permutation = Vec::with_capacity(self.bit_count());
        let mut offset = 1;
        for axis in num::range(0,self.bit_count()){
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
            &[Axis(target_qubit)])
        .permuted_axes(IxDyn(axis_permutation.as_slice()));

        
        self.amplitudes = Box::new(amplitude_result);

        if self.trace { println!("{}", self.to_string()) };
    }

    // Handling a measurement:
    // First determine if the target cbit is zero or one
    //      If the cbit is zero we do not need a hidden bit
    //      If the cbit is one we also do not need a hidden bit, but need to swap 0,1 or correct this later...?
    // If we need a hidden bit, we create one at the end, and then swap the index for the used cbit to the end
    //
    // At this point the bit at the cbit's index is zero
    //
    // We then proceed to do our 0 and 1 projection measurements, removing the cbit index appropriately
    // 
    // Then we paste in the projected amplitudes holding the cbit axis equal to whichever projection
    // we are working with
    fn measure(&mut self, from_qubit:usize, to_cbit:usize) {

        // Construct the index of the classical bit
        let cbit_index = to_cbit + self.qubit_count;

        // Construct projection tensors for zero and one values.
        let c_1 = Complex32::one();
        let c_0 = Complex32::zero();
        let zero_projector = array![
            [c_1,c_0],
            [c_0,c_0],
        ];
        let one_projector = array![
            [c_0,c_0],
            [c_0,c_1],
        ];

        // // Construct the permutation of axes we'll need to un-permute the result of the
        // // tensor dot products below. This unpermuting will help keep the tensor indexes
        // // that represent our qubits and cbits in a consistent order so we can address them.
        // let mut axis_permutation = Vec::with_capacity(self.bit_count());
        // let mut offset = 1;
        // for axis in num::range(0,self.bit_count()){
        //     axis_permutation.push(
        //         if axis == from_qubit { offset -= 1; 0}
        //         else {axis + offset }
        //     );
        // }

        let initial_cbit_value:usize = match self.cbit_values[to_cbit] {
            Some(true) => {1},
            Some(false) => {0},
            None => {
                // Make a new hidden bit and swap it into the spot where the cbit is
                self.hide_bit(cbit_index);
                0
            },
        };

        // Compute the new amplitudes by doing a tensor contraction that projects the state
        // into the world where the qubit is measured to be zero
        let mut amplitude_result_zero=
        tensordot(
            &zero_projector, 
            &self.amplitudes, 
            &[Axis(1)],
            &[Axis(from_qubit)]);

        amplitude_result_zero.swap_axes(from_qubit, 0); 
 
        // Compute the new amplitudes by doing a tensor contraction that projects the state
        // into the world where the qubit is measured to be one.
        let mut amplitude_result_one =
        tensordot(
            &one_projector, 
            &self.amplitudes, 
            &[Axis(1)],
            &[Axis(from_qubit)]);
        
        amplitude_result_one.swap_axes(from_qubit, 0); 
    

        // Some Debug information .
        // println!("Amplitudes when the qubit is measured as zero");
        // println!("{}",Self::to_string_internal(& amplitude_result_zero.view()));
        // println!("{:?}",amplitude_result_zero);
        // println!("Amplitudes when the qubit is measured as one");
        // println!("{}",Self::to_string_internal(& amplitude_result_one.view()));
        // println!("{:?}",amplitude_result_one);


        // Set the amplitude tensor to the projections we got. This effectively 
        // 1) "overwrites" whatever value the cbit had
        // 2) Collapses the state of the qubit to the basis state it was measured in
        //
        // When the cbit is zero, assign the amplitudes from the zero projection.
        self.amplitudes.index_axis_mut(Axis(cbit_index),0).assign(&amplitude_result_zero.index_axis_mut(Axis(cbit_index), initial_cbit_value));
        // When the cbit is one, assign the amplitudes from the one projection.
        self.amplitudes.index_axis_mut(Axis(cbit_index),1).assign(&amplitude_result_one.index_axis_mut(Axis(cbit_index), initial_cbit_value));


        // Set hints for the next measurement
        self.cbit_values[to_cbit] = 
            if abs_diff_eq!(self.bit_probability(cbit_index, false), 0.0, epsilon = 0.00001) { Some(true) }
            else if abs_diff_eq!(self.bit_probability(cbit_index, true), 0.0, epsilon = 0.00001) { Some(false) } 
            else { None };


        if self.trace { println!("{}", self.to_string()) };
    }



    // Adds a hidden bit at the end of amplitudes
    // This is used when state is "destroyed" in order to continue 
    // to be able to trace across all possibilities.
    fn hide_bit(&mut self, bit_index: usize) {

        //println!("Hiding bit {}", bit_index);
        //println!("total bits {}", self.bit_count());

        let hidden_bit_index = self.bit_count();

        // Create the new index we will need
        self.amplitudes.insert_axis_inplace(Axis(hidden_bit_index));


        // // Modify the dimensions of the new axis to accomodate 0 and 1 states
        let mut dimension = self.amplitudes.dim();
        let mut dimension_spec = dimension.as_array_view_mut();
        dimension_spec[hidden_bit_index] = 2;

        // Take ownership of amplitudes, because we're going to end up replacing it
        let amplitudes = self.amplitudes.to_owned();

        // Make a copy of our amplitudes by broadcasting the existing one to the larger shape that includes the hidden bit
        let mut new_amplitudes = amplitudes.broadcast(dimension_spec.to_vec()).unwrap().to_owned();

        // Reset the |1> state of the new hidden bit to zero
        let zeros = ArrayD::zeros(vec![2; hidden_bit_index - 1]);
        new_amplitudes.index_axis_mut(Axis(hidden_bit_index),1).assign(&zeros);
        
        // Swap our newly-zeroed hidden bit into place
        new_amplitudes.swap_axes(hidden_bit_index, bit_index);

        
        self.amplitudes = Box::new(new_amplitudes);
        self.hidden_qubit_count += 1;
    }



    // Despite the language spec talking about tracing over qubits, There isn't a way to implement 
    // reset as a projection in the formalism I've chosen. There is an implied measurement here, which means
    // I need to keep amplitudes around and defer that until the end such that probabilities for every possible
    // situation are fundamentally derivable from the amplitudes
    //
    // In the amplitude formalism the correct way to handle a reset is to introduce a new qubit in 
    // the |0> state and swap it into the position of the target qubit. Any other qubits
    // that are entangled with the target qubit will stay entangled with the now hidden qubit
    // and the squared norm logic in the probability() function will function correctly.
    //
    // This has the nice side effect that there's only one place in the whole program where we
    // actually take a trace, and that trace has the same meaning no matter how weird/crazy the 
    // circuit gets
    fn reset(&mut self, target_qubit:usize) {

        self.hide_bit(target_qubit);

        if self.trace { println!("{}", self.to_string()) };
    }


    // How to do conditionals
    // Collapse the relevant axes to the "only" values that make sense. i.e. those that satisfy the conditional
    // amplitudes still needs to have the same number of indicies so that indexing into it from all the other
    // handlers still works

    // We need to save off a handle to the original amplitudes somehow, so this can be undone
    fn start_conditional (&mut self, creg: usize, size: usize, value: u64) {

        // Save off the original amplitudes
        let amplitudes = self.amplitudes.to_owned();

        // Make a full copy for use during the condition
        let tensor_shape = vec![2; self.bit_count()];
        self.amplitudes = Box::new(ArrayD::<Complex32>::zeros(tensor_shape));
        self.amplitudes.assign(&amplitudes);

        // Make a list of which indicies are pinned and to what values
        let cbit_index = self.qubit_count + creg;
        let mut pinned_values:Vec<Option<bool>> = vec![None; cbit_index];
        for bit_index in 0..size {
            match (value >> bit_index) & 1 {
                0 => pinned_values.push(Some(false)),
                1 => pinned_values.push(Some(true)),
                _ => panic!("bitwise and of a number with 1 didn't give 1 or zero."),
            }
        }

        // Contract amplitudes along those values so gates can only access the elements
        // where the condition is true. Leave the axes we're pinning as length-1 so 
        // counting tensor axes still works.
        for (index, pinned_value) in pinned_values.iter().enumerate() {
            match pinned_value {
                Some(true) => {self.amplitudes.collapse_axis(Axis(index), 1);},
                Some(false) => {self.amplitudes.collapse_axis(Axis(index), 0);},
                None => (),
            }
        }

        // Save our settings in the object for when the condition is turned off
        self.amplitudes_save = Some(amplitudes);
        self.pinned_values = Some(pinned_values);

        if self.trace { println!("{}", self.to_string()) };
    }


    // Undo the amplitudes and put things back where they were
    pub fn end_conditional (&mut self) {
        // Recover the values we need to undo the condition.
        // Take ownership of the old amplitudes so we can put them in the right place.
        let pinned_values = match self.pinned_values.as_mut() {
            Some(pinned_values) => pinned_values,
            None => panic!("No pinned values, probably means a condition was ended without starting one."),
        };
        let mut amplitudes_saved = match self.amplitudes_save.as_mut() {
            Some(amplitudes_save) => amplitudes_save,
            None => panic!("No saved amplitudes, probably means a condition was ended without starting one."),
        }.to_owned();

        // Create a view of the saved amplitudes that has axes removed to the same slice the condition used.
        // At the same time remove the 1-length axes we left in self.amplitudes so that the two are the same shape.
        let mut amplitudes_saved_view = amplitudes_saved.view_mut();
        let mut amplitudes = *self.amplitudes.to_owned();

        let mut indexes_removed:usize = 0;
        for (index_iter, pinned_value) in pinned_values.iter().enumerate() {
            let index = index_iter - indexes_removed;
            match pinned_value {
                // take amplitudes (which was subject to conditional gates) and 
                // remove all its zero-length axes
                Some(true) => {
                    amplitudes.index_axis_inplace(Axis(index),0);
                    amplitudes_saved_view.index_axis_inplace(Axis(index), 1);
                    indexes_removed += 1;
                },
                Some(false) => {
                    amplitudes.index_axis_inplace(Axis(index), 0);
                    amplitudes_saved_view.index_axis_inplace(Axis(index), 0);
                    indexes_removed += 1;
                },
                // For none we leave the index alone and deal with the next index.
                None => (),
            }
        }

        // Now that the view on the saved amplitudes and self.amplitudes (which was altered during the condition) 
        // are the same shape. Assign the amplitudes to the slice of saved amplitudes, so that amplitudes_saved has 
        // the full set after the condition
        amplitudes_saved_view.assign(&amplitudes);

        // Put our saved amplitudes back now and return the object to normal mode
        self.amplitudes = amplitudes_saved;
        self.amplitudes_save = None;
        self.pinned_values = None;

        if self.trace { println!("{}", self.to_string()) };
    }
   

    // This is the core operation for converting amplitudes to probabilities 
    // it pins a bunch of c/qubit values, this part does not distinguish between c/q or hidden bits
    //
    // Caller must pass in the right order and no more than the total number of bits
    // function does a normsquared operation across all of them.
    //
    // Normalization is such that no values are pinned by caller, it should
    // always return 1
    pub fn probability(& self, pinned_values: &Vec<Option<bool>>) -> f32 {
        // Slice amplitudes down to only the relevant values
        let mut view = self.amplitudes.view();

        for (index, pinned_value) in pinned_values.iter().enumerate() {
            match pinned_value {
                // For true/false change our slice to be along that index being the desired value
                Some(true) => {view.collapse_axis(Axis(index), 1);},
                Some(false) => {view.collapse_axis(Axis(index), 0);},
                // For none we leave the index alone and deal with the next index
                None => (),
            }
            //println!("{:?}", view);
        }

        // Now view is an indexed object that lacks all the indicies that were pinned by caller
        // So we take the squared norm and return it.
        Self::norm_sqr(&view)
    }

    // Hand back the global probability that a single bit is set to the given value
    pub fn bit_probability(& self, bit_index: usize, value: bool) -> f32 {
        let mut pinned_values = vec![None; self.bit_count()];
        pinned_values[bit_index] = Some(value);
        self.probability(&pinned_values)
    }


    fn norm_sqr<T: Num + Clone> (amplitudes:&ArrayViewD<Complex<T>>) -> T {
        amplitudes.fold(T::zero(), |accum:T, elem:&Complex<T>| -> T{
            accum + elem.norm_sqr()
        })
    }    

    pub fn to_string(& self) -> String{
        self.to_string_internal(&(*self.amplitudes).view())
    }


    fn to_string_internal(&self, amplitudes:& ArrayViewD<Complex32>) -> String{
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
                    out.push_str(format!("{:.3} ", amplitude.re).as_str());
                } else if amplitude.re.is_zero() {
                    out.push_str(format!("{:.3} i ", amplitude.im).as_str());
                } else if amplitude.im.is_sign_negative(){
                    out.push_str(format!("({:.3} - {:.3} i)", amplitude.re, -amplitude.im).as_str());
                } else {
                    out.push_str(format!("({:.3} + {:.3} i)", amplitude.re, amplitude.im).as_str());
                }

                if self.debug {
                    out.push_str(Self::state_to_string(&indicies.as_array_view()).as_str());
                } else {
                    out.push_str(self.state_to_string_as_registers(&indicies.as_array_view()).as_str());
                }

                //println!("{:?} {:?}", indicies.as_array_view(), amplitude);
            }
        }

        if out.len() == 0 {
            out.push_str("0");
        }

        out
    }

    fn state_to_string (indicies: &ArrayView1<usize>) -> String{
        let mut out = String::new();
        out.push_str("|");
        for bit in indicies {
            out.push_str(format!("{:.3}", bit).as_str());
        }
        out.push_str(">");
        out
    }

    fn state_to_string_as_registers( &self, indicies: &ArrayView1<usize>) -> String {

        // Convert our indicies to Vec<Option<bool>> all with Some() so we can use the helper
        let pinned_bits = indicies.map( |bit| {
            match bit { 
                0 => Some(false), 
                1 => Some(true), 
                _ => panic!("Index specification had a value not one or zero.")}
        }).to_vec();

        let mut register_strings = Vec::new();

        register_strings.push(Self::register_printable_helper(&pinned_bits, 0, &self.qregister_extents));
        register_strings.push(Self::register_printable_helper(&pinned_bits, self.qubit_count, &self.cregister_extents));

        if self.hidden_qubit_count != 0 {
            let hidden_bit_extent: HashMap<String, (usize, usize)> = HashMap::from([(String::from("hidden"), (0, self.hidden_qubit_count))]);
            register_strings.push(Self::register_printable_helper(&pinned_bits, self.cbit_count + self.qubit_count, &hidden_bit_extent));
        }

        // concatenate everything and return
        format!("|{}>",register_strings.join(",")).to_string()
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
    // This is essentially a compile-time static proof that errors are unreachable 
    // in this dispatch code
    type Error = std::convert::Infallible;

    fn initialize(&mut self, qubits: &[Symbol], cbits: &[Symbol]) -> Result<(), Self::Error> {
        self.data.initialize(qubits.to_vec(), cbits.to_vec());
        Ok(())
    }

    fn write_cx(&mut self, copy: usize, xor: usize) -> Result<(), Self::Error> {
        if self.data.trace { println!("cx control qubit index: {copy} target qubit index: {xor}") };
        self.data.apply_cx(copy, xor);
        Ok(())
    }

    fn write_u(&mut self, theta: Value, phi: Value, lambda: Value, reg: usize) -> Result<(), Self::Error> {
        if self.data.trace { println!("u({theta}, {phi}, {lambda}) on qubit index {reg}") };
        self.data.apply_u(theta, phi, lambda, reg);
        Ok(())
    }

    fn write_opaque(&mut self, name: &Symbol, _: &[Value], _: &[usize]) -> Result<(), Self::Error> {
        if self.data.trace { println!("opaque gate {}", name) };
        Ok(())
    }

    fn write_barrier(&mut self, _: &[usize]) -> Result<(), Self::Error> {
        if self.data.trace { println!("barrier, does nothing") }
        Ok(())
    }

    fn write_measure(&mut self, from: usize, to: usize) -> Result<(), Self::Error> {
        if self.data.trace { println!("measure qubit index: {} into cbit index: {}", from, to) };
        self.data.measure(from, to);
        Ok(())
    }

    fn write_reset(&mut self, reg: usize) -> Result<(), Self::Error> {
        if self.data.trace { println!("reset qubit index: {reg}") };
        self.data.reset(reg);
        Ok(())
    }

    fn start_conditional(&mut self, reg: usize, count: usize, value: u64) -> Result<(), Self::Error> {
        if self.data.trace { println!("Begin Conditional segment cbit index {reg} through {} == {value}", reg+count)};
        self.data.start_conditional(reg, count, value);
        Ok(())
    }

    fn end_conditional(&mut self) -> Result<(), Self::Error> {
        if self.data.trace { println!("End Conditional segment") };
        self.data.end_conditional();
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

fn checked_div(lhs: Value, rhs: Rational64) -> Option<Value>{
    Some(Value {
        a: lhs.a.checked_div(&rhs)?,
        b: lhs.b.checked_div(&rhs)?,
    })
}