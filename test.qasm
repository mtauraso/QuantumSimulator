OPENQASM 2.0;

// Standard header for usual gate set
include "qelib1.inc";

// copied from standard header
//gate u2(phi, lambda) q { U(pi/2, phi, lambda) q; }

// Clifford gate: Hadamard (exists in standard header)
// gate h a { u2(0,pi) a; }


qreg q[2];
qreg j[2];
creg c[1];

h q[0];
CX q[0], q[1];
CX j, q;

measure q[1] -> c[0];