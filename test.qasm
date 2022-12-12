OPENQASM 2.0;

// Standard header for usual gate set
include "qelib1.inc";

// copied from standard header
//gate u2(phi, lambda) q { U(pi/2, phi, lambda) q; }

// Clifford gate: Hadamard (exists in standard header)
//gate h a { u2(0,pi) a; }


qreg q[2];
creg c[1];
h q;
cx q[0], q[1];

measure q[0] -> c[0];
reset q[1];

if (c == 0) h q[1];

//measure q[1] -> c[0];


// qreg r[2];
// qreg q[3];
// creg c[3];

// x q[0];

// h q[1];
// CX q[1], q[2];
// CX q[0], q[1];

// measure q -> c;

// h r[1];

// CX r[1], q[2];
