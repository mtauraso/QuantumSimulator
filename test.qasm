OPENQASM 2.0;

// Standard header for usual gate set
include "qelib1.inc";

qreg q[2];
creg c[1];

h q;
cx q[0], q[1];

measure q[0] -> c[0];
reset q[1];

if (c == 0) h q[1];

measure q[1] -> c[0];
