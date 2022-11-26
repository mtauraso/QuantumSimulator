# README

Goal here:  Write a program which takes a circuit specification in an existing quantum circuit language and runs it on a classical computer. Potential extension: Model noise within the simulation, possibly with some physical model of the noise. Several programs like this exist in the world with fairly large feature sets; however, for this idea I'd want to code my own from a blank slate. 

Deliverables would be program, source code, readme & getting started level documentation, screen-capture demo videos, and runnable sample circuits that the simulator can execute put up on a website. 

Decisions:
- Rust is the base language
- Not using QASM-rust because it doesn't support parsing all of QASM 2.0. Openqasm's AST is much nicer, and supports the expression language implicit in QASM
- Nalgebra for linal (and num/ndarray i think as deps?)
- Input will be OpenQASM 2.0 https://github.com/openqasm/openqasm/tree/OpenQASM2.x (and the paper in my zotero)
- Readme level docs and "how to run" will be on github repo, will directly link PDF
- Final paper with citations will be in this repo as a tex file. pdf will include hyperlink support


Program flow:
- Read in the qasm file
- Execute it: There may need to be some different execution modes 
- You see command-line output (Output format may be flaggable)

Open questions:
- Final docs: (Ideas: github.io style website, sequence of .md files in repo, A pdf/tex paper) md files and pdf/tex paper seem both easiest and most legible to academia
- Specifying initial input states
- How do you get the output? (Should be some form of amplitudes? or probability of measurements? or bloch spheres?)
- What is the noise model? 
- Do we allow specification of symbolic input states?
- Do I want particular pretty output?

Out of scope:
- Input circuit GUI

TODO:
- Separate out the visitor logic into separate files
- Use the visitor logic to create registers
- Data model the state vector of all quantum registers
- Use the visitor logic to evolve quantum registers for basic test program
- Look through spec and implement all (implementable) primitives

- Read the QASM spec, and consider user interface
	- An idea: Circuit is displayed graphically, with ability to expand/contract the gates. Circuit input is the openqasm file, graphic output is only a viewer of circuit and results 

- Write down the most cut down basic user interface that would work
- Write down 3-5 immediate extensions on that basic form. Ideally no more than 1 level of "A builds on B" type dependencies. Once basic is done, 3 ideas should be actionable