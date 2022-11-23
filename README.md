# README

Goal here:  Write a program which takes a circuit specification in an existing quantum circuit language and runs it on a classical computer. Potential extension: Model noise within the simulation, possibly with some physical model of the noise. Several programs like this exist in the world with fairly large feature sets; however, for this idea I'd want to code my own from a blank slate. 

Deliverables would be program, source code, readme & getting started level documentation, screen-capture demo videos, and runnable sample circuits that the simulator can execute put up on a website. 

Decisions:
- Rust is the base language
- I'm using QASM-rust https://github.com/libtangle/qasm-rust (Backup plan: Openqasm) for parsing so as not to waste time
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