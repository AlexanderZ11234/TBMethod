# TBMethod
 
This package consists of four subpackages responsible, respectively, for different group of functionalities.

- MDConstruct
	- Etymology: short for "MoDel Construction"
	- Major task: generation of matrix representation of the Hamiltonian for a lattice system in any dimension, any configuration (connection position and number) of terminals
		- Using [NNS (nearest neighbor search)](https://en.wikipedia.org/wiki/Nearest_neighbor_search) algorithm for matrix generation and adaptive partiton of the central scattering region
		- Supporting Slater-Koster method, both user-defined parameters and results from [DeePTB](https://github.com/deepmodeling/DeePTB)


- EigenSpect
	- Etymology: short for "Eigen Spectrum"
	- Major task: evaluation of the characteristic information for the given Hamiltonian matrices, e.g., bandstructures, energy spectrum, topological invariants (Chern number, WCC, Bott index, etc.), and so on

- LGFF
	- Etymology: short for "Lattice Green's Function Formalism"
	- Major task: manipulation of stuff related to nonequilibrium Green's function (NEGF), e.g., surface Green's function / self energy of a lead, Green's function of the central scattering region (CSR), transport coefficients (transmission and reflection), local density of states (LDOS) in either real space or the reciprical space, local current density vector (LCDV) field, etc.

- DataVisualization
	- Etymology: self-evident
	- Major task: visulization of the relavant data
