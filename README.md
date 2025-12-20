# TBMethod

## Installation & Uninstallation

### Installation

Programming environment version: the latest the best

#### Offline

1.  **Download** the latest `"TBMethod-<\*version #\*>.paclet"` file to one's local machine;

2.  **Run** `PacletInstall["<*path-to-download*>/TBMethod-<*version #*>.paclet"]`.

#### Online

- Outstanding

<!--
**Run** `PacletInstall["https://github.com/AlexanderZ11234/TBMethod/releases/download/0.2.1/TBMethod-0.2.1.paclet"]`
-->

### Installation Test (1 → 3 or 2 → 3)

1. For single kernel, load the package by

<!---->

    Needs["TBMethod`"]

2. For parallel computation, load it by

<!---->

    Needs["TBMethod`"]
    ParallelNeeds["TBMethod`"]

3. Check the installation by

<!---->

    Scan[Echo @* Information] @ {"TBMethod`MDConstruct`*", "TBMethod`EigenSpect`*", "TBMethod`LGFF`*", "TBMethod`DataVisualization`*"}

and four lists of functions should be indexed out.

### Uninstallation

- **Run** `PacletUninstall["TBMethod"]` for uninstallation or reinstallation.

## Functionality Highlights

- External degree of freedom (real-space coordinate): sufficient employment of the [NNS (nearest neighbor search)](https://en.wikipedia.org/wiki/Nearest_neighbor_search) algorithm, so that the total computation complexity tends to be fine as:
    - Model construction linear in system's size $ \text{O}(n) $:
        - Generation of Hamiltonian matrices,
        - Adaptive partition of central scattering region

    - Calculation of transport related quantities:
        - 5-terminal Hall calculation in $ \text{O}(n^{1.7}) $

- Internal degree of freedoms: spin, atomic orbital, (BdG) particle-hole, (Floquet) photon block, and lattice vibration polarization

- Workflow coordinated with [DeePTB](https://github.com/deepmodeling/DeePTB) on Slater-Koster model construction and transport calculation with nonidentity overlapping matrices

## Documentation

<a href="#" class="magic-button" title="Onsite testable"> _MMA-style_ </a> documentation under construction

Coorperation is highly welcome.

A tutorial in [Zhihu Column](https://www.zhihu.com/column/c_1954932040861450487) is also under compilation.


<details open>
<summary>

## Related Publications
</summary>

1. [npj Quant. Mater. **10**, 48 (2025)](https://www.nature.com/articles/s41535-025-00768-1).
1. [Phys. Rev. B **111**, 085137 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.085137).
1. [Phys. Rev. B **111**, 155303 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.155303).
1. [Phys. Rev. Lett. **133**, 246606 (2024)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.246606).
1. [Phys. Rev. Lett. **131**, 086601 (2023)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.131.086601).
1. [Phys. Rev. B **107**, 075303 (2023)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.075303).
1. [Phys. Rev. B (Letter) **106**, L201407 (2022)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.L201407).
1. [Front. Phys. **17**, 63503 (2022)](https://link.springer.com/article/10.1007/s11467-022-1185-y).
1. [Appl. Phys. Lett. **120**, 084002 (2022)](https://pubs.aip.org/aip/apl/article-abstract/120/8/084002/2833231/In-plane-magnetization-and-electronic-structures?redirectedFrom=fulltext).
1. [Chin. Phys. Lett. **39**, 017302 (2022)](https://cpl.iphy.ac.cn/article/doi/10.1088/0256-307X/39/1/017302).
1. [Phys. Rev. B **101**, 235432 (2020)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.235432).
1. [Phys. Rev. B **100**, 205408 (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.205408).
1. [Phys. Rev. B **95**, 045424 (2017)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.045424).
</details>


## Incomplete References

<details>
<summary>

### Topological Models & Characterization

</summary>

1. Bernevig, [_Topological Insulators and Topological Superconductors_](https://press.princeton.edu/books/hardcover/9780691151755/topological-insulators-and-topological-superconductors?srsltid=AfmBOop9JnAo53v7Hn3ErPpR2uf3vW0JLPykFSNWSK_QoP1xjsDuKoMG), PUP, 2013.
1. Shen, [_Topological Insulators: Dirac Equation in Condensed Matters_](https://link.springer.com/book/10.1007/978-981-10-4606-3), Springer, 2017.
1. [Phys. Rev. Lett. **61**, 2015 (1988)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.61.2015).
1. [Phys. Rev. Lett. **95**, 146802 (2005)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.146802).
1. [Phys. Rev. Lett. **95**, 226801 (2005)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.226801).
1. [Phys. Rev. B **82**, 161414(R) (2010)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.82.161414).
1. [Phys. Rev. B **84**, 075119 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.075119).
1. [Phys. Rev. Lett. **112**, 037001 (2014)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.112.037001).
1. [Phys. Rev. B **95**, 195102 (2017)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.195102).
1. [Phys. Rev. B **95**, 245433 (2017)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.245433).
1. [Phys. Rev. Lett. **124**, 136403 (2020)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.136403).
1. [Phys. Rev. Lett. **124**, 166804 (2020)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.166804).

</details>

<details>
<summary>

### Lattice Green's Function Formalism

</summary>

1. Datta, [_Electronic Transport in Mesoscopic Systems_](https://www.cambridge.org/core/books/electronic-transport-in-mesoscopic-systems/1E55DEF5978AA7B843FF70337C220D8B), CUP, 1995.
1. Datta, [_Quantum Transport: Atom to Transistor_](https://www.cambridge.org/core/books/quantum-transport/E96BE74AACD59A03A7D6A7F7DACDFB71), CUP, 2005.
1. Wimmer, [_Quantum transport in nanostructures: From computational concepts to spintronics in graphene and magnetic tunnel junctions_](https://epub.uni-regensburg.de/12142/), Ph.D. Dissertation, Universität Regensburg, 2008.
1. Qiao, [_Charge and Spin Transport in Two-Dimensional Mesoscopic Systems_](https://hub.hku.hk/handle/10722/55540), Ph.D. Dissertation, HKU, 2009.
1. Papior, [_Computational Tools and Studies of Graphene Nanostructures_](https://orbit.dtu.dk/en/publications/computational-tools-and-studies-of-graphene-nanostructures), Ph.D. Dissertation, TUD, 2016.
1. [J. Phys. F: Met. Phys. **14**, 1205 (1984)](https://iopscience.iop.org/article/10.1088/0305-4608/14/5/016).
1. [J. Phys. F: Met. Phys. **15**, 851 (1985)](https://iopscience.iop.org/article/10.1088/0305-4608/15/4/009).
1. [Phys. Rev. Lett. **97**, 066603 (2006)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.066603).
1. [Nanotechnology **18**, 435402 (2007)](https://iopscience.iop.org/article/10.1088/0957-4484/18/43/435402).
1. [Phys. Rev. B **83**, 085412 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.085412).
1. [Phys. Rev. B **91**, 125408 (2015)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.125408).
1. [Phys. Rev. B **97**, 165405 (2018)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.165405).
1. [Phys. Rev. B **100**, 195417 (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.195417).

</details>

<details>
<summary>

### Slater-Koster Method

</summary>

1. Saito, [_Physical Properties of Carbon Nanotubes_](https://www.worldscientific.com/worldscibooks/10.1142/p080?srsltid=AfmBOoosI-cgWaXJxEbpkiw1QPAPb82G87WuKIr6LAeeVNM8vWX1tifB#t=aboutBook), ICP, 1998.
1. [Phys. Rev. **94**, 1498 (1954)](https://journals.aps.org/pr/abstract/10.1103/PhysRev.94.1498).
1. [Phys. Rev. B **74**, 165310 (2006)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.74.165310).
1. [Phys. Rev. B **82**, 245412 (2010)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.82.245412).
1. [Nat. Commun. **15**, 6772 (2024)](https://www.nature.com/articles/s41467-024-51006-4).
1. [Phys. Rev. B **110**, 235130 (2024)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.110.235130). 
</details>

<details>
<summary>

### Other Theoretical Considerations
</summary>

1. [Z. Phys. **64**, 629 (1930)](https://link.springer.com/article/10.1007/BF01397213).
1. [Z. Phys. **80**, 763 (1933)](https://link.springer.com/article/10.1007/BF01342591).
1. [Phys. Rev. B **40**, 8169 (1989)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.40.8169).
1. [Phys. Rev. B **79**, 081406(R) (2009)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.79.081406).
1. [Phys. Rev. B **84**, 235108 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.235108).
1. [Phys. Rev. Lett. **114**, 056801 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.056801).
</details>