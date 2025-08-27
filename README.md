# TBMethod

## Installation & Uninstallation

### Installation

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

## Functionality Highlights:

- External degree of freedom (real-space coordinates): sufficient employment of the [NNS (nearest neighbor search)](https://en.wikipedia.org/wiki/Nearest_neighbor_search) algorithm, so that the total computation complexity tends to be fine as:
    - Model construction linear in system's size $ \text{O}(n) $:
        - Generation of Hamiltonian matrices,
        - Adaptive partition of central scattering region

    - Calculation of transport related quantities:
        - 5-terminal Hall calculation in $ \text{O}(n^{1.7}) $

- Internal degree of freedom: spin, atomic orbital, (BdG) particle-hole, (Floquet) photon block, and lattice vibration polarization

- Workflow coordinated with [DeePTB](https://github.com/deepmodeling/DeePTB) on Slater-Koster model construction and transport calculation with nonidentity overlapping matrices

## Documentation

<a href="#" class="magic-button" title="Onsite testable"> _MMA-style_ </a> documentation under construction

## Related Publications

1. npj Quant. Mater. 10, 48 (2025).
1. Phys. Rev. B 111, 085137 (2025).
1. Phys. Rev. Lett. 133, 246606 (2024).
1. Phys. Rev. Lett. 131, 086601 (2023).
1. Phys. Rev. B 107, 075303 (2023).
1. Phys. Rev. B (Letter) 106, L161413 (2022).
1. Front. Phys. 17, 63503 (2022).
1. Appl. Phys. Lett. 120, 084002 (2022).
1. Chin. Phys. Lett. 39, 017302 (2022).
1. Phys. Rev. B 101, 235432 (2020).
1. Phys. Rev. B 100, 205408 (2019).
1. Phys. Rev. B 95, 045424 (2017).



## Incomplete References

### Topological Models:
1. Shun-Qing Shen, _Topological Insulators_, Springer/WPC, 2006.
1. Phys. Rev. Lett. 61, 2015 (1988).
1. Phys. Rev. Lett. 95, 146802 (2005).
1. Phys. Rev. Lett. 95, 226801 (2005).
1. Phys. Rev. B 82, 161414(R) (2010).
1. Phys. Rev. Lett. 124, 136403 (2020).
1. Phys. Rev. Lett. 124, 166804 (2020).

### Lattice Green's Function Formalism:
1. Datta, _Electronic Transport in Mesoscopic Systems_, Cambridge/WPC, 2004.
1. Datta, _Quantum Transport: Atom to Transistor_, Cambridge/WPC, 2011.
1. Phys. Rev. Lett. 97, 066603 (2006).
1. Nanotechnology 18, 435402 (2007).
1. Phys. Rev. B 91, 125408 (2015).
1. Phys. Rev. B 97, 165405 (2018).

### Slater-Koster Method:
1. Saito and Dresselhauses, _Physical Properties of Carbon Nanotubes_, ICP/WPC, 2003.
1. Phys. Rev. 94, 1498 (1954).
1. Nat. Commun. 15, 6772 (2024).
1. Phys. Rev. B 74, 165310 (2006).
1. Phys. Rev. B 82, 245412 (2010).
1. Phys. Rev. B 110, 235130 (2024). 

