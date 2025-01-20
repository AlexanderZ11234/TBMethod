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

## Functionalities

- Highlights: sufficient employment of the [NNS (nearest neighbor search)](https://en.wikipedia.org/wiki/Nearest_neighbor_search) algorithm, so that the total computation complexity tends to be fine as:
    - Model construction linear in system's size $ \text{O}(n) $:
        - Generation of Hamiltonian matrices,
        - Adaptive partition of central scattering region

    - Calculation of transport related quantities:
        - 5-terminal Hall calculation in $ \text{O}(n^{1.7}) $

## Documentation

<a href="#" class="magic-button" title="Onsite testable"> _MMA-style_ </a> documentation under construction

## Incomplete References

- Topological Models:
    1. PRL 61, 2015 (2988).
    1. PRL 95, 146802 (2005).
    1. PRL 95, 226801 (2025).
    1. PRB(R) 82, 161414 (2010).
    1. PRL 124, 136403 (2020).
    1. PRL 124, 166804 (2020).

- Green's function:
    1. Datta, _Electronic Transport in Mesoscopic Systems_.
    1. Datta, _Quantum Transport: Atom to Transistor_.
    1. PRL 97, 066603 (2006).
    1. Nanotechnology 18, 435402 (2007).
    1. PRB 91, 125408 (2015).
    1. PRB 97, 165405 (2018).
