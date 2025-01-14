# TBMethod

## Installation & Uninstallation

### Installation

#### Offline

1.  **Download** the latest `"TBMethod-<\*version #\*>.paclet"` file to local machine;

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

- **Run** `PacletUninstall["TBMethod"]`

## Functionalities

- Highlights: fully employment of the [NNS (nearest neighbor search)](https://en.wikipedia.org/wiki/Nearest_neighbor_search) algorithm, so that the total computation complexity tends to be fine as:
    - Model construction linear in system's size $ \text{O}(n) $:
        - Generation of Hamiltonian matrices,
        - Adaptive partition of central scattering region

    - Calculation of transport related quantities:
        - 5-terminal Hall calculation in $ \text{O}(n^{1.7}) $

## Documentation

<a href="#" class="magic-button" title="Onsite testable"> _MMA-style_ </a> documentation under construction
