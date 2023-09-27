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

1.  For single kernel, load the package by

<!---->

    Needs["TBMethod`"]

1.  For parallel computation, load it by

<!---->

    Needs["TBMethod`"]
    ParallelNeeds["TBMethod`"]

1.  Check the installation by

<!---->

    Scan[Echo @* Information] @ {"TBMethod`MDConstruct`*", "TBMethod`EigenSpect`*", "TBMethod`LGFF`*", "TBMethod`DataVisualization`*"}

and four lists of functions should be indexed out.

### Uninstallation

- **Run** `PacletUninstall["TBMethod"]`

## Functionalities

-

- Highlights: fully employment of the nearest neighbor searching algorithm, so that the total computation complex tends to O(n) in
 - Generation of Hamiltonian matrices,
 - Adaptive partition of central scattering region, and
 - Calculation of transport related quantities.

## Documentation

Under construction
