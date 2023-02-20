# TBMethod

## Installation & Uninstallation
### Installation
1. **Download** the latest `"TBMethod-<\*version #\*>.paclet"` file to local machine;
2. **Run** `PacletInstall["<*path-to-download*>/TBMethod-<*version #*>.paclet"]`.

### Installation Test
1. For single kernel, load the package by 
```
Needs["TBMethod`"]
```
2. For parallel computation, load it by
```
Needs["TBMethod`"]
ParallelNeeds["TBMethod`"]
```	

_Choose one between Step 1 and Step 2_

3. Check the installation by
```
Scan[Echo @* Information] @ {"TBMethod`MDConstruct`*", "TBMethod`EigenSpect`*", "TBMethod`LGFF`*", "TBMethod`DataVisualization`*"}
```
and four lists of functions should be indexed out.

### Uninstallation
- **Run** `PacletUninstall["TBMethod"]`

## Functionalities
Under construction

## Documentation
Under construction
