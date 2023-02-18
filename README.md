# TBMethod

## Installation & Uninstallation
### Installation
1. **Download** the latest `"TBMethod-<\*version #\*>.paclet"` file to local machine;
2. **Run** `PacletInstall["<*path-to-download*>/TBMethod-<*version #*>.paclet"]`.

### Installation Test
1. Load the package by 
```
Needs["TBMethod`"]
```
3. For parallel computation, load it by
```
Needs["TBMethod`"]
ParallelNeeds["TBMethod`"]
```	
4. Check the installation by
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
