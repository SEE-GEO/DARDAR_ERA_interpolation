# Study to retrieve Ice Water Path

This repository contains code and documentation for a study aiming at retrieval
of IWP from GMI data. 

## Dependencies

For preparing the data required to run ARTS simulations, we are using the [pansat](https://github.com/SEE-MOF/pansat) package.
Install it manually as :

````
git clone https://github.com/SEE-MOF/pansat
cd pansat
pip install -e .
````

Note the `-e` flag for the `pip` command. This is allows the installed pansat package to be updated, when
pulling a new version of pansat

