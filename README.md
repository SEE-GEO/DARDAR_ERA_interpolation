# Interpolate ERA5 to DARDAR grid

This repository contains code and documentation for interpolating ERA5 data to DARDAR grid.
Preparing the data required to run ARTS simulations.

## Dependencies

We are using the [pansat](https://github.com/SEE-MOF/pansat) package.
Install it manually as :

````
git clone https://github.com/SEE-MOF/pansat
cd pansat
pip install -e .
````

Note the `-e` flag for the `pip` command. This is allows the installed pansat package to be updated, when
pulling a new version of pansat

To use pansat, Copy the identities.json file in pansat/test/unit_tests/test_data to ~/.config/pansat
When prompted for password, enter: not_a_secret
An environment variable PANSAT_PASSWORD can also be set with this password
