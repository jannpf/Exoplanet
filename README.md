## General

This branch contains modified code to handle gaussian quasiperiodic noise to account for stellar activity. Experimental results for nu oph are also present.

## Prerequisites
* Configured c++ environment (g++/gcc etc.). You also need to clone [DNest4](https://github.com/eggplantbren/DNest4), then build it by executing make in the `code` directory. `gsl` library also needs to be installed.
* Update Makefile with correct full paths to dnest4 and gsl libraries.
* Configured python environment. You might need to install some libraries, including dnest4 for python: `pip install -U dnest4`

## Run locally
1. Create orbits for orbit lookup `cd orbits && python3 orbit_lookup.py`. This also produces `kepler.txt` file for visualizations.
2. Compile by running `make`. If you need to update a dataset, you should update it in main.cpp and then recompile.
3. Run the compiled executable. You can also multithread, e.g. `./main -t 12`. Sampling options can be configured in the `OPTIONS` file.

## Available datasets
1. `nu_oph.txt`
2. `fake_data_like_nu_oph.txt`. True rv values are also available
3. `new_data.txt` (HD 142)
