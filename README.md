# SPMR-MATLAB

**Author: Nicholas Hu ([nicholas.hu@alumni.ubc.ca](mailto:nicholas.hu@alumni.ubc.ca), [@nick-hu](https://github.com/nick-hu))**

SPMR: A family of saddle-point minimal residual solvers [[Estrin and Greif, 2018]]

## Installation

Install the [SPMR-MATLAB toolbox](SPMR-MATLAB.mltbx?raw=true).

## Getting started

To familiarize yourself with the functions and examples available, type `help SPMR-MATLAB` in the MATLAB command window.
In each of the examples, the results will appear in a struct called `result` after the example is run. To access (for
instance) the approximate solution *x*, type `result.x`. Descriptions of the fields of the `result` struct can be found
by typing `help` followed by the name of the SPMR function used (e.g., `help spmr_sc`).

[Estrin and Greif, 2018]: https://epubs.siam.org/doi/10.1137/16M1102410