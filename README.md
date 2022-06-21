# ASCV Homotopy
## Homotopy Techniques for Analytic Combinatorics in Several Variables
ACSVHomotopy is a Julia package for studying the power series expansions of multivariate rational functions.
It is built using the functionality of [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/).

## Usage
To run a notebook with ACSVHomotopy:
- Clone this repository.
- Create a new julia notebook inside the `ACSVHomotopy` folder and put
```
using Pkg
Pkg.activate(".")
using ACSVHomotopy
```
inside the first cell.

The package reexports the interface of `HomotopyContinuation` to the user.
Polynomials can be constructed by declaring variables with `@polyvar` and including them into a Julia expression.
The Package exports the functions `find_min_crits`, `find_min_crits_comb` and `leading_asymptotics`.

### Example
Input:
```julia
using Pkg
Pkg.activate(".")
using ACSVHomotopy

@polyvar x y
h = 1-x-y
minimal = find_min_crits_comb(h)
leading_asymptotics(1, h, minimal)
```
Output:
```
(0.25 + 0.0im)^(-n) n^(-0.5) ((0.46065886596178074 + 0.0im))
```



