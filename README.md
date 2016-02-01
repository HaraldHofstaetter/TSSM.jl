# TSSM.jl
## Installation
```julia
Pkg.clone("https://github.com/HaraldHofstaetter/TSSM.jl")
Pkg.build("TSSM")
```
##Examples
To get easy access to the examples, make a symbolic link in the home directory:
```julia
symlink(joinpath(homedir(), ".julia/v0.4/TSSM/examples/"), joinpath(homedir(), "TSSM_examples"))
```
Then 'TSSM_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [GroundstateDemo.ipynb](https://github.com/HaraldHofstaetter/TSSM.jl/blob/master/examples/GroundstateDemo.ipynb)
+ [TimeStepperDemo.ipynb](https://github.com/HaraldHofstaetter/TSSM.jl/blob/master/examples/TimeStepperDemo.ipynb)
