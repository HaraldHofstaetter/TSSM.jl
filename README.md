# TSSM.jl
## Installation
```julia
Pkg.clone("https://github.com/HaraldHofstaetter/TSSM.jl")
Pkg.build("TSSM")
```
##Examples
To get easy access to the examples, copy them into the home directory:
```julia
cp(joinpath(homedir(), ".julia/v0.4/TSSM/examples/"), joinpath(homedir(), "TSSM_examples"), remove_destination=true )
```
Then 'TSSM_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [SchroedingerDemo.ipynb](https://github.com/HaraldHofstaetter/TSSM.jl/blob/master/examples/SchroedingerDemo.ipynb)
+ [GroundstateDemo.ipynb](https://github.com/HaraldHofstaetter/TSSM.jl/blob/master/examples/GroundstateDemo.ipynb)
+ [TimeStepperDemo.ipynb](https://github.com/HaraldHofstaetter/TSSM.jl/blob/master/examples/TimeStepperDemo.ipynb)
+ [SolitonDemo.ipynb](https://github.com/HaraldHofstaetter/TSSM.jl/blob/master/examples/SolitonDemo.ipynb)
