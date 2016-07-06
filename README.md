# SwitchTimeOpt.jl

[![Build Status](https://travis-ci.org/bstellato/SwitchTimeOpt.jl.svg?branch=master)](https://travis-ci.org/bstellato/SwitchTimeOpt.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/bstellato/SwitchTimeOpt.jl?branch=master&svg=true)](https://ci.appveyor.com/project/bstellato/switchtimeopt-jl/branch/master)



**SwitchTimeOpt.jl** is a [Julia](https://github.com/JuliaLang/julia) package to easily define and efficiently solve switching time optimization (STO) problems for linear and nonlinear systems. SwitchTimeOpt.jl supports a wide variety of nonlinear solvers through [MathProgBase.jl](https://github.com/JuliaOpt/MathProgBase.jl) interface such as [Ipopt](https://github.com/JuliaOpt/Ipopt.jl), [KNITRO](https://github.com/JuliaOpt/KNITRO.jl), [NLopt](https://github.com/JuliaOpt/NLopt.jl).


**Documentation** The complete interface documentation can be found [here](http://switchtimeoptjl.readthedocs.io/en/latest/).


## Todo List
- [ ] Solve NaN problem first-order method IPOPT inside Fishing Problem
- [X] Fix Linearization Grid for nonlinear systems
- [ ] Write State Constraints for Linear Systems within the Grid
- [ ] Write State Constraints for Nonlinear Systems within the Grid

<!-- ## Installation

You can install the package by running

    julia> Pkg.clone("git://github.com/bstellato/SwitchTimeOpt.jl.git")

This does not install any nonlinear solvers. If you don’t have a nonlinear solver installed already, you will want to install a solver such as [Ipopt](https://github.com/JuliaOpt/Ipopt.jl) by running:

    julia> Pkg.add("Ipopt")


## Usage

### Linear Switching Time Optimization


### Nonlinear Switching Time Optimization -->
