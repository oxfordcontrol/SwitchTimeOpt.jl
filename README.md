# SwitchTimeOpt.jl
[![SwitchTimeOpt](http://pkg.julialang.org/badges/SwitchTimeOpt_0.4.svg)](http://pkg.julialang.org/?pkg=SwitchTimeOpt&ver=0.4)
[![SwitchTimeOpt](http://pkg.julialang.org/badges/SwitchTimeOpt_0.5.svg)](http://pkg.julialang.org/?pkg=SwitchTimeOpt&ver=0.5)

[![Build Status](https://travis-ci.org/OxfordControl/SwitchTimeOpt.jl.svg?branch=master)](https://travis-ci.org/OxfordControl/SwitchTimeOpt.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/bstellato/SwitchTimeOpt.jl?branch=master&svg=true)](https://ci.appveyor.com/project/bstellato/switchtimeopt-jl/branch/master)



**SwitchTimeOpt.jl** is a [Julia](https://github.com/JuliaLang/julia) package to easily define and efficiently solve switching time optimization (STO) problems for linear and nonlinear systems. SwitchTimeOpt.jl supports a wide variety of nonlinear solvers through [MathProgBase.jl](https://github.com/JuliaOpt/MathProgBase.jl) interface such as [Ipopt](https://github.com/JuliaOpt/Ipopt.jl), [KNITRO](https://github.com/JuliaOpt/KNITRO.jl), [NLopt](https://github.com/JuliaOpt/NLopt.jl).


**Documentation** The complete interface documentation can be found [here](http://switchtimeoptjl.readthedocs.io/en/latest/).


## Citing this package

If you use SwitchTimeOpt.jl for published work, we encourage you to cite the following [paper](http://arxiv.org/abs/1608.08597):
```
@article{2016arXiv160808597S,
  author = {{Stellato}, B. and {Ober-Bl{\"o}baum}, S. and {Goulart}, P.~J.},
  title = "{Second-Order Switching Time Optimization for Switched Dynamical Systems}",
  journal = {ArXiv e-prints},
  archivePrefix = "arXiv",
  eprint = {1608.08597},
  primaryClass = "math.OC",
  keywords = {Mathematics - Optimization and Control},
  year = 2016,
  month = aug
}
```
