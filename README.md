# SwitchTimeOpt.jl
[![SwitchTimeOpt](http://pkg.julialang.org/badges/SwitchTimeOpt_0.6.svg)](http://pkg.julialang.org/?pkg=SwitchTimeOpt&ver=0.6)

[![Build Status](https://travis-ci.org/oxfordcontrol/SwitchTimeOpt.jl.svg?branch=master)](https://travis-ci.org/oxfordcontrol/SwitchTimeOpt.jl)



**SwitchTimeOpt.jl** is a [Julia](https://github.com/JuliaLang/julia) package to easily define and efficiently solve switching time optimization (STO) problems for linear and nonlinear systems. SwitchTimeOpt.jl supports a wide variety of nonlinear solvers through [MathProgBase.jl](https://github.com/JuliaOpt/MathProgBase.jl) interface such as [Ipopt](https://github.com/JuliaOpt/Ipopt.jl), [KNITRO](https://github.com/JuliaOpt/KNITRO.jl), [NLopt](https://github.com/JuliaOpt/NLopt.jl).


**Documentation** The complete interface documentation can be found [here](http://switchtimeoptjl.readthedocs.io/en/latest/).


## Citing this package

If you use SwitchTimeOpt.jl for published work, we encourage you to cite the following [paper](http://arxiv.org/abs/1608.08597):
```
@article{2016arXiv160808597S,
  title = {Second-Order Switching Time Optimization for Switched Dynamical Systems},
  author = {{Stellato}, B. and {Ober-Bl{\"o}baum}, S. and {Goulart}, P.},
  year = {2017},
  journal = {IEEE Transactions on Automatic Control (To appear)}
}
```
