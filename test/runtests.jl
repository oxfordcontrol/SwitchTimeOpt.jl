using SwitchTimeOpt
using Test
using SafeTestsets

@safetestset "Linear system" begin include("./linsys.jl") end
@safetestset "Nonlinear system" begin include("./nonlinsys.jl") end
