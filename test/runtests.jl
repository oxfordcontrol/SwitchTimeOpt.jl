using SwitchTimeOpt
using Test
using SafeTestsets

@safetestset "Linear system" begin include("./linsys.jl") end
@safetestset "Noninear system" begin include("./nonlinsys.jl") end
