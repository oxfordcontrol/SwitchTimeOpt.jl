using SwitchTimeOpt
using Test
using Ipopt

include("linsys.jl")
linsystest()

include("nonlinsys.jl")
nonlinsystest()
