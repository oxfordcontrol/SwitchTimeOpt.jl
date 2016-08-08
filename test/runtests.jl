using SwitchTimeOpt
using Base.Test
using Ipopt

include("linsys.jl")
linsystest()

include("nonlinsys.jl")
nonlinsystest()
