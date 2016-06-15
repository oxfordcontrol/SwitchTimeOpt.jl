module SwitchTimeOpt

  # Import Necessary Modules
  using MathProgBase
  using ODE


  export createsto,
         solve!,
         gettau,
         getdelta,
         gettaucomplete, getdeltacomplete,
         getobjval, getstat, getsoltime,
         simulatelinearized, simulateinput,
         simulate


  # Include Source Files
  include("types.jl")
  include("interface.jl")
  include("nlp.jl")
  include("simulation.jl")

end
