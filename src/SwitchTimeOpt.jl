module SwitchTimeOpt

  # Import Necessary Modules
  using MathProgBase
  using ODE


  export createsto,
         solve!,
         gettau, gettaucomplete, getobjval, getstat, getsoltime,
         simulate, simulatelinearized


  # Include Source Files
  include("types.jl")
  include("interface.jl")
  include("nlp.jl")
  include("simulation.jl")

end
