module SwitchTimeOpt

  # Import Necessary Modules
  using MathProgBase
  import OrdinaryDiffEq


  export stoproblem,
         setwarmstart!,setx0!,
         solve!,
         gettau,
         getdelta,
         gettaucomplete, getdeltacomplete,
         getobjval, getstat, getsoltime,
         getnobjeval, getngradeval, getnhesseval,
         simulatelinearized, simulateinput,
         simulate

  # Include Source Files
  include("types.jl")
  include("interface.jl")
  include("nlp.jl")
  include("simulation.jl")

end
