# Empty matrix and vector for function evaluation
const emptyfvec = Array(Float64, 0)
const emptyfmat = Array(Float64, 0, 0)


# Define Abstract NLP Evaluator for STO problem
abstract STOev <: MathProgBase.AbstractNLPEvaluator


# Linear Case
type linSTOev <: STOev
  # Parameters
  x0::Array{Float64,1}                    # Initial State x0
  nx::Int                                 # State dimension
  A::Array{Float64,3}                     # Dynamics Matrices
  N::Int                                  # Number of Switching Times
  t0::Float64                             # Initial Time
  tf::Float64                             # Final Time (from setup)
  tfdelta::Float64                        # Final Time (from current delta)
  Q::Array{Float64,2}                     # State Cost
  Qf::Array{Float64,2}                    # Final State Cost

  # Grid Variables
  ngrid::Int64                             # Number of grid points
  tgrid::Array{Float64,1}                  # Array of the grid
  tvec::Array{Float64,1}                   # Array of the complete grid with sw times
  tauIdx::Array{Int,1}                     # Array of the position of switching times in the complete grid vector tvec
  tgridIdx::Array{Int,1}                   # Array of the position of grid points in the complete grid vector tvec
  deltacomplete::Array{Float64,1}          # Array of the complete sw intervals
  ncons::Int64                             # Number of constraints per stage
  nconsf::Int64                            # Number of constraints at final stage

  # Precomputed Values
  V::Array{Complex{Float64}, 3}            # Dynamics Matrices decomp V
  invV::Array{Complex{Float64}, 3}         # Dynamics Matrices decomp V^-1
  D::Array{Complex{Float64}, 2}            # Dynamics Matrices decomp D
  isDiag::Array{Bool,1}                    # Vector of boolean variables stating if dynamics matrices are diagonalizable

  IndTril::Array{Int, 1}                   # Single Element of Lower Triangular Matrices (Hessian)
  Itril::Array{Int, 1}                     # Double Element Indeces of Lower Triangular Matrices (Hessian)
  Jtril::Array{Int, 1}                     # Double Element Indeces of Lower Triangular Matrices (Hessian)
  Ac::Array{Float64, 2}                    # Linear Constraints Matrix
  Acf::Array{Float64, 2}                    # Final Linear Constraints Matrix
  gsum::Array{Float64, 2}                  # Constraint sum of deltas
  Ig::Array{Int, 1}                        # Index Linear Constraints
  Jg::Array{Int, 1}                        # Index Linear Constraints
  Vg::Array{Float64, 1}                    # Value Linear Constraints
  # bg::Array{Float64, 1}                    # Constant Term Linear Constraints


  # Shared Data Between Functions
  deltafun_prev::Array{Float64, 1}            # Previous switching intervals for function evaluation
  deltagrad_prev::Array{Float64, 1}             # Previous switching intervals for gradient evaluation
  deltahess_prev::Array{Float64, 1}             # Previous switching intervals for hessian evaluation
  xpts::Array{Float64, 2}                  # States at Switching Times
  expMat::Array{Float64, 3}                # Matrix Exponentials
  Phi::Array{Float64, 4}                   # Matrix of State Transition Matrices
  M::Array{Float64, 3}                     # Integrals over Switching Intervals
  S::Array{Float64, 3}                     # S Matrices for each interval
  C::Array{Float64, 3}                     # C Matrices for each interval


  # Store Evaluations
  obj::Array{Float64, 1}                   # Store Objective Value
  deltaval::Array{Float64, 2}              # Store Iterates
  nobjeval::Int                            # Number of objective function evaluations
  ngradeval::Int                           # Number of gradient evaluations
  nhesseval::Int                           # Number of hessian evaluations
end

# Nonlinear Case
type nlinSTOev <: STOev
  # Parameters
  x0::Array{Float64,1}                     # Initial State x0
  nx::Int                                  # State dimension
  A::Array{Float64,3}                      # Linearized dynamics matrices
  N::Int                                   # Number of switching times
  t0::Float64                              # Initial Time
  tf::Float64                              # Final Time (from setup)
  tfdelta::Float64                         # Final Time (from current delta)
  Q::Array{Float64,2}                      # State Cost
  Qf::Array{Float64,2}                     # Final State Cost
  uvec::Array{Float64,2}                   # Array of switching inputs
  # nartsw::Int64                          # Number of artificial switches

  # Grid Variables
  ngrid::Int64                             # Number of grid points
  tgrid::Array{Float64,1}                  # Array of the grid
  tvec::Array{Float64,1}                   # Array of the complete grid with sw times
  tauIdx::Array{Int,1}                     # Array of the position of switching times in the complete grid vector tvec
  tgridIdx::Array{Int,1}                   # Array of the position of grid points in the complete grid vector tvec
  deltacomplete::Array{Float64,1}          # Array of the complete sw intervals including the grid ones


  # Nonlinear Dynamics and Derivatives Functions
  nonlin_dyn::Function
  nonlin_dyn_deriv::Function

  # Precomputed Values
  IndTril::Array{Int, 1}                  # Single Element Indeces of Lower Triangular Matrices (Hessian)
  Itril::Array{Int, 1}                    # Double Element Indeces of Lower Triangular Matrices (Hessian)
  Jtril::Array{Int, 1}                    # Double Element Indeces of Lower Triangular Matrices (Hessian)
  Ag::Array{Float64, 2}                   # Linear Constraints (Vector)
  Ig::Array{Int, 1}                       # Index Linear Constraint
  Jg::Array{Int, 1}                       # Index Linear Constraint
  Vg::Array{Float64, 1}                   # Value Linear Constraints
  bg::Array{Float64, 1}                   # Constant Term Linear Constraints

  # Shared Data Between Functions
  deltafun_prev::Array{Float64, 1}            # Previous switching intervals for function evaluation
  deltagrad_prev::Array{Float64, 1}             # Previous switching intervals for gradient evaluation
  deltahess_prev::Array{Float64, 1}             # Previous switching intervals for hessian evaluation
  xpts::Array{Float64, 2}                 # States at Switching Times
  expMat::Array{Float64, 3}               # Matrix Exponentials
  Phi::Array{Float64, 4}                  # Matrix of State Transition Matrices
  M::Array{Float64, 3}                    # Integrals over Switching Intervals
  S::Array{Float64, 3}                    # S Matrices for each interval
  C::Array{Float64, 3}                     # C Matrices for each interval

  # Store Evaluations
  obj::Array{Float64, 1}                   # Store Objective Value
  deltaval::Array{Float64, 2}              # Store Iterates
  nobjeval::Int                            # Number of objective function evaluations
  ngradeval::Int                           # Number of gradient evaluations
  nhesseval::Int                           # Number of hessian evaluations
end



# Create switching time optimization (STO) abstract type
abstract STO

type linSTO <: STO  # Linear STO type
  model::MathProgBase.AbstractNonlinearModel  # Nonlinear Program Model
  STOev::linSTOev                             # NLP Evaluator for linear STO

  # Data Obtained after Optimization
  delta::Array{Float64,1}                     # Optimal Switching Intervals
  tau::Array{Float64,1}                       # Optimal Switching Times
  objval::Float64                             # Optimal Value of Cost Function
  stat::Symbol                                # Status of Opt Problem
  soltime::Float64                            # Time Required to solve Opt


  # Inner Contructor for Incomplete Initialization
  linSTO(model, STOev, delta) = new(model, STOev, delta)

end

type nlinSTO <: STO  # Nonlinear STO type
  model::MathProgBase.AbstractNonlinearModel  # Nonlinear Program Model
  STOev::nlinSTOev                            # NLP Evaluator for nonlinear STO
  # nartsw::Int64                               # Number of artificial switches


  # Data Obtained After Optimization
  delta::Array{Float64,1}                     # Optimal Switching Intervals
  # deltacomplete::Array{Float64,1}                     # Optimal Switching Intervals (Complete)
  tau::Array{Float64,1}                       # Optimal Switching Times
  # taucomplete::Array{Float64,1}               # Optimal Switching Times (Complete)
  objval::Float64                             # Optimal Value of Cost Function
  stat::Symbol                                # Status of Opt Problem
  soltime::Float64                            # Time Required to solve Opt


  # Inner Contructor for Incomplete Initialization
  nlinSTO(model, STOev, delta) = new(model, STOev, delta)
end
