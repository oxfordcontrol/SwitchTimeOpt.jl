# Empty matrix and vector for function evaluation
const emptyfvec = Array(Float64, 0)
const emptyfmat = Array(Float64, 0, 0)


# Define Abstract NLP Evaluator for STO problem
abstract STOev <: MathProgBase.AbstractNLPEvaluator


# Linear Case
type linSTOev <: STOev
  # Parameters
  x0::Array{Float64,1}          # Initial State x0
  nx::Int                       # State dimension
  A::Array{Float64,3}           # Dynamics Matrices
  N::Int                        # Number of Switching Times
  t0::Float64                   # Initial Time
  tf::Float64                   # Final Time
  Q::Array{Float64,2}           # State Cost

  # Precomputed Values
  # V::Array{Float64, 3}            # Dynamics Matrices decomp V
  # invV::Array{Float64, 3}         # Dynamics Matrices decomp V^-1
  # D::Array{Float64, 2}            # Dynamics Matrices decomp D
  V::Array{Complex{Float64}, 3}            # Dynamics Matrices decomp V
  invV::Array{Complex{Float64}, 3}         # Dynamics Matrices decomp V^-1
  D::Array{Complex{Float64}, 2}            # Dynamics Matrices decomp D
  IndTril::Array{Int, 1}         # Indeces of Lower Triangular Matrices
  Itril::Array{Int, 1}           # Indeces of Lower Triangular Matrices
  Jtril::Array{Int, 1}           # Indeces of Lower Triangular Matrices
  Ag::Array{Float64, 2}          # Linear Constraints Matrix
  Ig::Array{Int, 1}              # Index Linear Constraint Matrix
  Jg::Array{Int, 1}              # Index Linear Constraint Matrix
  Vg::Array{Float64, 1}          # Index Linear Constraint Matrix vectors
  bg::Array{Float64, 1}          # Linear Constraint constant


  # Shared Data Between Functions
  prev_tau::Array{Float64, 1}    # Previous Switching Times Vector
  xpts::Array{Float64, 2}        # States at Switching Times
  expMat::Array{Float64, 3}      # Matrix Exponentials
  Phi::Array{Float64, 4}         # Matrix of State Transition Matrices
  M::Array{Float64, 3}           # Integrals over Switching Intervals
  P::Array{Float64, 3}           # P Matrices for each switching instant
end

# Nonlinear Case
type nlinSTOev <: STOev
  # Parameters
  x0::Array{Float64,1}          # Initial State x0
  nx::Int                       # State dimension
  A::Array{Float64,3}           # Linearized dynamics Matrices
  N::Int                        # Number of Switching Times
  t0::Float64                   # Initial Time
  tf::Float64                   # Final Time
  Q::Array{Float64,2}           # State Cost
  uvec::Array{Float64,2}        # Array of switching inputs
  # nartsw::Int64                 # Number of artificial switches

  # Nonlinear Dynamics and Derivatives Functions
  nonlin_dyn::Function
  nonlin_dyn_deriv::Function

  # Precomputed Values
  IndTril::Array{Int, 1}         # Indeces of Lower Triangular Matrices
  Itril::Array{Int, 1}           # Indeces of Lower Triangular Matrices
  Jtril::Array{Int, 1}           # Indeces of Lower Triangular Matrices
  Ag::Array{Float64, 2}          # Linear Constraints Matrix
  Ig::Array{Int, 1}              # Index Linear Constraint Matrix
  Jg::Array{Int, 1}              # Index Linear Constraint Matrix
  Vg::Array{Float64, 1}          # Index Linear Constraint Matrix vectors
  bg::Array{Float64, 1}          # Linear Constraint constant

  # Shared Data Between Functions
  prev_tau::Array{Float64, 1}    # Previous Switching Times Vector
  xpts::Array{Float64, 2}        # States at Switching Times
  expMat::Array{Float64, 3}      # Matrix Exponentials
  Phi::Array{Float64, 4}         # Matrix of State Transition Matrices
  M::Array{Float64, 3}           # Integrals over Switching Intervals
  P::Array{Float64, 3}           # P Matrices for each switching instant
end



# Create switching time optimization (STO) abstract type
abstract STO

type linSTO <: STO  # Linear STO type
  model::MathProgBase.AbstractNonlinearModel  # Nonlinear Program Model
  STOev::linSTOev                             # NLP Evaluator for linear STO

  # Data Obtained after Optimization
  tau::Array{Float64,1}                       # Optimal Switching Times
  objval::Float64                             # Optimal Value of Cost Function
  stat::Symbol                                # Status of Opt Problem
  soltime::Float64                            # Time Required to solve Opt


  # Inner Contructor for Incomplete Initialization
  linSTO(model, STOev) = new(model, STOev)
end

type nlinSTO <: STO  # Nonlinear STO type
  model::MathProgBase.AbstractNonlinearModel  # Nonlinear Program Model
  STOev::nlinSTOev                             # NLP Evaluator for nonlinear STO
  nartsw::Int64                               # Number of artificial switches


  # Data Obtained After Optimizatin
  tau::Array{Float64,1}                       # Optimal Switching Times
  taucomplete::Array{Float64,1}               # Optimal Switching Times
  objval::Float64                             # Optimal Value of Cost Function
  stat::Symbol                                # Status of Opt Problem
  soltime::Float64                            # Time Required to solve Opt


  # Inner Contructor for Incomplete Initialization
  nlinSTO(model, STOev, nartsw) = new(model, STOev, nartsw)
end
