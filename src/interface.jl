# TODO: Add all the checks for the varibles sizes!

# Create STO Problem
function createsto(x0::Array{Float64,1}, A::Array{Float64,3};
                            t0::Float64=0.0,  # Initial Time
                            tf::Float64=1.0,  # Final Time
                            Q::Array{Float64, 2}=emptyfmat,   # Cost Matrix
                            gl::Array{Float64, 1}=emptyfvec,  # Lower Bound on Intervals
                            gu::Array{Float64, 1}=emptyfvec,  # Upper Bound on Intervals
                            x0ws::Array{Float64,1}=emptyfvec, # Warm Start x0
                            solver::MathProgBase.AbstractMathProgSolver=IpoptSolver())

  # Get Dimensions
  nx = size(A, 1)      # State Dimension
  N = size(A, 3) - 1   # Get number of switches

  # Adjust variables which have not been initalized
  if isempty(Q)
    Q = eye(nx)
  end

  if isempty(gl)
    gl = zeros(N+1)
  end

  if isempty(gu)
    gu = Inf*ones(N+1)
  end

  if isempty(x0ws)
    x0ws = linspace(t0, tf, N+2)
    x0ws = x0ws[2:end-1]
  end


  # Define Bounds for switching times
  lb = t0*ones(N)
  ub = tf*ones(N)

  # Generate Model
  m = MathProgBase.NonlinearModel(solver)

  ### Initialize NLP Evaluator
  # Preallocate arrays
  prev_tau = Array(Float64, N)
  xpts = Array(Float64, nx, N+1)
  expMat = Array(Float64, nx, nx, N+1)
  Phi = Array(Float64, nx, nx, N+1, N+1)
  M = Array(Float64, nx, nx, N+1)
  P = Array(Float64, nx, nx, N+1)

  # Decompose Dynamics Matrices
  # V = Array(Float64, 2*nx, 2*nx, N+1)
  # invV = Array(Float64, 2*nx, 2*nx, N+1)
  # D =  Array(Float64, 2*nx, N+1)
  V = Array(Complex{Float64}, 2*nx, 2*nx, N+1)
  invV = Array(Complex{Float64}, 2*nx, 2*nx, N+1)
  D =  Array(Complex{Float64}, 2*nx, N+1)
  for i = 1:N+1
  D[:, i], V[:, :, i] = eig([-A[:, :, i]'  Q;
                           zeros(nx, nx) A[:, :, i]])
  invV[:, :, i] = inv(V[:, :, i])
  end



  # Construct Matrix of Indeces for lower triangular Matrix
  IndTril = find(tril(ones(N, N)))
  Itril, Jtril, _ = findnz(tril(ones(N, N)))

  # Construct Constraints Matrix
  Ag = [-[eye(N-1) zeros(N-1)] + [zeros(N-1) eye(N-1)]; 1 zeros(1,N-1); zeros(1,N-1) -1]
  Ig, Jg, Vg = findnz(Ag)
  bg = -[zeros(N-1); -t0; tf]

  # Construct NLPEvaluator
  STOev = linSTOev(x0, nx, A, N, t0, tf, Q, V, invV, D, IndTril, Jtril, Itril, Ag, Ig, Jg, Vg, bg, prev_tau, xpts, expMat, Phi, M, P)

  # Old Call
  # nlpSTO = STOev(x0, nx, A, N, t0, tf, Q)


  ### Load NLP Program into the model
  MathProgBase.loadproblem!(m, N, N+1, lb, ub, gl, gu, :Min, STOev)

  ### Add Warm Starting Point
  MathProgBase.setwarmstart!(m, x0ws)

  # Create STO
  STOproblem = linSTO(m, STOev)

  return STOproblem  # Return STO


end


# Create STO Problem for Nonlinear STO
function createsto(x0::Array{Float64,1},         # Initial State
                            nonlin_dyn::Function,             # Nonlinear Dynamics
                            nonlin_dyn_deriv::Function,       # Nonlinear Dynamics Derivative
                            uvec::Array{Float64, 2};           # Vector of integer Inputs per switching combination
                            t0::Float64=0.0,  # Initial Time
                            tf::Float64=1.0,  # Final Time
                            Q::Array{Float64, 2}=emptyfmat,   # Cost Matrix
                            gl::Array{Float64, 1}=emptyfvec,  # Lower Bound on Intervals
                            gu::Array{Float64, 1}=emptyfvec,  # Upper Bound on Intervals
                            x0ws::Array{Float64,1}=emptyfvec, # Warm Start x0
                            nartsw::Int64=6,                  # Number of Artificial Switches added at each switching time
                            solver::MathProgBase.AbstractMathProgSolver=IpoptSolver())

  # Get Dimensions
  nx = length(x0)         # State Dimension


  # Generate artificial switches repeating the input sequence
  uvec = hcat([repmat(uvec[:,i], 1,nartsw) for i = 1:size(uvec, 2)]...)
  N = size(uvec, 2) - 1   # Get total number of switches

  # Define Maximum distance between switching times
  maxSwDist = 2*(tf - t0)/(N+1)

  # Adjust variables which have not been initalized
  if isempty(Q)
    Q = eye(nx)
  end

  if isempty(gl)
    gl = zeros(N+1)
  end

  if isempty(gu)
    gu = maxSwDist*ones(N+1)
  end

  if isempty(x0ws)
    x0ws = linspace(t0, tf, N+2)
    x0ws = x0ws[2:end-1]
  end



  # Define Bounds for switching times
  lb = t0*ones(N)
  ub = tf*ones(N)

  # Generate Model
  m = MathProgBase.NonlinearModel(solver)

  ### Initialize NLP Evaluator
  # Extend System State To get Affine to Linear Dynamics
  nx = nx + 1

  # Extend Initial State and Cost Matrix
  x0 = [x0; 1]
  Q = full(blkdiag(sparse(Q), sparse([0.0])))


  # Define Required Matrices
  A = Array(Float64, nx, nx, N+1)

  prev_tau = Array(Float64, N)
  # xpts = Array(Float64, nx, N+1)
  xpts = repmat(x0, 1, N+1)        # Initialize states at switching instants
  expMat = Array(Float64, nx, nx, N+1)
  Phi = Array(Float64, nx, nx, N+1, N+1)
  M = Array(Float64, nx, nx, N+1)
  P = Array(Float64, nx, nx, N+1)


  # Construct Matrix of Indeces for lower triangular Matrix
  IndTril = find(tril(ones(N, N)))
  Itril, Jtril, _ = findnz(tril(ones(N, N)))

  # Construct Constraints Matrix
  Ag = [-[eye(N-1) zeros(N-1)] + [zeros(N-1) eye(N-1)]; 1 zeros(1,N-1); zeros(1,N-1) -1]
  Ig, Jg, Vg = findnz(Ag)
  bg = -[zeros(N-1); -t0; tf]

  # Construct NLPEvaluator
  STOev = nlinSTOev(x0, nx, A, N, t0, tf, Q, uvec, nonlin_dyn, nonlin_dyn_deriv, IndTril, Jtril, Itril, Ag, Ig, Jg, Vg, bg, prev_tau, xpts, expMat, Phi, M, P)

  # Old Stuff
  # nlpSTO = STOev(x0, nx, N, t0, tf, Q, uvec, nartsw, nldyn, nldyn_deriv)


  ### Load NLP Program into the model
  MathProgBase.loadproblem!(m, N, N+1, lb, ub, gl, gu, :Min, STOev)

  ### Add Warm Starting Point
  MathProgBase.setwarmstart!(m, x0ws)


  # Create STO
  STOproblem = nlinSTO(m, STOev, nartsw)

  return STOproblem  # Return STO


end


# Solve Optimization for Linear System
function solve!(m::linSTO)

  # Perform STO
  m.soltime = @elapsed MathProgBase.optimize!(m.model)
  m.stat = MathProgBase.status(m.model)
  m.tau = MathProgBase.getsolution(m.model)
  m.objval = MathProgBase.getobjval(m.model)

  # return tauopt, Jopt, solTime, stat
end

# Solve Optimization for Nonlinear System
function solve!(m::nlinSTO)

  # Perform STO
  m.soltime = @elapsed MathProgBase.optimize!(m.model)
  m.stat = MathProgBase.status(m.model)
  m.taucomplete = MathProgBase.getsolution(m.model)
  m.tau = m.taucomplete[m.nartsw:m.nartsw:end]  # N.B. nrep and not nrep+1 to start because uvec has one more element than tauopt
  m.objval = MathProgBase.getobjval(m.model)

  # return tauopt, Jopt, solTime, stat
end

# # Solve Optimization for Nonlinear System
# function solveSTO(m::MathProgBase.AbstractNonlinearModel)
#
#   # Perform STO
#   solTime = @elapsed MathProgBase.optimize!(m)
#   stat = MathProgBase.status(m)
#   tauopt = MathProgBase.getsolution(m)
#   Jopt = MathProgBase.getobjval(m)
#
#   return tauopt, Jopt, solTime, stat
#
# end


# Return Variables from STO
gettau(m::STO) = m.tau
gettaucomplete(m::nlinSTO) = m.taucomplete
getobjval(m::STO) = m.objval
getstat(m::STO) = m.stat
getsoltime(m::STO) = m.soltime
