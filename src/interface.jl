# Create STO Problem
function stoproblem(
  x0::Array{Float64,1},                 # Initial State
  A::Array{Float64,3};                  # Linear Dynamics
  ngrid::Int64=2,                       # Number of Linearization points in the fixed grid (2 for linear case. We do not need them by default)
  t0::Float64=0.0,                      # Initial Time
  tf::Float64=1.0,                      # Final Time
  Q::Array{Float64, 2}=emptyfmat,       # Cost Matrix
  E::Array{Float64, 2}=emptyfmat,      # Final Cost Matrix
  lb::Array{Float64, 1}=emptyfvec,      # Lower Bound on Intervals
  ub::Array{Float64, 1}=emptyfvec,      # Upper Bound on Intervals
  tau0ws::Array{Float64,1}=emptyfvec,   # Warm Start tau vector
  solver::MathProgBase.AbstractMathProgSolver=Ipopt.IpoptSolver(maxiter = 30))

  # Get Dimensions
  nx = size(A, 1)      # State Dimension
  N = size(A, 3) - 1   # Get number of switches

  # Adjust variables which have not been initalized
  if isempty(Q)
    Q = eye(nx)
  end

  if isempty(E)
    E = zeros(nx, nx)
  end

  if isempty(lb)
    lb = zeros(N+1)
  end

  if isempty(ub)
    ub = Inf*ones(N+1)
  end

  if isempty(tau0ws)
    tau0ws = collect(linspace(t0, tf, N+2))
    tau0ws = tau0ws[2:end-1]
  end

  # Define warm starting delta0
  delta0ws = tau2delta(tau0ws, t0, tf)

  # Create Discretization Grid
  tgrid = collect(linspace(t0, tf, ngrid))

  # Initialize time vectors
  tvec = Array(Float64, N + ngrid)    # Complete grid
  tauIdx = Array(Int, N + 2)      # Indeces of switching times in the complete grid
  tgridIdx = Array(Int, ngrid)      # Indeces of switching times in the complete grid
  deltacomplete = Array(Float64, N + ngrid - 1)   # Intervals over the whole grid


  ### Initialize NLP Evaluator
  # Preallocate arrays
  deltafun_prev = Array(Float64, N+1)
  deltagrad_prev = Array(Float64, N+1)
  deltahess_prev = Array(Float64, N+1)
  xpts = Array(Float64, nx, N+ngrid); xpts[:, 1] = x0   # Set Initial State
  expMat = Array(Float64, nx, nx, N+ngrid-1)
  Phi = Array(Float64, nx, nx, N+2, N+2)
  M = Array(Float64, nx, nx, N+ngrid-1)
  S = Array(Float64, nx, nx, N+ngrid)
  C = Array(Float64, nx, nx, N+1)


  # Decompose Dynamics Matrices
  V = Array(Complex{Float64}, 2*nx, 2*nx, N+1)
  invV = Array(Complex{Float64}, 2*nx, 2*nx, N+1)
  D =  Array(Complex{Float64}, 2*nx, N+1)
  isDiag = Array(Bool, N+1)

  for i = 1:N+1
    D[:, i], V[:, :, i] = eig([-A[:, :, i]'  Q;
                             zeros(nx, nx) A[:, :, i]])
    if cond(V[:, :, i]) == Inf  # Non diagonalizable matrix
      isDiag[i] = false
    else
      invV[:, :, i] = inv(V[:, :, i])
      isDiag[i] = true
    end
  end


  # Construct Matrix of Indeces for lower triangular Matrix (Hessian)
  IndTril = find(tril(ones(N+1, N+1)))
  Itril, Jtril, _ = findnz(tril(ones(N+1, N+1)))

  # Construct Constraints Matrix (Vector)
  Ag = ones(1, N+1)
  Ig, Jg, Vg = findnz(Ag)
  bg = [tf]   # Only one constraints for the sum of the switching intervals


  # Initialize objective evaluator
  obj = Array(Float64, 0)
  deltaval = Array(Float64, N+1, 0)
  nobjeval = 0                           # Number of objective function evaluations
  ngradeval = 0                           # Number of gradient evaluations
  nhesseval = 0                           # Number of hessian evaluations

  # Construct NLP evaluator
  STOev = linSTOev(x0, nx, A, N, t0, tf, tf, Q, E, ngrid, tgrid, tvec, tauIdx, tgridIdx, deltacomplete, V, invV, D, isDiag, IndTril, Itril, Jtril, Ag, Ig, Jg, Vg, bg, deltafun_prev, deltagrad_prev, deltahess_prev, xpts, expMat, Phi, M, S, C,
  obj, deltaval, nobjeval, ngradeval, nhesseval)


  # Generate Model
  m = MathProgBase.NonlinearModel(solver)


  ### Load NLP Program into the model
  MathProgBase.loadproblem!(m, N+1, length(bg), lb, ub, bg, bg, :Min, STOev)


  # Propagate dynamic for new switching times
  propagateDynamics!(STOev, delta0ws)

  ### Add Warm Starting Point
  MathProgBase.setwarmstart!(m, delta0ws)

  # Create STO
  STOproblem = linSTO(m, STOev, delta0ws)

  return STOproblem  # Return STO


end


# Create STO Problem for Nonlinear STO
function stoproblem(
  x0::Array{Float64,1},             # Initial State
  nonlin_dyn::Function,             # Nonlinear Dynamics
  nonlin_dyn_deriv::Function,       # Nonlinear Dynamics Derivative
  uvec::Array{Float64, 2};          # Vector of integer Inputs per switching combination
  ngrid::Int64=10,                  # Number of Linearization points in the fixed grid
  t0::Float64=0.0,                  # Initial Time
  tf::Float64=1.0,                  # Final Time
  Q::Array{Float64, 2}=emptyfmat,   # Cost Matrix
  E::Array{Float64, 2}=emptyfmat,  # Final Cost Matrix
  lb::Array{Float64, 1}=emptyfvec,  # Lower Bound on Intervals
  ub::Array{Float64, 1}=emptyfvec,  # Upper Bound on Intervals
  tau0ws::Array{Float64,1}=emptyfvec, # Warm Start tau vector
  solver::MathProgBase.AbstractMathProgSolver=Ipopt.IpoptSolver(maxiter = 30))

  # Get Dimensions
  nx = length(x0)         # State Dimension


  # Generate artificial switches repeating the input sequence
  N = size(uvec, 2) - 1   # Get total number of switches

  # Adjust variables which have not been initalized
  if isempty(Q)
    Q = eye(nx)
  end

  if isempty(E)
    E = zeros(nx, nx)
  end

  if isempty(lb)
    lb = zeros(N+1)
  end

  if isempty(ub)
    ub = Inf*ones(N+1)
  end

  if isempty(tau0ws)
    tau0ws = collect(linspace(t0, tf, N+2))  # Currently counting tau_0 and tau_{N+1}. They are removed below after xpts initialization.
    tau0ws = tau0ws[2:end-1]  # Include only switching instants

  end



  # Create Discretization grid
  tgrid = collect(linspace(t0, tf, ngrid))


  # Initialize time vectors
  tvec = Array(Float64, N + ngrid)    # Complete grid
  tauIdx = Array(Int, N + 2)      # Indeces of switching times in the complete grid
  tgridIdx = Array(Int, ngrid)      # Indeces of switching times in the complete grid
  deltacomplete = Array(Float64, N + ngrid - 1)   # Intervals over the whole grid

  ### Initialize NLP Evaluator
  # Extend System State To get Affine to Linear Dynamics
  nx = nx + 1

  # Extend Initial State and Cost Matrix
  x0 = [x0; 1]
  spz = Array(Float64, 1,1); spz[1,1] = 0.0; spz = sparse(spz)  # Sparse scalar
  Q = full(blkdiag(sparse(Q), spz))
  E = full(blkdiag(sparse(E), spz))


  # Define Required Matrices and Switching Instants
  A = Array(Float64, nx, nx, N+ngrid-1)

  # Initialize Switching Instants
  xpts = Array(Float64, nx, N+ngrid)
  xpts[:, 1] = x0



  # Define warm starting delta0ws from tau0ws
  delta0ws = tau2delta(tau0ws, t0, tf)


  deltafun_prev = Array(Float64, N+1)
  deltagrad_prev = Array(Float64, N+1)
  deltahess_prev = Array(Float64, N+1)
  expMat = Array(Float64, nx, nx, N+ngrid-1)
  Phi = Array(Float64, nx, nx, N+2, N+2)
  M = Array(Float64, nx, nx, N+ngrid-1)
  S = Array(Float64, nx, nx, N+ngrid)
  C = Array(Float64, nx, nx, N+1)  # Only at sw times (including initial time)


  # Construct Matrix of Indeces for lower triangular Matrix (Hessian)
  IndTril = find(tril(ones(N+1, N+1)))
  Itril, Jtril, _ = findnz(tril(ones(N+1, N+1)))


  # Construct Constraints Matrix (Vector)
  Ag = ones(1, N+1)
  Ig, Jg, Vg = findnz(Ag)
  bg = [tf]   # Only one constraints for the sum of the switching intervals

  # Initialize objective evaluator
  obj = Array(Float64, 0)
  deltaval = Array(Float64, N+1, 0)
  nobjeval = 0                           # Number of objective function evaluations
  ngradeval = 0                           # Number of gradient evaluations
  nhesseval = 0                           # Number of hessian evaluations

  # Construct NLPEvaluator
  STOev = nlinSTOev(x0, nx, A, N, t0, tf, tf, Q, E, uvec, ngrid, tgrid, tvec, tauIdx, tgridIdx, deltacomplete, nonlin_dyn, nonlin_dyn_deriv, IndTril, Itril, Jtril, Ag, Ig, Jg, Vg, bg, deltafun_prev, deltagrad_prev, deltahess_prev, xpts, expMat, Phi, M, S, C, obj, deltaval, nobjeval, ngradeval, nhesseval)


  # Propagate Dynamics to compute matrix exponentials and states at the switching times
  propagateDynamics!(STOev, delta0ws)

  # Generate Model
  m = MathProgBase.NonlinearModel(solver)

  ### Load NLP Program into the model
  MathProgBase.loadproblem!(m, N+1, length(bg), lb, ub, bg, bg, :Min, STOev)

  ### Add Warm Starting Point for the solver
  MathProgBase.setwarmstart!(m, delta0ws)


  # Create STO
  STOproblem = nlinSTO(m, STOev, delta0ws)

  return STOproblem  # Return STO


end


# Solve Optimization for Linear System
function solve!(m::STO)

  # Perform STO
  m.soltime = @elapsed MathProgBase.optimize!(m.model)
  m.stat = MathProgBase.status(m.model)
  m.delta = MathProgBase.getsolution(m.model)
  m.tau, _ = delta2tau(m.delta, m.STOev.t0)
  m.objval = MathProgBase.getobjval(m.model)

end

"Set warm starting point"
function setwarmstart!(m::STO, tau0ws::Array{Float64,1})

  # Define warm starting delta0
  delta0ws = tau2delta(tau0ws, m.STOev.t0, m.STOev.tf)

  # Set Warm Starting Point for Nonlinear Solver
  MathProgBase.setwarmstart!(m.model, delta0ws)

  # Reset Cost Function Iterates
  m.STOev.obj = Array(Float64,0)
  m.STOev.deltaval = Array(Float64, m.STOev.N+1, 0)

end


"Set initial state x0"
function setx0!(m::STO, x0::Array{Float64,1})

  # Define initial point
  m.STOev.x0 = x0

end


# Return Variables from STO
gettau(m::STO) = m.tau
getdelta(m::STO) = m.delta
getobjval(m::STO) = m.objval
getstat(m::STO) = m.stat
getsoltime(m::STO) = m.soltime
getnobjeval(m::STO) = m.STOev.nobjeval
getngradeval(m::STO) = m.STOev.ngradeval
getnhesseval(m::STO) = m.STOev.nhesseval


# Convert from Switching Times to Intervals
function tau2delta(tau::Array{Float64, 1}, t0::Float64, tf::Float64)

  # Extend vector of switching times with initial and final time
  tauexp = [t0; tau; tf]  # Extend tau vector to simplify numbering

  # Return Vector of differences
  return  diff(tauexp)

end


# Convert from Intervals to Switching Times
function delta2tau(delta::Array{Float64, 1}, t0::Float64)

  # Define tau vector
  tau = Array(Float64, length(delta)-1)

  # Initialize first tau element
  tau[1] = t0 + delta[1]

  # Iteratively add all deltas to create switching times vector
  for i = 2:length(tau)
    tau[i] = tau[i-1] + delta[i]
  end

  tfdelta = tau[end] + delta[end]

  return tau, tfdelta

end
