# TODO: Add all the checks for the varibles sizes!

# Create STO Problem
function createsto(
  x0::Array{Float64,1},                 # Initial State
  A::Array{Float64,3};                  # Linear Dynamics
  t0::Float64=0.0,                      # Initial Time
  tf::Float64=1.0,                      # Final Time
  Q::Array{Float64, 2}=emptyfmat,       # Cost Matrix
  Qf::Array{Float64, 2}=emptyfmat,      # Final Cost Matrix
  lb::Array{Float64, 1}=emptyfvec,      # Lower Bound on Intervals
  ub::Array{Float64, 1}=emptyfvec,      # Upper Bound on Intervals
  tau0ws::Array{Float64,1}=emptyfvec,   # Warm Start tau vector
  solver::MathProgBase.AbstractMathProgSolver=IpoptSolver())

  # Get Dimensions
  nx = size(A, 1)      # State Dimension
  N = size(A, 3) - 1   # Get number of switches

  # Adjust variables which have not been initalized
  if isempty(Q)
    Q = eye(nx)
  end

  if isempty(Qf)
    Qf = zeros(nx, nx)
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


  # # Define Bounds for switching times
  # lbtau = t0*ones(N)
  # ubtau = tf*ones(N)


  # Generate Model
  m = MathProgBase.NonlinearModel(solver)

  ### Initialize NLP Evaluator
  # Preallocate arrays
  prev_delta = Array(Float64, N+1)
  xpts = Array(Float64, nx, N+2)
  expMat = Array(Float64, nx, nx, N+1)
  Phi = Array(Float64, nx, nx, N+2, N+2)
  M = Array(Float64, nx, nx, N+1)
  S = Array(Float64, nx, nx, N+1)
  C = Array(Float64, nx, nx, N+1)


  # Decompose Dynamics Matrices
  V = Array(Complex{Float64}, 2*nx, 2*nx, N+1)
  invV = Array(Complex{Float64}, 2*nx, 2*nx, N+1)
  D =  Array(Complex{Float64}, 2*nx, N+1)
  for i = 1:N+1
  D[:, i], V[:, :, i] = eig([-A[:, :, i]'  Q;
                           zeros(nx, nx) A[:, :, i]])
  invV[:, :, i] = inv(V[:, :, i])
  end



  # Construct Matrix of Indeces for lower triangular Matrix (Hessian)
  IndTril = find(tril(ones(N+1, N+1)))
  Itril, Jtril, _ = findnz(tril(ones(N+1, N+1)))

  # Construct Constraints Matrix (Vector)
  Ag = ones(1, N+1)
  Ig, Jg, Vg = findnz(Ag)
  bg = [tf]   # Only one constraints for the sum of the switching intervals

  # Construct NLPEvaluator
  STOev = linSTOev(x0, nx, A, N, t0, tf, Q, Qf, V, invV, D, IndTril, Jtril, Itril, Ag, Ig, Jg, Vg, bg, prev_delta, xpts, expMat, Phi, M, S, C)

  ### Load NLP Program into the model
  MathProgBase.loadproblem!(m, N+1, 1, lb, ub, bg, bg, :Min, STOev)

  ### Add Warm Starting Point
  MathProgBase.setwarmstart!(m, delta0ws)

  # Create STO
  STOproblem = linSTO(m, STOev, delta0ws)

  return STOproblem  # Return STO


end


# Create STO Problem for Nonlinear STO
function createsto(
  x0::Array{Float64,1},             # Initial State
  nonlin_dyn::Function,             # Nonlinear Dynamics
  nonlin_dyn_deriv::Function,       # Nonlinear Dynamics Derivative
  uvec::Array{Float64, 2},          # Vector of integer Inputs per switching combination
  ngrid::Int64=10;                  # Number of Linearization points in the grid (Apart from sw times)
  t0::Float64=0.0,                  # Initial Time
  tf::Float64=1.0,                  # Final Time
  Q::Array{Float64, 2}=emptyfmat,   # Cost Matrix
  Qf::Array{Float64, 2}=emptyfmat,  # Final Cost Matrix
  lb::Array{Float64, 1}=emptyfvec,  # Lower Bound on Intervals
  ub::Array{Float64, 1}=emptyfvec,  # Upper Bound on Intervals
  tau0ws::Array{Float64,1}=emptyfvec, # Warm Start tau vector
  solver::MathProgBase.AbstractMathProgSolver=IpoptSolver())

  # Get Dimensions
  nx = length(x0)         # State Dimension


  # Generate artificial switches repeating the input sequence
  # uvec = hcat([repmat(uvec[:,i], 1,nartsw+1) for i = 1:size(uvec, 2)]...)
  N = size(uvec, 2) - 1   # Get total number of switches

  # Define Maximum distance between switching times
  # maxSwDist = 2*(tf - t0)/(N+1)
  maxSwDist = Inf

  # Adjust variables which have not been initalized
  if isempty(Q)
    Q = eye(nx)
  end

  if isempty(Qf)
    Qf = zeros(nx, nx)
  end

  if isempty(lb)
    lb = zeros(N+1)
  end

  if isempty(ub)
    ub = maxSwDist*ones(N+1)
  end

  if isempty(tau0ws)
    tau0ws = collect(linspace(t0, tf, N+2))  # Currently counting tau_0 and tau_{N+1}. They are removed below after xpts initialization.
  else
    tau0ws = [t0; tau0ws; tf]
  end


  # Create Discretization grid
  tgrid = collect(linspace(t0, tf, ngrid))

  ### GO ON FROM HERE

  # Generate Model
  m = MathProgBase.NonlinearModel(solver)

  ### Initialize NLP Evaluator
  # Extend System State To get Affine to Linear Dynamics
  nx = nx + 1

  # Extend Initial State and Cost Matrix
  x0 = [x0; 1]
  Q = full(blkdiag(sparse(Q), sparse([0.0])))
  Qf = full(blkdiag(sparse(Qf), sparse([0.0])))


  # Define Required Matrices and Switching Instants
  A = Array(Float64, nx, nx, N+1)

  # Initialize Switching Instants and A matrices by performing initial linearized simulation
  xpts = Array(Float64, nx, N+2)
  xpts[:, 1] = x0

  for i = 2:N+2

    # Generate Linearized Dynamics
    A[:,:,i-1] = linearizeDyn(nonlin_dyn, nonlin_dyn_deriv, xpts[1:end-1,i-1], uvec[:,i-1])

    # Compute Next Point in Simulation from warm starting definition
    xpts[:, i] = expm(A[:, :, i-1]*(tau0ws[i] - tau0ws[i-1]))*xpts[:, i-1]

  end

  # Reduce tau0ws dimension (remove tau_0 and tau_{N+1})
  tau0ws = tau0ws[2:end-1]

  # Define warm starting delta0ws from tau0ws
  delta0ws = tau2delta(tau0ws, t0, tf)


  prev_delta = Array(Float64, N+1)
  expMat = Array(Float64, nx, nx, N+1)
  Phi = Array(Float64, nx, nx, N+2, N+2)
  M = Array(Float64, nx, nx, N+1)
  S = Array(Float64, nx, nx, N+1)
  C = Array(Float64, nx, nx, N+1)


  # Construct Matrix of Indeces for lower triangular Matrix (Hessian)
  IndTril = find(tril(ones(N+1, N+1)))
  Itril, Jtril, _ = findnz(tril(ones(N+1, N+1)))

  # # Construct Constraints Matrix
  # Ag = [-[eye(N-1) zeros(N-1)] + [zeros(N-1) eye(N-1)]; 1 zeros(1,N-1); zeros(1,N-1) -1]
  # Ig, Jg, Vg = findnz(Ag)
  # bg = -[zeros(N-1); -t0; tf]

  # Construct Constraints Matrix (Vector)
  Ag = ones(1, N+1)
  Ig, Jg, Vg = findnz(Ag)
  bg = [tf]   # Only one constraints for the sum of the switching intervals

  # # Debug: Add constraint of equidistant discretization points
  # Norig = Int((N+1)/(nartsw+1)) - 1
  # gMatTemp = -nartsw*eye(nartsw+1) + ones(nartsw+1, nartsw+1) - eye(nartsw+1)
  # gMatTemp = gMatTemp[1:end-1,:]
  # gMatTemp = kron(eye(Norig+1), gMatTemp)
  # Ag = [Ag; gMatTemp]
  # # bg = [bg; zeros(N+1)]
  # bg = [bg; zeros(N+1 - 1*(Norig+1))]
  # Ig, Jg, Vg = findnz(Ag)
  #
  # display(Norig)
  # display(size(Ag))
  # display(size(bg))

  # Construct NLPEvaluator
  STOev = nlinSTOev(x0, nx, A, N, t0, tf, Q, Qf, uvec, nonlin_dyn, nonlin_dyn_deriv, IndTril, Jtril, Itril, Ag, Ig, Jg, Vg, bg, prev_delta, xpts, expMat, Phi, M, S, C)

  ### Load NLP Program into the model
  MathProgBase.loadproblem!(m, N+1, length(bg), lb, ub, bg, bg, :Min, STOev)

  ### Add Warm Starting Point
  MathProgBase.setwarmstart!(m, delta0ws)


  # Create STO
  STOproblem = nlinSTO(m, STOev, nartsw, delta0ws)

  return STOproblem  # Return STO


end


# Solve Optimization for Linear System
function solve!(m::linSTO)

  # Perform STO
  m.soltime = @elapsed MathProgBase.optimize!(m.model)
  m.stat = MathProgBase.status(m.model)
  m.delta = MathProgBase.getsolution(m.model)
  m.tau = delta2tau(m.delta, m.STOev.t0, m.STOev.tf)
  m.objval = MathProgBase.getobjval(m.model)

  # return tauopt, Jopt, solTime, stat
end

# Solve Optimization for Nonlinear System
function solve!(m::nlinSTO)

  # Perform STO
  m.soltime = @elapsed MathProgBase.optimize!(m.model)
  m.stat = MathProgBase.status(m.model)

  # Get delta and tau vectors (Complete and Not Complete)
  m.deltacomplete = MathProgBase.getsolution(m.model)
  m.taucomplete = delta2tau(m.deltacomplete, m.STOev.t0, m.STOev.tf)
  m.tau = m.taucomplete[m.nartsw+1:m.nartsw+1:end]  # N.B. nrep and not nrep+1 to start because uvec has one more element than tauopt
  m.delta = tau2delta(m.tau, m.STOev.t0, m.STOev.tf)

  m.objval = MathProgBase.getobjval(m.model)

end



# Return Variables from STO
gettau(m::STO) = m.tau
getdelta(m::STO) = m.delta
gettaucomplete(m::nlinSTO) = m.taucomplete
getdeltacomplete(m::nlinSTO) = m.deltacomplete
getobjval(m::STO) = m.objval
getstat(m::STO) = m.stat
getsoltime(m::STO) = m.soltime


# Convert from Switching Times to Intervals
function tau2delta(tau::Array{Float64, 1}, t0::Float64, tf::Float64)
  # Extend vector of switching times with initial and final time
  tauexp = [t0; tau; tf]  # Extend tau vector to simplify numbering

  # Return Vector of differences
  return  diff(tauexp)

end


# Convert from Intervals to Switching Times
function delta2tau(delta::Array{Float64, 1}, t0::Float64, tf::Float64)


  # Define tau vector
  tau = Array(Float64, length(delta)-1)

  # Initialize first tau element
  tau[1] = t0 + delta[1]

  # Iteratively add all deltas to create switching times vector
  for i = 2:length(tau)
    tau[i] = tau[i-1] + delta[i]
  end

  return tau

end
