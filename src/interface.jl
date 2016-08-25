# TODO: Add all the checks for the varibles sizes!

# Create STO Problem
function stoproblem(
  x0::Array{Float64,1},                 # Initial State
  A::Array{Float64,3};                  # Linear Dynamics
  ngrid::Int64=2,                       # Number of Linearization points in the fixed grid (2 for linear case. We do not need them by default)
  t0::Float64=0.0,                      # Initial Time
  tf::Float64=1.0,                      # Final Time
  Q::Array{Float64, 2}=emptyfmat,       # Cost Matrix
  Qf::Array{Float64, 2}=emptyfmat,      # Final Cost Matrix
  Ac::Array{Float64, 2}=emptyfmat,       # Stage Constraint Matrix
  bc::Array{Float64, 1}=emptyfvec,       # Stage Constraint Vector
  Acf::Array{Float64, 2}=emptyfmat,       # Final Stage Constraint Matrix
  bcf::Array{Float64, 1}=emptyfvec,       # Final Stage Constraint Vector
  lb::Array{Float64, 1}=emptyfvec,      # Lower Bound on Intervals
  ub::Array{Float64, 1}=emptyfvec,      # Upper Bound on Intervals
  tau0ws::Array{Float64,1}=emptyfvec,   # Warm Start tau vector
  solver::MathProgBase.AbstractMathProgSolver=Ipopt.IpoptSolver())

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





  # Create Discretization Grid
  tgrid = collect(linspace(t0, tf, ngrid))

  # # Create merged and sorted time vector with grid and switching times
  # tvec = sort(vcat(tgrid, tau0ws))
  #
  # # Create index of the tau vector elements inside tvec
  # tauIdx = Array(Int, N+2); tauIdx[1] = 1
  # for i = 1:N
  #   tauIdx[i+1] = findfirst(tvec, tau0ws[i])  # i+1 because tau0ws
  # end
  # tauIdx[end] = N + ngrid


  # Initialize time vectors
  tvec = Array(Float64, N + ngrid)    # Complete grid
  tauIdx = Array(Int, N + 2)      # Indeces of switching times in the complete grid
  tgridIdx = Array(Int, ngrid)      # Indeces of switching times in the complete grid
  deltacomplete = Array(Float64, N + ngrid - 1)   # Intervals over the whole grid


  # tvec, tauIdx = mergeSortFindIndex(tgrid, tau0ws)
  #
  # # Get complete delta vector with all intervals
  # deltacomplete = tau2delta(tvec[2:end-1], t0, tf)




  # # Define Bounds for switching times
  # lbtau = t0*ones(N)
  # ubtau = tf*ones(N)


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

  # # Construct Constraints Matrix (Vector)
  # Ag = ones(1, N+1)
  # Ig, Jg, Vg = findnz(Ag)
  # bg = [tf]   # Only one constraints for the sum of the switching intervals

  #-----------------------------------------------------------------------------
  # Constraints FIXME
  #-----------------------------------------------------------------------------

  # First line is sum constraint
  gsum = ones(1, N+1)
  bgsuml = tf; bgsumu = tf;

  # Final Stage Constraints
  if isempty(Acf)  # No last stage constraints
    if !isempty(bcf) # If provided only bcf and not Acf, throw error!
      error("Only a linear constraints vector, but no matrix has been specified!")
    end
    nconsf = 0
    # Acf = Array(Float64, 0, 0)
    bgfl = Array(Float64, 0)
    bgfu = Array(Float64, 0)

  else
    nconsf = size(Acf, 1)             # Number of constraints at last stage
    bgfu = bcf                        # Upper Bound
    bgfl = -Inf*ones(length(bcf))     # Lower bound
  end

  # Stagewise Constraints
  if isempty(Ac)  # No last stage constraints
    if !isempty(bc) # If provided only bc and not Ac, throw error!
      error("Only a linear constraints vector, but no matrix has been specified!")
    end
    ncons = 0
    bgcl = Array(Float64, 0)
    bgcu = Array(Float64, 0)
  else
    ncons = size(Ac, 1)              # Number of consraints per stage
    bgcu = repmat(bc, ngrid - 2)      # Upper bound for all the stages
    bgcl = -Inf*ones(length(bc)*(ngrid-2))     # Lower bound for all the stages
  end

  # Compute index matrices
  Ig, Jg, Vg = findnz(ones(1 + ncons*(ngrid-2) + nconsf, N+1))  # Compute indexed version of jacobian to define vectors (need to define it at every iteration stage)
  bgu = [tf; bgcu; bgfu]        # Constraints lower bound
  bgl = [tf; bgcl; bgfl]        # Constraints upper bound



  # show(bgl)
  # show(bgu)
  # # Old Constraints
  # # Final Stage Constraints
  # if isempty(Acf)  # No Constraints
  #   if !isempty(bcf) # If provided only bc and not Ac, throw error!
  #     error("Only a linear constraints vector, but no matrix has been specified!")
  #   end
  #
  #   ncons = 0
  #   Acf = Array(Float64, 0, 0)
  #   Ig, Jg, Vg = findnz(gsum)
  #   bgu = [tf]; bgl = [tf]
  #
  # else
  #   # Then Stage constraint (Start by last stage)
  #   ncons = size(Acf, 1)             # Number of Constraints per stage
  #   Ig, Jg, Vg = findnz(ones(ncons+1, N+1))         # Compute indexed version of jacobian to define vectors (need to define it at every iteration stage)
  #   bgu = [tf; bcf]                  # Constraints lower bound
  #   bgl = [tf; -Inf*ones(length(bcf))]   # Constraints upper bound
  # end



  # Initialize objective evaluator
  obj = Array(Float64, 0)
  deltaval = Array(Float64, N+1, 0)
  nobjeval = 0                           # Number of objective function evaluations
  ngradeval = 0                           # Number of gradient evaluations
  nhesseval = 0                           # Number of hessian evaluations

  #-----------------------------------------------------------------------------
  # Construct NLP evaluator
  #-----------------------------------------------------------------------------
  STOev = linSTOev(x0, nx, A, N, t0, tf, tf, Q, Qf, ngrid, tgrid, tvec, tauIdx, tgridIdx, deltacomplete, ncons, nconsf, V, invV, D, isDiag, IndTril, Itril, Jtril, Ac, Acf, gsum, Ig, Jg, Vg, deltafun_prev, deltagrad_prev, deltahess_prev, xpts, expMat, Phi, M, S, C,
  obj, deltaval, nobjeval, ngradeval, nhesseval)


  # Generate Model
  m = MathProgBase.NonlinearModel(solver)


  ### Load NLP Program into the model
  MathProgBase.loadproblem!(m, N+1, length(bgu), lb, ub, bgl, bgu, :Min, STOev)


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
  Qf::Array{Float64, 2}=emptyfmat,  # Final Cost Matrix
  lb::Array{Float64, 1}=emptyfvec,  # Lower Bound on Intervals
  ub::Array{Float64, 1}=emptyfvec,  # Upper Bound on Intervals
  tau0ws::Array{Float64,1}=emptyfvec, # Warm Start tau vector
  solver::MathProgBase.AbstractMathProgSolver=Ipopt.IpoptSolver())

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
    tau0ws = tau0ws[2:end-1]  # Include only switching instants
  # else
  #   tau0ws = [t0; tau0ws; tf]
  end

  # tau0full = [t0; tau0ws; tf]  # tau vector including initial and final time


  # Create Discretization grid
  tgrid = collect(linspace(t0, tf, ngrid))

  # # Create merged and sorted time vector with grid and switching times
  # tvec = sort(vcat(tgrid, tau0ws))
  #
  # # Create index of the tau vector elements inside tvec
  # tauIdx = Array(Int, N+2)
  # tauIdx[1] = 1
  #
  # for i = 1:N
  #   tauIdx[i+1] = findfirst(tvec, tau0ws[i])  # i+1 because tau0ws
  # end
  # tauIdx[end] = N + ngrid

  # Initialize time vectors
  tvec = Array(Float64, N + ngrid)    # Complete grid
  tauIdx = Array(Int, N + 2)      # Indeces of switching times in the complete grid
  tgridIdx = Array(Int, ngrid)      # Indeces of switching times in the complete grid
  deltacomplete = Array(Float64, N + ngrid - 1)   # Intervals over the whole grid

  # tvec, tauIdx = mergeSortFindIndex(tgrid, tau0ws)
  #
  # # Get complete delta vector with all intervals
  # deltacomplete = tau2delta(tvec[2:end-1], t0, tf)



  ### Initialize NLP Evaluator
  # Extend System State To get Affine to Linear Dynamics
  nx = nx + 1

  # Extend Initial State and Cost Matrix
  x0 = [x0; 1]
  spz = Array(Float64, 1,1); spz[1,1] = 0.0; spz = sparse(spz)  # Sparse scalar
  Q = full(blkdiag(sparse(Q), spz))
  Qf = full(blkdiag(sparse(Qf), spz))
  # Q = full(blkdiag(sparse(Q), sparse([0.0])))
  # Qf = full(blkdiag(sparse(Qf), sparse([0.0])))


  # Define Required Matrices and Switching Instants
  A = Array(Float64, nx, nx, N+ngrid-1)

  # Initialize Switching Instants
  xpts = Array(Float64, nx, N+ngrid)
  xpts[:, 1] = x0



  # uIdx = 1  # Initialize index for current u
  #
  # for i = 1:N+ngrid-1
  #
  #   # Verify which U input applies
  #   if uIdx<=N
  #     if i>= tauIdx[uIdx+1]
  #       uIdx += 1
  #     end
  #   end
  #
  #   # Generate Linearized Dynamics
  #   A[:,:,i] = linearizeDyn(nonlin_dyn, nonlin_dyn_deriv, xpts[1:end-1,i], uvec[:,uIdx])
  #
  #   # Compute Next Point in Simulation from warm starting definition
  #   xpts[:, i+1] = expm(A[:, :, i]*deltacomplete[i])*xpts[:, i]
  #
  # end



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

  # display(tauIdx)
  # display(delta0ws)

  # Initialize objective evaluator
  obj = Array(Float64, 0)
  deltaval = Array(Float64, N+1, 0)
  nobjeval = 0                           # Number of objective function evaluations
  ngradeval = 0                           # Number of gradient evaluations
  nhesseval = 0                           # Number of hessian evaluations

  # Construct NLPEvaluator
  STOev = nlinSTOev(x0, nx, A, N, t0, tf, tf, Q, Qf, uvec, ngrid, tgrid, tvec, tauIdx, tgridIdx, deltacomplete, nonlin_dyn, nonlin_dyn_deriv, IndTril, Itril, Jtril, Ag, Ig, Jg, Vg, bg, deltafun_prev, deltagrad_prev, deltahess_prev, xpts, expMat, Phi, M, S, C, obj, deltaval, nobjeval, ngradeval, nhesseval)


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

  # return tauopt, Jopt, solTime, stat
end

# # Solve Optimization for Nonlinear System
# function solve!(m::nlinSTO)
#
#   # Perform STO
#   m.soltime = @elapsed MathProgBase.optimize!(m.model)
#   m.stat = MathProgBase.status(m.model)
#
#   m.delta = MathProgBase.getsolution(m.model)
#   m.tau = delta2tau(m.delta, m.STOev.t0, m.STOev.tf)
#
#   # # Get delta and tau vectors (Complete and Not Complete)
#   # m.deltacomplete = MathProgBase.getsolution(m.model)
#   # m.taucomplete = delta2tau(m.deltacomplete, m.STOev.t0, m.STOev.tf)
#   # m.tau = m.taucomplete[m.nartsw+1:m.nartsw+1:end]  # N.B. nrep and not nrep+1 to start because uvec has one more element than tauopt
#   # m.delta = tau2delta(m.tau, m.STOev.t0, m.STOev.tf)
#
#   m.objval = MathProgBase.getobjval(m.model)
#
# end


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
# gettaucomplete(m::nlinSTO) = m.taucomplete
# getdeltacomplete(m::nlinSTO) = m.deltacomplete
getobjval(m::STO) = m.objval
getstat(m::STO) = m.stat
getsoltime(m::STO) = m.soltime
getnobjeval(m::STO) = m.STOev.nobjeval
getngradeval(m::STO) = m.STOev.ngradeval
getnhesseval(m::STO) = m.STOev.nhesseval


# Convert from Switching Times to Intervals
function tau2delta(tau::Array{Float64, 1}, t0::Float64, tf::Float64)
# function tau2delta(tau::Vector, t0, tf)

  # Extend vector of switching times with initial and final time
  tauexp = [t0; tau; tf]  # Extend tau vector to simplify numbering

  # Return Vector of differences
  return  diff(tauexp)

end


# Convert from Intervals to Switching Times
function delta2tau(delta::Array{Float64, 1}, t0::Float64)
# function delta2tau(delta::Vector, t0, tf)

  # Define tau vector
  tau = Array(Float64, length(delta)-1)
  # tau = Array(Float64, length(delta)-1)

  # Initialize first tau element
  tau[1] = t0 + delta[1]

  # Iteratively add all deltas to create switching times vector
  for i = 2:length(tau)
    tau[i] = tau[i-1] + delta[i]
  end

  tfdelta = tau[end] + delta[end]

  # tau = max(min(tau, tf-1e-10), t0+1e-10)  # Constrain vectors to be within the bounds of the grid

  return tau, tfdelta

  # Slower Version
  # Example
  #
  # [tau[1]]   [1 0 0] [t0 + delta[1]]
  # [tau[2]] = [1 1 0] [delta[2]]
  # [tau[3]]   [1 1 1] [delta[3]]
  #
  # N = length(delta)-1
  #
  # # Construct RHS
  # rhs = delta[1:end-1]
  # rhs[1] += t0
  #
  # return tril(ones(N, N))*rhs



end
