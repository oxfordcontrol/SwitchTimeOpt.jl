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
  solver::MathOptInterface.AbstractOptimizer=Ipopt.Optimizer())

  # Get Dimensions
  nx = size(A, 1)      # State Dimension
  N = size(A, 3) - 1   # Get number of switches

  # Adjust variables which have not been initalized
  if isempty(Q)
    Q = 1.0 * Matrix(I, nx, nx)
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
    tau0ws = collect(range(t0, stop=tf, length=N+2))
    tau0ws = tau0ws[2:end-1]
  end

  # Define warm starting delta0
  delta0ws = tau2delta(tau0ws, t0, tf)

  # Create Discretization Grid
  tgrid = collect(range(t0, stop=tf, length=ngrid))

  # Initialize time vectors
  tvec = Array{Float64}(undef, N + ngrid)    # Complete grid
  tauIdx = Array{Int}(undef, N + 2)      # Indeces of switching times in the complete grid
  tgridIdx = Array{Int}(undef, ngrid)      # Indeces of switching times in the complete grid
  deltacomplete = Array{Float64}(undef, N + ngrid - 1)   # Intervals over the whole grid


  ### Initialize NLP Evaluator
  # Preallocate arrays
  deltafun_prev = Array{Float64}(undef, N+1)
  deltagrad_prev = Array{Float64}(undef, N+1)
  deltahess_prev = Array{Float64}(undef, N+1)
  xpts = Array{Float64}(undef, nx, N+ngrid); xpts[:, 1] = x0   # Set Initial State
  expMat = Array{Float64}(undef, nx, nx, N+ngrid-1)
  Phi = Array{Float64}(undef, nx, nx, N+2, N+2)
  M = Array{Float64}(undef, nx, nx, N+ngrid-1)
  S = Array{Float64}(undef, nx, nx, N+ngrid)
  C = Array{Float64}(undef, nx, nx, N+1)

  # Decompose Dynamics Matrices
  V = Array{Complex{Float64}}(undef, 2*nx, 2*nx, N+1)
  invV = Array{Complex{Float64}}(undef, 2*nx, 2*nx, N+1)
  D =  Array{Complex{Float64}}(undef, 2*nx, N+1)
  isDiag = Array{Bool}(undef, N+1)

  for i = 1:N+1
    D[:, i], V[:, :, i] = eigen([-A[:, :, i]'  Q;
                             zeros(nx, nx) A[:, :, i]])
    if cond(V[:, :, i]) == Inf  # Non diagonalizable matrix
      isDiag[i] = false
    else
      invV[:, :, i] = inv(V[:, :, i])
      isDiag[i] = true
    end
  end


  # Construct Matrix of Indeces for lower triangular Matrix (Hessian)
  IndTril = (LinearIndices(tril(ones(N+1, N+1))))[findall(!iszero,tril(ones(N+1, N+1)))]
  Itril, Jtril, _ = begin
  I_temp = findall(!iszero, tril(ones(N+1, N+1)))
  (getindex.(I_temp, 1), getindex.(I_temp, 2), tril(ones(N+1, N+1))[I_temp])
  end

  # Construct Constraints Matrix (Vector)
  Ag = ones(1, N+1)
  Ig, Jg, Vg = begin
    I_temp = findall(!iszero, Ag)
    (getindex.(I_temp, 1), getindex.(I_temp, 2), Ag[I_temp])
  end
  bg = [tf]   # Only one constraint for the sum of the switching intervals


  # Initialize objective evaluator
  obj = Array{Float64}(undef, 0)
  deltaval = Array{Float64}(undef, N+1, 0)
  nobjeval = 0                           # Number of objective function evaluations
  ngradeval = 0                           # Number of gradient evaluations
  nhesseval = 0                           # Number of hessian evaluations

  # Construct NLP evaluator
  STOev = linSTOev(x0, nx, A, N, t0, tf, tf, Q, E, ngrid, tgrid, tvec,
          tauIdx, tgridIdx, deltacomplete, V, invV, D, isDiag, IndTril,
          Itril, Jtril, Ag, Ig, Jg, Vg, bg, deltafun_prev, deltagrad_prev, 
          deltahess_prev, xpts, expMat, Phi, M, S, C,
          obj, deltaval, nobjeval, ngradeval, nhesseval)


model = solver
MathOptInterface.empty!(model)
delta = MathOptInterface.add_variables(model, N+1)
for w_i in delta
  MathOptInterface.add_constraint(model, MathOptInterface.SingleVariable(w_i), MathOptInterface.GreaterThan(0.0))
  MathOptInterface.add_constraint(model, MathOptInterface.SingleVariable(w_i), MathOptInterface.LessThan(tf))
end

lbub = [MathOptInterface.NLPBoundsPair(tf, tf)]
MathOptInterface.set(model, MathOptInterface.VariablePrimalStart(), delta, delta0ws)
block_data = MathOptInterface.NLPBlockData(
  lbub,
  STOev,
  true
)
MathOptInterface.set(model, MathOptInterface.NLPBlock(), block_data)
MathOptInterface.set(model, MathOptInterface.ObjectiveSense(), MathOptInterface.MIN_SENSE)

# Propagate dynamic for new switching times
propagateDynamics!(STOev, delta0ws)

# Create STO
STOproblem = linSTO(model, STOev, delta0ws, delta)
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
  solver::MathOptInterface.AbstractOptimizer=Ipopt.Optimizer())

  # Get Dimensions
  nx = length(x0)         # State Dimension


  # Generate artificial switches repeating the input sequence
  N = size(uvec, 2) - 1   # Get total number of switches

  # Adjust variables which have not been initalized
  if isempty(Q)
    Q = 1.0 * Matrix(I, nx, nx)
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
    tau0ws = collect(range(t0, stop=tf, length=N+2))  # Currently counting tau_0 and tau_{N+1}. They are removed below after xpts initialization.
    tau0ws = tau0ws[2:end-1]  # Include only switching instants

  end



  # Create Discretization grid
  tgrid = collect(range(t0, stop=tf, length=ngrid))


  # Initialize time vectors
  tvec = Array{Float64}(undef, N + ngrid)    # Complete grid
  tauIdx = Array{Int}(undef, N + 2)      # Indeces of switching times in the complete grid
  tgridIdx = Array{Int}(undef, ngrid)      # Indeces of switching times in the complete grid
  deltacomplete = Array{Float64}(undef, N + ngrid - 1)   # Intervals over the whole grid

  ### Initialize NLP Evaluator
  # Extend System State To get Affine to Linear Dynamics
  nx = nx + 1

  # Extend Initial State and Cost Matrix
  x0 = [x0; 1]
  spz = Array{Float64}(undef, 1,1); spz[1,1] = 0.0; spz = sparse(spz)  # Sparse scalar
  Q = Matrix(blockdiag(sparse(Q), spz))
  E = Matrix(blockdiag(sparse(E), spz))


  # Define Required Matrices and Switching Instants
  A = Array{Float64}(undef, nx, nx, N+ngrid-1)

  # Initialize Switching Instants
  xpts = Array{Float64}(undef, nx, N+ngrid)
  xpts[:, 1] = x0



  # Define warm starting delta0ws from tau0ws
  delta0ws = tau2delta(tau0ws, t0, tf)


  deltafun_prev = Array{Float64}(undef, N+1)
  deltagrad_prev = Array{Float64}(undef, N+1)
  deltahess_prev = Array{Float64}(undef, N+1)
  expMat = Array{Float64}(undef, nx, nx, N+ngrid-1)
  Phi = Array{Float64}(undef, nx, nx, N+2, N+2)
  M = Array{Float64}(undef, nx, nx, N+ngrid-1)
  S = Array{Float64}(undef, nx, nx, N+ngrid)
  C = Array{Float64}(undef, nx, nx, N+1)  # Only at sw times (including initial time)


  # Construct Matrix of Indeces for lower triangular Matrix (Hessian)
  IndTril = (LinearIndices(tril(ones(N+1, N+1))))[findall(!iszero,tril(ones(N+1, N+1)))]
  Itril, Jtril, _ = begin
  I_temp = findall(!iszero, tril(ones(N+1, N+1)))
  (getindex.(I_temp, 1), getindex.(I_temp, 2), tril(ones(N+1, N+1))[I_temp])
  end


  # Construct Constraints Matrix (Vector)
  Ag = ones(1, N+1)
  Ig, Jg, Vg = begin
  I_temp = findall(!iszero, Ag)
  (getindex.(I_temp, 1), getindex.(I_temp, 2), Ag[I_temp])
  end
  bg = [tf]   # Only one constraints for the sum of the switching intervals

  # Initialize objective evaluator
  obj = Array{Float64}(undef, 0)
  deltaval = Array{Float64}(undef, N+1, 0)
  nobjeval = 0                           # Number of objective function evaluations
  ngradeval = 0                           # Number of gradient evaluations
  nhesseval = 0                           # Number of hessian evaluations

  # Construct NLPEvaluator
  STOev = nlinSTOev(x0, nx, A, N, t0, tf, tf, Q, E, uvec, ngrid, tgrid, tvec, tauIdx, tgridIdx, deltacomplete, nonlin_dyn,
                 nonlin_dyn_deriv, IndTril, Itril, Jtril, Ag, Ig, Jg, Vg, bg, deltafun_prev, deltagrad_prev, deltahess_prev,
                 xpts, expMat, Phi, M, S, C, obj, deltaval, nobjeval, ngradeval, nhesseval)


  # Propagate Dynamics to compute matrix exponentials and states at the switching times
  propagateDynamics!(STOev, delta0ws)

  # Generate Model
  model = solver
  MathOptInterface.empty!(model)
  delta = MathOptInterface.add_variables(model, N+1)
  for i=1:N+1
    MathOptInterface.add_constraint(model, MathOptInterface.SingleVariable(delta[i]), MathOptInterface.GreaterThan(lb[i]))
    MathOptInterface.add_constraint(model, MathOptInterface.SingleVariable(delta[i]), MathOptInterface.LessThan(ub[i]))
  end

  lbub = [MathOptInterface.NLPBoundsPair(tf, tf)]
  MathOptInterface.set(model, MathOptInterface.VariablePrimalStart(), delta, delta0ws)
  block_data = MathOptInterface.NLPBlockData(
    lbub,
    STOev,
    true
  )
  MathOptInterface.set(model, MathOptInterface.NLPBlock(), block_data)
  MathOptInterface.set(model, MathOptInterface.ObjectiveSense(), MathOptInterface.MIN_SENSE)

  # Create STO
  STOproblem = nlinSTO(model, STOev, delta0ws, delta)

  return STOproblem  # Return STO

end


# Solve Optimization for Linear System
function solve!(m::STO)

  # Perform STO
  m.soltime = @elapsed MathOptInterface.optimize!(m.model)
  m.stat = MathOptInterface.get(m.model, MathOptInterface.TerminationStatus())
  m.delta = MathOptInterface.get(m.model, MathOptInterface.VariablePrimal(), m.delta_MOI)
  m.objval = MathOptInterface.get(m.model, MathOptInterface.ObjectiveValue())
  m.tau, _ = delta2tau(m.delta, m.STOev.t0)

end

"Set warm starting point"
function setwarmstart!(m::STO, tau0ws::Array{Float64,1})

  # Define warm starting delta0
  delta0ws = tau2delta(tau0ws, m.STOev.t0, m.STOev.tf)

  # Set Warm Starting Point for Nonlinear Solver
  MathOptInterface.set(m.model, VariablePrimalStart, delta_MOI, delta0ws)

  # Reset Cost Function Iterates
  m.STOev.obj = Array{Float64}(undef, 0)
  m.STOev.deltaval = Array{Float64}(undef, m.STOev.N+1, 0)

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
  tau = Array{Float64}(undef, length(delta)-1)

  # Initialize first tau element
  tau[1] = t0 + delta[1]

  # Iteratively add all deltas to create switching times vector
  for i = 2:length(tau)
    tau[i] = tau[i-1] + delta[i]
  end

  tfdelta = tau[end] + delta[end]

  return tau, tfdelta

end


function plotSolution(d::nlinSTO, savename::String, title::String, show=true, save=false)

  tauopt = gettau(d)
  xsim, ~, ~, t = simulate(d, tauopt)
  usim, ~ = simulateinput(d, t)

  n_omega = size(usim)[1]
  nx = d.STOev.nx -1
  # Plot solution

  l = @layout [a;b]
  p1 = plot(title=title, titlefont=font(10))
  for i=1:nx
    plot!(t, xsim[i,:], label="x"*string(i))
  end
  plot!(grid=true, ylim=[(1-0.2*sign(minimum(xsim)))*minimum(xsim), (1+0.2*sign(maximum(xsim)))*maximum(xsim)])

  p2 = plot()
  for i=1:n_omega
    plot!(t,usim[i,:], label="u"*string(i), linetype=:steppre)
  end
  plot!(ylim=[(1-0.2*sign(minimum(xsim)))*minimum(usim), (1+0.2*sign(maximum(xsim)))*maximum(usim)], grid=true, xlabel="time")
  
  if show && save
    display(plot(p1, p2, layout=l))
    savefig(savename*".pdf")
  elseif save && !show
    plot(p1, p2, layout=l)
    savefig(savename*".pdf")
  elseif !save && show
    display(plot(p1, p2, layout=l))
  end
  return
end

function plotSolution(d::linSTO, u::Array{Float64,2}, savename::String, title::String, show=true, save=false)

  tauopt = gettau(d)
  xsim, ~, ~, t = simulate(d, tauopt)
  println(length(t))
  n_omega = size(u)[1]
  nx = d.STOev.nx
  # Plot solution
  l = @layout [a;b]
  p1 = plot(title=title, titlefont=font(10))
  for i=1:nx
    plot!(t, xsim[i,:], label="x"*string(i))
  end
  plot!(grid=true, ylim=[(1-0.2*sign(minimum(xsim)))*minimum(xsim), (1+0.2*sign(maximum(xsim)))*maximum(xsim)])

  p2 = plot()
  tSwitch = [d.STOev.t0; tauopt; d.STOev.tf]
  println(length(tSwitch))
  println(length([u[1,1]; u[1,:]]))
  for i=1:n_omega
    plot!(tSwitch, [u[i,1]; u[i,:]], label="A"*string(i), linetype=:steppre)
  end
  plot!(ylim=[(1-0.2*sign(minimum(xsim)))*minimum(u), (1+0.2*sign(maximum(xsim)))*maximum(u)], grid=true, xlabel="time")
  
  if show && save
    display(plot(p1, p2, layout=l))
    savefig(savename*".pdf")
  elseif save && !show
    plot(p1, p2, layout=l)
    savefig(savename*".pdf")
  elseif !save && show
    display(plot(p1, p2, layout=l))
  end
  return
end