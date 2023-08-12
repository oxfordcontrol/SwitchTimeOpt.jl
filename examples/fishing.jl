using SwitchTimeOpt
using MathOptInterface
const MOI = MathOptInterface
using Ipopt

# Define Solver options
maxiter = 100
maxtime = 1000.0
verbose = 5
tolerance = 1e-06

#Define solver
solver = Ipopt.Optimizer()
MOI.set(solver, MOI.RawOptimizerAttribute("tol"), tolerance)
MOI.set(solver, MOI.RawOptimizerAttribute("print_level"), verbose)
MOI.set(solver, MOI.RawOptimizerAttribute("max_cpu_time"), maxtime)
MOI.set(solver, MOI.RawOptimizerAttribute("max_iter"), maxiter)

### Define system parameters
# Time vector
t0 = 0.0
tf = 12.0

# Integer input
uvec = [repeat([0.0; 1.0], 4, 1); 0.0] 
uvec = copy(uvec')

# Number of switching times
N = size(uvec,2) - 1

# Cost funcction matrix
C = [1.0 0.0 -1.0;
     0.0 1.0 -1.0]
Q = C'*C

# Define initial state
x0 = [0.5; 0.7; 1]


### Define system dynamics
function nldyn(x::Array{Float64,1}, u::Array{Float64,1})
  n = length(x)
  f = zeros(n)
  f[1] = x[1] - x[1]*x[2] - 0.4*x[1]*u[1]
  f[2] = -x[2] + x[1]*x[2] - 0.2*x[2]*u[1]
  return f
end

function nldyn_deriv(x::Array{Float64,1}, u::Array{Float64,1})
  df = [1.0-x[2]-0.4*u[1]       -x[1]                   0;
        x[2]                     -1+x[1]-0.2*u[1]       0;
        0                        0                      0]
  return df
end


# Generate and solve problem
ngrid = 300

# Initialize the problem
m = stoproblem( 
  x0,                 # Initial state
  nldyn,              # Nonlinear dynamics
  nldyn_deriv,        # Nonlinear dynamics derivative
  uvec,               # Vector of integer inputs
  ngrid = ngrid,      # Number of grid points
  t0 = t0,            # InitialtTime
  tf = tf,            # Final time
  Q = Q,              # Cost matrix
  solver = solver)

solve!(m)
objlin = getobjval(m)
tauopt = gettau(m)
status = getstat(m)
xsim, xpts, objode45, t = simulate(m, tauopt)

println("Status: ", status)
println("Objective value: ", objlin)

showPlots = true
savePlots = false
plotSolution(m, "Lotka-Volterra", "Lotka-Volterra", showPlots, savePlots)