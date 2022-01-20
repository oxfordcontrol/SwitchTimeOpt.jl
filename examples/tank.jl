
using SwitchTimeOpt
using MathOptInterface
const MOI = MathOptInterface

using Plots
using Ipopt

# Define Solver options
maxiter = 300
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
#  Time vector
t0 = 0.0
tf = 10.0

# Integer inputs
uvec = repeat([1.0; 2.0], 8, 1)
uvec = copy(uvec')

# Number of switches
N = size(uvec,2) -1

# Cost function matrix
K = 1.0
C = sqrt(K)*[0.0 1.0 -1.0]
Q = C'*C

# Initial state
x0 = [2.0; 2.0; 3.0]


### Define system dynamics
function nldyn(x::Array{Float64,1}, u::Array{Float64,1})
  n = length(x)
  f = zeros(n)

  f[1] = -sqrt(x[1]) + u[1]
  f[2] = sqrt(x[1]) - sqrt(x[2])
  f[3] = -0.05
  return f

end

function nldyn_deriv(x::Array{Float64,1}, u::Array{Float64,1})
  df = [-1/(2*sqrt(x[1]))        0                      0;
        1/(2*sqrt(x[1]))        -1/(2*sqrt(x[2]))       0;
        0                        0                      0]
end




# Generate and solve problem
ngrid = 100

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
plotSolution(m, "Tank", "Tank", showPlots, savePlots)