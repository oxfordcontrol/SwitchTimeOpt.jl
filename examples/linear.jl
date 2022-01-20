using SwitchTimeOpt
using MathOptInterface
const MOI = MathOptInterface

using Plots
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
# Size of the state space
nx = 2;

# Number of switches
N = 5;

# Initial state
x0 = [1.0; 1.0]

# Define initial and final time
t0 = 0.0
tf = 1.0


### Define system dynamics
A = zeros(nx, nx, N+1)

A[:, :, 1] = [-1 0;
              1  2]
A[:, :, 2] = [1 1;
              1 -2]
for i = 3:N+1
  A[:, :, i] = A[:, :, mod(i+1, 2)+1]
end

### Define and solve switching time optimization problem
# m = stoproblem(x0, A, solver=solver)
#
m = stoproblem(x0, A, t0=t0, tf=tf, ngrid=100)

# Solve problem
solve!(m)

# Obtain values
tauopt = gettau(m)
Jopt = getobjval(m)
soltime = getsoltime(m)

# Simulate linear system
xsim, xpts, Jsim, t = simulate(m)

println("Objective: ", Jsim)

p1 = plot(title="Linear example", titlefont=font(10))
for i=1:nx
  plot!(t, xsim[i,:], label="x"*string(i))
end
plot!(grid=true, ylim=[(1-0.2*sign(minimum(xsim)))*minimum(xsim), (1+0.2*sign(maximum(xsim)))*maximum(xsim)])
display(plot(p1))

