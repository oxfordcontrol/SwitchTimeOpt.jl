# Test Linear System Optimizaiton
using SwitchTimeOpt
using MathOptInterface
using LinearAlgebra
const MOI = MathOptInterface
using Ipopt

# Define test options
objtol = 1e-04
primaltol = 1e-03

# Define Solver options
maxiter = 50
maxtime = 100.0
verbose = 0
tolerance = 1e-06

#Define solver
solver = Ipopt.Optimizer()
MOI.set(solver, MOI.RawOptimizerAttribute("tol"), tolerance)
MOI.set(solver, MOI.RawOptimizerAttribute("print_level"), verbose)
MOI.set(solver, MOI.RawOptimizerAttribute("max_cpu_time"), maxtime)
MOI.set(solver, MOI.RawOptimizerAttribute("max_iter"), maxiter)



# Size of the state space
nx = 2;

# Cost function Matrix
Q = 1.0 * Matrix(I, nx, nx)

# Initial State
x0 = [1.0; 1.0]


# Define initial and final time
t0 = 0.0
tf = 1.0

### Define System Dynamics
A = zeros(nx, nx, 2)
A[:, :, 1] = [-1 0;
              1  2]
A[:, :, 2] = [1 1;
              1 -2]

m = stoproblem(x0, A, solver=solver)
solve!(m)

@testset "Test optimal switching time" begin 
  @test isapprox(gettau(m)[1], 0.26486646235103123, atol=primaltol)
end

@testset "Test status and optimal objective value" begin 
  @test isapprox(getobjval(m), 5.2545429449272145, atol=objtol)
  @test string(getstat(m)) == "LOCALLY_SOLVED"
end
