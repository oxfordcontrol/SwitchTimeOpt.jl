# Test Linear System Optimizaiton
using SwitchTimeOpt
using MathOptInterface
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



# Define Time Interval
t0 = 0.0; tf = 12.0

uvec = [repeat([0.0; 1.0], 4, 1); 0.0]  # Input Vector
uvec = copy(uvec')

# Cost Funcction Matrix
C = [1.0 0.0 -1.0;
      0.0 1.0 -1.0]
Q = C'*C

# Define Initial State
x0 = [0.5; 0.7; 1]

### Define System Dynamics
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
end

m = stoproblem(x0, nldyn, nldyn_deriv, uvec, ngrid = 300, t0=t0, tf=tf, Q=Q, solver=solver)


solve!(m)

# Test Optimal Solution
@testset "Test optimal switching times" begin
  tauopt = gettau(m)

  @test isapprox(tauopt[1], 2.4435555103155107, atol=primaltol)
  @test isapprox(tauopt[2], 4.122622390312436, atol=primaltol)
  @test isapprox(tauopt[3], 4.433725972515174, atol=primaltol)
  @test isapprox(tauopt[4], 4.682336152264894,  atol=primaltol)
  @test isapprox(tauopt[5], 5.199622121419033, atol=primaltol)
  @test isapprox(tauopt[6], 5.367438990179279, atol=primaltol)
  @test isapprox(tauopt[7], 6.370950931911677, atol=primaltol)
  @test isapprox(tauopt[8], 6.4660389494043295, atol=primaltol)
end

@testset "Test status and optimal objective value" begin
  # Test Optimal Value
  @test isapprox(getobjval(m), 1.3454355602054182, atol=objtol)
  @test string(getstat(m)) == "LOCALLY_SOLVED"
end
