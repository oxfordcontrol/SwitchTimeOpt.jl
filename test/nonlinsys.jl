# Test Fishing Problem - Switched Nonlinear System
function nonlinsystest(solver=IpoptSolver(tol=1e-04); objtol = 1e-3, primaltol = 1e-3)
  println("Testing Switched Nonlinear Systems Optimization ", string(typeof(solver)))

  # Define Time Interval
  t0 = 0.0; tf = 12.0

  nartsw = 6  # Number of artificial switchings per switching
  uvec = [repmat([0.0; 1.0], 4, 1); 0.0]'  # Input Vector

  # Cost Funcction Matrix
  C = [1.0 0.0 -1.0 0.0;
       0.0 1.0 0.0  -1.0]
  Q = C'*C

  # Define Initial State
  x0 = [0.5; 0.7; 1; 1]

  ### Define System Dynamics
  function nldyn(x::Array{Float64,1}, u::Array{Float64,1})
    n = length(x)
    f = zeros(n)
    f[1] = x[1] - x[1]*x[2] - 0.4*x[1]*u[1]
    f[2] = -x[2] + x[1]*x[2] - 0.2*x[2]*u[1]
    return f

  end

  function nldyn_deriv(x::Array{Float64,1}, u::Array{Float64,1})
    df = [1.0-x[2]-0.4*u[1]       -x[1]                   0    0;
          x[2]                     -1+x[1]-0.2*u[1]       0    0;
          0                        0                      0    0;
          0                        0                      0    0]
  end

  m = stoproblem(x0, nldyn, nldyn_deriv, uvec, nartsw, t0=t0, tf=tf, Q=Q, solver=solver)


  solve!(m)
  # @test getstat(m) == :Optimal

  # Test Optimal Solution
  tauopt = gettau(m)
  @test_approx_eq_eps tauopt[1] 2.415555396457526 primaltol
  @test_approx_eq_eps tauopt[2] 4.229649407353941 primaltol
  @test_approx_eq_eps tauopt[3] 4.833725131673787 primaltol
  @test_approx_eq_eps tauopt[4] 5.221665532976622  primaltol
  @test_approx_eq_eps tauopt[5] 6.607137544331319 primaltol
  @test_approx_eq_eps tauopt[6] 6.737507093345905 primaltol
  @test_approx_eq_eps tauopt[7] 9.357656064178713 primaltol
  @test_approx_eq_eps tauopt[8] 9.377554759276379 primaltol

  # Test Optimal Value
  @test_approx_eq_eps getobjval(m) 0.6785595836871159 objtol

  println("Passed")
end
