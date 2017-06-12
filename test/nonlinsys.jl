# Test Fishing Problem - Switched Nonlinear System
function nonlinsystest(solver=IpoptSolver(tol=1e-04, max_iter = 25); objtol = 1e-3, primaltol = 1e-3)
  println("Testing Switched Nonlinear Systems Optimization ", string(typeof(solver)))

  # Define Time Interval
  t0 = 0.0; tf = 12.0

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

  m = stoproblem(x0, nldyn, nldyn_deriv, uvec, ngrid = 300, t0=t0, tf=tf, Q=Q, solver=solver)


  solve!(m)
  # @test getstat(m) == :Optimal

  # Test Optimal Solution
  tauopt = gettau(m)
  @test isapprox(tauopt[1], 2.4434718158206916, atol=primaltol)
  @test isapprox(tauopt[2], 4.122362522221089, atol=primaltol)
  @test isapprox(tauopt[3], 4.433072844726073, atol=primaltol)
  @test isapprox(tauopt[4], 4.68252002334405,  atol=primaltol)
  @test isapprox(tauopt[5], 5.201667827750571, atol=primaltol)
  @test isapprox(tauopt[6], 5.369908894975762, atol=primaltol)
  @test isapprox(tauopt[7], 6.376845190917198, atol=primaltol)
  @test isapprox(tauopt[8], 6.47160340206027, atol=primaltol)

  # Test Optimal Value
  @test isapprox(getobjval(m), 1.3454355602054182, atol=objtol)

  println("Passed")
end
