# Test Linear System Optimizaiton withh one switching time
function linsystest(solver=IpoptSolver(tol=1e-03); objtol = 1e-3, primaltol = 1e-3)
  println("Testing Switched Linear Systems Optimization ", string(typeof(solver)))

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
  @test isapprox(gettau(m)[1], 0.26486646235103123, atol=primaltol)
  @test isapprox(getobjval(m), 5.2545429449272145, atol=objtol)

  println("Passed")
end
