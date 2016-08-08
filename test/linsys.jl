# Test Linear System Optimizaiton withh one switching time
function linsystest(solver=IpoptSolver(tol=1e-03); objtol = 1e-3, primaltol = 1e36)
  println("Testing Switched Linear Systems Optimization ", string(typeof(solver)))

  # Size of the state space
  nx = 2;

  # Cost function Matrix
  Q = eye(2)

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
  # @test getstat(m) == :Optimal
  # @test_approx_eq_eps gettau(m) 0.2648660991124608 primaltol
  # @test_approx_eq_eps getobjval(m) 2.6272714724616115 objtol

  println("Passed")
end
