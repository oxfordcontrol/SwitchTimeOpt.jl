================================
Quick Example
================================

Consider the Switching time Optimization Problem in the form


.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2}\int_0^{T} x(t)^\top Q x(t)\; \mathrm{d}t \\
    \mbox{subject to} & \dot{x}(t) = \begin{cases}
    A_0 x(t) & t< \tau\\
    A_1 x(t) & t\geq \tau
    \end{cases}\\
    & 0\leq \tau \leq T
  \end{array}

with variable :math:`\tau\in \mathbb{R}` and dynamics defined by matrices :math:`A_0,A_1\in \mathbb{R}^{n\times n}`.


This problem can be solved by SwitchTimeOpt.jl as follows


.. code-block:: julia

  using SwitchTimeOpt
  using Ipopt

  # Time Interval
  t0 = 0.0; tf = 1.0

  # Cost function Matrix
  Q = eye(2)

  # Initial State
  x0 = [1.0; 1.0]

  # Dynamics
  A = Array(Float64, 2, 2, 2)
  A[:, :, 1] = randn(2, 2)  # A_0 matrix
  A[:, :, 2] = randn(2, 2)  # A_1 matrix

  # Create Problem
  m = createsto(x0, A, t0=t0, tf=tf, Q=Q)

  # Solve Problem
  solve!(m)

  # Get optimal Solution
  tauopt = gettau(m)

  # Get optimum value
  objval = getobjval(m)
