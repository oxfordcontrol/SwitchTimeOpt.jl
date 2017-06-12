================================
Quick Example
================================

Consider the switching time optimization problem in the form


.. math::
  \begin{array}{ll}
    \underset{\tau}{\mbox{minimize}} & \int_{t_0}^{t_f} \|x(t)\|_2^2\; \mathrm{d}t \\
    \mbox{subject to} & \dot{x}(t) = \begin{cases}
    A_0 x(t) & t< \tau\\
    A_1 x(t) & t\geq \tau
    \end{cases}\\
    & x(0) = x_0\\
    & 0\leq \tau \leq T
  \end{array}

with variable :math:`\tau\in \mathbb{R}`, the dynamics defined by matrices :math:`A_0,A_1\in \mathbb{R}^{n\times n}` and the initial state :math:`x_0\in \mathbb{R}^{n}`.


This problem can be solved by SwitchTimeOpt.jl as follows


::

  using SwitchTimeOpt
  using Ipopt

  # Time Interval
  t0 = 0.0; tf = 1.0

  # Initial State
  x0 = [1.0; 1.0]

  # Dynamics
  A = Array{Float64}(2, 2, 2)
  A[:, :, 1] = randn(2, 2)  # A_0 matrix
  A[:, :, 2] = randn(2, 2)  # A_1 matrix

  # Create Problem
  m = stoproblem(x0, A, t0=t0, tf=tf)

  # Solve Problem
  solve!(m)

  # Get optimal Solution
  tauopt = gettau(m)

  # Get optimum value
  objval = getobjval(m)
