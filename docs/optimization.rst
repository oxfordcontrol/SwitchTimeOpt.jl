===============================
Optimization
===============================
We now describe the interface for defining and solving switching time optimization problems.

Problem Definition
==================
This package allows us to define and solve problems in the form

.. math::
  \begin{array}{ll}
    \underset{\delta}{\mbox{minimize}} & \int_{t_0}^{T_\delta} x(t)^\top Q x(t)\; \mathrm{d}t + x(T_\delta)^\top Q x(T_\delta)\\
    \mbox{subject to} & \dot{x}(t) = f_i(x(t)) \quad t\in[\tau_i,\tau_{i+1}) \quad i = 0,\dots,N\\
    & x(0) = x_0\\
    & \delta \in \Delta
  \end{array}


where the decision variable is the vector of :math:`N+1` intervals :math:`\delta = \begin{bmatrix}\delta_0 & \dots & \delta_{N}\end{bmatrix}^\top\in \mathbb{R}^{N+1}` such that :math:`\delta_i = \tau_{i+1} - \tau_i`. Each interval :math:`\delta_i` defines how long the :math:`i`-th dynamics are active. The state trajectory is :math:`x(t) \in \mathbb{R}^{n}`. The value :math:`T_\delta` is defined as the final time when intervals :math:`\delta` are applied, i.e.

.. math::

  T_\delta = \sum_{i=0}^{N}\delta_i.


The set :math:`\Delta` defines the set of feasible intervals

.. math::
  \Delta = \left\{\delta \in \mathbb{R}^{N+1} \;\middle|\; T_\delta = T \wedge 0\leq lb_i \leq \delta_i \leq ub_i\; \forall i\right\}

The variable :math:`T` defines the desired final time of the interval. The scalars :math:`lb_i` and :math:`ub_i` define additional constraints on the interval in case we would like to have a minimum or a maximum time in which the :math:`i`-th dynamics are active.


Linear Dynamics
--------------------

In case when the dynamics are linear of the form

.. math::
  \dot{x}(t) = A_ix(t), \quad t\in [\tau_i,\tau_{i+1})

we can define a 3-dimensional matrix :code:`A` whose slices :code:`A[:,:,i]` represent dynamics :math:`A_{i-1}`:

::

  A = Array(Float64, n, n, N+1)
  A[:, :, i] = ...   # Dynamics A_{i-1}
  ...


The switching time optimization problem can be quickle defined as

::

  p = stoproblem(x0, A)

Where :code:`x0` is the initial state vector :math:`x_0` and :code:`A` is the 3-dimensional matrix defining the dynamics.


Noninear Dynamics
-------------------

Given a nonlinear system defined by dynamics

.. math::

  \dot{x}(t) = f(x(t), u(t))

where :math:`u(t)` the input vector assuming integer values :math:`u_i` between switching instants

.. math::

  u(t) = u_i \quad t\in [\tau_i, \tau_{i+1}),

we can define our switched nonlinear system as

.. math::

  \dot{x}(t) = f(x(t), u_i)  = f_i(x(t)) \quad t\in [\tau_i, \tau_{i+1}).

To create the optimization problem we need to define the nonlinear dynamics by means of an additional function

::

  function nldyn(x, ui)
    ...
  end

returning the vector of states derivatives. The variable :code:`x` is the state :math:`x(t)` and :code:`ui` is the input vector :math:`u_i`. Moreover, we need to define the jacobian of the switched dynamics with respect to the system states

.. math::

  J_{f_i} = \frac{\partial f_i (x(t))}{\partial x(t)}

by means of an additional function

::

  function jac_nldyn(x, ui)
    ...
  end


Note that the function :code:`jac_nldyn` returns a matrix having in each row the gradient of every component of the function :math:`f_i(x(t))` with respect to each state component. Last  element necessary to construct the matrix :code:`U` having a column each integer input vector :code:`ui`. Then, we can define the switching time optimization problem as:

::

  p = stoproblem(x0, nldyn, jac_nldyn, U)


.. note::
  The nonliner switched system optimization operates by introducing additional linearization points at an equally spaced linearization grid. To set the number of linearization points to :math:`100` for example, it is just necessary to add an extra argument to the previous function call as follows:
  ::

    p = stoproblem(x0, nldyn, jac_nldyn, U, ngrid = 100)

  where :code:`ngrid` defines the number of linearization points.

Optional Arguments
---------------------
There are many additional keyword arguments that can be be passed to the :code:`stoproblem(...)` function to customize the optimization problem.

+--------------------------+----------------------------------------+----------------------------------------------------+
|Parameter                 | Description                            | Default value                                      |
+==========================+========================================+====================================================+
|:code:`t0`                | Initial Time :math:`t_0`               | :code:`0.0`                                        |
+--------------------------+----------------------------------------+----------------------------------------------------+
|:code:`tf`                | Final Time :math:`t_0`                 | :code:`1.0`                                        |
+--------------------------+----------------------------------------+----------------------------------------------------+
|:code:`Q`                 | Cost matrix :math:`Q`                  | :code:`eye(n)`                                     |
+--------------------------+----------------------------------------+----------------------------------------------------+
|:code:`lb`                | Vector of lower bounds :math:`lb_i`    | :code:`zeros(N+1)`                                 |
+--------------------------+----------------------------------------+----------------------------------------------------+
|:code:`ub`                | Vector of lower bounds :math:`ub_i`    | :code:`Inf*ones(N+1)`                              |
+--------------------------+----------------------------------------+----------------------------------------------------+
|:code:`tau0ws`            | Warm starting initial switching times  | Equally spaced between :code:`t0` and :code:`tf`   |
+--------------------------+----------------------------------------+----------------------------------------------------+
|:code:`solver`            | MathProgbase.jl solver                 | :code:`IpoptSolver()`                              |
+--------------------------+----------------------------------------+----------------------------------------------------+



Problem Solution
======================

Once the problem is defined, it can be solved by simply running

::

  solve!(p)



Choosing Solver
-----------------

Any NLP solver supported by `JuliaOpt <http://www.juliaopt.org/>`_ may be used through `MathProgBase.jl <https://github.com/JuliaOpt/MathProgBase.jl/>`_ interface. The default solver is `Ipopt <https://github.com/JuliaOpt/Ipopt.jl/>`_. To use `KNITRO <https://github.com/JuliaOpt/KNITRO.jl/>`_ solver with the linear example, it is just necessary to specify an :code:`AbstractMathProgSolver` object (see `here <http://mathprogbasejl.readthedocs.io/en/latest/solvers.html>`_ for more details) when the problem is created

::

  using KNITRO
  p = stoproblem(x0, A, solver = KnitroSolver())

All the solver-specific options can be passed when creating the :code:`AbstractMathProgSolver` object: algorithm types (first/second order methods), tolerances, verbosity and so on.

Obtaining Results
-----------------

The optimal cost function and the optimal switching times and intervals can be obtained as follows:
::

  objval = getobjval(p)
  tauopt = gettau(p)
  deltaopt = getdelta(p)


We can get the execution time (including the time for the function calls) and the status of the solver by executing:

::

  stat = getstat(p)
  soltime = getsoltime(p)


Optimizing in a Loop
---------------------
The toolbox is suited for receeding horizon implementations. To run the optimization in a loop it is just necessary to update the value of the current state :code:`x0` and to update the warm starting point :code:`tau0ws` which is usually chosen as the optimal solution at the previous optimizaton.

To set the initial state at :code:`x0` it is just necessary to return

::

  setx0!(m, x0)


We can set the warm starting point at :code:`tau0ws` with

::

  setwarmstart!(m, tau0ws)
