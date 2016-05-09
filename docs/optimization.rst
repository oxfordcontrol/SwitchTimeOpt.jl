===============================
Optimization
===============================
We now describe the interface for defining and solving switching time optimization problems.

Problem Definition
==================
This package allows us to define and solve problems in the form

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2}\int_{t_0}^{t_f} x(t)^\top Q x(t)\; \mathrm{d}t \\
    \mbox{subject to} & \dot{x}(t) = f_i(x(t)) \quad t\in[\tau_i,\tau_{i+1}) \quad i = 0,\dots,N\\
    & x(0) = x_0\\
    & l_i \leq \tau_{i+1} - \tau_i \leq u_i,\quad i = 0,\dots,N\\
    &\tau_0 = t_0,\;\tau_{N+1} = t_f
  \end{array}


where the decision variable is the vector of :math:`N` switches :math:`\tau = \begin{bmatrix}\tau_1 & \dots & \tau_N\end{bmatrix}^\top\in \mathbb{R}^{N}`. Note that :math:`\tau_0` and :math:`\tau_{N+1}` are just used to simplify indexing and are not intended as optimizaiton variables. The parameters :math:`l_i` and :math:`u_i` define the limits of each switching interval where dynamics :math:`\dot{x}(t) = f_i(x(t))` are active.

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

  p = createsto(x0, A)

Where :code:`x0` is the initial state vector :math:`x_0` and :code:`A` is the 3-dimensional matrix defining the dynamics.


Noninear Dynamics
-------------------




Optional Arguments
---------------------
There are many additional keyword arguments which can be passed to customize the optimization problem.

+--------------------------+-------------------------------------+----------------------------------------------------+
|Parameter                 | Description                         | Default value                                      |
+==========================+=====================================+====================================================+
|:code:`t0`                | Initial Time :math:`t_0`            | :code:`0.0`                                        |
+--------------------------+-------------------------------------+----------------------------------------------------+
|:code:`tf`                | Final Time :math:`t_0`              | :code:`1.0`                                        |
+--------------------------+-------------------------------------+----------------------------------------------------+
|:code:`Q`                 | Cost matrix :math:`Q`               | :code:`eye(n)`                                     |
+--------------------------+-------------------------------------+----------------------------------------------------+
|:code:`l`                 | Vector of lower bounds :math:`l_i`  | :code:`zeros(N+1)`                                 |
+--------------------------+-------------------------------------+----------------------------------------------------+
|:code:`u`                 | Vector of lower bounds :math:`u_i`  | :code:`tf*ones(N+1)`                               |
+--------------------------+-------------------------------------+----------------------------------------------------+
|:code:`tau0ws`            | Warm starting initial solution      | Equally spaced between :code:`t0` and :code:`tf`   |
+--------------------------+-------------------------------------+----------------------------------------------------+
|:code:`solver`            | MathProgbase.jl solver              | :code:`IpoptSolver()`                              |
+--------------------------+-------------------------------------+----------------------------------------------------+

Note that any NLP solver supported by `JuliaOpt <http://www.juliaopt.org/>`_ may be used through `MathProgBase.jl <https://github.com/JuliaOpt/MathProgBase.jl/>`_ interface. For instance, in order to use `KNITRO <https://github.com/JuliaOpt/KNITRO.jl/>`_ solver with the linear example we can easily create the problem as

::

  using KNITRO
  p = createsto(x0, A, solver = KnitroSolver())


Problem Solution
======================

Once the problem is defined, it can be solved by simply running

::

  solve!(p)

The optimizer and the optimal cost function can be obtained as follows:
::

  tauopt = gettau(p)
  objval = getobjval(p)

We can get the execution time (including function calls) and the status of the solver by executing:

::

  stat = getstat(p)
  soltime = getsoltime(p)
