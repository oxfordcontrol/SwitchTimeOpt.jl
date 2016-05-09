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

Given a nonlinear system defined by dynamics

.. math::

  \dot{x}(t) = f(x(t), u(t))

where :math:`u(t)` the input vector assuming integer values :math:`u_i` between switching instants

.. math::

  u(t) = u_i \quad t\in [\tau_i, \tau_{i+1}),

we can define our switched nonlinear system as

.. math::

  \dot{x}(t) = f_i(x(t)) = f(x(t), u_i) \quad t\in [\tau_i, \tau_{i+1}).

To create the optimizaiton problem we beed to define the nonlinear dynamics by means of an additional function

::

  function nldyn(x, ui)
    ...
  end

returning the vector of states derivatives. The variable :code:`x` is the state :math:`x(t)` and :code:`ui` is the input vector :math:`u_i`. Moreover, we need to define the jacobian of the switched dynamics with respect to the system states

.. math::

  J_{f_i}(x(t)) = \frac{\partial f_i (x(t))}{\partial x}

by means of an additional function

::

  function jac_nldyn(x, ui)
    ...
  end


Note that the function :code:`jac_nldyn` returns a matrix having in each row the gradient of every component of the function :math:`f_i(x(t))` with respect to each state component. Last  element necessary to construct the matrix :code:`U` having a column each integer input vector :code:`ui`, i.e. :code:`U[:, i] = ui`. Then, we can define the switching time optimization problem as:

::

  p = createsto(x0, nldyn, jac_nldyn, U)


.. note::
  The nonliner switched system optimization operates by introducing additional linearization points between the switching intervals. To vary the number of linearization points per interval, it is just necessary to add an extra argument to the previous function call as follows:
  ::

    p = createsto(x0, nldyn, jac_nldyn, U, nlinpts)

  where :code:`nlinpts` defines the number of linearization points.

Optional Arguments
---------------------
There are many additional keyword arguments that can be be passed to the :code:`createsto(...)` function to customize the optimization problem.

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
