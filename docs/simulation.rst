===============================
Simulation
===============================

The system can be simulated with the obtained switching times by running

::

  x, xsw, optval, t = simulate(m)
  x, xsw, optval, t = simulate(m, tau)      # Specify switching time vector
  x, xsw, optval, t = simulate(m, t)        # Specify time vector
  x, xsw, optval, t = simulate(m, tau, t)   # Specify switching time and time vectors


The outputs of the simulation are

  * :code:`x` State trajectory. Each :math:`x(t)` can be obtained as :code:`x[:, i]`
  * :code:`xsw` States at each switching time. Each :math:`x(\tau_i)` can be obtained as :code:`xsw[:, i]`
  * :code:`optval` Simulated optimal value function
  * :code:`t` Time vector during the simulation


Additional Functions for Nonlinear Dynamics
--------------------------------------------

In case of nonlinear dynamics it is possible to simulate the system linearized at the linearization points obtained after the optimization

::

  x, xsw, optval, t = simulatelinearized(m)
  x, xsw, optval, t = simulatelinearized(m, t)  # Specify time vector


In addition, it is possible to easily obtain the input vector at each time instant by running

::

  u, t = simulateinput(m)
  u, t = simulateinput(m)  # Specify time vector

Each vector :math:`u(t)` can be obtained by slicing the output :code:`u[:, i]`.
