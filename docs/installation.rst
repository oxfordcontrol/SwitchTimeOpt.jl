===========================================
Installation
===========================================
You can easily install the package by running

::

  Pkg.clone("git://github.com/bstellato/SwitchTimeOpt.jl.git")

This does not install any nonlinear solvers. If you don’t have a nonlinear solver installed already, you will want to install a solver such as `Ipopt <https://github.com/JuliaOpt/Ipopt.jl/>`_ by running:

::

  Pkg.add("Ipopt")
