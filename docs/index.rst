.. SwitchTimeOpt.jl documentation master file, created by
   sphinx-quickstart on Mon May  9 09:22:48 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SwitchTimeOpt.jl - Switching Time Optimization in Julia
=======================================================

**SwitchTimeOpt.jl** is a `Julia <https://github.com/JuliaLang/julia/>`_ package to easily define and efficiently solve switching time optimization (STO) problems for linear and nonlinear systems. SwitchTimeOpt.jl supports a wide variety of nonlinear solvers through `MathProgBase.jl <https://github.com/JuliaOpt/MathProgBase.jl/>`_ interface such as `Ipopt <https://github.com/JuliaOpt/Ipopt.jl/>`_, `KNITRO <https://github.com/JuliaOpt/KNITRO.jl/>`_, `NLopt <https://github.com/JuliaOpt/NLopt.jl/>`_.

.. toctree::
   :maxdepth: 2

   Installation <installation>
   Quick Example <quick_example>
   Optimization <optimization>
   Simulation <simulation>


Citing this package
--------------------
If you use SwitchTimeOpt.jl for published work, we encourage you to cite the following `paper <http://arxiv.org/abs/1608.08597/>`_:

.. code-block:: latex

  @article{2016arXiv160808597S,
    author = {{Stellato}, B. and {Ober-Bl{\"o}baum}, S. and {Goulart}, P.~J.},
    title = "{Second-Order Switching Time Optimization for Switched Dynamical Systems}",
    journal = {ArXiv e-prints},
    archivePrefix = "arXiv",
    eprint = {1608.08597},
    primaryClass = "math.OC",
    keywords = {Mathematics - Optimization and Control},
    year = 2016,
    month = aug
  }
..
.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
