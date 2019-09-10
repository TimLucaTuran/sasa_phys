SASA - Semi Analytic Stacking Algorithm
============================================================

This Software allows you to calculate the optical behavior of stacked materials.
It requires you to already know the Jonas-Matrices of the complex layers in your
Stack and works out their interactions. Calculations are based on 4 by 4 composite
Jonas matrices and the Starproduct between them.

.. math::
   S =
   \left(
   \begin{matrix}
       T_f & R_f \\
       R_b & T_b \\
  \end{matrix}
  \right)

Symmetry operations can be applied directly to these Matrices 


.. toctree::
   :maxdepth: 2

   stack
   operations
   starproduct


References
^^^^^^^^^^

[1] J. Sperrhake, M. Decker, M. Falkner, S. Fasold, T. Kaiser, I. Staude, T. Pertsch,
    "Analyzing the polarization response of a chiral metasurface stack by semi-analytic modeling",
    Optics Express 1246, 2019
[2] C. Menzel, J. Sperrhake, T. Pertsch,
    "Efficient treatment of stacked metasurfaces for optimizing and enhancing the range of accessible optical functionalities",
    Physical Review A 93, 2016
[3] J. Sperrhake, T. Kaiser, M. Falkner, S. Fasold, T. Pertsch,
    "Interaction of reflection paths of light in metasurfaces stacks",
