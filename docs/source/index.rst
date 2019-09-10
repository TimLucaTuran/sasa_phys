SASA - Semi Analytic Stacking Algorithm
============================================================

This Software allows you to calculate the optical behavior of stacked materials.
It requires you to already know the Jonas-Matrices of the complex layers in your
Stack and works out their interactions. Calculations are based on Lx4x4 composite
Jonas matrices, called S-matrices, and the Starproduct between them. L represents
the wavelengths were you wish to calculate the behavior.

.. math::
   S =
   \left(
   \begin{matrix}
       T_f & R_f \\
       R_b & T_b \\
  \end{matrix}
  \right)

:math:`T_f:` Transmission Jonas matrix for light coming from the front

:math:`R_b:` Reflection for the back

Usage
^^^^^
The exact usage is described in example_usage.py. In general you have to define
multiple Layer-Objects::

    l1 = MetaLayer(s_mat, cladding, substrate)
    l2 = NonMetaLayer(n_vec, cladding, substrate)

These can be Meta-Layers where you need to provide a Lx4x4 S-matrix or Non-Meta-Layers
where you need to provide a vector of refractive indices's at the desired
wavelengths. Then you pass the layers to a stack object and build your result::

    s = Stack([l1,l2,...], wavelengths, cladding, substrate)
    result = s.build()

In the case of layers cladding and substrate represent the environment in which
s_mat or n_vec were measured. For Stack its what materials are blow/on-top
of the Stack.


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
