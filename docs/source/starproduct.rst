The Starproduct
===============
The Starproduct between two S-matrices is defined as follows:

.. math::
   \left( \begin{array}{cc}
   r_1 &u_1\\
   w_1 & s_1 \\
   \end{array}\right) \star
   \left( \begin{array}{cc}
   r_2 &u_2\\
   w_2 & s_2 \\
   \end{array}\right) = \left( \begin{array}{cc}
   r_2(1-u_1 w_2)^{-1}r_1 &u_2 + r_2  u_1 (a-w_2  u_1)^{-1}s_2\\
   w_1 + s_1 w_2 (1- u_1 w_2)^{-1} r_1  & s_1 (1-w_2 \cdot u_1)^{-1}s_2 \\
   \end{array}\right)

. The following functions just apply this definition once analytically and once by usage of a geometric matrix series. Thats possible, because

.. math::
   (1-u_1 w_2)^{-1}



could be easily written down in a series.

.. automodule:: star_product
    :members:
