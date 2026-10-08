.. _trusselement:

Truss Element
=============

Theory
------

A truss (bar) element in two or three spatial dimensions with 2 (linear) or 3 (quadratic; node order end, end,
mid) nodes. It carries only an axial force. The material sees a uniaxial stress state along the truss axis, found by
the uniaxial stress reductions of the material base classes, so any 3D material of Marmot can be used.

**Small strain** (``TR2D2``, ``TR2D3``, ``TR3D2``, ``TR3D3``; hypoelastic materials): the axial strain is
:math:`\varepsilon = \mathbf{T} \cdot \partial \mathbf{u} / \partial S` with the reference axis :math:`\mathbf{T}`
and arc length :math:`S`. The material is called through ``MarmotMaterialHypoElastic::computeUniaxialStress``.

**Finite strain** (``TR2D2FS``, ``TR2D3FS``, ``TR3D2FS``, ``TR3D3FS``; finite strain materials): the material is
evaluated in the co-rotated frame of the truss, :math:`\mathbf{F} = \mathrm{diag}(\lambda, \lambda_2, \lambda_3)`,
with the axial stretch :math:`\lambda = |\mathbf{g}| / |\mathbf{G}|` and the lateral stretches found by
``MarmotMaterialFiniteStrain::computeUniaxialStress`` such that :math:`\tau_{22} = \tau_{33} = 0`. Rigid body
rotations are exactly stress free. With the nominal stress :math:`P = \tau_{11} / \lambda` and the current axis
:math:`\mathbf{n}`,

.. math::

   \mathbf{f}_a = \int P\, \mathbf{n}\, \frac{\partial N_a}{\partial S}\, A_0\, \mathrm{d}S, \qquad
   \mathbf{K}_{ab} = \int \Big[ \frac{\mathrm{d}P}{\mathrm{d}\lambda}\, \mathbf{n} \otimes \mathbf{n}
     + \frac{P}{\lambda} (\mathbf{I} - \mathbf{n} \otimes \mathbf{n}) \Big]
     \frac{\partial N_a}{\partial S} \frac{\partial N_b}{\partial S}\, A_0\, \mathrm{d}S .

Element property: the cross section area :math:`A_0`.

Quadrature point states: ``stress`` (axial Cauchy stress), ``strain`` (small strain: axial strain; finite strain:
:math:`\ln\lambda`), ``normal force``, ``elastic energy``, ``dissipation``, ``lateral stretches``, ``kirchhoff
stress``, followed by the material state.

Implementation
--------------

.. doxygenclass:: Marmot::Elements::TrussElement
   :project: Marmot
   :members:
