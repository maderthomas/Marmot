.. _beamelement:

Beam Element
============

Theory
------

A two-node Euler-Bernoulli beam (frame) element, geometrically linear, in two (``B23``: nodal dofs
:math:`u_x, u_y, \theta_z`) or three (``B33``: :math:`u_x, u_y, u_z, \theta_x, \theta_y, \theta_z`) spatial
dimensions. The node fields are ``displacement`` and ``rotation`` (ordered node by node). The names follow the
Abaqus convention for cubic Euler-Bernoulli beams.

In the local frame :math:`(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3)` (:math:`\mathbf{e}_1` from node 1 to node 2,
section coordinates :math:`y, z`), the axial displacement :math:`u` and the twist :math:`\varphi` are linear, the
transverse displacements :math:`v, w` are cubic Hermite polynomials of the nodal displacements and rotations
:math:`\theta_z = v'`, :math:`\theta_y = -w'`. No shear deformation, hence no shear locking; for a linear elastic
material, the element reproduces the exact nodal values of a beam loaded at its nodes (e.g., a cantilever with a
tip force or moment with one element). The generalized strains are

.. math::

   \varepsilon_0 = u', \qquad \kappa_y = \theta_y' = -w'', \qquad \kappa_z = \theta_z' = v'', \qquad
   \chi = \varphi' \qquad \text{(2D: } \varepsilon_0, \ \kappa = \kappa_z \text{)} .

**Fiber section.** The section is integrated with fibers :math:`(y_f, z_f, A_f)` of the element's hypoelastic
material, at 3 Gauss points along the axis. A fiber has the strains

.. math::

   \varepsilon_{11} = \varepsilon_0 + z_f \kappa_y - y_f \kappa_z, \qquad
   \gamma_{12} = - s\, z_f \chi, \qquad \gamma_{13} = s\, y_f \chi ,

and its material is evaluated with :math:`\sigma_{22} = \sigma_{33} = \sigma_{23} = 0` (Newton iteration on the
free strains, condensed tangent; a cutback is requested if it does not converge). The section forces
:math:`(N, M_y, M_z, T)` are the work conjugates of :math:`(\varepsilon_0, \kappa_y, \kappa_z, \chi)`. The fiber
coordinates are scaled such that the fibers reproduce :math:`A, I_y, I_z` exactly, and
:math:`s = \sqrt{J / \sum_f (y_f^2 + z_f^2) A_f}` gives the torsional stiffness :math:`GJ`. For a linear elastic
material: :math:`N = EA\varepsilon_0`, :math:`M_y = EI_y\kappa_y`, :math:`M_z = EI_z\kappa_z`,
:math:`T = GJ\chi` exactly, for every profile:

* 0, generic: 2 fibers at :math:`y = \pm\sqrt{I/A}` (2D) or 4 fibers at :math:`(\pm\sqrt{I_z/A}, \pm\sqrt{I_y/A})`
  (3D) -- exact for elastic materials, an idealized sandwich section for inelastic ones.
* 1, rectangle: :math:`n` layers (2D) or :math:`n \times n` fibers (3D). With the scaling, the plastic moment of the
  fibers is :math:`f_y b h^2/4 \cdot n / \sqrt{n^2 - 1}`.
* 2, circle: :math:`\max(1, n/2)` rings of equal area with :math:`2n` (rounded up to a multiple of 4) sectors.

Element properties:

* ``B23``: :math:`[A, I, \text{profile}, n]` with :math:`I = \int y^2 \mathrm{d}A`;
* ``B33``: :math:`[A, I_y, I_z, J, v_x, v_y, v_z, \text{profile}, n]` with :math:`I_y = \int z^2 \mathrm{d}A`,
  :math:`I_z = \int y^2 \mathrm{d}A` and the orientation vector :math:`\mathbf{v}`, whose part normal to the axis
  is :math:`\mathbf{e}_2` (:math:`\mathbf{e}_3 = \mathbf{e}_1 \times \mathbf{e}_2`).

profile and :math:`n` are optional (generic, 8). A body force (per volume) gives the consistent nodal forces and
moments of the line load :math:`\mathbf{b} A`.

Quadrature point states: ``normal force``, ``bending moment`` (2D) / ``bending moments`` and ``torque`` (3D),
``axial strain``, ``curvature`` (2D) / ``curvatures`` and ``twist rate`` (3D), ``section forces`` and
``section strains`` (all components), ``elastic energy``, ``dissipation``; the state of a fiber :math:`k` as
``fiber k stress`` (:math:`\sigma_{11}, \sigma_{12}, \sigma_{13}`) or ``fiber k <material state>``.

Limitations: geometrically linear (small rotations; a finite rigid rotation strains the element), no shear
deformation, no distributed surface loads, no inertia. For nonlinear materials, the linear axial and cubic
transverse interpolation cannot represent a shifting neutral axis within an element (as in any displacement-based
fiber beam): refine the mesh where the beam yields.

Implementation
--------------

.. doxygenclass:: Marmot::Elements::BeamElement
   :project: Marmot
   :members:
