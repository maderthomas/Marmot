.. _bondsliplaws:

Bond-Slip Laws
==============

Bond-slip laws relate the slip of a reinforcement bar relative to its host to the bond stress on the bar surface,
in the local frame of the bar (tangential, normal 1, normal 2). They derive from ``MarmotBondSlipLaw``, live as
modules in the materials category and register in the ``MarmotBondSlipLawFactory``, from which the
:ref:`embedded bond element <embeddedbondelement>` creates them by name.

``LINEARELASTICBONDSLIP``
-------------------------

Properties ``[K_t, K_n]``: :math:`\tau_t = K_t s_t`, :math:`\sigma_n = K_n s_n`. Large stiffnesses give a penalty
formulation of perfect bond.

``MODELCODE2010BONDSLIP``
-------------------------

Properties ``[tauMax, s1, s2, s3, alpha, tauF, K0, Kn]``. Monotonic envelope of the fib Model Code 2010:

.. math::

   \tau_{\mathrm{env}}(s) = \begin{cases}
     K_0\, s & 0 \le s \le s_0 \\
     \tau_{\max} (s / s_1)^{\alpha} & s_0 < s \le s_1 \\
     \tau_{\max} & s_1 < s \le s_2 \\
     \tau_{\max} - (\tau_{\max} - \tau_f) \frac{s - s_2}{s_3 - s_2} & s_2 < s \le s_3 \\
     \tau_f & s_3 < s
   \end{cases}

The linear branch regularizes the infinite initial slope of the power law. Unloading and reloading are elastic with
:math:`K_0`, bounded by :math:`|\tau_t| \le \tau_{\mathrm{env}}(\kappa)`, :math:`\kappa = \max |s|`, beyond which the
slip is frictional (state variables ``plastic slip``, ``max slip``). The normal response is linear with :math:`K_n`.

Implementation
--------------

.. doxygenclass:: MarmotBondSlipLaw
   :project: Marmot
   :members:

.. doxygenclass:: Marmot::Materials::ModelCode2010BondSlip
   :project: Marmot
   :members:
