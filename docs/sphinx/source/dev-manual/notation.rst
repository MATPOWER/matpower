.. _sec_notation:

Notation
========

This section introduces and summarizes the mathematical notation used throughout this manual.

This notation is consistent with what was used in the MP-Element technical note, |TN5| [TN5]_ where you can find more detail.

**Styles**

.. list-table::
   :widths: 17 83
   :class: longtable

   * - :math:`x, \theta`
     - real scalars
   * - :math:`\cscal{x}, \cscal{\uptheta}`
     - complex scalars
   * - :math:`\rvec{x}, \rvec{\theta}`
     - real vectors
   * - :math:`\cvec{x}, \cvecG{\uptheta}`
     - complex vectors
   * - :math:`\rmat{X}, \rmatG{\Theta}`
     - real matrices
   * - :math:`\cmat{X}, \cmatG{\Theta}`
     - complex matrices
   * - :math:`x, \cscal{x}, \rvec{x}, \cvec{x}, \rmat{X}, \cmat{X}`
     - variables, functions
   * - :math:`\param{x}, \param{\cscal{x}}, \param{\rvec{x}}, \param{\cvec{x}}, \param{\rmat{X}}, \param{\cmat{X}}`
     - constants, parameters [#]_
   * - :math:`\hat{\rvec{x}}, \hat{\cvec{x}}, \hat{\rmat{X}}, \hat{\cmat{X}}`
     - selected rows of interest of :math:`\rvec{x}, \cvec{x}, \rmat{X}, \cmat{X}`, respectively [#]_

**Operators**

.. list-table::
   :widths: 17 83
   :class: longtable

   * - :math:`\diag{\cvec{a}}`
     - diagonal matrix with vector :math:`\cvec{a}` on the diagonal
   * - :math:`\trans{\cmat{A}}`
     - (non-conjugate) transpose of matrix :math:`\cmat{A}`
   * - :math:`\conj{\cscal{a}}`, :math:`\conj{\cvec{a}}`, :math:`\conj{\cmat{A}}`
     - complex conjugate of :math:`\cscal{a}`, :math:`\cvec{a}`, and :math:`\cmat{A}`, respectively
   * - :math:`\Re\{{\cvec{a}\}}`, :math:`\Im\{{\cvec{a}\}}`
     - real and imaginary parts of :math:`\cvec{a}`, respectively
   * - :math:`\cvec{a}^{n}`
     - element-wise exponent [#]_ for vector :math:`\cvec{a}`
   * - :math:`\cmat{A}^{n}`
     - matrix exponent [3]_ for matrix :math:`\cmat{A}`
   * - :math:`a^{\rvec{b}}`, :math:`a^\rmat{B}`
     - element-wise exponent [3]_ for vector :math:`\rvec{b}` and matrix :math:`\rmat{B}`, respectively
   * - :math:`\f(\x), \F(\x)`
     - scalar, vector functions of :math:`\x`, respectively
   * - :math:`\f_\x, \F_\x`
     - transpose of gradient of :math:`\f`, Jacobian of :math:`\F`, respectively, w.r.t. :math:`\x`
   * - :math:`\f_{\x\x}, \F_{\x\x}(\lam)`
     - Hessian of :math:`\f`, Jacobian of :math:`\trans{\F_\x} \lam`, respectively, w.r.t. :math:`\x`

**Constants and Dimensions**

.. list-table::
   :widths: 17 83
   :class: longtable

   * - :math:`e, j`
     - constants, :math:`e` is base of natural log (:math:`\approx 2.71828`), :math:`j` is :math:`\sqrt{-1}`
   * - :math:`n_k, n_n, n_p, n_p^k`
     - number of elements, nodes, ports, ports for element :math:`k`, respectively
   * - :math:`n_\X, n_\V, n_\Z`
     - dimension of vector :math:`\X`, :math:`\V`, :math:`\Z`, respectively.
   * - :math:`\ones{n}, \Id{n}`
     - :math:`n \times 1` vector of all ones, :math:`n \times n` identity matrix
   * - :math:`\zeros`
     - appropriately-sized vector or matrix of all zeros

**Variables**

.. list-table::
   :widths: 17 83
   :class: longtable

   * - :math:`\vvi{i}`
     - complex voltage at node/port :math:`i`
   * - :math:`\vri{i}, \vii{i}`
     - real and imaginary parts of voltage at node/port :math:`i`, :math:`\vvi{i} = \vri{i} + j \vii{i}`
   * - :math:`\vmi{i}, \vai{i}`
     - voltage magnitude and angle at node/port :math:`i`, :math:`\vvi{i} = \vmi{i} e^{j \vai{i}}`
   * - :math:`\V`
     - column vector of complex voltages :math:`\vvi{i}`
   * - :math:`\E`
     - column vector :math:`\V` with elements scaled to unit magnitude, :math:`\E = e^{j \Va}`
   * - :math:`\Vr, \Vi`
     - column vectors of real (:math:`\vri{i}`) and imaginary (:math:`\vii{i}`) parts of voltage, respectively, :math:`\V = \Vr + j \Vi`
   * - :math:`\Vm, \Va`
     - column vectors of voltage magnitudes :math:`\vmi{i}` and angles :math:`\vai{i}`, respectively, :math:`\V = \dVm \E = \dVm e^{j \Va}`
   * - :math:`\inV`
     - column vector of inverse of complex voltages :math:`\frac{1}{\vvi{i}}`, :math:`\inV = \V^{-1}`
   * - :math:`\z`
     - column vector of real non-voltage state variables :math:`z_i`
   * - :math:`\Z`
     - column vector of complex non-voltage state variables :math:`\cscal{z}_i`
   * - :math:`\Zr, \Zi`
     - column vectors of real and imaginary parts of :math:`\Z = \Zr + j \Zi`

**Parameters**

.. list-table::
   :widths: 17 83
   :class: longtable

   * - :math:`\J_\kk`
     - matrix formed by taking selected rows, indexed by vector :math:`\kk`, from an identity matrix [#]_
   * - :math:`\YY`
     - AC model admittance matrix
   * - :math:`\LL`
     - linear coefficient (of :math:`\Z`) for affine complex current injections
   * - :math:`\iv`
     - vector of constant complex current injections
   * - :math:`\MM`
     - linear coefficient (of :math:`\V`) for affine complex power injections
   * - :math:`\NN`
     - linear coefficient (of :math:`\Z`) for affine complex power injections
   * - :math:`\sv`
     - vector of constant complex power injections
   * - :math:`\BB`
     - DC model susceptance matrix
   * - :math:`\KK`
     - linear coefficient (of :math:`\z`) for affine active power injections
   * - :math:`\pv`
     - vector of constant active power injections
   * - :math:`\CC`
     - element-node incidence matrix for a given port
   * - :math:`\DD`
     - element-variable incidence matrix for a given state variable
   * - :math:`\Aa`
     - combined incidence matrix :math:`\Aa = \left[\begin{array}{ccc}\CC & \zeros \\ \zeros & \DD \end{array}\right]`

.. [#] Constants and parameters are underlined, with the following exceptions: constants :math:`e` and :math:`j`, :math:`p`, :math:`q`, :math:`m` and :math:`n` when used as dimensions, and :math:`i`, :math:`j`, and :math:`k` as indices.

.. [#] Obtained by multiplying by matrix :math:`\J` or :math:`\J_\kk`.

.. [#] Superscripts may also be used as indices, indicated by context.

.. [#] Often used simply as :math:`\J` without the subscript.

..
    Careful the 3rd footnote above is explicitly numbered as [3]_ in two
    references above (to avoid repeating the footnote itself).

