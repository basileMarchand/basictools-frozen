**********************
Solving linear systems
**********************

Finite element computation make use of a linear solvers to solve assembled linear system of equations.
The libraries Numpy [#numpyurl]_, Scipy [#scipyurl]_ and Eigen [#eigenurl]_ offer capabilities to solve sparse linear systems by different methods.
BasicTools propose a single class to wrap all this solvers into one unique simple interface.

Treatment of kinematic constraints
##################################

.. math::
    \begin{equation}
    \begin{matrix}
    \text{minimize } \UFE \text{ for:} & \frac{1}{2} \UFE^T \KFE \UFE - \UFE^T \FFE \\
    \end{matrix}
    \end{equation}
    :label: KUF

.. math::
    \begin{equation}
    \begin{matrix}
    \text{subject to} &\CFE \UFE = \GFE \\
    \end{matrix}
    \end{equation}
    :label: CUG

Equation :eq:`KUF` represent the assembled finite element system to be solved under the kinematic constraints defined the equation :eq:`CUG`.
This problem can be solved by four different approaches.

The constrains represented by the equation :eq:`CUG` are created from the Dirichlet like boundary conditions (like RBE2 and RBE3).
In some cases, this system can have linearly dependent equations, for example when a user specified twice some boundary condition in the same entity (e.g. the shared edge in the case of :math:`u_x=0` on 2 adjacent faces).

To be able to treat correctly the constraint defined by :eq:`CUG`, we must clean the linearly dependent equations, and find (in some cases) a reduced number of dof to eliminate .

Let extract from :math:`\CFE` only the non zero columns  (:math:`\CFE^*`) and define the augmented matrix of :math:`\CPFE` as :

.. math::
    \begin{equation}
    \CPFE :=  \left[ \begin{array}{c|c}
    \CFE^* &  \GFE
    \end{array} \right]
    \end{equation}

The :math:`\CPFE^T`  matrix can be then decompose using a QR decomposition to extract a the suitable (orthonormal) base.
This is done using the routine of EIGEN SpQR.
SpQR is only available if the c++ interface was compiled, a pure Python (using `scipy.linalg.qr`) backup algorithm is available.
This backup approach uses dense linear and can be expensive (in terms of memory and cpu time).
This algorithm is capable of finding the rank of the matrix :math:`\CPFE^T` (using a prescribed tolerance).

.. math::
    \begin{equation}
    \CPFE \mathbf{P} = \mathbf{Q} \mathbf{R}
    \end{equation}

With,
:math:`\mathbf{P}` is the column permutation,
:math:`\mathbf{Q}` is the orthogonal matrix represented as Householder reflectors (in the case using Eigen).
:math:`\mathbf{R}` is the sparse triangular factor.

In many cases the constraints are on single nodes or only in a reduced number of nodes.
This make the matrix :math:`\CPFE^T` very sparse, and we can compute the :math:`\mathbf{QR}` decomposition block by block.
For this the adjacency matrix of the matrix :math:`\CFE^*` is computed (:math:`abs(\CFE^*)abs(\CFE^*)^T`).
A not zero :math:`(i,j)` term of the adjacency matrix means the row :math:`i` and :math:`j` of :math:`\CPFE` share at least one degree of freedom in common.
Then we can calculate the connected components of the adjacency matrix and work on each block independently.
As each block is treated (computation of the QR decomposition), the rank is identified to remove redundant equations.
At the same time for each of the sub-matrices a subset of dofs equal to the rank is stored to created the a list of slave dofs for later use.

For example:

.. _two-springs:

.. figure:: images/Spring1.svg
    :width: 500px
    :align: center
    :alt: 2 Springs
    :figclass: align-center

    Two springs system with 4 degrees of freedom, A) initial state B) constrained solution.


For the system of 2 independent springs (4 degrees of freedom) represented in figure, the
tangent matrix and the right hand side member are,

.. math::
    \begin{equation}
   \KFE =
   \begin{bmatrix}
    1000 & -1000 &    0  &    0 \\
   -1000 &  1000 &    0  &    0 \\
       0 &     0 & 1000  &-1000 \\
       0 &     0 &-1000  & 1000
   \end{bmatrix}, \FFE =
   \begin{bmatrix}
   0 \\
   0 \\
   0 \\
   0
   \end{bmatrix}.
   \end{equation}


The constraint imposed to the problem are: 1) blockage of the first dof (:math:`u_0` to the value 0), 2) kinematic relation between :math:`u_1` and :math:`u_2`;  :math:`u_2 - u_1 = 1`, 3) prescribed solution on :math:`u_3` equal to 3. To demonstrate the treatment, each constraint was added 2 time to the matrix :math:`\CFE` :

.. math::
    \begin{equation}
    \CFE =
    \begin{bmatrix}
    1 & 0 &0 &0 \\
    1 & 0 &0 &0 \\
    0 &-1 &1 &0 \\
    0 &-2 &2 &0 \\
    0 & 0 &0 &1\\
    0 & 0 &0 &1
    \end{bmatrix}, \GFE =
    \begin{bmatrix}
    0 \\
    0 \\
    1 \\
    2 \\
    3 \\
    3
    \end{bmatrix}
    \end{equation}

In this case all degrees of freedom are present in matrix :math:`\CFE`, this means :math:`\CFE^* = \CFE`.
We start by building the adjacency matrix :math:`\CFE\CFE^T`, we are interested in only the non zero values so we really calculate :math:`abs(\CFE)abs(\CFE^T) != 0`.

.. math::
    \begin{equation}
    \mathbf{AM} := abs(\CFE)abs(\CFE^T)!=0  \to
    \begin{bmatrix}
    1 & 1 &0 &0 &0 &0 \\
    1 & 1 &0 &0 &0 &0 \\
    0 & 0 &1 &1 &0 &0\\
    0 & 0 &1 &1 &0 &0\\
    0 &0 &0 & 0 &1 &1\\
    0 &0 &0 & 0 &1 &1\\
    \end{bmatrix}
    \end{equation}

The computation of the connected components can be expressed by the solution vector:

.. math::
    \begin{equation}
    Connected Component Of(\mathbf{AM}) =
    \begin{bmatrix}
    0 & 0& 1 & 1 &2 &2
    \end{bmatrix}
    \end{equation}

We can see 3 connected components (rows of A).
Each of this component can be treated independent (they are orthogonal).
For the first component (line 1 and 2), it is evident that the lines are linearly dependent and we have only one \emph{real} constraint.
This can be detected by the SpQr routine and only the relevant base vector are kept (we also normalize the vectors).
We also store the rank first dofs of each component to define the slave indices.
At the end we obtain 1) a vector of slave dofs of the system :math:`[0,1,3]`, and a clean constants system defined by :

.. math::
    \begin{equation}
    \MFE = \begin{bmatrix}
    1 &0 &0 &0 \\
    0 &-\sqrt{\frac{1}{2}} &\sqrt{\frac{1}{2}} &0 \\
    0 &0 &0 &1 \\
    \end{bmatrix}, \VFE =
    \begin{bmatrix}
    0 \\
    \sqrt{\frac{1}{2}} \\
    -3 \\
    \end{bmatrix}
    \end{equation}


Penalization
############

The penalization approach is the simplest way to solve the contained system (in terms of programming).
The idea of this method is to add to the original system to be solve a penalization term scaled by a really large value (:math:`\alpha_p = 10^8`).
The perturbed system becomes :

.. math::
    \begin{equation}
    \KFE_p =  \KFE + maxdiag(\KFE)*\alpha_p*(\MFE^T \MFE) \\
    \end{equation}


.. math::
    \begin{equation}
    \FFE_p =  \FFE + maxdiag(\KFE)*\alpha_p*(\MFE^T \VFE) \\
    \end{equation}


Because the modified system has the same number of dofs and no base transformation was done the solution
of this system gives directly an original system (equations :eq:`KUF` and :eq:`CUG`) solution approximation


.. note: This technique has the disadvantage of altering heavily the condition number of the system to be solved, making it (in some case) impossible to solve using an iterative solver.


Lagrange Multipliers
####################

This method impose the constraints  by adding one extra dof for each row of equation :eq:`CUG`.
This yield to the following modified system:

.. math::
    \begin{equation}
    \KFE_{\lambda} = \begin{bmatrix}
    \KFE & \MFE \\
    \MFE & 0  \\
    \end{bmatrix},
    \FFE_{\lambda} = \begin{bmatrix}
    \FFE \\
    \VFE
    \end{bmatrix},
    \end{equation}

The resolution of this linear system gives the solution vector:

.. math::
    \begin{equation}
    \UFE_{\lambda} = \begin{bmatrix}
    \UFE \\
    \lambda
    \end{bmatrix},
    \end{equation}

The solution to the original system is then  contained in the first part of the vector :math:`\UFE_{\lambda}`.

In the case that we have an initial a solution :math:`\UFE^0`, and we want to solve the system :math:`\KFE_{\lambda}\UFE_{\lambda} = \FFE_{\lambda}` using this initial guess, we cand define the initial guess :math:`\UFE_{\lambda}`:

.. math::
    \begin{equation}
    \UFE_{\lambda} = \begin{bmatrix}
    \UFE^0 \\
    \MFE(\FFE-\KFE)
    \end{bmatrix}
    \end{equation}


.. note:: This technique has the disadvantage of losing the positive define properties of the original system (zeros on the diagonal of the operator :math:`\KFE_{\lambda}`)


Substitution
############

Let decompose the :math:`\MFE` (using the indices colleted durint the QR decomposition) into slave (:math:`s`) and master (:math:`m`) dofs :

.. math::
    \begin{equation}
    \begin{matrix}
    \MFE = \left[  \begin{array}{c|c}
    \MFE_s & \MFE_m
    \end{array}
    \right]
    \end{matrix}
    \end{equation}

The the original system becomes:

.. math::
    \begin{equation}
    \begin{bmatrix}%[l]
    \KFE_s & \KFE_{s,m}  \\
    \KFE_{m,s}   & \KFE_{m}
    \end{bmatrix}
    \begin{bmatrix}
    \UFE_{s} \\
    \UFE_{m}
    \end{bmatrix}
    =
    \begin{bmatrix}
    \FFE_{s} \\
    \FFE_{m}
    \end{bmatrix}
    \end{equation}

Then the constraint equation :eq:`CUG` can be rewritten to solve :math:`\UFE_s` :

.. math::
    \begin{equation}
    \begin{array}{l}
    %\MFE_{s}\UFE_s + \MFE_m \UFE_m = \VFE  \\
    \UFE_s  = \MFE_{s}^{-1} (\VFE - \MFE_m \UFE_m  )
    \end{array}
    \end{equation}

Now we can rewrite the :math:`\UFE` in the form of :

.. math::
    \begin{equation}
    \UFE = \begin{bmatrix}
    \UFE_s \\
    \UFE_m \\
    \end{bmatrix}=
    \begin{bmatrix}
    \MFE_{s}^{-1} (\VFE - \MFE_m \UFE_m  ) \\
    \UFE_m \\
    \end{bmatrix}=
    \begin{bmatrix}
    \MFE_{s}^{-1} (\VFE  ) \\
    0 \\
    \end{bmatrix}+
    \begin{bmatrix}
    -\MFE_{s}^{-1} ( \MFE_m   ) \\
    I \\
    \end{bmatrix}\UFE_m
    \end{equation}
    :label: Ulambda

By defying:

.. math::
    \begin{equation}
    \XFE = \begin{bmatrix}
    -\MFE_{s}^{-1} ( \MFE_m   ) \\
    I \\
    \end{bmatrix}, \DFE= \begin{bmatrix}
    \MFE_{s}^{-1} (\VFE  ) \\
    0 \\
    \end{bmatrix}
    \end{equation}

And injection this expression on the original system we obtain:

.. math::
    \begin{equation}
    \XFE^T \KFE \XFE \UFE_m =
    \XFE^T \left( \FFE - \KFE \DFE  \right)
    \end{equation}

.. math::
    \begin{equation}
    \begin{matrix}
    \underbrace{
    \begin{bmatrix}
    -\MFE_{s}^{-1} ( \MFE_m   ) \\
     I \\
    \end{bmatrix}^T }_{X^T}
    \begin{bmatrix}%[l]
    \KFE_s & \KFE_{s,m}  \\
    \KFE_{m,s}   & \KFE_{m}
    \end{bmatrix}
    \underbrace{
    \begin{bmatrix}
    -\MFE_{s}^{-1} ( \MFE_m   ) \\
    I \\
    \end{bmatrix} }_{X}
    \UFE_m = \dotsb \\
    \dotsb \underbrace{
    \begin{bmatrix}
    -\MFE_{s}^{-1} ( \MFE_m   ) \\
    I \\
    \end{bmatrix}^T }_{X^T}
    \left(\begin{bmatrix}
    \FFE_{s}\\
    \FFE_{m} \\
    \end{bmatrix} - \begin{bmatrix}%[l]
    \KFE_s & \KFE_{s,m}  \\
    \KFE_{m,s}   & \KFE_{m}
    \end{bmatrix}
    \underbrace{\begin{bmatrix}
    \MFE_{s}^{-1} (\VFE  ) \\
    0 \\
    \end{bmatrix} }_D \right) \\
	\end{matrix}
    \end{equation}

The modified system involves only :math:`\UFE_m` dofs, and can be solved using a standard solver.
The solution to the original system is then calculates using equation :eq:`Ulambda`.


Ainsworth Method
################

The Ainsworth metod [#ainsworth]_ propose a general approach to solve the original system of equations without chnging the number of degrees of freedom of the system.




.. rubric:: Footnotes
.. [#numpyurl] https://numpy.org/
.. [#scipyurl] https://www.scipy.org/
.. [#eigenurl] https://eigen.tuxfamily.org/dox-devel/classEigen\_1\_1SPQR.html
.. [#ainsworth] https://doi.org/10.1016/S0045-7825(01)00236-5