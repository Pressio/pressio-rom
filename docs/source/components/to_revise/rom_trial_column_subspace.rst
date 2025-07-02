.. include:: ../mydefs.rst

.. role:: cpp(code)
   :language: cpp

Trial Subspace
==============

Header: ``<pressio/rom_subspaces.hpp>``

API
---

.. code-block:: cpp

   namespace pressio{ namespace rom{

   template<class ReducedStateType, class BasisMatrixType, class FullStateType>
   /*impl defined*/ create_trial_column_subspace(BasisMatrixType && basisMatrix,
						 FullStateType && translation,
						 bool isAffine);

   }} //end namespace

Parameters
-----------

- ``basisMatrix``: the basis matrix

- ``translation``: the affine translation (potentially a "zero" vector)

- ``isAffine``: boolean to indicate the space is affine

Constraints
-----------

Let ``basis_matrix_type = std::remove_cv_t<BasisMatrixType>`` and
``full_state_type = std::remove_cv_t<FullStateType>``,
then all of the following must hold:

- :cpp:`pressio::is_vector_eigen<ReducedStateType>::value == true`

- ``basis_matrix_type`` is a rank-2 matrix data type already supported
  in pressio, i.e. an Eigen Matrix, a Kokkos rank-2 View, a Trilinos Tpetra multivector,
  or a Tpetra block multivector

- ``full_state_type`` is an Eigen vector, a Kokkos rank-1 View, a Trilinos Tpetra vector,
  or a Tpetra block vector

..
   - :cpp:`std::is_class<BasisMatrixType>::value == true &&
     std::is_copy_constructible<basis_matrix_type>::value == true &&
     std::is_pointer<basis_matrix_type>::value == false &&
     pressio::mpl::is_std_shared_ptr<basis_matrix_type>::value == false`

   - :cpp:`std::is_class<full_state_type>::value == true &&
     std::is_copy_constructible<full_state_type>::value == true &&
     std::is_pointer<full_state_type>::value == false &&
     pressio::mpl::is_std_shared_ptr<full_state_type>::value == false`

   And let ``M`` and ``S`` be instances of ``basis_matrix_type`` and ``full_state_type``,
   then the following operations must be defined:

   - ``pressio::ops::extent(M, /*int*/)``

   - ``pressio::ops::extent(M, /*int*/)``


..
   Mandates
   --------

   - :cpp:`pressio::all_have_traits_and_same_scalar<
     reduced_state_type, full_state_type, basis_matrix_type>::value == true`

Preconditions
-------------

- ``basisMatrix`` *represents a linear basis with full column rank*

- ``translation`` must be compatible with the full space, i.e. the following
  must hold ``pressio::ops::extent(translation, 0) == pressio::ops::extent(basisMatrix, 0)``

..
   - the following operation ``auto basisMatrixClone = pressio::ops::clone(basisMatrix)``
     conforms to the requirements/semantics `explained here <ops/clone.html>`__;
     in brief, ``basisMatrixClone`` is an object independent of ``basisMatrix`` and if any
     modification is applied to ``basisMatrixClone`` (or any data it references/uses),
     these do not affect the original ``basisMatrix`` object (or any data it references/uses),
     and viceversa. The same should hold for ``translation``.


Return value and effects
------------------------

Creates an instance of a trial subspace. The return type is implementation defined,
so just use `auto`. If passing lvalue arguments, the implementation will "clone" the
arguments so new allocations will occur.
If passing rvalues (temporaries), the code will use move semantics,
so the temporary objects are moved-constructed and no new memory allocation should occur.

Postconditions
--------------

The return type is guaranteed to model
the `PossiblyAffineTrialColumnSubspace concept <rom_concepts.html>`__.


Example
-------

.. code-block:: cpp

   int main(){
     using basis_t  = Eigen::MatrixXd;
     using translation_t = Eigen::VectorXd;

     basis_t basis1(56,5);
     // fill such that columns are linearly independent
     translation_t translation1(56);
     auto space = create_trial_column_subspace<Eigen::VectorXd>(basis1, translation1, false);

     basis_t basis2(56,5);
     translation_t translation2(156); // violates precondition

     basis_t basis3(56,3);
     // assume the following:
     basis3.col(0).setConstant(1.);
     basis3.col(1).setConstant(2.);
     basis3.col(3).setConstant(4.);
     // basis3 is NOT full rank since it does NOT contain
     // linearly indepdent columns so this violates the precondition
   }






..
   .. code-block:: cpp

      namespace pressio{ namespace rom{

      template <class ReducedStateType, class BasisMatrixType, class FullStateType>
      class TrialColumnSubspace
      {
	public:
	  // type aliases
	  using reduced_state_type = ReducedStateType;
	  using basis_matrix_type  = std::remove_cv_t<BasisMatrixType>;
	  using full_state_type    = std::remove_cv_t<FullStateType>;
	  using offset_type        = full_state_type;
	  using tangent_space_basis_matrix_type = basis_matrix_type;

	  // constructors
	  TrialColumnSubspace(const basis_matrix_type & basisMatrix,                         (1)
			      const full_state_type & offset,
			      bool isAffine);

	  TrialColumnSubspace(basis_matrix_type && basisMatrix,                              (2)
			      full_state_type && offset,
			      bool isAffine);

	  TrialColumnSubspace(const basis_matrix_type & basisMatrix,                         (3)
			      full_state_type && offset,
			      bool isAffine);

	  TrialColumnSubspace(basis_matrix_type && basisMatrix,                              (4)
			      const full_state_type & offset,
			      bool isAffine);

	  // copy and copy assignment
	  TrialColumnSubspace(const TrialColumnSubspace &);                                  (5)
	  TrialColumnSubspace & operator=(const TrialColumnSubspace &) = delete;             (6)

	  // accessors
	  reduced_state_type createReducedState() const;                                     (7)

	  full_state_type createFullState() const;                                           (8)

	  void mapFromReducedState(const reduced_state_type & reducedState,                  (9)
				   full_state_type & fullState) const;

	  full_state_type createFullStateFromReducedState(const reduced_state_type &)const; (10)

	  const basis_matrix_type & viewBasis() const;                                      (11)
	  std::size_t dimension() const;                                                    (12)
	  bool isColumnSpace() const;                                                       (13)
	  bool isRowSpace() const;                                                          (14)

	  const tangent_space_basis_matrix_type & viewTangentSpaceBasis() const;            (15)
	  const offset_type & viewAffineOffset() const;                                     (16)
	  bool isAffine() const;                                                            (17)
      };

      }} //end namespace

   Constraints
   -----------

   - :cpp:`pressio::is_vector_eigen<ReducedStateType>::value == true`

   - :cpp:`std::is_copy_constructible<basis_matrix_type>::value == true &&` \
     :cpp:`std::is_pointer<basis_matrix_type>::value == false &&
     pressio::mpl::is_std_shared_ptr<basis_matrix_type>::value == false`

   - :cpp:`std::is_copy_constructible<full_state_type>::value == true &&
     std::is_pointer<full_state_type>::value == false &&
     pressio::mpl::is_std_shared_ptr<full_state_type>::value == false`

   Mandates
   --------

   - :cpp:`std::is_same<typename pressio::Traits<basis_matrix_type>::scalar_type,
     typename pressio::Traits<full_state_type>::scalar_type>::value == true`

   Member functions
   ----------------

   Constructors
   ^^^^^^^^^^^^

   :memberfunction:`(1)`: constructor

     *Parameters*:

     - ``basisMatrix``: the basis matrix

     - ``offset``: the affine translation

     - ``isAffine``: boolean to indicate the space is affine

     *Preconditions*:

     - ``basisMatrix`` *represents a linear basis with full column rank*

     - the following operation ``auto basisMatrixClone = pressio::ops::clone(basisMatrix)``
       conforms to the requirements/semantics `explained here <ops/clone.html>`__;
       in brief, ``basisMatrixClone`` is an object independent of ``basisMatrix`` and if any
       modification is applied to ``basisMatrixClone`` (or any data it references/uses),
       these do not affect the original ``basisMatrix`` object (or any data it references/uses),
       and viceversa. The same should hold for ``offset``.

     - ``offset`` must be compatible with the full space, i.e. the following
       must hold ``pressio::ops::extent(offset, 0) == pressio::ops::extent(basisMatrix, 0)``

     *Effects*:

     - calls ``pressio::ops::clone(...)`` on both ``basismatrix`` and ``offset``,
       and invokes the copy-constructor on the results to initialize
       data members, therefore making new allocations

     *Postconditions*:

     - the constructed object is immutable and owns, until its destruction,
       data members that are clones (according to ``ops::clone()``)
       of the arguments originally passed to the constructor. Therefore, this class
       has a well-defined invariant established upon construction.


   :memberfunction:`(2)`: constructs an instance of this class from rvalues

     *Parameters*:

     - same as in (1)

     *Preconditions*:

     - ``basisMatrix`` *represents a linear basis with full column rank*

     - ``offset`` must be compatible with the full space, i.e. the following
       must hold ``pressio::ops::extent(offset, 0) == pressio::ops::extent(basisMatrix, 0)``

     *Effects*:

     - tries to use move semantics by move-constructing ``basisMatrix`` and ``offset``.
       If their types implement move semantics,
       the temporary objects are moved-constructed and no new memory allocation should occur.
       If move semantics are not implemented and the move operation
       falls back to a copy, a new allocation (potentially) occurs.

     *Postconditions*:

     - the constructed object is immutable and owns, until its destruction,
       data members that were either move constructed from the
       original arguments, or are just copies. As in (1), this class
       has a clear invariant that is established upon construction and cannot be
       changed until its destruction.


   :memberfunction:`(3)`: constructor

     - behaves like (1) for ``basisMatrix`` and like (2) for ``offset``

   :memberfunction:`(4)`: constructor

     - behaves like (2) for ``basisMatrix`` and like (1) for ``offset``

   :memberfunction:`(5)`: copy constructor

     *Parameters*:

     - ``other``: another subspace object to be used as source to copy from

     *Effects*:

     - constructs the object with a "clone" of the contents in ``other``,
       with semantics similar to constructor (1)

     *Postconditions*:

     - the new object will be identical to ``other``


   :memberfunction:`(6)`: copy assignment operator

      - this is explicitly marked as deleted because it would violate immutability


   Accessors
   ^^^^^^^^^

   :memberfunction:`(7)`: create an instance of the reduced state

     *Effects*: allocates a new instance of a reduced state,
     zero initializes all values and returns it

     *Postconditions*:

     - let ``s`` be an instance of the class, and ``auto result = s.createReducedState()``,
       then the following holds: ``pressio::ops::extent(result, 0) == s.dimension();``


   memberfunction:`(8)`: create an instance of the full state

     *Effects*: allocates a new instance of a full state,
     zero initializes all values and returns it

     *Postconditions*:

     - let ``s`` be an instance of the class, ``M`` the basis of the tangent space,
       and ``auto result = s.createFullState()``, then the following holds:
       ``pressio::ops::extent(result, 0) == pressio::ops::extent(M, 0);``

   :memberfunction:`(9)`: map a reduced state to the full state

     *Parameters*:

     - ``reducedState`` reduced state instance to map from

     - ``fullState`` the full state instance to overwrite

     *Preconditions*:

     - ``reducedState`` must match the dimensionality of the subspace

     - ``fullState`` must match the dimensionality of the *full space*

     *Effects*:

     Let ``y_r`` be a reduced state instance, ``y_f`` a full state instance,
     ``s`` be an instance of the class, ``M`` and ``o`` be the basis matrix and offset stored in ``s``,
     then calling ``s.mapFromReducedState(y_r, y_f)`` means that:

     - if ``isAffine == false``, the operation performed is ``y_f = M y_r``,
       which corresponds to calling: ``::pressio::ops::product(::pressio::nontranspose(), /*one*/,
       M, y_r, /*zero*/, y_f)``

     - if ``isAffine == true``, the operation performed is ``y_f = M y_r + o``,
       which corresponds to first calling ``::pressio::ops::product(::pressio::nontranspose(), /*one*/, M,
       y_r, /*zero*/, y_f)``, followed by:
       ``::pressio::ops::update(y_f, /*one*/, o, /*one*/);``


   :memberfunction:`(10)`: creates a full state instance and maps it from a reduced state

     *Effects*: equivalent to calling (8), then (9), and returning


   ..
      :memberfunction:`(11)`: returns a reference to the basis object

	 - Note: do NOT use ``const_cast`` on the returned value to modify
	   the referenced object, as doing so will break the invariant
	   and the immutability of the class.
	   If you need a different linear subspace, create a new object.

      :memberfunction:`(12)`: returns the dimension of the subspace

	*Effects*: let ``s`` be an instance of the class and ``M`` be the basis matrix stored in ``s``,
	  then the method returns ``pressio::ops::extent(M, 1)``

      :memberfunction:`(13)`: query if this represents an affine subspace

	*Effects*: returns ``true`` is the object was constructed as an
	affine space, ``false`` otherwise


      :memberfunction:`(14)`: query if the spanning set is the columns set

	*Effects*: always returns ``true`` since by constructor this is a column subspace

      :memberfunction:`(15)`: query if the spanning set is the rows set

	*Effects*: always returns ``false`` since by constructor this is a column subspace

   Notes
   ^^^^^

   - the immutability and invariance of this class strictly
     hinges on the assumption that users do not violate const
     correctness: please do not const_cast the return of (5)

   - the design and semantics of the class functions needs a clarification:
     the intent is to ensure the class *always* has value semantics because
     it nees to work even for template arguments with view semantics.
     Suppose that ``basisMatrix`` were an instance of a ``Tpetra::MultiVector``
     or a ``Kokkos::View``, which are both classes with view semantics.
     In such case, if constructor (1) were to simply copy-construct it,
     it would lead to a shallow copy. This would
     compromise the immutability and invariant of the class,
     since any changes to the original object would be reflected
     in the data member copy-constructed from it.
     Therefore, where needed, we rely on ``pressio::ops::clone()``,
     to ensure the correct semantics and achieve the desired immutability.
