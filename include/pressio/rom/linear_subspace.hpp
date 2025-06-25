
#ifndef PRESSIO_ROM_LINEAR_SUBSPACE_HPP_
#define PRESSIO_ROM_LINEAR_SUBSPACE_HPP_

namespace pressio{ namespace rom{

// note: these are inside the impl namespace,
// do not use them outside
namespace impl{

template<class ReducedStateType, class = void>
struct CreateReducedState;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class ReducedStateType>
struct CreateReducedState<
  ReducedStateType,
  std::enable_if_t< ::pressio::is_vector_eigen<ReducedStateType>::value >
  >
{
  template<class BasisType>
  ReducedStateType operator()(const BasisType & basis){
    return ReducedStateType(::pressio::ops::extent(basis, 1));
  }
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<class ReducedStateType>
struct CreateReducedState<
  ReducedStateType,
  std::enable_if_t< ::pressio::is_vector_kokkos<ReducedStateType>::value >
  >
{
  template<class BasisType>
  ReducedStateType operator()(const BasisType & basis){
    return ReducedStateType("tmp", ::pressio::ops::extent(basis, 1));
  }
};
#endif

} // end namespace impl

template <class BasisMatrixType, class FullStateType, class ReducedStateType>
class TrialColumnSubspace
{
public:
  using reduced_state_type = ReducedStateType;
  using basis_matrix_type  = std::remove_cv_t<BasisMatrixType>;
  using full_state_type    = std::remove_cv_t<FullStateType>;

private:
  basis_matrix_type basis_;
  full_state_type translation_;
  bool isAffine_;
  basis_matrix_type * dummy_ = nullptr;

public:
  TrialColumnSubspace() = default;

  TrialColumnSubspace(TrialColumnSubspace const &) = delete;
  TrialColumnSubspace& operator=(TrialColumnSubspace const&) = delete;
  TrialColumnSubspace(TrialColumnSubspace &&) = default;
  TrialColumnSubspace& operator=(TrialColumnSubspace &&) = default;

  TrialColumnSubspace(const basis_matrix_type & basis,
		      const full_state_type & translation,
		      bool isAffine)
    : basis_(basis),
      translation_(translation),
      isAffine_(isAffine)
  {
    setShiftToZeroIfNonAffine();
  }

  TrialColumnSubspace(basis_matrix_type && basis,
		      full_state_type && translation,
		      bool isAffine)
    : basis_(basis),
      translation_(std::move(translation)),
      isAffine_(isAffine)
  {
    setShiftToZeroIfNonAffine();
  }

  TrialColumnSubspace(const basis_matrix_type & basis,
		      full_state_type && translation,
		      bool isAffine)
    : basis_(basis),
      translation_(std::move(translation)),
      isAffine_(isAffine)
  {
    setShiftToZeroIfNonAffine();
  }

  TrialColumnSubspace(basis_matrix_type && basis,
		      const full_state_type & translation,
		      bool isAffine)
    : basis_(basis),
      translation_(translation),
      isAffine_(isAffine)
  {
    setShiftToZeroIfNonAffine();
  }

  ~TrialColumnSubspace() = default;

  //
  // methods
  //
  reduced_state_type createReducedState() const{
    auto result = impl::CreateReducedState<ReducedStateType>()(basis_);
    using sc_t = typename ::pressio::Traits<ReducedStateType>::scalar_type;
    ::pressio::ops::fill(result, sc_t(0));
    return result;
  }

  full_state_type createFullState() const
  {
    // we need to use clone here because full_state_type might
    // NOT have value semantics so we need to ensure a new object
    // is created every time
    auto result = ::pressio::ops::clone(translation_);
    using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    ::pressio::ops::fill(result, sc_t(0));
    return result;
  }

  void mapFromReducedStateWithoutTranslation(const reduced_state_type & latState,
					     full_state_type & fullState) const
  {
    using basis_sc_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    using full_state_sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    constexpr auto alpha = static_cast<basis_sc_t>(1);
    constexpr auto beta  = static_cast<full_state_sc_t>(0);
    ::pressio::ops::product(::pressio::nontranspose(), alpha,
			    basis_, latState, beta, fullState);
  }

  void mapFromReducedState(const reduced_state_type & latState,
			   full_state_type & fullState) const
  {
    // always do y = phi*latState
    mapFromReducedStateWithoutTranslation(latState, fullState);

    if (isAffine_){
      // update full state to account for translation
      using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
      constexpr auto one = static_cast<sc_t>(1);
      ::pressio::ops::update(fullState, one, translation_, one);
    }
  }

  full_state_type createFullStateFromReducedState(const reduced_state_type & latState) const
  {
    auto fomState = this->createFullState();
    this->mapFromReducedState(latState, fomState);
    return fomState;
  }

  bool isColumnSpace() const{ return true; }
  bool isRowSpace() const{ return false; }
  const full_state_type & translationVector() const{ return translation_; }

  std::size_t dimension() const{
    // this is a column subspace, so num of columns defines the dimension
    // of the subspace
    return ::pressio::ops::extent(basis_, 1);
  }

  const basis_matrix_type & basisOfTranslatedSpace() const{
    return basis_;
  }

  const basis_matrix_type & basis() const{
    if (isAffine_){
      return basis_;
    } else{
      return *dummy_;
    }
  }

private:
  void setShiftToZeroIfNonAffine(){
    if (!isAffine_){
      ::pressio::ops::fill(translation_, 0);
    }
  }
};


/*
  ------------------------------------------------------------------------------
  create_trial_column_subspace
  ------------------------------------------------------------------------------

  Factory function to construct a `TrialColumnSubspace` object, which defines a
  reduced-order trial subspace.

  Parameters:
    - basisMatrix: A matrix whose columns span the reduced trial subspace.
                   It is typically a tall matrix mapping a reduced state to full state.
    - offset:      A full-order state used as an offset in the affine case.
    - isAffine:    A boolean flag indicating whether the subspace is
                   affine (true) or linear (false).

  Returns:
    - A `TrialColumnSubspace<basis_matrix_type, full_state_type, ReducedStateType>` instance
      representing either an affine or linear trial subspace.

  A note on the use of universal references:
  The parameters `basisMatrix` and `offset` are declared as `BasisMatrixType&&`
  and `FullStateType&&` respectively, where `BasisMatrixType` and `FullStateType`
  are template parameters. This makes them *universal references* (also known as
  forwarding references), enabling perfect forwarding.
  Perfect forwarding allows:
    - Passing lvalues without unnecessary copies.
    - Passing rvalues while preserving their move semantics.
    - Supporting const, non-const, lvalue, and rvalue variants uniformly.

  This gives maximum flexibility and performance:
    - If the caller provides an lvalue, it is passed as a reference.
    - If the caller provides an rvalue (e.g., `std::move(basis)`), it is moved into
      the `TrialColumnSubspace`, avoiding a deep copy.
    - These references are then forwarded using `std::forward<...>` to the
      constructor of `TrialColumnSubspace`, preserving the original value category.
    - Enables generic and efficient construction.
    - Avoids unnecessary copies or moves.

  Important:
    - Not all parameters must be rvalues — the function accepts **any combination**
      of lvalues and rvalues.
      For example: Both lvalues, Both rvalues, One lvalue and one rvalue (in any order)
    - Each parameter is treated independently and forwarded accordingly.
*/

template<
  class ReducedStateType,
  class BasisMatrixType,
  class FullStateType
>
auto create_trial_column_subspace(BasisMatrixType && basisMatrix,
				  FullStateType && offset,
				  bool isAffine)
{
  // figure out the "raw" type
  using basis_matrix_type = mpl::remove_cvref_t<BasisMatrixType>;
  using full_state_type = mpl::remove_cvref_t<FullStateType>;

  using return_t = TrialColumnSubspace<basis_matrix_type, full_state_type, ReducedStateType>;
  return return_t(std::forward<BasisMatrixType>(basisMatrix),
		  std::forward<FullStateType>(offset),
		  isAffine);
}

}}
#endif  // PRESSIO_ROM_LINEAR_SUBSPACE_HPP_
