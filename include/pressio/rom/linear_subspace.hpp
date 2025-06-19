
#ifndef PRESSIO_ROM_LINEAR_SUBSPACE_HPP_
#define PRESSIO_ROM_LINEAR_SUBSPACE_HPP_

namespace pressio{ namespace rom{

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

}

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
  std::size_t dimension() const{ return ::pressio::ops::extent(basis_, 1); }

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


template<
  class ReducedStateType,
  class BasisMatrixType,
  class FullStateType
>
auto create_trial_column_subspace(BasisMatrixType && basisMatrix,
				  FullStateType && offset,
				  bool isAffine)
{
  using basis_matrix_type = mpl::remove_cvref_t<BasisMatrixType>;
  using full_state_type = mpl::remove_cvref_t<FullStateType>;

  using ret_t = TrialColumnSubspace<basis_matrix_type, full_state_type, ReducedStateType>;
  return ret_t(std::forward<BasisMatrixType>(basisMatrix),
	       std::forward<FullStateType>(offset),
	       isAffine);
}

}}
#endif  // PRESSIO_ROM_LINEAR_SUBSPACE_HPP_
