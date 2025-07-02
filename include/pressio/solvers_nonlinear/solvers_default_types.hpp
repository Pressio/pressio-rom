
#ifndef PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_DEFAULT_TYPES_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_DEFAULT_TYPES_HPP_

namespace pressio{
namespace nonlinearsolvers{

template<class T, class = void>
struct normal_eqs_default_types{
  using hessian_type  = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct normal_eqs_default_types<
  T, std::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
  using gradient_type = T;

  static hessian_type createHessian(const T & v){
    const auto ext = ::pressio::ops::extent(v, 0);
    return hessian_type(ext, ext);
  }
};
#endif

template<class T> using normal_eqs_default_hessian_t =
  typename normal_eqs_default_types<T>::hessian_type;
template<class T> using normal_eqs_default_gradient_t =
  typename normal_eqs_default_types<T>::gradient_type;

// ====================================================================

template<class T, class = void>
struct valid_state_for_least_squares : std::false_type{};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct valid_state_for_least_squares<
  T, std::enable_if_t< ::pressio::is_vector_eigen<T>::value >
  > : std::true_type{};
#endif

}}
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_DEFAULT_TYPES_HPP_
