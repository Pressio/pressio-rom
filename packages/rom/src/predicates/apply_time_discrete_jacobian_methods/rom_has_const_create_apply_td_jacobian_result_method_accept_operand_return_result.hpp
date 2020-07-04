
#ifndef ODE_HAS_CONST_CREATE_TD_JAC_OBJECT_METHOD_RETURN_RESULT_HPP_
#define ODE_HAS_CONST_CREATE_TD_JAC_OBJECT_METHOD_RETURN_RESULT_HPP_

namespace pressio{ namespace rom{ namespace meta {

template <
  typename T, typename operand_t, typename result_t,
  typename = void
  >
struct has_const_create_td_jacobian_method_return_result
  : std::false_type{};


template <typename T, typename operand_t, typename result_t>
struct has_const_create_td_jacobian_method_return_result<
  T, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    !std::is_void<result_t>::value and
    mpl::is_same<
      result_t,
      decltype(
	       std::declval<T const>().createApplyTimeDiscreteJacobianResult
          (
            std::declval<operand_t const & >()
          )
	       )
      >::value
    >
  > : std::true_type{};


}}} 
#endif
