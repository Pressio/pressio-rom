
#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_LEVMAR_DAMPING_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_LEVMAR_DAMPING_HPP_

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

template<class T>
class LevenbergMarquardtDamping
{
  static_assert(std::is_floating_point<T>::value, "");
  using value_type = T;
  value_type v_{1};
public:
  LevenbergMarquardtDamping & operator = (T v){ v_ = v; return *this; }
  LevenbergMarquardtDamping & operator *= (T v){ v_ *= v; return *this; }
  operator value_type () const { return v_; }
};

}}}
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_LEVMAR_DAMPING_HPP_
