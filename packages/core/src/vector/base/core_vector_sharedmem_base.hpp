
#ifndef CORE_VECTOR_BASE_VECTOR_SHAREDMEM_BASE_HPP_
#define CORE_VECTOR_BASE_VECTOR_SHAREDMEM_BASE_HPP_

#include "../core_vector_traits.hpp"

namespace rompp{ namespace core{

template<typename derived_type>
class VectorSharedMemBase
  : private core::details::CrtpBase<
     VectorSharedMemBase<derived_type>>
{

  static_assert(details::traits<derived_type>::is_shared_mem==1,
  "OOPS: distributed concrete vector inheriting from sharedMem base!");

  using this_t = VectorSharedMemBase<derived_type>;
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

public:
   ord_t size() const{
    return this->underlying().sizeImpl();
  };

  template <typename T,
  	    core::meta::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  void putScalar(T value) {
    this->underlying().putScalarImpl(value);}

private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;
  VectorSharedMemBase() = default;
  ~VectorSharedMemBase() = default;

};//end class

}}//end namespace rompp::core
#endif
