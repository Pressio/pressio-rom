
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

struct CustomVecA{};
struct CustomVecB{};
struct CustomMat{};

namespace pressio{

template<> struct Traits<CustomVecA>{
  static constexpr int rank = 1;
  using scalar_type = double;
};
template<> struct Traits<CustomVecB>{
  static constexpr int rank = 1;
  using scalar_type = double;
};
template<> struct Traits<CustomMat>{
  static constexpr int rank = 2;
  using scalar_type = double;
};

namespace ops{
double norm2(const CustomVecA &){ return {}; }
double norm2(const CustomVecB &){ return {}; }
void deep_copy(CustomVecA &, const CustomVecA &){}
void scale(CustomVecA &, double){}
void update(CustomVecA &, double,
	    const CustomVecA &, double,
	    const CustomVecA &, double){}
void update(CustomVecA &, double,
	    const CustomVecA &, double){}
}//end namespace ops
}//end namespace pressio

#include "pressio/solvers_nonlinear_newton.hpp"

struct MyProblem{
  using state_type = CustomVecA;
  using residual_type = CustomVecB;
  using jacobian_type = CustomMat;
  state_type createState() const { return CustomVecA{}; }
  residual_type createResidual() const { return CustomVecB{}; }
  jacobian_type createJacobian() const { return CustomMat{}; }
  void residualAndJacobian(const state_type& /*x*/,
			   residual_type& /*r*/,
			   std::optional<jacobian_type*> /*Jo*/) const{}
};

struct MyLinSolver{
  void solve(const CustomMat & /*A*/,
	     const CustomVecB & /*b*/,
	     CustomVecA & /*x*/){}
};

TEST(solvers_nonlinear, newton_compile_only)
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);
  {

    using namespace pressio;
    namespace pnonls = pressio::nonlinearsolvers;

    using problem_t  = MyProblem;
    using state_t    = typename problem_t::state_type;

    problem_t sys;
    state_t y;
    MyLinSolver ls{};
    auto nonLinSolver = pnonls::create_newton_solver(sys, ls);
    nonLinSolver.solve(y);
    (void)y;
  }

  PRESSIOLOG_FINALIZE();
}
