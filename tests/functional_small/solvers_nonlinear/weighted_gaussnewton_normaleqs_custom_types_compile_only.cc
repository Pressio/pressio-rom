
#include "pressio/type_traits.hpp"
#include "pressio/ops.hpp"
#include <optional>

struct CustomVecB{};
struct CustomMat{};

struct MyProblem{
  using state_type    = Eigen::VectorXd;
  using residual_type = CustomVecB;
  using jacobian_type = CustomMat;
  state_type createState()       const { return state_type{};    }
  residual_type createResidual() const { return residual_type{}; }
  jacobian_type createJacobian() const { return jacobian_type{}; }
  void residualAndJacobian(const state_type& /*x*/,
			   residual_type& /*r*/,
			   std::optional<jacobian_type*> /*Jo*/) const{}
};

struct Weigher{
  void operator()(const CustomVecB & /*operand*/, CustomVecB & /*result*/) const{}
  void operator()(const CustomMat  & /*operand*/, CustomMat  & /*result*/) const{}
};

using my_hessian_type  = Eigen::MatrixXd;
using my_gradient_type = Eigen::VectorXd;

namespace pressio{
template<> struct Traits<CustomVecB>{
  static constexpr int rank = 1;
  using scalar_type = double;
};
template<> struct Traits<CustomMat>{
  static constexpr int rank = 2;
  using scalar_type = double;
};

namespace ops{
double norm2(const CustomVecB &){ return {}; }
double dot(const CustomVecB &, const CustomVecB &){ return {}; }
void product(transpose, nontranspose, double, const CustomMat &, const CustomMat &, double, my_hessian_type &){}
void product(transpose, double, const CustomMat &, const CustomVecB &, double, my_gradient_type &){}
}//end namespace ops
}//end namespace pressio

#include "pressio/solvers_nonlinear_gaussnewton.hpp"

struct MyLinSolver{
  void solve(const my_hessian_type & /*A*/,
	     const my_gradient_type & /*b*/,
	     typename MyProblem::state_type & /*x*/){}
};

int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);
  {
    using namespace pressio;
    using problem_t  = MyProblem;
    using state_t    = typename problem_t::state_type;
    problem_t sys;
    state_t y;
    MyLinSolver ls{}; 
    Weigher W{};
    auto nonLinSolver = create_gauss_newton_solver(sys, ls, W);
    nonLinSolver.solve(y);
    (void)y;
    std::cout << "PASSED" << std::endl;
  }
  PRESSIOLOG_FINALIZE();
}
