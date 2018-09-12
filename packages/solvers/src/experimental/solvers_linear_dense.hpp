
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_DENSE_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_DENSE_HPP

#include "solvers_linear_base.hpp"


namespace solvers {

// Forward declarations
struct LinearSolvers;


/**
 * Linear dense solver
 *
 * @tparam SolverT linear dense solver
 * @tparam MatrixT matrix defining the linear system
 * @tparam PolicyT policy that implements the solution algorithm
 */
template<
  typename SolverT,
  typename MatrixT,
  typename PolicyT
>
class LinearDenseSolver
  : public LinearSolverBase<
      SolverT,
      MatrixT,
      PolicyT,
      LinearDenseSolver<
        SolverT,
        MatrixT,
        PolicyT
      >
    >
{

  private:

    friend LinearSolvers;
    typedef LinearSolverBase<SolverT, MatrixT, PolicyT, LinearDenseSolver<SolverT, MatrixT, PolicyT>> base_type;


  public:

    LinearDenseSolver(LinearDenseSolver&& other) : base_type(std::move(other)) {}


    template <typename T>
    auto _solve(const T& b) {
      auto solver = this->getSolver();
      return PolicyT::solve(solver, b);
    }


  protected:

    LinearDenseSolver() : base_type() {};

    LinearDenseSolver(std::shared_ptr<SolverT> solver) : base_type(solver) {};

};

} //end namespace solvers

#endif
