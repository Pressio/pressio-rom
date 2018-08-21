
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_BASE_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_BASE_HPP

#include <memory>
#include <type_traits>

#include "meta/core_meta_static_checks.hpp"
#include "matrix/core_matrix_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"


namespace solvers {


/**
 * @brief Base class for linear solver implemented through CRTP
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a linear solver. 
 * Objects of the class cannot be created directly. To create a solver,
 * use the factory class LinearSolvers.
 */
template<
  typename SolverT,
  typename MatrixT,
  typename PolicyT,
  typename Derived
>
class LinearSolverBase {

  public: 


    /**
     * @brief  Initialize a new linear solver
     *
     * @param  A Matrix representing the new linear system to solve
     */
    template <
      typename CompatibleMatrixT,
      typename std::enable_if<
        core::meta::are_matrix_compatible<
          MatrixT,
          CompatibleMatrixT
        >::value, MatrixT*
      >::type = nullptr
    >
    void resetLinearSystem(const CompatibleMatrixT& A) {
      PolicyT::resetLinearSystem(solver_, A);
    }


    /**
     * @brief  Solve the linear system
     *
     * @param  b is the RHS vector
     * @return Solution vector
     */
    template <
      typename VectorLT,
      typename std::enable_if<
        core::meta::are_vector_matrix_compatible<
          VectorLT,
          MatrixT
        >::value, MatrixT*
      >::type = nullptr
    >
    auto solve(const VectorLT& b) {
      return this->underlying()._solve(b);  
    }


    /**
     * @brief  Specify and solve the linear system
     *
     * @param  A is the system matrix
     * @param  b is the RHS vector
     * @return Solution vector
     */
    template <
      typename CompatibleMatrixT, 
      typename VectorLT,
      typename std::enable_if<
        core::details::matrix_traits<CompatibleMatrixT>::matrix_class != core::details::WrappedClass::Undefined,
        CompatibleMatrixT*
      >::type = nullptr
    >
    auto solve(const CompatibleMatrixT& A, const VectorLT& b) {
      this->resetLinearSystem(A);
      return this->solve(b);
    }


    /**
     * @brief  Solve the linear system
     *
     * @param  b is the RHS vector
     * @param  x is the solution vector
     * @return void
     */
    template <
      typename VectorLT, 
      typename VectorRT,
      typename std::enable_if<
        core::meta::are_vector_compatible<
          VectorLT, 
          VectorRT
        >::value &&
        core::details::vector_traits<VectorLT>::vector_class != core::details::WrappedClass::Undefined,
        VectorLT*
      >::type = nullptr
    >
    void solve(const VectorLT& b, VectorRT& x) {
      x = VectorRT(*this->solve(b).data());
    }


    /**
     * @brief  Specify and solve the linear system
     *
     * @param  A is the system matrix
     * @param  b is the RHS vector
     * @param  x is the solution vector
     * @return void
     */
    template <
      typename CompatibleMatrixT,
      typename VectorLT,
      typename VectorRT
    >
    void solve(const CompatibleMatrixT& A, const VectorLT& b, VectorRT& x) {
      this->resetLinearSystem(A);
      this->solve(b, x);
    }


  protected:

    LinearSolverBase() : solver_(nullptr) {};


    LinearSolverBase(std::shared_ptr<SolverT> solver) : solver_(solver) {}


    LinearSolverBase(LinearSolverBase&& other) : solver_(std::move(other.solver_)) {}


    LinearSolverBase(const LinearSolverBase&) = delete;


    ~LinearSolverBase() = default;


    std::shared_ptr<SolverT> getSolver() {
      return solver_;
    }


  private:

    Derived& underlying() {
      return static_cast<Derived&>(*this);
    }
  

    Derived const& underlying() const {
      return static_cast<Derived const&>(*this);
    }  


  private:

    std::shared_ptr<SolverT> solver_;
};

} //end namespace solvers

#endif
