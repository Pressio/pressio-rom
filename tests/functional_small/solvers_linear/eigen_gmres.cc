
#include <gtest/gtest.h>
#include "pressio/solvers_linear.hpp"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

class OperatorWrapper;

namespace Eigen {
  namespace internal {
    template<>
    struct traits<OperatorWrapper>
      : public Eigen::internal::traits<Eigen::Matrix<double,-1,-1> >
    {};
  }
}

class OperatorWrapper : public Eigen::EigenBase<OperatorWrapper> {
public:
  using Scalar = double;
  using RealScalar = double;
  using StorageIndex = int;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
  };

  OperatorWrapper() : mp_mat(0) {}

  Index rows() const { return mp_mat->rows(); }
  Index cols() const { return mp_mat->cols(); }

  template<typename Rhs>
  Eigen::Product<OperatorWrapper, Rhs, Eigen::AliasFreeProduct>
  operator*(const Eigen::MatrixBase<Rhs>& x) const
  {
    using r_t = Eigen::Product<OperatorWrapper,Rhs,Eigen::AliasFreeProduct>;
    return r_t(*this, x.derived());
  }

  void attachMyMatrix(const Eigen::MatrixXd &mat) {
    mp_mat = &mat;
  }

  const Eigen::MatrixXd my_matrix() const { return *mp_mat; }


  template<class OperandT, class ResultT>
  void applyAndAddTo(OperandT const & operand, ResultT & out) const {
    //
    // compute: out += operator * operand

    out += *mp_mat * operand;
  }

private:
  const Eigen::MatrixXd *mp_mat;
};


namespace Eigen {
  namespace internal {

    template<typename Rhs>
    struct generic_product_impl<
      OperatorWrapper, Rhs, DenseShape, DenseShape, GemvProduct
      >
      : generic_product_impl_base<
      OperatorWrapper, Rhs, generic_product_impl<OperatorWrapper,Rhs>
      >
    {
      typedef typename Product<OperatorWrapper,Rhs>::Scalar Scalar;

      template<typename Dest>
      static void scaleAndAddTo(
	 Dest& dst,
	 const OperatorWrapper& lhs,
	 const Rhs& rhs,
	 const Scalar& alpha)
      {
	// This method should implement "dst += alpha * lhs * rhs" inplace,
	// however, for iterative solvers, alpha is always equal to 1,
	// so let's not bother about it.
	assert(alpha==Scalar(1) && "scaling is not implemented");
	EIGEN_ONLY_USED_FOR_DEBUG(alpha);

	// // // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
	// // // but let's do something fancier (and less efficient):
	// for(Index i=0; i<lhs.cols(); ++i)
	//    dst += rhs(i) * lhs.my_matrix().col(i);

	lhs.applyAndAddTo(rhs, dst);
      }
    };
  }
}

TEST(solvers_linear_eigen, gmres)
{
  srand(3451677);

  int n = 10;
  Eigen::MatrixXd S = Eigen::MatrixXd::Random(n,n);
  S = S.transpose()*S;
  std::cout << S << "\n";

  OperatorWrapper A;
  A.attachMyMatrix(S);

  Eigen::VectorXd b(n), x(n);
  b.setRandom();
  x.setZero();


  Eigen::GMRES<OperatorWrapper, Eigen::IdentityPreconditioner> gmres;
  gmres.compute(A);
  gmres.setMaxIterations(8);
  x = gmres.solve(b);
  std::cout << "GMRES: #iterations: " << gmres.iterations()
	    << ", estimated error: " << gmres.error()
	    << std::endl;
  std::cout << "my solution \n " << x << std::endl;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

struct MyOperator{
  using scalar_type = double;

  MyOperator(int n){
    Eigen::MatrixXd S = Eigen::MatrixXd::Random(n,n);
    m_A = S.transpose()*S;
  }

  int rows() const { return m_A.rows(); }
  int cols() const { return m_A.cols(); }

  template<class OperandT, class ResultT>
  void applyAndAddTo(OperandT const & operand, ResultT & out) const {
    out += m_A * operand;
  }

private:
  Eigen::MatrixXd m_A;
};

TEST(solvers_linear_eigen, gmres2)
{
  srand(3451677);

  int n = 10;
  MyOperator A(n);
  Eigen::VectorXd b(n), x(n);
  b.setRandom();
  x.setZero();

  using operator_t = MyOperator;
  using tag = pressio::linsol::iterative::GMRES;
  using solver_t = pressio::linsol::Solver<tag, operator_t>;
  solver_t solver;
  solver.setMaxIterations(8);
  solver.solve(A, b, x);

  std::cout << "my error \n " << solver.finalError() << std::endl;
  std::cout << "my solution \n " << x << std::endl;
}
