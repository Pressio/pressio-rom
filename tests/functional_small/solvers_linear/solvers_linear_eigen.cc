
#include <gtest/gtest.h>
#include "pressio/solvers_linear.hpp"

TEST(solvers_linear_eigen, sparse_iterative_lscg)
{
  using matrix_t = Eigen::SparseMatrix<double>;
  using vector_t = Eigen::VectorXd;

  // Define linear system
  vector_t b(2);
  b << 7.0, 1.0;

  matrix_t A(2, 2);
  A.insert(0, 0) +=  3.0;
  A.insert(0, 1) += -1.0;
  A.insert(1, 0) +=  2.0;
  A.insert(1, 1) +=  3.0;

  vector_t y(2);
  using tag = pressio::linsol::iterative::LSCG;
  using solver_t = pressio::linsol::Solver<tag, matrix_t>;
  solver_t solver;
  solver.solve(A,b,y);

  const auto e1 = std::abs(y(0) - (2.0));
  const auto e2 = std::abs(y(1) - (-1.0));
  ASSERT_TRUE(e1 <= 1e-13 and e2 <= 1e-13);
}

#define PRESSIO_SOLVERS_LINEAR_EIGEN_DENSE_UTEST(TAGIN) \
  using matrix_t = Eigen::MatrixXd; \
  using vector_t = Eigen::VectorXd; \
  vector_t b(3); \
  b << 4.0, 1.0, 3.0; \
  matrix_t A(3, 2);\
  A(0, 0) = 1.0;\
  A(0, 1) = 3.0;\
  A(1, 0) = 2.0;\
  A(1, 1) = 4.0;\
  A(2, 0) = 1.0;\
  A(2, 1) = 6.0;\
  vector_t y(2);\
  using solver_t = pressio::linsol::Solver<TAGIN, matrix_t>;\
  solver_t solver;\
  solver.solve(A,b,y); \
  const auto e1 = std::abs(y(0) - (-0.377)); \
  const auto e2 = std::abs(y(1) - (0.662)); \
  ASSERT_TRUE(e1 <= 1e-3 and e2<= 1e-3);\


TEST(solvers_linear_eigen, dense_iterative_lscg)
{
  using tag = pressio::linsol::iterative::LSCG;
  PRESSIO_SOLVERS_LINEAR_EIGEN_DENSE_UTEST(tag);
}

TEST(solvers_linear_eigen, dense_direct_colpivqr)
{
  using tag = pressio::linsol::direct::ColPivHouseholderQR;
  PRESSIO_SOLVERS_LINEAR_EIGEN_DENSE_UTEST(tag);
}

TEST(solvers_linear_eigen, dense_direct_houseqr)
{
  using tag = pressio::linsol::direct::HouseholderQR;
  PRESSIO_SOLVERS_LINEAR_EIGEN_DENSE_UTEST(tag);
}


// Define all tags we want to test
using SolverTags = ::testing::Types<
  pressio::linsol::iterative::CG,
  pressio::linsol::iterative::LSCG,
  pressio::linsol::iterative::Bicgstab,
  pressio::linsol::direct::HouseholderQR,
  pressio::linsol::direct::ColPivHouseholderQR,
  pressio::linsol::direct::PartialPivLU
>;

// Define typed test fixture
template <typename T>
class LinearSolverDefaultConstructibilityTest : public ::testing::Test {};

// Register typed test suite
TYPED_TEST_SUITE(LinearSolverDefaultConstructibilityTest, SolverTags);

TYPED_TEST(LinearSolverDefaultConstructibilityTest, IsDefaultConstructible)
{
  using TagType = TypeParam;
  using Solver = pressio::linsol::Solver<TagType, Eigen::MatrixXd>;
  static_assert(std::is_default_constructible<Solver>::value,
		"Solver is not default constructible");
  Solver solver;  // actually default construct it
  (void)solver;   // avoid unused warning
}

TYPED_TEST(LinearSolverDefaultConstructibilityTest, IsDefaultConstructible2)
{
  using TagType = TypeParam;
  using Solver = pressio::linsol::Solver<TagType, Eigen::SparseMatrix<double>>;
  static_assert(std::is_default_constructible<Solver>::value,
		"Solver is not default constructible");
  Solver solver;  // actually default construct it
  (void)solver;   // avoid unused warning
}
