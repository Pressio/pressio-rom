
#include <gtest/gtest.h>
#include "pressio/solvers_linear.hpp"

TEST(solvers_linear_kokkos, dense_getrs)
{
  using d_layout = Kokkos::LayoutLeft;
  using exe_space = Kokkos::DefaultExecutionSpace;
  using k1d_d = Kokkos::View<double*, d_layout, exe_space>;
  using k1d_h = typename k1d_d::host_mirror_type;
  using k2d_d = Kokkos::View<double**, d_layout, exe_space>;
  using k2d_h = typename k2d_d::host_mirror_type;

   // size of problems, rows and cols
  constexpr int Nr = 4;
  constexpr int Nc = 4;

  // create/fill the matrix
  k2d_h A_h("Ah", Nr, Nc);
  A_h(0,0)=1.; A_h(0,1)=3.; A_h(0,2)=2.; A_h(0,3)=3.;
  A_h(1,0)=2.; A_h(1,1)=4.; A_h(1,2)=1.; A_h(1,3)=2.;
  A_h(2,0)=1.; A_h(2,1)=6.; A_h(2,2)=3.; A_h(2,3)=1.;
  A_h(3,0)=0.; A_h(3,1)=2.; A_h(3,2)=2.; A_h(3,3)=1.;
  k2d_d A_d("Ad", Nr, Nc);
  Kokkos::deep_copy(A_d, A_h);

  // Construct rhs
  k1d_h b_h("bh", Nr);
  b_h(0) = 4.;
  b_h(1) = 1.;
  b_h(2) = 3.0;
  b_h(3) = 2.;
  k1d_d b_d("bd", Nr);
  Kokkos::deep_copy(b_d, b_h);

  using solver_tag   = pressio::linsol::direct::PartialPivLU;
  using linear_solver_t = pressio::linsol::Solver<solver_tag, k2d_d>;
  linear_solver_t lsObj;
  k1d_d x_d("xd", Nr);
  lsObj.solve(A_d, b_d, x_d);

  // create host of solution
  k1d_h x_h("xh", Nc);
  Kokkos::deep_copy(x_h, x_d);

  std::cout << x_h(0) << " "
      << x_h(1) << " "
      << x_h(2) << " "
      << x_h(3) << std::endl;

  std::vector<double> sum(Nc, 0.0);
  for (int i=0; i<Nr; ++i)
    for (int j=0; j<Nc; ++j)
      sum[i] += A_h(i,j) * x_h(j);

  std::string strOutput  = "PASSED";
  for (int i=0; i<Nr; ++i){
    const auto err = sum[i] - b_h(i);
    EXPECT_TRUE( err <= 1e-12);
  }

}
