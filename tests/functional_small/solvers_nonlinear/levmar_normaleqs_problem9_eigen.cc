
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_levmarq.hpp"
#include "./problems/problem9.hpp"

#include <iomanip>

template <typename problem_t, typename solver_t>
void testC1(std::string & sentinel,
            problem_t & problem,
            solver_t & solver)
{
  using namespace pressio;
  auto x = problem.createState();
  x << 1.3, 6.5e-1, 6.5e-1, 7.0e-1,
    6.0e-1, 3.0, 5.0, 7.0, 2.0, 4.5, 5.5;

  solver.setUpdateCriterion(nonlinearsolvers::Update::LMSchedule1);
  solver.setStopTolerance(1e-6);
  solver.solve(problem, x);
  std::cout << std::setprecision(14) << x << " ";
  std::cout << std::endl;

  // this gold solution corresponds to objective = 2.006887e-02
  // which matches the one in http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
  // 2.006887e-02 = 0.5 * 0.2003440**2
  Eigen::VectorXd gold(11);
  gold << 1.3099771465073,  0.43155377520055, 0.63366169194236,
    0.59943051899443, 0.75418317304122, 0.90428874189697,
    1.3658117486587,  4.8236989280488,  2.3986848749122,
    4.5688745952807,  5.6753414737204;

  if (!gold.isApprox(x)){ sentinel = "FAILED"; }
}

template <typename problem_t, typename solver_t>
void testC2(std::string & sentinel,
            problem_t & problem,
            solver_t & solver)
{
  using namespace pressio;
  auto x = problem.createState();
  x << 1.3, 6.5e-1, 6.5e-1, 7.0e-1,
    6.0e-1, 3.0, 5.0, 7.0, 2.0, 4.5, 5.5;

  solver.setUpdateCriterion(nonlinearsolvers::Update::LMSchedule2);
  solver.setStopTolerance(1e-6);
  solver.solve(problem, x);
  std::cout << std::setprecision(14) << x << " ";
  std::cout << std::endl;

  // this gold solution corresponds to objective = 2.006887e-02
  // which matches the one in http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
  // 2.006887e-02 = 0.5 * 0.2003440**2
  Eigen::VectorXd gold(11);
  gold <<  1.3099771465073, 0.43155377520055, 0.63366169194236,
    0.59943051899443, 0.75418317304122, 0.90428874189697,
    1.3658117486587,  4.8236989280488, 2.3986848749122,
    4.5688745952807, 5.6753414737204;

  if (!gold.isApprox(x)){ sentinel = "FAILED"; }
}

int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::info);

  using namespace pressio;
  namespace pnonls = pressio::nonlinearsolvers;

  using problem_t = solvers::test::Problem9<double>;
  using hessian_t = Eigen::MatrixXd;

  problem_t problem;

  using lin_tag      = linearsolvers::direct::HouseholderQR;
  using lin_solver_t = linearsolvers::Solver<lin_tag, hessian_t>;
  lin_solver_t linSolver;
  std::string sentinel = "PASSED";
  auto solver = pnonls::create_levenberg_marquardt_solver(problem, linSolver);
  testC1(sentinel, problem, solver);
  std::cout << "\n" << std::endl;

  testC2(sentinel, problem, solver);
  std::cout << "\n" << std::endl;

  std::cout << sentinel << "\n";
  PRESSIOLOG_FINALIZE();
  return 0;
}
