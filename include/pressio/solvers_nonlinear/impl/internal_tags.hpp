
#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_INTERNAL_TAGS_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_INTERNAL_TAGS_HPP_

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

// this tag is inside the impl namespace because we do NOT want
// to expose it outside, this is an impl detail
struct SystemTag{};

// this represent ths Q^T *r for QR gauss newton
struct QTransposeResidualTag{};


struct NewtonTag{};
struct MatrixFreeNewtonTag{};
struct GaussNewtonNormalEqTag{};
struct WeightedGaussNewtonNormalEqTag{};
struct CompactWeightedGaussNewtonNormalEqTag{};
struct LevenbergMarquardtNormalEqTag{};
struct GaussNewtonQrTag{};

}}}
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_INTERNAL_TAGS_HPP_
