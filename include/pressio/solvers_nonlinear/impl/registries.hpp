/*
//@HEADER
// ************************************************************************
//
// registries.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_REGISTRIES_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_REGISTRIES_HPP_

#include "levmar_damping.hpp"
#include "instance_or_reference_wrapper.hpp"

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

#define GETMETHOD(N) \
  template<class Tag, std::enable_if_t< std::is_same<Tag, Tag##N >::value, int> = 0> \
  auto & get(){ return d##N##_; } \
  template<class Tag, std::enable_if_t< std::is_same<Tag, Tag##N >::value, int> = 0> \
  const auto & get() const { return d##N##_; }


template<class SystemType, class InnSolverType>
class RegistryNewton
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::JacobianTag;
  using Tag5 = nonlinearsolvers::InnerSolverTag;
  using Tag6 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  r_t d3_;
  j_t d4_;
  InnSolverType * d5_;
  SystemType const * d6_;

public:
  RegistryNewton(const SystemType & system, InnSolverType & innS)
    : d1_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(&innS),
      d6_(&system){}

  RegistryNewton() = delete;
  // non-copyable
  RegistryNewton(RegistryNewton const &) = delete;
  RegistryNewton& operator=(RegistryNewton const &) = delete;
  // movable
  RegistryNewton(RegistryNewton &&) = default;
  RegistryNewton& operator=(RegistryNewton &&) = default;

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	   Tag1, Tag3, Tag4, Tag5, Tag6>::value) < 5;
  }

  GETMETHOD(1)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
};

template<class SystemType, class LinearSolverTag>
class RegistryMatrixFreeNewtonKrylov
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  r_t d3_;
  SystemType const * d4_;

public:
  using linear_solver_tag = LinearSolverTag;

  RegistryMatrixFreeNewtonKrylov(const SystemType & system)
    : d1_(system.createState()),
      d3_(system.createResidual()),
      d4_(&system){}

  RegistryMatrixFreeNewtonKrylov() = delete;
  // non-copyable
  RegistryMatrixFreeNewtonKrylov(RegistryMatrixFreeNewtonKrylov const &) = delete;
  RegistryMatrixFreeNewtonKrylov& operator=(RegistryMatrixFreeNewtonKrylov const &) = delete;
  // movable
  RegistryMatrixFreeNewtonKrylov(RegistryMatrixFreeNewtonKrylov &&) = default;
  RegistryMatrixFreeNewtonKrylov& operator=(RegistryMatrixFreeNewtonKrylov &&) = default;

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag3, Tag4>::value) < 3;
  }

  GETMETHOD(1)
  GETMETHOD(3)
  GETMETHOD(4)
};


template<class SystemType, class InnSolverType>
class RegistryGaussNewtonNormalEqs
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::JacobianTag;
  using Tag5 = nonlinearsolvers::GradientTag;
  using Tag6 = nonlinearsolvers::HessianTag;
  using Tag7 = nonlinearsolvers::InnerSolverTag;
  using Tag8 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  r_t d3_;
  j_t d4_;
  gradient_t d5_;
  hessian_t d6_;
  InnSolverType * d7_;
  SystemType const * d8_;

public:
  RegistryGaussNewtonNormalEqs(const SystemType & system, InnSolverType & innS)
    : d1_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createState()),
      d6_( hg_default::createHessian(system.createState()) ),
      d7_(&innS),
      d8_(&system){}

  RegistryGaussNewtonNormalEqs() = delete;
  // non-copyable
  RegistryGaussNewtonNormalEqs(RegistryGaussNewtonNormalEqs const &) = delete;
  RegistryGaussNewtonNormalEqs& operator=(RegistryGaussNewtonNormalEqs const &) = delete;
  // movable
  RegistryGaussNewtonNormalEqs(RegistryGaussNewtonNormalEqs &&) = default;
  RegistryGaussNewtonNormalEqs& operator=(RegistryGaussNewtonNormalEqs &&) = default;

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	   Tag1, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8>::value) < 7;
  }

  GETMETHOD(1)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
};

template<class SystemType, class InnSolverType, class WeightingOpType>
class RegistryWeightedGaussNewtonNormalEqs
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;

  using Tag1  = nonlinearsolvers::CorrectionTag;
  using Tag3  = nonlinearsolvers::ResidualTag;
  using Tag4  = nonlinearsolvers::JacobianTag;
  using Tag5  = nonlinearsolvers::WeightedResidualTag;
  using Tag6  = nonlinearsolvers::WeightedJacobianTag;
  using Tag7  = nonlinearsolvers::GradientTag;
  using Tag8  = nonlinearsolvers::HessianTag;
  using Tag9  = nonlinearsolvers::InnerSolverTag;
  using Tag10 = nonlinearsolvers::WeightingOperatorTag;
  using Tag11 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  r_t d3_;
  j_t d4_;
  r_t d5_;
  j_t d6_;
  gradient_t d7_;
  hessian_t d8_;
  InnSolverType * d9_;
  WeightingOpType const * d10_;
  SystemType const * d11_;

public:
  RegistryWeightedGaussNewtonNormalEqs(const SystemType & system,
				       InnSolverType & innS,
				       const WeightingOpType & weigher)
    : d1_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createResidual()),
      d6_(system.createJacobian()),
      d7_(system.createState()),
      d8_( hg_default::createHessian(system.createState()) ),
      d9_(&innS),
      d10_(&weigher),
      d11_(&system){}

  RegistryWeightedGaussNewtonNormalEqs() = delete;
  // non-copyable
  RegistryWeightedGaussNewtonNormalEqs(RegistryWeightedGaussNewtonNormalEqs const &) = delete;
  RegistryWeightedGaussNewtonNormalEqs& operator=(RegistryWeightedGaussNewtonNormalEqs const &) = delete;
  // movable
  RegistryWeightedGaussNewtonNormalEqs(RegistryWeightedGaussNewtonNormalEqs &&) = default;
  RegistryWeightedGaussNewtonNormalEqs& operator=(RegistryWeightedGaussNewtonNormalEqs &&) = default;

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10, Tag11>::value) < 10;
  }

  GETMETHOD(1)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
  GETMETHOD(9)
  GETMETHOD(10)
  GETMETHOD(11)
};

template<class SystemType, class InnSolverType, class WeightingOpType>
class RegistryCompactWeightedGaussNewtonNormalEqs
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;

  using Tag1  = nonlinearsolvers::CorrectionTag;
  using Tag3  = nonlinearsolvers::ResidualTag;
  using Tag4  = nonlinearsolvers::JacobianTag;
  using Tag5  = nonlinearsolvers::WeightedResidualTag;
  using Tag6  = nonlinearsolvers::WeightedJacobianTag;
  using Tag7  = nonlinearsolvers::GradientTag;
  using Tag8  = nonlinearsolvers::HessianTag;
  using Tag9  = nonlinearsolvers::InnerSolverTag;
  using Tag10 = nonlinearsolvers::WeightingOperatorTag;
  using Tag11 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  r_t d3_;
  j_t d4_;
  r_t d5_;
  j_t d6_;
  gradient_t d7_;
  hessian_t d8_;
  InnSolverType * d9_;
  WeightingOpType const * d10_;
  SystemType const * d11_;

public:
  RegistryCompactWeightedGaussNewtonNormalEqs(const SystemType & system,
					      InnSolverType & innS,
					      const WeightingOpType & weigher)
    : d1_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createResidual()),
      d6_(system.createJacobian()),
      d7_(system.createState()),
      d8_( hg_default::createHessian(system.createState()) ),
      d9_(&innS),
      d10_(&weigher),
      d11_(&system){
        // resize Wr and WJ leading dimension according to weighing operator
        pressio::ops::resize(d5_, weigher.leadingDim());
        pressio::ops::resize(d6_, weigher.leadingDim(), pressio::ops::extent(d6_, 1));
      }

  RegistryCompactWeightedGaussNewtonNormalEqs() = delete;

  // non-copyable
  RegistryCompactWeightedGaussNewtonNormalEqs(RegistryCompactWeightedGaussNewtonNormalEqs const &) = delete;
  RegistryCompactWeightedGaussNewtonNormalEqs& operator=(RegistryCompactWeightedGaussNewtonNormalEqs const &) = delete;
  // movable
  RegistryCompactWeightedGaussNewtonNormalEqs(RegistryCompactWeightedGaussNewtonNormalEqs &&) = default;
  RegistryCompactWeightedGaussNewtonNormalEqs& operator=(RegistryCompactWeightedGaussNewtonNormalEqs &&) = default;

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10, Tag11>::value) < 10;
  }

  GETMETHOD(1)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
  GETMETHOD(9)
  GETMETHOD(10)
  GETMETHOD(11)
};

template<class SystemType, class QRSolverType>
class RegistryGaussNewtonQr
{
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using QTr_t      = state_t; // type of Q^T*r
  using gradient_t = state_t; // type of J^T r

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::JacobianTag;
  using Tag5 = nonlinearsolvers::GradientTag;
  using Tag6 = nonlinearsolvers::impl::QTransposeResidualTag;
  using Tag7 = nonlinearsolvers::InnerSolverTag;
  using Tag8 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  r_t d3_;
  j_t d4_;
  gradient_t d5_;
  QTr_t d6_;
  QRSolverType * d7_;
  SystemType const * d8_;

public:
  RegistryGaussNewtonQr(const SystemType & system, QRSolverType & qrs)
    : d1_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createState()),
      d6_(system.createState()),
      d7_(&qrs),
      d8_(&system){}

  RegistryGaussNewtonQr() = delete;
  // non-copyable
  RegistryGaussNewtonQr(RegistryGaussNewtonQr const &) = delete;
  RegistryGaussNewtonQr& operator=(RegistryGaussNewtonQr const &) = delete;
  // movable
  RegistryGaussNewtonQr(RegistryGaussNewtonQr &&) = default;
  RegistryGaussNewtonQr& operator=(RegistryGaussNewtonQr &&) = default;

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	   Tag1, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8>::value) < 7;
  }

  GETMETHOD(1)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
};

template<class SystemType, class InnSolverType>
class RegistryLevMarNormalEqs
{
  using scalar_t   = system_scalar_t<SystemType>;
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;
  using lm_damp_t  = LevenbergMarquardtDamping<scalar_t>;

  using Tag1 = nonlinearsolvers::CorrectionTag;
  using Tag3 = nonlinearsolvers::ResidualTag;
  using Tag4 = nonlinearsolvers::JacobianTag;
  using Tag5 = nonlinearsolvers::GradientTag;
  using Tag6 = nonlinearsolvers::HessianTag;
  using Tag7 = nonlinearsolvers::LevenbergMarquardtUndampedHessianTag;
  using Tag8 = nonlinearsolvers::LevenbergMarquardtDampingTag;
  using Tag9 = nonlinearsolvers::InnerSolverTag;
  using Tag10 = nonlinearsolvers::impl::SystemTag;

  state_t d1_;
  r_t d3_;
  j_t d4_;
  gradient_t d5_;
  hessian_t d6_;
  hessian_t d7_;
  lm_damp_t d8_;
  InnSolverType * d9_;
  SystemType const * d10_;

public:
  RegistryLevMarNormalEqs(const SystemType & system, InnSolverType & innS)
    : d1_(system.createState()),
      d3_(system.createResidual()),
      d4_(system.createJacobian()),
      d5_(system.createState()),
      d6_( hg_default::createHessian(system.createState()) ),
      d7_( hg_default::createHessian(system.createState()) ),
      d8_{},
      d9_(&innS),
      d10_(&system){}

  RegistryLevMarNormalEqs() = delete;
  // non-copyable
  RegistryLevMarNormalEqs(RegistryLevMarNormalEqs const &) = delete;
  RegistryLevMarNormalEqs& operator=(RegistryLevMarNormalEqs const &) = delete;
  // movable
  RegistryLevMarNormalEqs(RegistryLevMarNormalEqs &&) = default;
  RegistryLevMarNormalEqs& operator=(RegistryLevMarNormalEqs &&) = default;

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10>::value) < 9;
  }

  GETMETHOD(1)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
  GETMETHOD(9)
  GETMETHOD(10)
};

}}}
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_REGISTRIES_HPP_
