
#ifndef PRESSIOROM_ROM_LINEAR_ROM_HPP_
#define PRESSIOROM_ROM_LINEAR_ROM_HPP_

#include <fstream>

namespace pressio{ namespace rom{ namespace linear{

namespace impl{

template<class ScalarType>
struct LinearRomDefaultReducedOperatorsTraits
{
  using reduced_matrix_type = Eigen::Matrix<ScalarType, -1, -1>;
  using reduced_vector_type = Eigen::Matrix<ScalarType, -1, 1>;

  static auto createMatrix(std::size_t nr, std::size_t nc){
    return reduced_matrix_type(nr, nc);
  }

  static auto createVector(std::size_t n){
    return reduced_vector_type(n);
  }
};

template <typename T>
void write_matrix_to_ascii(const std::string fileName,
			   const Eigen::DenseBase<T> & A,
			   bool writeSize = false)
{
  const int m = A.rows();
  const int n = A.cols();
  std::ofstream file; file.open(fileName);
  if (writeSize){
    file << m << " " << n << std::endl;
  }

  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      file << A(i,j) << " ";
    }
    file << std::endl;
  }
  file.close();
}

template <
  class FomMaxOperatorType, class LeftBasisType,
  class RightBasisType, class ShiftVecType>
class RomLinearThingy;

template <class FomVecType, class LeftBasisType>
class RomLinearThingyConstVector;
}

template <class ...Ts>
void export_ascii(const impl::RomLinearThingy<Ts...> & o,
		  const std::string & matFile);

template <class ...Ts>
void export_ascii(const impl::RomLinearThingy<Ts...> & o,
		  const std::string & matFile,
		  const std::string & shiftFile);

template <class ...Ts>
void export_ascii(const impl::RomLinearThingy<Ts...> & o,
		  const std::string & matFile,
		  const std::string & shiftFile,
		  const std::string & vecFile);

template <class ...Ts>
void export_ascii(const impl::RomLinearThingyConstVector<Ts...> & o,
		  const std::string & fileout);

namespace impl{

template <
  class FomMaxOperatorType, class LeftBasisType,
  class RightBasisType, class ShiftVecType>
class RomLinearThingy{

  // we need to make this friend to access private methods
  template <class ...Ts> friend
  void ::pressio::rom::linear::export_ascii(const RomLinearThingy<Ts...> &,
		    const std::string &);

  template <class ...Ts> friend
  void ::pressio::rom::linear::export_ascii(const RomLinearThingy<Ts...> &,
		    const std::string &,
		    const std::string &);

  // we need to make this friend to access private methods
  template <class ...Ts> friend
  void ::pressio::rom::linear::export_ascii(const RomLinearThingy<Ts...> &,
		    const std::string &,
		    const std::string &,
		    const std::string &);


  using fomop_type = pressio::mpl::remove_cvref_t<FomMaxOperatorType>;
  using lb_type    = pressio::mpl::remove_cvref_t<LeftBasisType>;
  using rb_type    = pressio::mpl::remove_cvref_t<RightBasisType>;
  using shift_type = pressio::mpl::remove_cvref_t<ShiftVecType>;
  using reduced_operators = LinearRomDefaultReducedOperatorsTraits<
    typename pressio::Traits<lb_type>::scalar_type>;

  const fomop_type * fop_ = nullptr;
  const lb_type * lb_ = nullptr;
  const rb_type * rb_ = nullptr;
  const shift_type * shift_ = nullptr;

  rb_type fomATimesRightBasis_;  // fomATimesRightBasis_ = fomA * rb
  shift_type fomATimeShift_;     // fomATimeShift_ = fomA * shift

  // this is the hiddentVector when the fom is a linear operator that evaluates as:
  //    result = fomA * operand + hiddenVector
  shift_type fomHiddenForcing_;

  typename reduced_operators::reduced_matrix_type redMat_;
  typename reduced_operators::reduced_vector_type redShift_;
  typename reduced_operators::reduced_vector_type redForcing_;
  bool solve_;

public:
  explicit RomLinearThingy(const fomop_type & fop, const lb_type & lb,
			   const rb_type & rb, const shift_type & shift,
			   bool solve)
    : fop_(&fop), lb_(&lb), rb_(&rb), shift_(&shift)
    , fomATimesRightBasis_(fop.createResultOfActionOn(*rb_))
    , fomATimeShift_(fop.createResultOfActionOn(*shift_))
    , fomHiddenForcing_(pressio::ops::clone(*shift_))
    , redMat_(reduced_operators::createMatrix(pressio::ops::extent(*lb_, 1),
					      pressio::ops::extent(*rb_, 1)) )
    , redShift_(pressio::ops::extent(*lb_, 1))
    , redForcing_(pressio::ops::extent(*lb_, 1))
    , solve_(solve)
  {
    if (solve_){
      this->compute2();
    }
    else{
      this->compute();
    }
  }

  auto operator()() const{
    // return tie bc we are return referencs to lvalues
    return std::tie(std::as_const(redMat_),
		    std::as_const(redShift_),
		    std::as_const(redForcing_));
  }

private:
  void compute2()
  {
    constexpr auto pT  = ::pressio::transpose();
    constexpr auto pNT = ::pressio::nontranspose();

    // in this case, the fom apply method on a vector (v)
    // computes:   result = fomA * v + hiddenVector
    // so we don't know how to split these yet, therefore we can do as follows

    // 1.
    auto zerovec = pressio::ops::clone(*shift_);
    pressio::ops::fill(zerovec, 0);
    fop_->apply(zerovec, fomHiddenForcing_);
    pressio::ops::product(pT, 1., *lb_, fomHiddenForcing_, 0., redForcing_);

    // 2.
    const std::size_t nCols = pressio::ops::extent(*rb_, 1);
    for (std::size_t j=0; j<nCols; ++j){
      auto rightBasisCol_j = pressio::column(*rb_, j);
      auto currCol = pressio::column(fomATimesRightBasis_, j);
      fop_->apply(rightBasisCol_j, currCol);
      // remove the forcing
      pressio::ops::update(currCol, 1., fomHiddenForcing_, -1.);
    }
    pressio::ops::product(pT, pNT, 1., *lb_, fomATimesRightBasis_, 0., redMat_);

    // 3.
    fop_->apply(*shift_, fomATimeShift_);
    pressio::ops::update(fomATimeShift_, 1., fomHiddenForcing_, -1.);
    pressio::ops::product(pT, 1., *lb_, fomATimeShift_, 0., redShift_);
  }

  void compute()
  {
    const std::size_t nCols = pressio::ops::extent(*rb_, 1);
    for (std::size_t j=0; j<nCols; ++j){
      auto rightBasisCol_j = pressio::column(*rb_, j);
      auto currCol = pressio::column(fomATimesRightBasis_, j);
      fop_->apply(rightBasisCol_j, currCol);
    }

    constexpr auto pT  = ::pressio::transpose();
    constexpr auto pNT = ::pressio::nontranspose();
    pressio::ops::product(pT, pNT, 1., *lb_, fomATimesRightBasis_, 0., redMat_);

    fop_->apply(*shift_, fomATimeShift_);
    pressio::ops::product(pT, 1., *lb_, fomATimeShift_, 0., redShift_);
  }
};

template <class FomVecType, class LeftBasisType>
class RomLinearThingyConstVector
{

  // we need to make this friend to access private methods
  template <class ...Ts> friend
  void ::pressio::rom::linear::export_ascii(const RomLinearThingyConstVector<Ts...> &,
					    const std::string &);

  using fomop_type = pressio::mpl::remove_cvref_t<FomVecType>;
  using lb_type    = pressio::mpl::remove_cvref_t<LeftBasisType>;
  using reduced_operators = LinearRomDefaultReducedOperatorsTraits<
    typename pressio::Traits<lb_type>::scalar_type>;

  const fomop_type * fop_ = nullptr;
  const lb_type * lb_ = nullptr;
  typename reduced_operators::reduced_vector_type result_;

public:
  explicit RomLinearThingyConstVector(const fomop_type & fop,
				      const lb_type & lb)
    : fop_(&fop), lb_(&lb)
    , result_(reduced_operators::createVector(pressio::ops::extent(*lb_, 1)) )
  {
    this->compute();
  }

  const auto & operator()() const{
    return std::as_const(result_);
    // // return tie bc we are return referencs to lvalues
    // return std::tie(std::as_const(result_));
  }

private:
  void compute()
  {
    const auto fomF = fop_->evaluate();
    constexpr auto pT = ::pressio::transpose();
    pressio::ops::product(pT, 1., *lb_, fomF, 0., result_);
  }
};

} //end anonym namespace

template <
  class FomMatOperatorType,
  class LeftBasisType,
  class RightBasisType,
  class ShiftType
  >
auto create_reduced_matrix_operator(const FomMatOperatorType & fomOp,
				    const LeftBasisType & lb,
				    const RightBasisType & rb,
				    const ShiftType & shift)
{
  using ret_t = impl::RomLinearThingy<
    FomMatOperatorType, LeftBasisType, RightBasisType, ShiftType>;
  return ret_t(fomOp, lb, rb, shift, false);
}

template <class FomConstVecOperatorType, class LeftBasisType>
auto create_reduced_vector_operator(const FomConstVecOperatorType & fomOp,
				    const LeftBasisType & lb)
{
  using ret_t = impl::RomLinearThingyConstVector<
    FomConstVecOperatorType, LeftBasisType>;
  return ret_t(fomOp, lb);
}

template <
  class FomLinearOperatorType,
  class LeftBasisType,
  class RightBasisType,
  class ShiftType
  >
auto create_reduced_linear_operator(const FomLinearOperatorType & fomOp,
				    const LeftBasisType & lb,
				    const RightBasisType & rb,
				    const ShiftType & shift)
{
  using ret_t = impl::RomLinearThingy<
    FomLinearOperatorType, LeftBasisType, RightBasisType, ShiftType>;
  return ret_t(fomOp, lb, rb, shift, true);
}

template <class ...Ts>
void export_ascii(const impl::RomLinearThingy<Ts...> & o,
		  const std::string & matFile)
{
  impl::write_matrix_to_ascii(matFile, o.redMat_);
}

template <class ...Ts>
void export_ascii(const impl::RomLinearThingy<Ts...> & o,
		  const std::string & matFile,
		  const std::string & shiftFile)
{
  impl::write_matrix_to_ascii(matFile, o.redMat_);
  impl::write_matrix_to_ascii(shiftFile, o.redShift_);
}

template <class ...Ts>
void export_ascii(const impl::RomLinearThingy<Ts...> & o,
		  const std::string & matFile,
		  const std::string & shiftFile,
		  const std::string & vecFile)
{
  impl::write_matrix_to_ascii(matFile, o.redMat_);
  impl::write_matrix_to_ascii(shiftFile, o.redShift_);
  impl::write_matrix_to_ascii(vecFile, o.redForcing_);
}

template <class ...Ts>
void export_ascii(const impl::RomLinearThingyConstVector<Ts...> & o,
		  const std::string & fileout)
{
  impl::write_matrix_to_ascii(fileout, o.result_);
}

}}}
#endif  // PRESSIOROM_ROM_LINEAR_ROM_HPP_
