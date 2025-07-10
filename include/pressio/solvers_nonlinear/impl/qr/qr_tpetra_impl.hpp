/*
//@HEADER
// ************************************************************************
//
// qr_tpetra_impl.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_TPETRA_IMPL_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_TPETRA_IMPL_HPP_

#include "Tpetra_TsqrAdaptor.hpp"
#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

namespace pressio{ namespace qr{ namespace impl{


#ifdef PRESSIO_ENABLE_TPL_EIGEN
template< typename VectorType, typename R_type>
std::enable_if_t<
  ::pressio::is_vector_eigen<VectorType>::value and
  ::pressio::is_dense_matrix_eigen<R_type>::value
>
solve(const VectorType & rhs,
     const R_type & Rmatrix,
     VectorType & y)
{
  y = Rmatrix.template triangularView<Eigen::Upper>().solve(rhs);
}
#endif

template<typename VectorType, typename R_type>
std::enable_if_t<
  ::pressio::is_dense_vector_teuchos<VectorType>::value and
  ::pressio::is_dense_matrix_teuchos_rcp<R_type>::value
>
solve(const VectorType & rhs, R_type Rmatrix, VectorType & y)
{
  using ord_t  = typename R_type::element_type::ordinalType;
  using sc_t   = typename R_type::element_type::scalarType;
  using solver_t = Teuchos::SerialDenseSolver<ord_t, sc_t>;

  solver_t My_Solver;
  My_Solver.setMatrix(Rmatrix);
  My_Solver.setVectors( Teuchos::rcpFromRef(y),
    Teuchos::rcpFromRef(const_cast<VectorType&>(rhs)));

  int info = My_Solver.factor();
  if (info != 0){
    std::cout << "SerialDenseSolver::factor() returned : "
        << info << std::endl;
  }

  info = My_Solver.solve();
  if (info != 0){
    std::cout << "SerialDenseSolver::solve() returned : "
        << info << std::endl;
  }
}

template<typename VectorType, typename R_type>
std::enable_if_t<
  !::pressio::is_dense_vector_teuchos<VectorType>::value and
  ::pressio::is_dense_matrix_teuchos_rcp<R_type>::value
>
solve(const VectorType & rhs, R_type Rmatrix, VectorType & y)
{
  // todo: maybe here we should use directly backward substitution but it does
  // not seem to be available directly from teuchos

  using ord_t   = typename R_type::element_type::ordinalType;
  using sc_t    = typename R_type::element_type::scalarType;
  using tservec_t = Teuchos::SerialDenseVector<ord_t, sc_t>;

  auto vecSize = ::pressio::ops::extent(rhs,0);
  tservec_t rhsTV(Teuchos::View, const_cast<sc_t*>(rhs.data()), vecSize);
  tservec_t yTV(Teuchos::View, y.data(), vecSize);

  Teuchos::SerialQRDenseSolver<ord_t, sc_t> My_Solver;
  My_Solver.setMatrix(Rmatrix);
  My_Solver.setVectors(Teuchos::rcpFromRef(yTV), Teuchos::rcpFromRef(rhsTV));
  int info = My_Solver.factor();
  if (info != 0){
    std::cout << "SerialDenseSolver::factor() returned : " << info << std::endl;
  }

  info = My_Solver.solve();
  if (info != 0){
    std::cout << "SerialDenseSolver::solve() returned : " << info << std::endl;
  }
}

template<typename MatrixType, typename R_t>
class TpetraBlockMVTSQR
{

public:
  using int_t      = int;
  using sc_t       = typename ::pressio::Traits<MatrixType>::scalar_type;
  using serden_mat_t = Teuchos::SerialDenseMatrix<int_t, sc_t>;
  using trcp_mat     = Teuchos::RCP<serden_mat_t>;

  using lo_t     = typename MatrixType::local_ordinal_type;
  using go_t     = typename MatrixType::global_ordinal_type;
  using node_t   = typename MatrixType::node_type;

  using Q_type            = Tpetra::BlockMultiVector<sc_t, lo_t, go_t, node_t>;
  using tpetra_mv_t       = Tpetra::MultiVector<sc_t, lo_t, go_t, node_t>;
  using tsqr_adaptor_type = Tpetra::TsqrAdaptor<tpetra_mv_t>;

public:
  TpetraBlockMVTSQR() = default;
  ~TpetraBlockMVTSQR() = default;

  void computeThinOutOfPlace(const MatrixType & Ain)
  {
    auto & A = const_cast<MatrixType &>(Ain);

    auto nVecs     = ::pressio::ops::extent(A,1);
    auto blockSize = A.getBlockSize();
    createLocalRIfNeeded(nVecs);

    // this is the row map of the block MV
    auto ArowMap = A.getMap();
    createQIfNeeded(*ArowMap, blockSize, nVecs);

    // get the multivector
    auto mv = A.getMultiVectorView();
    auto Qv = Qmat_->getMultiVectorView();
    tsqrAdaptor_.factorExplicit(mv, Qv, *localR_.get(), false);
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const
  {
    qr::impl::solve<VectorType, trcp_mat>(rhs, this->localR_, y);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyQTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const
  {
    constexpr auto beta  = static_cast<sc_t>(0);
    constexpr auto alpha = static_cast<sc_t>(1);
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyRTranspose(const VectorInType & vecIn, VectorOutType & y) const
  {
    constexpr auto beta  = static_cast<sc_t>(0);
    constexpr auto alpha = static_cast<sc_t>(1);
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->localR_, vecIn, beta, y);
  }

  template <typename T = R_t>
  std::enable_if_t<
    !::pressio::is_dense_matrix_teuchos<T>::value and
    !std::is_void<T>::value,
    const T &
  >
  RFactor() const
  {
    this->Rmat_ = std::make_shared<T>(this->localR_->values());
    return *this->Rmat_;
  }

  template <typename T = R_t>
  std::enable_if_t<
    ::pressio::is_dense_matrix_teuchos<T>::value and
    !std::is_void<T>::value,
    const T &
  >
  RFactor() const
  {
    this->Rmat_ = std::make_shared<T>(*this->localR_, Teuchos::View);
    return *this->Rmat_;
  }

  const Q_type & QFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(int newsize)
  {
    if (localR_.is_null() or
      (localR_->numRows()!=newsize and localR_->numCols()!=newsize)){
      localR_ = Teuchos::rcp(new serden_mat_t(newsize, newsize) );
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int blockSize, int numVecs)
  {
    if (!Qmat_ or !Qmat_->getMap()->isSameAs(map) ){
     Qmat_ = std::make_shared<Q_type>(map, blockSize, numVecs);
    }
  }

private:
  tsqr_adaptor_type tsqrAdaptor_;
  trcp_mat localR_      = {};
  int computedRank_     = {};
  mutable std::shared_ptr<Q_type> Qmat_ = nullptr;
  mutable std::shared_ptr<R_t> Rmat_  = nullptr;
};


template<typename MatrixType, typename R_t>
class TpetraMVHouseholderUsingEigen
{

public:
  using sc_t = typename ::pressio::Traits<MatrixType>::scalar_type;
  using lo_t         = typename MatrixType::local_ordinal_type;
  using go_t         = typename MatrixType::global_ordinal_type;
  using node_t       = typename MatrixType::node_type;
  using Q_type       = Tpetra::MultiVector<sc_t, lo_t, go_t, node_t>;

  using eig_dyn_mat = Eigen::Matrix<sc_t, -1, -1>;
  using help_impl_t = QRHouseholderDenseEigenMatrix<eig_dyn_mat, R_t>;
  help_impl_t myImpl_ = {};

public:
  TpetraMVHouseholderUsingEigen() = default;
  ~TpetraMVHouseholderUsingEigen() = default;

  template < typename vector_in_t, typename vector_out_t>
  void applyQTranspose(const vector_in_t & vecIn, vector_out_t & vecOut) const
  {
    constexpr auto beta  = static_cast<sc_t>(0);
    constexpr auto alpha = static_cast<sc_t>(1);
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const{
    myImpl_.template doLinSolve<VectorType>(rhs, y);
  }

  const Q_type & QFactor() const {
    return *this->Qmat_;
  }

  void computeThinOutOfPlace(const MatrixType & Ain)
  {
    auto & A = const_cast<MatrixType &>(Ain);

    size_t rows = ::pressio::ops::extent(A,0);
    size_t cols = ::pressio::ops::extent(A,1);
    auto ArowMap = A.getMap();
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
      Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));

    // convert it to replicated eptra matrix
    using local_map_t = Tpetra::Map<lo_t, go_t, node_t>;
    using rcp_local_map_t = Teuchos::RCP<const local_map_t>;
    rcp_local_map_t rcp_local_map = Teuchos::rcp( new local_map_t(rows, 0, comm) );

    using import_t = Tpetra::Import<lo_t, go_t, node_t>;
    import_t importer(ArowMap, rcp_local_map);
    MatrixType A2(rcp_local_map, cols);
    A2.doImport(A, importer, Tpetra::INSERT);

    // store it into an Eigen matrix
    Eigen::Matrix<sc_t, -1, -1> eA2W(rows,cols);
    for (size_t j=0;j<cols;j++)
    {
      auto colData = A2.getData(j);
      for (size_t i=0;i<rows;i++){
       eA2W(i,j) = colData[i];
      }
    }

    myImpl_.computeThinOutOfPlace(eA2W);

    // store Q into replicated Tpetra::Multivector
    const auto & Q2 = myImpl_.QFactor();
    Q_type locQ( rcp_local_map, Q2.cols() );
    // // auto trilD = locQ.data();
    // locQ.template sync<Kokkos::HostSpace>();

    auto v2d = locQ.template getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWriteStruct());
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    // // //we are going to change the host view
    // locQ.template modify<Kokkos::HostSpace>();
    for (size_t i=0;i< (size_t)Q2.rows();i++)
      for (size_t j=0;j< (size_t)Q2.cols();j++)
      v2d(i,j) = Q2(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_type>(ArowMap, Q2.cols());
    import_t importer2(rcp_local_map, ArowMap);
    Qmat_->doImport(locQ, importer2, Tpetra::INSERT);
  }

private:
  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_type> Qmat_ = nullptr;
  mutable std::shared_ptr<R_t> Rmat_  = nullptr;
};



template<typename MatrixType, typename R_t>
class TpetraMVTSQR
{

public:
  using int_t        = int;
  using sc_t         = typename ::pressio::Traits<MatrixType>::scalar_type;
  using serden_mat_t = Teuchos::SerialDenseMatrix<int_t, sc_t>;
  using trcp_mat     = Teuchos::RCP<serden_mat_t>;

  using lo_t         = typename MatrixType::local_ordinal_type;
  using go_t         = typename MatrixType::global_ordinal_type;
  using node_t       = typename MatrixType::node_type;
  using Q_type       = Tpetra::MultiVector<sc_t, lo_t, go_t, node_t>;
  using tsqr_adaptor_type = Tpetra::TsqrAdaptor<MatrixType>;

public:
  TpetraMVTSQR() = default;
  ~TpetraMVTSQR() = default;

  void computeThinOutOfPlace(const MatrixType & A)
  {
    auto nVecs = ::pressio::ops::extent(A,1);
    auto ArowMap = A.getMap();
    createQIfNeeded(ArowMap, nVecs);
    createLocalRIfNeeded(nVecs);
    tsqrAdaptor_.factorExplicit(const_cast<MatrixType &>(A), *Qmat_, *localR_.get(), false);
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const {
      qr::impl::solve<VectorType, trcp_mat>(rhs, this->localR_, y);
  }


  template < typename VectorInType, typename VectorOutType>
  void applyQTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const
  {
    constexpr auto beta  = static_cast<sc_t>(0);
    constexpr auto alpha = static_cast<sc_t>(1);
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyRTranspose(const VectorInType & vecIn, VectorOutType & y) const
  {
    constexpr auto beta  = static_cast<sc_t>(0);
    constexpr auto alpha = static_cast<sc_t>(1);
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->localR_, vecIn, beta, y);
  }

  template <typename T = R_t>
  std::enable_if_t<
    !::pressio::is_dense_matrix_teuchos<T>::value and
    !std::is_void<T>::value,
    const T &
  >
  RFactor() const {
    this->Rmat_ = std::make_shared<T>(this->localR_->values());
    return *this->Rmat_;
  }

  template <typename T = R_t>
  std::enable_if_t<
    ::pressio::is_dense_matrix_teuchos<T>::value and
    !std::is_void<T>::value,
    const T &
  >
  RFactor() const {
    this->Rmat_ = std::make_shared<T>(*this->localR_, Teuchos::View);
    return *this->Rmat_;
  }

  const Q_type & QFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(int newsize)
  {
    if (localR_.is_null() or
      (localR_->numRows()!=newsize and localR_->numCols()!=newsize)){
      localR_ = Teuchos::rcp(new serden_mat_t(newsize, newsize) );
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int cols)
  {
    if (!Qmat_ or !Qmat_->getMap()->isSameAs(*map) ){
      Qmat_ = std::make_shared<Q_type>(map, cols);
    }
  }

private:
  tsqr_adaptor_type tsqrAdaptor_;
  trcp_mat localR_      = {};
  int computedRank_     = {};
  mutable std::shared_ptr<Q_type> Qmat_ = nullptr;
  mutable std::shared_ptr<R_t> Rmat_  = nullptr;
};




template<typename MatrixType, typename R_t>
class ModGramSchmidtMVTpetra
{

public:
  using int_t      = int;
  using sc_t        = typename ::pressio::Traits<MatrixType>::scalar_type;
  using lo_t        = typename MatrixType::local_ordinal_type;
  using go_t        = typename MatrixType::global_ordinal_type;
  using node_t      = typename MatrixType::node_type;
  using Q_type      = Tpetra::MultiVector<sc_t, lo_t, go_t, node_t>;
  using R_nat_t     = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic>;

public:
  ModGramSchmidtMVTpetra() = default;
  ~ModGramSchmidtMVTpetra() = default;

  void computeThinOutOfPlace(const MatrixType & Ain)
  {
    auto & A = const_cast<MatrixType &>(Ain);

    size_t nVecs = ::pressio::ops::extent(A,1);
    auto ArowMap = A.getMap();
    createQIfNeeded(ArowMap, nVecs);
    createLocalRIfNeeded(nVecs);

    sc_t rkkInv = {};
    for (size_t k=0; k<nVecs; k++)
    {
      auto ak = A.getVector(k);
      localR_(k,k) = ak->norm2();
      rkkInv = static_cast<sc_t>(1)/localR_(k,k);

      auto qk = Qmat_->getVectorNonConst(k);
      qk->update( rkkInv, *ak, static_cast<sc_t>(0) );

      for (size_t j=k+1; j<nVecs; j++){
        auto aj = A.getVectorNonConst(j);
        localR_(k,j) = qk->dot(*aj);
        aj->update(-localR_(k,j), *qk, static_cast<sc_t>(1));
      }
    }
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const {
    //auto vecSize = y.size();
    auto & Rm = localR_.template triangularView<Eigen::Upper>();
    y = Rm.solve(rhs);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyQTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const
  {
    constexpr auto beta  = static_cast<sc_t>(0);
    constexpr auto alpha = static_cast<sc_t>(1);
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  const Q_type & QFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(std::size_t newsize){
    const std::size_t locRext0 = ::pressio::ops::extent(localR_, 0);
    const std::size_t locRext1 = ::pressio::ops::extent(localR_, 1);
    if (locRext0!=newsize or locRext1!=newsize){
      localR_ = R_nat_t(newsize, newsize);
      ::pressio::ops::set_zero(localR_);
    }
  }

  template <typename MapType>
  void createQIfNeeded(const MapType & map, std::size_t cols){
    if (!Qmat_ or !Qmat_->getMap()->isSameAs(*map) )
      Qmat_ = std::make_shared<Q_type>(map, cols);
  }

private:
  R_nat_t localR_     = {};
  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_type> Qmat_ = nullptr;
  mutable std::shared_ptr<R_t> Rmat_  = nullptr;
};

}}} // end namespace pressio::qr::impl
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_TPETRA_IMPL_HPP_
