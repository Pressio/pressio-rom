
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "Epetra_MpiComm.h"

struct core_matrix_vec_product_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;

  core::Matrix<Epetra_CrsMatrix> * A_;
  core::Vector<Epetra_Vector> * b_;
  int nRowsSM_;
  Epetra_Map * smMap_;
  int nRowsV_;
  Epetra_Map * vMap_;

  core::Matrix<Epetra_MultiVector> * Ad_;
  core::Vector<Epetra_Vector> * b2_;
  int nRowsDM_;
  Epetra_Map * dmMap_;
  Epetra_Map * b2Map_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    assert( NumProc_ == 3 );

    // setup CRS and vector for test 1
    nRowsSM_ = 9;
    smMap_ = new Epetra_Map(nRowsSM_, 0, *Comm_);
    A_ = new core::Matrix<Epetra_CrsMatrix>(*smMap_, 5);
    nRowsV_ = 13;
    vMap_ = new Epetra_Map(nRowsV_, 0, *Comm_);
    vMap_->Print(std::cout);
    b_ = new core::Vector<Epetra_Vector>(*vMap_);
    //------------------------------------------

    // setup dense and vector for test 2
    nRowsDM_ = 9;
    dmMap_ = new Epetra_Map(nRowsDM_, 0, *Comm_);
    Ad_ = new core::Matrix<Epetra_MultiVector>(*dmMap_, 5);
    b2Map_ = new Epetra_Map(5, 0, *Comm_);
    b2Map_->Print(std::cout);
    b2_ = new core::Vector<Epetra_Vector>(*b2Map_);
    //------------------------------------------    
  }
  
  virtual void TearDown(){
    delete Comm_;
    delete vMap_;
    delete smMap_;
    delete A_;
    delete b_;
    delete b2Map_;
    delete dmMap_;
    delete Ad_;
    delete b2_;
  }
};



TEST_F(core_matrix_vec_product_distributed_epetraFix, CRSMatTimesVector)
{

  //-----------
  // FILL A
  //-----------
  {
    int myNR = smMap_->NumMyElements();

    std::vector<int> mygid(myNR);
    smMap_->MyGlobalElements( mygid.data() );
    std::array<double,5> vals;
    std::array<int,5> colind;
    for (auto const & it : mygid){
      if(it == 0){
      	vals = {1., 1., 1., 3., 3.};
      	colind = {0, 4, 7, 11, 12};
      	A_->insertGlobalValues(it, 5, vals.data(), colind.data());
      }
      if(it == 6){
      	vals = {1., 3., 6.};
      	colind = {1, 4, 7};
      	A_->insertGlobalValues(it, 3, vals.data(), colind.data());
      }

    }
    A_->fillingIsCompleted(*vMap_, *smMap_);
    A_->data()->Print(std::cout);
  }

  //-----------
  // FILL b
  //-----------
  {
    if (MyPID_ == 0){
      (*b_)[0] = 1.;
      (*b_)[4] = 3.;
    }
    if (MyPID_ == 1){
      (*b_)[2] = 1.;
    }
    if (MyPID_ == 2){
      (*b_)[2] = 2.;
      (*b_)[3] = 2.;
    }
    b_->data()->Print(std::cout);
  }
  
  //------------------
  // product: b = A b
  //------------------
  auto c = core::mat_ops::product(*A_, *b_);
  c.data()->Print(std::cout);

  assert( c.globalSize() == 9);
  static_assert( std::is_same<decltype(c),
		 core::Vector<Epetra_Vector>>::value, "" );
  if (MyPID_ == 0){
    EXPECT_DOUBLE_EQ( c[0], 17. );
  }
  if (MyPID_ == 2)
    EXPECT_DOUBLE_EQ( c[0], 15. );
  
}
///////////////////////////////////////////
///////////////////////////////////////////



TEST_F(core_matrix_vec_product_distributed_epetraFix, DenseMatTimesVector)
{

  //-----------
  // FILL A
  //-----------
  {
    if (MyPID_ == 0){
      (*Ad_)(0,0) = 1.0;
      (*Ad_)(0,2) = 2.0;
      (*Ad_)(0,3) = 3.0;
      (*Ad_)(0,4) = 1.0;
    }
    if (MyPID_ == 1){
      (*Ad_)(1,0) = 3.0;
      (*Ad_)(1,2) = 2.0;
      (*Ad_)(1,3) = 3.0;
      (*Ad_)(1,4) = 4.0;
    }

    Ad_->data()->Print(std::cout);
  }

  //-----------
  // FILL b
  //-----------
  {
    if (MyPID_ == 0){
      (*b2_)[0] = 1.;
      (*b2_)[1] = 2.;
    }
    if (MyPID_ == 1){
      (*b2_)[0] = 1.;
    }
    if (MyPID_ == 2){
      (*b2_)[0] = 1.;
    }
    b2_->data()->Print(std::cout);
  }
  
  //------------------
  // product: b = A b
  //------------------
  auto c = core::mat_ops::product( *Ad_, *b2_ );  
  c.data()->Print(std::cout);

  assert( c.globalSize() == 9);
  static_assert( std::is_same<decltype(c),
  		 core::Vector<Epetra_Vector>>::value, "" );
  if (MyPID_ == 0){
    EXPECT_DOUBLE_EQ( c[0], 4. );
  }
  if (MyPID_ == 1)
    EXPECT_DOUBLE_EQ( c[1], 9. );
  
}
