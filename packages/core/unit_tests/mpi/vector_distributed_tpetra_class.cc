
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraVectorGlobSize15Fixture, Constructor){
  using namespace rompp;

  std::cout << rank_ << std::endl;
  using myvec_t = core::Vector<typename tpetraVectorGlobSize15Fixture::vec_t>;
  myvec_t a( *x_ );
  myvec_t a3( a );
}


TEST_F(tpetraVectorGlobSize15Fixture,
       localSize){
  using namespace rompp;
  using myvec_t = core::Vector<typename tpetraVectorGlobSize15Fixture::vec_t>;
  myvec_t v1( *x_ );
  EXPECT_EQ( v1.localSize(), 5);
}


TEST_F(tpetraVectorGlobSize15Fixture,
       globalSize){
  using namespace rompp;
  using myvec_t = core::Vector<typename tpetraVectorGlobSize15Fixture::vec_t>;
  myvec_t v1( *x_ );
  EXPECT_EQ( v1.globalSize(), 15);
}


TEST_F(tpetraVectorGlobSize15Fixture,
       isGloballyDist){
  using namespace rompp;
  using myvec_t = core::Vector<typename tpetraVectorGlobSize15Fixture::vec_t>;
  myvec_t v1( *x_ );
  EXPECT_TRUE( v1.isDistributedGlobally() );
}


TEST_F(tpetraVectorGlobSize15Fixture,
       SetScalar){
  using namespace rompp;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = core::Vector<typename tpetraVectorGlobSize15Fixture::vec_t>;
  myvec_t v1( *x_ );
  v1.putScalar(43.3);
  
  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 43.3 );
  }
}


TEST_F(tpetraVectorGlobSize15Fixture,
       SetZero){
  using namespace rompp;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = core::Vector<typename tpetraVectorGlobSize15Fixture::vec_t>;
  myvec_t v1( *x_ );
  v1.setZero();
  
  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 0.0 );
  }
}


TEST_F(tpetraVectorGlobSize15Fixture,
       QueryWrappedData){
  using namespace rompp;
  using nvec_t = typename tpetraVectorGlobSize15Fixture::vec_t;
  
  using myvec_t = core::Vector<nvec_t>;
  myvec_t v1( *x_ );
  ::testing::StaticAssertTypeEq<decltype(v1.data()),
  				nvec_t * >(); 
  const myvec_t v2( *x_ );
  ::testing::StaticAssertTypeEq< decltype(v2.data()),
  				 const nvec_t * >();
}


TEST_F(tpetraVectorGlobSize15Fixture,
       empty){
  using namespace rompp;
  using nvec_t = typename tpetraVectorGlobSize15Fixture::vec_t;  
  using myvec_t = core::Vector<nvec_t>;
  myvec_t v1( *x_ );
  EXPECT_FALSE(v1.empty());
}


TEST_F(tpetraVectorGlobSize15Fixture,
       getMap){
  using namespace rompp;
  using nvec_t = typename tpetraVectorGlobSize15Fixture::vec_t;  
  using myvec_t = core::Vector<nvec_t>;
  myvec_t v1( *x_ );
  auto const & mapO = v1.getDataMap();

  ::testing::StaticAssertTypeEq<decltype(mapO),
  				const typename tpetraVectorGlobSize15Fixture::map_t & >(); 
  EXPECT_TRUE(mapO.isContiguous());
}



// TEST_F(tpetraVectorGlobSize15Fixture,
//        SubscriptOperatorConstSqBrack){
//   using namespace rompp;
//   using nvec_t = typename tpetraVectorGlobSize15Fixture::vec_t;  
//   using myvec_t = core::Vector<nvec_t>;
//   myvec_t v1( *x_ );
//   v1.putScalar(11.2);
  
//   const myvec_t v2(v1);
//   for (int i=0; i<v2.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v2[i], 11.2 );
//   }
// }


// TEST_F(tpetraVectorGlobSize15Fixture,
//        SubscriptOperatorConstParenth){
//   using namespace rompp;
//   using nvec_t = typename tpetraVectorGlobSize15Fixture::vec_t;  
//   using myvec_t = core::Vector<nvec_t>;
//   myvec_t v1( *x_ );
//   v1.putScalar(11.2);
  
//   const myvec_t v2(v1);
//   for (int i=0; i<v2.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v2(i), 11.2 );
//   }
// }


// TEST_F(tpetraVectorGlobSize15Fixture,
//        SubscriptOperatorNonConstSqBrack){
//   using namespace rompp;
//   using nvec_t = typename tpetraVectorGlobSize15Fixture::vec_t;  
//   using myvec_t = core::Vector<nvec_t>;
//   myvec_t v1( *x_ );
//   for (int i=0; i<v1.localSize(); i++){
//     v1[i] = 11.2;
//   }
//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], 11.2 );
//   }
// }


// TEST_F(tpetraVectorGlobSize15Fixture,
//        SubscriptOperatorNonConstParenth){
//   using namespace rompp;
//   using nvec_t = typename tpetraVectorGlobSize15Fixture::vec_t;  
//   using myvec_t = core::Vector<nvec_t>;
//   myvec_t v1( *x_ );
//   for (int i=0; i<v1.localSize(); i++){
//     v1(i) = 11.2;
//   }
//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1(i), 11.2 );
//   }
// }


// TEST_F(tpetraVectorGlobSize15Fixture,
//        compAssignPlus){
// This gives a linking issue for daxpy in kokkos:blas
//   using namespace rompp;
//   using nvec_t = typename tpetraVectorGlobSize15Fixture::vec_t;  
//   using myvec_t = core::Vector<nvec_t>;

//   myvec_t v1( *x_ );
//   v1.putScalar(1.2);
//   myvec_t v2( v1 );
//   v2.putScalar(2.);

//   v1 += v2;
//   for (int i=0; i<v1.localSize(); i++){
//     v1[i] = 3.2;
//   }
// }
