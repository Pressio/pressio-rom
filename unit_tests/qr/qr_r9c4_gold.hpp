
#ifndef QR_GOLD_SOLUTIONS_HPP_
#define QR_GOLD_SOLUTIONS_HPP_

#include "Eigen/Dense"

namespace pressio{ namespace qr{ namespace test{

template <typename T = double>
struct qrGoldr9c4Sol{

  Eigen::Matrix<T, 9, 9> trueQ_;
  Eigen::Matrix<T, 4, 4> trueR_;

  Eigen::Matrix<T, 9, 1> colDotOnes_;
  Eigen::Matrix<T, 4, 1> trueYForRSolve_;

  qrGoldr9c4Sol(){
    fillTrueQ();
    fillTrueR();
    fillTrueYForRSolve();
    // compute dot products between each column of Q and a vector of ones
    colDotOnes_ = trueQ_.colwise().sum();
  }

  template <typename mat_t>
  void checkR(const mat_t & R){
    EXPECT_EQ( R.cols(), 4 );
    EXPECT_EQ( R.rows(), 4 );
    for (auto i=0; i<4; i++)
      for (auto j=0; j<4; j++)
	EXPECT_NEAR( std::abs(R(i,j)), std::abs(trueR_(i,j)), 1e-6);
  }

  template <typename vec_t>
  void checkYForRsolve(const vec_t & y){
    EXPECT_EQ(y.extent(0), 4);
    for (auto i=0; i<y.extent(0); i++)
      EXPECT_NEAR( std::abs(trueYForRSolve_[i]),
		   std::abs(y[i]), 1e-12 );
  }

private:

  void fillTrueQ(){
    trueQ_ <<
      -0.897235446547271, -0.039024431200404,  0.173692309541681,
      -0.034843851055998,  0.06844323641866,  -0.162666577780091,
      -0.355064595330458, -0.069618789726194,  0.,
      //-------------------------------------------------------
      -0.336463292455227, -0.074296513246923, -0.7624257912921,
      -0.286577009249518, -0.15470306042946,   0.249902247740785,
       0.353018563741718,  0.082744998721086,  0.,
      //-------------------------------------------------------
      -0.,                 0.530332013749077, -0.155702484341873,
       0.126884420222574, -0.42395658892426,  -0.401481903098228,
      -0.154500472643382,  0.560006903183078,  0.,
      //-------------------------------------------------------
      -0.,                 0.530332013749077, -0.155702484341873,
      -0.139879089512256,  0.326913228407908, -0.468149250960295,
       0.260689856536577, -0.530040895116606,  0.,
      //-------------------------------------------------------
      -0.,                 0.,                -0.32818521819755,
       0.421198101926496,  0.751758650006291,  0.022892491611334,
      -0.138183613491497,  0.360730891969511,  0.,
      //-------------------------------------------------------
      -0.056077215409204,  0.650532264978526,  0.126820727560628,
      -0.069095509456248,  0.063878165111811,  0.719648521641246,
      -0.169043632892447, -0.04748965432464,   0.,
      //-------------------------------------------------------
      -0.280386077046022,  0.083927542741894,  0.333731413505014,
       0.469211836369866, -0.046150317046723,  0.0767206472791,
       0.746393155145893,  0.132984059523448,  0.,
      //-------------------------------------------------------
      -0.,                 0.014729971681881,  0.323860581694954,
      -0.691846753372529,  0.336778411839413, -0.08020531631897,
       0.229574011030821,  0.493648258802424,  0.,
      //-------------------------------------------------------
      -0., 0., 0.,-0., 0., 0., 0., 0.,  1.;
      //-------------------------------------------------------
  }//fillQ

  void fillTrueR(){
    trueR_ <<	-3.566510900025401, -1.665493297653372,
		-0.563576014862505,  0.841158231138066,
		//---------------------------------------------
		0.,  7.542444914314701,
		0.8945995630306  , -1.224066825632553,
		//---------------------------------------------
		0.,  0.,
		3.047059844718702,  2.566115091128629,
		//---------------------------------------------
		0., 0, 0., -7.497277277495907;
  }//fillR

  void fillTrueYForRSolve(){
    trueYForRSolve_ << -0.3569028575307429, -0.249377911416946,
			0.1687065812379429, -0.0273363044398516;
  }
};

}}} //end namespace
#endif
