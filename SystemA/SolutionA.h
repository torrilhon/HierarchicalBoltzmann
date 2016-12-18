/*************************************************************************\
 SolutionA.h  - implements analytical solutions for system (A) 
\*************************************************************************/
#include <boost/math/special_functions/bessel.hpp>

// struct for analytical solution up to integration constants
struct solution{

    double
        tau, zeta,
        A0, A1, A2,
        C1, C2, C3, C4, K1, K2,
        lambda1;
    solution() { };
    solution(
        double tau, double zeta,
        double A0, double A1, double A2,
        double C1, double C2, double C3, double C4, double K1, double K2) :
        tau(tau), zeta(zeta),
        A0(A0), A1(A1), A2(A2),
        C1(C1), C2(C2), C3(C3), C4(C4), K1(K1), K2(K2) {
        lambda1 = sqrt(2*zeta);
    }
    double theta(double x, double y) {
      double r = sqrt(x*x + y*y), phi = atan2(y,x);
      r = r/tau;
      return C3 + C4*zeta*log(r) - (A0*tau*zeta*pow(r,2))/4. + 
      cos(phi)*(C2*r*zeta + C1*zeta*pow(r,-1) - (A1*tau*(-2 + zeta*pow(r,2)))/3.) + 
      (2*A2*pow(r,2)*pow(tau,3))/3. - (A2*zeta*pow(r,4)*pow(tau,3))/16.;
    }
    double s_r(double r, double phi) {
      r = r/tau;
      return (A0*r*tau)/2. - C4*pow(r,-1) + cos(phi)*(-C2 + (2*A1*r*tau)/3. + C1*pow(r,-2) - 
	    K2*BI(1,r*lambda1)*pow(2,0.5)*pow(r,-1) + K1*BK(1,r*lambda1)*pow(2,0.5)*pow(r,-1))
	    + (A2*pow(r,3)*pow(tau,3))/4.;
    }
    double s_phi(double r, double phi) {
      r = r/tau;
      return (C2 - (A1*r*tau)/3. + C1*pow(r,-2) + 
      K2*(BI(0,lambda1*r)*pow(zeta,0.5) + BI(2,lambda1*r)*pow(zeta,0.5)) + 
      K1*(BK(0,lambda1*r)*pow(zeta,0.5) + BK(2,lambda1*r)*pow(zeta,0.5)))*sin(phi);
    }
    double s_x(double x, double y) {
      double r = sqrt(x*x + y*y), phi = atan2(y,x);
      return cos(phi) * s_r(r,phi) - sin(phi) * s_phi(r,phi);
    }    
    double BI( int n, double x ) { return boost::math::cyl_bessel_i(n,x); }
    double BK( int n, double x ) { return boost::math::cyl_bessel_k(n,x); }
};

// defines the integration constants for various cases 
void System::setupExact()
{
  if( SystemText.find("mixed") != -1 ) {
    if( fabs(kappa - 1.0) < 1e-5 ) {  // angular/radial, L2-BC
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-4.01158541255238,0.00912181074518269,-5.460256952965152,5.288321376634141,3056.3171171825607,-2.049209850510533e-16);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-0.8660700736646773,0.1469559929405082,2.612851381692125,0.23377258203784068,1.108165701182322,-0.012889353089531019);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-0.09513392903290223,3.370128868300453,14.441104663250647,0.023497664861756657,0.09540038379989395,-3.019422627148662);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-0.009938852672236868,42.0083289706534,139.93580191710387,0.002482757542509108,0.009939126238290338,-41.544679750449916);
    };
    if( fabs(kappa - 0.0) < 1e-5 ) {  // angular/radial, simple-BC
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-1.5344890193223897,0.0034708867647638745,-6.000111913805922,5.505007359332947,2978.3396962069764,-3.913905997606115e-17);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-0.2328090854784152,0.026142383393355065,3.138888983089352,0.8614005885635808,0.48106821824548013,-0.0014311991485566613);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-0.003086096578496753,-0.13269978080269537,24.046961217231864,0.9496081671066655,0.0033882467285893123,0.174022697668903);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.00023937059975533113,-2.8818956238064173,221.34696422106651,0.9943960536822419,-0.0002390612975664888,2.929045071741498);
    };
  };
  if( SystemText.find("radial0") != -1 ) {
    if( fabs(kappa - 1.0) < 1e-5 ) {  // radial sym, no source, L2-BC
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,1.4090818873107707,-0.29872247846830097,0.,0.);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,0.5617542446609948,-0.06547678446574094,0.,0.);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,0.473081495038552,-0.0012456827924228005,0.,0.);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,0.4674990475333622,-0.000013244546732703282,0.,0.);
    };
    if( fabs(kappa - 0.0) < 1e-5 ) {  // radial sym, no source, simple-BC
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,1.430679817129441,-0.30556849175829015,0.,0.);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,0.6535070519228661,-0.12865726410284523,0.,0.);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,0.5642485458769436,-0.018949231489539822,0.,0.);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,0.5916676406557906,-0.001988970803960152,0.,0.);
    };
  };
  if( SystemText.find("radial") != -1 ) {
    if( fabs(kappa - 1.0) < 1e-5 ) {  // radial sym, with source, L2-BC
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,-5.460256952965152,5.288321376634141,0.,0.);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,2.612851381692125,0.23377258203784068,0.,0.);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,14.441104663250647,0.023497664861756657,0.,0.);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,139.93580191710387,0.002482757542509108,0.,0.);
    };
    if( fabs(kappa - 0.0) < 1e-5 ) {  // radial sym, with source, simple-BC
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,-6.000111913805922,5.505007359332947,0.,0.);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,3.138888983089352,0.8614005885635808,0.,0.);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,24.046961217231864,0.9496081671066655,0.,0.);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,0.,0.,221.34696422106651,0.9943960536822419,0.,0.);
    };
  };  
};

