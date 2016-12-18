/*************************************************************************\
 SolutionB.h  - implements analytical solutions for system (B) 
\*************************************************************************/
#include <boost/math/special_functions/bessel.hpp>

// struct for analytical solution up to integration constants
struct solution{

    double
        tau, zeta,
        A0, A1, A2,
        C0, C1, C2, C3, K1, K2, K3, K4, K7, K8, 
        lambda1, lambda2, lambda3;
    solution() { };
    solution(
        double tau, double zeta,
        double A0, double A1, double A2,
	      double C0, double C1, double C2, double C3, 
	      double K1, double K2, double K3, double K4, double K7, double K8) :
        tau(tau), zeta(zeta),
        A0(A0), A1(A1), A2(A2),
        C0(C0), C1(C1), C2(C2), C3(C3), K1(K1), K2(K2), K3(K3), K4(K4), K7(K7), K8(K8) {
        lambda1 = sqrt(30.0/47.0);
	      lambda2 = sqrt(5.0/6.0);
	      lambda3 = sqrt(3.0/2.0);
    }
    double theta(double x, double y) {
      double r = sqrt(x*x + y*y), phi = atan2(y,x);
      r = r/tau;
      return C0 + K4*BI(0,lambda2*r) + K3*BK(0,lambda2*r) + C3*log(r) - (A0*tau*pow(r,2))/4. + 
      + A2*pow(r,2)*(2.0/3 - 1.0/16*pow(r,2))*pow(tau,3) +
      
      cos(phi)*(C2*r + K1*BI(1,lambda2*r) + K2*BK(1,lambda2*r) + C1*pow(r,-1) - 
      A1*tau*(5.0/27 + (1.0/3 - 1.0/54*pow(r,2)))*pow(r,2));
     }
    double s_r(double r, double phi) {
      r = r/tau;
      return (A0*r*tau)/2. - C3*pow(r,-1) + (A2*pow(r,3)*pow(tau,3))/4. + 
      cos(phi)*(K7*BI(1,lambda1*r)*pow(r,-1) + K8*BK(1,lambda1*r)*pow(r,-1) - 
      A1*r*tau*(-2.0/3 + 2.0/27*pow(r,2)) - C2 + C1*pow(r,-2) );
     }
    double s_phi(double r, double phi) {
      r = r/tau;
      return  -( (K7*lambda1*(BI(0,lambda1*r) + BI(2,lambda1*r)))/2. - 
      (K8*lambda1*(BK(0,lambda1*r) + BK(2,lambda1*r)))/2. - 
      A1*r*tau*(-1.0/3 + 1.0/54*pow(r,2)) - C2 - C1*pow(r,-2) )*sin(phi);
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
  if( SystemText.find("angular") != -1 ) {
    if( fabs(kappa - 2.0) < 1e-5 ) {  // angular sym, L2-BC
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2, 0.,-111.13916942083414,6.660569535327411,0.,-1.5467495693523e-7,-68.42914074629492,0.,0.,-1.6653724783169088e-8,75.25683558006557);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2, 0.,-0.12059077270413276,-0.227493425227386,0.,-0.22480575559819507,0.13226674757262266,0.,0.,0.0018817968308445336,0.11292256050210862);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2, 0.,-0.0028464275340363515,-32.99259427482727,0.,70.11677186241039,0.0013720775414666753,0.,0.,-81.1826112603669,0.001575527948114299);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2, 0.,-0.00003515261705927927,-5300.849059964678,0.,11607.28163008687,0.000013936541710291672,0.,0.,-13268.053747129135,0.00001753818745435866);
    };
    if( fabs(kappa - 0.0) < 1e-5 ) {  // angular sym, simple-BC
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2, 0.,-113.60118968229399,6.665977737329732,0.,-1.5436749612392012e-7,-7.145508969223169,0.,0.,-3.9934390840255894e-8,
   141.25586158213324);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2, 0.,-1.400871510085252,-0.46971858927299437,0.,0.2515347159703147,1.499822124011616,0.,0.,-0.5601225099345,
   1.2867933473653226);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2, 0.,0.073630583565882,586.2937877216784,0.,-1302.419430933237,-0.05656057881949366,0.,0.,1458.6984809917412,
   -0.051431506218773196);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2, 0.,0.0004466366390853809,41378.31361359346,0.,-90825.26266760116,-0.00032425996966908746,0.,0.,103577.24660976982,
   -0.00029760236956955037); 
    };
  };
  if( SystemText.find("radial") != -1 ) { 
    if( fabs(kappa - 0.0) < 1e-5 ) {  // radially sym, L2-BC (kappa = 0.0)
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-5.728296791046461,0.,0.,5.391875027366691,0.,0.,-18.85128610299785,9.110050427771062e-9,0.,0.);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-0.31387486744277676,0.,0.,0.3127881502526124,0.,0.,-0.17388051396380483,2.041573398301043,0.,0.);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-3141.1176770795305,0.,0.,0.01681998730891037,0.,0.,-0.00787996610778202,3142.10302135665,0.,0.);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-3.195005048246991e6,0.,0.,0.0014373430258577003,0.,0.,-0.0006126816730987822,3.1950059051543605e6,0.,0.);
    };
    if( fabs(kappa - (-1.0)) < 1e-5 ) {  // radially sym, old-BC (kappa = -1.0)
      if( fabs(tau - 0.1) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-5.54017004655196,0.,0.,5.312471905917583,0.,0.,-18.60437569020495,9.121250709283261e-9,0.,0.);
      if( fabs(tau - 1.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-0.4434225663711735,0.,0.,0.26403046038004896,0.,0.,-0.15253260569832836,2.0392253876596267,0.,0.);
      if(  fabs(tau - 10.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-3139.190396368959,0.,0.,0.011183361968163205,0.,0.,-0.005413431007865912,3139.961757627532,0.,0.);
      if(  fabs(tau - 100.0) < 1e-5 )
      solX = new solution(tau, zeta, A0, A1, A2,-3.1947583612112333e6,0.,0.,0.000871180699385234,0.,0.,-0.00036386886858405244,3.1947590075608473e6,0.,0.);
    };
  };
};

