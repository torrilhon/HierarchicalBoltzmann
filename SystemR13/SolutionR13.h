/*************************************************************************\
 SolutionR13.h  - implements analytical solutions for R13 system 
\*************************************************************************/
#include <boost/math/special_functions/bessel.hpp>

// struct for analytical solution up to integration constants
struct solution{

    double
      tau, zeta,
      A0, A1, A2,
      C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10,
      lambda1, lambda2, lambda3; 
    solution() { };
    solution(
      double tau, 
      double A0, double A1, double A2,
      double C0, double C1, double C2, double C3, double C4, double C5, double C6, double C7, double C8, double C9,
      double K1, double K2, double K3, double K4, double K5, double K6, double K7, double K8, double K9, double K10) :
      tau(tau), 
      A0(A0), A1(A1), A2(A2),
      C0(C0), C1(C1), C2(C2), C3(C3), C4(C4), C5(C5), C6(C6), C7(C7), C8(C8), C9(C9), 
      K1(K1), K2(K2), K3(K3), K4(K4), K5(K5), K6(K6), K7(K7), K8(K8), K9(K9), K10(K10) {
  lambda1 = sqrt(5.0/9.0);
	lambda2 = sqrt(5.0/6.0);
	lambda3 = sqrt(3.0/2.0);
    }
    double theta(double x, double y) {
      double r = sqrt(x*x + y*y), phi = atan2(y,x);
      r = r/tau;
      return C8 + (2*K4*BI(0,lambda2*r))/5. + (2*K3*BK(0,lambda2*r))/5. - (4*C6*log(r))/15. + 
      cos(phi)*((-4*C2*r)/15. + (8*C4*r)/5. + (2*K10*BI(1,lambda2*r))/5. + (2*K9*BK(1,lambda2*r))/5. - 
      (2*C1*pow(r,-1))/15. - (4*C5*pow(r,-1))/5.);
    }
    double s_r(double r, double phi) {
      r = r/tau;
      return C6*pow(r,-1) + cos(phi)*(C2 - (C1*pow(r,-2))/2. - 
      (5*(2*K1*BI(1,lambda1*r)*pow(lambda1,-1)*pow(r,-1) + 2*K2*BK(1,lambda1*r)*pow(lambda1,-1)*pow(r,-1)))/2.);
    }
    double s_phi(double r, double phi) {
      r = r/tau;
      return (-C2 + (5*(K1*(BI(0,lambda1*r) + BI(2,lambda1*r)) - K2*(BK(0,lambda1*r) + BK(2,lambda1*r))))/2. - (C1*pow(r,-2))/2.)*sin(phi);
    }
    double u_r(double r, double phi) {
      r = r/tau;
      return C7*pow(r,-1) + cos(phi)*(C0 + C5*(-0.5 + log(r)) - (C3*pow(r,-2))/2. + 
      2*K1*BI(1,lambda1*r)*pow(lambda1,-1)*pow(r,-1) + 2*K2*BK(1,lambda1*r)*pow(lambda1,-1)*pow(r,-1) + (C4*pow(r,2))/2.);
    }
    double u_phi(double r, double phi) {
      r = r/tau;
      return -(C0 + K1*(BI(0,lambda1*r) + BI(2,lambda1*r)) - K2*(BK(0,lambda1*r) + BK(2,lambda1*r)) + C5*(0.5 + log(r)) + 
       (C3*pow(r,-2))/2. + (3*C4*pow(r,2))/2.)*sin(phi);
    }
    double sig_rr(double r, double phi) {
      r = r/tau;
      return (4*C6*pow(r,-2))/5. + 2*C7*pow(r,-2) + K6*BI(1,lambda3*r)*pow(lambda3,-1)*pow(r,-1) + 
     K5*BK(1,lambda3*r)*pow(lambda3,-1)*pow(r,-1) + K4*(-BI(2,lambda2*r) - (BI(1,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2.) + 
     K3*(-BK(2,lambda2*r) + (BK(1,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2.) + 
     cos(phi)*(-2*C4*r - (4*C1*pow(r,-3))/5. - 2*(C3 - (64*C5)/15.)*pow(r,-3) - 2*C5*pow(r,-1) + 
      K7*BI(2,lambda3*r)*pow(lambda3,-1)*pow(r,-1) + K8*BK(2,lambda3*r)*pow(lambda3,-1)*pow(r,-1) + 
      K10*(-BI(1,lambda2*r) + (3*BI(2,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2.) + 
      K9*(-BK(1,lambda2*r) - (3*BK(2,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2.));
    }
    double sig_rp(double r, double phi) {
      r = r/tau;
      return (2*C4*r - (4*C1*pow(r,-3))/5. - 2*(C3 - (64*C5)/15.)*pow(r,-3) + (3*K10*BI(2,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2. - 
     (3*K9*BK(2,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2. + K7*BI(2,lambda3*r)*pow(lambda3,-1)*pow(r,-1) + 
     K8*BK(2,lambda3*r)*pow(lambda3,-1)*pow(r,-1))*sin(phi);
    }
    double sig_pp(double r, double phi) {
      r = r/tau;
      return (-4*C6*pow(r,-2))/5. - 2*C7*pow(r,-2) - K4*(-BI(0,lambda2*r)/2. + (3*BI(1,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2.) - 
       K3*(-BK(0,lambda2*r)/2. - (3*BK(1,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2.) - 
       K6*(-BI(0,lambda3*r) + BI(1,lambda3*r)*pow(lambda3,-1)*pow(r,-1)) - 
       K5*(BK(0,lambda3*r) + BK(1,lambda3*r)*pow(lambda3,-1)*pow(r,-1)) - 
       cos(phi)*(-2*C4*r - (4*C1*pow(r,-3))/5. - 2*(C3 - (64*C5)/15.)*pow(r,-3) - 2*C5*pow(r,-1) + 
      K10*(-BI(3,lambda2*r)/2. - (BI(2,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2.) + 
      K9*(-BK(3,lambda2*r)/2. + (BK(2,lambda2*r)*pow(lambda2,-1)*pow(r,-1))/2.) + 
      K7*(-BI(1,lambda3*r) + BI(2,lambda3*r)*pow(lambda3,-1)*pow(r,-1)) + 
      K8*(BK(1,lambda3*r) + BK(2,lambda3*r)*pow(lambda3,-1)*pow(r,-1)));
    }
    double u_x(double x, double y) {
      double r = sqrt(x*x + y*y), phi = atan2(y,x);
      return cos(phi) * u_r(r,phi) - sin(phi) * u_phi(r,phi);
    }    
    double s_x(double x, double y) {
      double r = sqrt(x*x + y*y), phi = atan2(y,x);
      return cos(phi) * s_r(r,phi) - sin(phi) * s_phi(r,phi);
    }    
    double sig_xy(double x, double y) {
      double r = sqrt(x*x + y*y), phi = atan2(y,x);
      return 0.5*(2*cos(2*phi)*sig_rp(r,phi) - sig_pp(r,phi)*sin(2*phi) + sig_rr(r,phi)*sin(2*phi));
    }    
    
    double BI( int n, double x ) { return boost::math::cyl_bessel_i(n,x); }
    double BK( int n, double x ) { return boost::math::cyl_bessel_k(n,x); }    
};

// defines the integration constants 
void System::setupExact()
{ // tau = 1.0
  solX = new solution(tau, A0, A1, A2,1.223903717080178,-0.5446547711684675,-0.14232449420977364,0.721586884988257,-0.06514216869392382,0.11496182185955707,-0.18670207167189568,3.7168702001158503e-7,
   1.9632445277597919,-0.015364900162247102,-0.025702982175726672,0.03063181277382055,-0.0722115557892615,0.012633722193122355,0.03846904408740305,0.003998658010456889,
   0.07955172795371973,0.04765241431972761,-0.0016312206698076686,0.27142027373457595); 
};


