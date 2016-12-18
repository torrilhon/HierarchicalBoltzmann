/*************************************************************************\
 SolutionNSF.h  - implements analytical solutions for NSF system 
\*************************************************************************/

// struct for analytical solution up to integration constants
struct solution{

    double
      tau, zeta,
      A0, A1, A2,
      C1, C2, C3, C4, C5, C6, C7, C8, D1, D2, D3, D4; 
    solution() { };
    solution(
      double tau, 
      double A0, double A1, double A2,
      double C1, double C2, double C3, double C4, double C5, double C6, double C7, double C8,
      double D1, double D2, double D3, double D4) :
      tau(tau), 
      A0(A0), A1(A1), A2(A2),
      C1(C1), C2(C2), C3(C3), C4(C4), C5(C5), C6(C6), C7(C7), C8(C8), D1(D1), D2(D2), D3(D3), D4(D4) {
    }
    double theta(double x, double y) {
      double r = sqrt(x*x + y*y), phi = atan2(y,x);
      r = r/tau;
      return 4.0/15.0*(D3 + D4*log(r) - (A0*tau*pow(r,2))/4. + cos(phi)*(D2*r + D1*pow(r,-1) - (A1*tau*pow(r,2))/3.) - (A2*pow(r,4)*pow(tau,3))/16.);
    }
    double s_r(double r, double phi) {
      r = r/tau;
      return (A0*r*tau)/2. + cos(phi)*(-D2 + (2*A1*r*tau)/3. + D1*pow(r,-2)) - D4*pow(r,-1) + (A2*pow(r,3)*pow(tau,3))/4.;
    }
    double s_phi(double r, double phi) {
      r = r/tau;
      return (D2 - (A1*r*tau)/3. + D1*pow(r,-2))*sin(phi);
    }
    double u_r(double r, double phi) {
      r = r/tau;
      return C1*pow(r,-1) + cos(phi)*(C7 + C8*log(r) - (C5*pow(r,-2))/2. + (C6*pow(r,2))/2.);
    }
    double u_phi(double r, double phi) {
      r = r/tau;
      return (C4*r + (C2*pow(r,-1))/2.) + (-C7 - C8 - C8*log(r) - (C5*pow(r,-2))/2. - (3*C6*pow(r,2))/2.)*sin(phi);
    }
    double sig_rr(double r, double phi) {
      r = r/tau;
      return 2*C1*pow(r,-2) + cos(phi)*(-2*C6*r - 2*C5*pow(r,-3) - 2*C8*pow(r,-1));
    }
    double sig_rp(double r, double phi) {
      r = r/tau;
      return C2*pow(r,-2) + (2*C6*r - 2*C5*pow(r,-3))*sin(phi);
    }
    double sig_pp(double r, double phi) {
      r = r/tau;
      return -2*C1*pow(r,-2) + cos(phi)*(2*C6*r + 2*C5*pow(r,-3) + 2*C8*pow(r,-1));
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
    
};

// defines the integration constants for various cases 
void System::setupExact()
{ // tau = 0.1
  solX = new solution(tau, A0, A1, A2,
  //0.,0.,0.,0.,-15.817784503771081,-0.0013209096987176169,-1.7573540614888068,0.9056047316086425,0.,0.,-3.1503195038973986,7.665265179237176); // A0=1.0
  //0.,0.,0.,0.,-15.817784503771081,-0.0013209096987176169,-1.7573540614888068,0.9056047316086425,0.,0.,1.2604086038272742,2.8240450660347496); // A0 = 0.2
  0.,0.,0.,0.,-15.817784503771081,-0.0013209096987176169,-1.7573540614888068,0.9056047316086425,0.,0.,2.3630906307584425,1.6137400377341427); // A0=0
};


