/*************************************************************************\
 Tools.h  - declaration of some useful tools. Definitions are in 
            Tools.cpp.
\*************************************************************************/

// string printf
string strprintf( const char* format, ...);

// Time Stamp Strings
string TimeStamp();
string CompleteTimeStamp();

// simple time measuring tools
void startTimer( string s );
void stopTimer();
double GlobalTime();

// simple command line access
int exists_arg(int argc,const char **argv,const char *keystr);
float get_float_arg(int argc,const char **argv,const char *keystr);
string get_string_arg(int argc,const char **argv,const char *keystr);

// sort/ordering
struct Ordering
{
  Ordering( VectorXd &Array );
  bool operator() (int i,int j);    
  VectorXi index;  
  VectorXd array;
};


