/*************************************************************************\
 Tools.cpp  - some useful tools: time stamp, memory usage, time 
              measurements, command line access, sort/ordering
\*************************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "EigenSetup.h"
#include "Tools.h"

// string printf
string strprintf( const char* format, ...)
{
  char tmp[300];
  va_list args;
  va_start(args, format);
  vsprintf(tmp, format, args);
  return(string(tmp));
}

// get time stamp suitable as file name
string TimeStamp()
{
  char output[300];
  time_t rawtime;
  struct tm * timeinfo;
  
  time( &rawtime );
  timeinfo = localtime ( &rawtime );
  
  sprintf(output, "[%d-%d-%d_%d-", timeinfo->tm_mday, timeinfo->tm_mon + 1, timeinfo->tm_year + 1900, timeinfo->tm_hour );
	  
  if( timeinfo->tm_min < 10 ) 
    sprintf(output+strlen(output), "0%d", timeinfo->tm_min );
  else
    sprintf(output+strlen(output), "%d", timeinfo->tm_min );
	    
  if( timeinfo->tm_sec < 10 ) 
    sprintf(output+strlen(output), "-0%d]", timeinfo->tm_sec );
  else
    sprintf(output+strlen(output), "-%d]", timeinfo->tm_sec );

  return(string(output));
}

// Store time stamp with run time
string CompleteTimeStamp()
{
  return( strprintf( "%s, RunTime = %.2f", TimeStamp().c_str(), GlobalTime() ) );
}


/***********************************************************************
 * Two functions to query the memory usage (resident set size) of the 
 * current process: size_t getCurrentRSS()
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 */
#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

// Returns the current resident set size (physical memory use) measured in bytes
size_t getCurrentRSS()
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
	/* OSX ------------------------------------------------------ */
	struct mach_task_basic_info info;
	mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
	if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
		(task_info_t)&info, &infoCount ) != KERN_SUCCESS )
		return (size_t)0L;		/* Can't access? */
	return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
	/* Linux ---------------------------------------------------- */
	long rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
	/* AIX, BSD, Solaris, and Unknown OS ------------------------ */
	return (size_t)0L;			/* Unsupported. */
#endif
}
/***********************************************************************/


// simple time measuring tools

#include <sys/time.h>

static struct timeval time0;
static struct timeval time1;
static int globalTimer = 0;

static inline unsigned int GetmSec( timeval timeR )
{
  struct timeval time;
  gettimeofday(&time, NULL );  
  return( (unsigned int)(1.0e3*(time.tv_sec-timeR.tv_sec)+1.0e-3*(time.tv_usec-timeR.tv_usec)) );
}

void startTimer( string s )
{
  if( globalTimer == 0 ) {
    gettimeofday(&time1,NULL);
    globalTimer = 1;
  };
  cout << s;
  gettimeofday(&time0,NULL);
};

void stopTimer()
{
  cout << " ...done after " << 1.0*GetmSec(time0)/1000 << " seconds. (current size: " 
       << 1.0*getCurrentRSS()/(1024*1024) << " MB) \n";
};

double GlobalTime()
{
  if( globalTimer == 1 ) return( 1.0*GetmSec(time1)/1000 ); // in seconds
  else return( -1 );
};
/***********************************************************************/


// simple command line access

int exists_arg(int argc,const char **argv,const char *keystr)
{
  int i;

  for (i=1;i<argc;i++) if (strstr(argv[i],keystr)) return 1;
  return 0;
}

float get_float_arg(int argc,const char **argv,const char *keystr)
{
  int i;
  float result=0.0;

  for (i=1;i<argc;i++) {
    if (strstr(argv[i],keystr)) {
      result=atof(argv[i+1]);
      break;
    }
  }
  if (i == argc) {
    printf("Error: get_float_argument(%s)\n",keystr);
    return (float) -999999999;
  }
  else return result;
}
  
string get_string_arg(int argc,const char **argv,const char *keystr)
{
  int i;

  for (i=1;i<argc;i++) 
    if (strstr(argv[i],keystr)) 
      return(strprintf("%s",argv[i+1]));
  
  if (i == argc) 
    fprintf(stderr,"Error: get_string_argument(%s)\n",keystr);
  
  return(string("Error"));
}
/***********************************************************************/

// returns the ordering of a vector in an array

Ordering::Ordering( VectorXd &Array )
{
  int n = Array.size();
  array = Array;
  
  vector<int> sorted(n);
  for( int i; i<n; i++ ) sorted[i] = i;
  sort(sorted.begin(),sorted.end(),*this);
  index = Map<VectorXi>(sorted.data(),n);
};

bool Ordering::operator() (int i,int j) 
{ 
  return( array[i]<array[j] );
};
/***********************************************************************/


