/*************************************************************************\
 main.cpp  - main program routine and logfile output
\*************************************************************************/
#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "EigenSetup.h"
#include "Mesh.h"
#include "System.h"
#include "Numerics.h"
#include "Tools.h"


int main(int argc, const char **argv)
{
  string MeshText, SysText, NumText, PostText, TimeStamp;
  string MasterLogFile("../output/MasterLog.txt");
  
  startTimer(" *************************************\n program start "); 
  stopTimer();
     
  // construct mesh
  MeshCore Mesh( argc, argv ); 
  MeshText = Mesh.InfoText();
  
  // initialize system
  System Sys( argc, argv );    
  SysText = Sys.InfoText();
  
  // numerical setup + solve
  Numerics Num( &Mesh, &Sys, argc, argv );
  Num.Assemblation();
  Num.LinearSolve();  
  Num.AttachSolution();
  NumText = Num.InfoText();
  
  // time stamp and total runtime
  TimeStamp = CompleteTimeStamp();
   
  // error computation and/or output 
  PostText = Sys.PostProcess( &Mesh );

  // write logfile
  if(exists_arg(argc,argv,"-mfile")) MasterLogFile = get_string_arg(argc,argv,"-mfile");
  FILE *fp;
  fp = fopen( MasterLogFile.c_str(), "a+" );
  fprintf( fp, " %s\n%s\n%s\n%s\n%s\n************\n", 
           TimeStamp.c_str(), SysText.c_str(), MeshText.c_str(), NumText.c_str(), PostText.c_str() );
  fclose( fp );
  cout << endl;
  return(0);
}  

