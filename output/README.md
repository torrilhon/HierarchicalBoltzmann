# Output File Folder

To accommodate log-files after program runs and final simulation data output.

Three pre-computed logfiles with empirical error data exist

    LogTwoCirclesLargeA.txt      // results for system (A)
    LogTwoCirclesLargeB.txt      // results for system (B)
    LogTwoCirclesLargeFlow.txt   // results for the NSF and R13 system
 
These files can be explored with the Mathematica srcipt `ErrorEvaluation.nb` in
the evaluation folder.

An example simulation is docmented in

    LogComplexChannel.txt        // Boltzmann simulation for curved channel

and the result can be visualized with the Mathematica script `Visualization.nb` 
in the evaluation folder. The result data files follow a time-stamped 
naming convention `result[dd-mm-yyyy_hh-mm-ss].dat` and are documented
in the logfiles.
