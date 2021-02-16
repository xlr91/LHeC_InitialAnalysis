# Introduction

The code "Process.cxx" takes as input a file from the Delphes simulation framework (https://cp3.irmp.ucl.ac.be/projects/delphes) and allows you to implement a basic physics analysis event selection and fill some histograms.

# How to compile the code

Set up the necessary environment:

cc7

source /disk/moose/general/asc/Build_cc7_Delphes342_Py8243/SetupEnv.sh


Compile the code:

make

# How to run the code

The program takes two arguments, the first is the input simulation file and the second is the name of a new output ROOT file which will contain the histograms produced by the code. To run do:

./Process [path/to/inputfile] [output.root]

(where you will need to substitute these dummy arguments with the appropriate path to the input file and your chosen name for the output file)

the path to the first input file is 
/disk/moose/general/asc/Y4_LHeC_2021/LHeC_CC_H4l.Delphes.root

the second one is
/disk/moose/general/asc/Y4_LHeC_2021/LHeC_CC_H4l.Had.Delphes.root

V2 is
/disk/moose/general/asc/Y4_LHeC_2021/LHeC_CC_H4l.Had.Delphes.v2.root

Theres another one, this one is the lighter one
/disk/moose/general/asc/Y4_LHeC_2021/LHeC_CC_H4l.Parton.LHE.Delphes.root

The sherpa backup
/disk/moose/general/asc/Y4_LHeC_2021/LHeC_CC_H4l.Sherpa.root

The Fixed file is
/disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_H4l.Delphes.root

The background sample is
/disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_Zll.100k.Delphes.root

The updated background sample is
/disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_Zll.100k.v2.Delphes.root
# Making changes to the code

Whenever you update the C++ code, make sure you re-compile the code before re-running it by issuing the following command again:

make

# Whenever you open a new terminal

Repeat the steps to setup the environment whenever you open a new terminal window, before you attempt to compile or run the code:

cc7

source /disk/moose/general/asc/Build_cc7_Delphes342_Py8243/SetupEnv.sh

