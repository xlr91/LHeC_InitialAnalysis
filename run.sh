make > make.out; 
./Process /disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_H4l.Delphes.root signal.root s > signal.log 2>&1
./Process /disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_Zll.100k.v2.Delphes.root bakgnd.root b > bakgnd.log 2>&1
rm make.out;
