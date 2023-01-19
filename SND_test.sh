#!/bin/bash
source /cvmfs/sndlhc.cern.ch/SNDLHC-2022/June10/setUp.sh
source  /afs/cern.ch/user/o/onur/SNDLHC_BETA/config.sh
python /afs/cern.ch/work/o/onur/scripts/US_scripts/main.py -r $2  -p /eos/experiment/sndlhc/convertedData/physics/2022/ -g geofile_sndlhc_TI18_V4_10August2022.root -t DS -npj 250000 -j $(expr $1 % 500)  -platform HTCondor
xrdcp *.root /eos/user/o/onur/sndlhc_recalibrated_data/run00$2
#sndlhc_MuFi_data/run00$2/
