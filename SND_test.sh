#!/bin/bash
source /cvmfs/sndlhc.cern.ch/SNDLHC-2022/June10/setUp.sh
source  /afs/cern.ch/user/o/onur/SNDLHC_BETA/config.sh
DIR=/eos/user/o/onur/sndlhc_veto_data/run00$2/
echo "$DIR"
FILE="$DIR"histograms_"$1"_TI18_run00"$2"_DS.root
echo "$FILE"
if [ ! -f "$FILE" ]; then
python /afs/cern.ch/work/o/onur/scripts/US_scripts/main.py -r $2  -p /eos/experiment/sndlhc/convertedData/physics/2022/  -t DS -npj 250000 -j $(expr $1 % 500)  -platform HTCondor
xrdcp *.root "$DIR"
fi
