#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Help
usage()
{
  echo "ERROR - Usage:"
  echo
  echo "      sh $(basename $0) OPTIONSTRINGS ..."
  echo
  echo "Options with arguments:"
  echo "  -h    Help                     (Display this message)"
  echo "  -y    year                     (e.g. -y 2016)"
  echo "  -t    Ntuple type              (e.g. -t WVZMVA, WVZ, or Trilep etc.)"
  echo "  -v    Ntuple version           (e.g. -v 0.0.9, 0.1.0, or etc.)"
  echo "  -T    tag                      (e.g. -T test1)"
  echo "  -f    filter cuts (comma sep.) (e.g. -f FiveLeptonsMT5th,SixLeptonsSumPtCut,ChannelEMuBDT)"
  echo "  -s    do syst                  (e.g. -s)"
  echo
  exit
}

# Command-line opts
while getopts ":y:t:v:T:f:skh" OPTION; do
  case $OPTION in
    y) YEAR=${OPTARG};;
    t) NTUPLETYPE=${OPTARG};;
    v) NTUPLEVERSION=${OPTARG};;
    T) TAG=${OPTARG};;
    f) FILTERCUTS=${OPTARG};;
    s) SYST=SYST;;
    k) SKIM=SKIM;;
    h) usage;;
    :) usage;;
  esac
done

if [ -z ${YEAR} ]; then usage; fi
if [ -z ${NTUPLETYPE}  ]; then usage; fi
if [ -z ${NTUPLEVERSION}  ]; then usage; fi
if [ -z ${TAG}  ]; then usage; fi

# to shift away the parsed options
shift $(($OPTIND - 1))

# Verbose
date
echo "================================================"
echo "$(basename $0) $*"
echo "$(basename $0) $*" >> $DIR/.$(basename $0).history
echo "------------------------------------------------"
echo "YEAR           : ${YEAR}"
echo "NTUPLEVERSION  : ${NTUPLEVERSION}"
echo "NTUPLETYPE     : ${NTUPLETYPE}"
echo "TAG            : ${TAG}"
echo "FILTERCUTS     : ${FILTERCUTS}"
echo "SYST           : ${SYST}"
echo "SKIM           : ${SKIM}"
echo "================================================"

# Sanity check that the sample referred exists
if [ ! -d ${WVZ_DATA_PATH}/babies/${NTUPLETYPE}${YEAR}_${NTUPLEVERSION}/ ]; then echo "Asked to run over the ${NTUPLEVERSION} of ${NTUPLETYPE} but ${WVZ_DATA_PATH}/babies/${NTUPLETYPE}${YEAR}_${NTUPLEVERSION}/ does not exists!"; exit; fi

echo " Will only show progress for ZZ to 4L sample. since that is the bottle neck processing"

rm -f .jobs_${YEAR}_${NTUPLEVERSION}_${NTUPLETYPE}_${TAG}_${SYST}_${SKIM}_${FILTERCUTS}.txt
for i in $(ls -r ${WVZ_DATA_PATH}/babies/${NTUPLETYPE}${YEAR}_${NTUPLEVERSION}/); do
    # # if [[ $i == *"ttbar_"* ]] || [[ $i == *"dy_"* ]] || [[ $i == *"wz_3l"* ]] || [[ $i == *"wz_2l"* ]] || [[ $i == *"zg_"* ]]; then
    # # if [[ $i == *"ttz_"* ]]; then
    # if [[ $i == *"zzz_4l"* ]]; then
    #     :
    # else
    #     continue
    # fi
    if [[ $i == *"zz_4l_powheg"* ]]; then
        echo ./Analysis.exe ${i} ${NTUPLETYPE}${YEAR}_${NTUPLEVERSION} ${TAG} ${SYST} ${SKIM} ${FILTERCUTS}" | tee outputs/${NTUPLETYPE}${YEAR}_${NTUPLEVERSION}/${TAG}/${SYST}${SKIM}${i}.log"  >> .jobs_${YEAR}_${NTUPLEVERSION}_${NTUPLETYPE}_${TAG}_${SYST}_${SKIM}_${FILTERCUTS}.txt
    else
        echo ./Analysis.exe ${i} ${NTUPLETYPE}${YEAR}_${NTUPLEVERSION} ${TAG} ${SYST} ${SKIM} ${FILTERCUTS}" > outputs/${NTUPLETYPE}${YEAR}_${NTUPLEVERSION}/${TAG}/${SYST}${SKIM}${i}.log 2>&1"  >> .jobs_${YEAR}_${NTUPLEVERSION}_${NTUPLETYPE}_${TAG}_${SYST}_${SKIM}_${FILTERCUTS}.txt
    fi
done

mkdir -p outputs/${NTUPLETYPE}${YEAR}_${NTUPLEVERSION}/${TAG}
sh rooutil/xargs.sh -n 9 .jobs_${YEAR}_${NTUPLEVERSION}_${NTUPLETYPE}_${TAG}_${SYST}_${SKIM}_${FILTERCUTS}.txt
