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
  echo "  -h    Help                   (Display this message)"
  echo "  -s    sample_set             (e.g. -s WVZ2016_v0.1.6)"
  echo "  -t    tag                    (e.g. -t y2016_test)"
  echo "  -b    do bdt"
  echo "  -w    wwz only"
  echo "  -s    split vvv v vh"
  echo
  exit
}

# Command-line opts
while getopts ":s:t:bSwh" OPTION; do
  case $OPTION in
    s) SAMPLESET=${OPTARG};;
    t) TAG=${OPTARG};;
    b) DOBDT=true;;
    w) WWZONLY="_wwzonly";;
    S) SPLITVH="_split";;
    h) usage;;
    :) usage;;
  esac
done

if [ -z ${SAMPLESET} ]; then usage; fi
if [ -z ${TAG}  ]; then usage; fi

# to shift away the parsed options
shift $(($OPTIND - 1))

# Verbose
date
echo "================================================"
echo "$(basename $0) $*"
echo "$(basename $0) $*" >> $DIR/.$(basename $0).history
echo "------------------------------------------------"
echo "SAMPLESET      : ${SAMPLESET}"
echo "TAG            : ${TAG}"
echo "================================================"

if [ $DOBDT ]; then
combineCards.py \
    emu1=stats/${SAMPLESET}/${TAG}/emu_bdt_datacard${WWZONLY}_bin1${SPLITVH}.txt \
    emu2=stats/${SAMPLESET}/${TAG}/emu_bdt_datacard${WWZONLY}_bin2${SPLITVH}.txt \
    emu3=stats/${SAMPLESET}/${TAG}/emu_bdt_datacard${WWZONLY}_bin3${SPLITVH}.txt \
    emu4=stats/${SAMPLESET}/${TAG}/emu_bdt_datacard${WWZONLY}_bin4${SPLITVH}.txt \
    emu5=stats/${SAMPLESET}/${TAG}/emu_bdt_datacard${WWZONLY}_bin5${SPLITVH}.txt \
    offzA=stats/${SAMPLESET}/${TAG}/offz_bdt_datacard${WWZONLY}_bin1${SPLITVH}.txt \
    offzB=stats/${SAMPLESET}/${TAG}/offz_bdt_datacard${WWZONLY}_bin2${SPLITVH}.txt  > stats/${SAMPLESET}/${TAG}/stat.txt
else
combineCards.py \
    emu1=stats/${SAMPLESET}/${TAG}/emu_datacard${WWZONLY}_bin1${SPLITVH}.txt \
    emu2=stats/${SAMPLESET}/${TAG}/emu_datacard${WWZONLY}_bin2${SPLITVH}.txt \
    emu3=stats/${SAMPLESET}/${TAG}/emu_datacard${WWZONLY}_bin3${SPLITVH}.txt \
    emu4=stats/${SAMPLESET}/${TAG}/emu_datacard${WWZONLY}_bin4${SPLITVH}.txt \
    offzA=stats/${SAMPLESET}/${TAG}/offz_datacard${WWZONLY}_bin1${SPLITVH}.txt \
    offzB=stats/${SAMPLESET}/${TAG}/offz_datacard${WWZONLY}_bin2${SPLITVH}.txt \
    offzC=stats/${SAMPLESET}/${TAG}/offz_datacard${WWZONLY}_bin3${SPLITVH}.txt  > stats/${SAMPLESET}/${TAG}/stat.txt
fi
