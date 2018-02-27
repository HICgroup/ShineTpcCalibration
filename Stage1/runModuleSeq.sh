#!/bin/bash

### Get settings from command line argument.
INPUTFILE=$1

INTERNALOUTPUTS="tpcPhases.root mhtdcDiffs.root"

### Check input arguments.
if [[ "$INPUTFILE" == "" ]]  ; then echo "No input file name given!" 1>&2  ; exit 1 ; fi
if ! [[ -f $INPUTFILE ]] ; then echo "Input file $INPUTFILE not existing as a local file!" 1>&2 ; exit 1 ; fi

### Get prefix of input file.
INPUTFILEBASE=`basename $INPUTFILE`
ND=`echo "$INPUTFILEBASE" | tr -d -c '.' | wc -c`
if (( $ND == 0 )) ; then ND=1 ; fi
INPUTPREFIX=`echo "$INPUTFILEBASE" | tr -d '\n' | cut -f -$ND -d '.'`

### Compile silently if not compiled.
make -s

### Export input file so that it can be read from EventFileReader.xml
export INPUTFILE

### Run the Shine module sequence.
ipcs -l
ShineOffline -b bootstrap.xml
EXITCODE=$?
ds kill > /dev/null 2>&1
if (( $EXITCODE != 0 )) ; then echo "Shine exited with exit code $EXITCODE !" 1>&2 ; exit $EXITCODE ; fi

### Save output under appropriate name.
for INTERNALOUTPUT in $INTERNALOUTPUTS ; do
  ND=`echo "$INTERNALOUTPUT" | tr -d -c '.' | wc -c`
  if (( $ND == 0 )) ; then ND=1 ; fi
  PREFIX=`echo "$INTERNALOUTPUT" | tr -d '\n' | cut -f -$ND -d '.'`
  SUFFIX=`echo "$INTERNALOUTPUT" | tr -d '\n' | cut -f $(( $ND + 1 ))- -d '.'`
  mv $INTERNALOUTPUT ${PREFIX}_${INPUTPREFIX}.${SUFFIX}
done
