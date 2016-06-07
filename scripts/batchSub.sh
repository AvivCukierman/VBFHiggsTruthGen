#!/bin/bash

[ "$USER" == "acukierm" ] && WorkDir=/u/at/acukierm/nfs/VBF/VBFHiggsTruthGen

SubFileLoc=`pwd`/_submitSingleJob.sh
DateSuffix=`date +%Y%m%d_%Hh%Mmin`

export BSUB_QUIET=

echo '#!/bin/bash
echo CD to $1
echo CMD is $4

cd $1
cmd=$4

echo MAKING TEMP DIR $2
JOBFILEDIR=$2
mkdir $JOBFILEDIR
REALOUT=$3
echo MADE TEMP DIR $JOBFILEDIR
echo WILL COPY TO $REALOUT

shift 4
echo Calling $cmd $*
$cmd $*
cp -r $JOBFILEDIR/*.root $REALOUT
echo COPYING to $REALOUT
rm -rf $JOBFILEDIR
' > $SubFileLoc
chmod u+x $SubFileLoc

#----------------
# OPTIONS

Process=2
ms="30"
mineta=0.
maxeta=5.0
Queue=short
nevents=1000
njobs=100

OutDirFinal=`pwd`/files/${DateSuffix}
mkdir -p $OutDirFinal
echo
echo "Submitting $njobs jobs each with $nevents events to $Queue"

echo "Job Parameter Grid:"
echo "M  - "$ms

for m in $ms ; do
    LogPrefix=`pwd`/logs/${DateSuffix}/${DateSuffix}_${m}_${nevents}
    mkdir -p `dirname $LogPrefix`
    echo $LogPrefix
    
    for (( ii=1; ii<=$njobs; ii++ )) ;  do
	OutDir=/scratch/${DateSuffix}_${ii}/
			
	bsub -q ${Queue} -R rhel60 -W60 -o $LogPrefix${ii}.log \
	    $SubFileLoc ${WorkDir} ${OutDir} ${OutDirFinal} \
	    setup_batch.sh  \
            --OutDir ${OutDir} \
	    --OutFile Sample_mu-${mu}_${nevents}_job-${ii} \
	    --Proc ${Process} \
	    --nevents ${nevents} \
            --ScalarMass ${m} \
	    --MinEta ${mineta} --MaxEta ${maxeta}
    done
done
