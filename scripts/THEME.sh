#!/bin/bash

THEME_EXE=/nfs/data/cwng/archive/cvEM.64/THEME_edit.py

OPT_SPEC='
{
"NAME": "THEME.sh",
"DESC": "Run old THEME version",
"ARGS": ["FG_FASTA","BG_FASTA","HYP_FN","MARKOV"],
"OPTS": {
    "CV":{"LONG":"--cv","DEFAULT":5,"TYPE":"int","HELP":"number of cross validation folds [default:%default]"},
    "NOREFINE":{"LONG":"--no-refine","ACTION":"store_true","HELP":"do not run with refinement"},
    "BETA":{"LONG":"--beta","DEFAULT":0.7,"TYPE":"float","HELP":"beta parameter to use [default:%default]"},
    "DELTA":{"LONG":"--delta","DEFAULT":0.001,"TYPE":"float","HELP":"delta parameter to use [default:%default]"},
    "RANDOMIZE":{"LONG":"--randomization","ACTION":"store_true","HELP":"run randomization"},
    "MOTIF_FN":{"LONG":"--motif-file","DEFAULT":"dummy.out","HELP":"filename to write motif results to [default:%default]"},
    "OUTPUT_FN":{"LONG":"--output-filename","DEFAULT":"dummy.txt","HELP":"filename to write motif results to [default:%default]"},
    "RANDOM_FN":{"LONG":"--random-output","DEFAULT":"random.txt","HELP":"filename to write motif results to [default:%default]"},
    "DUMP":{"LONG":"--dump","ACTION":"store_true","HELP":"dump categtories to file"},
    "REM_COM":{"LONG":"--remove-common","ACTION":"store_true","HELP":"remove common sequences from analysis"},
    "NOPARALLEL":{"LONG":"--no-parallelize","ACTION":"store_true","HELP":"do not use wqsub.py for parallelization"},
    "INTERACTIVE":{"LONG":"--interactive","ACTION":"store_true","HELP":"run the script interactively"},
    "HYP_INDS":{"LONG":"--hyp-indices","DEFAULT":"ALL","HELP":"0-based indices of hypotheses to run [default: %default]"},
    "VERBOSE":{"SHORT":"-v","LONG":"--verbose","ACTION":"store_true","HELP":"print out the commands that are being run"},
    "TRIALS":{"LONG":"--trials","HELP":"this option is here only for backwards compatibility with THEME.py"}
    }
}'
OUTPUT=$(echo $OPT_SPEC | getopts.py -- $@)
GETOPTS_RET=$?
if [ $GETOPTS_RET -ne 0 ]; then
    exit 1
fi
$OUTPUT

INTERACTIVE_FLAG="--auto"
if [ $INTERACTIVE != "None" ]; then
    INTERACTIVE_FLAG=
fi

eval "$(steplist.py $INTERACTIVE_FLAG -t "Run THEME" THEME "Wait for jobs" "Combine results")"

# run THEME
OUTDIR=THEME_data
test \! -e $OUTDIR && mkdir $OUTDIR

WQSUB_EXE="wqsub.py"
if [ $NOPARALLEL != "None" ]; then
    WQSUB_EXE=
fi

RANDOMIZE_FLAG=
if [ $RANDOMIZE != "None" ]; then
    RANDOMIZE_FLAG="-randomization"
fi

RC=
if [ $RC ]; then
    RC='-rc'
fi

if [ $HYP_INDS != "ALL" ]; then
    HYP_INDS=$(parse_steplist.py $HYP_INDS)
    HYP_INDS_STATUS=$?
    if [ $HYP_INDS_STATUS != 0 ]; then
        echo "Incorrectly formatted argument to --hyp-indices option, aborting"
        exit $HYP_INDS_STATUS
    fi
else
    NUM_HYPS=`grep -c '^Source' $HYP_FN`
    NUM_HYPS=$(($NUM_HYPS-1))
    HYP_INDS=$(seq 0 $NUM_HYPS)
fi

JOBIDS=
next_step && \
for i in $HYP_INDS
do

    WQSUB=
    REDIRECT=
    if [ ! -z $WQSUB_EXE ]; then
        WQSUB="$WQSUB_EXE --wqsub-name=THEME_$i"
    fi

    OUTPRE=$OUTDIR/$i

    CMD="$WQSUB python $THEME_EXE $FG_FASTA $BG_FASTA $i \
        -fse $HYP_FN -markov $MARKOV -cv $CV -beta $BETA \
        -delta $DELTA -motif_file $OUTPRE.tamo -out_file $OUTPRE.txt \
        $RC"
    JOBID=$($WQSUB $CMD)
    JOBIDS="$JOBID $JOBIDS"
    if [ $VERBOSE != "None" ]; then
        echo $WQSUB $CMD
    fi

    if [ $RANDOMIZE != "None" ]; then

        WQSUB="$WQSUB_EXE --wqsub-name=THEME_rand_$i"

        CMD="$WQSUB python $THEME_EXE $FG_FASTA $BG_FASTA $i \
            -fse $HYP_FN -markov $MARKOV -cv $CV -beta $BETA \
            -delta $DELTA -out_file ${OUTPRE}_rand_output.txt \
            -random_file ${OUTPRE}_rand.txt $RC -randomization"

        JOBID=$($WQSUB $CMD)
        JOBIDS="$JOBID $JOBIDS"

        if [ $VERBOSE != "None" ]; then
            echo $WQSUB $CMD -randomization
        fi
    fi

done


# wait for jobs
next_step && wait_for_jobid.py $JOBIDS

# compile results
next_step
DO_COMPILE=$?
if [ $DO_COMPILE == 0 ]; then

    rm -f $MOTIF_FN && touch $MOTIF_FN
    (
        cd $OUTDIR
        ls *.tamo | sort -n | xargs -n1 -I{} -t cat {} >> ../$MOTIF_FN
    )

    if [ $NOPARALLEL == "None" ]; then
        mv -f *.{err,out} THEME_data
    fi

    if [ $RANDOMIZE != "None" ]; then
        rm -f $RANDOM_FN && touch $RANDOM_FN
        (
            cd $OUTDIR
            for ind in $HYP_INDS
            do
                out_fn="${ind}_rand.txt"
                echo "Consolidating $out_fn"
                python >> ../$RANDOM_FN << EOF
import re
import sys

from TAMO.MotifTools import load

ind = re.match('(\d+)',"$out_fn").group(1)

motif = load("$HYP_FN")[int(ind)]

src = motif.source.split()
if len(src) == 0 :
    print 'Got weird motif source: %s\n'%src
src = src[0]+'_%s'%ind

cverrs = []
for l in open("$out_fn") :
    m = re.match("trial: \d+ mean test error: (\d+\.\d+)$",l)
    if m is not None :
         cverrs.append(float(m.group(1)))

print "\t".join([src,str(sum(cverrs)/len(cverrs)),repr(cverrs)])
sys.stdout.flush()

EOF
            done

        )
    fi
fi
