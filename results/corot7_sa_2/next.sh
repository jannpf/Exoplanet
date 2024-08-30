#!/bin/bash

# sleep 10800

runname=corot7_sa_2
mkdir -p $runname

# make
# ./main -t 4

cp corot7.txt $runname/
cp hip14810.txt $runname/
cp hd191939.txt $runname/
cp levels.txt $runname/
cp log_prior_weights.txt $runname/
cp posterior_sample.txt $runname/
cp sample.txt $runname/
cp sample_info.txt $runname/
cp sampler_state.txt $runname/
cp weights.txt $runname/
cp OPTIONS $runname/

cp display.py $runname/
cp showresults.py $runname/

cp next.sh $runname/
