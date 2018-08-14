#!/bin/sh
lam=0.0;
for dlam in 0.1 0.2
do
 for ran in 100 145
 do
    DIR=LAM"$lam"_DLAM"$dlam"_SEED"$ran"
    echo $lam, $dlam, $ran
    echo $DIR
    sed "s/LAM/$lam/g" < origQ > Q1
    sed "s/dR/$dlam/g" < Q1 > Q2
    sed "s/rSD/$ran/g" < Q2 > QUEUE
    rm -rf Q1 Q2
    mkdir $DIR
    cd $DIR
    cp ../2du1link .
    cp ../QUEUE .
    ./2du1link > out
    cd ..  
 done
done
