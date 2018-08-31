#!/bin/bash
docker build --rm -t smfretlikelihood .
mkdir -p ../dist
cp smfretlikelihood ../dist
cp install*.sh ../dist
if [ "$#" -ne 0 ]; then
    if [ "$1" = "tar" ]; then
        rm ../dist/smfretlikelihood.tar
        docker save smfretlikelihood > ../dist/smfretlikelihood.tar
    fi    
fi

