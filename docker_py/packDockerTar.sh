#!/bin/bash
rm -rf algo ui
cp -a ../algo .
cp -a ../ui .
docker build --rm -t smfretlikelihood .
# docker save smfretlikelihood > ../data/smfretlikelihood.tar