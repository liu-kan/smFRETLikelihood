#!/bin/bash
docker build --rm -t smfretlikelihood .
mkdir ../dist
docker save smfretlikelihood > ../dist/smfretlikelihood.tar
cp smfretlikelihood ../dist
cp install*.sh ../dist
