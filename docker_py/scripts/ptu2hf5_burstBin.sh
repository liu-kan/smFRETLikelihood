#!/bin/bash
cd $HOME
cd data/smFRETLikelihood
python3 untils/ptu2hdf.py -i $HOME/labdata/$1 -o $HOME/labdata/tmp.hf5
python3 algo/burstBin.py -i $HOME/labdata/tmp.hf5 -o $HOME/labdata/$2
rm $HOME/labdata/tmp.hf5