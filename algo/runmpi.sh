#!/bin/bash
if [[ $# -ge 1 ]];then
  #echo $1
  mpiexec -f host -n 184 python maxLikMPI.py $1
else
  #echo "lala"
  mpiexec -f host -n 184 python maxLikMPI.py
fi
