#!/bin/bash
if [[ $# -ge 1 ]];then

  #python ../test/echo.py $2 $1
  mpiexec $1 python maxLikMPI.py $2
else
  echo "lala"
  mpiexec -f host -n 184 python maxLikMPI.py
fi
