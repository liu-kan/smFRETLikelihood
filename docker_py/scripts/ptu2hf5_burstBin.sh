cd $HOME
cd data/smFRETLikelihood
python3 untils/ptu2hdf.py -i $HOME/labdata/$1 -o $HOME/labdata/tmp.hf5
python3 algo/burstBin.py - $HOME/labdata/tmp.hf5 -i $HOME/labdata/$2