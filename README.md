# smFRETLikelihood
The gSMFRETda branch of smFRETLikelihood is a set of Python toolbox to help you preparing hdf5 data file and serving as a  parameter server for running [gSMFRETda](https://github.com/liu-kan/gSMFRETda)

## Install
```bash
#if ubuntu
sudo apt install build-essential libhdf5-dev pkg-config protobuf-compiler libprotobuf-dev libnanomsg-dev
#elif centos or redhat
sudo dnf install protobuf-devel hdf5-devel nanomsg-devel python3-protobuf python3-devel 
#endif
pip3 install --user -r requirements.txt
```
##  Prepare data for gSMFRETda

```bash
git clone --single-branch --branch gSMFRETda https://github.com/liu-kan/smFRETLikelihood.git --depth=1
cd smFRETLikelihood
python3 untils/ptu2hdf.py
python3 untils/arrivalTimePDAdata.py
```

## Trouble shooting
If you meet errors like
```bash
Traceback (most recent call last):
  File "untils/arrivalTimePDAdata.py", line 10, in <module>
    from data import  prepData
  File "/home/liuk/data/build/g/smFRETLikelihood/untils/data/prepData.py", line 17, in <module>
    from fretbursts import *
  File "/home/liuk/miniconda3/envs/gSMFRETda/lib/python3.7/site-packages/fretbursts/__init__.py", line 144, in <module>
    from . import burst_plot as bpl
  File "/home/liuk/miniconda3/envs/gSMFRETda/lib/python3.7/site-packages/fretbursts/burst_plot.py", line 43, in <module>
    from matplotlib.mlab import normpdf
ImportError: cannot import name 'normpdf' from 'matplotlib.mlab' (/home/liuk/miniconda3/envs/gSMFRETda/lib/python3.7/site-packages/matplotlib/mlab.py)
```
Build FRETBursts from my repository.
```bash
pip3 uninstall fretbursts # or conda uninstall fretbursts
git clone --depth=1 https://github.com/liu-kan/FRETBursts.git
cd FRETBursts
python3 setup.py build
python3 setup.py install
```
