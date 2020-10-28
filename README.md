# Install

    # fedroa
    pip3 install --user numpy scipy pandas sklearn  numba fretbursts h5py
    # Use mpi4py from your disto.

#  Prepare data for gSMFRETda

```bash
git clone https://github.com/liu-kan/smFRETLikelihood.git --depth=1
cd smFRETLikelihood
python3 untils/ptu2hdf.py
python3 untils/arrivalTimePDAdata.py
```