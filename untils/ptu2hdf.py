# coding: utf-8
import numpy as np
import phconvert as phc
def ptu2hdf5(f,fsave):
    d=phc.pqreader.load_ptu(f)

    import json
    # print(json.dumps(d[3],indent=4))
    # print(int(d[3]["tags"]["TTResult_SyncRate"]['value']))
    # print(d[3]["timestamps_unit"]+1)

    # import matplotlib
    # matplotlib.pyplot.hist(d[2],101)

    for det, count in zip(*np.unique(d[1], return_counts=True)):
        print("%8d   %8d" % (det, count))

    nanotimes = d[2]
    detectors = d[1]
    timestamps = d[0]

    not_overflow = d[2] != 0

    detectors = detectors[not_overflow]
    timestamps = timestamps[not_overflow]
    nanotimes = nanotimes[not_overflow]
    print("Detector    Counts")
    print("--------   --------")
    for det, count in zip(*np.unique(detectors, return_counts=True)):
        print("%8d   %8d" % (det, count))


    measurement_specs = dict(
        measurement_type = 'smFRET-nsALEX',
        laser_repetition_rate=int(d[3]["tags"]["TTResult_SyncRate"]['value']),
        detectors_specs = {'spectral_ch1': [1],  # list of donor's detector IDs
                        'spectral_ch2': [0]},  # list of acceptor's detector IDs
        alex_excitation_period1=[10,1240],
        alex_excitation_period2=[1260,2490]
        )
    tcspc_num_bins=max(d[2])
    nanotimes_specs=dict(
        tcspc_unit=float(d[3]["tags"]["MeasDesc_Resolution"]['value']),
        tcspc_num_bins=tcspc_num_bins,
        tcspc_range=float(d[3]["tags"]["MeasDesc_Resolution"]['value'])*tcspc_num_bins
    )
    photon_data = dict(
        timestamps=timestamps,
        detectors=detectors,
        timestamps_specs={'timestamps_unit': d[3]["timestamps_unit"]},
        measurement_specs=measurement_specs,
        nanotimes_specs=nanotimes_specs,
        nanotimes=nanotimes
    )
    setup = dict(
        ## Mandatory fields
        num_pixels = 2,                   # using 2 detectors
        num_spots = 1,                    # a single confoca excitation
        num_spectral_ch = 2,              # donor and acceptor detection 
        num_polarization_ch = 1,          # no polarization selection 
        num_split_ch = 1,                 # no beam splitter
        excitation_cw=[False,False],
        modulated_excitation = True,     # CW excitation, no modulation 
        excitation_alternated = [True,True],  # CW excitation, no modulation 
        lifetime = True                 #  TCSPC in detection
    )

    description = 'This is a fake dataset which mimics smFRET data.'

    author = 'kk'
    author_affiliation = 'Name of Research Institution'
    identity = dict(
        author=author,
        author_affiliation=author_affiliation)
    data = dict(
        description=description,
        photon_data = photon_data,
        setup=setup,
        identity=identity,
    )
    phc.hdf5.save_photon_hdf5(data, h5_fname=fsave, overwrite=True)

# timestamps
# float(d[3]["tags"]["MeasDesc_Resolution"]['value'])

def usage():  
    print("Usage:%s -i|--ptu inputfilename.ptu -o|--hdf outputfilename.hf5" % sys.argv[0])

if __name__ == '__main__':
    import sys,getopt
    dbname=''
    savefn=''
    try:  
        opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["ptu=", "hdf="])  
        for o, v in opts: 
            if o in ("-i", "--ptu"):
                dbname=v
            if o in ("-o", "--hdf"):
                savefn = v
    except getopt.GetoptError:  
        # print help information and exit:  
        print("getopt error!")    
        usage()    
        sys.exit(1)
    if len(dbname)<1 or len(savefn)<1:
        usage()    
        sys.exit(1)
    ptu2hdf5(dbname,savefn)
