import numpy as np
import matplotlib.pyplot as plt
import readligo as rl

strain, time, chan_dict = rl.loaddata('H-H1_LOSC_4_V1-815411200-4096.hdf5',H1)

slice_list = rl.dq_channel_to_seglist(chan_dict['Data'])
for slice in slice_list:
	time_seg = time[slice]
	strain_seg = strain[slice]
	#action with strain segments here

start = 842656000
stop = 842670000
segList = rl.getsegs(start, stop, 'H1', flag='CBLOW_CAT2')

N = 10000
for (begin, end) in segList:
    #Uses the getstrain() method to load the data
    strain, meta, dq = rl.getstrain(begin, end, 'H1')

    #Creates a plot
    plt.figure()
    ts = meta['dt']
    rel_time = np.arange(0, end-begin, meta['dt'])
    plt.plot(rel_time[0:N], strain[0:N])
    plt.xlabel('Seconds since GPS ' + str(begin) )
plt.show()