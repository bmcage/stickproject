import os
import fipy
import matplotlib.pyplot as plt
import numpy as np

treshold = 2e-6
FIGFILEEXT = '.png'

times = []
xpos = None
for file in sorted(os.listdir('.')):
    if file[-3:] == '.gz':
        data = fipy.tools.dump.read(file)
        times.append(data['time'])
        if xpos is None:
            xpos = data['xpos']
            sol = [None] * len(xpos)
            for i in range(len(sol)):
                sol[i] = []
        conc = data['concentration']
        for ind, C in enumerate(conc):
            sol[ind].append(C)

plt.ion()
for ind, pos in enumerate(xpos):
    plt.rc("font", family="serif")
    plt.rc("font", size=10)
    width = 4.5  #width in inches
    height = 1.4 #height in inches
    plt.rc("figure.subplot", left=(50/72.27)/width)
    plt.rc("figure.subplot", right=(width-10/72.27)/width)
    plt.rc("figure.subplot", bottom=(14/72.27)/height)
    plt.rc("figure.subplot", top=(height-7/72.27)/height)
    plt.figure(ind)
    plt.gca().set_xlabel('Time [s]')
    plt.gca().set_ylabel('Concentration [$\mu$g/mm$^3$]')
    #plt.gca().yaxis.set_major_formatter(pylab.FormatStrFormatter('%e'))
    plt.title('Concentration at position %g mm' % pos)
    plt.plot(times, sol[ind])
    #plt.ylim(0, maxv*1.1)
    plt.plot(times, np.ones(len(times)) * treshold, 'b--')
    plt.show()
    plt.savefig('proc_AIconc_%03.1f_mm' % pos + FIGFILEEXT)
        
raw_input('Finished, press key to stop')
