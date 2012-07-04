import os
import const
import lib.utils.utils as utils
from bednet.bednetmodel import Bednet

#------------------------------------------------------------------------------

OUTPUTDIR = const.DATA_DIR

def read_bednet_sol(self,times):
        times = Bednet.times
        print times
        #for nr, time in enumerate(times):
        file = open(utils.OUTPUTDIR + os.sep + 'bednet_sol_%08d.gz'%str(6),'r')
        sols = file.read()
        print 'Time', sols['time'], 'Concentration', sols['concentration']  

read_bednet_sol(times)
