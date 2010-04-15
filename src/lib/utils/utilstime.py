#!/usr/bin env python

# Copyright (C) 2009  B. Malengier
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#-------------------------------------------------------------------------
#
# local modules
#
#-------------------------------------------------------------------------

import lib.utils.utilsbm as UtilsBm

#-------------------------------------------------------------------------
#
# Time settings
#
#-------------------------------------------------------------------------
class TimeParameters(object):
    def __init__(self, config, t_start = 0., t_end = 0., t_step = 0., 
                 exptime=None):
        if config != None:
            if config.get('init.initfromfile'):
                self.t_start = exptime[0]
            else:
                self.t_start = config.get('time.t0')
            ttype = config.get('time.type')
            if ttype == 'fixstep':
                self.t_step = config.get('time.dt')
                self.t_end = self.t_start \
                    + config.get('time.nrsteps')*self.t_step
            elif ttype == 'toend':
                self.t_end = config.get('time.tend')
                self.t_step = ((self.t_end - self.t_start) / 
                               float( config.get('time.nrsteps')))
            else:
                raise NotImplementedError, \
                    'time.type must be fixstep or toend' \
                    ', recieved ' + ttype

            print 'start %f, end %f' \
                % (self.t_start, self.t_end)
            #compute where output is wanted:
            self.t_out = [self.t_start]
            out_every = config.get('time.tout_every')
            if not out_every:
                self.t_out.append(self.t_end)
            else:
                print ' output calculated every %f' % out_every
                tout = self.t_start + out_every
                while tout < self.t_end:
                    self.t_out.append(tout)
                    tout += out_every
                if (not (self.t_out[-1] == self.t_end) 
                    and self.t_out[-1] < self.t_end):
                    self.t_out.append(self.t_end)
        else:
            self.t_start = t_start
            self.t_end   = t_end
            self.t_step  = t_step
            self.t_out   = [t_start, t_step]

        self.t_exp = []
        if exptime:
            self.t_exp = exptime
        self.t_out_all = sorted(UtilsBm.merge_nodup(self.t_out, self.t_exp))