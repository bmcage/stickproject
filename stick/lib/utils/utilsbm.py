#!/usr/bin env python

# Copyright (C) 2000-2006  B. Malengier
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
# Some Utility functions
#
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# linear interpolation of a function in x, given 2 values of the
#    function, namely at x1 and x2.    1e order Lagrange interpolant  
#-------------------------------------------------------------------------
def inter2(x1,v1,x2,v2,x):
    return (v1) * ((x) - (x2)) / ((x1) - (x2))     \
            + (v2) * ((x) - (x1)) / ((x2) - (x1))

#-------------------------------------------------------------------------
# interpolation of a function in x, given 3 values of the
#    function, namely at x1, x2 and x3. 2e order Lagrange interpolant  
#-------------------------------------------------------------------------
def inter3(x1,v1,x2,v2,x3,v3,x):
    return (v1) * ((x) - (x2)) * ((x) - (x3)) / (((x1) - (x2)) * ((x1) - (x3))) \
        + (v2) * ((x) - (x1)) * ((x) - (x3)) / (((x2) - (x1)) * ((x2) - (x3))) \
        + (v3) * ((x) - (x1)) * ((x) - (x2)) / (((x3) - (x1)) * ((x3) - (x2)))

#-------------------------------------------------------------------------
# interpolation of a function in x, given 4 values of the
#    function, namely at x1, x2, x3 and x4. 3e order Lagrange interpolant  
#-------------------------------------------------------------------------
def inter4(x1,v1,x2,v2,x3,v3,x4,v4,x):
    return (v1) * ((x) - (x2)) * ((x) - (x3)) *((x) - (x4)) / (((x1) - (x2)) *      \
            ((x1) - (x3))* ((x1) - (x4))) \
            + (v2) * ((x) - (x1)) * ((x) - (x3)) * ((x) - (x4))/ (((x2) - (x1)) *   \
            ((x2) - (x3))*((x2) - (x4)))  \
            + (v3) * ((x) - (x1)) * ((x) - (x2)) * ((x) - (x4))/ (((x3) - (x1)) *   \
            ((x3) - (x2))* ((x3) - (x4))) \
            + (v4) * ((x) - (x1)) * ((x) - (x2)) * ((x) - (x3))/ (((x4) - (x1)) *   \
            ((x4) - (x2))* ((x4) - (x3)))

#-------------------------------------------------------------------------
# derivative at a border of 2 cells, with left value at x1 v1, and right
#    value at x2 v2
#-------------------------------------------------------------------------
def deriv12(x1,v1,x2,v2):
    return (((v2)-(v1))/((x2)-(x1)))

#-------------------------------------------------------------------------
# derivative based on 2e order Lagrange interpolant in point x2,
#    so deriv132=subs(x=x2, diff(inter3 , x))
#-------------------------------------------------------------------------
def deriv132(x1,v1,x2,v2,x3,v3):
    return (  (v1) * ((x2) - (x3)) / (((x1) - (x2)) * ((x1) - (x3)))              \
    + (v2) * (2 * (x2) - (x1) - (x3)) / (((x2) - (x1)) * ((x2) - (x3)))   \
    + (v3) * ((x2) - (x1)) / (((x3) - (x1)) * ((x3) - (x2))))

#-------------------------------------------------------------------------
# derivative based on 2e order Lagrange interpolant in the left point x1,
#   so deriv131=subs(x=x1, diff(inter3 , x))
#-------------------------------------------------------------------------
def deriv131(x1,v1,x2,v2,x3,v3):
    return deriv132(x2,v2,x1,v1,x3,v3)

#-------------------------------------------------------------------------
# derivative based on 2e order Lagrange interpolant in the right point x3,
#   so deriv131=subs(x=x3, diff(inter3 , x))
#-------------------------------------------------------------------------
def deriv133(x1,v1,x2,v2,x3,v3):
    return deriv132(x1,v1,x3,v3,x2,v2)

#-------------------------------------------------------------------------
# second derivative of the 2e order Lagrange interpolant
#   Note: this is a constant
#-------------------------------------------------------------------------
def deriv23(x1,v1,x2,v2,x3,v3) :
    return (2 * ((v1) / (((x1) - (x2)) * ((x1) - (x3)))                           \
        + (v2) / (((x2) - (x1)) * ((x2) - (x3)))                              \
        + (v3) / (((x3) - (x1)) * ((x3) - (x2)))))

#-------------------------------------------------------------------------
# third order derivative based on the 3e order Lagrange interpolant
#   through the points x1,x2,x3,x4 . This is a constant  
#-------------------------------------------------------------------------
def deriv34(x1,v1,x2,v2,x3,v3,x4,v4):
    return (6 * (  (v1) / (((x1) - (x2)) * ((x1) - (x3)) * ((x1) - (x4)))         \
      + (v2) / (((x2) - (x1)) * ((x2) - (x3)) * ((x2) - (x4)))         \
      + (v3) / (((x3) - (x1)) * ((x3) - (x2)) * ((x3) - (x4)))         \
      + (v4) / (((x4) - (x1)) * ((x4) - (x2)) * ((x4) - (x3)))))

#-------------------------------------------------------------------------
# second order derivative based on the 3e order Lagrange interpolant
#   in the point x1
#-------------------------------------------------------------------------
def deriv241(x1,v1,x2,v2,x3,v3,x4,v4) :
    return ( 2*(3*(x1)-(x4)-(x3)-(x2))/((x1)-(x2))/((x1)-(x3))/((x1)-(x4))*(v1)   \
    -2*(2*(x1)-(x4)-(x3))/((x1)-(x2))/((x2)-(x3))/((x2)-(x4))*(v2)        \
    +2*(2*(x1)-(x4)-(x2))/((x1)-(x3))/((x2)-(x3))/((x3)-(x4))*(v3)        \
    -2*(2*(x1)-(x3)-(x2))/((x1)-(x4))/((x2)-(x4))/((x3)-(x4))*(v4))

#-------------------------------------------------------------------------
# first order derivative based on the 3e order Lagrange interpolant
#   in the point x1
#-------------------------------------------------------------------------
def deriv141(x1,v1,x2,v2,x3,v3,x4,v4) :
    return ((3*(x1)*(x1)-2*(x1)*(x4)-2*(x3)*(x1)                                  \
      +(x3)*(x4)-2*(x2)*(x1)+(x2)*(x4)+(x2)*(x3))                          \
      /((x1)-(x2))/((x1)-(x3))/((x1)-(x4))*(v1)                            \
      -((x1)-(x3))*((x1)-(x4))/((x1)-(x2))/((x2)-(x3))/((x2)-(x4))*(v2)    \
      +((x1)-(x2))*((x1)-(x4))/((x1)-(x3))/((x2)-(x3))/((x3)-(x4))*(v3)    \
      -((x1)-(x2))*((x1)-(x3))/((x1)-(x4))/((x2)-(x4))/((x3)-(x4))*(v4))
    
def xmerge(*ln):
    '''
    from http://aspn.activestate.com/ASPN/Cookbook/Python/
    General version of merge, on sorted sequences. 
    '''
    from itertools import chain
    from heapq import heappop

    heap = []
    for i in chain(*ln):
        heap.append(i)

    while heap:
        yield heappop(heap)
    
def merge(*ln):
     """ Merge several sorted sequences into a sorted list.
     
     Assuming l1, l2, l3...ln sorted sequences, return a list that contain
     all the items of l1, l2, l3...ln in ascending order.
     Input values doesn't need to be lists: any iterable sequence can be used.
     """
     return list(xmerge(*ln))

def merge_nodup(ln1, ln2):
    '''Merge two lists, remove duplicates
    '''
    workingSet = set(ln1)
    workingSet =workingSet.union(ln2)
    return list(workingSet)
    
def piecewise3(x, first, func1, second, func2, func3):
    ''' define a piecewise function
        if x< first then func1(x)
        elif x< second then func2(x) 
        else func3(x)
        '''
    if x < first : 
        return func1(x) 
    elif x<second: 
        return func2(x) 
    else :
        return func3(x)