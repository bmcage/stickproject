# Low level utilities in stickproject under BSD license!
#
# Copyright (C) 2011  B. Malengier
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#  a. Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#  b. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#  c. Neither the name of the Enthought nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#

import numpy as np

def fullcompare_array(a, b=None, func=None, funcdata_a=None,
                      funcdata_b=None, chunk_size=None):
    """
    Compare every element of a with every element from b
    
    Given two 2d arrays (m,n) interpreted as m points of dimension n, compare
    every point with each other using a given vectorized function func that 
    uses user data funcdata_a and funcdata_b.
    
    Warning: this method needs a lot of memory ma*na*mb+mb*nb*ma. 
    You can limit the memory usage by giving a chunk_size, indicating how many
    elements of a (ma) are compared with b at the same time. 
    
    Parameters
    ----------
    a : 2D array
        `a` represents ma points of dimension na
    b : 2D array, optional
        `b` represents mb points of dimension nb. nb=na is required
        if b is not given, a will be compared with itself
    func : a callable, optional
        `func` must have the signature func(x, y, data_x, data_y) where 
        `x` will be a slice of `a` and `y` a slice of `b`; `data_x` is the 
        corresponding slice of `funcdata_a`, and `data_y` is `funcdata_b`. The 
        data of `data_y` relevant for `y[i]` is `data_y[i%len(a)]`
        If `func` is not given, the return value is where the arrays are equal
        The return format of `func` should be a tuple `res1`, `res2` where 
        `res1` is a bool array of where the comparison has a result of which the
        output must be kept, and `res2` is a list of arrays with computed 
        results.
    funcdata_a : array-like, optional
        Data that will be passed to func for the points being compared
    funcdata_b : array-like, optional
        Data that will be passed to func for the points being compared
    chunk_size : int > 0, optional
        To reduce memory footprint, a can be divided in chunks to do the 
        comparison with b.
    
    Returns
    -------
    ind : list of integer arrays
        ind[m] contains the indices in b that func returns True with for a[m]
    res : results of comparison
        `res` has the format of the result as given by `func`. This result is
        only returned for those values given in `ind`
        
    Example
    -------
    >>> x = np.array([[3],[2],[1],[3],[5],[6]])
    
    Now use fullcompare_array to find that first and fourth are equal
    >>> ind, res = fullcompare_array(x)

    """
    #some argument checking
    if len(a.shape) != 2:
        raise ValueError, "a should be (m,n) matrix, m number of points, dim n"
    if b is None:
        b = a
        funcdata_b = funcdata_a
    if len(b.shape) != 2:
        raise ValueError, "b should be (m,n) matrix, m number of points, dim n"
    if b.shape[-1] != a.shape[-1]:
        raise ValueError, "dimension of points in b should be the same as a"
    
    nrptna = a.shape[0]
    nrptnb = b.shape[0]

    ind = [None] * nrptna
    res = None

    if func is None:
        def _equal(a, b, dataa, datab):
            c = np.sum(a-b, axis=1) #square of distance
            return c == 0, [c,]
        func = lambda x1,x2,x3,x4: _equal(x1,x2,x3,x4)
    if chunk_size is None:
        chunk_size = a.shape[0]
        nr_chunks = 1
    else:
        nr_chunks = int(round(nrptna / chunk_size + 0.5-1e-10))
    tmp_funcdata_a = None
    indices = np.arange(nrptna)
    #repeat b as many times as chunk_size
    tmp_b_repeat = np.repeat(b.T, chunk_size, axis=0).reshape(
                    (b.shape[1], chunk_size*b.shape[0])).T
    tmp_funcdata_a = None
    for chunk in range(nr_chunks):
        chunk_indices = indices[chunk*chunk_size:(chunk+1)*chunk_size]
        stchk = chunk*chunk_size
        if chunk_size != len(chunk_indices):  #last length may be different
            chunk_size = len(chunk_indices)
            #repeat b as many times as chunk_size
            tmp_b_repeat = np.repeat(b.T, chunk_size, axis=0).reshape(
                    (b.shape[1], chunk_size*b.shape[0])).T
        tmp_a = a.T[..., chunk_indices]
        if funcdata_a is not None:
            tmp_funcdata_a = funcdata_a[chunk_indices]
        #repeat tmp_a as many times as nrptnb 
        tmp_a_repeat  = np.repeat(tmp_a, nrptnb, axis=1).T
        tind, tres = func(tmp_a_repeat, tmp_b_repeat, tmp_funcdata_a, funcdata_b)
        tind = tind.reshape((chunk_size, nrptnb))
        for i,j in enumerate(tres):
            tres[i] = j.reshape((chunk_size, nrptnb))
        if res is None:
            #first time, size as output of func
            res = [None] * len(tres)
            for i in range(len(tres)):
                res[i] = [None] * nrptna
        for chk in np.arange(chunk_size):
            ind[stchk + chk] = indices[tind[chk]]
            for i,j in enumerate(tres):
                res[i][stchk + chk] = j[chk][tind[chk]]
    return ind, res

def circledist(a, b, rada, radb):
    """
    Compute distance between a and b and return those that do not overlap. 
    a, b are arrays of length len(rada) * len(radb). a is repeated radb times
    """
    #substract both and calc distance of all pairs
    tmp = a - b
    tmp = np.sum(tmp*tmp, axis=1) #square of distance
    tmp = tmp.reshape((len(rada), len(radb)))
    distreqsqr = np.empty(tmp.shape, float)
    #print 'tmp', tmp
    for i in range(len(rada)):
        distreqsqr[i,:] = np.power(rada[i] + radb[:] + 0.001, 2)
    #now determine all indices that overlap
    result = (tmp < distreqsqr) # point itself is returned too!
    return result.flatten(), [np.sqrt(tmp), np.sqrt(distreqsqr)]

def test():
    x = np.array([[3],[2],[1],[3],[5],[6]])
    #Now use fullcompare_array to find that first and fourth are equal
    ind, res = fullcompare_array(x)
    print ind
    raw_input("Continue")
    
    a = np.zeros((3,1), float)
    a[0,0] = 1; a[1,0]=2; a[2,0]=1
    print fullcompare_array(a)
    raw_input("Continue")
    
    a = np.zeros((4,2), float)
    a[0,0] = 1; a[1,0]=2; a[2,0]=1; a[3,0]=2
    a[0,1] = 4; a[1,1]=5; a[2,1]=4; a[3,1]=6
    print fullcompare_array(a)
    raw_input("Continue")
    

if __name__ == '__main__': 
    test()

