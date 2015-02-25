from __future__ import print_function
import numpy

def ring_select(logic):

    weven = numpy.arange(0,logic.size,2)
    wodd  = numpy.arange(1,logic.size,2)

    wboth, = numpy.where(  logic[wodd] & logic[weven] )

    w=numpy.concatenate( (wodd[wboth], weven[wboth] ) )

    w.sort()
    return w
