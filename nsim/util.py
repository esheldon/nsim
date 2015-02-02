from __future__ import print_function
import numpy

def ring_select(logic):
    ii = numpy.arange(logic.size)
    wodd,  = numpy.where( (ii % 2) != 0)
    weven, = numpy.where( (ii % 2) == 0)

    wboth, = numpy.where(  logic[wodd] & logic[weven] )

    w=numpy.concatenate( (wodd[wboth], weven[wboth] ) )

    return w
