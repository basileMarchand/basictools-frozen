#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

from numpy import int64, float
from numpy cimport int64_t, float64_t

int_DTYPE   = int64
float_DTYPE = float

ctypedef int64_t     int_DTYPE_t
ctypedef float64_t float_DTYPE_t