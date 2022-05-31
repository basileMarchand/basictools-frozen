#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import sys
if sys.platform == "win32":
    from numpy import int32 as PBasicIndexType
else:
    from numpy import int64 as PBasicIndexType

from numpy import float64 as PBasicFloatType

from numpy.typing import ArrayLike
