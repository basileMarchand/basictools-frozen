#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

IF UNAME_SYSNAME == "Linux":
    from numpy cimport int64_t as CBasicIndexType
ELIF UNAME_SYSNAME == "Darwin":
    from numpy cimport int64_t as CBasicIndexType
ELSE:
    from numpy cimport int32_t as CBasicIndexType

from numpy cimport float64_t as CBasicFloatType
