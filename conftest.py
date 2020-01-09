# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

import pytest


class fail_if_not_ok:
    def __init__(self, f):
        self.f = f

    def __call__(self):
        res = self.f()
        assert res.lower() == "ok"

def pytest_pycollect_makeitem(collector, name, obj):
    if name == "CheckIntegrity" and hasattr(obj, "__call__"):
        # Collect CheckIntegrity functions
        # Decorate them to detect failing return values
        return pytest.Function(name, parent=collector, callobj=fail_if_not_ok(obj))
    else:
        # Fallback to default handler
        return None
