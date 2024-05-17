# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

import os, sys
from setuptools import setup
import configparser

__config = configparser.ConfigParser()
__config.read('setup.cfg')

if __name__ == '__main__':
    setup(
        data_files=[("ParaViewPlugins",["extras/BasicToolsParaViewBridge.py"])]
)
