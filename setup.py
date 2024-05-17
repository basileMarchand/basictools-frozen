# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

import os, sys
from setuptools import setup
from distutils.command.build import build
import configparser

'''
PREFIX : Set this variable to point to the external libraries (if the mkl or eigen are installed with pip install --user for example)
'''



__config = configparser.ConfigParser()
__config.read('setup.cfg')

if __name__ == '__main__':
    setup(
        data_files=[("ParaViewPlugins",["extras/BasicToolsParaViewBridge.py"])]
)
