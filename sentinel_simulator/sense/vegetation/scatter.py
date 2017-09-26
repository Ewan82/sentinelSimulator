#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Filename: scatter.py
# Author: "Thomas Weiß"
# Date:   2017-07-26 13:46:39
# Last Modified by:   "Thomas Weiß"
# Last Modified time: 2017-07-26 13:48:55

"""
Major vegetation scatter class
"""

class VegetationScatterWaterCloudModel(object):
    def __init__(self, theta, V=None, V1=None, V2=None, A=None, B=None):
        self.theta = theta
        self.V = V
        self.V1 = V1
        self.V2 = V2
        self.A = A
        self.B = B

