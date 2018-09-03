# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 15:39:53 2018

@author: catar
"""

from scipy.io import loadmat
data = loadmat('InstanciaArtigo_II.mat')
C = data['C'][0][0]
R = data['R'][0][0]
