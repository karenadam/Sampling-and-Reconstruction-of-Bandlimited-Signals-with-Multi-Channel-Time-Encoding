import numpy as np
import SimulationSettings


##### IMPORTING PATHS #####
import sys
import os

sys.path.insert(0, os.path.split(os.path.realpath(__file__))[0] + "/../Multi-Channel-Time-Encoding/Source")
Figure_Path = os.path.split(os.path.realpath(__file__))[0] + "/../Figures/"
Data_Path = os.path.split(os.path.realpath(__file__))[0] + "/../Data/"



##### IMPORT FILES #####
from Time_Encoder import *
from Signal import *

##### MATPLOTLIB SETTINGS #####
if SimulationSettings.graphical_import:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rc

    if SimulationSettings.To_Svg:
        plt.rc('text', usetex=False)
        plt.rc('text.latex', unicode = False)
        plt.rc('svg',fonttype = 'none')
    else:
        matplotlib.rc("text", usetex=True)
        matplotlib.rc("font", family="serif")
        matplotlib.rc("font", size=7)
        matplotlib.rc("text.latex", preamble=r"\usepackage{amsmath}\usepackage{amssymb}")
    from matplotlib.colors import LogNorm   

    import seaborn as sns

import time
import pickle
import datetime
import math
import multiprocessing
import itertools
import csv


def getdatetime():
    return "{date:%Y-%m-%d %H-%M-%S}".format(date=datetime.datetime.now())


def getmostrecentfile(logdir):
    return max([f for f in os.listdir(logdir)])
