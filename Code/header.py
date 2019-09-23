import numpy as np
import SimulationSettings


##### IMPORTING PATHS #####
import sys

sys.path.insert(0, "../Multi-Channel-Time-Encoding/Source")

##### IMPORT FILES #####
from Time_Encoder import *
from Signal import *

##### MATPLOTLIB SETTINGS #####
if SimulationSettings.graphical_import:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rc

    matplotlib.rc("text", usetex=True)
    matplotlib.rc("font", family="serif")
    matplotlib.rc("font", size=7)
    matplotlib.rc("text.latex", preamble=r"\usepackage{amsmath}\usepackage{amssymb}")

    # import matplotlib2tikz

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
