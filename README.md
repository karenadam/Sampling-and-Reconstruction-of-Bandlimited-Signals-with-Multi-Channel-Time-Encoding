# Sampling and Reconstruction of Bandlimited Signals with Multi-Channel Time Encoding
This repository holds the code for reproducing the figures in the paper "Sampling and Reconstruction of Bandlimited Signals with Multi-Channel Time Encoding", by Karen Adam, Adam Scholefield and Martin Vetterli

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/karenadam/Sampling-and-Reconstruction-of-Bandlimited-Signals-with-Multi-Channel-Time-Encoding/master?filepath=Code%2FGenerate%20Paper%20Figures.ipynb)

# Abstract
    Sampling is classically performed by recording the amplitude of the input at given time instants; however, sampling and reconstructing a signal using multiple devices in parallel becomes a more difficult problem to solve when the devices have an unknown shift in their clocks.
    Alternatively, one can record the times at which a signal (or its integral) crosses given thresholds. This can model integrate-and-fire neurons, for example, and has been studied by Lazar and Tòth under the name of "Time Encoding Machines". This sampling method is closer to what is found in nature.
    In this paper, we show that, when using time encoding machines, reconstruction from multiple channels has a more intuitive solution, and does not require the knowledge of the shifts between machines. We show that, if single-channel time encoding can sample and perfectly reconstruct a 2Ω-bandlimited signal, then M-channel time encoding can sample and perfectly reconstruct a signal with M times the bandwidth.
    Furthermore, we present an algorithm to perform this reconstruction and prove that it converges to the correct unique solution, in the noiseless case, without knowledge of the relative shifts between the machines. This is quite unlike classical multi-channel sampling, where unknown shifts between sampling devices pose a problem for perfect reconstruction.

# Reference
https://arxiv.org/abs/1907.05673
