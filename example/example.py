# -*- coding: utf-8 -*-
"""
Created on Sat May 22 2021

@author: tourniert
"""
import matplotlib.pyplot as plt
from TestRankHist import TestRankHist


# Definition of histograms
hist_slope  = [3, 6, 3, 7, 4, 6, 4, 4, 2, 4, 1, 2, 3, 6, 2, 3]
hist_linear = [1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 6, 6, 6, 7]
hist_convex = [2, 2, 3, 3, 4, 4, 6, 7, 6, 6, 4, 4, 3, 3, 2, 1]
hist_wave   = [4, 4, 6, 7, 6, 6, 4, 4, 3, 3, 2, 1, 2, 2, 3, 3]


# Plotting figure
fig, axs = plt.subplots(2, 2, figsize=(14, 14))
print(axs[0, 0])
# Plotting with special x_label
TestRankHist(hist_slope).plot(x=[str(_i) for _i in range(16)],
                                   ax=axs[0, 0])
TestRankHist(hist_linear).plot(ax=axs[0, 1])
TestRankHist(hist_convex).plot(ax=axs[1, 0])
TestRankHist(hist_wave).plot(ax=axs[1, 1])