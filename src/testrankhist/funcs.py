# -*- coding: utf-8 -*-
"""
Created on Thu May 20 2021

@author: tourniert
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2


class TestRankHist:
    def __init__(self, histogram):
        self.hist = histogram

        if len(np.shape(self.hist)) > 1:
            raise ValueError('histogram must only have one dimension, but ' +
                             'have shape {}'.format(np.shape(self.hist)))
        elif np.shape(self.hist)[0] == 0:
            raise ValueError('histogram must not be empty')

        self.__shape = len(self.hist)
        self.__even = (self.__shape % 2 == 0)
        self.__h = self.__shape // 2

        self.__E = (np.ones(self.__shape)*np.sum(self.hist)) / self.__shape
        self.__X = (self.hist - self.__E) / np.sqrt(self.__E)

        # Calculation of chi2
        self.chisquare = ResultsJolliffePrimo(np.sum(self.__X**2),
                                              self.__shape-1)

        # Creation of a non-orthonormal vector base
        self.__basis = [self.__vec_linear(),
                        self.__vec_Ushape(),
                        self.__vec_wave()]

        # Apply Gram-Schmidt algorithm
        self.__ortho_linear, self.__ortho_Ushape, self.__ortho_wave = self.__gramSchmidt()

        # Calculation of statistics of chi2 component
        self.__U_linear = np.sum(self.__ortho_linear * self.__X)**2
        self.__U_Ushape = np.sum(self.__ortho_Ushape * self.__X)**2
        self.__U_wave = np.sum(self.__ortho_wave * self.__X)**2

        # Calculation of final results
        self.linear = ResultsJolliffePrimo(self.__U_linear, 1)
        self.Ushape = ResultsJolliffePrimo(self.__U_Ushape, 1)
        self.wave = ResultsJolliffePrimo(self.__U_wave, 1)

    def __vec_linear(self):
        """
        Apply the Chi² test on linear vector.

        Returns
        -------
        ResultsJolliffePrimo
            Object with 2 attributes: statistic and p_value.

        """
        if self.__even:
            a = -(2*self.__h - 1)
            b = 2
        else:
            a = -self.__h
            b = 1

        I_lin = [a+_i*b for _i in range(self.__shape)]
        I_lin = self.__normalize(I_lin)
        return I_lin

    def __vec_Ushape(self):
        """
        Apply the Chi² test on Ushape vector.

        Returns
        -------
        ResultsJolliffePrimo
            Object with 2 attributes: statistic and p_value.

        """
        if self.__even:
            b = (4*(self.__h**2) - 1) / 3
            I_Ushape = [-b for _i in range(self.__shape//2)]
            for _i in range(self.__shape//2):
                I_Ushape[_i] += (2*self.__h - (2*_i+1))**2
            I_Ushape_rev = list(np.flip(I_Ushape))
            I_Ushape = I_Ushape + I_Ushape_rev
        else:
            b = (self.__h*(self.__h + 1)) / 3
            I_Ushape = [-b for _i in range(self.__shape//2 + 1)]
            for _i in range(self.__shape//2 + 1):
                if _i // self.__h == 0:
                    I_Ushape[_i] += (self.__h - _i)**2
            I_Ushape_rev = list(np.flip(I_Ushape[:-1]))
            I_Ushape = I_Ushape + I_Ushape_rev

        I_Ushape = self.__normalize(I_Ushape)
        return I_Ushape

    def __vec_wave(self):
        """
        Apply the Chi² test on wave vector.

        Returns
        -------
        ResultsJolliffePrimo
            Object with 2 attributes: statistic and p_value.

        """
        I_wave = [np.sin(2*np.pi * (_i/(self.__shape-1)))
                  for _i in range(self.__shape)]

        I_wave = self.__normalize(I_wave)
        return I_wave
        # U_wave = np.sum(I_wave * self.__X)**2

        # return ResultsJolliffePrimo(U_wave, 1-chi2.cdf(U_wave, df=1))

    def plot(self, x=None, ax=None, **fig_kw):
        """
        Plot the rank histogram with the scores for each test.

        Parameters
        ----------
        x : list, optional
            list of string for plotting on x axis.
            The default is ['1', '2', ..., 'len(hist)'].
        ax : matplotlib.pyplot.AxesSubplot, optional
            Axis on which we plot the histogram. The default is to plot on a
            new figure.
        **fig_kw : dictionnary
            Arguments to give to matplotlib.pyplot.subplots().

        Returns
        -------
        Axis or figure with ranks histogram.

        """
        if ax is None:
            fig, ax = plt.subplots(**fig_kw)

        if x is None:
            x = [str(_i) for _i in range(1, self.__shape+1)]

        # Plotting
        ax.bar(x, self.hist, color='grey')
        ax.axhline(y=self.__E[0], color='k', ls='--')
        ax.set_ylabel('Counts')
        ax.set_xlabel('Rank')
        ax.set_xticks(x[::2])
        ax.set_title('$\chi^2$: {}, '.format(self.chisquare.p_value) +
                     'Linear: {}, '.format(self.linear.p_value) +
                     '\n Ushape: {}, '.format(self.Ushape.p_value) +
                     'Wave: {}'.format(self.wave.p_value))

    def __normalize(self, vec):
        vec = np.asarray(vec)
        sum_squared = np.sum(vec**2)
        return vec / np.sqrt(sum_squared)

    def __gramSchmidt(self):
        basis = []
        for _i, _vec in enumerate(self.__basis):
            sum_proj = 0
            for _j in range(_i):
                sum_proj += np.dot(_vec, basis[_j])*basis[_j]
            e = (_vec - sum_proj) / np.linalg.norm(_vec - sum_proj)
            basis.append(e)
        return basis


class ResultsJolliffePrimo:
    def __init__(self, U, df):

        self.statistic = round(U, 3)
        self.p_value = round(1-chi2.cdf(self.statistic, df=df), 3)
