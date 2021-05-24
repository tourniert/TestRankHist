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

        self.chisquare = ResultsJolliffePrimo(np.sum(self.__X**2),
                                              1-chi2.cdf(np.sum(self.__X**2),
                                                         df=self.__shape-1))


    def linear(self):
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
        I_lin = self.__normalizeI(I_lin)
        U_linear = np.sum(I_lin * self.__X)**2

        return ResultsJolliffePrimo(U_linear, 1-chi2.cdf(U_linear, df=1))


    def ends(self):
        """
        Apply the Chi² test on ends vector.

        Returns
        -------
        ResultsJolliffePrimo
            Object with 2 attributes: statistic and p_value.

        """
        if self.__even:
            a = self.__h - 1
            b = 1
        else:
            a = 2*self.__h - 1
            b = 2

        I_ends = [-b for _i in range(self.__shape)]
        I_ends[0] = a
        I_ends[-1] = a
        I_ends = self.__normalizeI(I_ends)
        U_ends = np.sum(I_ends * self.__X)**2

        return ResultsJolliffePrimo(U_ends, 1-chi2.cdf(U_ends, df=1))


    def Vshape(self):
        """
        Apply the Chi² test on Vshape vector.

        Returns
        -------
        ResultsJolliffePrimo
            Object with 2 attributes: statistic and p_value.

        """
        if self.__even:
            a = self.__h - 1
            b = 2
            I_Vshape = [a for _i in range(self.__shape//2)]
            for _i in range(self.__shape//2):
                I_Vshape[_i] -= _i*b
            I_Vshape_rev = list(np.flip(I_Vshape))
            I_Vshape = I_Vshape + I_Vshape_rev
        else:
            a = self.__h**2
            b = 2*self.__h + 1
            I_Vshape = [a for _i in range(self.__shape//2 + 1)]
            for _i in range(self.__shape//2 + 1):
                I_Vshape[_i] -= _i*b
            I_Vshape_rev = list(np.flip(I_Vshape[:-1]))
            I_Vshape = I_Vshape + I_Vshape_rev

        I_Vshape = self.__normalizeI(I_Vshape)
        U_Vshape = np.sum(I_Vshape * self.__X)**2

        return ResultsJolliffePrimo(U_Vshape, 1-chi2.cdf(U_Vshape, df=1))


    def Ushape(self):
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

        I_Ushape = self.__normalizeI(I_Ushape)
        U_Ushape = np.sum(I_Ushape * self.__X)**2

        return ResultsJolliffePrimo(U_Ushape, 1-chi2.cdf(U_Ushape, df=1))


    def wave(self):
        """
        Apply the Chi² test on wave vector.

        Returns
        -------
        ResultsJolliffePrimo
            Object with 2 attributes: statistic and p_value.

        """
        I_wave = [np.sin(2*np.pi * (_i/(self.__shape-1))) 
                  for _i in range(self.__shape)]

        I_wave = self.__normalizeI(I_wave)
        U_wave = np.sum(I_wave * self.__X)**2

        return ResultsJolliffePrimo(U_wave, 1-chi2.cdf(U_wave, df=1))


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
        p_chisq = round(self.chisquare.p_value, 3)
        p_linear = round(self.linear().p_value, 3)
        p_ends = round(self.ends().p_value, 3)
        p_Vshape = round(self.Vshape().p_value, 3)
        p_Ushape = round(self.Ushape().p_value, 3)
        p_wave = round(self.wave().p_value, 3)

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
        # ax.set_title('{}'.format(p_chisq))
        ax.set_title('$\chi^2$: {}, Linear: {}, '.format(p_chisq, p_linear) +
                     'Ends : {}, \n Vshape: {}, '.format(p_ends, p_Vshape) +
                     'Ushape: {}, Wave: {}'.format(p_Ushape, p_wave))


    def __normalizeI(self, I):
        I = np.asarray(I)
        sum_squared = np.sum(I**2)

        return I / np.sqrt(sum_squared)


class ResultsJolliffePrimo:
    def __init__(self, statistic, p_value):
        self.statistic = statistic
        self.p_value = p_value
