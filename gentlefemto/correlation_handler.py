"""
Module containing the class used to compute the correlation function 
from same-event and mixed-event distributions.
"""

from ROOT import TH1F
import logging
import numpy as np


class CorrelationHandler:
    """
    Class used to compute the correlation function 
    from same-event and mixed-event k* distributions.
    Performs rebinning, normalisation and advanced operations.

    Parameters:
    ____________________________________________________

    same_event = same-event k* distribution (TH1F)

    mixed_event = same-event k* distribution (TH1F)

    """

    def __init__(self, same_event=None, mixed_event=None) -> None:
        self.same_event = same_event
        self.mixed_event = mixed_event
        self.correlation_function = None

    def self_normalise(self) -> None:
        if self.same_event:
            n_same = self.same_event.GetNbinsX()
            integral_same = self.same_event.Integral(1, n_same+1)
            self.same_event.Scale(1/integral_same)
        else:
            logging.info('Same-event not set.')
        if self.mixed_event:
            n_mixed = self.mixed_event.GetNbinsX()
            integral_mixed = self.mixed_event.Integral(1, n_mixed+1)
            self.mixed_event.Scale(1/integral_mixed)
        else:
            logging.info('Mixed-event not set.')

    def set_same_event(self, same_event=None) -> None:
        self.same_event = same_event

    def set_mixed_event(self, mixed_event=None) -> None:
        self.mixed_event = mixed_event

    def rebin(self, n_rebin=2) -> None:
        self.same_event.Rebin(n_rebin)
        self.mixed_event.Rebin(n_rebin)

    def make_correlation_function(self, name = 'hCF'):
        self.correlation_function = self.same_event.Clone(name)
        self.correlation_function.Reset()
        self.correlation_function.GetYaxis().SetTitle('C(k*)')
        for i in range(1, self.correlation_function.GetNbinsX()):
            vSame = self.same_event.GetBinContent(i)
            eSame = self.same_event.GetBinError(i)
            vMixed = self.mixed_event.GetBinContent(i)
            eMixed = self.mixed_event.GetBinError(i)
            if vSame < 1.e-12 or vMixed < 1.e-12:
                continue
            else:
                vCF = vSame/vMixed
                eCF = np.sqrt((eMixed/vMixed)*(eMixed/vMixed) +
                              (eSame/vSame)*(eSame/vSame)) * vCF
                self.correlation_function.SetBinContent(i, vCF)
                self.correlation_function.SetBinError(i, eCF)

    def get_correlation_function(self):
        return self.correlation_function

    def get_same_event(self):
        return self.same_event

    def get_mixed_event(self):
        return self.mixed_event

    def normalise(self, low=0.6, high=1.) -> None:
        if self.correlation_function:
            low_bin = self.correlation_function.FindBin(low)
            low_edge = self.correlation_function.GetBinLowEdge(low_bin)
            high_bin = self.correlation_function.FindBin(high)
            high_edge = self.correlation_function.GetBinLowEdge(high_bin+1)
            cf_integral = self.correlation_function.Integral(
                low_bin, high_bin, 'width')
            self.correlation_function.Scale((high_edge-low_edge)/cf_integral)
        else:
            logging.error('Compute the correlation function first')
