#!/usr/bin/env python
from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TLorentzVector, AddressOf, gROOT, TNamed, gStyle
import ROOT as ROOT
import os
import sys, optparse
from array import array
import math
import numpy as numpy_

gStyle.SetOptStat(0)
h_depthFrac = TH1F("EnergyFraction","EnergyFraction",7,0,7)

arr = [0.6526692509651184, 0.20300354063510895, 0.1383916437625885, 0.005935580935329199, 0.0, 0.0, 0.0] 
lab = ['Depth1','Depth2','Depth3','Depth4','Depth5','Depth6','Depth7','Depth7']

for i in range(len(arr)):
    h_depthFrac.SetBinContent(i+1,arr[i])
    h_depthFrac.GetXaxis().SetBinLabel(i+1,lab[i])

#h_depthFrac.SetXTitle("MET [GeV]")
h_depthFrac.SetYTitle("Energy fraction per depth")

c=TCanvas()

h_depthFrac.Draw('hist')
c.SaveAs('depthFraction.png')
c.SaveAs('depthFraction.pdf')
