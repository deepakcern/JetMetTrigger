#!/usr/bin/env python
from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TLorentzVector, AddressOf, gROOT, TNamed,gStyle
import ROOT as ROOT
import os
import sys, optparse
from array import array
import math
import numpy as numpy_


## This setup has been tested only for slc6 for centos7 please contact, work is in progress. 
ntuple = TChain("pfRaw/events")
ntuple.Add("hltJetMetNtuple.root")
NEntries = ntuple.GetEntries()

depth1=[]
depth2=[]
depth3=[]
depth4=[]
depth5=[]
depth6=[]
depth7 = []
for ievent in range(NEntries):
    if ievent%100==0: print "Processed %d of %d events..." %(ievent,NEntries)
    ntuple.GetEntry(ievent)

    depthFractions  = ntuple.__getattr__('depthFractions')
    #print "lengh of condidates", len(depthFractions)
    for cnd in range(len(depthFractions)):
        fractions = depthFractions[cnd]
        fractions_ = []
        fractions_ = [i for i in fractions]
        if 1.0 in fractions_ or sum(fractions_)==0.0: continue

        #print fractions_
        dep1=fractions_[0]
        dep2=fractions_[1]
        dep3=fractions_[2]
        dep4=fractions_[3]
        dep5=fractions_[4]
        dep6=fractions_[5]
        dep7=fractions_[6]

        #print dep1,dep2,dep3,dep4
        depth1.append(dep1)
        depth2.append(dep2)
        depth3.append(dep3)
        depth4.append(dep4)
        depth5.append(dep5)
        depth6.append(dep6)
        depth7.append(dep7)
        #print "depth=================+'\n'"
        #print depth1,depth2
'''
print ""
print "depth1",depth1
print "depth2",depth2
print "depth3",depth3
print "depth4",depth4
print "depth5",depth5
print "depth6",depth6
print "depth7",depth7

print "sum(depth1)",sum(depth1)
print "sum(depth2)",sum(depth2)
print "sum(depth3)",sum(depth3)
print "sum(depth4)",sum(depth4)
print "sum(depth5)",sum(depth5)
print "sum(depth6)",sum(depth6)
print "sum(depth7)",sum(depth7)
'''
frac1 = sum(depth1)/len(depth1)
frac2 = sum(depth2)/len(depth2)
frac3 = sum(depth3)/len(depth3)
frac4 = sum(depth4)/len(depth4)
frac5 = sum(depth5)/len(depth5)
frac6 = sum(depth6)/len(depth6)
frac7 = sum(depth7)/len(depth7)


fractions = [frac1,frac2,frac3,frac4,frac5,frac6,frac7]
print fractions
gStyle.SetOptStat(0)
h_depthFrac = TH1F(" Average EnergyFraction"," Avarage EnergyFraction",7,0,7)

lab = ['Depth1','Depth2','Depth3','Depth4','Depth5','Depth6','Depth7','Depth7']

for i in range(len(fractions)):
    h_depthFrac.SetBinContent(i+1,fractions[i])
    h_depthFrac.GetXaxis().SetBinLabel(i+1,lab[i])

#h_depthFrac.SetXTitle("MET [GeV]")
h_depthFrac.SetYTitle("Energy fraction per depth")

c=TCanvas()

h_depthFrac.Draw('hist')
c.SaveAs('depthFraction_average.png')
c.SaveAs('depthFraction_average.pdf')


