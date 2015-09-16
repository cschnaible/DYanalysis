'''
Usage:python plot.py RootFile.root label[optional]
Script to make migration plots of RECO vs GEN.
Author: C. Schnaible, UC Los Angeles
'''

from subprocess import Popen
from sys import argv, exit, stdout, stderr
import math
from array import array

import ROOT

# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)

######## File #########
if len(argv) < 2:
   print 'Usage:python plot.py RootFile.root label[optional]'
   exit()

option = argv[1]
if option == "25ns":
    DY_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DYJetsToLL_M-25_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_25ns/DYJetsToLL_25ns.root")
    TTa_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_25ns/TTJets_nlo_25ns.root")
    TTb_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MC_25ns/TTJets_MLM_25ns.root")
    WW_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WW_TuneCUETP8M1_13TeV-pythia8/MC_25ns/WW_25ns.root")
    WZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WZ_TuneCUETP8M1_13TeV-pythia8/MC_25ns/WZ_25ns.root")
    ZZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/ZZ_TuneCUETP8M1_13TeV-pythia8/MC_25ns/ZZ_25ns.root")
    #data_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/Run2015B/DoubleMuon_2015RunB.root")
if option == "50ns":
    DY_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_50ns/DYJetsToLL_50ns.root")
    TTa_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_50ns/tt_nlo_50ns.root")
    TTb_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MC_50ns/tt_MLM_50ns.root")
    WW_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WW_TuneCUETP8M1_13TeV-pythia8/MC_50ns/WW_50ns.root")
    WZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WZ_TuneCUETP8M1_13TeV-pythia8/MC_50ns/WZ_50ns.root")
    ZZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/ZZ_TuneCUETP8M1_13TeV-pythia8/MC_50ns/ZZ_50ns.root")
    data_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/Run2015B/DoubleMuon_2015RunB.root")
outputName = '~/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/DYmigration.root'
outputRoot = ROOT.TFile(outputName,"recreate")

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/CMSSW_7_4_6_patch6/src/DYanalysis/plots/'+argv[2]+'_'
else:
   saveWhere='~/CMSSW_7_4_6_patch6/src/DYanalysis/plots/'



#####################################
#Get Effi NTUPLE                 #
#####################################

DYtree = DY_file.Get("analyzer/events")
datatree = data_file.Get("analyzer/events")

canvas = ROOT.TCanvas()

def make_migration(tree, recoVar, genVar, recosel, gensel, binning, xaxis, yaxis, title):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s:%s>>htemp(%s)" % (recoVar, genVar, ", ".join(str(x) for x in binning2))
    tree.Draw(draw_string, gensel+"||"+recosel, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.GetYaxis().SetTitle(yaxis)
    output_histo.SetTitle(title)
    return output_histo

def make_plot(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo
    
def produce_migration(tree,recoVar,genVar,recosel,gensel,binning1,binning2,xaxis,yaxis,title):
    migPlot = make_migration(tree,recoVar,genVar,recosel,gensel,binning2,xaxis,yaxis,title)
    genPlot = make_plot(tree,genVar,gensel,binning1)
    for i in range(genPlot.GetNbinsX()):
        G = genPlot.GetBinContent(i+1)
        for j in range(genPlot.GetNbinsX()):
            M = migPlot.GetBinContent(i+1,j+1)
            mig = M/G
            migPlot.SetBinContent(i+1,j+1,mig)
    return migPlot

def save_hist(tree,recoVar,genVar,recosel,gensel,binning1,binning2, xaxis, yaxis, filename):
    #frame = ROOT.TH2F("frame","frame",30,60,120,30,60,120)
    #frame.GetXaxis().SetTitle(xaxis)
    #frame.GetYaxis().SetTitle(yaxis)
    #frame.SetTitle(title)
    #frame.Draw()
    migPlot = produce_migration(tree,recoVar,genVar,recosel,gensel,binning1,binning2,xaxis,yaxis,title)
    migPlot.Draw("colz")
    migPlot.SetName(filename)
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.Write(filename)
    canvas.SaveAs(saveas)

recoVar = "llM"
genVar  = "llMGEN"
recosel = "(lepPt[0]>20 && lepPt[1]>20 && fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4 && nLeptons>1 && llchid==-4 && lepId[0]==1 && lepId[1]==1 && lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15 && isTriggered==1)"
gensel = "(lepPtGEN[0]>20 && lepPtGEN[1]>10 && fabs(lepEtaGEN[0])<2.4 && fabs(lepEtaGEN[1])<2.4 && fabs(lepPdgIdGEN[0])==13 && fabs(lepPdgIdGEN[1])==13 && nGenLeptons>1)"
xaxis = "Generated M(#mu#mu) [GeV]"
yaxis = "Reconstructed M(#mu#mu) [GeV]"
title = "Mass Migration"
filename = "massMigrationPlot"
binning1 = [30,60,120]
binning2 = [30,60,120,30,60,120]
save_hist(DYtree,recoVar,genVar,recosel,gensel,binning1,binning2,xaxis,yaxis,filename)
