'''
Usage:python plot.py RootFile.root label[optional]
Script to plot N-1 and cut-by-cut efficiency plots
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
ROOT.gStyle.SetLegendBorderSize(0)


######## File #########
if len(argv) < 2:
   print 'Usage:python plot.py RootFile.root label[optional]'
   exit()

option = argv[1]
#if option == "25ns":
    #DY_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DYJetsToLL_M-25_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_25ns/DYJetsToLL_25ns.root")
    #TTa_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_25ns/TTJets_nlo_25ns.root")
    #TTb_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MC_25ns/TTJets_MLM_25ns.root")
    #WW_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WW_TuneCUETP8M1_13TeV-pythia8/MC_25ns/WW_25ns.root")
    #WZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WZ_TuneCUETP8M1_13TeV-pythia8/MC_25ns/WZ_25ns.root")
    #ZZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/ZZ_TuneCUETP8M1_13TeV-pythia8/MC_25ns/ZZ_25ns.root")
    #data_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/Run2015B/DoubleMuon_2015RunB.root")
if option == "50ns":
    DY_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_50ns/DYJetsToLL_50ns.root")
    #TTa_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_50ns/tt_nlo_50ns.root")
    #TTb_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MC_50ns/tt_MLM_50ns.root")
    #WW_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WW_TuneCUETP8M1_13TeV-pythia8/MC_50ns/WW_50ns.root")
    #WZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WZ_TuneCUETP8M1_13TeV-pythia8/MC_50ns/WZ_50ns.root")
    #ZZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/ZZ_TuneCUETP8M1_13TeV-pythia8/MC_50ns/ZZ_50ns.root")
    singleMu_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/Run2015B/DoubleMuon_Run2015B.root")
    doubleMu_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/Run2015B/DoubleMuon_Run2015B.root")
outputName = '~/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/DYNm1andSeq.root'
outputRoot = ROOT.TFile(outputName,"update")

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/CMSSW_7_4_6_patch6/src/DYanalysis/plots/'+argv[2]+'_'
else:
   saveWhere='~/CMSSW_7_4_6_patch6/src/DYanalysis/plots/'



#####################################
#Get Effi NTUPLE                 #
#####################################

DYtree = DY_file.Get("analyzer/events")
singleMuTree = singleMu_file.Get("analyzer/events")
doubleMuTree = doubleMu_file.Get("analyzer/events")

canvas = ROOT.TCanvas()#, 800, 800)

def getNentries(tree,selection):
    N = tree.GetEntries(selection)
    return N

def make_efficiency(denom, num):
    ''' Make an efficiency/acceptance graph with the style '''
    eff = ROOT.TGraphAsymmErrors(num, denom)
    eff.SetMarkerStyle(20)
    eff.SetMarkerSize(1.0)
    return eff

def make_denom(tree,denomSel):
    denom = ROOT.TH1F("denom","",nVar,0,nVar)
    for i in range(denom.GetNbinsX()):
        ibin = i+1
        numD = getNentries(tree,denomSel[i])
        denom.SetBinContent(ibin,numD)
    return denom

def make_num(tree,numSel):
    num = ROOT.TH1F("num","",nVar,0,nVar)
    numN = getNentries(tree,numSel)
    for i in range(num.GetNbinsX()):
        ibin = i+1
        num.SetBinContent(ibin,numN)
    return num

def produce_Nm1(tree, numSel, denomSel,color):
    denom = make_denom(tree,  denomSel)
    num = make_num(tree, numSel)
    l1 = make_efficiency(denom,num)
    l1.SetMarkerColor(color)
    l1.SetLineColor(color)
    return l1

def save_plot(tree1,tree2,numSel,denomSel,labels,filename,legend):
    MC   = produce_Nm1(tree1,numSel,denomSel,ROOT.kRed)
    data = produce_Nm1(tree2,numSel,denomSel,ROOT.kBlack)
    frame = ROOT.TH1F("frame","N-1 Efficiency",nVar,0,nVar)
    frame.GetYaxis().SetRangeUser(0.75,1.05)
    for i in range(frame.GetNbinsX()):
        ibin = i+1
        frame.GetXaxis().SetBinLabel(ibin,labels[i]) 
    frame.GetXaxis().LabelsOption("v")
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    legend = ROOT.TLegend(0.5,0.1,0.8,0.4,"","brNDC")
    legend.SetFillStyle(0)
    legend.AddEntry(MC,"Drell-Yan MC","pe1")
    legend.AddEntry(data,legend,"pe1")
    legend.Draw()
    line = ROOT.TLine(0,1,13,1)
    line.SetLineStyle(7)
    line.Draw()
    MC.Draw("pe1")
    data.Draw("pe1same")
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    canvas.Write(filename)
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

xaxis='Cut'
yaxis='N-1 Efficiency'
ptCut =         "lepPt[0]>20 && lepPt[1]>20"
isTriggered =   "isTriggered==1" # 1 for MuMu only
singleMuTrig =  "singleMuTrig==1" # 1 for Mu only
globalMu =      "lepMuType[0]%2!=0          && lepMuType[1]%2!=0" # odd numbered if muon is global
PFMu =          "lepMuType[0]>8             && lepMuType[1]>8" # greater than 8 if muon is PF
chi2ndof =      "lepR9orChi2ndof[0]<10      && lepR9orChi2ndof[1]<10"
globalMuHits =  "globalMuonHits[0]>0        && globalMuonHits[1]>0"
matchStations = "matchedMuStations[0]>1     && matchedMuStations[0]>1"
dxy =           "fabs(lepMuonDxyPV[0])<0.2  && fabs(lepMuonDxyPV[1])<0.2"
dz =            "fabs(lepMuonDzPV[0])<0.5   && fabs(lepMuonDzPV[1])<0.5"
nPixelHits =    "lepMuonPixelHits[0]>0      && lepMuonPixelHits[1]>0"
trackerLayers = "lepMuonTrackerLayers[0]>5  && lepMuonTrackerLayers[1]>5"
isoCut =        "lepPFIsoDBCor[0]<0.12      && lepPFIsoDBCor[1]<0.12"
isTight =       "lepId[0]==1                && lepId[1]==1"
Sel =           "fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4 && llchid==-4 && nLeptons>1 && llM>60 && llM<120"
numSel1 = Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+isTight
denomSel1 = [
Sel+"&&"           +isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"            +singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"                 +globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"              +PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"          +chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"              +globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"                  +matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"                   +dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"         +dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"        +nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"                +trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig
]
labels1 = ["P_{T}","isoCut","HLT Mu Fired","GlobalMuon","PFMuon","#chi^{2}/ndof<10","N(Muon Hits)>0","N(matched stations)>1","|dxy|<0.2","|dz|<0.5","N(pixel hits)>0","N(tracker layers)>5","isTight"]
numSel2 = Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+isTight
denomSel2 = [
Sel+"&&"           +isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"            +isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"                 +globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"              +PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"          +chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"              +globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"                  +matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"                   +dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"         +dz+"&&"+nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"        +nPixelHits+"&&"+trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"                +trackerLayers,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered
]
labels2 = ["P_{T}","isoCut","HLT MuMu Fired","GlobalMuon","PFMuon","#chi^{2}/ndof<10","N(Muon Hits)>0","N(matched stations)>1","|dxy|<0.2","|dz|<0.5","N(pixel hits)>0","N(tracker layers)>5","isTight"]
nVar = len(denomSel1)
filename = "Nm1Eff"

save_plot(DYtree,singleMuTree,numSel1,denomSel1,labels1,"Nm1EffSingleMu","Single Mu Data")
save_plot(DYtree,doubleMuTree,numSel2,denomSel2,labels2,"Nm1EffDoubleMu","Double Mu Data")
