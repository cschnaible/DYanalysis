'''
Usage:python plot.py RootFile.root label[optional]
Script to plot cut-by-cut efficiency plots
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
        #print labels[i], numD
        denom.SetBinContent(ibin,numD)
    return denom

def make_num(tree,numSel):
    num = ROOT.TH1F("num","",nVar,0,nVar)
    numN = getNentries(tree,numSel)
    #print "numerator ", numN
    for i in range(num.GetNbinsX()):
        ibin = i+1
        num.SetBinContent(ibin,numN)
    return num

def produce_Ceff(tree, numSel, denomSel,color):
    denom = make_denom(tree,  denomSel)
    num = make_num(tree, numSel)
    l1 = make_efficiency(denom,num)
    for i in range(nVar):
        ibin=i+1
        #print "Upper: ", l1.GetErrorYhigh(ibin), "Lower: ", l1.GetErrorYlow(ibin)
    l1.SetMarkerColor(color)
    l1.SetLineColor(color)
    return l1

def save_plot(MCtree,datatree,numSel,denomSel,labels,filename,leg):
    MC   = produce_Ceff(MCtree,  numSel,denomSel,ROOT.kRed)
    data = produce_Ceff(datatree,numSel,denomSel,ROOT.kBlack)
    frame = ROOT.TH1F("frame","Cumulative Efficiency",nVar,0,nVar)
    #frame.GetYaxis().SetRangeUser(0.5,1.05)
    for i in range(frame.GetNbinsX()):
        ibin = i+1
        frame.GetXaxis().SetBinLabel(ibin,labels[i]) 
    frame.GetXaxis().LabelsOption("v")
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    legend = ROOT.TLegend(0.5,0.1,0.8,0.4,"","brNDC")
    legend.SetFillStyle(0)
    legend.AddEntry(MC,"Drell-Yan MC","f")
    legend.AddEntry(data,"%s"%(leg),"pe1")
    legend.Draw()
    line = ROOT.TLine(0,1,13,1)
    line.SetLineStyle(7)
    line.Draw()
    MC.Draw("pe1")
    data.Draw("pe1same")
    data.SetMarkerStyle(8)
    data.SetMarkerSize(0.8)
    canvas.Write(filename)
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

xaxis='Cut'
yaxis='Cumulative Efficiency'
ptCut =         "lepPt[0]>20 && lepPt[1]>20"
isTriggered =   "isTriggered==1" # 1 for MuMu only
singleMuTrig =  "singleMuTrig==1" # 1 for Mu only
globalMu =      "lepMuType[0]%2!=0          && lepMuType[1]%2!=0" # odd numbered if muon is global
PFMu =          "lepMuType[0]>8             && lepMuType[1]>8" # greater than 8 if muon is PF
chi2ndof =      "lepR9orChi2ndof[0]<10      && lepR9orChi2ndof[1]<10"
globalMuHits =  "globalMuonHits[0]>0        && globalMuonHits[1]>0"
matchStations = "matchedMuStations[0]>1     && matchedMuStations[0]>1"
dxy =           "fabs(lepDxy[0])<0.2        && fabs(lepDxy[1])<0.2"
dz =            "fabs(lepDz[0])<0.5         && fabs(lepDz[1])<0.5"
nPixelHits =    "lepMuonPixelHits[0]>0      && lepMuonPixelHits[1]>0"
trackerLayers = "lepMuonTrackerLayers[0]>5  && lepMuonTrackerLayers[1]>5"
isoCut =        "lepPFIsoDBCor[0]<0.12      && lepPFIsoDBCor[1]<0.12"
Sel =           "fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4 && nLeptons>1 && llM>60 && llM<120 && llchid==-4"
numSel1 = Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers
denomSel1 = [
Sel,
Sel+"&&"+ptCut,
Sel+"&&"+ptCut+"&&"+isoCut,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+singleMuTrig+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers
]
labels1 = ["|#eta|<2.4, N(#mu)>1","P_{T}","isoCut","HLT Mu Fired","GlobalMuon","PFMuon","#chi^{2}/ndof<10","N(Muon Hits)>0","N(matched stations)>1","|dxy|<0.2","|dz|<0.5","N(pixel hits)>0","N(tracker layers)>5"]
numSel2 = Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers
denomSel2 = [
Sel,
Sel+"&&"+ptCut,
Sel+"&&"+ptCut+"&&"+isoCut,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits,
Sel+"&&"+ptCut+"&&"+isoCut+"&&"+isTriggered+"&&"+globalMu+"&&"+PFMu+"&&"+chi2ndof+"&&"+globalMuHits+"&&"+matchStations+"&&"+dxy+"&&"+dz+"&&"+nPixelHits+"&&"+trackerLayers
]
labels2 = ["|#eta|<2.4, N(#mu)>1","P_{T}","isoCut","HLT MuMu Fired","GlobalMuon","PFMuon","#chi^{2}/ndof<10","N(Muon Hits)>0","N(matched stations)>1","|dxy|<0.2","|dz|<0.5","N(pixel hits)>0","N(tracker layers)>5"]
nVar = len(denomSel1)
filename = "CumulativeEff_test"

save_plot(DYtree,singleMuTree,numSel1,denomSel1,labels1,"sequentialCutsSingleMu","Single Mu Data")
save_plot(DYtree,doubleMuTree,numSel2,denomSel2,labels2,"sequentialCutsDoubleMu","Double Mu Data")
