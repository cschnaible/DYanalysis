'''
Usage:python plot.py RootFile.root label[optional]
Script to make some quick efficiency plots to test ntuple generation.
Author: L. Dodd, UW Madison
Edited: C. Schnaible, UC Los Angeles
Edited to plot acceptance for DY analysis
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
    singleMuon_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/SingleMuon/Run2015B/SingleMuon_Run2015B.root")
    doubleMuon_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/Run2015B/DoubleMuon_Run2015B.root")
outputName = '~/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/DYIsoEff.root'
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
singleMuonTree = singleMuon_file.Get("analyzer/events")
doubleMuonTree = doubleMuon_file.Get("analyzer/events")

canvas = ROOT.TCanvas()

def make_plot(tree, variable, selection, binning, filename, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    if filename in ("leadMuonEtaIsoEffSingleMu_rebin","subleadMuonEtaIsoEffSingleMu_rebin","leadMuonEtaIsoEffDoubleMu_rebin","subleadMuonEtaIsoEffDoubleMu_rebin"):
       output_histo = output_histo.Rebin(4,"htemp",EtaBinning2) 
    if filename in ("zMIsoEffSingleMu_rebin","zMIsoEffDoubleMu_rebin"):
       output_histo = output_histo.Rebin(13,"htemp",MBinning2) 
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    C = output_histo.Integral()
    return output_histo,C

def make_efficiency(denom, num):
    ''' Make an efficiency/acceptance graph with the style '''
    eff = ROOT.TGraphAsymmErrors(num, denom)
    eff.SetMarkerStyle(20)
    eff.SetMarkerSize(1.0)
    return eff

def make_num(ntuple, variable, numSel,binning,filename):
    num,N = make_plot(
        ntuple, variable,numSel,
        binning,filename
    )
    return num,N
def make_denom(ntuple, variable, denomSel,binning,filename):
    denom,D = make_plot(
        ntuple, variable,denomSel,
        binning,filename
    )
    return denom,D

def produce_efficiency(ntuple, variable, numSel, denomSel, binning, filename,color):
    denom,D = make_denom(ntuple, variable, denomSel, binning,filename)
    num,N = make_num(ntuple, variable, numSel, binning,filename)
    l1 = make_efficiency(denom,num)
    l1.SetMarkerColor(color)
    l1.SetLineColor(color)
    eff = N/D
    return l1,eff

def produce_acceptance(MCtree,datatree, variable, numSel, denomSel, binning, filename, title, xaxis, yaxis,leg):
    if filename in ("leadMuonEtaIsoEffSingleMu_rebin","subleadMuonEtaIsoEffSingleMu_rebin","leadMuonEtaIsoEffDoubleMu_rebin","subleadMuonEtaIsoEffDoubleMu_rebin"):
        frame = ROOT.TH1F("frame","frame",4,EtaBinning2)
    if filename in ("zMIsoEffSingleMu_rebin","zMIsoEffDoubleMu_rebin"):
        frame = ROOT.TH1F("frame","frame",13,MBinning2)
    else:
        frame = ROOT.TH1F("frame","frame",*binning)
    MC,mcEff = produce_efficiency(MCtree, variable, numSel, denomSel, binning, filename, ROOT.kRed )
    data,dataEff = produce_efficiency(datatree, variable, numSel, denomSel, binning, filename, ROOT.kBlack )
    frame.SetMaximum(1.05)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    line = ROOT.TLine(binning[1],1,binning[2],1)
    line.SetLineStyle(7)
    line.Draw()
    legend = ROOT.TLegend(0.5,0.1,0.8,0.4,"","brNDC")
    legend.SetFillStyle(0)
    legend.AddEntry(MC,"#splitline{Drell-Yan MC}{#varepsilon_{DY} = %s}"%(mcEff),"pe1")
    legend.AddEntry(data,"#splitline{%s}{#varepsilon_{data} = %s}"%(leg,dataEff),"pe1")
    legend.Draw()
    MC.Draw('pe1')
    data.Draw("pe1same")
    canvas.Write(filename)
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)




################################################################################

################################################################################

PtBinning = [100,0,500]
muoPtBinning = [50,0,250]
EtaBinning1 = [48,-2.4,2.4]
EtaBinning2 = array('d',[-2.4,-1.5,0,1.5,2.4])
#MBinning = [15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200]
MBinning1 = [60,60,120]
MBinning2 = array('d',[60, 64, 68, 72, 76, 81, 86, 91,96, 101, 106, 110, 115, 120])
YBinning = [26,-2.6,2.6]
PhiBinning = [50,-3.14159,3.14159]
DPhiBinning = [25,0,3.14159]

ptSel = "lepPt[0]>20 && lepPt[1]>20"
etaSel = "fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4"
muoSel = "llchid==-4 && nLeptons>1 && llM>60 && llM<120"
isoSel = "lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15"

produce_acceptance(DYtree, singleMuonTree, "llM", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      MBinning1, #binning
                     'zMIsoEffSingleMu', #filename
                     'DiMuon Candidate Isolation Efficiency', #title
                     'M(#mu#mu) [GeV]', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "llM", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      MBinning1, #binning
                     'zMIsoEffSingleMu_rebin', #filename
                     'DiMuon Candidiate Isolation Efficiency', #title
                     'M(#mu#mu) [GeV]', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "llY", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      YBinning, #binning
                     'zYIsoEffSingleMu', #filename
                     'DiMuon Candidiate Isolation Efficiency', #title
                     'Y(#mu#mu)', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "llDPhi", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      DPhiBinning, #binning
                     'zDPhiIsoEffSingleMu', #filename
                     'DiMuon Candidiate Isolation Efficiency', #title
                     '#Delta#phi(#mu#mu)', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "llPt", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      PtBinning, #binning
                     'zPtIsoEffSingleMu', #filename
                     'DiMuon Candidiate Isolation Efficiency', #title
                     'p_{T}(#mu#mu)', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "lepPt[0]", etaSel+"&&"+muoSel+"&&"+isoSel, etaSel+"&&"+muoSel, 
                      muoPtBinning, #binning
                     'leadMuonPtIsoEffSingleMu', #filename
                     'Leading Muon Isolation Efficiency', #title
                     'p_{T}(#mu_{1})', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "lepPt[1]", etaSel+"&&"+muoSel+"&&"+isoSel, etaSel+"&&"+muoSel, 
                      muoPtBinning, #binning
                     'subleadMuonPtIsoEffSingleMu', #filename
                     'Subleading Muon Isolation Efficiency', #title
                     'p_{T}(#mu_{2})', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "lepEta[0]", ptSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+muoSel, 
                      EtaBinning1, #binning
                     'leadMuonEtaIsoEffSingleMu', #filename
                     'Leading Muon Isolation Efficiency', #title
                     '#eta(#mu_{1})', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "lepEta[1]", ptSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+muoSel, 
                      EtaBinning1, #binning
                     'subleadMuonEtaIsoEffSingleMu', #filename
                     'Subleading Muon Isolation Efficiency', #title
                     '#eta(#mu_{2})', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "lepEta[0]", ptSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+muoSel, 
                      EtaBinning1, #binning
                     'leadMuonEtaIsoEffSingleMu_rebin', #filename
                     'Leading Muon Isolation Efficiency', #title
                     '#eta(#mu_{1})', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "lepEta[1]", ptSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+muoSel, 
                      EtaBinning1, #binning
                     'subleadMuonEtaIsoEffSingleMu_rebin', #filename
                     'Subleading Muon Isolation Efficiency', #title
                     '#eta(#mu_{2})', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "lepPhi[0]", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      PhiBinning, #binning
                     'leadMuonPhiIsoEffSingleMu', #filename
                     'Leading Muon Isolation Efficiency', #title
                     '#phi(#mu_{1})', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, singleMuonTree, "lepPhi[1]", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      PhiBinning, #binning
                     'subleadMuonPhiIsoEffSingleMu', #filename
                     'Subleading Muon Isolation Efficiency', #title
                     '#phi(#mu_{2})', #xaxis
                     'Efficiency',
                     'Single Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "llM", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      MBinning1, #binning
                     'zMIsoEffDoubleMu', #filename
                     'DiMuon Candidate Isolation Efficiency', #title
                     'M(#mu#mu) [GeV]', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "llM", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      MBinning1, #binning
                     'zMIsoEffDoubleMu_rebin', #filename
                     'DiMuon Candidiate Isolation Efficiency', #title
                     'M(#mu#mu) [GeV]', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "llY", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      YBinning, #binning
                     'zYIsoEffDoubleMu', #filename
                     'DiMuon Candidiate Isolation Efficiency', #title
                     'Y(#mu#mu)', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "llDPhi", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      DPhiBinning, #binning
                     'zDPhiIsoEffDoubleMu', #filename
                     'DiMuon Candidiate Isolation Efficiency', #title
                     '#Delta#phi(#mu#mu)', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "llPt", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      PtBinning, #binning
                     'zPtIsoEffDoubleMu', #filename
                     'DiMuon Candidiate Isolation Efficiency', #title
                     'p_{T}(#mu#mu)', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "lepPt[0]", etaSel+"&&"+muoSel+"&&"+isoSel, etaSel+"&&"+muoSel, 
                      muoPtBinning, #binning
                     'leadMuonPtIsoEffDoubleMu', #filename
                     'Leading Muon Isolation Efficiency', #title
                     'p_{T}(#mu_{1})', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "lepPt[1]", etaSel+"&&"+muoSel+"&&"+isoSel, etaSel+"&&"+muoSel, 
                      muoPtBinning, #binning
                     'subleadMuonPtIsoEffDoubleMu', #filename
                     'Subleading Muon Isolation Efficiency', #title
                     'p_{T}(#mu_{2})', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "lepEta[0]", ptSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+muoSel, 
                      EtaBinning1, #binning
                     'leadMuonEtaIsoEffDoubleMu', #filename
                     'Leading Muon Isolation Efficiency', #title
                     '#eta(#mu_{1})', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "lepEta[1]", ptSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+muoSel, 
                      EtaBinning1, #binning
                     'subleadMuonEtaIsoEffDoubleMu', #filename
                     'Subleading Muon Isolation Efficiency', #title
                     '#eta(#mu_{2})', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "lepEta[0]", ptSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+muoSel, 
                      EtaBinning1, #binning
                     'leadMuonEtaIsoEffDoubleMu_rebin', #filename
                     'Leading Muon Isolation Efficiency', #title
                     '#eta(#mu_{1})', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "lepEta[1]", ptSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+muoSel, 
                      EtaBinning1, #binning
                     'subleadMuonEtaIsoEffDoubleMu_rebin', #filename
                     'Subleading Muon Isolation Efficiency', #title
                     '#eta(#mu_{2})', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "lepPhi[0]", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      PhiBinning, #binning
                     'leadMuonPhiIsoEffDoubleMu', #filename
                     'Leading Muon Isolation Efficiency', #title
                     '#phi(#mu_{1})', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
produce_acceptance(DYtree, doubleMuonTree, "lepPhi[1]", ptSel+"&&"+etaSel+"&&"+muoSel+"&&"+isoSel, ptSel+"&&"+etaSel+"&&"+muoSel, 
                      PhiBinning, #binning
                     'subleadMuonPhiIsoEffDoubleMu', #filename
                     'Subleading Muon Isolation Efficiency', #title
                     '#phi(#mu_{2})', #xaxis
                     'Efficiency',
                     'Double Muon Data'
)
print " "
outputRoot.Close()
