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
outputName = '~/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/AccPlots.root'
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

canvas = ROOT.TCanvas()

def make_plot(tree, variable, selection, binning, filename, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    if filename in ("leadMuonAccEta_lowBinning","subleadMuonAccEta_lowBinning"):
       output_histo = output_histo.Rebin(4,"htemp",EtaBinning2) 
    if filename == "zAccM_rebin":
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
    acc = N/D
    return l1,N,D,acc

def produce_acceptance(ntuple, variable, numSel, denomSel, binning, filename, title="", xaxis="", yaxis=""):
    if filename in ("leadMuonAccEta_lowBinning","subleadMuonAccEta_lowBinning"):
        frame = ROOT.TH1F("frame","frame",4,EtaBinning2)
    elif filename == "zAccM_rebin":
        frame = ROOT.TH1F("frame","frame",13,MBinning2)
    else:
        frame = ROOT.TH1F("frame","frame",*binning)
    a1,N,D,acc = produce_efficiency(ntuple, variable, numSel, denomSel, binning, filename, ROOT.kBlack )
    frame.SetMaximum(1.05)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    legend = ROOT.TLegend(0.5,0.2,0.7,0.3,"","brNDC")
    legend.SetFillStyle(0)
    legend.AddEntry(a1,"#splitline{Drell-Yan MC}{#alpha = %s}"%(acc),"pe1")
    legend.Draw()
    a1.Draw('pe')
    a1.SetName(filename)
    canvas.Write(filename)
    #outputRoot.Write()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)




################################################################################

################################################################################

#numSel = "(lepPtGEN>20 && fabs(lepEtaGEN)<2.4 && llMGEN>60 && llMGEN<120)"
#denomSel = "(llMGEN>60 && llMGEN<120)"
ptBinning = [100,0,500]
muoPtBinning = [50,0,250]
EtaBinning1 = [48,-2.4,2.4]
EtaBinning2 = array('d',[-2.4,-1.5,0,1.5,2.4])
#MBinning = [15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200]
MBinning1 = [60,60,120]
MBinning2 = array('d',[60, 64, 68, 72, 76, 81, 86, 91,96, 101, 106, 110, 115, 120])
YBinning = [60,-3,3]
PhiBinning = [50,-3.14159,3.14159]
DPhiBinning = [25,0,3.14159]
llnumSel = "fabs(lepEtaGEN[0])<2.4 && fabs(lepEtaGEN[1])<2.4 && lepPtGEN[0]>20 && lepPtGEN[1]>20 && lepPdgIdGEN[0]*lepPdgIdGEN[1]==-169 && nGenLeptons>1 && llMGEN>60 && llMGEN<120"
lldenomSel = "lepPdgIdGEN[0]*lepPdgIdGEN[1]==-169 && nGenLeptons>1 && llMGEN>60 && llMGEN<120"
numPtSel = "lepPdgIdGEN[0]*lepPdgIdGEN[1]==-169 && nGenLeptons>1 && llMGEN>60 && llMGEN<120 && fabs(lepEtaGEN[0])<2.4 && fabs(lepEtaGEN[1])<2.4"
numEtaSel = "lepPtGEN[0]>20 && lepPtGEN[1]>20 && lepPdgIdGEN[0]*lepPdgIdGEN[1]==-169 && nGenLeptons>1 && llMGEN>60 && llMGEN<120"

produce_acceptance(DYtree, "llMGEN", llnumSel,lldenomSel, 
                      MBinning1, #binning
                     'zAccM', #filename
                     'Z candidate acceptance vs M', #title
                     'gen M(#mu#mu) [GeV]', #xaxis
                     'A(#mu#mu)'
)
print " "
produce_acceptance(DYtree, "llMGEN", llnumSel,lldenomSel,
                      MBinning1, #binning
                     'zAccM_rebin', #filename
                     'Z candidate acceptance vs M', #title
                     'gen M(#mu#mu) [GeV]', #xaxis
                     'A(Z)'
)
print " "
produce_acceptance(DYtree, "llYGEN", llnumSel,lldenomSel,
                     YBinning, #binning
                     'zAccY', #filename
                     'Z candidate acceptance vs Y', #title,
                     'Y(#mu#mu)', #xaxis
                     'A(Z)'
)
print " "
produce_acceptance(DYtree, "llPtGEN", llnumSel,lldenomSel,
                      ptBinning, #binning
                     'zAccPt', #filename
                     'Z candidate acceptance vs DiMuon p_{T}', #title
                     'gen #Delta#phi(#mu#mu)', #xaxis
                     'A(Z)'
)

print " "
produce_acceptance(DYtree, "llDPhiGEN", llnumSel,lldenomSel,
                      DPhiBinning, #binning
                     'zAccDPhi', #filename
                     'Z candidate acceptance vs muon #Delta#phi', #title
                     'gen #Delta#phi(#mu#mu)', #xaxis
                     'A(Z)'
)

print " "
produce_acceptance(DYtree, "lepPtGEN[0]", numPtSel,lldenomSel,
                     muoPtBinning, #binning
                     'leadMuonAccPt', #filename
                     'Leading Muon Acceptance for |#eta|<2.4 vs muon p_{T}', #title
                     'gen p_{T}(#mu_{1}) [GeV]', #xaxis
                     'A(#mu)'
)
print " "
produce_acceptance(DYtree, "lepPtGEN[1]", numPtSel,lldenomSel,
                     muoPtBinning, #binning
                     'subleadMuonAccPt', #filename
                     'Subleading Muon Acceptance for |#eta|<2.4 vs muon p_{T}', #title
                     'gen p_{T}(#mu_{2}) [GeV]', #xaxis
                     'A(#mu)'
)
print " "
produce_acceptance(DYtree, "lepEtaGEN[0]", numEtaSel,lldenomSel,
                      EtaBinning1, #binning
                     'leadMuonAccEta_lowBinning', #filename
                     'Leading Muon Acceptance for p_{T}>20GeV vs #eta', #title
                     'gen #eta(#mu_{1})', #xaxis
                     'A(#mu)'
)
print " "
produce_acceptance(DYtree, "lepEtaGEN[1]", numEtaSel,lldenomSel,
                      EtaBinning1, #binning
                     'subleadMuonAccEta_lowBinning', #filename
                     'Subleading Muon Acceptance for p_{T}>20GeV vs #eta', #title
                     'gen #eta(#mu_{2})', #xaxis
                     'A(#mu)'
)
print " "
produce_acceptance(DYtree, "lepEtaGEN[0]", numEtaSel,lldenomSel,
                      EtaBinning1, #binning
                     'leadMuonAccEta_highBinning', #filename
                     'Leading Muon Acceptance for p_{T}>20GeV vs #eta', #title
                     'gen #mu #eta(#mu_{1})', #xaxis
                     'A(#mu)'
)
print " "
produce_acceptance(DYtree, "lepEtaGEN[1]", numEtaSel,lldenomSel,
                      EtaBinning1, #binning
                     'subleadMuonAccEta_highBinning', #filename
                     'Subleading Muon Acceptance for p_{T}>20GeV vs #eta', #title
                     'gen #eta(#mu_{2})', #xaxis
                     'A(#mu)'
)
print " "
produce_acceptance(DYtree, "lepPhiGEN[0]",llnumSel,lldenomSel,
                      PhiBinning, #binning
                     'leadMuonAccPhi', #filename
                     'Leading Muon Acceptance for p_{T}>20GeV && |#eta|<2.4 vs #phi', #title
                     'gen #phi(#mu_{2})', #xaxis
                     'A(#mu)'
)
print " "
produce_acceptance(DYtree, "lepPhiGEN[1]",llnumSel,lldenomSel,
                      PhiBinning, #binning
                     'subleadMuonAccPhi', #filename
                     'Subleading Muon Acceptance for p_{T}>20GeV && |#eta|<2.4 vs #phi', #title
                     'gen #phi(#mu_{2})', #xaxis
                     'A(#mu)'
)

outputRoot.Close()
