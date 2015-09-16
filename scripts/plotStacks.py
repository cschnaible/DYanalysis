'''
Usage:python plot.py option (25ns or 50ns) 
Script to make stack plots of event variables
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
#ROOT.gStyle.SetTitleSize(0.01)
#ROOT.gStyle.SetTitleFont(41)
#ROOT.gStyle.SetTitleFontSize(0.05)

######## File #########
if len(argv) < 2:
   print 'Usage:python plot.py RootFile.root label[optional]'
   exit()

option = argv[1]

if option == "25ns":
    DY_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_25ns/")
    #TTa_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_25ns/")
    #TTb_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MC_25ns/")
    WW_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WW_TuneCUETP8M1_13TeV-pythia8/MC_25ns/")
    WZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WZ_TuneCUETP8M1_13TeV-pythia8/MC_25ns/")
    ZZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/ZZ_TuneCUETP8M1_13TeV-pythia8/MC_25ns/")
    data_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/")
if option == "50ns":
    DY_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_50ns/DYJetsToLL_50ns.root")
    TT_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MC_50ns/TT_50ns.root")
    #TTa_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_50ns/tt_nlo_50ns.root")
    #TTb_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MC_50ns/tt_MLM_50ns.root")
    WW_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WW_TuneCUETP8M1_13TeV-pythia8/MC_50ns/WW_50ns.root")
    WZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/WZ_TuneCUETP8M1_13TeV-pythia8/MC_50ns/WZ_50ns.root")
    ZZ_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/ZZ_TuneCUETP8M1_13TeV-pythia8/MC_50ns/ZZ_50ns.root")
    data_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/Run2015B/DoubleMuon_Run2015B.root")

outputName = '~/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/DYstackplots.root'
outputRoot = ROOT.TFile(outputName,"recreate")
#outputRoot = ROOT.TFile(outputName,"update")

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/CMSSW_7_4_6_patch6/src/DYanalysis/plots/'+argv[2]+'_'
else:
   saveWhere='~/CMSSW_7_4_6_patch6/src/DYanalysis/plots/'

# Get Tree
DYtree   = DY_file.Get("analyzer/events")
TTtree  = TT_file.Get("analyzer/events")
#TTatree  = TTa_file.Get("analyzer/events")
#TTbtree  = TTb_file.Get("analyzer/events")
WWtree   = WW_file.Get("analyzer/events")
WZtree   = WZ_file.Get("analyzer/events")
ZZtree   = ZZ_file.Get("analyzer/events")
datatree = data_file.Get("analyzer/events")
# Cross sections
DYxs   = 6025.2# = 2008.4*3
TTxs  = 831.76
#TTaxs  = 831.76
#TTbxs  = 831.76
WWxs   = 118.7
WZxs   = 66
ZZxs   = 15.4
# Luminosity
Lumi = 41.856481082

canvas = ROOT.TCanvas("asdf", "adsf")

def make_hist(tree,variable,weight,selection,binning,color):
    ''' Plot a distribution for a specific sample '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, weight+"*"+selection, "goff")
    hist = ROOT.gDirectory.Get("htemp").Clone()
    hist.SetFillColor(color)
    hist.SetLineColor(color)
    nInt = hist.Integral()
    #hist.Sumw2()
    return hist,nInt

def make_data(tree,variable,selection,binning):
    ''' Plot a distribution for data '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    data = ROOT.gDirectory.Get("htemp").Clone()
    #data.Sumw2()
    data.SetMarkerStyle(8)
    data.SetMarkerSize(0.8)
    data.SetLineWidth(2)
    data.SetMarkerColor(ROOT.kBlack)
    data.SetLineColor(ROOT.kBlack)
    return data

def getWsum(file_name):
    Wsum = file_name.Get("analyzer/WEvents").GetBinContent(5) # Sum of weights for every event
    #Nevt = file_name.Get("analyzer/WEvents").GetBinContent(6) # Number of events
    return Wsum#,Nevt
    
def make_stack(variable,selection,binning):
    stack = ROOT.THStack("stack","")
    DYw = getWsum(DY_file)   # weight = weight*sigma*lumi/sum(weights)
    TTw = getWsum(TT_file)   # weight = weight*sigma*lumi/sum(weights)
    #TTaw = getWsum(TTa_file) # weight = weight*sigma*lumi/sum(weights)
    #TTbw = getWsum(TTb_file) # weight = weight*sigma*lumi/sum(weights)
    WWw = getWsum(WW_file)   # weight = weight*sigma*lumi/sum(weights)
    WZw = getWsum(WZ_file)   # weight = weight*sigma*lumi/sum(weights)
    ZZw = getWsum(ZZ_file)   # weight = weight*sigma*lumi/sum(weights)
    DY,DYint    = make_hist(DYtree,variable, "(mcWeight*%s*%s/%s)"%(DYxs,Lumi,DYw)  ,selection,binning,ROOT.kBlue)
    TT,TTint    = make_hist(TTtree,variable,"(mcWeight*%s*%s/%s)"%(TTxs,Lumi,TTw)  ,selection,binning,ROOT.kGreen)
    #TTa,TTaint  = make_hist(TTatree,variable,"(mcWeight*%s*%s/%s)"%(TTaxs,Lumi,TTaw),selection,binning,ROOT.kGreen)
    #TTb,TTbint  = make_hist(TTbtree,variable,"(mcWeight*%s*%s/%s)"%(TTbxs,Lumi,TTbw),selection,binning,ROOT.kGreen)
    WW,WWint    = make_hist(WWtree,variable, "(mcWeight*%s*%s/%s)"%(WWxs,Lumi,WWw)  ,selection,binning,ROOT.kYellow)
    WZ,WZint    = make_hist(WZtree,variable, "(mcWeight*%s*%s/%s)"%(WZxs,Lumi,WZw)  ,selection,binning,ROOT.kYellow)
    ZZ,ZZint    = make_hist(ZZtree,variable, "(mcWeight*%s*%s/%s)"%(ZZxs,Lumi,ZZw)  ,selection,binning,ROOT.kYellow)
    stack.Add(WW)
    stack.Add(WZ)
    stack.Add(ZZ)
    stack.Add(TT)
    #stack.Add(TTa)
    #stack.Add(TTb)
    stack.Add(DY)
    legend = ROOT.TLegend(0.65, 0.7, 0.9, 0.9, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(DY,"Drell-Yan", "f")
    legend.AddEntry(TT,"t#bar{t}", "f")
    legend.AddEntry(WW,"DiBoson", "f")
    print "INT - DY : ",DYint," | TT : ",TTint," | DiBoson : ",WWint+WZint+ZZint
    #print "INT - DY : ",DYint," | TT : ",TTaint+TTbint," | DiBoson : ",WWint+WZint+ZZint
    Int = DYint + TTint + WWint + WZint + ZZint
    #Int = DYint + TTaint + TTbint + WWint + WZint + ZZint
    return stack,legend,Int

def make_plot(variable,selection,binning,filename,xaxis,yaxis,leg):
    stack,legend,Int = make_stack(variable,selection,binning)
    data = make_data(datatree,variable,selection,binning)
    legend.AddEntry(data,"%s"%(leg), "pe1")
    stack.Draw("hist")
    data.Draw("pe1same")
    D = data.Integral()
    print "MC    : ",Int
    print "Data  : ",D
    R = D/Int
    #err = R*math.sqrt((1/D)+(1/Int))
    print "Ratio : ",R#," +/- ",err
    legend.Draw()
    stack.GetXaxis().SetTitle(xaxis)
    stack.GetYaxis().SetTitle(yaxis)
    #stack.SetTitle("CMS Preliminary \\sqrt{s}=13TeV \\int\\mathscr{L}dt = 42.8 pb^{-1}")
    canvas.SetName(filename)
    canvas.Write()
    canvas.SetLogy()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)


################################################################################

################################################################################

PtBinning = [80,0,400]
muoPtBinning = [55,0,275]
EtaBinning1 = [46,-2.8,2.8]
EtaBinning2 = array('d',[-2.4,-1.5,0,1.5,2.4])
MBinning1 = [60,60,120]
MBinning2 = array('d',[60, 64, 68, 72, 76, 81, 86, 91,96, 101, 106, 110, 115, 120])
YBinning = [30,-3,3]
PhiBinning = [50,-3.14159,3.14159]
DPhiBinning = [25,0,3.14159]
chi2ndofBin = [50,0,25]
nTckBin = [20,0,20]
nPixBin = [11,0,11]
nMuoHitsBin = [27,0,54]
nSegMatchBin = [7,0,7]
DxyBin = [50,-1,1]
DzBin = [50,-1,1]

#isTight = "(fabs(lepMuonDxyPV[0])<0.2 && fabs(lepMuonDxyPV[1])<0.2 && matchedMuStations[0]>1 && matchedMuStations[1]>1 && fabs(lepMuonDzPV[0])<0.5 && fabs(lepMuonDzPV[1])<0.5 && lepMuonTrackerLayers[0]>5 && lepMuonTrackerLayers[1]>5 && lepMuonPixelHits[0]>0 && lepMuonPixelHits[1]>0 && lepR9orChi2ndof[0]<10 && lepR9orChi2ndof[1]<10 && globalMuonHits[0]>0 && globalMuonHits[1]>0 && lepMuType[0]%2!=0 && lepMuType[1]%2!=0 && lepMuType[0]>8 && lepMuType[1]>8)"
#isoSel = "(lepPFIsoDBCor[0]<0.12 && lepPFIsoDBCor[1]<0.12)"
#Muo = "(lepChId[0]*lepChId[1]==-4 && llM>60 && llM<120 && isTriggered==1)"
#pt2020 = "(lepPt[0]>20 && lepPt[1]>20)"
#pt2010 = "(lepPt[0]>20 && lepPt[1]>10)"
#etaBB = "(fabs(lepEta[0])<0.8 && fabs(lepEta[1])<0.8)"
#etaEE = "(fabs(lepEta[0])<2.4 && fabs(lepEta[0])>1.3 && fabs(lepEta[1])<2.4 && fabs(lepEta[1])>1.3)"
#etaMM = "(fabs(lepEta[0])<1.3 && fabs(lepEta[0])>0.8 && fabs(lepEta[1])<1.3 && fabs(lepEta[1])>0.8)"
#etaBE = "((fabs(lepEta[0])<0.8 && fabs(lepEta[1])>0.8 && fabs(lepEta[1]<2.4))||(fabs(lepEta[1])<0.8 && fabs(lepEta[0])>0.8 && fabs(lepEta[0]<2.4)))"

IdIsoTrigMu = "(lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15 && lepId[0]==1 && lepId[1]==1 && singleMuTrig==1 && llchid==-4 && llM>60 && llM<120 && nLeptons>1)"
IdIsoTrigMuMu = "(lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15 && lepId[0]==1 && lepId[1]==1 && isTriggered==1 && llchid==-4 && llM>60 && llM<120 && nLeptons>1)"
IsoTrigMu = "(lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15 && singleMuTrig==1 && llchid==-4 && llM>60 && llM<120 && nLeptons>1)"
IsoTrigMuMu = "(lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15 && isTriggered==1 && llchid==-4 && llM>60 && llM<120 && nLeptons>1)"
kinGeo = "(lepPt[0]>20 && lepPt[1]>20 && fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4)"
kin = "(lepPt[0]>20 && lepPt[1]>20)"
etaBB = "(fabs(lepEta[0])<0.8 && fabs(lepEta[1])<0.8)"
etaEE = "(fabs(lepEta[0])>0.8 && fabs(lepEta[1])>0.8 && fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4)"
#etaEB = "( (fabs(lepEta[0])<0.8 fabs(lepEta[1])>0.8 && fabs(lepEta[1])<2.4) || (fabs(lepEta[1])<0.8 fabs(lepEta[0])>0.8 && fabs(lepEta[0])<2.4) )"

selectionMu = "("+kinGeo+"&&"+IdIsoTrigMu+")"
selectionIDMu = "("+kinGeo+"&&"+IsoTrigMu+")"
selectionBBMu = "("+kin+"&&"+etaBB+"&&"+IdIsoTrigMu+")"
selectionEEMu = "("+kin+"&&"+etaEE+"&&"+IdIsoTrigMu+")"
selectionMuMu = "("+kinGeo+"&&"+IdIsoTrigMuMu+")"
selectionIDMuMu = "("+kinGeo+"&&"+IsoTrigMuMu+")"
selectionBBMuMu = "("+kin+"&&"+etaBB+"&&"+IdIsoTrigMuMu+")"
selectionEEMuMu = "("+kin+"&&"+etaEE+"&&"+IdIsoTrigMuMu+")"
#selectionEB = "("+kin+"&&"+etaEB+"&&"+IdIsoTrig+")"
#selectionBB2010 = "("+pt2010+"&&"+etaBB+"&&"+Muo+"&&"+isoSel+"&&"+isTight+")"
#selectionEE2010 = "("+pt2010+"&&"+etaEE+"&&"+Muo+"&&"+isoSel+"&&"+isTight+")"

make_plot("llM",                 selectionMu,MBinning1,     "llM_SingleMu",             "M(#mu#mu) [GeV]",     "Events/Gev","Single Muon Data")
make_plot("llY",                 selectionMu,YBinning,      "llY_SingleMu",               "Y(#mu#mu)",           "Events","Single Muon Data")
make_plot("llDPhi",              selectionMu,DPhiBinning,   "llDPhi_SingleMu",            "#Delta#phi(#mu#mu)",  "Events","Single Muon Data")
make_plot("llPt",                selectionMu,PtBinning,     "llPt_SingleMu",              "p_{T}(#mu#mu) [GeV]", "Events/5 GeV","Single Muon Data")
make_plot("llM",                 selectionBBMu,MBinning1,   "llM_BB_SingleMu",                  "M(#mu#mu) [GeV]",     "Events/GeV","Single Muon Data")
make_plot("llM",                 selectionEEMu,MBinning1,   "llM_EE_SingleMu",                  "M(#mu#mu) [GeV]",     "Events/GeV","Single Muon Data")
make_plot("lepPt[0]",            selectionMu,muoPtBinning,  "leadMuoPt_SingleMu",         "p_{T}(#mu_{1}) [GeV]","Events/5 Gev","Single Muon Data")
make_plot("lepPt[1]",            selectionMu,muoPtBinning,  "subleadMuoPt_SingleMu",      "p_{T}(#mu_{2}) [GeV]","Events/5 Gev","Single Muon Data")
make_plot("lepEta[0]",           selectionMu,EtaBinning1,   "leadMuoEta_SingleMu",        "#eta(#mu_{1})",       "Events","Single Muon Data")
make_plot("lepEta[1]",           selectionMu,EtaBinning1,   "subleadMuoEta_SingleMu",     "#eta(#mu_{2})",       "Events","Single Muon Data")
make_plot("lepR9orChi2ndof",     selectionIDMu,chi2ndofBin, "chi2ndof_SingleMu",          "#chi^{2}/ndof(#mu)",  "Events","Single Muon Data")
make_plot("lepMuonTrackerLayers",selectionIDMu,nTckBin,     "nTrackerLayers_SingleMu",    "N_{Tk. layers}",      "Events","Single Muon Data")
make_plot("lepMuonPixelHits",    selectionIDMu,nPixBin,     "nPixelHits_SingleMu",        "N_{pixel} Hits",      "Events","Single Muon Data")
make_plot("globalMuonHits",      selectionIDMu,nMuoHitsBin, "nGlobalMuonHits_SingleMu",   "N_{muon} Hits",       "Events","Single Muon Data")
make_plot("matchedMuStations",   selectionIDMu,nSegMatchBin,"nMatchedMuStations_SingleMu","N_{matched stations}","Events","Single Muon Data")
make_plot("lepMuonDxyPV",        selectionIDMu,DxyBin,      "Dxy_SingleMu",               "d_{xy}(PV) [cm]",     "Events","Single Muon Data")
make_plot("lepMuonDzPV",         selectionIDMu,DzBin,       "Dz_SingleMu",                "d_{z}(PV) [cm]",      "Events","Single Muon Data")

make_plot("llM",                 selectionMuMu,MBinning1,     "llM_DoubleMu",             "M(#mu#mu) [GeV]",     "Events/Gev","Double Muon Data")
make_plot("llY",                 selectionMuMu,YBinning,      "llY_DoubleMu",               "Y(#mu#mu)",           "Events","Double Muon Data")
make_plot("llDPhi",              selectionMuMu,DPhiBinning,   "llDPhi_DoubleMu",            "#Delta#phi(#mu#mu)",  "Events","Double Muon Data")
make_plot("llPt",                selectionMuMu,PtBinning,     "llPt_DoubleMu",              "p_{T}(#mu#mu) [GeV]", "Events/5 GeV","Double Muon Data")
make_plot("llM",                 selectionBBMuMu,MBinning1,   "llM_BB_DoubleMu",                  "M(#mu#mu) [GeV]",     "Events/GeV","Double Muon Data")
make_plot("llM",                 selectionEEMuMu,MBinning1,   "llM_EE_DoubleMu",                  "M(#mu#mu) [GeV]",     "Events/GeV","Double Muon Data")
make_plot("lepPt[0]",            selectionMuMu,muoPtBinning,  "leadMuoPt_DoubleMu",         "p_{T}(#mu_{1}) [GeV]","Events/5 Gev","Double Muon Data")
make_plot("lepPt[1]",            selectionMuMu,muoPtBinning,  "subleadMuoPt_DoubleMu",      "p_{T}(#mu_{2}) [GeV]","Events/5 Gev","Double Muon Data")
make_plot("lepEta[0]",           selectionMuMu,EtaBinning1,   "leadMuoEta_DoubleMu",        "#eta(#mu_{1})",       "Events","Double Muon Data")
make_plot("lepEta[1]",           selectionMuMu,EtaBinning1,   "subleadMuoEta_DoubleMu",     "#eta(#mu_{2})",       "Events","Double Muon Data")
make_plot("lepR9orChi2ndof",     selectionIDMuMu,chi2ndofBin, "chi2ndof_DoubleMu",          "#chi^{2}/ndof(#mu)",  "Events","Double Muon Data")
make_plot("lepMuonTrackerLayers",selectionIDMuMu,nTckBin,     "nTrackerLayers_DoubleMu",    "N_{Tk. layers}",      "Events","Double Muon Data")
make_plot("lepMuonPixelHits",    selectionIDMuMu,nPixBin,     "nPixelHits_DoubleMu",        "N_{pixel} Hits",      "Events","Double Muon Data")
make_plot("globalMuonHits",      selectionIDMuMu,nMuoHitsBin, "nGlobalMuonHits_DoubleMu",   "N_{muon} Hits",       "Events","Double Muon Data")
make_plot("matchedMuStations",   selectionIDMuMu,nSegMatchBin,"nMatchedMuStations_DoubleMu","N_{matched stations}","Events","Double Muon Data")
make_plot("lepMuonDxyPV",        selectionIDMuMu,DxyBin,      "Dxy_DoubleMu",               "d_{xy}(PV) [cm]",     "Events","Double Muon Data")
make_plot("lepMuonDzPV",         selectionIDMuMu,DzBin,       "Dz_DoubleMu",                "d_{z}(PV) [cm]",      "Events","Double Muon Data")

outputRoot.Close()
#make_plot("llM",selectionBB2020,MBinning1,"llM_BB2020","M(#mu#mu) [GeV]","Events/Gev")
#make_plot("llM",selectionBB2010,MBinning1,"llM_BB2010","M(#mu#mu) [GeV]","Events/Gev")
#make_plot("llM",selectionEE2010,MBinning1,"llM_EE2010","M(#mu#mu) [GeV]","Events/Gev")
