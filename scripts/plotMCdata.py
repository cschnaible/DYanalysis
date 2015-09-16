'''
Script to plot some simple data/MC comparisions
    Data & MC Single/Double Muon
    Data / MC Single/Double Muon
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


#####################################
#Get Effi NTUPLE                 #
#####################################

DY_file = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_50ns/DYJetsToLL_50ns.root")
DYtree = DY_file.Get("analyzer/events")
singleMuFile = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/SingleMuon/Run2015B/SingleMuon_Run2015B.root")
singleMu = singleMuFile.Get("analyzer/events")
doubleMuFile = ROOT.TFile("/afs/cern.ch/work/c/cschnaib/DYntuples/DoubleMuon/Run2015B/DoubleMuon_Run2015B.root")
doubleMu = doubleMuFile.Get("analyzer/events")


######## File #########

outputName = '~/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/DYdataComparision_noScale.root'
#outputName = '~/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/DYdataComparision_Scale.root'
outputRoot = ROOT.TFile(outputName,"recreate")

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/CMSSW_7_4_6_patch6/src/DYanalysis/plots/'+argv[2]+'_'
else:
   saveWhere='~/CMSSW_7_4_6_patch6/src/DYanalysis/plots/'


DYxs = 6025.2
Lumi = 41.856481082


canvas = ROOT.TCanvas()



def make_histW(tree,variable,weight,selection,binning,color):
    ''' Plot a distribution for a specific sample '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, weight+"*"+selection, "goff")
    hist = ROOT.gDirectory.Get("htemp").Clone()
    nInt = hist.Integral()
    return hist,nInt
def make_data(tree,variable,selection,binning):
    ''' Plot a quantity using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string,selection,"goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    nInt=output_histo.Integral()
    return output_histo,nInt

def make_ratio(data,MC):
    ''' Make a ratio plot '''
    ratio = data.Clone()
    ratio.Sumw2()
    #ratio.SetMinimum(0.5)
    #ratio.SetMaximum(3)
    ratio.Divide(MC)
    return ratio

def getWsum(file_name):
    Wsum = file_name.Get("analyzer/WEvents").GetBinContent(5) # Sum of weights for every event
    return Wsum#,Nevt

def make_compP(variable,dataTree,legendData,MCtree,legendMC,selection,binning,filename,xaxis,yaxis):
    data,dInt = make_data(dataTree,variable,selection,binning)
    data.SetLineColor(ROOT.kBlack)
    data.SetLineWidth(2)
    DYw = getWsum(DY_file)
    DY,DYint = make_histW(MCtree,variable, "(mcWeight*%s*%s/%s)"%(DYxs,Lumi,DYw)  ,selection,binning,ROOT.kBlue)
    #DY,DYint = make_data(MCtree,variable,selection,binning)
    #scale = dInt/DYint
    #DY.Scale(scale)
    DY.SetLineColor(ROOT.kBlue)
    DY.SetLineWidth(2)
    DY.GetXaxis().SetTitle(xaxis)
    DY.GetYaxis().SetTitle(yaxis)
    DY.Draw("hist")
    data.Draw("pe1same")
    legend = ROOT.TLegend(0.65, 0.7, 0.9, 0.9, "", "brNDC")
    legend.SetFillStyle(0)
    legend.AddEntry(DY,legendMC,"l")
    legend.AddEntry(data,legendData,"pe1")
    legend.Draw()
    canvas.Write(filename)
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def make_compR(variable,dataTree,legendData,MCtree,legendMC,selection,binning,filename,xaxis,yaxis):
    data,dInt = make_data(dataTree,variable,selection,binning)
    DYw = getWsum(DY_file)
    DY,DYint = make_histW(MCtree,variable, "(mcWeight*%s*%s/%s)"%(DYxs,Lumi,DYw)  ,selection,binning,ROOT.kBlue)
    #DY,DYint = make_data(MCtree,variable, selection,binning)
    #scale = dInt/DYint
    #DY.Scale(scale)
    ratio = make_ratio(data,DY)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(2)
    ratio.Draw("pe1")
    ratio.GetXaxis().SetTitle(xaxis)
    ratio.GetYaxis().SetTitle(yaxis)
    line = ROOT.TLine(binning[1],1,binning[2],1)
    line.SetLineStyle(7)
    line.Draw()
    canvas.Write(filename)
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def make_dataP(variable,dataTree,legendData,MCtree,legendMC,selection,binning,filename,xaxis,yaxis):
    data1,data1int = make_data(dataTree,variable,selection,binning)
    data1.SetLineColor(ROOT.kBlack)
    data1.SetLineWidth(2)
    data2,data2int = make_data(MCtree,variable,selection,binning)
    data2.SetLineColor(ROOT.kBlue)
    data2.SetLineWidth(2)
    data1.Draw()
    data2.Draw("same")
    data1.GetXaxis().SetTitle(xaxis)
    data1.GetYaxis().SetTitle(yaxis)
    legend = ROOT.TLegend(0.65, 0.7, 0.9, 0.9, "", "brNDC")
    legend.SetFillStyle(0)
    legend.AddEntry(data1,legendData,"l")
    legend.AddEntry(data2,legendMC,"l")
    legend.Draw()
    canvas.Write(filename)
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def make_dataR(variable,dataTree,legendData,MCtree,legendMC,selection,binning,filename,xaxis,yaxis):
    data1,int1 = make_data(dataTree,variable,selection,binning)
    data2,int2 = make_data(MCtree,variable,selection,binning)
    ratio = make_ratio(data1,data2)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(2)
    ratio.Draw("pe1")
    #ratio.SetMinimum(0)
    #ratio.SetMaximum(3)
    ratio.GetXaxis().SetTitle(xaxis)
    ratio.GetYaxis().SetTitle(yaxis)
    line = ROOT.TLine(binning[1],1,binning[2],1)
    line.SetLineStyle(7)
    line.Draw()
    canvas.Write(filename)
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

selection20single = "(llM>60 && llM<120 && lepId[0]==1 && lepId[1]==1 && lepPt[0]>20 && lepPt[1]>20 && fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4 && singleMuTrig==1 && lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15)"
selection48single = "(llM>60 && llM<120 && lepId[0]==1 && lepId[1]==1 && lepPt[0]>48 && lepPt[1]>48 && fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4 && singleMuTrig==1 && lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15)"
selection20double = "(llM>60 && llM<120 && lepId[0]==1 && lepId[1]==1 && lepPt[0]>20 && lepPt[1]>20 && fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4 && isTriggered==1 && lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15)"
selection48double = "(llM>60 && llM<120 && lepId[0]==1 && lepId[1]==1 && lepPt[0]>48 && lepPt[1]>48 && fabs(lepEta[0])<2.4 && fabs(lepEta[1])<2.4 && isTriggered==1 && lepPFIsoDBCor[0]<0.15 && lepPFIsoDBCor[1]<0.15)"

binning = [30,60,120]

make_compP("llM",singleMu,"Single Muon Data",DYtree,"Drell-Yan MC",selection20single,binning,"SingleMuonMCplot20noScale","M(#mu#mu)[GeV]","Events")
make_compR("llM",singleMu,"Single Muon Data",DYtree,"Drell-Yan MC",selection20single,binning,"SingleMuonMCratio20noScale","M(#mu#mu)[GeV]","Ratio")
make_compP("llM",singleMu,"Single Muon Data",DYtree,"Drell-Yan MC",selection48single,binning,"SingleMuonMCplot48noScale","M(#mu#mu)[GeV]","Events")
make_compR("llM",singleMu,"Single Muon Data",DYtree,"Drell-Yan MC",selection48single,binning,"SingleMuonMCratio48noScale","M(#mu#mu)[GeV]","Ratio")
    
make_compP("llM",doubleMu,"Double Muon Data",DYtree,"Drell-Yan MC",selection20double,binning,"DoubleMuonMCplot20noScale","M(#mu#mu)[GeV]","Events")
make_compR("llM",doubleMu,"Double Muon Data",DYtree,"Drell-Yan MC",selection20double,binning,"DoubleMuonMCratio20noScale","M(#mu#mu)[GeV]","Data/MC")
make_compP("llM",doubleMu,"Double Muon Data",DYtree,"Drell-Yan MC",selection48double,binning,"DoubleMuonMCplot48noScale","M(#mu#mu)[GeV]","Events")
make_compR("llM",doubleMu,"Double Muon Data",DYtree,"Drell-Yan MC",selection48double,binning,"DoubleMuonMCratio48noScale","M(#mu#mu)[GeV]","Data/MC")

make_dataP("llM",singleMu,"Single Muon Data",doubleMu,"Double Muon Data",selection20double,binning,"SingleDoubleMuonplot20noScale","M(#mu#mu)[GeV]","Events")
make_dataR("llM",singleMu,"Single Muon Data",doubleMu,"Double Muon Data",selection20double,binning,"SingleDoubleMuonratio20noScale","M(#mu#mu)[GeV]","Single/Double")
make_dataP("llM",singleMu,"Single Muon Data",doubleMu,"Double Muon Data",selection48double,binning,"SingleDoubleMuonplot48noScale","M(#mu#mu)[GeV]","Events")
make_dataR("llM",singleMu,"Single Muon Data",doubleMu,"Double Muon Data",selection48double,binning,"SingleDoubleMuonratio48noScale","M(#mu#mu)[GeV]","Single/Double")

outputRoot.Close()
