import ROOT as rt
from array import array

outFile = rt.TFile.Open("tmp.root","RECREATE")
dataFile = rt.TFile.Open("outData_ABCD_jetVeto10.root","READ")
sigFile = rt.TFile.Open("outSig_ABCD_jetVeto10_test.root","READ")
dataHistos = {}
sigHistos = {}
years = ["2016","2017","2018"]
mX = ['15','40','55']
ctau = '10000'
colors = [1,2,3,4,6,7,38,41]
styles = [1,7]
vars = ["dtRechitClusterNSegmentStation2","dtRechitClusterNSegmentStation3","dtRechitClusterNSegmentStation4","dtRechitClusterMaxStation","dtRechitClusterMaxStationRatio","dtRechitClusterNStation","dtRechitClusterMaxChamber","dtRechitClusterNChamber","dtRechitClusterMaxChamberRatio","dtRechitClusterX","dtRechitClusterY","dtRechitClusterZ","dtRechitClusterEta","dtRechitClusterPhi","dtRechitClusterXSpread","dtRechitClusterYSpread","dtRechitClusterZSpread","dtRechitClusterEtaSpread","dtRechitClusterPhiSpread","dtRechitClusterTime","dtRechitClusterTimeSpread","dtRechitClusterMajorAxis","dtRechitClusterMinorAxis"]
vars = ["dtRechitClustersDR","dtRechitClustersVetoDR"]
vars = ["nMB1MatchAdjacent0p8","nMB1MatchAdjacent0p8Pi2"]
vars = ["leadJetPt","leadJetPtMET"]
vars = ["muonOverlap10GeVLoose","muonOverlap10GeVTight"]
vars = ['nDtSegsStation_50hits','nDtSegsStation_100hits','nDtSegsStation_150hits','nDtSegsWheel_50hits','nDtSegsWheel_100hits','nDtSegsWheel_150hits']
vars = ['nStations1Seg','nWheels1Seg','nStations5Seg','nWheels5Seg','nStations10Seg','nWheels10Seg']
vars = ['minSegmentDR','nMatchedSegments','nMatchedSegmentsOtherStations','nAlignedSegments','nMatchedSegmentsInnerStation','nMatchedSegmentsOuterStation']
#vars = ['nMatchedSegmentsInnerStation','nMatchedSegmentsOuterStation','MB1VetoOutcomes']
#vars = ['segmentAlignmentDeltaPhi','segmentAlignmentDeltaEta']
#vars = ['nMatchedSegmentsClusterStation','nMatchedSegmentsInnerAndAdjacentStation']
vars = ['MB1VetoOutcomesPassRecoMuon','MB1VetoOutcomesPassLooseMuon','MB1VetoOutcomesPassSegMuon','MB1VetoOutcomesPassOneSegMuon','AdjacentMB1VetoOutcomes','nMatchedSegmentsInnerStationPassRecoMuon','nMatchedSegmentsInnerStationPassLooseMuon','nMatchedSegmentsInnerStationPassSegMuon']
#vars = ['matchedSegmentTimeMeanClusterStation','matchedSegmentTimeMeanInnerStation','matchedSegmentTimeMeanOuterStation']
vars = ['muonLooseIDVetoPt']
vetoes = ['MB1CR','MB2withMB1CR','MB1HitsCR','MB2CR','SR']
vetoes = ['50hits_clusterMETCR','100hits_clusterMETCR','150hits_clusterMETCR']
vetoes = ['MB2_invertedShowerVetoes','MB3_invertedShowerVetoes','MB4_invertedShowerVetoes']
#vetoes = ['invertedShowerVetoes']
#vetoes = ['invertedJetVeto']
vetoes = ['MB2_invertedShowerVetoes','MB3_invertedShowerVetoes','MB4_invertedShowerVetoes','MB2_invertedJetVeto','MB3_invertedJetVeto','MB4_invertedJetVeto']
vetoes = ['Nminus1_clusterMETCR']
decay = ''
sigEff = {'noVeto':0.120958, 'muonVeto':0.0846705, 'jetVeto':0.0888669, 'MB1Veto':0.0540607, 'fullVeto':0.0409775}
bkgEff = {'noVeto':604688.0/9808483.0, 'muonVeto':90527.0/9808483.0, 'jetVeto':42786.0/9808483.0, 'MB1Veto':8477.0/9808483.0, 'fullVeto':1142.0/9808483.0}
sigEvents = {'15':5190000,'40':5135998,'55':5106227}
#sigEvents = {'15':1768000+1177000+2245000,'40':1699999+1173999+2262000,'55':1716228+1182999+2207000}

canvas = rt.TCanvas("canvas")
legend = rt.TLegend(0.55,0.62,0.85,0.85)
for veto in vetoes:
    nveto=0
    canvas.Clear()
    for v in vars:
        print(v+'_'+veto)
        nfile = 0
        canvas.cd()
        canvas.Clear()
        legend.Clear()
        rt.gStyle.SetOptStat(0)
        dataStack = rt.THStack("dataStack","")
        sigStack = rt.THStack("sigStack","")
        for y in years:
            dataHistos[y+v+veto] = dataFile.Get('h_'+v+'_'+veto+'_'+y)
            #dataHistos[y+v+veto] = dataFile.Get('h_'+v+'_'+y)
            dataHistos[y+v+veto].SetStats(0)
            dataHistos[y+v+veto].SetLineWidth(2)
            dataHistos[y+v+veto].SetLineColor(colors[nfile])
            dataHistos[y+v+veto].SetLineStyle(styles[nveto])
            #if(dataHistos[y+v+veto].Integral()>0.0):
            #    dataHistos[y+v+veto].Scale(1/dataHistos[y+v+veto].Integral())
            dataStack.Add(dataHistos[y+v+veto],"hist")
            if(y=='2018'):
                legend.AddEntry(dataHistos[y+v+veto],"Background","L")
        nfile+=1
        for m in mX:
            print('h_'+v+'_'+veto+'_'+decay+'_'+m+ctau)
            sigHistos[veto] = sigFile.Get('h_'+v+'_'+veto+'_'+m+'_'+ctau)
            sigHistos[veto].SetStats(0)
            sigHistos[veto].SetLineWidth(2)
            sigHistos[veto].SetLineColor(colors[nfile])
            if(sigHistos[veto].Integral()>0.0):
                sigHistos[veto].Scale(1/sigHistos[veto].Integral())
                #sigHistos[veto].Scale(1/sigHistos[veto].GetBinContent(5))
            #if(m!='450'):
            #    sigHistos[veto].Scale(48580.0*137.0/sigEvents[m]*0.01)
            #else:
            #    sigHistos[veto].Scale(184.5*137.0/49984)
            sigStack.Add(sigHistos[veto],"hist")
            legend.AddEntry(sigHistos[veto],'mX='+m+'GeV c#tau='+ctau,'L')
            nfile+=1
        
        #outFile.cd()
        #canvas.Write(v)
        #canvas.SaveAs(v+"_DataAll.png")
        #canvas.SaveAs(v+".C")

        canvas.SetLogy(1)
        canvas.SetLogx(0)
        dataStack.GetStack().Last().Draw('hist')
        if(dataStack.GetStack().Last().Integral()>0):
            dataStack.GetStack().Last().Scale(1/dataStack.GetStack().Last().Integral())
            #dataStack.GetStack().Last().Scale(1/dataStack.GetStack().Last().GetBinContent(5))
        #dataStack.Draw('nostack')
        dataStack.GetStack().Last().GetXaxis().SetTitle(v)
        #dataStack.GetStack().Last().SetMinimum(0.004)
        #dataStack.GetXaxis().SetTitle(v+'_'+veto)
        canvas.Update()
        sigStack.Draw("nostacksame")
        #sigStack.GetXaxis().SetTitle(v+'_'+veto)
        legend.Draw()
        canvas.Update()
        canvas.WaitPrimitive()
        #outFile.cd()
        #canvas.Write(v)
        canvas.SaveAs("plots/"+v+veto+".png")

outFile.Close()
dataFile.Close()
sigFile.Close()
