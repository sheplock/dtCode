import ROOT as rt
from array import array

outFile = rt.TFile.Open("tmp.root","RECREATE")
dataFile = rt.TFile.Open("outData_beamHalo.root","READ")
sigFile = rt.TFile.Open("outSig_central_V1p17_altTest4.root","READ")
dataHistos = {}
sigHistos = {}
years = ["2016","2017","2018"]
mX = ['15','40','55']
ctau = '10000'
colors = [1,2,3,4,6,7,38,41]
vars = ['nRpcRechits']
vetoes = ['fullVeto']
decay = ''

canvas = rt.TCanvas("canvas")
legend = rt.TLegend(0.55,0.62,0.85,0.85)
for v in vars:
    for veto in vetoes:
        print(v+'_'+veto)
        nfile = 0
        canvas.cd()
        canvas.Clear()
        legend.Clear()
        rt.gStyle.SetOptStat(0)
        dataStack = rt.THStack("dataStack","")
        sigStack = rt.THStack("sigStack","")
        for y in years:
            dataHistos[y+veto] = dataFile.Get('h_'+v+'_'+veto+'_'+y)
            #dataHistos[y+veto] = dataFile.Get('h_'+v+'_'+y)
            dataHistos[y+veto].SetStats(0)
            dataHistos[y+veto].SetLineColor(colors[nfile])
            dataStack.Add(dataHistos[y+veto],"hist")
            if(y=='2018'):
                legend.AddEntry(dataHistos[y+veto],"Background","L")
        nfile+=1
        for m in mX:
            sigHistos[veto] = sigFile.Get('h_'+v+decay+'_'+veto+'_'+m+'_'+ctau)
            sigHistos[veto].SetStats(0)
            sigHistos[veto].SetLineColor(colors[nfile])
            if(sigHistos[veto].Integral()>0.0):
                sigHistos[veto].Scale(1/sigHistos[veto].Integral())
            sigStack.Add(sigHistos[veto],"hist")
            legend.AddEntry(sigHistos[veto],'mX='+m+'GeV c#tau='+ctau,'L')
            nfile+=1
        
        canvas.SetLogy(0)
        canvas.SetLogx(0)
        dataStack.GetStack().Last().Draw('hist')
        if(dataStack.GetStack().Last().Integral()>0):
            dataStack.GetStack().Last().Scale(1/dataStack.GetStack().Last().Integral())
        #dataStack.Draw('nostack')
        dataStack.GetStack().Last().GetXaxis().SetTitle(v)
        dataStack.GetStack().Last().SetMinimum(0.004)
        canvas.Update()
        legend.Draw()
        canvas.Update()
        #canvas.WaitPrimitive()
        outFile.cd()
        canvas.Write(v)
        canvas.SaveAs("plots/"+v+"_"+veto+"_central_beamHalo.png")

outFile.Close()
dataFile.Close()
sigFile.Close()
