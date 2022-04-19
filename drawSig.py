import ROOT as rt
from array import array

sigFile = rt.TFile.Open("outSig_2DT.root","READ")
sigHistos = {}
years = ["2016","2017","2018"]
mX = ['15','40','55']
ctau = ['1000','10000']
colors = [1,2,3,4,6,7,38,41]
vars = ["MB1Ratio"]
vetoes = ["MB1Cluster","nonMB1Cluster"]

canvas = rt.TCanvas("canvas")
legend = rt.TLegend(0.55,0.6,0.8,0.8)
for veto in vetoes:
    canvas.Clear()
    for v in vars:
        for c in ctau:
            canvas.cd()
            canvas.Clear()
            legend.Clear()
            rt.gStyle.SetOptStat(0)
            nfile = 0
            sigStack = rt.THStack("sigStack",";"+v+";")
            for m in mX:
                sigHistos[m] = sigFile.Get('h_'+v+'_'+veto+'_'+m+'_'+c)
                sigHistos[m].SetLineWidth(2)
                sigHistos[m].SetLineColor(colors[nfile])
                if(sigHistos[m].Integral()>0.0):
                    sigHistos[m].Scale(1/sigHistos[m].Integral())
                    sigStack.Add(sigHistos[m],"hist")
                    legend.AddEntry(sigHistos[m],"mX="+m+" GeV","L")
                    
                nfile+=1
                
            sigStack.Draw("nostack")
            legend.Draw()
            sigStack.GetStack().Last().GetXaxis().SetTitle(v)
            canvas.Update()
            canvas.WaitPrimitive()
            canvas.SaveAs("plots/"+v+'_'+veto+'_ctau'+c+'sig.png')

sigFile.Close()
