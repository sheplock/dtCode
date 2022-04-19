import ROOT as rt
from array import array
import numpy as np

dataFile = rt.TFile.Open("outData_ABCD_looseMuonVetoNoJetMET.root","READ")
years = ["2016","2017","2018"]
vars = ["matchedJetPt_nMB1Match"]
vetoes = ["invertedJetVetoNoClusterMET",
          "invertedJetVetoLooseMuonNoClusterMET",
          "invertedJetVetoAllMuonNoClusterMET"]
bins = [5,9,13,15];
j=15
while(j<=100):
    j=bins[-1]+2
    bins.append(j)


canvas = rt.TCanvas("canvas")
for veto in vetoes:
    canvas.Clear()
    for v in vars:
        canvas.cd()
        canvas.Clear()
        rt.gStyle.SetOptStat(0)
        dataHistos = []
        itr_year = 0
        for y in years:
            dataHistos.append(dataFile.Get("h_"+v+"_"+veto+"_"+y))
            itr_year+=1
        for i in range(1,len(dataHistos)):
            dataHistos[0].Add(dataHistos[i])
        py = []
        for j in range(len(bins)-1):
            py.append(dataHistos[0].ProjectionY("py_"+str(j),bins[j],bins[j+1]-1))
        x = np.zeros(100)
        y = np.zeros(100)
        ex = np.zeros(100)
        ey = np.zeros(100)
        for k in range(len(py)):
            x[k] = 0.5*(bins[k]+bins[k+1])*5
            ex[k] = 0.5*(bins[k+1]-bins[k])*5
            if(py[k].Integral(0,2)>0  and py[k].Integral(3,-1)>0):
                y[k] = py[k].Integral(0,2)/py[k].Integral(3,-1)
                ey[k] = y[k]*np.sqrt(1.0/py[k].Integral(0,2)+1.0/py[k].Integral(3,-1))
            elif(py[k].Integral(3,-1)>0):
                y[k] = py[k].Integral(0,2)/py[k].Integral(3,-1)
                ey[k] = y[k]*np.sqrt(1.0/py[k].Integral(3,-1)+1.0/1.3)
            elif(py[k].Integral(0,2)>0):
                y[k] = -1
                ey[k] = y[k]*np.sqrt(1.0/1.3+1.0/py[k].Integral(0,-2))
            else:
                y[k] = -1
                ey[k] = 1
        gr = rt.TGraphErrors(100,x,y,ex,ey)
        gr.GetXaxis().SetTitle("matched jet pT [GeV]")
        gr.GetYaxis().SetTitle("N(pass MB1 veto) / N(fail MB1 veto)")
        gr.SetTitle("")
        gr.Draw("ap")
        canvas.SaveAs("plots/"+v+"_"+veto+".C")

dataFile.Close()
