import ROOT as rt

privateDir = "/storage/af/user/mcitron/hvHaddedFilesV3/"
centralDir = "/storage/af/user/mcitron/hvCentralProduction/"

years = ["2016","2017","2018"]
hierarchies = [("1","1"), ("2p5","1"), ("2p5", "2p5")]
portals = ["higgs","gluon","photon","darkphoton","vector"]
variables = ["nJets","jetPt","nDtRechitClusters","dtRechitClusterSize"]

colors = {"2016": 1, "2017": 2, "2018":3} 
tunes = {"2016": "_TuneCUETP8M1_13TeV-powheg-pythia8_RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_GEN-SIM-RECO_v2", "2017": "_TuneCP5_13TeV-powheg-pythia8_RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11_ext1-v1_GEN-SIM-RECO_v2", "2018": "_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18DRPremix-102X_upgrade2018_realistic_v15_ext1-v1_GEN-SIM-RECO_v2"}

canvas = rt.TCanvas("canvas")
legend = rt.TLegend(0.55,0.6,0.8,0.8)
h1 = rt.TH1D(); h1.SetLineColor(1); h1.SetLineStyle(7); h1.SetMarkerColor(1)
h2 = rt.TH1D(); h2.SetLineColor(2); h2.SetLineStyle(7); h2.SetMarkerColor(2)
h3 = rt.TH1D(); h3.SetLineColor(3); h3.SetLineStyle(7); h3.SetMarkerColor(3)
h4 = rt.TH1D(); h4.SetLineColor(1)
h5 = rt.TH1D(); h5.SetLineColor(2)
h6 = rt.TH1D(); h6.SetLineColor(3)
legend.AddEntry(h1,"2016 private","lep")
legend.AddEntry(h2,"2017 private","lep")
legend.AddEntry(h3,"2018 private","lep")
legend.AddEntry(h4,"2016 central","l")
legend.AddEntry(h5,"2016 central","l")
legend.AddEntry(h6,"2016 central","l")

firstPlot = True

for portal in portals:
    for xiO, xiL in hierarchies:
        if portal=="vector" and xiO!="1":
            continue
        for var in variables:
            firstPlot = True
            for year in years:
                privateFile = rt.TFile.Open(privateDir + year + "/HV_params_" + portal + "_m_15_ctau_5000mm_xiO_" + xiO + "_xiL_" + xiL + "_" + year + "_LLPNTUPLE_v3_filter_updatedKnapenCode_benchmarks.root")
                privateTree = privateFile.Get("ntuples/llp")

                privateTree.SetLineColor(colors[year])
                privateTree.SetLineStyle(7)
                privateTree.SetMarkerColor(colors[year])
                #h = rt.TH1D()
                #h.GetXaxis().SetTitle(var); h.SetTitle(portal+" "+year)
                if firstPlot:
                    rt.gStyle.SetOptStat(0)
                    privateTree.Draw(var,"","norm")
                    firstPlot = False
                else:
                    privateTree.Draw(var,"","same norm")

                h = rt.gPad.GetListOfPrimitives().At(0)
                h.GetXaxis().SetTitle(var); h.SetTitle(portal+" ("+xiO+", "+xiL+") ")

                centralFile = rt.TFile.Open(centralDir + "DarkShowerHiggs_" + portal + "_M15_CTAU5000_XIOMEGA" + xiO + "_XILAMBDA" + xiL + tunes[year] + "_generationForDT.root")
                centralTree = centralFile.Get("ntuples/llp")

                centralTree.SetLineColor(colors[year]);
                centralTree.Draw(var,"","same hist norm")
                
                privateFile.Close()
                centralFile.Close()

            canvas.SetTitle(portal+" " +var)
            legend.Draw()
            canvas.Update()
            canvas.SaveAs("plots/darkShowerComp/"+portal+"_"+xiO+"_"+xiL+"_"+var+".png")
            canvas.SetLogy()
            canvas.Update()
            canvas.SaveAs("plots/darkShowerComp/"+portal+"_"+xiO+"_"+xiL+"_"+var+"_logy.png")
            canvas.Clear()
