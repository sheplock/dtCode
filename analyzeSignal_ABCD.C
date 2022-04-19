#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector> 
#include <map> 
#include <fstream>      // std::ofstream
//#include "fastjet/ClusterSequence.hh"
using namespace std;
R__LOAD_LIBRARY(libTreePlayer)

int getStation(float hitX, float hitY){
  float hitR = sqrt(pow(hitX,2)+pow(hitY,2));
  if(hitR > 400. && hitR < 480.){ return 1; }
  else if(hitR > 485. && hitR < 560.){ return 2; }
  else if(hitR > 590. && hitR < 650.){ return 3; }
  else if(hitR > 690. && hitR < 800.){ return 4; }
  else{ return -1; }
}

int getWheel(float hitZ){
  if(hitZ > 0){
    if(hitZ < 127.){ return 0; }
    else if(hitZ < 395.){ return 1; }
    else if(hitZ < 661.){ return 2; }
    else{ return -99; }
  }
  else{
    return -1*getWheel(-1.0*hitZ);
  }
}

int getRPCLayer(float hitX, float hitY){
  float hitR = sqrt(pow(hitX,2)+pow(hitY,2));
  if(hitR > 410. && hitR < 440.){ return 1; }
  else if(hitR > 445. && hitR < 475.){ return 2; }
  else if(hitR > 490. && hitR < 520.){ return 3; }
  else if(hitR > 525. && hitR < 555.){ return 4; }
  else if(hitR > 600. && hitR < 630.){ return 5; }
  else if(hitR > 700. && hitR < 770.){ return 6; }
  else{ return -1; }
}

void analyzeSignal_ABCD(){

  TString name;
  std::vector<TString> mX = {"15","40","55"};
  //std::vector<TString> mX = {"55"};
  //std::vector<TString> mX = {"vector_m_2","vector_m_5","vector_m_10","vector_m_15","vector_m_20"};
  //std::vector<TString> mX = {"photon_m_2","photon_m_5","photon_m_10","photon_m_15","photon_m_20"};
  //std::vector<TString> mX = {"higgs_m_4","higgs_m_5","higgs_m_10","higgs_m_15","higgs_m_20"};
  //std::vector<TString> mX = {"gluon_m_3","gluon_m_5","gluon_m_10","gluon_m_15","gluon_m_20"};
  //std::vector<TString> mX = {"darkphoton_m_2","darkphoton_m_5","darkphoton_m_10","darkphoton_m_15","darkphoton_m_20"};
  //std::vector<TString> mX = {"photon_m_15","higgs_m_15","gluon_m_15","darkphoton_m_15","vector_m_15"};
  //char mX[2][10] = {"450"};
  //char ctau[2][20] = {"1m","10m"};
  //std::vector<TString> ctau = {"100","500","1000","2000","5000","10000","100000"};
  //std::vector<TString> ctau = {"3","10","30","100","300","1000","3000","10000","30000","100000"};
  //std::vector<TString> ctau = {"100000"};
  //std::vector<TString> ctau = {"1000","3000","10000","30000","100000"};
  std::vector<TString> ctau = {"100","1000","10000","100000"};
  //std::vector<TString> ctau = {"1000","10000"};
  //std::vector<TString> ctau = {"500mm_xiO_1_xiL_1","1000mm_xiO_1_xiL_1","5000mm_xiO_1_xiL_1","10000mm_xiO_1_xiL_1"};
  //std::vector<TString> ctau = {"500mm_xiO_1_xiL_1","500mm_xiO_2p5_xiL_1","500mm_xiO_2p5_xiL_2p5","1000mm_xiO_1_xiL_1","1000mm_xiO_2p5_xiL_1","1000mm_xiO_2p5_xiL_2p5",
  //			       "5000mm_xiO_1_xiL_1","5000mm_xiO_2p5_xiL_1","5000mm_xiO_2p5_xiL_2p5","10000mm_xiO_1_xiL_1","10000mm_xiO_2p5_xiL_1","10000mm_xiO_2p5_xiL_2p5"};
  //std::vector<TString> ctau = {"500mm_xi_1","500mm_xi_2p5","1000mm_xi_1","1000mm_xi_2p5","5000mm_xi_1","5000mm_xi_2p5","10000mm_xi_1","10000mm_xi_2p5"};
  //std::vector<TString> ctau = {"500mm_xiO_1_xiL_1","500mm_xiO_2p5_xiL_1","500mm_xiO_2p5_xiL_2p5"};
  //Float_t lifetime = 10000;
  std::vector<TString> years = {"MC_Fall18","MC_Fall17","MC_Summer16"};
  //std::vector<TString> years = {"2018","2017","2016"};
  //char years[3][20] = {"2018","2017","2016"};
  std::map<TString,Float_t> lumi;
  //lumi["2018"] = 59.74;
  //lumi["2017"] = 41.53;
  //lumi["2016"] = 35.92;
  lumi["MC_Fall18"] = 59.74;
  lumi["MC_Fall17"] = 41.53;
  lumi["MC_Summer16"] = 35.92;
  // Float_t lumi[3] = {59.74,41.53,35.92};
  Double_t weight = 1.0;
  Float_t decay1 = 0.0;
  Float_t decay2 = 0.0;
  Float_t ctau1 = 0.0;
  Float_t ctau2 = 0.0;

  //TString dir("/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/");
  //TString dir("/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/driftTube/V1p17/");
  //TString dir("/mnt/hadoop/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/driftTube/V1p15/");
  //TString dir("/storage/af/user/mcitron/hvHaddedFilesV3/");
  TString dir("/storage/af/user/mcitron/skims/signal/");
  TFile *_ofile = TFile::Open("outSig_ABCD_jetVeto10_test.root","RECREATE");

  TH1D *h_dtRechitClusterSize_signalRegion[5][12];
  TH1D *h_dtRechitClusterSize_MB2_signalRegion[5][12];
  TH1D *h_dtRechitClusterSize_MB3_signalRegion[5][12];
  TH1D *h_dtRechitClusterSize_MB4_signalRegion[5][12];
  TH1D *h_dtRechitClusterSize_signalRegionUnweighted[5][12];
  TH1D *h_dtRechitClusterSize_signalRegionNew[5][12];
  TH1D *h_dtRechitClusterSizeTotal_signalRegionNew[5][12];
  TH1D *h_dtRechitClusterSize_fullSelection_rpcCR[5][12];
  TH1D *h_dtRechitClusterSize_fullSelection_clusterMETCR[5][12];

  TH1D *h_nRPCMatched_fullVeto_clusterMETCR[5][12];
  TH1D *h_rpcSpread_fullVeto_clusterMETCR[5][12];
  TH1D *h_rpcBx_fullVeto_clusterMETCR[5][12];
  TH1D *h_dPhiJetMET_fullVeto_clusterMETCR[5][12];
  TH1D *h_dtRechitClusterMaxStation_fullVeto_clusterMETCR[5][12];

  TH1D *h_nRPCMatched_Nminus1_clusterMETCR[5][12];
  TH1D *h_rpcSpread_Nminus1_clusterMETCR[5][12];
  TH1D *h_rpcBx_Nminus1_clusterMETCR[5][12];
  TH1D *h_dPhiJetMET_Nminus1_clusterMETCR[5][12];
  TH1D *h_dtRechitClusterMaxStation_Nminus1_clusterMETCR[5][12];
  TH1D *h_jetVetoPt_Nminus1_clusterMETCR[5][12];
  TH1D *h_muonVetoPt_Nminus1_clusterMETCR[5][12];
  TH1D *h_muonLooseIDVetoPt_Nminus1_clusterMETCR[5][12];
  
  TH1D *h_dPhiClusterMET_fullVeto_rpcCR[5][12];
  TH1D *h_dPhiJetMET_fullVeto_rpcCR[5][12];
  TH1D *h_dtRechitClusterMaxStation_fullVeto_rpcCR[5][12];

  TH1D *h_dPhiClusterMET_Nminus1_rpcCR[5][12];
  TH1D *h_dPhiJetMET_Nminus1_rpcCR[5][12];
  TH1D *h_dtRechitClusterMaxStation_Nminus1_rpcCR[5][12];

  TH1D *h_nStations1_50hits_clusterMETCR[5][12];
  TH1D *h_nStations1_100hits_clusterMETCR[5][12];
  TH1D *h_nStations1_150hits_clusterMETCR[5][12];
  TH1D *h_nStations25_50hits_clusterMETCR[5][12];
  TH1D *h_nStations25_100hits_clusterMETCR[5][12];
  TH1D *h_nStations25_150hits_clusterMETCR[5][12];
  TH1D *h_nStations50_50hits_clusterMETCR[5][12];
  TH1D *h_nStations50_100hits_clusterMETCR[5][12];
  TH1D *h_nStations50_150hits_clusterMETCR[5][12];

  TH1D *h_nWheels1_50hits_clusterMETCR[5][12];
  TH1D *h_nWheels25_50hits_clusterMETCR[5][12];
  TH1D *h_nWheels50_50hits_clusterMETCR[5][12];
  TH1D *h_nWheels1_100hits_clusterMETCR[5][12];
  TH1D *h_nWheels25_100hits_clusterMETCR[5][12];
  TH1D *h_nWheels50_100hits_clusterMETCR[5][12];
  TH1D *h_nWheels1_150hits_clusterMETCR[5][12];
  TH1D *h_nWheels25_150hits_clusterMETCR[5][12];
  TH1D *h_nWheels50_150hits_clusterMETCR[5][12];

  TH1D *h_nDtSegsStation_50hits_clusterMETCR[5][12];
  TH1D *h_nDtSegsStation_100hits_clusterMETCR[5][12];
  TH1D *h_nDtSegsStation_150hits_clusterMETCR[5][12];
  TH1D *h_nDtSegsWheel_50hits_clusterMETCR[5][12];
  TH1D *h_nDtSegsWheel_100hits_clusterMETCR[5][12];
  TH1D *h_nDtSegsWheel_150hits_clusterMETCR[5][12];
  TH1D *h_nDtSegsChamber_50hits_clusterMETCR[5][12];
  TH1D *h_nDtSegsChamber_100hits_clusterMETCR[5][12];
  TH1D *h_nDtSegsChamber_150hits_clusterMETCR[5][12];
  TH1D *h_nDtSegs_50hits_clusterMETCR[5][12];
  TH1D *h_nDtSegs_100hits_clusterMETCR[5][12];
  TH1D *h_nDtSegs_150hits_clusterMETCR[5][12];

  TH1D *h_nStations1Seg_50hits_clusterMETCR[5][12];
  TH1D *h_nStations5Seg_50hits_clusterMETCR[5][12];
  TH1D *h_nStations10Seg_50hits_clusterMETCR[5][12];
  TH1D *h_nStations1Seg_100hits_clusterMETCR[5][12];
  TH1D *h_nStations5Seg_100hits_clusterMETCR[5][12];
  TH1D *h_nStations10Seg_100hits_clusterMETCR[5][12];
  TH1D *h_nStations1Seg_150hits_clusterMETCR[5][12];
  TH1D *h_nStations5Seg_150hits_clusterMETCR[5][12];
  TH1D *h_nStations10Seg_150hits_clusterMETCR[5][12];

  TH1D *h_nWheels1Seg_50hits_clusterMETCR[5][12];
  TH1D *h_nWheels5Seg_50hits_clusterMETCR[5][12];
  TH1D *h_nWheels10Seg_50hits_clusterMETCR[5][12];
  TH1D *h_nWheels1Seg_100hits_clusterMETCR[5][12];
  TH1D *h_nWheels5Seg_100hits_clusterMETCR[5][12];
  TH1D *h_nWheels10Seg_100hits_clusterMETCR[5][12];
  TH1D *h_nWheels1Seg_150hits_clusterMETCR[5][12];
  TH1D *h_nWheels5Seg_150hits_clusterMETCR[5][12];
  TH1D *h_nWheels10Seg_150hits_clusterMETCR[5][12];


  TH1D *h_nRPCStations1_50hits_clusterMETCR[5][12];
  TH1D *h_nRPCStations5_50hits_clusterMETCR[5][12];
  TH1D *h_nRPCStations10_50hits_clusterMETCR[5][12];
  TH1D *h_nRPCStations1_100hits_clusterMETCR[5][12];
  TH1D *h_nRPCStations5_100hits_clusterMETCR[5][12];
  TH1D *h_nRPCStations10_100hits_clusterMETCR[5][12];
  TH1D *h_nRPCStations1_150hits_clusterMETCR[5][12];
  TH1D *h_nRPCStations5_150hits_clusterMETCR[5][12];
  TH1D *h_nRPCStations10_150hits_clusterMETCR[5][12];

  TH1D *h_nRPCWheels1_50hits_clusterMETCR[5][12];
  TH1D *h_nRPCWheels5_50hits_clusterMETCR[5][12];
  TH1D *h_nRPCWheels10_50hits_clusterMETCR[5][12];
  TH1D *h_nRPCWheels1_100hits_clusterMETCR[5][12];
  TH1D *h_nRPCWheels5_100hits_clusterMETCR[5][12];
  TH1D *h_nRPCWheels10_100hits_clusterMETCR[5][12];
  TH1D *h_nRPCWheels1_150hits_clusterMETCR[5][12];
  TH1D *h_nRPCWheels5_150hits_clusterMETCR[5][12];
  TH1D *h_nRPCWheels10_150hits_clusterMETCR[5][12];

  TH1D *h_minSegmentDR_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegments_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsClusterStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsClusterStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsClusterStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsOtherStations_invertedShowerVetoes[5][12];
  TH1D *h_nAlignedSegments_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStation_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsOuterStation_invertedShowerVetoes[5][12];
  TH1D *h_segmentAlignmentDeltaPhi_invertedShowerVetoes[5][12];
  TH1D *h_segmentAlignmentDeltaEta_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsOuterStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsOuterStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsOuterStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedShowerVetoes[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedShowerVetoes[5][12];
  TH1D *h_muonVetoOutcomes_invertedShowerVetoes[5][12];
  TH1D *h_muonVetoOutcomesMB2_invertedShowerVetoes[5][12];
  TH1D *h_muonVetoOutcomesMB3_invertedShowerVetoes[5][12];
  TH1D *h_muonVetoOutcomesMB4_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomes_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesMB2_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesMB3_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesMB4_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassSegMuon_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassOneSegMuon_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[5][12];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[5][12];
  TH1D *h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[5][12];
  TH1D *h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[5][12];
  TH1D *h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[5][12];
  TH1D *h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimes_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMean_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesClusterStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesClusterStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesClusterStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesInnerStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesInnerStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesInnerStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesOuterStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesOuterStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimesOuterStationMB4_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB2_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB3_invertedShowerVetoes[5][12];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB4_invertedShowerVetoes[5][12];
  
  TH1D *h_minSegmentDR_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegments_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsClusterStationMB2_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsClusterStationMB3_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsClusterStationMB4_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsOtherStations_invertedJetVeto[5][12];
  TH1D *h_nAlignedSegments_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStation_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsOuterStation_invertedJetVeto[5][12];
  TH1D *h_segmentAlignmentDeltaPhi_invertedJetVeto[5][12];
  TH1D *h_segmentAlignmentDeltaEta_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB2_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB3_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB4_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB2_invertedMB1_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB3_invertedMB1_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationMB4_invertedMB1_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsOuterStationMB2_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsOuterStationMB3_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsOuterStationMB4_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedJetVeto[5][12];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedJetVeto[5][12];
  TH1D *h_muonVetoOutcomes_invertedJetVeto[5][12];
  TH1D *h_muonVetoOutcomesMB2_invertedJetVeto[5][12];
  TH1D *h_muonVetoOutcomesMB3_invertedJetVeto[5][12];
  TH1D *h_muonVetoOutcomesMB4_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomes_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesMB2_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesMB3_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesMB4_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassSegMuon_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassOneSegMuon_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[5][12];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[5][12];
  TH1D *h_AdjacentMB1VetoOutcomes_invertedJetVeto[5][12];
  TH1D *h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[5][12];
  TH1D *h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[5][12];
  TH1D *h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimes_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMean_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesClusterStationMB2_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesClusterStationMB3_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesClusterStationMB4_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB2_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB3_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB4_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesInnerStationMB2_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesInnerStationMB3_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesInnerStationMB4_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB2_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB3_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB4_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesOuterStationMB2_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesOuterStationMB3_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimesOuterStationMB4_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB2_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB3_invertedJetVeto[5][12];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB4_invertedJetVeto[5][12];
  
  TH1D *h_dtRechitClusterJetVetoPt[5][12];
  TH1D *h_dtRechitClusterMuonVetoPt[5][12];
  TH1D *h_dtRechitClusterMB1Veto[5][12];

  TH1D *h_efficiency[5][12];
  TH1D *h_efficiency_MB1CR[5][12];
  
  TH1D *h_nDtRechitClusters_dPhiJetMET[5][12];
  TH1D *h_nDtRechitClustersVeto_dPhiJetMET[5][12];
  TH1D *h_dtRechitClustersDR_dPhiJetMET[5][12];
  TH1D *h_dtRechitClustersVetoDR_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterNSegmentStation2_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterNSegmentStation3_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterNSegmentStation4_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterMaxStationRatio_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterNStation_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterMaxStation_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterMaxChamberRatio_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterNChamber_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterMaxChamber_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterX_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterY_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterZ_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterEta_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterPhi_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterTime_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterXSpread_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterYSpread_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterZSpread_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterEtaSpread_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterPhiSpread_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterTimeSpread_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterMajorAxis_dPhiJetMET[5][12];
  TH1D *h_dtRechitClusterMinorAxis_dPhiJetMET[5][12];
  
  TH1D *h_nRB1Match_dPhiJetMET[5][12];
  TH1D *h_nRB1Match_MB1Veto_dPhiJetMET[5][12];
  TH1D *h_nMB1MatchAdjacent_dPhiJetMET[5][12];
  TH1D *h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[5][12];
  TH1D *h_nMB1MatchAdjacent_dPhiClusterMET[5][12];
  TH1D *h_nMB1MatchAdjacent_MB1Veto_dPhiClusterMET[5][12];
  TH1D *h_nMB1MatchPi2_dPhiClusterMET[5][12];
  TH1D *h_nMB1MatchAdjacentPi2_dPhiClusterMET[5][12];
  TH1D *h_nMB1MatchAdjacentPi2_MB1Veto_dPhiClusterMET[5][12];
  TH1D *h_nMB1MatchAdjacent0p8_dPhiClusterMET[5][12];
  TH1D *h_nMB1MatchAdjacent0p8_MB1Veto_dPhiClusterMET[5][12];
  TH1D *h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[5][12];
  TH1D *h_nMB1MatchAdjacent0p8Pi2_MB1Veto_dPhiClusterMET[5][12];
  
  TH1D *h_jetChargedHadronicEnergyFraction_SR[5][12];
  TH1D *h_jetNeutralHadronicEnergyFraction_SR[5][12];
  TH1D *h_jetNeutralEMEnergyFraction_SR[5][12];
  TH1D *h_jetChargedEMEnergyFraction_SR[5][12];
  TH1D *h_leadingJetChargedHadronicEnergyFraction_SR[5][12];
  TH1D *h_leadingJetNeutralHadronicEnergyFraction_SR[5][12];
  TH1D *h_leadingJetNeutralEMEnergyFraction_SR[5][12];
  TH1D *h_leadingJetChargedEMEnergyFraction_SR[5][12];

  TH1D *h_decayVertexRadius_noVeto[5][12];
  TH1D *h_decayVertexRadius_clusterReco[5][12];
  TH1D *h_decayVertexRadius_clusterReco_signalRegionEq2[5][12];
  TH1D *h_decayVertexRadius_clusterReco_signalRegionGt2[5][12];

  TH1D *h_decayVertexZ_noVeto[5][12];
  TH1D *h_decayVertexZ_clusterReco[5][12];
  TH1D *h_decayVertexZ_clusterReco_signalRegionEq2[5][12];
  TH1D *h_decayVertexZ_clusterReco_signalRegionGt2[5][12];
  
  TH1D *h_leadJetPt_MB1HitsCR[5][12];
  TH1D *h_leadJetPt_MB2CR[5][12];
  TH1D *h_leadJetPt_MB1CR[5][12];
  TH1D *h_leadJetPt_MB2withMB1CR[5][12];
  TH1D *h_leadJetPt_SR[5][12];
  TH1D *h_leadJetPtMET_MB1HitsCR[5][12];
  TH1D *h_leadJetPtMET_MB2CR[5][12];
  TH1D *h_leadJetPtMET_MB1CR[5][12];
  TH1D *h_leadJetPtMET_MB2withMB1CR[5][12];
  TH1D *h_leadJetPtMET_SR[5][12];
  TH1D *h_leadJetEta_SR[5][12];

  TH1D *h_matchedJetNHF[5][12];
  TH1D *h_matchedJetNEMF[5][12];
  TH1D *h_matchedJetCHF[5][12];
  TH1D *h_matchedJetCEMF[5][12];
  TH1D *h_matchedJetNConst[5][12];
  TH1D *h_matchedJetChargedMult[5][12];

  TH1D *h_nCAClusters[5][12];
  TH1D *h_sizeCACluster[5][12];
  TH1D *h_radiusCACluster[5][12];
  TH1D *h_zCACluster[5][12];
  TH1D *h_phiCACluster[5][12];

  TH1D *h_nDtRechitClusters[5][12];
  TH1D *h_dtRechitClusterSize[5][12];
  TH1D *h_dtRechitClusterR[5][12];
  TH1D *h_dtRechitClusterZ[5][12];
  TH1D *h_dtRechitClusterPhi[5][12];

  TH1D *h_jetOverlap10GeV[5][12];
  TH1D *h_jetOverlap20GeV[5][12];
  TH1D *h_muonOverlap10GeVLoose[5][12];
  TH1D *h_muonOverlap10GeVTight[5][12];

  TH1D *h_MB1Ratio[5][12];
  TH1D *h_MB1Ratio_100hits[5][12];

  Bool_t overlapJet10GeV = false;
  Bool_t overlapJet20GeV = false;
  Bool_t overlapMuon10GeVLoose = false;
  Bool_t overlapMuon10GeVTight = false;

  //vector<fastjet::PseudoJet> dtPoints;
  //vector<fastjet::PseudoJet> clustersCA;
  //vector<fastjet::PseudoJet> constituents;
  Int_t nCAClusters = 0;
  Int_t SRyield_CA = 0;
  Int_t SRyieldMB2_CA = 0;
  Bool_t CAclusterJetVeto = false;
  Bool_t CAclusterMuonVeto = false;
  Int_t CAclusterMB1Match = 0;
  Int_t CAclusterAdjacentMB1MatchPlus = 0;
  Int_t CAclusterAdjacentMB1MatchMinus = 0;
  Int_t CAclusterRPCMatch = 0;
  Double_t CAclusterdPhiMET = 0.0;
  Int_t CAclusterSize = 0;
  Int_t CAclusterStation = 0;
  Int_t CAclusterWheel = 0;
  Int_t CAclusterStationHits[4] = {0,0,0,0};
  Int_t CAclusterWheelHits[5] = {0,0,0,0,0};

  Double_t chargedHadFraction_mindPhi = 0.0;
  Double_t chargedEMFraction_mindPhi = 0.0;
  Double_t neutralHadFraction_mindPhi = 0.0;
  Double_t neutralEMFraction_mindPhi = 0.0;

  Double_t dPhi_tmp = 0.0;
  Double_t dPhi_min = 0.0;
  Double_t dPhiClusterRPC = 0.0;
  Double_t dZClusterRPC = 0.0;
  Double_t dPhiClusterMET_max = 0.0;
  Double_t dPhiClusterMET = 0.0;
  vector<Int_t> rpcBx = {};
  Int_t rpcSpread = 0;
  Double_t rpcMedian = 0;

  vector<Double_t> clusterEta = {};
  vector<Double_t> clusterPhi = {};
  vector<Int_t> clusterSize = {};
  Int_t clusterSizeTotal = 0;
  Int_t nClustersVeto_dPhiJetMET = 0;

  Int_t evtNum = 0;
  Int_t totalNum = 0;
  Bool_t HLT = false;

  Bool_t passMuon = false;
  Bool_t passMuonLoose = false;
  Bool_t passMuon_alt = false;
  Bool_t passMB1 = false;
  Bool_t passJet = false;
  Bool_t passJetTightId = false;
  Double_t hoMatchedEnergy = 0.0;
  
  Bool_t passClusterCR = false;
  Bool_t passNoVeto_clusterCR = false;
  Bool_t passFullVeto_clusterCR = false;
  Bool_t passRPCMatch_clusterCR = false;
  Bool_t passRPCSpread_clusterCR = false;
  Bool_t passRPCBx_clusterCR = false;
  Bool_t passMaxStation_clusterCR = false;
  Bool_t passLepton_clusterCR = false;
  Bool_t pass50Hits_clusterCR = false;
  Bool_t pass25Hits_clusterCR = false;

  Bool_t passNoVeto = false;  
  Bool_t passFullVeto_rpcCR = false;
  Bool_t passRPCCR = false;
  Bool_t passClusterMET_rpcCR = false;
  Bool_t passMaxStation_rpcCR = false;
  Bool_t passJetMET_rpcCR = false;
  Bool_t passLepton_rpcCR = false;
  Bool_t pass50Hits_rpcCR = false;
  Bool_t pass25Hits_rpcCR = false;
  Bool_t passFullPlus = false;

  Bool_t passMET = false;
  Bool_t passOneJet = false;
  Bool_t passNHFJet = false;
  Bool_t passJetMET = false;
  Bool_t passStations25 = false;
  Bool_t passWheels25 = false;
  Bool_t passJetVeto = false;
  Bool_t passJetTightIdVeto = false;
  Bool_t passMuonVeto = false;
  Bool_t passMB1Veto = false;
  Bool_t passMaxStation = false;
  Bool_t passClusterMET = false;
  Bool_t passRPCMatch = false;
  Bool_t passRPCSpread = false;
  Bool_t passRPCBx = false;
  Bool_t passNoVetoCluster = false;
  Bool_t passClusterSize = false;
  Bool_t passClusterSizeA = false;
  Bool_t passClusterSizeB = false;
  Bool_t passClusterSizeC = false;
  Bool_t passClusterMETC = false;
  Bool_t passClusterMETB = false;
  Bool_t passClusterMETA = false;
  Bool_t passAdjacentMB1 = false;
  Bool_t passAdjacent0p8MB1 = false;
  Bool_t passOtherStations = false;
  Bool_t passCscCluster = false;
  Bool_t passLargeCscCluster = false;
  Bool_t passCscClusterME11Veto = false;
  Bool_t passLargeCscClusterME11Veto = false;

  Bool_t passMB1CR = false;
  Bool_t passJetVetoMB1CR = false;
  Bool_t passMuonVetoMB1CR = false;
  Bool_t passMuonVetoLooseMB1CR = false;
  Bool_t passRpcMatchMB1CR = false;
  Bool_t passClusterMETMB1CR = false;
  Bool_t passClusterSizeMB1CR = false;

  Bool_t passMB2Cluster = false;
  Bool_t passMB2CR = false;
  Bool_t passMB2CRwithAdjacent = false;
  Bool_t passMB2CRwithAdjacent0p8 = false;
  Bool_t passMB2CRwithOther = false;
  Bool_t passMB2CRwithNHF = false;
  Bool_t passMB2CRwithNHFnoRPC = false;
  Bool_t passSRwithNHFnoRPC = false;
  Double_t nPassMB2CR = 0;
  Double_t nPassMB2CRwithAdjacent = 0;
  Double_t nPassMB2CRwithAdjacent0p8 = 0;
  Double_t nPassMB2CRwithOther = 0;
  Double_t nPassMB2CRwithNHF = 0;
  Double_t nPassMB2CRwithNHFnoRPC = 0;
  Double_t nPassSRwithNHFnoRPC = 0;

  Bool_t passSignalRegion = false;
  Double_t nPassSignalRegion = 0;

  Double_t nPassClusterCR = 0;
  Double_t nPassNoVeto_clusterCR = 0;
  Double_t nPassFullVeto_clusterCR = 0;
  Double_t nPassRPCMatch_clusterCR = 0;
  Double_t nPassRPCSpread_clusterCR = 0;
  Double_t nPassRPCBx_clusterCR = 0;
  Double_t nPassMaxStation_clusterCR = 0;
  Double_t nPassLepton_clusterCR = 0;
  Double_t nPass50Hits_clusterCR = 0;
  Double_t nPass25Hits_clusterCR = 0;

  Double_t nPassNoVeto = 0;  
  Double_t nPassFullVeto_rpcCR = 0;
  Double_t nPassRPCCR = 0;
  Double_t nPassClusterMET_rpcCR = 0;
  Double_t nPassMaxStation_rpcCR = 0;
  Double_t nPassFullPlus = 0;
  Double_t nPassJetMET_rpcCR = 0;
  Double_t nPassLepton_rpcCR = 0;
  Double_t nPass50Hits_rpcCR = 0;
  Double_t nPass25Hits_rpcCR = 0;

  Int_t nStations1 = 0;
  Int_t nStations25 = 0;
  Int_t nStations50 = 0;
  Int_t nStations1Seg = 0;
  Int_t nStations5Seg = 0;
  Int_t nStations10Seg = 0;
  Int_t hitStation1 = 0;
  Int_t hitStation2 = 0;
  Int_t hitStation3 = 0;
  Int_t hitStation4 = 0;
  Int_t segStation1 = 0;
  Int_t segStation2 = 0;
  Int_t segStation3 = 0;
  Int_t segStation4 = 0;

  Int_t nWheels1 = 0;
  Int_t nWheels25 = 0;
  Int_t nWheels50 = 0;
  Int_t nWheels1Seg = 0;
  Int_t nWheels5Seg = 0;
  Int_t nWheels10Seg = 0;
  Int_t hitWheelm1 = 0;
  Int_t hitWheelm2 = 0;
  Int_t hitWheel0 = 0;
  Int_t hitWheel1 = 0;
  Int_t hitWheel2 = 0;
  Int_t segWheelm1 = 0;
  Int_t segWheelm2 = 0;
  Int_t segWheel0 = 0;
  Int_t segWheel1 = 0;
  Int_t segWheel2 = 0;
  Int_t segChambers[4][5] = {0};

  Int_t nRPCStations1 = 0;
  Int_t nRPCStations5 = 0;
  Int_t nRPCStations10 = 0;
  Int_t hitRPCStation1 = 0;
  Int_t hitRPCStation2 = 0;
  Int_t hitRPCStation3 = 0;
  Int_t hitRPCStation4 = 0;

  Int_t nRPCWheels1 = 0;
  Int_t nRPCWheels5 = 0;
  Int_t nRPCWheels10 = 0;
  Int_t hitRPCWheelm1 = 0;
  Int_t hitRPCWheelm2 = 0;
  Int_t hitRPCWheel0 = 0;
  Int_t hitRPCWheel1 = 0;
  Int_t hitRPCWheel2 = 0;

  Double_t minMB1SegmentDR = 999.;
  Double_t minSegmentDR = 999.;
  Int_t nAlignedSegments = 0;
  Int_t nMB1SegMatchCluster = 0;
  Int_t nMB2SegMatchCluster = 0;
  Int_t nMB3SegMatchCluster = 0;
  Int_t nMB4SegMatchCluster = 0;
  Double_t meanSegTime = -999.;
  Double_t meanMB1SegTime = -999.;
  Double_t meanMB2SegTime = -999.;
  Double_t meanMB3SegTime = -999.;
  Double_t meanMB4SegTime = -999.;

  Float_t dtR = 0.0;
  Int_t rpcStation = 99;
  Int_t rpcWheel = 99;
  Int_t dtStation = 99;
  Int_t dtWheel = 99;
  Int_t maxClusterSize = 0;
  Int_t maxClusterSizeSR = 0;
  Int_t maxClusterSizeMB1CR = 0;
  Int_t hitsMB1 = 0;
  Int_t nMB1MatchClusterAdjacentPlus = 0;
  Int_t nMB1MatchClusterAdjacentMinus = 0;
  Int_t nMB1MatchPi2AdjacentPlus = 0;
  Int_t nMB1MatchPi2 = 0;
  Int_t nMB1MatchPi2AdjacentMinus = 0;
  Int_t nMB1MatchClusterAdjacent0p8Plus = 0;
  Int_t nMB1MatchClusterAdjacent0p8Minus = 0;
  Int_t nMB1MatchPi2Adjacent0p8Plus = 0;
  Int_t nMB1MatchPi2Adjacent0p8Minus = 0;
  Int_t nRB1MatchCluster = 0;
  Int_t nMB1SegMatchClusterAdjacent0p8Plus = 0;
  Int_t nMB1SegMatchClusterAdjacent0p8Minus = 0;

  Int_t nMB1MatchCluster = 0;
  Int_t nMB2MatchCluster = 0;
  Int_t nMB3MatchCluster = 0;
  Int_t nMB4MatchCluster = 0;

  Double_t matchedLooseIDMuonPt = 0;

  TRandom3 *rand = new TRandom3();
  Int_t pmRand = 0;

  ofstream myout;
  myout.open("test.csv");
  myout << "mX,ctau,MET,OneJet,JetMET,Stations25,Wheels25,Cluster,NHFJet,JetVeto,MuonVeto,MB1Veto,MaxStation,RPCMatch,RPCSpread,RPCBx,AdjacentMB1,ClusterMET,Size,AdjacentMB1,NHFJet,MB2,MB2Unc,MB34,MB34Unc,C1,mc_unc_C1,C2,mc_unc_C2,B1,mc_unc_B1,B2,mc_unc_B2,A1,mc_unc_A1,A2,mc_unc_A2\n";
  
  for(Int_t itr_mX=0; itr_mX<mX.size(); itr_mX++){
    for(Int_t itr_ctau=0; itr_ctau<ctau.size(); itr_ctau++){
      name = "h_nDtRechitClusters_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nDtRechitClusters_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
      
      name = "h_nDtRechitClustersVeto_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nDtRechitClustersVeto_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
      
      name = "h_dtRechitClustersDR_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClustersDR_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,2);
      
      name = "h_dtRechitClustersVetoDR_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClustersVetoDR_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,2);
      
      name = "h_dtRechitClusterNSegmentStation2_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterNSegmentStation2_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,500);
      
      name = "h_dtRechitClusterNSegmentStation3_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterNSegmentStation3_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,500);
      
      name = "h_dtRechitClusterNSegmentStation4_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterNSegmentStation4_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,500);

      name = "h_dtRechitClusterMaxStation_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterMaxStation_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

      name = "h_dtRechitClusterNStation_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterNStation_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

      name = "h_dtRechitClusterMaxStationRatio_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterMaxStationRatio_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

      name = "h_dtRechitClusterMaxChamberRatio_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterMaxChamberRatio_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);
      
      name = "h_dtRechitClusterNChamber_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterNChamber_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

      name = "h_dtRechitClusterMaxChamber_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterMaxChamber_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

      name = "h_dtRechitClusterX_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterX_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",100,-800,800);

      name = "h_dtRechitClusterY_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterY_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",100,-800,800);

      name = "h_dtRechitClusterZ_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterZ_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",100,-600,600);

      name = "h_dtRechitClusterEta_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterEta_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,-1.5,1.5);

      name = "h_dtRechitClusterPhi_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterPhi_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",100,-3.5,3.5);

      name = "h_dtRechitClusterTime_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterTime_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,300,800);

      name = "h_dtRechitClusterXSpread_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterXSpread_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",100,0,500);

      name = "h_dtRechitClusterYSpread_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterYSpread_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",100,0,500);

      name = "h_dtRechitClusterZSpread_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterZSpread_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",100,0,140);

      name = "h_dtRechitClusterEtaSpread_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterEtaSpread_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,0.5);

      name = "h_dtRechitClusterPhiSpread_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterPhiSpread_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

      name = "h_dtRechitClusterTimeSpread_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterTimeSpread_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",100,0,150);

      name = "h_dtRechitClusterMajorAxis_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterMajorAxis_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

      name = "h_dtRechitClusterMinorAxis_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterMinorAxis_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);
      
      name = "h_dtRechitClusterSize_signalRegion_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterSize_signalRegion[itr_mX][itr_ctau] = new TH1D(name,"",250,0,500);

      name = "h_dtRechitClusterSize_MB2_signalRegion_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterSize_MB2_signalRegion[itr_mX][itr_ctau] = new TH1D(name,"",100,1,501);

      name = "h_dtRechitClusterSize_MB3_signalRegion_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterSize_MB3_signalRegion[itr_mX][itr_ctau] = new TH1D(name,"",100,1,501);

      name = "h_dtRechitClusterSize_MB4_signalRegion_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterSize_MB4_signalRegion[itr_mX][itr_ctau] = new TH1D(name,"",100,1,501);

      name = "h_dtRechitClusterSize_signalRegionUnweighted_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterSize_signalRegionUnweighted[itr_mX][itr_ctau] = new TH1D(name,"",250,0,500);

      name = "h_dtRechitClusterSize_signalRegionNew_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterSize_signalRegionNew[itr_mX][itr_ctau] = new TH1D(name,"",250,0,500);

      name = "h_dtRechitClusterSizeTotal_signalRegionNew_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterSizeTotal_signalRegionNew[itr_mX][itr_ctau] = new TH1D(name,"",250,0,500);
      
      name = "h_dtRechitClusterSize_fullSelection_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_dtRechitClusterSize_fullSelection_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",250,0,500);


      name = "h_nRB1Match_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nRB1Match_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",30,0,30);
      
      name = "h_nRB1Match_MB1Veto_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nRB1Match_MB1Veto_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",30,0,30);
      
      name = "h_nMB1MatchAdjacent_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacent_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,50);
      
      name = "h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_mX][itr_ctau] = new TH1D(name,"",30,0,30);

      name = "h_nMB1MatchAdjacent_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacent_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,50);
      
      name = "h_nMB1MatchAdjacent_MB1Veto_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacent_MB1Veto_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",30,0,30);

      name = "h_nMB1MatchAdjacentPi2_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,50);

      name = "h_nMB1MatchPi2_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchPi2_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",3,0,3);
      
      name = "h_nMB1MatchAdjacentPi2_MB1Veto_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacentPi2_MB1Veto_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",30,0,30);

      name = "h_nMB1MatchAdjacent0p8_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,50);
      
      name = "h_nMB1MatchAdjacent0p8_MB1Veto_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacent0p8_MB1Veto_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",30,0,30);

      name = "h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",50,0,50);
      
      name = "h_nMB1MatchAdjacent0p8Pi2_MB1Veto_dPhiClusterMET_"+mX[itr_mX]+"_"+ctau[itr_ctau];
      h_nMB1MatchAdjacent0p8Pi2_MB1Veto_dPhiClusterMET[itr_mX][itr_ctau] = new TH1D(name,"",30,0,30);


    name = "h_dtRechitClusterSize_fullSelection_rpcCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dtRechitClusterSize_fullSelection_rpcCR[itr_mX][itr_ctau] = new TH1D(name,"",250,0,500);

    name = "h_nRPCMatched_fullVeto_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCMatched_fullVeto_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_rpcSpread_fullVeto_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_rpcSpread_fullVeto_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",10,0,10);

    name = "h_rpcBx_fullVeto_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_rpcBx_fullVeto_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_dPhiJetMET_fullVeto_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dPhiJetMET_fullVeto_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",35,0,3.5);

    name = "h_dtRechitClusterMaxStation_fullVeto_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dtRechitClusterMaxStation_fullVeto_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCMatched_Nminus1_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCMatched_Nminus1_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_rpcSpread_Nminus1_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_rpcSpread_Nminus1_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",10,0,10);

    name = "h_rpcBx_Nminus1_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_rpcBx_Nminus1_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_jetVetoPt_Nminus1_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_jetVetoPt_Nminus1_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",20,0,200);

    name = "h_muonVetoPt_Nminus1_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoPt_Nminus1_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",20,0,200);

    name = "h_muonLooseIDVetoPt_Nminus1_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonLooseIDVetoPt_Nminus1_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",20,0,200);

    name = "h_dPhiJetMET_Nminus1_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dPhiJetMET_Nminus1_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",35,0,3.5);

    name = "h_dtRechitClusterMaxStation_Nminus1_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dtRechitClusterMaxStation_Nminus1_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_dPhiClusterMET_fullVeto_rpcCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dPhiClusterMET_fullVeto_rpcCR[itr_mX][itr_ctau] = new TH1D(name,"",35,0,3.5);

    name = "h_dPhiJetMET_fullVeto_rpcCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dPhiJetMET_fullVeto_rpcCR[itr_mX][itr_ctau] = new TH1D(name,"",35,0,3.5);

    name = "h_dtRechitClusterMaxStation_fullVeto_rpcCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dtRechitClusterMaxStation_fullVeto_rpcCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_dPhiClusterMET_Nminus1_rpcCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dPhiClusterMET_Nminus1_rpcCR[itr_mX][itr_ctau] = new TH1D(name,"",35,0,3.5);

    name = "h_dPhiJetMET_Nminus1_rpcCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dPhiJetMET_Nminus1_rpcCR[itr_mX][itr_ctau] = new TH1D(name,"",35,0,3.5);

    name = "h_dtRechitClusterMaxStation_Nminus1_rpcCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dtRechitClusterMaxStation_Nminus1_rpcCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);



    name = "h_nStations1_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations1_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations1_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations1_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations1_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations1_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations25_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations25_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations25_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations25_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations25_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations25_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations50_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations50_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations50_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations50_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations50_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations50_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nWheels1_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels1_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels1_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels1_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels1_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels1_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels25_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels25_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels25_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels25_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels25_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels25_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels50_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels50_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels50_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels50_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels50_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels50_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);


    name = "h_nDtSegsStation_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsStation_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsStation_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsStation_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsStation_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsStation_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsWheel_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsWheel_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsWheel_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsWheel_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsWheel_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsWheel_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsChamber_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsChamber_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsChamber_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsChamber_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsChamber_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegsChamber_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegs_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegs_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,50);

    name = "h_nDtSegs_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegs_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,50);

    name = "h_nDtSegs_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nDtSegs_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,50);

    name = "h_nStations1Seg_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations1Seg_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations1Seg_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations1Seg_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations1Seg_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations1Seg_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations5Seg_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations5Seg_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations5Seg_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations5Seg_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations5Seg_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations5Seg_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations10Seg_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations10Seg_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations10Seg_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations10Seg_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nStations10Seg_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nStations10Seg_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nWheels1Seg_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels1Seg_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels1Seg_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels1Seg_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels1Seg_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels1Seg_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels5Seg_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels5Seg_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels5Seg_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels5Seg_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels5Seg_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels5Seg_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels10Seg_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels10Seg_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels10Seg_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels10Seg_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nWheels10Seg_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nWheels10Seg_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);



    name = "h_nRPCStations1_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations1_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations1_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations1_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations1_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations1_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations5_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations5_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations5_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations5_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations5_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations5_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations10_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations10_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations10_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations10_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations10_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCStations10_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nRPCWheels1_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels1_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels1_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels1_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels1_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels1_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels5_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels5_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels5_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels5_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels5_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels5_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels10_50hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels10_50hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels10_100hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels10_100hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels10_150hits_clusterMETCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nRPCWheels10_150hits_clusterMETCR[itr_mX][itr_ctau] = new TH1D(name,"",6,0,6);

    
    name = "h_minSegmentDR_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_minSegmentDR_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,0,2);

    name = "h_nMatchedSegments_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegments_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsClusterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsClusterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsClusterStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOtherStations_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOtherStations_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nAlignedSegments_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nAlignedSegments_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStation_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStation_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStation_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOuterStation_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOuterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOuterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOuterStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_segmentAlignmentDeltaPhi_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_segmentAlignmentDeltaPhi_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",70,0,3.5);
    
    name = "h_segmentAlignmentDeltaEta_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_segmentAlignmentDeltaEta_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",70,0,3.5);

    
    name = "h_MB1VetoOutcomes_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_muonVetoOutcomes_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuon_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassSegMuon_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuon_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassOneSegMuon_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomes_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_matchedSegmentTimes_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimes_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMean_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMean_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesClusterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanClusterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesClusterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanClusterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesClusterStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanClusterStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesInnerStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanInnerStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesInnerStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanInnerStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesInnerStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanInnerStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesOuterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB2_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanOuterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesOuterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB3_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanOuterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesOuterStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB4_invertedShowerVetoes_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanOuterStationMB4_invertedShowerVetoes[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);



    name = "h_minSegmentDR_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_minSegmentDR_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,0,2);

    name = "h_nMatchedSegments_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegments_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsClusterStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsClusterStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsClusterStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOtherStations_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOtherStations_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nAlignedSegments_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nAlignedSegments_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStation_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStation_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationMB2_invertedMB1_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB2_invertedMB1_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationMB3_invertedMB1_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB3_invertedMB1_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationMB4_invertedMB1_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationMB4_invertedMB1_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStation_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOuterStation_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOuterStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOuterStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsOuterStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_segmentAlignmentDeltaPhi_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_segmentAlignmentDeltaPhi_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",70,0,3.5);
    
    name = "h_segmentAlignmentDeltaEta_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_segmentAlignmentDeltaEta_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",70,0,3.5);

    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_muonVetoOutcomes_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoOutcomes_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_muonVetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomes_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuon_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassSegMuon_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuon_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassOneSegMuon_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomes_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);
    
    name = "h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_matchedSegmentTimes_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimes_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMean_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMean_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesClusterStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanClusterStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesClusterStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanClusterStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesClusterStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanClusterStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesInnerStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanInnerStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesInnerStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanInnerStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesInnerStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanInnerStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesOuterStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB2_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanOuterStationMB2_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesOuterStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB3_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanOuterStationMB3_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimesOuterStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB4_invertedJetVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedSegmentTimeMeanOuterStationMB4_invertedJetVeto[itr_mX][itr_ctau] = new TH1D(name,"",50,-100,100);

    name = "h_dtRechitClusterJetVetoPt_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dtRechitClusterJetVetoPt[itr_mX][itr_ctau] = new TH1D(name,"",20,0,200);

    name = "h_dtRechitClusterMuonVetoPt_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dtRechitClusterMuonVetoPt[itr_mX][itr_ctau] = new TH1D(name,"",20,0,200);

    name = "h_dtRechitClusterMB1Veto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_dtRechitClusterMB1Veto[itr_mX][itr_ctau] = new TH1D(name,"",60,0,60);


    name = "h_jetNeutralEMEnergyFraction_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_jetNeutralEMEnergyFraction_SR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedEMEnergyFraction_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_jetChargedEMEnergyFraction_SR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

    name = "h_jetNeutralHadronicEnergyFraction_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_jetNeutralHadronicEnergyFraction_SR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedHadronicEnergyFraction_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_jetChargedHadronicEnergyFraction_SR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetNeutralEMEnergyFraction_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadingJetNeutralEMEnergyFraction_SR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetChargedEMEnergyFraction_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadingJetChargedEMEnergyFraction_SR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetNeutralHadronicEnergyFraction_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadingJetNeutralHadronicEnergyFraction_SR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetChargedHadronicEnergyFraction_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadingJetChargedHadronicEnergyFraction_SR[itr_mX][itr_ctau] = new TH1D(name,"",50,0,1);

    name = "h_matchedJetNHF_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedJetNHF[itr_mX][itr_ctau] = new TH1D(name,"",100,0,1);

    name = "h_matchedJetNEMF_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedJetNEMF[itr_mX][itr_ctau] = new TH1D(name,"",100,0,1);

    name = "h_matchedJetCHF_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedJetCHF[itr_mX][itr_ctau] = new TH1D(name,"",100,0,1);

    name = "h_matchedJetCEMF_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedJetCEMF[itr_mX][itr_ctau] = new TH1D(name,"",100,0,1);

    name = "h_matchedJetNConst_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedJetNConst[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);

    name = "h_matchedJetChargedMult_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_matchedJetChargedMult[itr_mX][itr_ctau] = new TH1D(name,"",20,0,20);
  

    name = "h_efficiency_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_efficiency[itr_mX][itr_ctau] = new TH1D(name,"",100,0,100);

    name = "h_efficiency_MB1CR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_efficiency_MB1CR[itr_mX][itr_ctau] = new TH1D(name,"",100,0,100);

    name = "h_decayVertexRadius_noVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_decayVertexRadius_noVeto[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    name = "h_decayVertexRadius_clusterReco_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_decayVertexRadius_clusterReco[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    name = "h_decayVertexRadius_clusterReco_signalRegionGt2_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_decayVertexRadius_clusterReco_signalRegionGt2[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    name = "h_decayVertexRadius_clusterReco_signalRegionEq2_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_decayVertexRadius_clusterReco_signalRegionEq2[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    // name = "h_decayVertexRadius_clusterReco_noMax_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    // h_decayVertexRadius_clusterReco_noMax[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);
    //
    // name = "h_decayVertexRadius_clusterReco_noMB1_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    // h_decayVertexRadius_clusterReco_noMB1[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);
    //
    // name = "h_decayVertexRadius_clusterReco_noMB1_noMax_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    // h_decayVertexRadius_clusterReco_noMB1_noMax[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    name = "h_decayVertexZ_noVeto_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_decayVertexZ_noVeto[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    name = "h_decayVertexZ_clusterReco_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_decayVertexZ_clusterReco[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    name = "h_decayVertexZ_clusterReco_signalRegionGt2_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_decayVertexZ_clusterReco_signalRegionGt2[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    name = "h_decayVertexZ_clusterReco_signalRegionEq2_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_decayVertexZ_clusterReco_signalRegionEq2[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    // name = "h_decayVertexZ_clusterReco_noMax_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    // h_decayVertexZ_clusterReco_noMax[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);
    //
    // name = "h_decayVertexZ_clusterReco_noMB1_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    // h_decayVertexZ_clusterReco_noMB1[itr_mX][itr_ctau] = new TH1D(name,"",80,0,800);

    name = "h_leadJetPt_MB1CR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPt_MB1CR[itr_mX][itr_ctau] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_MB2CR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPt_MB2CR[itr_mX][itr_ctau] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_MB2withMB1CR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPt_MB2withMB1CR[itr_mX][itr_ctau] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_MB1HitsCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPt_MB1HitsCR[itr_mX][itr_ctau] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPt_SR[itr_mX][itr_ctau] = new TH1D(name,"",200,0,2000);
    
    name = "h_leadJetPtMET_MB1CR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPtMET_MB1CR[itr_mX][itr_ctau] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_MB2CR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPtMET_MB2CR[itr_mX][itr_ctau] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_MB2withMB1CR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPtMET_MB2withMB1CR[itr_mX][itr_ctau] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_MB1HitsCR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPtMET_MB1HitsCR[itr_mX][itr_ctau] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetPtMET_SR[itr_mX][itr_ctau] = new TH1D(name,"",500,0,10);

    name = "h_leadJetEta_SR_"+mX[itr_mX]+"_"+ctau[itr_ctau];
    h_leadJetEta_SR[itr_mX][itr_ctau] = new TH1D(name,"",200,-5,5);
    

    name = "h_nCAClusters_"+mX[itr_mX]+ctau[itr_ctau];
    h_nCAClusters[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_sizeCACluster_"+mX[itr_mX]+ctau[itr_ctau];
    h_sizeCACluster[itr_mX][itr_ctau] = new TH1D(name,"",50,0,500);

    name = "h_radiusCACluster_"+mX[itr_mX]+ctau[itr_ctau];
    h_radiusCACluster[itr_mX][itr_ctau] = new TH1D(name,"",50,300,800);

    name = "h_zCACluster_"+mX[itr_mX]+ctau[itr_ctau];
    h_zCACluster[itr_mX][itr_ctau] = new TH1D(name,"",120,-600,600);

    name = "h_phiCACluster_"+mX[itr_mX]+ctau[itr_ctau];
    h_phiCACluster[itr_mX][itr_ctau] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_nDtRechitClusters_"+mX[itr_mX]+ctau[itr_ctau];
    h_nDtRechitClusters[itr_mX][itr_ctau] = new TH1D(name,"",5,0,5);

    name = "h_dtRechitClusterSize_"+mX[itr_mX]+ctau[itr_ctau];
    h_dtRechitClusterSize[itr_mX][itr_ctau] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterR_"+mX[itr_mX]+ctau[itr_ctau];
    h_dtRechitClusterR[itr_mX][itr_ctau] = new TH1D(name,"",50,300,800);

    name = "h_dtRechitClusterZ_"+mX[itr_mX]+ctau[itr_ctau];
    h_dtRechitClusterZ[itr_mX][itr_ctau] = new TH1D(name,"",120,-600,600);

    name = "h_dtRechitClusterPhi_"+mX[itr_mX]+ctau[itr_ctau];
    h_dtRechitClusterPhi[itr_mX][itr_ctau] = new TH1D(name,"",70,-3.5,3.5);


    name = "h_jetOverlap10GeV_"+mX[itr_mX]+ctau[itr_ctau];
    h_jetOverlap10GeV[itr_mX][itr_ctau] = new TH1D(name,"",2,0,2);

    name = "h_jetOverlap20GeV_"+mX[itr_mX]+ctau[itr_ctau];
    h_jetOverlap20GeV[itr_mX][itr_ctau] = new TH1D(name,"",2,0,2);

    name = "h_muonOverlap10GeVLoose_"+mX[itr_mX]+ctau[itr_ctau];
    h_muonOverlap10GeVLoose[itr_mX][itr_ctau] = new TH1D(name,"",2,0,2);

    name = "h_muonOverlap10GeVTight_"+mX[itr_mX]+ctau[itr_ctau];
    h_muonOverlap10GeVTight[itr_mX][itr_ctau] = new TH1D(name,"",2,0,2);


    nPassClusterCR = 0;
    nPassNoVeto_clusterCR = 0;
    nPassFullVeto_clusterCR = 0;
    nPassRPCMatch_clusterCR = 0;
    nPassRPCSpread_clusterCR = 0;
    nPassRPCBx_clusterCR = 0;
    nPassMaxStation_clusterCR = 0;
    nPassLepton_clusterCR = 0;
    nPass50Hits_clusterCR = 0;
    nPass25Hits_clusterCR = 0;
  
    nPassNoVeto = 0;  
    nPassFullVeto_rpcCR = 0;
    nPassRPCCR = 0;
    nPassClusterMET_rpcCR = 0;
    nPassMaxStation_rpcCR = 0;
    nPassFullPlus = 0;
    nPassJetMET_rpcCR = 0;
    nPassLepton_rpcCR = 0;
    nPass50Hits_rpcCR = 0;
    nPass25Hits_rpcCR = 0;
   
    nPassMB2CR = 0;
    nPassMB2CRwithAdjacent = 0;
    nPassMB2CRwithAdjacent0p8 = 0;
    nPassMB2CRwithOther = 0;
    nPassMB2CRwithNHF = 0;
    nPassMB2CRwithNHFnoRPC = 0;
    nPassSRwithNHFnoRPC = 0;
    nPassSignalRegion = 0;

    evtNum = 0;
    totalNum = 0;

    Int_t lifetime = ctau[itr_ctau].Atoi();

    cout << mX[itr_mX] << "_" << ctau[itr_ctau] << endl;
    for(Int_t itr_year = 0; itr_year<years.size(); itr_year++){
      cout << "  " << years[itr_year] << endl;

      evtNum = 0;
      
      TFile *_file;
      
      if(years[itr_year]=="MC_Summer16"){
	if(lifetime%3==0){
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v4/normalized/ggH_HToSSTo4Tau_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau+1]+"_TuneCUETP8M1_13TeV-powheg-pythia8_1pb_weighted.root");
	  _file = TFile::Open(dir+years[itr_year]+"/ggH_HToSSTobbbb_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau+1]+"_TuneCUETP8M1_13TeV-powheg-pythia8_1pb_weighted.root");
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v4/normalized/VBFH_HToSSTo4b_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau+1]+"_TuneCUETP8M1_13TeV-powheg-pythia8_1pb_weighted.root");
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v5/normalized/ggH_HToSSTobbbb_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCUETP8M1_13TeV-powheg-pythia8_1pb_weighted.root");
	}	  
	else{
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v4/normalized/ggH_HToSSTo4Tau_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCUETP8M1_13TeV-powheg-pythia8_1pb_weighted.root");
	  _file = TFile::Open(dir+years[itr_year]+"/ggH_HToSSTobbbb_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCUETP8M1_13TeV-powheg-pythia8_1pb_weighted.root");
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v4/normalized/VBFH_HToSSTo4b_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCUETP8M1_13TeV-powheg-pythia8_1pb_weighted.root");
	  //_file = TFile::Open(dir+years[itr_year]+"/v1/v124/normalized/ggH_HToSSTobbbb_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCUETP8M1_13TeV-powheg-pythia8_35920pb_weighted.root");
	}
      }
      else{
	if(lifetime%3==0){
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v4/normalized/ggH_HToSSTo4Tau_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau+1]+"_TuneCP5_13TeV-powheg-pythia8_1pb_weighted.root");
	  _file = TFile::Open(dir+years[itr_year]+"/ggH_HToSSTobbbb_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau+1]+"_TuneCP5_13TeV-powheg-pythia8_1pb_weighted.root");
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v4/normalized/VBFH_HToSSTo4b_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau+1]+"_TuneCP5_13TeV-powheg-pythia8_1pb_weighted.root");
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v5/normalized/ggH_HToSSTobbbb_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCP5_13TeV-powheg-pythia8_1pb_weighted.root");
	  //_file = TFile::Open(dir+years[itr_year]+"/HV_params_"+mX[itr_mX]+"_ctau_"+ctau[itr_ctau]+"_LLPNTUPLE_v0_filter.root");
	}
	else{
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v4/normalized/ggH_HToSSTo4Tau_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCP5_13TeV-powheg-pythia8_1pb_weighted.root");
	  _file = TFile::Open(dir+years[itr_year]+"/ggH_HToSSTobbbb_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCP5_13TeV-powheg-pythia8_1pb_weighted.root");
	  //_file = TFile::Open(dir+years[itr_year]+"/v3/v4/normalized/VBFH_HToSSTo4b_MH-125_MS-"+mX[itr_mX]+"_ctau-"+ctau[itr_ctau]+"_TuneCP5_13TeV-powheg-pythia8_1pb_weighted.root");
	  //if(years[itr_year]=="2016"){ _file = TFile::Open(dir+years[itr_year]+"/HV_params_"+mX[itr_mX]+"_ctau_"+ctau[itr_ctau]+"_"+years[itr_year]+"_LLPNTUPLE_v3_filter_updatedKnapenCode_benchmarks.root"); }
	  //_file = TFile::Open(dir+years[itr_year]+"/HV_params_"+mX[itr_mX]+"_ctau_"+ctau[itr_ctau]+"_"+years[itr_year]+"_LLPNTUPLE_v3_filter_updatedKnapenCode_benchmarks.root");
	}
      }
      //_file = TFile::Open(dir+"signalPointsGGHWithRPCWithFlags/signal_1000_"+mX[itr_mX]+"_"+ctau[itr_ctau]+".root");
      //_file = TFile::Open(dir+"signalPointsGenFilter/HiddenValleyGenFilter_"+mX[itr_mX]+"_ctau_"+ctau[itr_ctau]+".py_privateMC_102X_LLPNTUPLE_v2_generationGenFilter_forHV_2018_MS.root");
      
      TTreeReader treeReader("MuonSystem",_file);
      //TTreeReader treeReader("ntuples/llp",_file);

      //TTreeReaderValue<bool> Flag2_all(treeReader,"Flag2_all");
      bool* Flag2_all = new bool(false);
      TTreeReaderValue<bool> Flag2_HBHENoiseFilter(treeReader,"Flag2_HBHENoiseFilter");
      TTreeReaderValue<bool> Flag2_HBHEIsoNoiseFilter(treeReader,"Flag2_HBHEIsoNoiseFilter");
      TTreeReaderValue<bool> Flag2_BadPFMuonFilter(treeReader,"Flag2_BadPFMuonFilter");
      TTreeReaderValue<bool> Flag2_globalSuperTightHalo2016Filter(treeReader,"Flag2_globalSuperTightHalo2016Filter");
      TTreeReaderValue<bool> Flag2_EcalDeadCellTriggerPrimitiveFilter(treeReader,"Flag2_EcalDeadCellTriggerPrimitiveFilter");
      TTreeReaderValue<bool> Flag2_ecalBadCalibFilter(treeReader,"Flag2_ecalBadCalibFilter");
      TTreeReaderArray<bool> HLTDecision(treeReader,"HLTDecision");

      TTreeReaderValue<unsigned int> runNum(treeReader,"runNum");
      //TTreeReaderValue<unsigned int> lumiSec(treeReader,"lumiSec");
      //TTreeReaderValue<unsigned int> eventNum(treeReader,"evtNum");
      
      TTreeReaderValue<float> MET(treeReader,"met");
      TTreeReaderValue<float> METphi(treeReader,"metPhi");
     
      //TTreeReaderValue<float> MET(treeReader,"metType1Pt");
      //TTreeReaderValue<float> METphi(treeReader,"metType1Phi");
 
      //TTreeReaderArray<float> gLLP_ctau(treeReader,"gLLP_ctau");
      TTreeReaderArray<float> gLLP_beta(treeReader,"gLLP_beta");
      TTreeReaderArray<float> gLLP_decay_vertex_x(treeReader,"gLLP_decay_vertex_x");
      TTreeReaderArray<float> gLLP_decay_vertex_y(treeReader,"gLLP_decay_vertex_y");
      TTreeReaderArray<float> gLLP_decay_vertex_z(treeReader,"gLLP_decay_vertex_z");

      TTreeReaderValue<int> nDtRechitClusters(treeReader,"nDtRechitClusters");
      TTreeReaderArray<float> dtRechitClusterX(treeReader,"dtRechitClusterX");
      TTreeReaderArray<float> dtRechitClusterY(treeReader,"dtRechitClusterY");
      TTreeReaderArray<float> dtRechitClusterZ(treeReader,"dtRechitClusterZ");
      TTreeReaderArray<float> dtRechitClusterPhi(treeReader,"dtRechitClusterPhi");
      TTreeReaderArray<float> dtRechitClusterEta(treeReader,"dtRechitClusterEta");
      TTreeReaderArray<float> dtRechitClusterTime(treeReader,"dtRechitClusterTime");
      TTreeReaderArray<float> dtRechitClusterXSpread(treeReader,"dtRechitClusterXSpread");
      TTreeReaderArray<float> dtRechitClusterYSpread(treeReader,"dtRechitClusterYSpread");
      TTreeReaderArray<float> dtRechitClusterZSpread(treeReader,"dtRechitClusterZSpread");
      TTreeReaderArray<float> dtRechitClusterPhiSpread(treeReader,"dtRechitClusterPhiSpread");
      TTreeReaderArray<float> dtRechitClusterEtaSpread(treeReader,"dtRechitClusterEtaSpread");
      TTreeReaderArray<float> dtRechitClusterTimeSpread(treeReader,"dtRechitClusterTimeSpread");
      TTreeReaderArray<float> dtRechitClusterJetVetoPt(treeReader,"dtRechitClusterJetVetoPt");
      TTreeReaderArray<float> dtRechitClusterMuonVetoPt(treeReader,"dtRechitClusterMuonVetoPt");
      TTreeReaderArray<int> dtRechitClusterSize(treeReader,"dtRechitClusterSize");
      TTreeReaderArray<int> dtRechitClusterNSegmentStation1(treeReader,"dtRechitClusterNSegmentStation1");
      TTreeReaderArray<int> dtRechitClusterNSegmentStation2(treeReader,"dtRechitClusterNSegmentStation2");
      TTreeReaderArray<int> dtRechitClusterNSegmentStation3(treeReader,"dtRechitClusterNSegmentStation3");
      TTreeReaderArray<int> dtRechitClusterNSegmentStation4(treeReader,"dtRechitClusterNSegmentStation4");
      TTreeReaderArray<int> dtRechitClusterMaxStation(treeReader,"dtRechitClusterMaxStation");
      TTreeReaderArray<float> dtRechitClusterMaxStationRatio(treeReader,"dtRechitClusterMaxStationRatio");
      TTreeReaderArray<int> dtRechitClusterNStation(treeReader,"dtRechitClusterNStation");
      TTreeReaderArray<int> dtRechitClusterMaxChamber(treeReader,"dtRechitClusterMaxChamber");
      TTreeReaderArray<float> dtRechitClusterMaxChamberRatio(treeReader,"dtRechitClusterMaxChamberRatio");
      TTreeReaderArray<int> dtRechitClusterNChamber(treeReader,"dtRechitClusterNChamber");
      TTreeReaderArray<float> dtRechitClusterMajorAxis(treeReader,"dtRechitClusterMajorAxis");
      TTreeReaderArray<float> dtRechitClusterMinorAxis(treeReader,"dtRechitClusterMinorAxis");
      
      /*TTreeReaderValue<int> nCscRechitClusters(treeReader,"nCscRechitClusters");
      TTreeReaderArray<int> cscRechitClusterSize(treeReader,"cscRechitClusterSize");
      TTreeReaderArray<float> cscRechitClusterMe11Ratio(treeReader,"cscRechitClusterMe11Ratio");
      TTreeReaderArray<float> cscRechitClusterMe12Ratio(treeReader,"cscRechitClusterMe12Ratio");
      */
      TTreeReaderValue<int> nDtRechits(treeReader,"nDtRechits");
      TTreeReaderArray<float> dtRechitX(treeReader,"dtRechitX");
      TTreeReaderArray<float> dtRechitY(treeReader,"dtRechitY");
      TTreeReaderArray<float> dtRechitZ(treeReader,"dtRechitZ");
      TTreeReaderArray<float> dtRechitEta(treeReader,"dtRechitEta");
      TTreeReaderArray<float> dtRechitPhi(treeReader,"dtRechitPhi");
      TTreeReaderArray<int> dtRechitStation(treeReader,"dtRechitStation");
      TTreeReaderArray<int> dtRechitWheel(treeReader,"dtRechitWheel");

      /*TTreeReaderValue<int> nDtSeg(treeReader,"nDtSeg");
      TTreeReaderArray<float> dtSegX(treeReader,"dtSegX");
      TTreeReaderArray<float> dtSegY(treeReader,"dtSegY");
      TTreeReaderArray<float> dtSegZ(treeReader,"dtSegZ");
      TTreeReaderArray<float> dtSegEta(treeReader,"dtSegEta");
      TTreeReaderArray<float> dtSegPhi(treeReader,"dtSegPhi");
      TTreeReaderArray<int> dtSegStation(treeReader,"dtSegStation");
      TTreeReaderArray<int> dtSegWheel(treeReader,"dtSegWheel");
      TTreeReaderArray<float> dtSegTime(treeReader,"dtSegTime");
      TTreeReaderArray<float> dtSegTimeError(treeReader,"dtSegTimeError");
      */
      TTreeReaderValue<int> nJets(treeReader,"nJets");
      TTreeReaderArray<float> jetPt(treeReader,"jetPt");
      TTreeReaderArray<float> jetEta(treeReader,"jetEta");
      TTreeReaderArray<float> jetPhi(treeReader,"jetPhi");
      TTreeReaderArray<bool> jetTightPassId(treeReader,"jetTightPassId");
      //TTreeReaderArray<bool> jetTightPassId(treeReader,"jetPassIDTight");
      //TTreeReaderArray<float> jetRechitT(treeReader,"jetRechitT");
      //TTreeReaderArray<float> jetRechitT_rms(treeReader,"jetRechitT_rms");
      TTreeReaderArray<float> jetElectronEnergyFraction(treeReader,"jetElectronEnergyFraction");
      TTreeReaderArray<float> jetPhotonEnergyFraction(treeReader,"jetPhotonEnergyFraction");
      //TTreeReaderArray<float> jetElectronEnergyFraction(treeReader,"jetChargedEMEnergyFraction");
      //TTreeReaderArray<float> jetPhotonEnergyFraction(treeReader,"jetNeutralEMEnergyFraction");
      TTreeReaderArray<float> jetNeutralHadronEnergyFraction(treeReader,"jetNeutralHadronEnergyFraction");
      TTreeReaderArray<float> jetChargedHadronEnergyFraction(treeReader,"jetChargedHadronEnergyFraction");
      //TTreeReaderArray<float> jetMuonEnergyFraction(treeReader,"jetMuonEnergyFraction");
      /*TTreeReaderArray<int> jetNPFCands(treeReader,"jetNPFCands");
      TTreeReaderArray<int> jetElectronMultiplicity(treeReader,"jetElectronMultiplicity");
      TTreeReaderArray<int> jetMuonMultiplicity(treeReader,"jetMuonMultiplicity");
      TTreeReaderArray<int> jetPhotonMultiplicity(treeReader,"jetPhotonMultiplicity");
      TTreeReaderArray<int> jetNeutralHadronMultiplicity(treeReader,"jetNeutralHadronMultiplicity");
      TTreeReaderArray<int> jetChargedHadronMultiplicity(treeReader,"jetChargedHadronMultiplicity");
      */
      TTreeReaderValue<int> nRPCRechits(treeReader,"nRpc");
      TTreeReaderArray<float> RPCRechitX(treeReader,"rpcX");
      TTreeReaderArray<float> RPCRechitY(treeReader,"rpcY");
      TTreeReaderArray<float> RPCRechitZ(treeReader,"rpcZ");
      TTreeReaderArray<float> RPCRechitPhi(treeReader,"rpcPhi");
      TTreeReaderArray<float> RPCRechitEta(treeReader,"rpcEta");
      TTreeReaderArray<int> RPCRechitBx(treeReader,"rpcBx");
      
      /*TTreeReaderValue<int> nHORechits(treeReader,"nHORechits");
      TTreeReaderArray<float> hoRechitX(treeReader,"hoRechit_X");
      TTreeReaderArray<float> hoRechitY(treeReader,"hoRechit_Y");
      TTreeReaderArray<float> hoRechitZ(treeReader,"hoRechit_Z");
      TTreeReaderArray<float> hoRechitEta(treeReader,"hoRechit_Eta");
      TTreeReaderArray<float> hoRechitPhi(treeReader,"hoRechit_Phi");
      TTreeReaderArray<float> hoRechitE(treeReader,"hoRechit_E");
      TTreeReaderArray<float> hoRechitT(treeReader,"hoRechit_T");
      
      TTreeReaderValue<int> nMuons(treeReader,"nMuons");
      TTreeReaderArray<float> muonPhi(treeReader,"muonPhi");
      TTreeReaderArray<float> muonEta(treeReader,"muonEta");
      TTreeReaderArray<float> muonPt(treeReader,"muonPt");
      TTreeReaderArray<float> muon_ip3dSignificance(treeReader,"muon_ip3dSignificance");
      TTreeReaderArray<float> muon_chargedIso(treeReader,"muon_chargedIso");
      TTreeReaderArray<float> muon_photonIso(treeReader,"muon_photonIso");
      TTreeReaderArray<float> muon_neutralHadIso(treeReader,"muon_neutralHadIso");
      TTreeReaderArray<float> muon_pileupIso(treeReader,"muon_pileupIso");
      TTreeReaderArray<bool> muonIsLoose(treeReader,"muonIsLoose");
      
      TTreeReaderValue<int> nElectrons(treeReader,"nElectrons");
      TTreeReaderArray<float> elePhi(treeReader,"elePhi");
      TTreeReaderArray<float> eleEta(treeReader,"eleEta");
      TTreeReaderArray<float> elePt(treeReader,"elePt");
      TTreeReaderArray<float> ele_d0(treeReader,"ele_d0");
      TTreeReaderArray<float> ele_dZ(treeReader,"ele_dZ");
      TTreeReaderArray<bool> eleIsLoose(treeReader,"ele_passCutBasedIDLoose");
      */
      TTreeReaderValue<int> nLeptons(treeReader,"nLeptons");
      TTreeReaderArray<int> lepPdgId(treeReader,"lepPdgId");
      TTreeReaderArray<float> lepPhi(treeReader,"lepPhi");
      TTreeReaderArray<float> lepEta(treeReader,"lepEta");
      TTreeReaderArray<float> lepPt(treeReader,"lepPt");
      TTreeReaderArray<bool> lepPassId(treeReader,"lepPassId");
      
      TTreeReaderValue<float> higgsPtWeight(treeReader,"higgsPtWeight");
      TTreeReaderValue<float> sf_facScaleUp(treeReader,"sf_facScaleUp");
      TTreeReaderValue<float> sf_facScaleDown(treeReader,"sf_facScaleDown");
      TTreeReaderValue<float> sf_renScaleUp(treeReader,"sf_renScaleUp");
      TTreeReaderValue<float> sf_renScaleDown(treeReader,"sf_renScaleDown");
      TTreeReaderValue<float> sf_facRenScaleUp(treeReader,"sf_facRenScaleUp");
      TTreeReaderValue<float> sf_facRenScaleDown(treeReader,"sf_facRenScaleDown");
      TTreeReaderValue<float> pileupWeight(treeReader,"pileupWeight");
      TTreeReaderValue<float> pileupWeightUp(treeReader,"pileupWeightUp");
      TTreeReaderValue<float> pileupWeightDown(treeReader,"pileupWeightDown");
      TTreeReaderValue<float> metJESUp(treeReader,"metJESUp");
      TTreeReaderValue<float> metJESDown(treeReader,"metJESDown");
      
      
      TH1F *hEvents;
      if(lifetime%3==0){
	hEvents = (TH1F*)_file->Get("NEvents"+mX[itr_mX]+ctau[itr_ctau+1]);
      }
      else{
	hEvents = (TH1F*)_file->Get("NEvents"+mX[itr_mX]+ctau[itr_ctau]);
      }
      
      _ofile->cd();
      totalNum += treeReader.GetEntries(1);
      std::vector<bool> gLLP_plotted; 
      std::vector<bool> gLLP_plotted_id_sr1; 
      std::vector<bool> gLLP_plotted_id_sr2; 
      //weight = 48.58*1000*0.01*lumi[years[itr_year]]/treeReader.GetEntries(1);
      //weight = 48.58*1000*0.01*lumi[itr_year]/1E6;
      //weight = 48.58*1000*0.01*137/500000;
      //weight = 0.1845*1000*1.00*137/treeReader.GetEntries(1);
      while(treeReader.Next()){
	//weight = 48.58*1000*0.01*lumi[years[itr_year]]/1E6;
	//weight = 48.58*1000*0.01*lumi[years[itr_year]]/treeReader.GetEntries(1);
	//weight = 48.58*1000*0.01/treeReader.GetEntries(1);
	weight = 48.58*1000*0.01*lumi[years[itr_year]]/hEvents->GetEntries();
	//weight = 3.782*1000*0.01*lumi[years[itr_year]]/treeReader.GetEntries(1);
       	if(lifetime%3==0){
	  //weight = 100.;
	  //weight = weight*exp((gLLP_ctau[0]+gLLP_ctau[1])*(10./lifetime - 100./lifetime));
	  //weight = weight*lumi[itr_year]/treeReader.GetEntries(1);
	  
	  decay1 = sqrt(pow(gLLP_decay_vertex_x[0],2)+pow(gLLP_decay_vertex_y[0],2)+pow(gLLP_decay_vertex_z[0],2));
	  decay2 = sqrt(pow(gLLP_decay_vertex_x[1],2)+pow(gLLP_decay_vertex_y[1],2)+pow(gLLP_decay_vertex_z[1],2));
	  ctau1 = decay1 / (gLLP_beta[0]*(1.0/sqrt(1-gLLP_beta[0]*gLLP_beta[0])));
	  ctau2 = decay2 / (gLLP_beta[1]*(1.0/sqrt(1-gLLP_beta[1]*gLLP_beta[1])));
	  Int_t upLifetime = ctau[itr_ctau+1].Atoi();
	  weight = weight*pow(upLifetime,2)/pow(lifetime,2)*exp((ctau1+ctau2)*(10./upLifetime - 10./lifetime));
	  
	}
	weight = weight*(*higgsPtWeight)*(*pileupWeight);

	gLLP_plotted_id_sr1.clear();
	gLLP_plotted_id_sr2.clear();
	gLLP_plotted.clear();
	if(evtNum%100000==0){ cout << evtNum << " of " << treeReader.GetEntries(1) << endl; }
	passFullVeto_clusterCR = false;
	passRPCMatch_clusterCR = false;
	passRPCSpread_clusterCR = false;
	passRPCBx_clusterCR = false;
	passMaxStation_clusterCR = false;
	passLepton_clusterCR = false;
	pass50Hits_clusterCR = false;
	pass25Hits_clusterCR = false;
	passFullVeto_rpcCR = false;
	passRPCCR = false;
	passClusterMET_rpcCR = false;
	passMaxStation_rpcCR = false;
	passJetMET_rpcCR = false;
	passLepton_rpcCR = false;
	pass50Hits_rpcCR = false;
	pass25Hits_rpcCR = false;
	passSignalRegion = false;
	passMB2CR = false;
	passMB2Cluster = false;
	passMB2CRwithAdjacent = false;
	passMB2CRwithAdjacent0p8 = false;
	passMB2CRwithOther = false;
	passMB2CRwithNHF = false;
	passMB2CRwithNHFnoRPC = false;
	passSRwithNHFnoRPC = false;
	nWheels1=0;
	nWheels25=0;
	nWheels50=0;
	nStations1=0;
	nStations25=0;
	nStations50=0;
	hitStation1=0;
	hitStation2=0;
	hitStation3=0;
	hitStation4=0;
	hitWheelm2=0;
	hitWheelm1=0;
	hitWheel0=0;
	hitWheel1=0;
	hitWheel2=0;
	nWheels1Seg=0;
	nWheels5Seg=0;
	nWheels10Seg=0;
	nStations1Seg=0;
	nStations5Seg=0;
	nStations10Seg=0;
	segStation1=0;
	segStation2=0;
	segStation3=0;
	segStation4=0;
	segWheelm2=0;
	segWheelm1=0;
	segWheel0=0;
	segWheel1=0;
	segWheel2=0;
	for(int i=0; i<4; i++){
	  for(int j=0; j<5; j++){
	    segChambers[i][j]=0;
	  }
	}

	nRPCWheels1=0;
	nRPCWheels5=0;
	nRPCWheels10=0;
	nRPCStations1=0;
	nRPCStations5=0;
	nRPCStations10=0;
	hitRPCStation1=0;
	hitRPCStation2=0;
	hitRPCStation3=0;
	hitRPCStation4=0;
	hitRPCWheelm2=0;
	hitRPCWheelm1=0;
	hitRPCWheel0=0;
	hitRPCWheel1=0;
	hitRPCWheel2=0;
	
	maxClusterSize=0;
	maxClusterSizeSR=0;
	maxClusterSizeMB1CR=0;

	passMET = false;
	passOneJet = false;
	passNHFJet = true;
	passJetMET = false;
	passStations25 = false;
	passWheels25 = false;
	passJetVeto = false;
	passJetTightIdVeto = false;
	passMuonVeto = false;
	passMB1Veto = false;
	passMaxStation = false;
	passClusterMET = false;
	passRPCMatch = false;
	passRPCSpread = false;
	passRPCBx = false;
	passNoVetoCluster = false;
	passClusterSize = false;
	passClusterSizeA = false;
	passClusterMETA = false;
	passClusterMETB = false;
	passClusterMETC = false;
	passClusterSizeB = false;
	passClusterSizeC = false;
	passAdjacentMB1 = false;
	passAdjacent0p8MB1 = false;
	passOtherStations = false;
	passCscCluster = false;
	passLargeCscCluster = false;
	passCscClusterME11Veto = false;
	passLargeCscClusterME11Veto = false;

	passMB1CR = false;
	passJetVetoMB1CR = false;
	passMuonVetoMB1CR = false;
	passMuonVetoLooseMB1CR = false;
	passRpcMatchMB1CR = false;
	passClusterMETMB1CR = false;
	passClusterSizeMB1CR = false;

	nClustersVeto_dPhiJetMET = 0;
	clusterEta.clear();
	clusterPhi.clear();
	clusterSize.clear();
	clusterSizeTotal = 0;

	nMB1MatchClusterAdjacentPlus = 0;
	nMB1MatchClusterAdjacentMinus = 0;
	nMB1MatchPi2AdjacentPlus = 0;
	nMB1MatchPi2 = 0;
	nMB1MatchPi2AdjacentMinus = 0;
	nMB1MatchClusterAdjacent0p8Plus = 0;
	nMB1MatchClusterAdjacent0p8Minus = 0;
	nMB1MatchPi2Adjacent0p8Plus = 0;
	nMB1MatchPi2Adjacent0p8Minus = 0;
	nRB1MatchCluster = 0;
	
	//dtPoints.clear();
	//clustersCA.clear();
	//constituents.clear();
	nCAClusters = 0;
	CAclusterJetVeto = false;
	CAclusterMuonVeto = false;
	CAclusterMB1Match = 0;
	CAclusterAdjacentMB1MatchPlus = 0;
	CAclusterAdjacentMB1MatchMinus = 0;
	CAclusterRPCMatch = 0;
	CAclusterdPhiMET = 0.0;
	CAclusterSize = 0;
	CAclusterStation = 0;
	CAclusterWheel = 0;
	
	if(rand->Uniform()<0.5){ pmRand = -1; }
	else{ pmRand = 1; }
	
	*Flag2_all = *Flag2_HBHENoiseFilter && *Flag2_HBHEIsoNoiseFilter && *Flag2_BadPFMuonFilter && *Flag2_globalSuperTightHalo2016Filter && *Flag2_EcalDeadCellTriggerPrimitiveFilter;
	HLT = false;
	if(years[itr_year]=="MC_Summer16" || years[itr_year]=="2016"){
	  if(HLTDecision[310] || HLTDecision[467]){
	    HLT = true;
	  }
	}
	else{
	  if(HLTDecision[310] || HLTDecision[467] || HLTDecision[703] || HLTDecision[717] || HLTDecision[710] || HLTDecision[709]){
	    HLT = true;
	  }
	  *Flag2_all = *Flag2_all && *Flag2_ecalBadCalibFilter;
	}
	/*
	for(Int_t itr_lep = 0; itr_lep<*nMuons; itr_lep++){
	  if(muonIsLoose[itr_lep] && muonPt[itr_lep]>25. && fabs(muonEta[itr_lep])<2.4){
	    if(fabs(muon_ip3dSignificance[itr_lep]) < 4){
	      if((muon_chargedIso[itr_lep] + fmax(0.0,  muon_photonIso[itr_lep] + muon_neutralHadIso[itr_lep] - 0.5*muon_pileupIso[itr_lep])) / muonPt[itr_lep] < 0.25){
		HLT = false;
	      }
	    }
	  }
	}
	for(Int_t itr_lep = 0; itr_lep<*nElectrons; itr_lep++){
	  if(eleIsLoose[itr_lep] && elePt[itr_lep]>35. && fabs(muonEta[itr_lep])<2.5 && fabs(ele_d0[itr_lep])<0.1 && fabs(ele_dZ[itr_lep])<0.2){
	    HLT = false;
	  }
	  }*/
	//if(*nLeptons>0){ HLT = false; }
	
	if(*MET > 200 && *Flag2_all & HLT){
	  passMET = true;
	  for(Int_t itr_llp=0; itr_llp<2; itr_llp++){
	      float decayR = sqrt(pow(gLLP_decay_vertex_x[itr_llp],2)+pow(gLLP_decay_vertex_y[itr_llp],2));
	      if(decayR > 300 && decayR < 800 && gLLP_decay_vertex_z[itr_llp] < 650.){
		  h_decayVertexRadius_noVeto[itr_mX][itr_ctau]->Fill(sqrt(pow(gLLP_decay_vertex_x[itr_llp],2)+pow(gLLP_decay_vertex_y[itr_llp],2)));
		  h_decayVertexZ_noVeto[itr_mX][itr_ctau]->Fill(fabs(gLLP_decay_vertex_z[itr_llp]));
	      }
	      // float phi = TMath::ATan2(gLLP_decay_vertex_y[itr_llp],gLLP_decay_vertex_x[itr_llp]);
	      // float theta = TMath::ACos(gLLP_decay_vertex_z[itr_llp]/TMath::Sqrt(gLLP_decay_vertex_x[itr_llp]*gLLP_decay_vertex_x[itr_llp]+gLLP_decay_vertex_y[itr_llp]*gLLP_decay_vertex_y[itr_llp]+gLLP_decay_vertex_z[itr_llp]*gLLP_decay_vertex_z[itr_llp]));
	      // float eta = -TMath::Log(TMath::Tan(theta/2.));
	      gLLP_plotted.push_back(false); 
	      gLLP_plotted_id_sr1.push_back(false); 
	      gLLP_plotted_id_sr2.push_back(false); 
	      // gLLP_plotted_id_noMax.push_back(false); 
	      // gLLP_plotted_id_noMB1.push_back(false); 
	      // gLLP_plotted_id_noMB1_noMax.push_back(false); 
	  }
	  dPhi_min = 999.;
	  dPhiClusterMET = 0.0;
	  dPhiClusterMET_max = 0.0;
	  chargedHadFraction_mindPhi = -1.0;
	  chargedEMFraction_mindPhi = -1.0;
	  neutralHadFraction_mindPhi = -1.0;
	  neutralEMFraction_mindPhi = -1.0;
	  if(*nDtRechitClusters>0){
	    nPassNoVeto+=1;
	    for(Int_t itr_clust = 0; itr_clust<*nDtRechitClusters; itr_clust++){
	      dPhiClusterMET = dtRechitClusterPhi[itr_clust] - *METphi;
	      if(dPhiClusterMET > TMath::Pi()){ dPhiClusterMET -= 2*TMath::Pi(); }
	      if(dPhiClusterMET < -1.0*TMath::Pi()){ dPhiClusterMET += 2*TMath::Pi(); }
	      if(fabs(dPhiClusterMET)>dPhiClusterMET_max){ dPhiClusterMET_max=fabs(dPhiClusterMET); }
	      if(dtRechitClusterSize[itr_clust]>maxClusterSize){ maxClusterSize = dtRechitClusterSize[itr_clust]; }
	      if(nStations25<3 && nWheels25<3 && dtRechitClusterSize[itr_clust]>50){
		h_dtRechitClusterSize[itr_mX][itr_ctau]->Fill(dtRechitClusterSize[itr_clust],weight);
		h_dtRechitClusterR[itr_mX][itr_ctau]->Fill(sqrt(pow(dtRechitClusterX[itr_clust],2)+pow(dtRechitClusterY[itr_clust],2)),weight);
		h_dtRechitClusterZ[itr_mX][itr_ctau]->Fill(dtRechitClusterZ[itr_clust],weight);
		h_dtRechitClusterPhi[itr_mX][itr_ctau]->Fill(dtRechitClusterPhi[itr_clust],weight);
	      }
	    }
	  }
	  if(fabs(dPhiClusterMET)<1.0){ nPassClusterCR+=1; }

	  /*if(*nCscRechitClusters>0){ passCscCluster = true; }
	  for(Int_t itr_csc=0; itr_csc<*nCscRechitClusters; itr_csc++){
	    if(cscRechitClusterMe11Ratio[itr_csc]==0 && cscRechitClusterMe12Ratio[itr_csc]==0){
	      passCscClusterME11Veto = true;
	      if(cscRechitClusterSize[itr_csc]>130.){
		passLargeCscClusterME11Veto = true;
	      }
	    }
	    if(cscRechitClusterSize[itr_csc]>130.){
	      passLargeCscCluster = true;
	    }
	    }*/

	  if(fabs(jetEta[0])<=2.4 && jetTightPassId[0]){ passNHFJet = true; }
	  for(Int_t itr_jet = 0; itr_jet<*nJets; itr_jet++){
	    if(fabs(jetEta[itr_jet])<3.0 && jetPt[itr_jet]>30.0){
	      passOneJet = true;
	      //if(jetNeutralHadronEnergyFraction[itr_jet]>=0.9 && fabs(jetEta[itr_jet])<2.4 && jetPt[itr_jet]>50.0){ passNHFJet = false; }
	      dPhi_tmp = jetPhi[itr_jet] - *METphi;
	      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	      if(fabs(dPhi_tmp) < dPhi_min){ 
		dPhi_min = fabs(dPhi_tmp); 
		chargedHadFraction_mindPhi = jetChargedHadronEnergyFraction[itr_jet];
		chargedEMFraction_mindPhi = jetElectronEnergyFraction[itr_jet];
		neutralHadFraction_mindPhi = jetNeutralHadronEnergyFraction[itr_jet];
		neutralEMFraction_mindPhi = jetPhotonEnergyFraction[itr_jet];
	      }
	    }
	  }
	  if(fabs(dPhi_min)>0.6 && passOneJet){ 
	    passJetMET = true; 
	    h_nDtRechitClusters_dPhiJetMET[itr_mX][itr_ctau]->Fill(*nDtRechitClusters,weight);
	  }
	  

	  for(Int_t itr_dt = 0; itr_dt<*nDtRechits; itr_dt++){
	    //dtPoints.push_back( fastjet::PseudoJet( dtRechitX[itr_dt], dtRechitY[itr_dt], dtRechitZ[itr_dt], 0) );
	    //dtR = sqrt(pow(dtRechitX[itr_dt],2)+pow(dtRechitY[itr_dt],2));
	    //if(dtR>400. && dtR<480.){ hitStation1+=1; }
	    //else if(dtR>485. && dtR<560.){ hitStation2+=1; }
	    //else if(dtR>590. && dtR<650.){ hitStation3+=1; }
	    //else if(dtR>690. && dtR<800.){ hitStation4+=1; }
	    if(dtRechitStation[itr_dt]==1){ hitStation1+=1; }
	    else if(dtRechitStation[itr_dt]==2){ hitStation2+=1; }
	    else if(dtRechitStation[itr_dt]==3){ hitStation3+=1; }
	    else if(dtRechitStation[itr_dt]==4){ hitStation4+=1; }
	    //if(dtRechitZ[itr_dt]>-661. && dtRechitZ[itr_dt]<-395.){ hitWheelm2+=1; }
	    //else if(dtRechitZ[itr_dt]>-395. && dtRechitZ[itr_dt]<-127.){ hitWheelm1+=1; }
	    //else if(fabs(dtRechitZ[itr_dt]<127.)){ hitWheel0+=1; }
	    //else if(dtRechitZ[itr_dt]<395.){ hitWheel1+=1; }
	    //else if(dtRechitZ[itr_dt]<661.){ hitWheel2+=1; }
	    if(dtRechitWheel[itr_dt]==-2){ hitWheelm2+=1; }
	    else if(dtRechitWheel[itr_dt]==-1){ hitWheelm1+=1; }
	    else if(dtRechitWheel[itr_dt]==0){ hitWheel0+=1; }
	    else if(dtRechitWheel[itr_dt]==1){ hitWheel1+=1; }
	    else if(dtRechitWheel[itr_dt]==2){ hitWheel2+=1; }
	  }
	  /*for(Int_t itr_dt = 0; itr_dt<*nDtSeg; itr_dt++){
	    if(dtSegStation[itr_dt]==1){ segStation1+=1; }
	    else if(dtSegStation[itr_dt]==2){ segStation2+=1; }
	    else if(dtSegStation[itr_dt]==3){ segStation3+=1; }
	    else if(dtSegStation[itr_dt]==4){ segStation4+=1; }
	    if(dtSegWheel[itr_dt]==-2){ segWheelm2+=1; }
	    else if(dtSegWheel[itr_dt]==-1){ segWheelm1+=1; }
	    else if(dtSegWheel[itr_dt]==0){ segWheel0+=1; }
	    else if(dtSegWheel[itr_dt]==1){ segWheel1+=1; }
	    else if(dtSegWheel[itr_dt]==2){ segWheel2+=1; }
	    if(dtSegStation[itr_dt]>=1&&dtSegStation[itr_dt]<=4&&dtSegWheel[itr_dt]>=-2&&dtSegWheel[itr_dt]<=2){
	      segChambers[dtSegStation[itr_dt]-1][dtSegWheel[itr_dt]+2]+=1;
	    }
	    }*/
	  
	  if(hitStation1>25){ nStations25+=1; }
	  if(hitStation2>25){ nStations25+=1; }
	  if(hitStation3>25){ nStations25+=1; }
	  if(hitStation4>25){ nStations25+=1; }
	  if(hitWheelm2>25){ nWheels25+=1; }
	  if(hitWheelm1>25){ nWheels25+=1; }
	  if(hitWheel0>25){ nWheels25+=1; }
	  if(hitWheel1>25){ nWheels25+=1; }
	  if(hitWheel2>25){ nWheels25+=1; }
	  
	  if(segStation1>0){ nStations1Seg+=1; }
	  if(segStation2>0){ nStations1Seg+=1; }
	  if(segStation3>0){ nStations1Seg+=1; }
	  if(segStation4>0){ nStations1Seg+=1; }
	  if(segWheelm2>0){ nWheels1Seg+=1; }
	  if(segWheelm1>0){ nWheels1Seg+=1; }
	  if(segWheel0>0){ nWheels1Seg+=1; }
	  if(segWheel1>0){ nWheels1Seg+=1; }
	  if(segWheel2>0){ nWheels1Seg+=1; }
	  if(segStation1>=5){ nStations5Seg+=1; }
	  if(segStation2>=5){ nStations5Seg+=1; }
	  if(segStation3>=5){ nStations5Seg+=1; }
	  if(segStation4>=5){ nStations5Seg+=1; }
	  if(segWheelm2>=5){ nWheels5Seg+=1; }
	  if(segWheelm1>=5){ nWheels5Seg+=1; }
	  if(segWheel0>=5){ nWheels5Seg+=1; }
	  if(segWheel1>=5){ nWheels5Seg+=1; }
	  if(segWheel2>=5){ nWheels5Seg+=1; }
	  if(segStation1>=10){ nStations10Seg+=1; }
	  if(segStation2>=10){ nStations10Seg+=1; }
	  if(segStation3>=10){ nStations10Seg+=1; }
	  if(segStation4>=10){ nStations10Seg+=1; }
	  if(segWheelm2>=10){ nWheels10Seg+=1; }
	  if(segWheelm1>=10){ nWheels10Seg+=1; }
	  if(segWheel0>=10){ nWheels10Seg+=1; }
	  if(segWheel1>=10){ nWheels10Seg+=1; }
	  if(segWheel2>=10){ nWheels10Seg+=1; }
	  
	  if(nStations25<3){ passStations25=true; }
	  if(nWheels25<3){ passWheels25=true; }
	  if(*nDtRechitClusters>0){ passNoVetoCluster=true; }

	  /*if(nStations25<3 && nWheels25<3){
	    fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);
	    fastjet::ClusterSequence cs(dtPoints, jet_def);
	    clustersCA = cs.inclusive_jets();
	    for(int i=0; i<clustersCA.size(); i++){
	      constituents.clear();
	      constituents = clustersCA[i].constituents();
	      if(constituents.size() > 50){ 
		nCAClusters+=1;
		h_sizeCACluster[itr_mX][itr_ctau]->Fill(constituents.size(),weight);
		h_radiusCACluster[itr_mX][itr_ctau]->Fill(sqrt(pow(clustersCA[i].px(),2)+pow(clustersCA[i].py(),2)),weight);
		h_zCACluster[itr_mX][itr_ctau]->Fill(clustersCA[i].pz(),weight);
		h_phiCACluster[itr_mX][itr_ctau]->Fill(clustersCA[i].phi_std(),weight);
		
		fill(CAclusterStationHits, CAclusterStationHits+4, 0);
		fill(CAclusterWheelHits, CAclusterWheelHits+5, 0);
		for(int j=0; j<constituents.size(); j++){
		  float hitR = sqrt(pow(constituents[j].px(),2)+pow(constituents[j].py(),2));
		  if(hitR > 402. && hitR < 449.){ CAclusterStationHits[0]+=1; }
		  else if(hitR > 490.5  && hitR < 533.5){ CAclusterStationHits[1]+=1; }
		  else if(hitR > 597.5 && hitR < 636.){ CAclusterStationHits[2]+=1; }
		  else if(hitR > 700. && hitR < 738.){ CAclusterStationHits[3]+=1; }
		  
		  if(constituents[j].pz() < -395.){ CAclusterWheelHits[0]+=1; }
		  else if(constituents[j].pz() < -127.){ CAclusterWheelHits[1]+=1; }
		  else if(constituents[j].pz() < 127.){ CAclusterWheelHits[2]+=1; }
		  else if(constituents[j].pz() < 395.){ CAclusterWheelHits[3]+=1; }
		  else if(constituents[j].pz() < 661.){ CAclusterWheelHits[4]+=1; }
		}
		CAclusterStation = distance(CAclusterStationHits, max_element(CAclusterStationHits, CAclusterStationHits+4)) + 1;
		CAclusterWheel = distance(CAclusterWheelHits, max_element(CAclusterWheelHits, CAclusterWheelHits+5)) - 2;
		
		CAclusterJetVeto = true;
		for(Int_t itr_jet = 0; itr_jet<*nJets; itr_jet++){
		  if(fabs(jetEta[itr_jet])<3.0){
		    dPhi_tmp = jetPhi[itr_jet] - clustersCA[i].phi_std();
		    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		    if(sqrt(pow(dPhi_tmp,2)+pow(jetEta[itr_jet] - clustersCA[i].eta(),2))<0.4){
		      if(jetPt[itr_jet]>=20.0){
			CAclusterJetVeto = false;
			break;
		      }
		    }
		  }
		}
		if(CAclusterJetVeto){
		  
		  CAclusterMuonVeto = true;
		  for(Int_t itr_lep=0; itr_lep<*nLeptons; itr_lep++){
		    if(abs(lepPdgId[itr_lep])==13){
		      dPhi_tmp = lepPhi[itr_lep] - clustersCA[i].phi_std();
		      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		      if(sqrt(pow(dPhi_tmp,2)+pow(lepEta[itr_lep]-clustersCA[i].eta(),2))<0.4){
			if(lepPt[itr_lep]>=10.0){
			  CAclusterMuonVeto=false;
			  break;
			}
		      }
		    }
		  }
		  if(CAclusterMuonVeto){
		    
		    CAclusterMB1Match = 0;
		    CAclusterAdjacentMB1MatchPlus = 0;
		    CAclusterAdjacentMB1MatchMinus = 0;
		    for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
		      if(dtRechitStation[itr_dt]==1){
			dPhi_tmp = clustersCA[i].phi_std() - dtRechitPhi[itr_dt];
			if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			if(sqrt(pow(dPhi_tmp,2)+pow(clustersCA[i].eta() - dtRechitEta[itr_dt],2))<0.4){
			  CAclusterMB1Match += 1;
			}
			if(dtRechitWheel[itr_dt]==CAclusterWheel+1 && fabs(dPhi_tmp)<TMath::Pi()/4.0){
			  CAclusterAdjacentMB1MatchPlus += 1;
			}
			else if(dtRechitWheel[itr_dt]==CAclusterWheel-1 && fabs(dPhi_tmp)<TMath::Pi()/4.0){
			  CAclusterAdjacentMB1MatchMinus += 1;
			}
		      }
		    }
		    if(CAclusterMB1Match<=1 && CAclusterAdjacentMB1MatchPlus<=8 && CAclusterAdjacentMB1MatchMinus<=8){
		      
		      CAclusterRPCMatch = 0;
		      for(Int_t itr_rpc=0; itr_rpc<*nRPCRechits; itr_rpc++){
			dPhi_tmp = clustersCA[i].phi_std() - RPCRechitPhi[itr_rpc];
			if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			if(fabs(RPCRechitZ[itr_rpc] - clustersCA[i].pz())<5. && fabs(dPhi_tmp)<0.4){
			  CAclusterRPCMatch += 1;
			}
		      }
		      if(CAclusterRPCMatch>0){
			CAclusterdPhiMET = *METphi - clustersCA[i].phi_std();
			if(CAclusterdPhiMET > TMath::Pi()){ CAclusterdPhiMET -= 2*TMath::Pi(); }
			if(CAclusterdPhiMET < -1.0*TMath::Pi()){ CAclusterdPhiMET += 2*TMath::Pi(); }
			
			if(constituents.size()>CAclusterSize){ CAclusterSize = constituents.size(); }
			
		      }
		    }
		  }
		}
	      }
	    }
	    h_nCAClusters[itr_mX][itr_ctau]->Fill(nCAClusters,weight);
	    h_nDtRechitClusters[itr_mX][itr_ctau]->Fill(*nDtRechitClusters,weight);
	    
	    if(passJetMET && passNHFJet){
	      if(fabs(CAclusterdPhiMET)<1.0 && CAclusterSize>=100){
		if(CAclusterStation==2){ SRyieldMB2_CA+=1; }
		if(CAclusterStation>2){ SRyield_CA+=1; }
	      }
	    }
	  }
	    */
	  for(Int_t itr_clust=0; itr_clust<*nDtRechitClusters; itr_clust++){
	    if(dtRechitClusterSize[itr_clust]>50){  
		TVector3 dtRechitClusterVec = TVector3();
		dtRechitClusterVec.SetPtEtaPhi(1,dtRechitClusterEta[itr_clust],dtRechitClusterPhi[itr_clust]);
		for(Int_t itr_llp=0; itr_llp<2; itr_llp++){
		    float decayR = sqrt(pow(gLLP_decay_vertex_x[itr_llp],2)+pow(gLLP_decay_vertex_y[itr_llp],2));
		    if(decayR > 300 && decayR < 800 && fabs(gLLP_decay_vertex_z[itr_llp]) < 650. && !gLLP_plotted[itr_llp]){
			TVector3 llpDecay = TVector3(gLLP_decay_vertex_x[itr_llp],gLLP_decay_vertex_y[itr_llp],gLLP_decay_vertex_z[itr_llp]);
			// if (llpDecay.DeltaR(dtRechitClusterVec) < 0.5)
			{
			    h_decayVertexRadius_clusterReco[itr_mX][itr_ctau]->Fill(sqrt(pow(gLLP_decay_vertex_x[itr_llp],2)+pow(gLLP_decay_vertex_y[itr_llp],2)));
			    h_decayVertexZ_clusterReco[itr_mX][itr_ctau]->Fill(fabs(gLLP_decay_vertex_z[itr_llp]));
			    gLLP_plotted[itr_llp] = true;
			}
		    }
		}
	      passMuon=false;
	      passMuonLoose=false;
	      passMuon_alt=false;
	      passJet=false;
	      passJetTightId=false;
	      rpcBx.clear();
	      rpcSpread = 99;
	      rpcMedian = 99;
	      dPhiClusterRPC = -0.1;
	      dZClusterRPC = -1.;
	      //nStations25 = 0;
	      //nStations50 = 0;
	      hoMatchedEnergy = 0.;
	      hitsMB1 = 0;
	      
	      dPhiClusterMET = dtRechitClusterPhi[itr_clust] - *METphi;
	      if(dPhiClusterMET > TMath::Pi()){ dPhiClusterMET -= 2*TMath::Pi(); }
	      if(dPhiClusterMET < -1.0*TMath::Pi()){ dPhiClusterMET += 2*TMath::Pi(); }
	      
	      h_dtRechitClusterJetVetoPt[itr_mX][itr_ctau]->Fill(dtRechitClusterJetVetoPt[itr_clust],weight);
	      h_dtRechitClusterMuonVetoPt[itr_mX][itr_ctau]->Fill(dtRechitClusterMuonVetoPt[itr_clust],weight);
	      if(dtRechitClusterJetVetoPt[itr_clust]<10.){ passJet = true; } 
	      passJetTightId = true;
	      overlapJet10GeV = false;
	      overlapJet20GeV = false;
	      for(Int_t itr_jet = 0; itr_jet<*nJets; itr_jet++){
		//if((jetChargedHadronEnergyFraction[itr_jet]+jetElectronEnergyFraction[itr_jet])>0.1){
		dPhi_tmp = dtRechitClusterPhi[itr_clust] - jetPhi[itr_jet] + pmRand*TMath::Pi()/2.0;
		//dPhi_tmp = jetPhi[itr_jet] - dtRechitClusterPhi[itr_clust];
		if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust]-jetEta[itr_jet],2))<0.4){
		  if(jetPt[itr_jet]>10.){
		    overlapJet10GeV = true;
		    if(jetPt[itr_jet]>20.){
		      //passJet = false;
		      //break;
		      overlapJet20GeV = true;
		    }
		  }
		}
		dPhi_tmp = jetPhi[itr_jet] - dtRechitClusterPhi[itr_clust];
		if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust]-jetEta[itr_jet],2))<0.4){
		  if(jetPt[itr_jet]>20.){
		    h_matchedJetNHF[itr_mX][itr_ctau]->Fill(jetPhotonEnergyFraction[itr_jet],weight);
		    h_matchedJetNEMF[itr_mX][itr_ctau]->Fill(jetElectronEnergyFraction[itr_jet],weight);
		    h_matchedJetCHF[itr_mX][itr_ctau]->Fill(jetChargedHadronEnergyFraction[itr_jet],weight);
		    h_matchedJetCEMF[itr_mX][itr_ctau]->Fill(jetNeutralHadronEnergyFraction[itr_jet],weight);
		    //h_matchedJetNConst[itr_mX][itr_ctau]->Fill(jetMuonMultiplicity[itr_jet]+jetPhotonMultiplicity[itr_jet]+jetElectronMultiplicity[itr_jet]+jetNeutralHadronMultiplicity[itr_jet]+jetChargedHadronMultiplicity[itr_year],weight);
		    //h_matchedJetChargedMult[itr_mX][itr_ctau]->Fill(jetChargedHadronMultiplicity[itr_jet]+jetElectronMultiplicity[itr_jet],weight);
		    //if(jetTightPassId[itr_jet]){
		    if(jetChargedHadronEnergyFraction[itr_jet]>0 || jetElectronEnergyFraction[itr_jet]>0){
		      passJetTightId = false;
		    }
		  }
		}
	      }
	      if(overlapJet10GeV){ h_jetOverlap10GeV[itr_mX][itr_ctau]->Fill(1); }
	      else{ h_jetOverlap10GeV[itr_mX][itr_ctau]->Fill(0.0); }
	      if(overlapJet20GeV){ h_jetOverlap20GeV[itr_mX][itr_ctau]->Fill(1); }
	      else{ h_jetOverlap20GeV[itr_mX][itr_ctau]->Fill(0.0); }
	      //passJet = passJetTightId;
	      /*passJet = true;
	      for(Int_t itr_jet = 0; itr_jet<*nJets; itr_jet++){
		if((jetChargedHadronEnergyFraction[itr_jet]+jetElectronEnergyFraction[itr_jet])>0.1){
		  dPhi_tmp = jetPhi[itr_jet] - dtRechitClusterPhi[itr_clust];
		  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		  if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust]-jetEta[itr_jet],2))<0.4){
		    if(jetPt[itr_jet]>20.){
		      passJet = false;
		      break;
		    }
		  }
		}
		}
	      */
	      passMuonLoose=true;
	      matchedLooseIDMuonPt = 0.0;
	      overlapMuon10GeVLoose=false;
	      overlapMuon10GeVTight=false;
	      for(Int_t itr_lep=0; itr_lep<*nLeptons; itr_lep++){
		if(abs(lepPdgId[itr_lep])==13){
		  dPhi_tmp = lepPhi[itr_lep] - dtRechitClusterPhi[itr_clust];
		  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		  if(sqrt(pow(dPhi_tmp,2)+pow(lepEta[itr_lep]-dtRechitClusterEta[itr_clust],2))<0.4){
		    if(lepPt[itr_lep]>matchedLooseIDMuonPt){
		      matchedLooseIDMuonPt = lepPt[itr_lep];
		    }
		    if(lepPt[itr_lep]>=10.0){
		      passMuonLoose=false;
		      //break;
		    }
		  }
		  dPhi_tmp += pmRand*TMath::Pi()/2.0;
		  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		  if(sqrt(pow(dPhi_tmp,2)+pow(lepEta[itr_lep]-dtRechitClusterEta[itr_clust],2))<0.4){
		    if(lepPt[itr_lep]>=10.0){
		      overlapMuon10GeVLoose=true;
		      if(lepPassId[itr_lep]){
			overlapMuon10GeVTight=true;
		      }
		    }
		  }
		}
	      }
	      if(overlapMuon10GeVLoose){ h_muonOverlap10GeVLoose[itr_mX][itr_ctau]->Fill(1); }
	      else{ h_muonOverlap10GeVLoose[itr_mX][itr_ctau]->Fill(0.0); }
	      if(overlapMuon10GeVTight){ h_muonOverlap10GeVTight[itr_mX][itr_ctau]->Fill(1); }
	      else{ h_muonOverlap10GeVTight[itr_mX][itr_ctau]->Fill(0.0); }
	      /*
	      for(Int_t itr_mu=0; itr_mu<*nMuons; itr_mu++){
		if(muonIsLoose[itr_mu]){
		  dPhi_tmp = muonPhi[itr_mu] - dtRechitClusterPhi[itr_clust];
		  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		  if(sqrt(pow(dPhi_tmp,2)+pow(muonEta[itr_mu]-dtRechitClusterEta[itr_clust],2))<0.4){
		    if(muonPt[itr_mu]>10.0){
		      passMuonLoose=false;
		      break;
		    }
		  }
		}
	      }
	      */

	      //if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ passMuon = true; } 
	      if(passMuonLoose){ passMuon = true; }
	      else{ passMuon = false; }
	      //if(*nLeptons==0){ passMuon_alt = true; }
	    
	      /*passMB1 = false;
	      for(Int_t itr_ho = 0; itr_ho<*nHORechits; itr_ho++){
		dPhi_tmp = hoRechitPhi[itr_ho] - dtRechitClusterPhi[itr_clust];
		if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		if(sqrt(pow(dPhi_tmp,2)+pow(hoRechitEta[itr_ho]-dtRechitClusterEta[itr_clust],2))<0.5){
		  hoMatchedEnergy+=hoRechitE[itr_ho];
		}
	      }
	      if(hoMatchedEnergy>40.0){ passMB1 = true; }
	      */
	      
	      //cout << "doing rpc" << endl;
	      for(Int_t itr_rpc=0; itr_rpc<*nRPCRechits; itr_rpc++){
		dPhi_tmp = RPCRechitPhi[itr_rpc] - dtRechitClusterPhi[itr_clust];
		if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		if(fabs(dPhi_tmp)<dPhiClusterRPC || dPhiClusterRPC==-0.1){ dPhiClusterRPC=fabs(dPhi_tmp); }
		if(fabs(RPCRechitZ[itr_rpc] - dtRechitClusterZ[itr_clust])<dZClusterRPC || dZClusterRPC==-1.){ dZClusterRPC=fabs(RPCRechitZ[itr_rpc]-dtRechitClusterZ[itr_rpc]); }
		if(fabs(RPCRechitZ[itr_rpc] - dtRechitClusterZ[itr_clust])<5. && fabs(dPhi_tmp)<0.4){
		  rpcBx.push_back(RPCRechitBx[itr_rpc]);
		}
		/*if(itr_clust==0 && *nRPCRechits<500){
		  rpcStation=getRPCLayer(RPCRechitX[itr_rpc],RPCRechitY[itr_rpc]);
		  rpcWheel=getWheel(RPCRechitZ[itr_rpc]);
		  if(abs(rpcWheel)<=2){
		    if(rpcStation==1 || rpcStation==2){ hitRPCStation1+=1; }
		    else if(rpcStation==3 || rpcStation==4){ hitRPCStation2+=1; }
		    else if(rpcStation==5){ hitRPCStation3+=1; }
		    else if(rpcStation==6){ hitRPCStation4+=1; }
		    if(rpcWheel==-2){ hitRPCWheelm2+=1; }
		    else if(rpcWheel==-1){ hitRPCWheelm1+=1; }
		    else if(rpcWheel==0){ hitRPCWheel0+=1; }
		    else if(rpcWheel==1){ hitRPCWheel1+=1; }
		    else if(rpcWheel==2){ hitRPCWheel2+=1; }
		  }
		  }*/
	      }
	      if(!rpcBx.empty()){
		rpcSpread = max_element(rpcBx.begin(), rpcBx.end()) - min_element(rpcBx.begin(), rpcBx.end());
		if(rpcBx.size()%2 == 0){ rpcMedian = float(rpcBx[rpcBx.size()/2 - 1] + rpcBx[rpcBx.size()/2]) / 2.0; }
		else{ rpcMedian = rpcBx[rpcBx.size()/2]; }
	      }

	      nMB1MatchCluster=0;
	      nMB2MatchCluster=0;
	      nMB3MatchCluster=0;
	      nMB4MatchCluster=0;
	      passMB1 = true;
	      for(Int_t itr_dt = 0; itr_dt<*nDtRechits; itr_dt++){
		dPhi_tmp = dtRechitPhi[itr_dt] - dtRechitClusterPhi[itr_clust];
		if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitEta[itr_dt]-dtRechitClusterEta[itr_clust],2))<0.4){
		  if(dtRechitStation[itr_dt]==1){ nMB1MatchCluster+=1; }
		  if(dtRechitStation[itr_dt]==2){ nMB2MatchCluster+=1; }
		  if(dtRechitStation[itr_dt]==3){ nMB3MatchCluster+=1; }
		  if(dtRechitStation[itr_dt]==4){ nMB4MatchCluster+=1; }
		    //break;
		}
		/*if(itr_clust==0 && *nDtRechits<750){
		  dtStation=getStation(dtRechitX[itr_dt],dtRechitY[itr_dt]);
		  dtWheel=getWheel(dtRechitZ[itr_dt]);
		  if(dtStation==1){ hitStation1+=1; }
		  else if(dtStation==2){ hitStation2+=1; }
		  else if(dtStation==3){ hitStation3+=1; }
		  else if(dtStation==4){ hitStation4+=1; }
		  if(dtWheel==-2){ hitWheelm2+=1; }
		  else if(dtWheel==-1){ hitWheelm1+=1; }
		  else if(dtWheel==0){ hitWheel0+=1; }
		  else if(dtWheel==1){ hitWheel1+=1; }
		  else if(dtWheel==2){ hitWheel2+=1; }
		  }*/
	      }
	      /*if(itr_clust==0){
		if(hitStation1>0){
		  nStations1+=1;
		  if(hitStation1>25){
		    nStations25+=1;
		    if(hitStation1>50){
		      nStations50+=1;
		    }
		  }
		}
		if(hitStation2>0){
		  nStations1+=1;
		  if(hitStation2>25){
		    nStations25+=1;
		    if(hitStation2>50){
		      nStations50+=1;
		    }
		  }
		}
		if(hitStation3>0){
		  nStations1+=1;
		  if(hitStation3>25){
		    nStations25+=1;
		    if(hitStation3>50){
		      nStations50+=1;
		    }
		  }
		}
		if(hitStation4>0){
		  nStations1+=1;
		  if(hitStation4>25){
		    nStations25+=1;
		    if(hitStation4>50){
		      nStations50+=1;
		    }
		  }
		}
		if(hitWheel1>0){
		  nWheels1+=1;
		  if(hitWheel1>25){
		    nWheels25+=1;
		    if(hitWheel1>50){
		      nWheels50+=1;
		    }
		  }
		}
		if(hitWheel2>0){
		  nWheels1+=1;
		  if(hitWheel2>25){
		    nWheels25+=1;
		    if(hitWheel2>50){
		      nWheels50+=1;
		    }
		  }
		}
		if(hitWheel0>0){
		  nWheels1+=1;
		  if(hitWheel0>25){
		    nWheels25+=1;
		    if(hitWheel0>50){
		      nWheels50+=1;
		    }
		  }
		}
		if(hitWheelm1>0){
		  nWheels1+=1;
		  if(hitWheelm1>25){
		    nWheels25+=1;
		    if(hitWheelm1>50){
		      nWheels50+=1;
		    }
		  }
		}
		if(hitWheelm2>0){
		  nWheels1+=1;
		  if(hitWheelm2>25){
		    nWheels25+=1;
		    if(hitWheelm2>50){
		      nWheels50+=1;
		    }
		  }
		}
	      }
	      if(*nDtRechits>=750){
		nStations1=4;
		nStations25=4;
		nStations50=4;
		nWheels1=5;
		nWheels25=5;
		nWheels50=5;
		}*/
	      hitsMB1 = nMB1MatchCluster;
	      //if(hitsMB1>1 || dtRechitClusterNSegmentStation1[itr_clust]>0){ passMB1 = false; } 
	      h_dtRechitClusterMB1Veto[itr_mX][itr_ctau]->Fill(hitsMB1,weight);
	      
	      nMB1SegMatchCluster=0;
	      nMB2SegMatchCluster=0;
	      nMB3SegMatchCluster=0;
	      nMB4SegMatchCluster=0;
	      minMB1SegmentDR = 4.9;
	      minSegmentDR = 4.9;
	      meanSegTime = 0;
	      meanMB1SegTime = 0;
	      meanMB2SegTime = 0;
	      meanMB3SegTime = 0;
	      meanMB4SegTime = 0;
	      /*for(Int_t itr_dt=0; itr_dt<*nDtSeg; itr_dt++){
		dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtSegPhi[itr_dt];
		if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtSegEta[itr_dt],2))<0.4){
		  if(dtSegStation[itr_dt]==1){ 
		    nMB1SegMatchCluster+=1; 
		    meanMB1SegTime+=dtSegTime[itr_dt];
		  }
		  if(dtSegStation[itr_dt]==2){ 
		    nMB2SegMatchCluster+=1;
		    meanMB2SegTime+=dtSegTime[itr_dt];
		  }
		  if(dtSegStation[itr_dt]==3){ 
		    nMB3SegMatchCluster+=1;
		    meanMB3SegTime+=dtSegTime[itr_dt]; 
		  }
		  if(dtSegStation[itr_dt]==4){ 
		    nMB4SegMatchCluster+=1;
		    meanMB4SegTime+=dtSegTime[itr_dt];
		  }
		  meanSegTime+=dtSegTime[itr_dt];
		}
		if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtSegEta[itr_dt],2))<minSegmentDR){
		  minSegmentDR = sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtSegEta[itr_dt],2));
		}
		if(dtSegStation[itr_dt]==1){
		  if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtSegEta[itr_dt],2))<minMB1SegmentDR){
		    minMB1SegmentDR = sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtSegEta[itr_dt],2));
		  }
		} 
		}*/
	      if(nMB1SegMatchCluster>0){ meanMB1SegTime = meanMB1SegTime/nMB1SegMatchCluster; }
	      else{ meanMB1SegTime = -999.; }
	      if(nMB2SegMatchCluster>0){ meanMB2SegTime = meanMB2SegTime/nMB2SegMatchCluster; }
	      else{ meanMB2SegTime = -999.; }
	      if(nMB3SegMatchCluster>0){ meanMB3SegTime = meanMB3SegTime/nMB3SegMatchCluster; }
	      else{ meanMB3SegTime = -999.; }
	      if(nMB4SegMatchCluster>0){ meanMB4SegTime = meanMB4SegTime/nMB4SegMatchCluster; }
	      else{ meanMB4SegTime = -999.; }
	      if(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster>0){
		meanSegTime = meanSegTime/(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster);
	      }
	      else{ meanSegTime = -999.; }

	      //if(hitsMB1>1 || nMB1SegMatchCluster>0){ passMB1 = false; }
	      if(hitsMB1>1 || dtRechitClusterNSegmentStation1[itr_clust]>0){ passMB1 = false; }

	      /*passMuon = true;
	      if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster>0){ passMuon = false; }
	      if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster>0){ passMuon = false; }
	      if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster>0){ passMuon = false; }
	      */
	      
	      if(passStations25 && passWheels25 && passJetMET && passJet && !rpcBx.empty() && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJet){
		if(passMB1){ 
		  h_MB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		}
		else if(nMB1MatchCluster>1 && nMB1SegMatchCluster>0){ 
		  h_MB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		}
		else if(nMB1MatchCluster>1 && nMB1SegMatchCluster==0){ 
		  h_MB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		}
		else if(nMB1MatchCluster<=1 && nMB1SegMatchCluster>0){
		  h_MB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		}
		else{
		  h_MB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		}
	      }
	      if(passStations25 && passWheels25 && passJetMET && passJet && passMB1 && !rpcBx.empty() && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJet){
		h_minSegmentDR_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(minSegmentDR);
		h_nMatchedSegments_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster);
		h_matchedSegmentTimeMean_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanSegTime);
		nMB1SegMatchClusterAdjacent0p8Plus=0;
		nMB1SegMatchClusterAdjacent0p8Minus=0;
		/*for(Int_t itr_dt=0; itr_dt<*nDtSeg; itr_dt++){
		  if(dtSegStation[itr_dt]==dtRechitClusterMaxStation[itr_clust]-1){
		    dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtSegPhi[itr_dt];
		    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		    if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
		      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1SegMatchClusterAdjacent0p8Plus+=1; }
		      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1SegMatchClusterAdjacent0p8Minus+=1; }
		    }
		  }
		  dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtSegPhi[itr_dt];
		  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		  if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtSegEta[itr_dt],2))<0.4){
		    nAlignedSegments = 0;
		    for(Int_t j=0; j<*nDtSeg; j++){
		      if(j!=itr_dt){
			dPhi_tmp = dtSegPhi[j] - dtSegPhi[itr_dt];
			if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			h_segmentAlignmentDeltaPhi_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(fabs(dPhi_tmp));
			h_segmentAlignmentDeltaEta_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(fabs(dtSegEta[itr_dt]-dtSegEta[j]));
			if(fabs(dPhi_tmp) < 0.28 && fabs(dtSegEta[itr_dt] - dtSegEta[itr_dt]) < 0.28){
			  nAlignedSegments+=1;
			}
		      }
		    }
		    h_nAlignedSegments_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nAlignedSegments);
		  }
		  }*/
		nMB1MatchClusterAdjacent0p8Plus = 0;
		nMB1MatchClusterAdjacent0p8Minus = 0;
		for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
		  if(dtRechitStation[itr_dt]==1){
		    dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
		    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		    if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
		      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacent0p8Plus+=1; }
		      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacent0p8Minus+=1; }
		    }
		  }
		}
		if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
		  h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); } 
		}
		else if((nMB1MatchClusterAdjacent0p8Plus>=8 && nMB1SegMatchClusterAdjacent0p8Plus>0) || (nMB1MatchClusterAdjacent0p8Minus>=8 && nMB1SegMatchClusterAdjacent0p8Minus>0)){
		  h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); } 
		}
		else if((nMB1MatchClusterAdjacent0p8Plus>=8 && nMB1SegMatchClusterAdjacent0p8Plus==0) || (nMB1MatchClusterAdjacent0p8Minus>=8 && nMB1SegMatchClusterAdjacent0p8Minus==0)){
		  h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); } 
		}
		else if((nMB1MatchClusterAdjacent0p8Plus<8 && nMB1SegMatchClusterAdjacent0p8Plus>0) || (nMB1MatchClusterAdjacent0p8Minus<8 && nMB1SegMatchClusterAdjacent0p8Minus>0)){
		  h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); } 
		}
		else{
		  h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2){ 
		  h_nMatchedSegmentsOtherStations_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster); 
		  h_nMatchedSegmentsClusterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster);
		  h_nMatchedSegmentsOuterStation_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsOuterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsInnerStationMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster);
		  h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		  h_matchedSegmentTimeMeanClusterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanMB2SegTime);
		  h_matchedSegmentTimeMeanInnerStationMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanMB1SegTime);
		  h_matchedSegmentTimeMeanOuterStationMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanMB3SegTime);
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		    h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster); 
		    h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0);
		  }
		  if(passMuonLoose){ 
		    h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster); 
		    h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0);
		  }
		  if(nMB1SegMatchCluster!=1){ 
		    h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster); 
		    h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0);
		  }
		  if(nMB1SegMatchCluster==0){
		    h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0);
		  }
		  h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0);
		}
		if(dtRechitClusterMaxStation[itr_clust]==3){ 
		  h_nMatchedSegmentsClusterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsOtherStations_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB4SegMatchCluster); 
		  h_nMatchedSegmentsInnerStation_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster);
		  h_nMatchedSegmentsOuterStation_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB4SegMatchCluster);
		  h_nMatchedSegmentsInnerStationMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster);
		  h_nMatchedSegmentsOuterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB4SegMatchCluster);
		  h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		  h_matchedSegmentTimeMeanClusterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanMB3SegTime);
		  h_matchedSegmentTimeMeanInnerStationMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanMB2SegTime);
		  h_matchedSegmentTimeMeanOuterStationMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanMB4SegTime);
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		    h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster); 
		    h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0);
		  }
		  if(passMuonLoose){ 
		    h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster); 
		    h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0);
		  }
		  if(nMB2SegMatchCluster!=1){ 
		    h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster); 
		    h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0);
		  }
		  if(nMB2SegMatchCluster==0){
		    h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0);
		  }
		  h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0);
		}
		if(dtRechitClusterMaxStation[itr_clust]==4){ 
		  h_nMatchedSegmentsClusterStationMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB4SegMatchCluster);
		  h_nMatchedSegmentsOtherStations_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster); 
		  h_nMatchedSegmentsInnerStation_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsInnerStationMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		  if(nMB3SegMatchCluster==0 && nMB2SegMatchCluster==0){ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0); }
		  else if(nMB3SegMatchCluster>0 && nMB2SegMatchCluster==0){ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0); }
		  else if(nMB3SegMatchCluster==0 && nMB2SegMatchCluster>0){ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0); }
		  else if(nMB2SegMatchCluster>0 && nMB3SegMatchCluster>0){ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0); }
		  else{ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0); }
		  h_matchedSegmentTimeMeanClusterStationMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanMB4SegTime);
		  h_matchedSegmentTimeMeanInnerStationMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(meanMB3SegTime);
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		    h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster); 
		    h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(0.0);
		  }
		  if(passMuonLoose){ 
		    h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster); 
		    h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(1.0);
		  }
		  if(nMB3SegMatchCluster!=1){ 
		    h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster); 
		    h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(2.0);
		  }
		  if(nMB3SegMatchCluster==0){
		    h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(3.0);
		  }
		  h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_mX][itr_ctau]->Fill(4.0);
		}
	      }

	      if(passStations25 && passWheels25 && passJetMET && passJet && !rpcBx.empty() && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJet){
		if(passMB1){ 
		  h_MB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		}
		else if(nMB1MatchCluster>1 && nMB1SegMatchCluster>0){ 
		  h_MB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		}
		else if(nMB1MatchCluster>1 && nMB1SegMatchCluster==0){ 
		  h_MB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		}
		else if(nMB1MatchCluster<=1 && nMB1SegMatchCluster>0){ 
		  h_MB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		}
		else{
		  h_MB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0);
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		    h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  }
		  if(passMuonLoose){
		    h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); 
		    if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		    if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  }
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		}
		if(passMB1){
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_nMatchedSegmentsInnerStationMB2_invertedMB1_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster); }
		  else if(dtRechitClusterMaxStation[itr_clust]==3){ h_nMatchedSegmentsInnerStationMB3_invertedMB1_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster); }
		  else if(dtRechitClusterMaxStation[itr_clust]==4){ h_nMatchedSegmentsInnerStationMB4_invertedMB1_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster); }
		}
	      }
	      if(passStations25 && passWheels25 && passJetMET && passJet && passMB1 && !rpcBx.empty() && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJet){
		h_minSegmentDR_invertedJetVeto[itr_mX][itr_ctau]->Fill(minSegmentDR);
		h_nMatchedSegments_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster);
		h_matchedSegmentTimeMean_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanSegTime);
		nMB1SegMatchClusterAdjacent0p8Plus=0;
		nMB1SegMatchClusterAdjacent0p8Minus=0;
		/*for(Int_t itr_dt=0; itr_dt<*nDtSeg; itr_dt++){
		  if(dtSegStation[itr_dt]==dtRechitClusterMaxStation[itr_clust]-1){
		    dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtSegPhi[itr_dt];
		    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		    if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
		      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1SegMatchClusterAdjacent0p8Plus+=1; }
		      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1SegMatchClusterAdjacent0p8Minus+=1; }
		    }
		  }
		  dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtSegPhi[itr_dt];
		  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		  if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtSegEta[itr_dt],2))<0.4){
		    nAlignedSegments = 0;
		    for(Int_t j=0; j<*nDtSeg; j++){
		      if(j!=itr_dt){
			dPhi_tmp = dtSegPhi[j] - dtSegPhi[itr_dt];
			if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			h_segmentAlignmentDeltaPhi_invertedJetVeto[itr_mX][itr_ctau]->Fill(fabs(dPhi_tmp));
			h_segmentAlignmentDeltaEta_invertedJetVeto[itr_mX][itr_ctau]->Fill(fabs(dtSegEta[itr_dt]-dtSegEta[j]));
			if(fabs(dPhi_tmp) < 0.28 && fabs(dtSegEta[itr_dt] - dtSegEta[itr_dt]) < 0.28){
			  nAlignedSegments+=1;
			}
		      }
		    }
		    h_nAlignedSegments_invertedJetVeto[itr_mX][itr_ctau]->Fill(nAlignedSegments);
		  }
		  }*/
		nMB1MatchClusterAdjacent0p8Plus = 0;
		nMB1MatchClusterAdjacent0p8Minus = 0;
		for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
		  if(dtRechitStation[itr_dt]==1){
		    dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
		    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		    if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
		      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacent0p8Plus+=1; }
		      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacent0p8Minus+=1; }
		    }
		  }
		}
		if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
		  h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); } 
		}
		else if((nMB1MatchClusterAdjacent0p8Plus>=8 && nMB1SegMatchClusterAdjacent0p8Plus>0) || (nMB1MatchClusterAdjacent0p8Minus>=8 && nMB1SegMatchClusterAdjacent0p8Minus>0)){
		  h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); } 
		}
		else if((nMB1MatchClusterAdjacent0p8Plus>=8 && nMB1SegMatchClusterAdjacent0p8Plus==0) || (nMB1MatchClusterAdjacent0p8Minus>=8 && nMB1SegMatchClusterAdjacent0p8Minus==0)){
		  h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); } 
		}
		else if((nMB1MatchClusterAdjacent0p8Plus<8 && nMB1SegMatchClusterAdjacent0p8Plus>0) || (nMB1MatchClusterAdjacent0p8Minus<8 && nMB1SegMatchClusterAdjacent0p8Minus>0)){
		  h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); } 
		}
		else{
		  h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2){ 
		  h_nMatchedSegmentsClusterStationMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster);
		  h_nMatchedSegmentsOtherStations_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster); 
		  h_nMatchedSegmentsOuterStation_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsOuterStationMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsInnerStationMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster);
		  h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		  h_matchedSegmentTimeMeanClusterStationMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanMB2SegTime);
		  h_matchedSegmentTimeMeanInnerStationMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanMB1SegTime);
		  h_matchedSegmentTimeMeanOuterStationMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanMB3SegTime);
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		    h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster); 
		    h_muonVetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0);
		  }
		  if(passMuonLoose){ 
		    h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster); 
		    h_muonVetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0);
		  }
		  if(nMB1SegMatchCluster!=1){ 
		    h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster); 
		    h_muonVetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0);
		  }
		  if(nMB1SegMatchCluster==0){
		    h_muonVetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0);
		  }		
		  h_muonVetoOutcomesMB2_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0);
		}
		if(dtRechitClusterMaxStation[itr_clust]==3){ 
		  h_nMatchedSegmentsClusterStationMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsOtherStations_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB4SegMatchCluster); 
		  h_nMatchedSegmentsInnerStation_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster);
		  h_nMatchedSegmentsOuterStation_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB4SegMatchCluster);
		  h_nMatchedSegmentsInnerStationMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster);
		  h_nMatchedSegmentsOuterStationMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB4SegMatchCluster);
		  h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		  h_matchedSegmentTimeMeanClusterStationMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanMB3SegTime);
		  h_matchedSegmentTimeMeanInnerStationMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanMB2SegTime);
		  h_matchedSegmentTimeMeanOuterStationMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanMB4SegTime);
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		    h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster); 
		    h_muonVetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0);
		  }
		  if(passMuonLoose){ 
		    h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster); 
		    h_muonVetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0);
		  }
		  if(nMB2SegMatchCluster!=1){ 
		    h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB2SegMatchCluster); 
		    h_muonVetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0);
		  }
		  if(nMB2SegMatchCluster==0){
		    h_muonVetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0);
		  }
		  h_muonVetoOutcomesMB3_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0);
		}
		if(dtRechitClusterMaxStation[itr_clust]==4){
		  h_nMatchedSegmentsClusterStationMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB4SegMatchCluster);
		  h_nMatchedSegmentsOtherStations_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster); 
		  h_nMatchedSegmentsInnerStation_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsInnerStationMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster);
		  h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		  if(nMB3SegMatchCluster==0 && nMB2SegMatchCluster==0){ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0); }
		  else if(nMB3SegMatchCluster>0 && nMB2SegMatchCluster==0){ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0); }
		  else if(nMB3SegMatchCluster==0 && nMB2SegMatchCluster>0){ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0); }
		  else if(nMB2SegMatchCluster>0 && nMB3SegMatchCluster>0){ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0); }
		  else{ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0); }
		  h_matchedSegmentTimeMeanClusterStationMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanMB4SegTime);
		  h_matchedSegmentTimeMeanInnerStationMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(meanMB3SegTime);
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		    h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster); 
		    h_muonVetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(0.0);
		  }
		  if(passMuonLoose){ 
		    h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster); 
		    h_muonVetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(1.0);
		  }
		  if(nMB3SegMatchCluster!=1){ 
		    h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(nMB3SegMatchCluster); 
		    h_muonVetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(2.0);
		  }
		  if(nMB3SegMatchCluster==0){
		    h_muonVetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(3.0);
		  }
		  h_muonVetoOutcomesMB4_invertedJetVeto[itr_mX][itr_ctau]->Fill(4.0);
		}
	      }
	      
	      if(fabs(dPhi_min)>0.6 && dtRechitClusterJetVetoPt[itr_clust]<20.0 && dtRechitClusterMuonVetoPt[itr_clust]<10.0){
		nClustersVeto_dPhiJetMET+=1;
		h_dtRechitClusterNSegmentStation2_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterNSegmentStation2[itr_clust],weight);
		h_dtRechitClusterNSegmentStation3_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterNSegmentStation3[itr_clust],weight);
		h_dtRechitClusterNSegmentStation4_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterNSegmentStation4[itr_clust],weight);
		h_dtRechitClusterMaxStationRatio_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterMaxStationRatio[itr_clust],weight);
		h_dtRechitClusterMaxStation_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterMaxStation[itr_clust],weight);
		h_dtRechitClusterNStation_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterNStation[itr_clust],weight);
		h_dtRechitClusterMaxChamberRatio_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterMaxChamberRatio[itr_clust],weight);
		h_dtRechitClusterNChamber_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterNChamber[itr_clust],weight);
		h_dtRechitClusterMaxChamber_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterMaxChamber[itr_clust],weight);
		h_dtRechitClusterX_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterX[itr_clust],weight);
		h_dtRechitClusterY_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterY[itr_clust],weight);
		h_dtRechitClusterZ_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterZ[itr_clust],weight);
		h_dtRechitClusterEta_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterEta[itr_clust],weight);
		h_dtRechitClusterPhi_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterPhi[itr_clust],weight);
		h_dtRechitClusterTime_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterTime[itr_clust],weight);
		h_dtRechitClusterXSpread_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterXSpread[itr_clust],weight);
		h_dtRechitClusterYSpread_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterYSpread[itr_clust],weight);
		h_dtRechitClusterZSpread_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterZSpread[itr_clust],weight);
		h_dtRechitClusterEtaSpread_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterEtaSpread[itr_clust],weight);
		h_dtRechitClusterPhiSpread_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterPhiSpread[itr_clust],weight);
		h_dtRechitClusterTimeSpread_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterTimeSpread[itr_clust],weight);
		h_dtRechitClusterMajorAxis_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterMajorAxis[itr_clust],weight);
		h_dtRechitClusterMinorAxis_dPhiJetMET[itr_mX][itr_ctau]->Fill(dtRechitClusterMinorAxis[itr_clust],weight);
		if(passStations25 && passWheels25 && dtRechitClusterMaxStation[itr_clust]>2){
		  nRB1MatchCluster=0;
		  for(Int_t itr_rpc=0; itr_rpc<*nRPCRechits; itr_rpc++){
		    if(sqrt(pow(RPCRechitX[itr_rpc],2)+pow(RPCRechitY[itr_rpc],2))<475.){
		      dPhi_tmp = dtRechitClusterPhi[itr_clust] - RPCRechitPhi[itr_rpc];
		      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		      if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - RPCRechitEta[itr_rpc],2))<0.4){
			nRB1MatchCluster+=1;
		      }
		    }
		  }
		  nMB1MatchClusterAdjacentPlus = 0;
		  nMB1MatchClusterAdjacentMinus = 0;
		  for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
		    if(dtRechitStation[itr_dt]==1){
		      dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
		      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		      if(fabs(dPhi_tmp)<0.4){
			if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacentPlus+=1; }
			if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacentMinus+=1; }
		      }
		    }
		  }
		  h_nRB1Match_dPhiJetMET[itr_mX][itr_ctau]->Fill(nRB1MatchCluster);
		  if(dtRechitClusterMaxChamber[itr_clust]==-2){ h_nMB1MatchAdjacent_dPhiJetMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentPlus); }
		  else if(dtRechitClusterMaxChamber[itr_clust]==2){ h_nMB1MatchAdjacent_dPhiJetMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentMinus); }
		  else{
		    h_nMB1MatchAdjacent_dPhiJetMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentPlus);
		    h_nMB1MatchAdjacent_dPhiJetMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentMinus); 
		  }
		  if(passMB1){
		    h_nRB1Match_MB1Veto_dPhiJetMET[itr_mX][itr_ctau]->Fill(nRB1MatchCluster);
		    if(dtRechitClusterMaxChamber[itr_clust]==-2){ h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentPlus); }
		    else if(dtRechitClusterMaxChamber[itr_clust]==2){ h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentMinus); }
		    else{
		      h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentPlus);
		      h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentMinus); 
		    }
		  }
		}
	      }
	      
	      
	      if(!rpcBx.empty() && rpcSpread>0){ passRPCCR=true; }
	      /*if(itr_clust==0){
		if(hitRPCStation1>0){
		  nRPCStations1+=1;
		  if(hitRPCStation1>5){
		    nRPCStations5+=1;
		    if(hitRPCStation1>10){
		      nRPCStations10+=1;
		    }
		  }
		}
		if(hitRPCStation2>0){
		  nRPCStations1+=1;
		  if(hitRPCStation2>5){
		    nRPCStations5+=1;
		    if(hitRPCStation2>10){
		      nRPCStations10+=1;
		    }
		  }
		}
		if(hitRPCStation3>0){
		  nRPCStations1+=1;
		  if(hitRPCStation3>5){
		    nRPCStations5+=1;
		    if(hitRPCStation3>10){
		      nRPCStations10+=1;
		    }
		  }
		}
		if(hitRPCStation4>0){
		  nRPCStations1+=1;
		  if(hitRPCStation4>5){
		    nRPCStations5+=1;
		    if(hitRPCStation4>10){
		      nRPCStations10+=1;
		    }
		  }
		}
		if(hitRPCWheel1>0){
		  nRPCWheels1+=1;
		  if(hitRPCWheel1>5){
		    nRPCWheels5+=1;
		    if(hitRPCWheel1>10){
		      nRPCWheels10+=1;
		    }
		  }
		}
		if(hitRPCWheel2>0){
		  nRPCWheels1+=1;
		  if(hitRPCWheel2>5){
		    nRPCWheels5+=1;
		    if(hitRPCWheel2>10){
		      nRPCWheels10+=1;
		    }
		  }
		}
		if(hitRPCWheel0>0){
		  nRPCWheels1+=1;
		  if(hitRPCWheel0>5){
		    nRPCWheels5+=1;
		    if(hitRPCWheel0>10){
		      nRPCWheels10+=1;
		    }
		  }
		}
		if(hitRPCWheelm1>0){
		  nRPCWheels1+=1;
		  if(hitRPCWheelm1>5){
		    nRPCWheels5+=1;
		    if(hitRPCWheelm1>10){
		      nRPCWheels10+=1;
		    }
		  }
		}
		if(hitRPCWheelm2>0){
		  nRPCWheels1+=1;
		  if(hitRPCWheelm2>5){
		    nRPCWheels5+=1;
		    if(hitRPCWheelm2>10){
		      nRPCWheels10+=1;
		    }
		  }
		}
	      }
	      if(*nRPCRechits>750){
		nRPCStations1=4;
		nRPCStations5=4;
		nRPCStations10=4;
		nRPCWheels1=5;
		nRPCWheels5=5;
		nRPCWheels10=5;
		}*/
	      
	      if(passJet && passMB1){
		
		if(fabs(dPhiClusterMET)<1.0 && passMuon){
		  passFullVeto_clusterCR=true;
		  h_nRPCMatched_fullVeto_clusterMETCR[itr_mX][itr_ctau]->Fill(rpcBx.size(),weight);
		  h_rpcSpread_fullVeto_clusterMETCR[itr_mX][itr_ctau]->Fill(rpcSpread,weight);
		  h_rpcBx_fullVeto_clusterMETCR[itr_mX][itr_ctau]->Fill(rpcMedian,weight);
		  h_dPhiJetMET_fullVeto_clusterMETCR[itr_mX][itr_ctau]->Fill(fabs(dPhi_min),weight);
		  h_dtRechitClusterMaxStation_fullVeto_clusterMETCR[itr_mX][itr_ctau]->Fill(dtRechitClusterMaxStation[itr_clust],weight);
		  
		  if(dtRechitClusterMaxStation[itr_clust]>2 && !rpcBx.empty() && nStations25<3 && nWheels25<3){
		    nMB1MatchClusterAdjacentPlus = 0;
		    nMB1MatchClusterAdjacentMinus = 0;
		    nMB1MatchPi2 = 0;
		    nMB1MatchPi2AdjacentPlus = 0;
		    nMB1MatchPi2AdjacentMinus = 0;
		    nMB1MatchClusterAdjacent0p8Plus = 0;
		    nMB1MatchClusterAdjacent0p8Minus = 0;
		    nMB1MatchPi2Adjacent0p8Plus = 0;
		    nMB1MatchPi2Adjacent0p8Minus = 0;
		    for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
		      if(dtRechitStation[itr_dt]==1){
			dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
			if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			if(fabs(dPhi_tmp)<0.4){
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacentPlus+=1; }
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacentMinus+=1; }
			}
			if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacent0p8Plus+=1; }
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacent0p8Minus+=1; }
			}
			dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt] + pmRand*TMath::Pi()/2.0;
			if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			if(fabs(dPhi_tmp)<0.4){
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchPi2AdjacentPlus+=1; }
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchPi2AdjacentMinus+=1; }
			}
			if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchPi2Adjacent0p8Plus+=1; }
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchPi2Adjacent0p8Minus+=1; }
			}
			if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitEta[itr_dt]-dtRechitClusterEta[itr_clust],2))<0.4){
			  nMB1MatchPi2+=1; 
			}
		      }
		    }
		    if(nMB1MatchPi2<2){ h_nMB1MatchPi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(0.0); }
		    else{ h_nMB1MatchPi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(1.0); }
		    if(dtRechitClusterMaxChamber[itr_clust]==-2){ 
		      h_nMB1MatchAdjacent_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentPlus); 
		      h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchPi2AdjacentPlus); 
		      h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacent0p8Plus); 
		      h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchPi2Adjacent0p8Plus); 
		    }
		    else if(dtRechitClusterMaxChamber[itr_clust]==2){ 
		      h_nMB1MatchAdjacent_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentMinus); 
		      h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchPi2AdjacentMinus); 
		      h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacent0p8Minus); 
		      h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchPi2Adjacent0p8Minus); 
		    }
		    else{
		      h_nMB1MatchAdjacent_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentPlus);
		      h_nMB1MatchAdjacent_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacentMinus); 
		      h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchPi2AdjacentPlus);
		      h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchPi2AdjacentMinus); 
		      h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacent0p8Plus);
		      h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchClusterAdjacent0p8Minus); 
		      h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchPi2Adjacent0p8Plus);
		      h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_mX][itr_ctau]->Fill(nMB1MatchPi2Adjacent0p8Minus); 
		    }
		  }
		  
		  if(dtRechitClusterMaxStation[itr_clust]>2){
		    passMaxStation_clusterCR=true;
		    if(!rpcBx.empty()){
		      passRPCMatch_clusterCR=true;
		      if(rpcSpread==0){
			passRPCSpread_clusterCR=true;
			if(rpcMedian>=0.){
			  passRPCBx_clusterCR=true;
			  if(0==0){
			    passLepton_clusterCR=true;
			    if(nStations50<3 && nWheels50<3){
			      pass50Hits_clusterCR=true;
			      if(nStations25<3 && nWheels25<3){
				pass25Hits_clusterCR=true;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		  
		  if(!rpcBx.empty() && rpcSpread==0 && rpcMedian>=0.){ h_dtRechitClusterMaxStation_Nminus1_clusterMETCR[itr_mX][itr_ctau]->Fill(dtRechitClusterMaxStation[itr_clust],weight); }
		  if(!rpcBx.empty() && rpcSpread==0 && rpcMedian>=0.){ h_dPhiJetMET_Nminus1_clusterMETCR[itr_mX][itr_ctau]->Fill(fabs(dPhi_min),weight); }
		  if(!rpcBx.empty() && rpcSpread==0 && dtRechitClusterMaxStation[itr_clust]>2){ h_rpcBx_Nminus1_clusterMETCR[itr_mX][itr_ctau]->Fill(rpcMedian,weight); }
		  if(!rpcBx.empty() && rpcMedian>=0. && dtRechitClusterMaxStation[itr_clust]>2){ h_rpcSpread_Nminus1_clusterMETCR[itr_mX][itr_ctau]->Fill(rpcSpread,weight); }
		  if(dtRechitClusterMaxStation[itr_clust]>2){ h_nRPCMatched_Nminus1_clusterMETCR[itr_mX][itr_ctau]->Fill(rpcBx.size(),weight); }
		}

		if(!rpcBx.empty() && rpcSpread==0 && passMuon){
		  passFullVeto_rpcCR=true;
		  h_dPhiClusterMET_fullVeto_rpcCR[itr_mX][itr_ctau]->Fill(fabs(dPhiClusterMET),weight);
		  h_dPhiJetMET_fullVeto_rpcCR[itr_mX][itr_ctau]->Fill(fabs(dPhi_min),weight);
		  h_dtRechitClusterMaxStation_fullVeto_rpcCR[itr_mX][itr_ctau]->Fill(dtRechitClusterMaxStation[itr_clust],weight);
		  
		  if(dtRechitClusterMaxStation[itr_clust]>2){
		    passMaxStation_rpcCR=true;
		    h_dPhiClusterMET_Nminus1_rpcCR[itr_mX][itr_ctau]->Fill(fabs(dPhiClusterMET),weight);
		    h_dPhiJetMET_Nminus1_rpcCR[itr_mX][itr_ctau]->Fill(fabs(dPhi_min),weight);
		    if(fabs(dPhiClusterMET)<1.0){ 
		      passClusterMET_rpcCR=true;
		      if(fabs(dPhi_min)>0.6){
			passJetMET_rpcCR=true;
			if(0==0){
			  passLepton_rpcCR=true;
			  if(nStations50<3 && nWheels50<3){
			    pass50Hits_rpcCR=true;
			    if(nStations25<3 && nWheels25<3){
			      pass25Hits_rpcCR=true;
			    }
			  }
			}
		      }
		    }
		  }
		  if(fabs(dPhiClusterMET)<1.0){
		    h_dtRechitClusterMaxStation_Nminus1_rpcCR[itr_mX][itr_ctau]->Fill(dtRechitClusterMaxStation[itr_clust],weight); 
		  }
		}
		
		if(passMuon){
		  if(!rpcBx.empty() && rpcSpread==0){
		    if(rpcMedian>=0.){
		      if(dtRechitClusterMaxStation[itr_clust]>2){
			if(fabs(dPhi_min)>0.6){
			  if(fabs(dPhiClusterMET)<1.0){
			    if(0==0){
			      if(nStations50<3 && nWheels50<3){
				if(nStations25<3 && nWheels25<3){
				  h_dtRechitClusterSize_fullSelection_rpcCR[itr_mX][itr_ctau]->Fill(dtRechitClusterSize[itr_clust],weight);
				  h_dtRechitClusterSize_fullSelection_clusterMETCR[itr_mX][itr_ctau]->Fill(dtRechitClusterSize[itr_clust],weight);
				  if(dtRechitClusterSize[itr_clust]>100){
				    passSignalRegion=true;
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}

		if(passMuon){
		  if(dtRechitClusterMaxStation[itr_clust]==2 && fabs(dPhi_min)>0.6 && passMET && passOneJet && passJetMET && nStations25<3 && nWheels25<3 && fabs(dPhiClusterMET)<1.0 && dtRechitClusterSize[itr_clust]>=100){
		    nMB1MatchClusterAdjacentPlus = 0;
		    nMB1MatchClusterAdjacentMinus = 0;
		    nMB1MatchClusterAdjacent0p8Plus = 0;
		    nMB1MatchClusterAdjacent0p8Minus = 0;
		    for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
		      if(dtRechitStation[itr_dt]==1){
			dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
			if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			if(fabs(dPhi_tmp)<0.4){
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacentPlus+=1; }
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacentMinus+=1; }
			}
			if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacent0p8Plus+=1; }
			  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacent0p8Minus+=1; }
			}
		      }
		    }
		    if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
		      if(passNHFJet){
			passMB2CRwithNHFnoRPC = true;
		      }
		    }
		  }
		  if(!rpcBx.empty()){
		    if(dtRechitClusterMaxStation[itr_clust]==2){
		      if(fabs(dPhi_min)>0.6){
			if(passMET && passOneJet && passJetMET && nStations25<3 && nWheels25<3){
			  if(dtRechitClusterSize[itr_clust]>=100){
			    if(fabs(dPhiClusterMET)<1.0){
			      passMB2CR = true;
			      nMB1MatchClusterAdjacentPlus = 0;
			      nMB1MatchClusterAdjacentMinus = 0;
			      nMB1MatchClusterAdjacent0p8Plus = 0;
			      nMB1MatchClusterAdjacent0p8Minus = 0;
			      for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
				if(dtRechitStation[itr_dt]==1){
				  dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
				  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
				  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
				  if(fabs(dPhi_tmp)<0.4){
				    if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacentPlus+=1; }
				    if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacentMinus+=1; }
				  }
				  if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
				    if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacent0p8Plus+=1; }
				    if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacent0p8Minus+=1; }
				  }
				}
			      }
			      if(nMB1MatchClusterAdjacentPlus<5 && nMB1MatchClusterAdjacentMinus<5){
				passMB2CRwithAdjacent = true;
			      }
			      if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
				passMB2CRwithAdjacent0p8 = true;
				if (passNHFJet){
				  
				  //passMB2CRwithNHF = true;
				  
				  for(Int_t itr_llp=0; itr_llp<2; itr_llp++){
				    float decayR = sqrt(pow(gLLP_decay_vertex_x[itr_llp],2)+pow(gLLP_decay_vertex_y[itr_llp],2));
				    dtRechitClusterVec.SetPtEtaPhi(1,dtRechitClusterEta[itr_clust],dtRechitClusterPhi[itr_clust]);
				    if(decayR > 300 && decayR < 800 && fabs(gLLP_decay_vertex_z[itr_llp]) < 650. && !gLLP_plotted_id_sr2[itr_llp] && dtRechitClusterSize[itr_clust] > 100){
				      TVector3 llpDecay = TVector3(gLLP_decay_vertex_x[itr_llp],gLLP_decay_vertex_y[itr_llp],gLLP_decay_vertex_z[itr_llp]);
				      if (llpDecay.DeltaR(dtRechitClusterVec) < 0.5)
					{
					  h_decayVertexRadius_clusterReco_signalRegionEq2[itr_mX][itr_ctau]->Fill(sqrt(pow(gLLP_decay_vertex_x[itr_llp],2)+pow(gLLP_decay_vertex_y[itr_llp],2)));
					  h_decayVertexZ_clusterReco_signalRegionEq2[itr_mX][itr_ctau]->Fill(fabs(gLLP_decay_vertex_z[itr_llp]));
					  gLLP_plotted_id_sr2[itr_llp] = true;
					}
				    }
				  }
				  h_leadJetPt_MB2CR[itr_mX][itr_ctau]->Fill(jetPt[0],weight);
				  h_leadJetPtMET_MB2CR[itr_mX][itr_ctau]->Fill(jetPt[0]/(*MET),weight);
				  h_leadJetPt_MB2withMB1CR[itr_mX][itr_ctau]->Fill(jetPt[0],weight);
				  h_leadJetPtMET_MB2withMB1CR[itr_mX][itr_ctau]->Fill(jetPt[0]/(*MET),weight);
				  h_leadJetPt_MB1CR[itr_mX][itr_ctau]->Fill(jetPt[0],weight);
				  h_leadJetPtMET_MB1CR[itr_mX][itr_ctau]->Fill(jetPt[0]/(*MET),weight);
				}
			      }
			      if(nMB3MatchCluster<5 && nMB4MatchCluster<5){
				passMB2CRwithOther = true;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
		
	      }
	    
	      if(passMET && passOneJet && passJetMET && passStations25 && passWheels25 && passNHFJet && passMuon && passMB1 && fabs(dPhiClusterMET)<1.0 && !rpcBx.empty()){
		h_jetVetoPt_Nminus1_clusterMETCR[itr_mX][itr_ctau]->Fill(dtRechitClusterJetVetoPt[itr_clust],weight);
	      }
	      if(passMET && passOneJet && passJetMET && passStations25 && passWheels25 && passNHFJet && passJet && passMB1 && fabs(dPhiClusterMET)<1.0 && !rpcBx.empty()){
		h_muonVetoPt_Nminus1_clusterMETCR[itr_mX][itr_ctau]->Fill(dtRechitClusterMuonVetoPt[itr_clust],weight);
		h_muonLooseIDVetoPt_Nminus1_clusterMETCR[itr_mX][itr_ctau]->Fill(matchedLooseIDMuonPt);
	      }
	      
	     

	      if(passJetTightId && passMET && passOneJet && passJetMET && passStations25 && passWheels25 && passNHFJet){ passJetTightIdVeto = true; }
	      if(passJet && passMET && passOneJet && passJetMET && passStations25 && passWheels25 && passNHFJet){
		passJetVeto = true;
		if(passMuon){
		  passMuonVeto = true;
		  if(passMB1){
		    passMB1Veto = true;
		    if(dtRechitClusterMaxStation[itr_clust]>1){
		      passMaxStation = true;
		      if(fabs(dPhiClusterMET)<1.0 && dtRechitClusterSize[itr_clust]>=100){
			nMB1MatchClusterAdjacentPlus = 0;
			nMB1MatchClusterAdjacentMinus = 0;
			nMB1MatchClusterAdjacent0p8Plus = 0;
			nMB1MatchClusterAdjacent0p8Minus = 0;
			for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
			  if(dtRechitStation[itr_dt]==1){
			    dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
			    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			    if(fabs(dPhi_tmp)<0.4){
			      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacentPlus+=1; }
			      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacentMinus+=1; }
			    }
			    if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
			      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacent0p8Plus+=1; }
			      if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacent0p8Minus+=1; }
			    }
			  }
			}
			if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
			  if(passNHFJet){
			    passSRwithNHFnoRPC = true;
			  }
			}
		      }
		      if(!rpcBx.empty()){
			passRPCMatch = true;
			if(0==0){
			  passRPCSpread = true;
			  if(1>=0){
			    passRPCBx = true;
			    nMB1MatchClusterAdjacentPlus = 0;
			    nMB1MatchClusterAdjacentMinus = 0;
			    nMB1MatchClusterAdjacent0p8Plus = 0;
			    nMB1MatchClusterAdjacent0p8Minus = 0;
			    for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
			      if(dtRechitStation[itr_dt]==1){
				dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
				if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
				if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
				if(fabs(dPhi_tmp)<0.4){
				  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacentPlus+=1; }
				  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacentMinus+=1; }
				}
				if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
				  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1MatchClusterAdjacent0p8Plus+=1; }
				  if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1MatchClusterAdjacent0p8Minus+=1; }
				}
			      }
			    }
			    if(nMB1MatchClusterAdjacentPlus<5 && nMB1MatchClusterAdjacentMinus<5){
			      passAdjacentMB1 = true;
			    }
			    if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
			      passAdjacent0p8MB1 = true;
			    }
			    if(fabs(dPhiClusterMET)>=1.0){
			      passClusterMETA = true;
			      passClusterMETB = true;
			    }
			    if(dtRechitClusterSize[itr_clust]>=100){
			      passClusterSizeB = true;
			    }
			    else{
			      passClusterSizeA = true;
			      passClusterSizeC = true;
			    }	
			    if(dtRechitClusterMaxStation[itr_clust]==2){ passMB2Cluster = true; }
			    if(fabs(dPhiClusterMET)<1.0){
			      passClusterMET = true;
			      passClusterMETC = true;
			      if(dtRechitClusterSize[itr_clust]>maxClusterSizeSR){ maxClusterSizeSR = dtRechitClusterSize[itr_clust]; }
			      if(dtRechitClusterMaxStation[itr_clust]==2){ h_dtRechitClusterSize_MB2_signalRegion[itr_mX][itr_ctau]->Fill(dtRechitClusterSize[itr_clust],weight); }
			      if(dtRechitClusterMaxStation[itr_clust]==3){ h_dtRechitClusterSize_MB3_signalRegion[itr_mX][itr_ctau]->Fill(dtRechitClusterSize[itr_clust],weight); }
			      if(dtRechitClusterMaxStation[itr_clust]==4){ h_dtRechitClusterSize_MB4_signalRegion[itr_mX][itr_ctau]->Fill(dtRechitClusterSize[itr_clust],weight); }

			      if(dtRechitClusterSize[itr_clust]>=100){ 
				passClusterSize = true;
				
				if(dtRechitClusterMaxStation[itr_clust]==2){
				  passMB2CRwithNHF = true;
				}
				if(dtRechitClusterMaxStation[itr_clust]>2){
				  
				  h_jetChargedHadronicEnergyFraction_SR[itr_mX][itr_ctau]->Fill(chargedHadFraction_mindPhi);
				  h_jetChargedEMEnergyFraction_SR[itr_mX][itr_ctau]->Fill(chargedEMFraction_mindPhi);
				  h_jetNeutralHadronicEnergyFraction_SR[itr_mX][itr_ctau]->Fill(neutralHadFraction_mindPhi);
				  h_jetNeutralEMEnergyFraction_SR[itr_mX][itr_ctau]->Fill(neutralEMFraction_mindPhi);
				  h_leadingJetChargedHadronicEnergyFraction_SR[itr_mX][itr_ctau]->Fill(jetChargedHadronEnergyFraction[0]);
				  h_leadingJetChargedEMEnergyFraction_SR[itr_mX][itr_ctau]->Fill(jetElectronEnergyFraction[0]);
				  h_leadingJetNeutralHadronicEnergyFraction_SR[itr_mX][itr_ctau]->Fill(jetNeutralHadronEnergyFraction[0]);
				  h_leadingJetNeutralEMEnergyFraction_SR[itr_mX][itr_ctau]->Fill(jetPhotonEnergyFraction[0]);
				  
				  if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
				    if (passNHFJet){
				      for(Int_t itr_llp=0; itr_llp<2; itr_llp++){
					float decayR = sqrt(pow(gLLP_decay_vertex_x[itr_llp],2)+pow(gLLP_decay_vertex_y[itr_llp],2));
					dtRechitClusterVec.SetPtEtaPhi(1,dtRechitClusterEta[itr_clust],dtRechitClusterPhi[itr_clust]);
					if(decayR > 300 && decayR < 800 && fabs(gLLP_decay_vertex_z[itr_llp]) < 650. && !gLLP_plotted_id_sr1[itr_llp] && dtRechitClusterSize[itr_clust] > 100){
					  TVector3 llpDecay = TVector3(gLLP_decay_vertex_x[itr_llp],gLLP_decay_vertex_y[itr_llp],gLLP_decay_vertex_z[itr_llp]);
					  if (llpDecay.DeltaR(dtRechitClusterVec) < 0.5)
					    {
					      h_decayVertexRadius_clusterReco_signalRegionGt2[itr_mX][itr_ctau]->Fill(sqrt(pow(gLLP_decay_vertex_x[itr_llp],2)+pow(gLLP_decay_vertex_y[itr_llp],2)));
					      h_decayVertexZ_clusterReco_signalRegionGt2[itr_mX][itr_ctau]->Fill(fabs(gLLP_decay_vertex_z[itr_llp]));
					      gLLP_plotted_id_sr1[itr_llp] = true;
					    }
					}
				      }
				      h_leadJetPt_SR[itr_mX][itr_ctau]->Fill(jetPt[0],weight);
				      h_leadJetPtMET_SR[itr_mX][itr_ctau]->Fill(jetPt[0]/(*MET),weight);
				      h_leadJetPt_MB1HitsCR[itr_mX][itr_ctau]->Fill(jetPt[0],weight);
				      h_leadJetPtMET_MB1HitsCR[itr_mX][itr_ctau]->Fill(jetPt[0]/(*MET),weight);
				      h_leadJetPt_MB1CR[itr_mX][itr_ctau]->Fill(jetPt[0],weight);
				      h_leadJetPtMET_MB1CR[itr_mX][itr_ctau]->Fill(jetPt[0]/(*MET),weight);
				      h_leadJetEta_SR[itr_mX][itr_ctau]->Fill(jetEta[0],weight);
				    }
				  }
				  if((dtRechitClusterMaxStation[itr_clust]==3 && nMB2MatchCluster<5 && nMB4MatchCluster<5) || (dtRechitClusterMaxStation[itr_clust]==4 && nMB2MatchCluster<5 && nMB3MatchCluster<5)){
				    passOtherStations = true;
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	      
	      if(dtRechitClusterNSegmentStation1[itr_clust]==0 && dtRechitClusterMaxStation[itr_clust]>2 && passMET && passOneJet && passJetMET && passStations25 && passWheels25){
		passMB1CR = true;
		if(passJet){
		  passJetVetoMB1CR = true;
		  if(passMuon){
		    passMuonVetoMB1CR = true;
		    if(passMuonLoose){ passMuonVetoLooseMB1CR = true; }
		    if(!rpcBx.empty()){
		      passRpcMatchMB1CR = true;
		      if(fabs(dPhiClusterMET)<1.0){
			passClusterMETMB1CR = true;
			clusterEta.push_back(dtRechitClusterEta[itr_clust]);
			clusterPhi.push_back(dtRechitClusterPhi[itr_clust]);
			clusterSize.push_back(dtRechitClusterSize[itr_clust]);
			if(dtRechitClusterSize[itr_clust]>maxClusterSizeMB1CR){
			    maxClusterSizeMB1CR=dtRechitClusterSize[itr_clust];
			}
			if(dtRechitClusterSize[itr_clust]>100){
			  passClusterSizeMB1CR = true;
			}
		      }
		    }
		  }
		}
	      }
	      
	    }
	  }
	}
	if(passFullVeto_clusterCR){ nPassFullVeto_clusterCR+=weight; }
	if(passRPCMatch_clusterCR){ nPassRPCMatch_clusterCR+=weight; }
	if(passRPCSpread_clusterCR){ nPassRPCSpread_clusterCR+=weight; }
	if(passRPCBx_clusterCR){ nPassRPCBx_clusterCR+=weight; }
	if(passMaxStation_clusterCR){ 
	  nPassMaxStation_clusterCR+=weight; 
	  if(maxClusterSize>150){
	    h_nDtSegsStation_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation1);
	    h_nDtSegsStation_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation2);
	    h_nDtSegsStation_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation3);
	    h_nDtSegsStation_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation4);
	    
	    h_nDtSegsWheel_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel2);
	    h_nDtSegsWheel_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel1);
	    h_nDtSegsWheel_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel0);
	    h_nDtSegsWheel_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheelm1);
	    h_nDtSegsWheel_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheelm2);
	    
	    for(int i=0; i<4; i++){
	      for(int j=0; j<5; j++){
		h_nDtSegsChamber_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segChambers[i][j]);
	      }
	    }
	    
	    //h_nDtSegs_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(*nDtSeg);
	  

	    h_nStations1_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations1,weight);
	    h_nStations25_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations25,weight);
	    h_nStations50_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations50,weight);
	    h_nWheels1_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels1,weight);
	    h_nWheels25_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels25,weight);
	    h_nWheels50_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels50,weight);
	    
	    h_nStations1Seg_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations1Seg);
	    h_nStations5Seg_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations5Seg);
	    h_nStations10Seg_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations10Seg);
	    h_nWheels1Seg_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels1Seg);
	    h_nWheels5Seg_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels5Seg);
	    h_nWheels10Seg_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels10Seg);

	    
	    h_nRPCStations1_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations1,weight);
	    h_nRPCStations5_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations5,weight);
	    h_nRPCStations10_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations10,weight);
	    h_nRPCWheels1_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels1,weight);
	    h_nRPCWheels5_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels5,weight);
	    h_nRPCWheels10_150hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels10,weight);
	  }
	  else if(maxClusterSize>100){
	    h_nDtSegsStation_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation1);
	    h_nDtSegsStation_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation2);
	    h_nDtSegsStation_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation3);
	    h_nDtSegsStation_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation4);
	    
	    h_nDtSegsWheel_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel2);
	    h_nDtSegsWheel_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel1);
	    h_nDtSegsWheel_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel0);
	    h_nDtSegsWheel_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheelm1);
	    h_nDtSegsWheel_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheelm2);

	    for(int i=0; i<4; i++){
	      for(int j=0; j<5; j++){
		h_nDtSegsChamber_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segChambers[i][j]);
	      }
	    }
	    
	    //h_nDtSegs_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(*nDtSeg);
	  
	    h_nStations1_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations1,weight);
	    h_nStations25_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations25,weight);
	    h_nStations50_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations50,weight);
	    h_nWheels1_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels1,weight);
	    h_nWheels25_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels25,weight);
	    h_nWheels50_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels50,weight);
	    
	    h_nStations1Seg_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations1Seg);
	    h_nStations5Seg_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations5Seg);
	    h_nStations10Seg_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations10Seg);
	    h_nWheels1Seg_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels1Seg);
	    h_nWheels5Seg_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels5Seg);
	    h_nWheels10Seg_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels10Seg);

	    h_nRPCStations1_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations1,weight);
	    h_nRPCStations5_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations5,weight);
	    h_nRPCStations10_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations10,weight);
	    h_nRPCWheels1_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels1,weight);
	    h_nRPCWheels5_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels5,weight);
	    h_nRPCWheels10_100hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels10,weight);
	  }
	  else{
	    h_nDtSegsStation_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation1);
	    h_nDtSegsStation_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation2);
	    h_nDtSegsStation_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation3);
	    h_nDtSegsStation_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segStation4);
	    
	    h_nDtSegsWheel_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel2);
	    h_nDtSegsWheel_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel1);
	    h_nDtSegsWheel_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheel0);
	    h_nDtSegsWheel_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheelm1);
	    h_nDtSegsWheel_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segWheelm2);

	    for(int i=0; i<4; i++){
	      for(int j=0; j<5; j++){
		h_nDtSegsChamber_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(segChambers[i][j]);
	      }
	    }

	    //h_nDtSegs_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(*nDtSeg);
	  
	    h_nStations1_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations1,weight);
	    h_nStations25_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations25,weight);
	    h_nStations50_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations50,weight);
	    h_nWheels1_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels1,weight);
	    h_nWheels25_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels25,weight);
	    h_nWheels50_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels50,weight);
	    
	    h_nStations1Seg_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations1Seg);
	    h_nStations5Seg_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations5Seg);
	    h_nStations10Seg_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nStations10Seg);
	    h_nWheels1Seg_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels1Seg);
	    h_nWheels5Seg_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels5Seg);
	    h_nWheels10Seg_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nWheels10Seg);

	    h_nRPCStations1_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations1,weight);
	    h_nRPCStations5_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations5,weight);
	    h_nRPCStations10_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCStations10,weight);
	    h_nRPCWheels1_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels1,weight);
	    h_nRPCWheels5_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels5,weight);
	    h_nRPCWheels10_50hits_clusterMETCR[itr_mX][itr_ctau]->Fill(nRPCWheels10,weight);
	  }
	}
	if(passLepton_clusterCR){ nPassLepton_clusterCR+=weight; }
	if(pass50Hits_clusterCR){ nPass50Hits_clusterCR+=weight; }
	if(pass25Hits_clusterCR){ nPass25Hits_clusterCR+=weight; }
	
	if(passRPCCR){ nPassRPCCR+=weight; }
	if(passFullVeto_rpcCR){ nPassFullVeto_rpcCR+=weight; }
	if(passClusterMET_rpcCR){ nPassClusterMET_rpcCR+=weight; }
	if(passMaxStation_rpcCR){ nPassMaxStation_rpcCR+=weight; }
	if(passJetMET_rpcCR){ nPassJetMET_rpcCR+=weight; }
	if(passLepton_rpcCR){ nPassLepton_rpcCR+=weight; }
	if(pass50Hits_rpcCR){ nPass50Hits_rpcCR+=weight; }
	if(pass25Hits_rpcCR){ nPass25Hits_rpcCR+=weight; }
	if(passSignalRegion){ nPassSignalRegion+=weight; }
	if(passMB2CR){ nPassMB2CR+=weight; }
	if(passMB2CRwithAdjacent){ nPassMB2CRwithAdjacent+=weight; }
	if(passMB2CRwithAdjacent0p8){ nPassMB2CRwithAdjacent0p8+=weight; }
	if(passMB2CRwithOther){ nPassMB2CRwithOther+=weight; }
	if(passMB2CRwithNHF){ nPassMB2CRwithNHF+=weight; }
	if(passMB2CRwithNHFnoRPC){ nPassMB2CRwithNHFnoRPC+=weight; }
	if(passSRwithNHFnoRPC){ nPassSRwithNHFnoRPC+=weight; }

	h_nDtRechitClustersVeto_dPhiJetMET[itr_mX][itr_ctau]->Fill(nClustersVeto_dPhiJetMET,weight);
	if(clusterPhi.size()>1){
	  for(Int_t i=0; i<clusterPhi.size()-1; i++){
	    for(Int_t j=i+1; j<clusterPhi.size(); j++){
	      dPhi_tmp = clusterPhi[i] - clusterPhi[j];
	      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	      h_dtRechitClustersVetoDR_dPhiJetMET[itr_mX][itr_ctau]->Fill(sqrt(pow(dPhi_tmp,2)+pow(clusterEta[i]-clusterEta[j],2)),weight);
	      if(sqrt(pow(dPhi_tmp,2)+pow(clusterEta[i]-clusterEta[j],2))<0.6){
		clusterSizeTotal+=clusterSize[i];
		clusterSizeTotal+=clusterSize[j];
		break;
	      }
	    }
	  }
	}
	
	if(passMET){
	  h_efficiency[itr_mX][itr_ctau]->Fill(0.0,weight);
	  h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(0.0,weight);
	  if(passOneJet){
	    h_efficiency[itr_mX][itr_ctau]->Fill(1,weight);
	    h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(1,weight);
	    if(passJetMET){
	      h_efficiency[itr_mX][itr_ctau]->Fill(2,weight);
	      h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(2,weight);
	      if(passStations25){
		h_efficiency[itr_mX][itr_ctau]->Fill(3,weight);
		h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(3,weight);
		if(passWheels25){
		  h_efficiency[itr_mX][itr_ctau]->Fill(4,weight);
		  h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(4,weight);
		  if(passNoVetoCluster){
		    h_efficiency[itr_mX][itr_ctau]->Fill(5,weight);
		    h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(5,weight);
		    if(passNHFJet){
		      h_efficiency[itr_mX][itr_ctau]->Fill(6,weight);
		    }
		  }
		}
	      }
	    }
	  }
	}
	if(passJetVeto){ h_efficiency[itr_mX][itr_ctau]->Fill(7,weight); }
	if(passJetTightIdVeto){ h_efficiency[itr_mX][itr_ctau]->Fill(55,weight); }
	if(passMuonVeto){ h_efficiency[itr_mX][itr_ctau]->Fill(8,weight); }
	if(passMB1Veto){ h_efficiency[itr_mX][itr_ctau]->Fill(9,weight); }
	if(passMaxStation){ h_efficiency[itr_mX][itr_ctau]->Fill(10,weight); }
	if(passRPCMatch){ h_efficiency[itr_mX][itr_ctau]->Fill(11,weight); }
	if(passRPCSpread){ h_efficiency[itr_mX][itr_ctau]->Fill(12,weight); }
	if(passRPCBx){ 
	  h_efficiency[itr_mX][itr_ctau]->Fill(13,weight); 
	  if(passAdjacent0p8MB1){
	    h_efficiency[itr_mX][itr_ctau]->Fill(14,weight); 
	    if(passClusterMET){ 
	      h_efficiency[itr_mX][itr_ctau]->Fill(15,weight); 
	      h_dtRechitClusterSize_signalRegion[itr_mX][itr_ctau]->Fill(maxClusterSizeSR,weight);
	      h_dtRechitClusterSize_signalRegionUnweighted[itr_mX][itr_ctau]->Fill(maxClusterSizeSR);
	    }
	    if(passClusterSize){ 
	      h_efficiency[itr_mX][itr_ctau]->Fill(16,weight); 
	      if(passAdjacent0p8MB1){
		h_efficiency[itr_mX][itr_ctau]->Fill(17,weight); 
		if(passNHFJet){ 
		  h_efficiency[itr_mX][itr_ctau]->Fill(18,weight); 
		  if(passMB2Cluster){ 
		    h_efficiency[itr_mX][itr_ctau]->Fill(19,weight); 
		    h_efficiency[itr_mX][itr_ctau]->Fill(20,weight*weight); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(21,weight*(*pileupWeightUp)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(22,weight*(*pileupWeightDown));
		    //h_efficiency[itr_mX][itr_ctau]->Fill(23,weight*(*sf_facScaleUp)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(24,weight*(*sf_facScaleDown)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(25,weight*(*sf_renScaleUp)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(26,weight*(*sf_renScaleDown));
		    //h_efficiency[itr_mX][itr_ctau]->Fill(27,weight*(*sf_facRenScaleUp)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(28,weight*(*sf_facRenScaleDown)); 
		  }
		  else{
		    h_efficiency[itr_mX][itr_ctau]->Fill(30,weight); 
		    h_efficiency[itr_mX][itr_ctau]->Fill(31,weight*weight); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(32,weight*(*pileupWeightUp)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(33,weight*(*pileupWeightDown));
		    //h_efficiency[itr_mX][itr_ctau]->Fill(34,weight*(*sf_facScaleUp)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(35,weight*(*sf_facScaleDown)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(36,weight*(*sf_renScaleUp)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(37,weight*(*sf_renScaleDown));
		    //h_efficiency[itr_mX][itr_ctau]->Fill(38,weight*(*sf_facRenScaleUp)); 
		    //h_efficiency[itr_mX][itr_ctau]->Fill(39,weight*(*sf_facRenScaleDown)); 
		  }
		}
	      }
	    }
	    if(passAdjacent0p8MB1 && passNHFJet){
	      if(passClusterMETC && passClusterSizeC){
		if(passMB2Cluster){ 
		  h_efficiency[itr_mX][itr_ctau]->Fill(40,weight); 
		  h_efficiency[itr_mX][itr_ctau]->Fill(41,weight*weight);
		}
		else{
		  h_efficiency[itr_mX][itr_ctau]->Fill(42,weight); 
		  h_efficiency[itr_mX][itr_ctau]->Fill(43,weight*weight);
		}
	      }
	      if(passClusterMETB && passClusterSizeB){
		if(passMB2Cluster){ 
		  h_efficiency[itr_mX][itr_ctau]->Fill(45,weight); 
		  h_efficiency[itr_mX][itr_ctau]->Fill(46,weight*weight);
		}
		else{
		  h_efficiency[itr_mX][itr_ctau]->Fill(47,weight); 
		  h_efficiency[itr_mX][itr_ctau]->Fill(48,weight*weight);
		}
	      }
	      if(passClusterMETA && passClusterSizeA){
		if(passMB2Cluster){ 
		  h_efficiency[itr_mX][itr_ctau]->Fill(50,weight); 
		  h_efficiency[itr_mX][itr_ctau]->Fill(51,weight*weight);
		}
		else{
		  h_efficiency[itr_mX][itr_ctau]->Fill(52,weight); 
		  h_efficiency[itr_mX][itr_ctau]->Fill(53,weight*weight);
		}
	      }
	    }
	  }
	}
	
	if(passMB1CR){ h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(6,weight); }
	if(passJetVetoMB1CR){ h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(7,weight); }
	if(passMuonVetoMB1CR){ h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(8,weight); }
	if(passRpcMatchMB1CR){ h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(9,weight); }
	if(passClusterMETMB1CR){ 
	  h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(10,weight); 
	  h_dtRechitClusterSize_signalRegionNew[itr_mX][itr_ctau]->Fill(maxClusterSizeMB1CR,weight);
	  h_dtRechitClusterSizeTotal_signalRegionNew[itr_mX][itr_ctau]->Fill(max(maxClusterSizeMB1CR,clusterSizeTotal),weight);
	}
	if(passClusterSizeMB1CR){ h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(11,weight); }
	if(passMuonVetoLooseMB1CR){ h_efficiency_MB1CR[itr_mX][itr_ctau]->Fill(12,weight); }

	evtNum+=1;	
      }
      //_file->Close();
    }
    
    cout << setprecision(3);
    myout << setprecision(3);
    myout << mX[itr_mX] << "," << ctau[itr_ctau];
    for(int i=1; i<57; i++){
      if(i > 1){
	if(h_efficiency[itr_mX][itr_ctau]->GetBinContent(i-1)>0){
	  cout << h_efficiency[itr_mX][itr_ctau]->GetBinContent(i) << " (" << h_efficiency[itr_mX][itr_ctau]->GetBinContent(i) / h_efficiency[itr_mX][itr_ctau]->GetBinContent(i-1) << ")" << endl;
	}
	else{
	  cout << h_efficiency[itr_mX][itr_ctau]->GetBinContent(i) << endl;
	}
	if(i<=19){
	  myout << "," << 100.0 * h_efficiency[itr_mX][itr_ctau]->GetBinContent(i) / h_efficiency[itr_mX][itr_ctau]->GetBinContent(i-1) << " (" << 100.0 * h_efficiency[itr_mX][itr_ctau]->GetBinContent(i) / h_efficiency[itr_mX][itr_ctau]->GetBinContent(1) << ") ";
	}
	if(i==20 || i==21 || i==31 || i==32 || (i>=41 && i<=44) || (i>=46 && i<=49) || (i>=51 && i<=54)){
	  myout << setprecision(6);
	  myout << "," << h_efficiency[itr_mX][itr_ctau]->GetBinContent(i);
	}
      }
      else{
	cout << h_efficiency[itr_mX][itr_ctau]->GetBinContent(i) << " (" << h_efficiency[itr_mX][itr_ctau]->GetBinContent(i) / (1000*0.01*48.58*(59.74+41.53+35.92)) << ")" << endl;
	myout << "," << 100.0 * h_efficiency[itr_mX][itr_ctau]->GetBinContent(i) / (1000*0.01*48.58*(59.74+41.53+35.92));
      }
    }
    myout << "\n";
    cout << " " << endl;
    cout << "nPassMB2CR: " << nPassMB2CR << endl;
    cout << "nPassMB2CRwithAdjacent: " << nPassMB2CRwithAdjacent << endl;
    cout << "nPassMB2CRwithAdjacent0p8: " << nPassMB2CRwithAdjacent0p8 << endl;
    cout << "nPassMB2CRwithOther: " << nPassMB2CRwithOther << endl;
    cout << "nPassMB2CRwithNHF: " << nPassMB2CRwithNHF << endl;
    cout << "nPassMB2CRwithNHFnoRPC: " << nPassMB2CRwithNHFnoRPC << endl;
    cout << "nPassSRwithNHFnoRPC: " << nPassSRwithNHFnoRPC << endl;
    cout << " " << endl;
    /*cout << "nPassNoVeto: " << nPassNoVeto << endl;
    cout << " " << endl;
    cout << "nPassClusterCR: " << nPassClusterCR << endl;
    cout << "nPassFullVeto_clusterCR: " << nPassFullVeto_clusterCR << endl;
    cout << "nPassMaxStation_clusterCR: " << nPassMaxStation_clusterCR << endl;
    cout << "nPassRPCMatch_clusterCR: " << nPassRPCMatch_clusterCR << endl;
    cout << "nPassRPCSpread_clusterCR: " << nPassRPCSpread_clusterCR << endl;
    cout << "nPassRPCBx_clusterCR: " << nPassRPCBx_clusterCR << endl;
    cout << "nPassLepton_clusterCR: " << nPassLepton_clusterCR << endl;
    cout << "nPass50Hits_clusterCR: " << nPass50Hits_clusterCR << endl;
    cout << "nPass25Hits_clusterCR: " << nPass25Hits_clusterCR << endl;
    cout << " " << endl;
    cout << "nPassRPCCR: " << nPassRPCCR << endl;
    cout << "nPassFullVeto_rpcCR: " << nPassFullVeto_rpcCR << endl;
    cout << "nPassMaxStation_rpcCR: " << nPassMaxStation_rpcCR << endl;
    cout << "nPassClusterMET_rpcCR: " << nPassClusterMET_rpcCR << endl;
    cout << "nPassJetMET_rpcCR: " << nPassJetMET_rpcCR << endl;
    cout << "nPassLepton_rpcCR: " << nPassLepton_rpcCR << endl;
    cout << "nPass50Hits_rpcCR: " << nPass50Hits_rpcCR << endl;
    cout << "nPass25Hits_rpcCR: " << nPass25Hits_rpcCR << endl;
    cout << "nPassSignalRegion: " << nPassSignalRegion << endl;
    */
    }
  }

  _ofile->Write();
  _ofile->Close();

}
    

  
  
