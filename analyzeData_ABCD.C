#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
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

void analyzeData_ABCD(){

  TString name;
  TString years[3] = {"2016","2017","2018"};
  TString runNames[3] = {"Run2016","Run2017","17Sept2018_Run2018"};
  TString dates[20] = {"07Aug17","17Nov2017","17Sep2018"};
  //TString dir("/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/driftTube/V1p17/Data");
  TString dir("/storage/af/user/mcitron/skims/v3/");
  TFile *_ofile = TFile::Open("outData_ABCD_jetVeto20.root","RECREATE");

  TH1D *h_dtRechitClusterSize_dPhiJetMETLow_rpcCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiJetMETHigh_rpcCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiJetMETLow_rpcCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_dPhiJetMETHigh_rpcCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETLow_rpcCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETHigh_rpcCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETLow_rpcSpreadCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETHigh_rpcSpreadCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETLow_rpcCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETHigh_rpcCRnoLepton[4];
  
  TH1D *h_dtRechitClusterSize_rpcMatch_jetMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcNoMatch_jetMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETLow_jetMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETHigh_jetMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcMatch_jetMETCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_rpcNoMatch_jetMETCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETLow_jetMETCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_dPhiClusterMETHigh_jetMETCRnoLepton[4];

  TH1D *h_dtRechitClusterSize_rpcGood_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcBad_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcMatch_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcNoMatch_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcSpread_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcNoSpread_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcNegativeSpread_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcNegativeNoSpread_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiJetMETLow_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_dPhiJetMETHigh_clusterMETCRmuonVeto[4];
  TH1D *h_dtRechitClusterSize_rpcMatch_clusterMETCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_rpcNoMatch_clusterMETCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_dPhiJetMETLow_clusterMETCRnoLepton[4];
  TH1D *h_dtRechitClusterSize_dPhiJetMETHigh_clusterMETCRnoLepton[4];

  TH1D *h_dtRechitClusterPhi_fullSelectionMB2_rpcCR[4];
  TH1D *h_dtRechitClusterPhi_fullSelectionMB34_rpcCR[4];
  TH1D *h_metPhi_fullSelectionMB2_rpcCR[4];
  TH1D *h_metPhi_fullSelectionMB34_rpcCR[4];

  TH1D *h_dtRechitClusterSize_fullSelection_rpcCR[4];
  TH1D *h_dtRechitClusterSize_lowClusterMET_fullSelection_rpcCR[4];
  TH1D *h_dtRechitClusterSize_highClusterMET_fullSelection_rpcCR[4];
  TH1D *h_dtRechitClusterSize_fullSelection_clusterMETCR[4];
  TH1D *h_dtRechitClusterSize_goodRPC_fullSelection_clusterMETCR[4];
  TH1D *h_dtRechitClusterSize_badRPC_fullSelection_clusterMETCR[4];
  TH1D *h_dtRechitClusterSize_lowClusterMET_lowClusterSize_SR[4];
  TH1D *h_dtRechitClusterSize_highClusterMET_lowClusterSize_SR[4];
  TH1D *h_dtRechitClusterSize_highClusterMET_highClusterSize_SR[4];

  TH1D *h_dPhiClusterRPC_fullVeto[4];
  TH1D *h_dZClusterRPC_fullVeto[4];
  TH1D *h_rpcSpread_invertedJetVeto[4];
  TH1D *h_rpcSpread_invertedJetVeto_muonVeto[4];
  TH1D *h_rpcSpread_fullVeto[4];
  TH1D *h_rpcSpread_fullVeto_negBx[4];

  TH1D *h_nRPCMatched_fullVeto_clusterMETCR[4];
  TH1D *h_rpcSpread_fullVeto_clusterMETCR[4];
  TH1D *h_rpcBx_fullVeto_clusterMETCR[4];
  TH1D *h_dPhiJetMET_fullVeto_clusterMETCR[4];
  TH1D *h_dtRechitClusterMaxStation_fullVeto_clusterMETCR[4];

  TH1D *h_nRPCMatched_Nminus1_clusterMETCR[4];
  TH1D *h_rpcSpread_Nminus1_clusterMETCR[4];
  TH1D *h_rpcBx_Nminus1_clusterMETCR[4];
  TH1D *h_dPhiJetMET_Nminus1_clusterMETCR[4];
  TH1D *h_dtRechitClusterMaxStation_Nminus1_clusterMETCR[4];
  TH1D *h_nMB1Matched_Nminus1_clusterMETCR[4];
  TH1D *h_nMB1MatchedAdjacent_Nminus1_clusterMETCR[4];
  TH1D *h_jetVetoPt_Nminus1_clusterMETCR[4];
  TH1D *h_muonVetoPt_Nminus1_clusterMETCR[4];
  TH1D *h_muonLooseIDVetoPt_Nminus1_clusterMETCR[4];
  
  TH1D *h_nRPCMatched_Nminus1_MB1CR[4];
  TH1D *h_nMB1MatchedAdjacent_Nminus1_MB1CR[4];
  TH1D *h_dPhiClusterMET_Nminus1_MB1CR[4];
  TH1D *h_jetVetoPt_Nminus1_MB1CR[4];
  TH1D *h_muonVetoPt_Nminus1_MB1CR[4];
  TH1D *h_muonLooseIDVetoPt_Nminus1_MB1CR[4];

  TH1D *h_dPhiClusterMET_fullVeto_rpcCR[4];
  TH1D *h_dPhiJetMET_fullVeto_rpcCR[4];
  TH1D *h_dtRechitClusterMaxStation_fullVeto_rpcCR[4];

  TH1D *h_dPhiClusterMET_Nminus1_rpcCR[4];
  TH1D *h_dPhiJetMET_Nminus1_rpcCR[4];
  TH1D *h_dtRechitClusterMaxStation_Nminus1_rpcCR[4];

  TH1D *h_nStations1_50hits_clusterMETCR[4];
  TH1D *h_nStations25_50hits_clusterMETCR[4];
  TH1D *h_nStations50_50hits_clusterMETCR[4];
  TH1D *h_nStations1_100hits_clusterMETCR[4];
  TH1D *h_nStations25_100hits_clusterMETCR[4];
  TH1D *h_nStations50_100hits_clusterMETCR[4];
  TH1D *h_nStations1_150hits_clusterMETCR[4];
  TH1D *h_nStations25_150hits_clusterMETCR[4];
  TH1D *h_nStations50_150hits_clusterMETCR[4];

  TH1D *h_nWheels1_50hits_clusterMETCR[4];
  TH1D *h_nWheels25_50hits_clusterMETCR[4];
  TH1D *h_nWheels50_50hits_clusterMETCR[4];
  TH1D *h_nWheels1_100hits_clusterMETCR[4];
  TH1D *h_nWheels25_100hits_clusterMETCR[4];
  TH1D *h_nWheels50_100hits_clusterMETCR[4];
  TH1D *h_nWheels1_150hits_clusterMETCR[4];
  TH1D *h_nWheels25_150hits_clusterMETCR[4];
  TH1D *h_nWheels50_150hits_clusterMETCR[4];

  TH1D *h_nDtSegsStation_50hits_clusterMETCR[4];
  TH1D *h_nDtSegsStation_100hits_clusterMETCR[4];
  TH1D *h_nDtSegsStation_150hits_clusterMETCR[4];
  TH1D *h_nDtSegsWheel_50hits_clusterMETCR[4];
  TH1D *h_nDtSegsWheel_100hits_clusterMETCR[4];
  TH1D *h_nDtSegsWheel_150hits_clusterMETCR[4];
  TH1D *h_nDtSegsChamber_50hits_clusterMETCR[4];
  TH1D *h_nDtSegsChamber_100hits_clusterMETCR[4];
  TH1D *h_nDtSegsChamber_150hits_clusterMETCR[4];
  TH1D *h_nDtSegs_50hits_clusterMETCR[4];
  TH1D *h_nDtSegs_100hits_clusterMETCR[4];
  TH1D *h_nDtSegs_150hits_clusterMETCR[4];

  TH1D *h_nStations1Seg_50hits_clusterMETCR[4];
  TH1D *h_nStations5Seg_50hits_clusterMETCR[4];
  TH1D *h_nStations10Seg_50hits_clusterMETCR[4];
  TH1D *h_nStations1Seg_100hits_clusterMETCR[4];
  TH1D *h_nStations5Seg_100hits_clusterMETCR[4];
  TH1D *h_nStations10Seg_100hits_clusterMETCR[4];
  TH1D *h_nStations1Seg_150hits_clusterMETCR[4];
  TH1D *h_nStations5Seg_150hits_clusterMETCR[4];
  TH1D *h_nStations10Seg_150hits_clusterMETCR[4];

  TH1D *h_nWheels1Seg_50hits_clusterMETCR[4];
  TH1D *h_nWheels5Seg_50hits_clusterMETCR[4];
  TH1D *h_nWheels10Seg_50hits_clusterMETCR[4];
  TH1D *h_nWheels1Seg_100hits_clusterMETCR[4];
  TH1D *h_nWheels5Seg_100hits_clusterMETCR[4];
  TH1D *h_nWheels10Seg_100hits_clusterMETCR[4];
  TH1D *h_nWheels1Seg_150hits_clusterMETCR[4];
  TH1D *h_nWheels5Seg_150hits_clusterMETCR[4];
  TH1D *h_nWheels10Seg_150hits_clusterMETCR[4];


  TH1D *h_nRPCStations1_50hits_clusterMETCR[4];
  TH1D *h_nRPCStations5_50hits_clusterMETCR[4];
  TH1D *h_nRPCStations10_50hits_clusterMETCR[4];
  TH1D *h_nRPCStations1_100hits_clusterMETCR[4];
  TH1D *h_nRPCStations5_100hits_clusterMETCR[4];
  TH1D *h_nRPCStations10_100hits_clusterMETCR[4];
  TH1D *h_nRPCStations1_150hits_clusterMETCR[4];
  TH1D *h_nRPCStations5_150hits_clusterMETCR[4];
  TH1D *h_nRPCStations10_150hits_clusterMETCR[4];

  TH1D *h_nRPCWheels1_50hits_clusterMETCR[4];
  TH1D *h_nRPCWheels5_50hits_clusterMETCR[4];
  TH1D *h_nRPCWheels10_50hits_clusterMETCR[4];
  TH1D *h_nRPCWheels1_100hits_clusterMETCR[4];
  TH1D *h_nRPCWheels5_100hits_clusterMETCR[4];
  TH1D *h_nRPCWheels10_100hits_clusterMETCR[4];
  TH1D *h_nRPCWheels1_150hits_clusterMETCR[4];
  TH1D *h_nRPCWheels5_150hits_clusterMETCR[4];
  TH1D *h_nRPCWheels10_150hits_clusterMETCR[4];

  TH1D *h_matchedRPCStation_MB3Cluster_clusterMETCR[4];
  TH1D *h_matchedRPCStation_MB4Cluster_clusterMETCR[4];
  TH1D *h_matchedRPCStation_MB4Cluster_spreadRPC_clusterMETCR[4];
  TH1D *h_matchedRPCStation_MB3Cluster_spreadRPC_clusterMETCR[4];

  TH1D *h_minSegmentDR_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegments_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsClusterStationMB2_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsClusterStationMB3_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsClusterStationMB4_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsOtherStations_invertedShowerVetoes[4];
  TH1D *h_nAlignedSegments_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStation_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsOuterStation_invertedShowerVetoes[4];
  TH1D *h_segmentAlignmentDeltaPhi_invertedShowerVetoes[4];
  TH1D *h_segmentAlignmentDeltaEta_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationMB2_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationMB3_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationMB4_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsOuterStationMB2_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsOuterStationMB3_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsOuterStationMB4_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedShowerVetoes[4];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedShowerVetoes[4];
  TH1D *h_muonVetoOutcomes_invertedShowerVetoes[4];
  TH1D *h_muonVetoOutcomesMB2_invertedShowerVetoes[4];
  TH1D *h_muonVetoOutcomesMB3_invertedShowerVetoes[4];
  TH1D *h_muonVetoOutcomesMB4_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomes_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesMB2_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesMB3_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesMB4_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassSegMuon_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassOneSegMuon_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[4];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[4];
  TH1D *h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[4];
  TH1D *h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[4];
  TH1D *h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[4];
  TH1D *h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimes_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMean_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesClusterStationMB2_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesClusterStationMB3_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesClusterStationMB4_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB2_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB3_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB4_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesInnerStationMB2_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesInnerStationMB3_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesInnerStationMB4_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB2_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB3_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB4_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesOuterStationMB2_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesOuterStationMB3_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimesOuterStationMB4_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB2_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB3_invertedShowerVetoes[4];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB4_invertedShowerVetoes[4];
  
  TH1D *h_minSegmentDR_invertedJetVeto[4];
  TH1D *h_nMatchedSegments_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsClusterStationMB2_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsClusterStationMB3_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsClusterStationMB4_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsOtherStations_invertedJetVeto[4];
  TH1D *h_nAlignedSegments_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStation_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsOuterStation_invertedJetVeto[4];
  TH1D *h_segmentAlignmentDeltaPhi_invertedJetVeto[4];
  TH1D *h_segmentAlignmentDeltaEta_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationMB2_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationMB3_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationMB4_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationMB2_invertedMB1_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationMB3_invertedMB1_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationMB4_invertedMB1_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsOuterStationMB2_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsOuterStationMB3_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsOuterStationMB4_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedJetVeto[4];
  TH1D *h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedJetVeto[4];
  TH1D *h_muonVetoOutcomes_invertedJetVeto[4];
  TH1D *h_muonVetoOutcomesMB2_invertedJetVeto[4];
  TH1D *h_muonVetoOutcomesMB3_invertedJetVeto[4];
  TH1D *h_muonVetoOutcomesMB4_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomes_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesMB2_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesMB3_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesMB4_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassSegMuon_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassOneSegMuon_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[4];
  TH1D *h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[4];
  TH1D *h_AdjacentMB1VetoOutcomes_invertedJetVeto[4];
  TH1D *h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[4];
  TH1D *h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[4];
  TH1D *h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimes_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMean_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesClusterStationMB2_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesClusterStationMB3_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesClusterStationMB4_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB2_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB3_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanClusterStationMB4_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesInnerStationMB2_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesInnerStationMB3_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesInnerStationMB4_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB2_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB3_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanInnerStationMB4_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesOuterStationMB2_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesOuterStationMB3_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimesOuterStationMB4_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB2_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB3_invertedJetVeto[4];
  TH1D *h_matchedSegmentTimeMeanOuterStationMB4_invertedJetVeto[4];
  
  TH1D *h_dPhiClusterMET_lowClusterSize_fullSelection_rpcCR[4];
  TH1D *h_dPhiClusterMET_highClusterSize_fullSelection_rpcCR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_rpcCR[4];

  TH1D *h_dtRechitClusterJetVetoPt[4];
  TH1D *h_dtRechitClusterMuonVetoPt[4];
  TH1D *h_dtRechitClusterMB1Veto[4];

  TH1D *h_efficiency[4];

  TH1D *h_dtRechitClusterNSegmentStation2_dPhiJetMET[4];
  TH1D *h_dtRechitClusterNSegmentStation3_dPhiJetMET[4];
  TH1D *h_dtRechitClusterNSegmentStation4_dPhiJetMET[4];
  TH1D *h_dtRechitClusterNStation_dPhiJetMET[4];
  TH1D *h_nDtRechitClusters_dPhiJetMET[4];
  TH1D *h_nDtRechitClustersVeto_dPhiJetMET[4];
  TH1D *h_dtRechitClustersDR_dPhiJetMET[4];
  TH1D *h_dtRechitClustersVetoDR_dPhiJetMET[4];
  TH1D *h_dtRechitClusterMaxStation_dPhiJetMET[4];
  TH1D *h_dtRechitClusterMaxStationRatio_dPhiJetMET[4];
  TH1D *h_dtRechitClusterMaxChamberRatio_dPhiJetMET[4];
  TH1D *h_dtRechitClusterNChamber_dPhiJetMET[4];
  TH1D *h_dtRechitClusterMaxChamber_dPhiJetMET[4];
  TH1D *h_dtRechitClusterX_dPhiJetMET[4];
  TH1D *h_dtRechitClusterY_dPhiJetMET[4];
  TH1D *h_dtRechitClusterZ_dPhiJetMET[4];
  TH1D *h_dtRechitClusterEta_dPhiJetMET[4];
  TH1D *h_dtRechitClusterPhi_dPhiJetMET[4];
  TH1D *h_dtRechitClusterTime_dPhiJetMET[4];
  TH1D *h_dtRechitClusterXSpread_dPhiJetMET[4];
  TH1D *h_dtRechitClusterYSpread_dPhiJetMET[4];
  TH1D *h_dtRechitClusterZSpread_dPhiJetMET[4];
  TH1D *h_dtRechitClusterEtaSpread_dPhiJetMET[4];
  TH1D *h_dtRechitClusterPhiSpread_dPhiJetMET[4];
  TH1D *h_dtRechitClusterTimeSpread_dPhiJetMET[4];
  TH1D *h_dtRechitClusterMajorAxis_dPhiJetMET[4];
  TH1D *h_dtRechitClusterMinorAxis_dPhiJetMET[4];

  TH1D *h_nRB1Match_dPhiJetMET[4];
  TH1D *h_nRB1Match_MB1Veto_dPhiJetMET[4];
  TH1D *h_nMB1MatchAdjacent_dPhiJetMET[4];
  TH1D *h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[4];

  TH1D *h_nMB1MatchPi2_dPhiClusterMET[4];

  TH1D *h_nRB1Match_dPhiClusterMET[4];
  TH1D *h_nRB1Match_MB1Veto_dPhiClusterMET[4];
  TH1D *h_nMB1MatchAdjacent_dPhiClusterMET[4];
  TH1D *h_nMB1MatchAdjacent_MB1Veto_dPhiClusterMET[4];
  TH1D *h_nMB1MatchAdjacentPi2_dPhiClusterMET[4];
  TH1D *h_nMB1MatchAdjacentPi2_MB1Veto_dPhiClusterMET[4];
  TH1D *h_nMB1MatchAdjacent0p8_dPhiClusterMET[4];
  TH1D *h_nMB1SegMatchAdjacent0p8_dPhiClusterMET[4];
  TH1D *h_nMB1MatchAdjacent0p8_MB1Veto_dPhiClusterMET[4];
  TH1D *h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[4];
  TH1D *h_nMB1MatchAdjacent0p8Pi2_MB1Veto_dPhiClusterMET[4];

  TH1D *h_nMB1MatchJet_invertedJetVeto[4];
  TH1D *h_nMB1MatchJet_MB2_invertedJetVeto[4];
  TH1D *h_nMB1MatchJet_MB3_invertedJetVeto[4];
  TH1D *h_nMB1MatchJet_MB4_invertedJetVeto[4];
  TH1D *h_matchedJetPt_invertedJetVeto[4];
  TH1D *h_clusterPhi_passMB1Veto_invertedJetVeto[4];
  TH1D *h_clusterEta_passMB1Veto_invertedJetVeto[4];
  TH2D *h_clusterEtaPhi_passMB1Veto_invertedJetVeto[4];
  TH2D *h_clusterEtaPhi_MB2_passMB1Veto_invertedJetVeto[4];
  TH2D *h_clusterEtaPhi_MB3_passMB1Veto_invertedJetVeto[4];
  TH2D *h_clusterEtaPhi_MB4_passMB1Veto_invertedJetVeto[4];
  TH2D *h_jetEtaPhi_passMB1Veto_invertedJetVeto[4];
  TH2D *h_jetEtaPhi_MB2_passMB1Veto_invertedJetVeto[4];
  TH2D *h_jetEtaPhi_MB3_passMB1Veto_invertedJetVeto[4];
  TH2D *h_jetEtaPhi_MB4_passMB1Veto_invertedJetVeto[4];
  TH1D *h_clusterPhi_MB2_passMB1Veto_invertedJetVeto[4];
  TH1D *h_clusterEta_MB2_passMB1Veto_invertedJetVeto[4];
  TH1D *h_clusterPhi_MB3_passMB1Veto_invertedJetVeto[4];
  TH1D *h_clusterEta_MB3_passMB1Veto_invertedJetVeto[4];
  TH1D *h_clusterPhi_MB4_passMB1Veto_invertedJetVeto[4];
  TH1D *h_clusterEta_MB4_passMB1Veto_invertedJetVeto[4];
  TH2D *h_matchedJetPt_nMB1Match_invertedJetVeto[4];
  TH2D *h_matchedJetPt_nMB1Match_invertedJetVetoNoJetMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB2_invertedJetVetoNoJetMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB3_invertedJetVetoNoJetMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB4_invertedJetVetoNoJetMET[4];
  TH2D *h_matchedJetPt_nMB1Match_invertedJetVetoNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_invertedJetVetoAllMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB2_invertedJetVetoNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB2_invertedJetVetoAllMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB3_invertedJetVetoNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB3_invertedJetVetoAllMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB4_invertedJetVetoNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB4_invertedJetVetoAllMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB2_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB3_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB4_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonLowClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonLowSizeNoClusterMET[4];
  TH2D *h_matchedJetPt_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMETNoSize[4];
  TH2D *h_MET_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonLowClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonLowSizeNoClusterMET[4];
  TH2D *h_matchedJetPt_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMETNoSize[4];
  TH2D *h_MET_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonLowClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonLowSizeNoClusterMET[4];
  TH2D *h_matchedJetPt_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMETNoSize[4];
  TH2D *h_MET_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMET[4];
  TH2D *h_matchedJetPt_nMB1Match_invertedJetVetoNoClusterMETNoJetMET[4];
  TH2D *h_matchedJetPt_nMB1Match_invertedJetVetoLooseMuonNoClusterMETNoJetMET[4];
  TH2D *h_matchedJetPt_nMB1Match_invertedJetVetoAllMuonNoClusterMETNoJetMET[4];
  TH2D *h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonNoClusterMETCluster80[4];
  TH2D *h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonNoClusterMETCluster80[4];
  TH2D *h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonNoClusterMETCluster80[4];

  TH1D *h_nMB1Match_veryLowdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[4];
  TH1D *h_nMB1Match_lowdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[4];
  TH1D *h_nMB1Match_highdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[4];
  TH1D *h_nMB1Match_veryLowdPhiClusterMET_highClusterSize_fullSelection_MB1CR[4];
  TH1D *h_nMB1Match_lowdPhiClusterMET_highClusterSize_fullSelection_MB1CR[4];
  TH1D *h_nMB1Match_highdPhiClusterMET_highClusterSize_fullSelection_MB1CR[4];

  TH1D *h_nMB1Match_veryLowdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[4];
  TH1D *h_nMB1Match_lowdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[4];
  TH1D *h_nMB1Match_highdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[4];
  TH1D *h_nMB1Match_veryLowdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[4];
  TH1D *h_nMB1Match_lowdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[4];
  TH1D *h_nMB1Match_highdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[4];
  TH2D *h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1HitsCR[4];
  TH2D *h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1Hits30CR[4];
  TH2D *h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[4];
  TH2D *h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1HitsCR[4];
  TH2D *h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1Hits30CR[4];
  TH2D *h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[4];

  TH1D *h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1CR[4];
  TH1D *h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1CR[4];
  TH1D *h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1CR[4];
  TH1D *h_dPhiClusterMET_lowClusterSize_fullSelection_MB1CR[4];
  TH1D *h_dPhiClusterMET_highClusterSize_fullSelection_MB1CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB1CR[4];

  TH1D *h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1HitsCR[4];
  TH1D *h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1HitsCR[4];
  TH1D *h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1HitsCR[4];
  TH1D *h_dPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[4];
  TH1D *h_dPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1HitsCR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB1HitsCR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB1HitsCR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB1HitsCR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB1HitsCR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB1HitsCR[4];

  TH1D *h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1Hits30CR[4];
  TH1D *h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1Hits30CR[4];
  TH1D *h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1Hits30CR[4];
  TH1D *h_dPhiClusterMET_lowClusterSize_fullSelection_MB1Hits30CR[4];
  TH1D *h_dPhiClusterMET_highClusterSize_fullSelection_MB1Hits30CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB3CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB4CR[4];

  TH1D *h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[4];
  TH1D *h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[4];
  TH1D *h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[4];
  TH1D *h_dPhiClusterMET_lowClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[4];
  TH1D *h_dPhiClusterMET_highClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[4];

  TH1D *h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB2CR[4];
  TH1D *h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB2CR[4];
  TH1D *h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB2CR[4];
  TH1D *h_dPhiClusterMET_lowClusterSize_fullSelection_MB2CR[4];
  TH1D *h_dPhiClusterMET_highClusterSize_fullSelection_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacent0p8MB1Cut_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCutNoRPC_MB2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCutNoRPC_MB2CR[4];

  TH1D *h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB2withMB1CR[4];
  TH1D *h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB2withMB1CR[4];
  TH1D *h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB2withMB1CR[4];
  TH1D *h_dPhiClusterMET_lowClusterSize_fullSelection_MB2withMB1CR[4];
  TH1D *h_dPhiClusterMET_highClusterSize_fullSelection_MB2withMB1CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB2withMB1CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB2withMB1CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB2withMB1CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB2withMB1CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB2withMB1CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB2withMB1CR[4];

  TH1D *h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1or2CR[4];
  TH1D *h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1or2CR[4];
  TH1D *h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1or2CR[4];
  TH1D *h_dPhiClusterMET_lowClusterSize_fullSelection_MB1or2CR[4];
  TH1D *h_dPhiClusterMET_highClusterSize_fullSelection_MB1or2CR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1or2CR[4];

  TH1D *h_nMatchedHitsMB2_fullSelection_SRMB3[4];
  TH1D *h_nMatchedHitsMB4_fullSelection_SRMB3[4];
  TH1D *h_nMatchedHitsMB2and4_lowClusterSize_fullSelection_SRMB3[4];
  TH1D *h_nMatchedHitsMB2and4_highClusterSize_fullSelection_SRMB3[4];
  TH1D *h_nMatchedHitsMB2_fullSelection_SRMB4[4];
  TH1D *h_nMatchedHitsMB3_fullSelection_SRMB4[4];
  TH1D *h_nMatchedHitsMB2and3_lowClusterSize_fullSelection_SRMB4[4];
  TH1D *h_nMatchedHitsMB2and3_highClusterSize_fullSelection_SRMB4[4];
  TH1D *h_nMatchedHitsMB3_fullSelection_MB2CR[4];
  TH1D *h_nMatchedHitsMB4_fullSelection_MB2CR[4];
  TH1D *h_nMatchedHitsMB3and4_lowClusterSize_fullSelection_MB2CR[4];
  TH1D *h_nMatchedHitsMB3and4_highClusterSize_fullSelection_MB2CR[4];

  TH1D *h_nMatchedHitsMB2_fullSelectionWithAdjacentMB1Cut_SRMB3[4];
  TH1D *h_nMatchedHitsMB4_fullSelectionWithAdjacentMB1Cut_SRMB3[4];
  TH1D *h_nMatchedHitsMB2and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3[4];
  TH1D *h_nMatchedHitsMB2and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3[4];
  TH1D *h_nMatchedHitsMB2_fullSelectionWithAdjacentMB1Cut_SRMB4[4];
  TH1D *h_nMatchedHitsMB3_fullSelectionWithAdjacentMB1Cut_SRMB4[4];
  TH1D *h_nMatchedHitsMB2and3_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4[4];
  TH1D *h_nMatchedHitsMB2and3_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4[4];
  TH1D *h_nMatchedHitsMB3_fullSelectionWithAdjacentMB1Cut_MB2CR[4];
  TH1D *h_nMatchedHitsMB4_fullSelectionWithAdjacentMB1Cut_MB2CR[4];
  TH1D *h_nMatchedHitsMB3and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR[4];
  TH1D *h_nMatchedHitsMB3and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR[4];

  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SRMB3[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SRMB4[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacent0p8MB1Cut_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SRMB3[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SRMB4[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SRMB3[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SRMB4[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCutNoRPC_SR[4];
  TH2D *h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCutNoRPC_SR[4];

  TH1D *h_dtRechitClusterSize_SRhighNHF[4];
  TH1D *h_MET_SRhighNHF[4];
  TH1D *h_dPhiClusterMET_SRhighNHF[4];
  TH1D *h_dtRechitClusterSize_MB2CRhighNHF[4];
  TH1D *h_MET_MB2CRhighNHF[4];
  TH1D *h_dPhiClusterMET_MB2CRhighNHF[4];
  TH1D *h_dtRechitClusterSize_MB1CRhighNHF[4];
  TH1D *h_MET_MB1CRhighNHF[4];
  TH1D *h_dPhiClusterMET_MB1CRhighNHF[4];
  TH1D *h_dtRechitClusterSize_MB2withMB1CRhighNHF[4];
  TH1D *h_MET_MB2withMB1CRhighNHF[4];
  TH1D *h_dPhiClusterMET_MB2withMB1CRhighNHF[4];  
  TH1D *h_dtRechitClusterSize_MB1HitsCRhighNHF[4];
  TH1D *h_MET_MB1HitsCRhighNHF[4];
  TH1D *h_dPhiClusterMET_MB1HitsCRhighNHF[4];

  TH1D *h_dtRechitClusterSize_SRlowNHF[4];
  TH1D *h_MET_SRlowNHF[4];
  TH1D *h_dPhiClusterMET_SRlowNHF[4];
  TH1D *h_dtRechitClusterSize_MB2CRlowNHF[4];
  TH1D *h_MET_MB2CRlowNHF[4];
  TH1D *h_dPhiClusterMET_MB2CRlowNHF[4];
  TH1D *h_dtRechitClusterSize_MB1CRlowNHF[4];
  TH1D *h_MET_MB1CRlowNHF[4];
  TH1D *h_dPhiClusterMET_MB1CRlowNHF[4];
  TH1D *h_dtRechitClusterSize_MB2withMB1CRlowNHF[4];
  TH1D *h_MET_MB2withMB1CRlowNHF[4];
  TH1D *h_dPhiClusterMET_MB2withMB1CRlowNHF[4];  
  TH1D *h_dtRechitClusterSize_MB1HitsCRlowNHF[4];
  TH1D *h_MET_MB1HitsCRlowNHF[4];
  TH1D *h_dPhiClusterMET_MB1HitsCRlowNHF[4];
  
  TH1D *h_nCAClusters[4];
  TH1D *h_sizeCACluster[4];
  TH1D *h_radiusCACluster[4];
  TH1D *h_zCACluster[4];
  TH1D *h_phiCACluster[4];

  TH1D *h_matchedJetNHF[4];
  TH1D *h_matchedJetNEMF[4];
  TH1D *h_matchedJetCHF[4];
  TH1D *h_matchedJetCEMF[4];
  TH1D *h_matchedJetNConst[4];
  TH1D *h_matchedJetChargedMult[4];

  TH1D *h_jetOverlap10GeV[4];
  TH1D *h_jetOverlap20GeV[4];
  TH1D *h_muonOverlap10GeVLoose[4];
  TH1D *h_muonOverlap10GeVTight[4];

  TH1D *h_nDtRechitClusters[4];
  TH1D *h_dtRechitClusterSize[4];
  TH1D *h_dtRechitClusterR[4];
  TH1D *h_dtRechitClusterZ[4];
  TH1D *h_dtRechitClusterPhi[4];

  std::vector<TString> selsStrings = {"oneCluster_lowNHF","oneVetoCluster_lowNHF","oneCluster_highNHF","oneVetoCluster_highNHF"};
  std::map<const TString, TH1D*> runNumHists;
  std::map<const TString, TH1D*> etaClusterHists;
  std::map<const TString, TH1D*> phiClusterHists;
  std::map<const TString, TH2D*> etaPhiClusterHists;
      

  TH1D *h_MET_highNHF[4];
  TH1D *h_MET_oneCluster_highNHF[4];
  TH1D *h_MET_oneVetoCluster_highNHF[4];
  TH1D *h_MET_lowNHF[4];
  TH1D *h_MET_oneCluster_lowNHF[4];
  TH1D *h_MET_oneVetoCluster_lowNHF[4];

  TH1D *h_leadJetPt_highNHF[4];
  TH1D *h_leadJetPt_oneCluster_highNHF[4];
  TH1D *h_leadJetPt_lowNHF[4];
  TH1D *h_leadJetPt_oneCluster_lowNHF[4];
  TH1D *h_leadJetPt_oneVetoCluster_highNHF[4];
  TH1D *h_leadJetPt_oneVetoCluster_lowNHF[4];
  TH1D *h_leadJetPt_MB1CR[4];
  TH1D *h_leadJetPt_MB2withMB1CR[4];
  TH1D *h_leadJetPt_MB1HitsCR[4];
  TH1D *h_leadJetPt_MB2CR[4];
  TH1D *h_leadJetPt_SR[4];
  
  TH1D *h_leadJetPtMET_highNHF[4];
  TH1D *h_leadJetPtMET_oneCluster_highNHF[4];
  TH1D *h_leadJetPtMET_lowNHF[4];
  TH1D *h_leadJetPtMET_oneCluster_lowNHF[4];
  TH1D *h_leadJetPtMET_oneVetoCluster_highNHF[4];
  TH1D *h_leadJetPtMET_oneVetoCluster_lowNHF[4];
  TH1D *h_leadJetPtMET_MB1CR[4];
  TH1D *h_leadJetPtMET_MB2withMB1CR[4];
  TH1D *h_leadJetPtMET_MB1HitsCR[4];
  TH1D *h_leadJetPtMET_MB2CR[4];
  TH1D *h_leadJetPtMET_SR[4];

  TH1D *h_rpcBxMedian_MB1CR[4];
  TH1D *h_rpcBxMedian_MB2CR[4];
  TH1D *h_rpcBxMedian_MB2withMB1CR[4];
  TH1D *h_rpcBxMedian_SR[4];
  TH1D *h_rpcBxMedian_MB1HitsCR[4];

  TH1D *h_jetNeutralHadronicEnergyFraction_passJetMET[4];
  TH1D *h_jetNeutralHadronicEnergyFraction_passStationsWheels[4];
  TH1D *h_jetNeutralHadronicEnergyFraction_SR[4];
  TH1D *h_jetChargedHadronicEnergyFraction_passJetMET[4];
  TH1D *h_jetChargedHadronicEnergyFraction_passStationsWheels[4];
  TH1D *h_jetChargedHadronicEnergyFraction_SR[4];
  TH1D *h_jetNeutralEMEnergyFraction_passJetMET[4];
  TH1D *h_jetNeutralEMEnergyFraction_passStationsWheels[4];
  TH1D *h_jetNeutralEMEnergyFraction_SR[4];
  TH1D *h_jetChargedEMEnergyFraction_passJetMET[4];
  TH1D *h_jetChargedEMEnergyFraction_passStationsWheels[4];
  TH1D *h_jetChargedEMEnergyFraction_SR[4];
  TH1D *h_jetChargedEMEnergyFraction_MB2CR[4];
  TH1D *h_jetChargedEMEnergyFraction_MB2withMB1CR[4];
  TH1D *h_jetChargedEMEnergyFraction_MB1CR[4];
  TH1D *h_jetChargedEMEnergyFraction_MB1HitsCR[4];
  TH1D *h_jetNeutralEMEnergyFraction_MB2CR[4];
  TH1D *h_jetNeutralEMEnergyFraction_MB2withMB1CR[4];
  TH1D *h_jetNeutralEMEnergyFraction_MB1CR[4];
  TH1D *h_jetNeutralEMEnergyFraction_MB1HitsCR[4];
  TH1D *h_jetChargedHadronicEnergyFraction_MB2CR[4];
  TH1D *h_jetChargedHadronicEnergyFraction_MB2withMB1CR[4];
  TH1D *h_jetChargedHadronicEnergyFraction_MB1CR[4];
  TH1D *h_jetChargedHadronicEnergyFraction_MB1HitsCR[4];
  TH1D *h_jetNeutralHadronicEnergyFraction_MB2CR[4];
  TH1D *h_jetNeutralHadronicEnergyFraction_MB2withMB1CR[4];
  TH1D *h_jetNeutralHadronicEnergyFraction_MB1CR[4];
  TH1D *h_jetNeutralHadronicEnergyFraction_MB1HitsCR[4];
  TH1D *h_leadingJetChargedEMEnergyFraction_passStationsWheels[4];
  TH1D *h_leadingJetChargedEMEnergyFraction_SR[4];
  TH1D *h_leadingJetNeutralEMEnergyFraction_passStationsWheels[4];
  TH1D *h_leadingJetNeutralEMEnergyFraction_SR[4];
  TH1D *h_leadingJetChargedHadronicEnergyFraction_passStationsWheels[4];
  TH1D *h_leadingJetChargedHadronicEnergyFraction_SR[4];
  TH1D *h_leadingJetNeutralHadronicEnergyFraction_passStationsWheels[4];
  TH1D *h_leadingJetNeutralHadronicEnergyFraction_SR[4];
  
  TH1D *h_invertedLeadJetVeto_failedMETFilters_SR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_MB2CR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR[4];
  TH1D *h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR[4];

  TH2D *h_clusterSize_looseMuonVetoResult_ABCD[4];
  TH2D *h_clusterSize_recoMuonVetoResult_ABCD[4];

  TH1D *h_clusterSize_looseMuonVeto_D[4];
  TH1D *h_clusterSize_looseMuonVeto_C[4];
  TH1D *h_clusterSizeMuonMB1_looseMuonVeto_C[4];
  TH1D *h_clusterSize_looseMuonVeto_B[4];
  TH1D *h_clusterSize_looseMuonVeto_A[4];
  TH1D *h_matchedMB1Hits_looseMuonVeto_D[4];
  TH1D *h_matchedMB1Hits_looseMuonVeto_C[4];
  TH1D *h_matchedMB1Hits_looseMuonVeto_B[4];
  TH1D *h_matchedMB1Hits_looseMuonVeto_A[4];
  TH1D *h_matchedMB1Hits_invertedLooseMuonVeto_D[4];
  TH1D *h_matchedMB1Hits_invertedLooseMuonVeto_C[4];
  TH1D *h_matchedMB1Hits_invertedLooseMuonVeto_B[4];
  TH1D *h_matchedMB1Hits_invertedLooseMuonVeto_A[4];
  TH1D *h_matchedMB1Hits_invertedRecoMuonVeto_D[4];
  TH1D *h_matchedMB1Hits_invertedRecoMuonVeto_C[4];
  TH1D *h_matchedMB1Hits_invertedRecoMuonVeto_B[4];
  TH1D *h_matchedMB1Hits_invertedRecoMuonVeto_A[4];
  TH1D *h_matchedMB1Segments_looseMuonVeto_D[4];
  TH1D *h_matchedMB1Segments_looseMuonVeto_C[4];
  TH1D *h_matchedMB1Segments_looseMuonVeto_B[4];
  TH1D *h_matchedMB1Segments_looseMuonVeto_A[4];
  TH1D *h_minMB1SegmentDR_looseMuonVeto_D[4];
  TH1D *h_minMB1SegmentDR_looseMuonVeto_C[4];
  TH1D *h_minMB1SegmentDR_looseMuonVeto_B[4];
  TH1D *h_minMB1SegmentDR_looseMuonVeto_A[4];
  TH1D *h_matchedMB1HitsRatio_looseMuonVeto_D[4];
  TH1D *h_matchedMB1HitsRatio_looseMuonVeto_C[4];
  TH1D *h_matchedMB1HitsRatio_looseMuonVeto_B[4];
  TH1D *h_matchedMB1HitsRatio_looseMuonVeto_A[4];
  TH1D *h_matchedMB1Hits_looseMuonVetoDoubleInverted_D[4];
  TH1D *h_matchedMB1Hits_looseMuonVetoDoubleInverted_C[4];
  TH1D *h_matchedMB1Hits_looseMuonVetoDoubleInverted_B[4];
  TH1D *h_matchedMB1Hits_looseMuonVetoDoubleInverted_A[4];
  TH1D *h_matchedMB1Segments_looseMuonVetoDoubleInverted_D[4];
  TH1D *h_matchedMB1Segments_looseMuonVetoDoubleInverted_C[4];
  TH1D *h_matchedMB1Segments_looseMuonVetoDoubleInverted_B[4];
  TH1D *h_matchedMB1Segments_looseMuonVetoDoubleInverted_A[4];
  TH1D *h_minMB1SegmentDR_looseMuonVetoDoubleInverted_D[4];
  TH1D *h_minMB1SegmentDR_looseMuonVetoDoubleInverted_C[4];
  TH1D *h_minMB1SegmentDR_looseMuonVetoDoubleInverted_B[4];
  TH1D *h_minMB1SegmentDR_looseMuonVetoDoubleInverted_A[4];
  TH1D *h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_D[4];
  TH1D *h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_C[4];
  TH1D *h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_B[4];
  TH1D *h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_A[4];

  TH1D *h_dtRechitClusterPhi_invertedMB1AC_MB3_10GeVJet[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1AC_MB3_jetVeto20[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1AC_MB3_jetVeto10[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1AC_MB4_10GeVJet[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1AC_MB4_jetVeto20[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1AC_MB4_jetVeto10[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1AC_10GeVJet[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1AC_jetVeto20[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1AC_jetVeto10[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1A_10GeVJet[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1A_jetVeto20[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1A_jetVeto10[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1C_10GeVJet[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1C_jetVeto20[4];
  TH1D *h_dtRechitClusterPhi_invertedMB1C_jetVeto10[4];
  TH1D *h_dtRechitClusterPhi_AC_10GeVJet[4];
  TH1D *h_dtRechitClusterPhi_AC_jetVeto20[4];
  TH1D *h_dtRechitClusterPhi_AC_jetVeto10[4];
  

  Double_t dtR = 0.0;
  Double_t dPhi_tmp = 0.0;
  Double_t dPhi_min = 0.0;
  Double_t dPhi_min_invertedJet = 0.0;
  Double_t dPhiClusterRPC = 0.0;
  Double_t dZClusterRPC = 0.0;
  Double_t dPhiClusterMET_max = 0.0;
  Double_t dPhiClusterMET = 0.0;
  vector<Int_t> rpcBx = {};
  vector<Int_t> rpcMatchStation = {};
  Int_t rpcSpread = 0;
  Double_t rpcMedian = 0;
  Double_t chargedHadFraction_mindPhi = 0.0;
  Double_t chargedEMFraction_mindPhi = 0.0;
  Double_t neutralHadFraction_mindPhi = 0.0;
  Double_t neutralEMFraction_mindPhi = 0.0;

  Double_t matchedLooseIDMuonPt = 0.0;

  Double_t matchedJetPhi = -999.;
  Double_t matchedJetEta = -999.;
  Double_t matchedJetNHF = -999.;
  Double_t matchedJetCHF = -999.;
  Double_t matchedJetNEF = -999.;
  Double_t matchedJetCEF = -999.;
  
  vector<Double_t> clusterEta = {};
  vector<Double_t> clusterPhi = {};
  vector<Int_t> clusterSize = {};
  Int_t clusterSizeTotal = 0;
  Int_t nClustersVeto_dPhiJetMET = 0;

  //vector<fastjet::PseudoJet> dtPoints;
  //vector<fastjet::PseudoJet> clustersCA;
  //vector<fastjet::PseudoJet> constituents;
  Int_t nCAClusters = 0;
  Int_t SRyieldA_CA = 0;
  Int_t SRyieldB_CA = 0;
  Int_t SRyieldC_CA = 0;
  Int_t SRyieldMB2A_CA = 0;
  Int_t SRyieldMB2B_CA = 0;
  Int_t SRyieldMB2C_CA = 0;
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

  Int_t evtNum = 0;
  Bool_t HLT = false;

  Bool_t goodInvertedJet = false;

  Bool_t passMuon = false;
  Bool_t passMuonLoose = false;
  Bool_t passMuon_alt = false;
  Bool_t passMB1 = false;
  Bool_t passJet = false;

  Bool_t overlapJet10GeV = false;
  Bool_t overlapJet20GeV = false;
  Bool_t overlapMuon10GeVLoose = false;
  Bool_t overlapMuon10GeVTight = false;
  
  Bool_t passMET = false;
  Bool_t passOneJet = false;
  Bool_t passJetMET = false;
  Bool_t passNHFJet50 = false;
  Bool_t passNHFJetLead = false;
  Bool_t passStations25 = false;
  Bool_t passWheels25 = false;
  Bool_t passJetVeto = false;
  Bool_t passMuonVeto = false;
  Bool_t passMuonVetoLoose = false;
  Bool_t passRpcMatch = false;
  Bool_t passRpcMatchMB1CRwithNHFLead = false;
  Bool_t passRpcMatchMB1HitsCR = false;
  Bool_t passRpcMatchMB1HitsCRwithAdjacent = false;
  Bool_t passRpcMatchMB1HitsCRwithAdjacent0p8 = false;
  Bool_t passRpcMatchMB1HitsCRwithNHF50 = false;
  Bool_t passRpcMatchMB1HitsCRinvertedNHF50 = false;
  Bool_t passRpcMatchMB1HitsCRwithNHFLead = false;
  Bool_t passRpcMatchMB1HitsCRinvertedNHFLead = false;
  Bool_t passRpcMatchMB1Hits30CR = false;
  Bool_t passRpcMatchMB1Hits30MaxMB2CR = false;
  Bool_t passRpcMatchMB1Hits30MaxMB3CR = false;
  Bool_t passRpcMatchMB1Hits30MaxMB4CR = false;
  Bool_t passRpcMatchMB1HitsNoMB1ClusterCR = false;
  Bool_t passRpcMatchMB1or2CR = false;
  Bool_t passRpcMatchMB2CR = false;
  Bool_t passRpcMatchMB2CRwithAdjacent = false;
  Bool_t passRpcMatchMB2CRwithAdjacent0p8 = false;
  Bool_t passRpcMatchMB2CRwithOther = false;
  Bool_t passRpcMatchMB2CRwithNHF50 = false;
  Bool_t passRpcMatchMB2CRinvertedNHF50 = false;
  Bool_t passRpcMatchMB2CRwithNHFLead = false;
  Bool_t passRpcMatchMB2CRinvertedNHFLead = false;
  Bool_t noRpcMatchMB2CRwithNHFLead = false;
  Bool_t passRpcMatchMB2withMB1CR = false;
  Bool_t passRpcMatchMB2withMB1CRwithAdjacent = false;
  Bool_t passRpcMatchMB2withMB1CRwithNHF50 = false;
  Bool_t passRpcMatchMB2withMB1CRinvertedNHF50 = false;
  Bool_t passRpcMatchMB2withMB1CRwithNHFLead = false;
  Bool_t passRpcMatchMB2withMB1CRinvertedNHFLead = false;
  Bool_t passMB2CRwithInvertedNHFLeadNoRPC = false;
  Bool_t passMB1Veto = false;
  Bool_t passMaxStation = false;
  Bool_t passClusterMET = false;
  Bool_t passRPCMatch = false;
  Bool_t passRPCSpread = false;
  Bool_t passRPCBx = false;
  Bool_t passNoVetoCluster = false;
  Bool_t passClusterSize = false;
  Bool_t passMB1CR = false;
  Int_t clusterSizeMB1CR = 0;
  Double_t dPhiClusterMETMB1CR = 0.0;
  Int_t clusterSizeMB1CRwithNHFLead = 0;
  Double_t dPhiClusterMETMB1CRwithNHFLead = 0.0;
  Int_t clusterSizeMB1HitsCR = 0;
  Double_t dPhiClusterMETMB1HitsCR = 0.0;
  Int_t clusterSizeMB1HitsCRwithAdjacent = 0;
  Double_t dPhiClusterMETMB1HitsCRwithAdjacent = 0.0;
  Int_t clusterSizeMB1HitsCRwithNHF50 = 0;
  Double_t dPhiClusterMETMB1HitsCRwithNHF50 = 0.0;
  Int_t clusterSizeMB1HitsCRinvertedNHF50 = 0;
  Double_t dPhiClusterMETMB1HitsCRinvertedNHF50 = 0.0;
  Int_t clusterSizeMB1HitsCRwithNHFLead = 0;
  Double_t dPhiClusterMETMB1HitsCRwithNHFLead = 0.0;
  Int_t clusterSizeMB1HitsCRinvertedNHFLead = 0;
  Double_t dPhiClusterMETMB1HitsCRinvertedNHFLead = 0.0;
  Int_t clusterSizeMB1Hits30CR = 0;
  Double_t dPhiClusterMETMB1Hits30CR = 0.0;
  Int_t clusterSizeMB1Hits30MaxMB2CR = 0;
  Double_t dPhiClusterMETMB1Hits30MaxMB2CR = 0.0;
  Int_t clusterSizeMB1Hits30MaxMB3CR = 0;
  Double_t dPhiClusterMETMB1Hits30MaxMB3CR = 0.0;
  Int_t clusterSizeMB1Hits30MaxMB4CR = 0;
  Double_t dPhiClusterMETMB1Hits30MaxMB4CR = 0.0;
  Int_t clusterSizeMB1HitsNoMB1ClusterCR = 0;
  Double_t dPhiClusterMETMB1HitsNoMB1ClusterCR = 0.0;
  Int_t clusterSizeMB1or2CR = 0;
  Double_t dPhiClusterMETMB1or2CR = 0.0;
  Int_t clusterSizeMB2CR = 0;
  Double_t dPhiClusterMETMB2CR = 0.0;
  Int_t clusterSizeMB2CRwithAdjacent = 0;
  Double_t dPhiClusterMETMB2CRwithAdjacent = 0.0;
  Int_t clusterSizeMB2CRwithAdjacent0p8 = 0;
  Double_t dPhiClusterMETMB2CRwithAdjacent0p8 = 0.0;
  Int_t clusterSizeMB2CRwithOther = 0;
  Double_t dPhiClusterMETMB2CRwithOther = 0.0;
  Int_t clusterSizeMB2CRwithNHF50 = 0;
  Double_t dPhiClusterMETMB2CRwithNHF50 = 0.0;
  Int_t clusterSizeMB2CRinvertedNHF50 = 0;
  Double_t dPhiClusterMETMB2CRinvertedNHF50 = 0.0;
  Int_t clusterSizeMB2CRwithNHFLead = 0;
  Double_t dPhiClusterMETMB2CRwithNHFLead = 0.0;
  Int_t clusterSizeMB2CRinvertedNHFLead = 0;
  Double_t dPhiClusterMETMB2CRinvertedNHFLead = 0.0;
  Int_t clusterSizeMB2CRwithNHFLeadNoRPC = 0;
  Double_t dPhiClusterMETMB2CRwithNHFLeadNoRPC = 0.0;
  Int_t clusterSizeMB2CRwithInvertedNHFLeadNoRPC = 0;
  Double_t dPhiClusterMETMB2CRwithInvertedNHFLeadNoRPC = 0.0;
  Int_t clusterSizeMB2withMB1CR = 0;
  Double_t dPhiClusterMETMB2withMB1CR = 0.0;
  Int_t clusterSizeMB2withMB1CRwithAdjacent = 0;
  Double_t dPhiClusterMETMB2withMB1CRwithAdjacent = 0.0;
  Int_t clusterSizeMB2withMB1CRwithNHF50 = 0;
  Double_t dPhiClusterMETMB2withMB1CRwithNHF50 = 0.0;
  Int_t clusterSizeMB2withMB1CRinvertedNHF50 = 0;
  Double_t dPhiClusterMETMB2withMB1CRinvertedNHF50 = 0.0;
  Int_t clusterSizeMB2withMB1CRwithNHFLead = 0;
  Double_t dPhiClusterMETMB2withMB1CRwithNHFLead = 0.0;
  Int_t clusterSizeMB2withMB1CRinvertedNHFLead = 0;
  Double_t dPhiClusterMETMB2withMB1CRinvertedNHFLead = 0.0;

  Int_t nPassMET = 0;
  Int_t nPassOneJet = 0;
  Int_t nPassNoVetoCluster = 0;
  Int_t nPassJetMET = 0;
  Int_t nPassStations25 = 0;
  Int_t nPassWheels25 = 0;
  Int_t nPassMB1CR = 0;
  Int_t nPassJetVeto = 0;
  Int_t nPassMuonVeto = 0;
  Int_t nPassMuonVetoLoose = 0;
  Int_t nPassRpcMatch = 0;
  Int_t nPassHighClusterMETMB1CR = 0;
  Int_t nPassHighClusterMETHighClusterSizeMB1CR = 0;
  Int_t nPassHighClusterMETLowClusterSizeMB1CR = 0;
  Int_t nPassLowClusterMETMB1CR = 0;
  Int_t nPassLowClusterMETHighClusterSizeMB1CR = 0;
  Int_t nPassLowClusterMETLowClusterSizeMB1CR = 0;
  Int_t nPassHighClusterMETMB1HitsCR = 0;
  Int_t nPassHighClusterMETHighClusterSizeMB1HitsCR = 0;
  Int_t nPassHighClusterMETLowClusterSizeMB1HitsCR = 0;
  Int_t nPassLowClusterMETMB1HitsCR = 0;
  Int_t nPassLowClusterMETHighClusterSizeMB1HitsCR = 0;
  Int_t nPassLowClusterMETLowClusterSizeMB1HitsCR = 0;

  Int_t nPassCluster_JetMET_StationsWheels_InvertedJet = 0;
  Int_t nPassMaxStation_InvertedJet = 0;
  Int_t nPassMuonVeto_InvertedJet = 0;
  Int_t nPassRpcMatch_InvertedJet = 0;
  Int_t nPassInvertedJetVeto_InvertedJet = 0;
  Int_t nPassdPhiClusterMET_InvertedJet = 0;
  Int_t nPassClusterSize_InvertedJet = 0;
  Int_t nPassLooseMuonVeto_InvertedJet = 0;
  Int_t nPassRpcMatch_InvertedJetLooseMuon = 0;
  Int_t nPassInvertedJetVeto_InvertedJetLooseMuon = 0;
  Int_t nPassdPhiClusterMET_InvertedJetLooseMuon = 0;
  Int_t nPassClusterSize_InvertedJetLooseMuon = 0;

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
  Bool_t passJetMET_rpcCR = false;
  Bool_t passMaxStation_rpcCR = false;
  Bool_t passLepton_rpcCR = false;
  Bool_t pass50Hits_rpcCR = false;
  Bool_t pass25Hits_rpcCR = false;
  
  Bool_t passFullPlus = false;

  Int_t nPassClusterCR = 0;
  Int_t nPassNoVeto_clusterCR = 0;
  Int_t nPassFullVeto_clusterCR = 0;
  Int_t nPassRPCMatch_clusterCR = 0;
  Int_t nPassRPCSpread_clusterCR = 0;
  Int_t nPassRPCBx_clusterCR = 0;
  Int_t nPassMaxStation_clusterCR = 0;
  Int_t nPassLepton_clusterCR = 0;
  Int_t nPass50Hits_clusterCR = 0;
  Int_t nPass25Hits_clusterCR = 0;

  Int_t nPassNoVeto = 0;  
  Int_t nPassFullVeto_rpcCR = 0;
  Int_t nPassRPCCR = 0;
  Int_t nPassClusterMET_rpcCR = 0;
  Int_t nPassJetMET_rpcCR = 0;
  Int_t nPassMaxStation_rpcCR = 0;
  Int_t nPassLepton_rpcCR = 0;
  Int_t nPass50Hits_rpcCR = 0;
  Int_t nPass25Hits_rpcCR = 0;
  Int_t nPassFullPlus = 0;

  Int_t passLowClusterMET_rpcCR = false;
  Int_t passHighClusterMET_rpcCR = false;
  Int_t passLowClusterMET_LowClusterSize = false;
  Int_t passHighClusterMET_LowClusterSize = false;
  Int_t passHighClusterMET_HighClusterSize = false;

  Int_t SRyieldA = 0;
  Int_t SRyieldB = 0;
  Int_t SRyieldC = 0;
  Bool_t passABCD = false;
  Bool_t passABCDwithAdjacent = false;
  Bool_t passABCDwithAdjacent0p8 = false;
  Bool_t passABCDwithOther = false;
  Bool_t passABCDwithNHF50 = false;
  Bool_t passABCDinvertedNHF50 = false;
  Bool_t passABCDwithNHFLead = false;
  Bool_t passABCDinvertedNHFLead = false;
  Bool_t passABCDwithNHFLeadNoRPC = false;
  Bool_t passABCDwithInvertedNHFLeadNoRPC = false;
  Int_t clusterSizeSR = 0;
  Double_t dPhiClusterMETSR = 0.0;
  Int_t clusterSizeSRwithAdjacent = 0;
  Double_t dPhiClusterMETSRwithAdjacent = 0.0;
  Int_t clusterSizeSRwithAdjacent0p8 = 0;
  Double_t dPhiClusterMETSRwithAdjacent0p8 = 0.0;
  Int_t clusterSizeSRwithOther = 0;
  Double_t dPhiClusterMETSRwithOther = 0.0;
  Int_t clusterSizeSRwithNHF50 = 0;
  Double_t dPhiClusterMETSRwithNHF50 = 0.0;
  Int_t clusterSizeSRinvertedNHF50 = 0;
  Double_t dPhiClusterMETSRinvertedNHF50 = 0.0;
  Int_t clusterSizeSRwithNHFLead = 0;
  Double_t dPhiClusterMETSRwithNHFLead = 0.0;
  Int_t clusterSizeSRinvertedNHFLead = 0;
  Double_t dPhiClusterMETSRinvertedNHFLead = 0.0;
  Int_t clusterSizeSRwithNHFLeadNoRPC = 0;
  Double_t dPhiClusterMETSRwithNHFLeadNoRPC = 0.0;
  Int_t clusterSizeSRwithInvertedNHFLeadNoRPC = 0;
  Double_t dPhiClusterMETSRwithInvertedNHFLeadNoRPC = 0.0;
  Int_t maxStationSR = 0;
  Int_t maxStationSRwithAdjacent = 0;
  Int_t maxStationSRwithOther = 0;

  Int_t invertedMB1_MB2Cluster_5MB1 = 0;
  
  ofstream eventListInvertedMB1_phiSpike;
  eventListInvertedMB1_phiSpike.open("events/invertedMB1_phiSpike.txt");

  ofstream eventListInvertedMB1_10GeVJetsA;
  eventListInvertedMB1_10GeVJetsA.open("events/invertedMB1A_10GeVJetEvents.txt");

  ofstream eventListInvertedMB1_10GeVJetsC;
  eventListInvertedMB1_10GeVJetsC.open("events/invertedMB1C_10GeVJetEvents.txt");

  ofstream eventListABCD_10GeVJetsA;
  eventListABCD_10GeVJetsA.open("events/A_10GeVJetEvents.txt");

  ofstream eventListABCD_10GeVJetsC;
  eventListABCD_10GeVJetsC.open("events/C_10GeVJetEvents.txt");

  ofstream eventListInvertedJet;
  eventListInvertedJet.open("events/invertedJet_events.txt");
  /*
  ofstream eventListInvertedJetPassMB1;
  eventListInvertedJetPassMB1.open("events/invertedJetPassMB1_events.txt");

  ofstream eventListMB1;
  eventListMB1.open("events/MB1CR_tailEvents.txt");

  ofstream eventListRPC;
  eventListRPC.open("events/rpcCR_tailEvents.txt");
  ofstream eventListNoRPC;
  eventListNoRPC.open("events/noRpcMatch_events.txt");
  */
  Int_t event100_rpc = 0;
  Int_t event150_rpc = 0;
  //ofstream eventListClusterMET;
  //eventListClusterMET.open("events/clusterMETCR_tailEvents.txt");
  Int_t event100_clusterMET = 0;
  Int_t event150_clusterMET = 0;
  /*
  ofstream eventListSRB;
  eventListSRB.open("events/regionB_events.txt");
  ofstream eventListABC;
  eventListABC.open("events/regionABC_events.txt");
  ofstream eventListMB2CR;
  eventListMB2CR.open("events/MB2CR_highClusterSize.txt");
  ofstream eventListAdjacentMB1;
  eventListAdjacentMB1.open("events/AdjacentMB1Hits_events.txt");

  ofstream eventListInvertedJetIDSR;
  eventListInvertedJetIDSR.open("events/invertedJetID_SRevents.txt");
  ofstream eventListInvertedJetIDABC;
  eventListInvertedJetIDABC.open("events/invertedJetID_ABCevents.txt");

  ofstream eventListInvertedMB1LooseMuon;
  eventListInvertedMB1LooseMuon.open("events/invertedMB1_looseMuonVeto.txt");
  */
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

  Double_t minMB1SegmentDR = 999.;
  Double_t minSegmentDR = 999.;
  Int_t nAlignedSegments = 0;
  Double_t meanSegTime = -999.;
  Double_t meanMB1SegTime = -999.;
  Double_t meanMB2SegTime = -999.;
  Double_t meanMB3SegTime = -999.;
  Double_t meanMB4SegTime = -999.;

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

  Int_t rpcStation = 99;
  Int_t rpcWheel = 99;
  Int_t dtStation = 99;
  Int_t dtWheel = 99;
  Int_t maxClusterSize = 0;
  Int_t maxGoodClusterSize = 0;
  Int_t hitsMB1 = 0;
  Int_t nMB1MatchCluster = 0;
  Int_t nMB2MatchCluster = 0;
  Int_t nMB3MatchCluster = 0;
  Int_t nMB4MatchCluster = 0;
  Int_t nMB1SegMatchCluster = 0;
  Int_t nMB2SegMatchCluster = 0;
  Int_t nMB3SegMatchCluster = 0;
  Int_t nMB4SegMatchCluster = 0;
  Int_t nMB1MatchClusterMB1CR = 0;
  Int_t nMB1MatchClusterMB1HitsCR = 0;
  Int_t nMB1MatchClusterMB1Hits30CR = 0;
  Int_t nMB1MatchClusterMB1HitsNoMB1ClusterCR = 0;
  Int_t nMB1MatchClusterMB1or2CR = 0;
  Int_t nMB1MatchClusterMB2CR = 0;
  Int_t nMB1MatchClusterMB2withMB1CR = 0;
  Int_t nMB1MatchClusterMB2withMB1CRwithNHFLead = 0;
  Int_t nMB1MatchClusterAdjacentPlus = 0;
  Int_t nMB1MatchClusterAdjacentMinus = 0;
  Int_t nMB1MatchClusterAdjacent0p8Plus = 0;
  Int_t nMB1MatchClusterAdjacent0p8Minus = 0;
  Int_t nMB1SegMatchClusterAdjacent0p8Plus = 0;
  Int_t nMB1SegMatchClusterAdjacent0p8Minus = 0;
  Int_t nMB1MatchPi2 = 0;
  Int_t nMB1MatchPi2AdjacentPlus = 0;
  Int_t nMB1MatchPi2AdjacentMinus = 0;
  Int_t nMB1MatchPi2Adjacent0p8Plus = 0;
  Int_t nMB1MatchPi2Adjacent0p8Minus = 0;
  Int_t nRB1MatchCluster = 0;
  Int_t invertedJetMB1 = 0;

  TRandom3 *rand = new TRandom3();
  Int_t pmRand = 0;

  for(Int_t itr_year=0; itr_year<3; itr_year++){
    TString year = years[itr_year];

    name = "h_dtRechitClusterNSegmentStation2_dPhiJetMET_"+year;
    h_dtRechitClusterNSegmentStation2_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,500);
    
    name = "h_dtRechitClusterNSegmentStation3_dPhiJetMET_"+year;
    h_dtRechitClusterNSegmentStation3_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,500);
    
    name = "h_dtRechitClusterNSegmentStation4_dPhiJetMET_"+year;
    h_dtRechitClusterNSegmentStation4_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterMaxStation_dPhiJetMET_"+year;
    h_dtRechitClusterMaxStation_dPhiJetMET[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_dtRechitClusterNStation_dPhiJetMET_"+year;
    h_dtRechitClusterNStation_dPhiJetMET[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nDtRechitClusters_dPhiJetMET_"+year;
    h_nDtRechitClusters_dPhiJetMET[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nDtRechitClustersVeto_dPhiJetMET_"+year;
    h_nDtRechitClustersVeto_dPhiJetMET[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_dtRechitClustersDR_dPhiJetMET_"+year;
    h_dtRechitClustersDR_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,2);

    name = "h_dtRechitClustersVetoDR_dPhiJetMET_"+year;
    h_dtRechitClustersVetoDR_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,2);

    name = "h_dtRechitClusterMaxStationRatio_dPhiJetMET_"+year;
    h_dtRechitClusterMaxStationRatio_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_dtRechitClusterMaxChamberRatio_dPhiJetMET_"+year;
    h_dtRechitClusterMaxChamberRatio_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_dtRechitClusterNChamber_dPhiJetMET_"+year;
    h_dtRechitClusterNChamber_dPhiJetMET[itr_year] = new TH1D(name,"",6,0,6);
    
    name = "h_dtRechitClusterMaxChamber_dPhiJetMET_"+year;
    h_dtRechitClusterMaxChamber_dPhiJetMET[itr_year] = new TH1D(name,"",6,0,6);
    
    name = "h_dtRechitClusterX_dPhiJetMET_"+year;
    h_dtRechitClusterX_dPhiJetMET[itr_year] = new TH1D(name,"",100,-800,800);
    
    name = "h_dtRechitClusterY_dPhiJetMET_"+year;
    h_dtRechitClusterY_dPhiJetMET[itr_year] = new TH1D(name,"",100,-800,800);
    
    name = "h_dtRechitClusterZ_dPhiJetMET_"+year;
    h_dtRechitClusterZ_dPhiJetMET[itr_year] = new TH1D(name,"",100,-600,600);
    
    name = "h_dtRechitClusterEta_dPhiJetMET_"+year;
    h_dtRechitClusterEta_dPhiJetMET[itr_year] = new TH1D(name,"",50,-1.5,1.5);
    
    name = "h_dtRechitClusterPhi_dPhiJetMET_"+year;
    h_dtRechitClusterPhi_dPhiJetMET[itr_year] = new TH1D(name,"",100,-3.5,3.5);
    
    name = "h_dtRechitClusterTime_dPhiJetMET_"+year;
    h_dtRechitClusterTime_dPhiJetMET[itr_year] = new TH1D(name,"",50,300,800);
    
    name = "h_dtRechitClusterXSpread_dPhiJetMET_"+year;
    h_dtRechitClusterXSpread_dPhiJetMET[itr_year] = new TH1D(name,"",100,0,500);
    
    name = "h_dtRechitClusterYSpread_dPhiJetMET_"+year;
    h_dtRechitClusterYSpread_dPhiJetMET[itr_year] = new TH1D(name,"",100,0,500);
    
    name = "h_dtRechitClusterZSpread_dPhiJetMET_"+year;
    h_dtRechitClusterZSpread_dPhiJetMET[itr_year] = new TH1D(name,"",100,0,140);
    
    name = "h_dtRechitClusterEtaSpread_dPhiJetMET_"+year;
    h_dtRechitClusterEtaSpread_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,0.5);
    
    name = "h_dtRechitClusterPhiSpread_dPhiJetMET_"+year;
    h_dtRechitClusterPhiSpread_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_dtRechitClusterTimeSpread_dPhiJetMET_"+year;
    h_dtRechitClusterTimeSpread_dPhiJetMET[itr_year] = new TH1D(name,"",100,0,150);
    
    name = "h_dtRechitClusterMajorAxis_dPhiJetMET_"+year;
    h_dtRechitClusterMajorAxis_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_dtRechitClusterMinorAxis_dPhiJetMET_"+year;
    h_dtRechitClusterMinorAxis_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,1);
    

    name = "h_dtRechitClusterPhi_fullSelectionMB2_rpcCR_"+year;
    h_dtRechitClusterPhi_fullSelectionMB2_rpcCR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterPhi_fullSelectionMB34_rpcCR_"+year;
    h_dtRechitClusterPhi_fullSelectionMB34_rpcCR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_metPhi_fullSelectionMB2_rpcCR_"+year;
    h_metPhi_fullSelectionMB2_rpcCR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_metPhi_fullSelectionMB34_rpcCR_"+year;
    h_metPhi_fullSelectionMB34_rpcCR[itr_year] = new TH1D(name,"",70,0,3.5);


    name = "h_dtRechitClusterSize_dPhiJetMETLow_rpcCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiJetMETLow_rpcCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);
    
    name = "h_dtRechitClusterSize_dPhiJetMETHigh_rpcCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiJetMETHigh_rpcCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);
    
    name = "h_dtRechitClusterSize_dPhiClusterMETLow_rpcCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiClusterMETLow_rpcCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiClusterMETHigh_rpcCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiClusterMETHigh_rpcCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiClusterMETLow_rpcSpreadCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiClusterMETLow_rpcSpreadCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiClusterMETHigh_rpcSpreadCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiClusterMETHigh_rpcSpreadCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiJetMETLow_rpcCRnoLepton_"+year;
    h_dtRechitClusterSize_dPhiJetMETLow_rpcCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiJetMETHigh_rpcCRnoLepton_"+year;
    h_dtRechitClusterSize_dPhiJetMETHigh_rpcCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiClusterMETLow_rpcCRnoLepton_"+year;
    h_dtRechitClusterSize_dPhiClusterMETLow_rpcCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiClusterMETHigh_rpcCRnoLepton_"+year;
    h_dtRechitClusterSize_dPhiClusterMETHigh_rpcCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);



    name = "h_dtRechitClusterSize_dPhiClusterMETLow_jetMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiClusterMETLow_jetMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiClusterMETHigh_jetMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiClusterMETHigh_jetMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcMatch_jetMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcMatch_jetMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcNoMatch_jetMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcNoMatch_jetMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiClusterMETLow_jetMETCRnoLepton_"+year;
    h_dtRechitClusterSize_dPhiClusterMETLow_jetMETCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiClusterMETHigh_jetMETCRnoLepton_"+year;
    h_dtRechitClusterSize_dPhiClusterMETHigh_jetMETCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcMatch_jetMETCRnoLepton_"+year;
    h_dtRechitClusterSize_rpcMatch_jetMETCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcNoMatch_jetMETCRnoLepton_"+year;
    h_dtRechitClusterSize_rpcNoMatch_jetMETCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);



    name = "h_dtRechitClusterSize_dPhiJetMETLow_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiJetMETLow_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiJetMETHigh_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_dPhiJetMETHigh_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcGood_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcGood_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcBad_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcBad_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcMatch_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcMatch_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcNoMatch_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcNoMatch_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcSpread_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcSpread_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcNoSpread_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcNoSpread_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcNegativeSpread_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcNegativeSpread_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcNegativeNoSpread_clusterMETCRmuonVeto_"+year;
    h_dtRechitClusterSize_rpcNegativeNoSpread_clusterMETCRmuonVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiJetMETLow_clusterMETCRnoLepton_"+year;
    h_dtRechitClusterSize_dPhiJetMETLow_clusterMETCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_dPhiJetMETHigh_clusterMETCRnoLepton_"+year;
    h_dtRechitClusterSize_dPhiJetMETHigh_clusterMETCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcMatch_clusterMETCRnoLepton_"+year;
    h_dtRechitClusterSize_rpcMatch_clusterMETCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_rpcNoMatch_clusterMETCRnoLepton_"+year;
    h_dtRechitClusterSize_rpcNoMatch_clusterMETCRnoLepton[itr_year] = new TH1D(name,"",50,0,500);


    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_rpcCR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_rpcCR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);


    name = "h_dtRechitClusterSize_fullSelection_rpcCR_"+year;
    h_dtRechitClusterSize_fullSelection_rpcCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_lowClusterMET_fullSelection_rpcCR_"+year;
    h_dtRechitClusterSize_lowClusterMET_fullSelection_rpcCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highClusterMET_fullSelection_rpcCR_"+year;
    h_dtRechitClusterSize_highClusterMET_fullSelection_rpcCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_fullSelection_clusterMETCR_"+year;
    h_dtRechitClusterSize_fullSelection_clusterMETCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_goodRPC_fullSelection_clusterMETCR_"+year;
    h_dtRechitClusterSize_goodRPC_fullSelection_clusterMETCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_badRPC_fullSelection_clusterMETCR_"+year;
    h_dtRechitClusterSize_badRPC_fullSelection_clusterMETCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_lowClusterMET_lowClusterSize_SR_"+year;
    h_dtRechitClusterSize_lowClusterMET_lowClusterSize_SR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highClusterMET_lowClusterSize_SR_"+year;
    h_dtRechitClusterSize_highClusterMET_lowClusterSize_SR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highClusterMET_highClusterSize_SR_"+year;
    h_dtRechitClusterSize_highClusterMET_highClusterSize_SR[itr_year] = new TH1D(name,"",50,0,500);



    name = "h_dPhiClusterRPC_fullVeto_"+year;
    h_dPhiClusterRPC_fullVeto[itr_year] = new TH1D(name,"",73,-0.15,3.5);

    name = "h_dZClusterRPC_fullVeto_"+year;
    h_dZClusterRPC_fullVeto[itr_year] = new TH1D(name,"",101,-5,500);

    name = "h_rpcSpread_invertedJetVeto_"+year;
    h_rpcSpread_invertedJetVeto[itr_year] = new TH1D(name,"",10,0,10);

    name = "h_rpcSpread_invertedJetVeto_muonVeto_"+year;
    h_rpcSpread_invertedJetVeto_muonVeto[itr_year] = new TH1D(name,"",10,0,10);

    name = "h_rpcSpread_fullVeto_"+year;
    h_rpcSpread_fullVeto[itr_year] = new TH1D(name,"",10,0,10);

    name = "h_rpcSpread_fullVeto_negBx_"+year;
    h_rpcSpread_fullVeto_negBx[itr_year] = new TH1D(name,"",10,0,10);
    
    

    name = "h_nRPCMatched_fullVeto_clusterMETCR_"+year;
    h_nRPCMatched_fullVeto_clusterMETCR[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_rpcSpread_fullVeto_clusterMETCR_"+year;
    h_rpcSpread_fullVeto_clusterMETCR[itr_year] = new TH1D(name,"",10,0,10);

    name = "h_rpcBx_fullVeto_clusterMETCR_"+year;
    h_rpcBx_fullVeto_clusterMETCR[itr_year] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_dPhiJetMET_fullVeto_clusterMETCR_"+year;
    h_dPhiJetMET_fullVeto_clusterMETCR[itr_year] = new TH1D(name,"",35,0,3.5);

    name = "h_dtRechitClusterMaxStation_fullVeto_clusterMETCR_"+year;
    h_dtRechitClusterMaxStation_fullVeto_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCMatched_Nminus1_clusterMETCR_"+year;
    h_nRPCMatched_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_rpcSpread_Nminus1_clusterMETCR_"+year;
    h_rpcSpread_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",10,0,10);

    name = "h_rpcBx_Nminus1_clusterMETCR_"+year;
    h_rpcBx_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_dPhiJetMET_Nminus1_clusterMETCR_"+year;
    h_dPhiJetMET_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",35,0,3.5);

    name = "h_dtRechitClusterMaxStation_Nminus1_clusterMETCR_"+year;
    h_dtRechitClusterMaxStation_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_nMB1Matched_Nminus1_clusterMETCR_"+year;
    h_nMB1Matched_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",60,0,60);

    name = "h_jetVetoPt_Nminus1_clusterMETCR_"+year;
    h_jetVetoPt_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",20,0,200);

    name = "h_muonVetoPt_Nminus1_clusterMETCR_"+year;
    h_muonVetoPt_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",20,0,200);

    name = "h_muonLooseIDVetoPt_Nminus1_clusterMETCR_"+year;
    h_muonLooseIDVetoPt_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",20,0,200);

    name = "h_jetVetoPt_Nminus1_MB1CR_"+year;
    h_jetVetoPt_Nminus1_MB1CR[itr_year] = new TH1D(name,"",20,0,200);

    name = "h_muonVetoPt_Nminus1_MB1CR_"+year;
    h_muonVetoPt_Nminus1_MB1CR[itr_year] = new TH1D(name,"",20,0,200);

    name = "h_muonLooseIDVetoPt_Nminus1_MB1CR_"+year;
    h_muonLooseIDVetoPt_Nminus1_MB1CR[itr_year] = new TH1D(name,"",20,0,200);

    name = "h_nMB1MatchedAdjacent_Nminus1_MB1CR_"+year;
    h_nMB1MatchedAdjacent_Nminus1_MB1CR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMB1MatchedAdjacent_Nminus1_clusterMETCR_"+year;
    h_nMB1MatchedAdjacent_Nminus1_clusterMETCR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nRPCMatched_Nminus1_MB1CR_"+year;
    h_nRPCMatched_Nminus1_MB1CR[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_dPhiClusterMET_Nminus1_MB1CR_"+year;
    h_dPhiClusterMET_Nminus1_MB1CR[itr_year] = new TH1D(name,"",35,0,3.5);

    name = "h_dPhiClusterMET_fullVeto_rpcCR_"+year;
    h_dPhiClusterMET_fullVeto_rpcCR[itr_year] = new TH1D(name,"",35,0,3.5);

    name = "h_dPhiJetMET_fullVeto_rpcCR_"+year;
    h_dPhiJetMET_fullVeto_rpcCR[itr_year] = new TH1D(name,"",35,0,3.5);

    name = "h_dtRechitClusterMaxStation_fullVeto_rpcCR_"+year;
    h_dtRechitClusterMaxStation_fullVeto_rpcCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_dPhiClusterMET_Nminus1_rpcCR_"+year;
    h_dPhiClusterMET_Nminus1_rpcCR[itr_year] = new TH1D(name,"",35,0,3.5);

    name = "h_dPhiJetMET_Nminus1_rpcCR_"+year;
    h_dPhiJetMET_Nminus1_rpcCR[itr_year] = new TH1D(name,"",35,0,3.5);

    name = "h_dtRechitClusterMaxStation_Nminus1_rpcCR_"+year;
    h_dtRechitClusterMaxStation_Nminus1_rpcCR[itr_year] = new TH1D(name,"",5,0,5);



    name = "h_nStations1_50hits_clusterMETCR_"+year;
    h_nStations1_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations1_100hits_clusterMETCR_"+year;
    h_nStations1_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations1_150hits_clusterMETCR_"+year;
    h_nStations1_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations25_50hits_clusterMETCR_"+year;
    h_nStations25_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations25_100hits_clusterMETCR_"+year;
    h_nStations25_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations25_150hits_clusterMETCR_"+year;
    h_nStations25_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations50_50hits_clusterMETCR_"+year;
    h_nStations50_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations50_100hits_clusterMETCR_"+year;
    h_nStations50_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations50_150hits_clusterMETCR_"+year;
    h_nStations50_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nWheels1_50hits_clusterMETCR_"+year;
    h_nWheels1_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels1_100hits_clusterMETCR_"+year;
    h_nWheels1_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels1_150hits_clusterMETCR_"+year;
    h_nWheels1_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels25_50hits_clusterMETCR_"+year;
    h_nWheels25_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels25_100hits_clusterMETCR_"+year;
    h_nWheels25_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels25_150hits_clusterMETCR_"+year;
    h_nWheels25_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels50_50hits_clusterMETCR_"+year;
    h_nWheels50_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels50_100hits_clusterMETCR_"+year;
    h_nWheels50_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels50_150hits_clusterMETCR_"+year;
    h_nWheels50_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);


    name = "h_nDtSegsStation_50hits_clusterMETCR_"+year;
    h_nDtSegsStation_50hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsStation_100hits_clusterMETCR_"+year;
    h_nDtSegsStation_100hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsStation_150hits_clusterMETCR_"+year;
    h_nDtSegsStation_150hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsWheel_50hits_clusterMETCR_"+year;
    h_nDtSegsWheel_50hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsWheel_100hits_clusterMETCR_"+year;
    h_nDtSegsWheel_100hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsWheel_150hits_clusterMETCR_"+year;
    h_nDtSegsWheel_150hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsChamber_50hits_clusterMETCR_"+year;
    h_nDtSegsChamber_50hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsChamber_100hits_clusterMETCR_"+year;
    h_nDtSegsChamber_100hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegsChamber_150hits_clusterMETCR_"+year;
    h_nDtSegsChamber_150hits_clusterMETCR[itr_year] = new TH1D(name,"",15,0,15);

    name = "h_nDtSegs_50hits_clusterMETCR_"+year;
    h_nDtSegs_50hits_clusterMETCR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nDtSegs_100hits_clusterMETCR_"+year;
    h_nDtSegs_100hits_clusterMETCR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nDtSegs_150hits_clusterMETCR_"+year;
    h_nDtSegs_150hits_clusterMETCR[itr_year] = new TH1D(name,"",50,0,50);


    name = "h_nStations1Seg_50hits_clusterMETCR_"+year;
    h_nStations1Seg_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations1Seg_100hits_clusterMETCR_"+year;
    h_nStations1Seg_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations1Seg_150hits_clusterMETCR_"+year;
    h_nStations1Seg_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations5Seg_50hits_clusterMETCR_"+year;
    h_nStations5Seg_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations5Seg_100hits_clusterMETCR_"+year;
    h_nStations5Seg_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations5Seg_150hits_clusterMETCR_"+year;
    h_nStations5Seg_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations10Seg_50hits_clusterMETCR_"+year;
    h_nStations10Seg_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations10Seg_100hits_clusterMETCR_"+year;
    h_nStations10Seg_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nStations10Seg_150hits_clusterMETCR_"+year;
    h_nStations10Seg_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nWheels1Seg_50hits_clusterMETCR_"+year;
    h_nWheels1Seg_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels1Seg_100hits_clusterMETCR_"+year;
    h_nWheels1Seg_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels1Seg_150hits_clusterMETCR_"+year;
    h_nWheels1Seg_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels5Seg_50hits_clusterMETCR_"+year;
    h_nWheels5Seg_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels5Seg_100hits_clusterMETCR_"+year;
    h_nWheels5Seg_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels5Seg_150hits_clusterMETCR_"+year;
    h_nWheels5Seg_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels10Seg_50hits_clusterMETCR_"+year;
    h_nWheels10Seg_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels10Seg_100hits_clusterMETCR_"+year;
    h_nWheels10Seg_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nWheels10Seg_150hits_clusterMETCR_"+year;
    h_nWheels10Seg_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);


    name = "h_nRPCStations1_50hits_clusterMETCR_"+year;
    h_nRPCStations1_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations1_100hits_clusterMETCR_"+year;
    h_nRPCStations1_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations1_150hits_clusterMETCR_"+year;
    h_nRPCStations1_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations5_50hits_clusterMETCR_"+year;
    h_nRPCStations5_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations5_100hits_clusterMETCR_"+year;
    h_nRPCStations5_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations5_150hits_clusterMETCR_"+year;
    h_nRPCStations5_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations10_50hits_clusterMETCR_"+year;
    h_nRPCStations10_50hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations10_100hits_clusterMETCR_"+year;
    h_nRPCStations10_100hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCStations10_150hits_clusterMETCR_"+year;
    h_nRPCStations10_150hits_clusterMETCR[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nRPCWheels1_50hits_clusterMETCR_"+year;
    h_nRPCWheels1_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels1_100hits_clusterMETCR_"+year;
    h_nRPCWheels1_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels1_150hits_clusterMETCR_"+year;
    h_nRPCWheels1_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels5_50hits_clusterMETCR_"+year;
    h_nRPCWheels5_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels5_100hits_clusterMETCR_"+year;
    h_nRPCWheels5_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels5_150hits_clusterMETCR_"+year;
    h_nRPCWheels5_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels10_50hits_clusterMETCR_"+year;
    h_nRPCWheels10_50hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels10_100hits_clusterMETCR_"+year;
    h_nRPCWheels10_100hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_nRPCWheels10_150hits_clusterMETCR_"+year;
    h_nRPCWheels10_150hits_clusterMETCR[itr_year] = new TH1D(name,"",6,0,6);


    name = "h_matchedRPCStation_MB3Cluster_clusterMETCR_"+year;
    h_matchedRPCStation_MB3Cluster_clusterMETCR[itr_year] = new TH1D(name,"",7,0,7);

    name = "h_matchedRPCStation_MB4Cluster_clusterMETCR_"+year;
    h_matchedRPCStation_MB4Cluster_clusterMETCR[itr_year] = new TH1D(name,"",7,0,7);

    name = "h_matchedRPCStation_MB3Cluster_spreadRPC_clusterMETCR_"+year;
    h_matchedRPCStation_MB3Cluster_spreadRPC_clusterMETCR[itr_year] = new TH1D(name,"",7,0,7);

    name = "h_matchedRPCStation_MB4Cluster_spreadRPC_clusterMETCR_"+year;
    h_matchedRPCStation_MB4Cluster_spreadRPC_clusterMETCR[itr_year] = new TH1D(name,"",7,0,7);


    name = "h_minSegmentDR_invertedShowerVetoes_"+year;
    h_minSegmentDR_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,0,2);

    name = "h_nMatchedSegments_invertedShowerVetoes_"+year;
    h_nMatchedSegments_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB2_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsClusterStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB3_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsClusterStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsClusterStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOtherStations_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsOtherStations_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nAlignedSegments_invertedShowerVetoes_"+year;
    h_nAlignedSegments_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStation_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStation_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationMB2_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationMB3_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStation_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsOuterStation_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB2_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsOuterStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB3_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsOuterStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsOuterStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_segmentAlignmentDeltaPhi_invertedShowerVetoes_"+year;
    h_segmentAlignmentDeltaPhi_invertedShowerVetoes[itr_year] = new TH1D(name,"",70,0,3.5);
    
    name = "h_segmentAlignmentDeltaEta_invertedShowerVetoes_"+year;
    h_segmentAlignmentDeltaEta_invertedShowerVetoes[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedShowerVetoes_"+year;
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_muonVetoOutcomes_invertedShowerVetoes_"+year;
    h_muonVetoOutcomes_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB2_invertedShowerVetoes_"+year;
    h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB3_invertedShowerVetoes_"+year;
    h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB4_invertedShowerVetoes_"+year;
    h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomes_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomes_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesMB2_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB3_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB4_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuon_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassSegMuon_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuon_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassOneSegMuon_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes_"+year;
    h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomes_invertedShowerVetoes_"+year;
    h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes_"+year;
    h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes_"+year;
    h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes_"+year;
    h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_matchedSegmentTimes_invertedShowerVetoes_"+year;
    h_matchedSegmentTimes_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMean_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMean_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB2_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesClusterStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB2_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanClusterStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB3_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesClusterStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB3_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanClusterStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB4_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesClusterStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB4_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanClusterStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB2_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesInnerStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB2_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanInnerStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB3_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesInnerStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB3_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanInnerStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB4_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesInnerStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB4_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanInnerStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB2_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesOuterStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB2_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanOuterStationMB2_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB3_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesOuterStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB3_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanOuterStationMB3_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB4_invertedShowerVetoes_"+year;
    h_matchedSegmentTimesOuterStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB4_invertedShowerVetoes_"+year;
    h_matchedSegmentTimeMeanOuterStationMB4_invertedShowerVetoes[itr_year] = new TH1D(name,"",50,-100,100);


    
    name = "h_minSegmentDR_invertedJetVeto_"+year;
    h_minSegmentDR_invertedJetVeto[itr_year] = new TH1D(name,"",100,0,2);

    name = "h_nMatchedSegments_invertedJetVeto_"+year;
    h_nMatchedSegments_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsClusterStationMB2_invertedJetVeto_"+year;
    h_nMatchedSegmentsClusterStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsClusterStationMB3_invertedJetVeto_"+year;
    h_nMatchedSegmentsClusterStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsClusterStationMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsClusterStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOtherStations_invertedJetVeto_"+year;
    h_nMatchedSegmentsOtherStations_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nAlignedSegments_invertedJetVeto_"+year;
    h_nAlignedSegments_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStation_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStation_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationMB2_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationMB3_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationMB2_invertedMB1_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationMB2_invertedMB1_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationMB3_invertedMB1_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationMB3_invertedMB1_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationMB4_invertedMB1_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationMB4_invertedMB1_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStation_invertedJetVeto_"+year;
    h_nMatchedSegmentsOuterStation_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB2_invertedJetVeto_"+year;
    h_nMatchedSegmentsOuterStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB3_invertedJetVeto_"+year;
    h_nMatchedSegmentsOuterStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsOuterStationMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsOuterStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_segmentAlignmentDeltaPhi_invertedJetVeto_"+year;
    h_segmentAlignmentDeltaPhi_invertedJetVeto[itr_year] = new TH1D(name,"",70,0,3.5);
    
    name = "h_segmentAlignmentDeltaEta_invertedJetVeto_"+year;
    h_segmentAlignmentDeltaEta_invertedJetVeto[itr_year] = new TH1D(name,"",70,0,3.5);

name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassSegMuonMB2_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassSegMuonMB3_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassSegMuonMB4_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
  
    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);
    
    name = "h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedJetVeto_"+year;
    h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedJetVeto[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_muonVetoOutcomes_invertedJetVeto_"+year;
    h_muonVetoOutcomes_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB2_invertedJetVeto_"+year;
    h_muonVetoOutcomesMB2_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB3_invertedJetVeto_"+year;
    h_muonVetoOutcomesMB3_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_muonVetoOutcomesMB4_invertedJetVeto_"+year;
    h_muonVetoOutcomesMB4_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomes_invertedJetVeto_"+year;
    h_MB1VetoOutcomes_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB2_invertedJetVeto_"+year;
    h_MB1VetoOutcomesMB2_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB3_invertedJetVeto_"+year;
    h_MB1VetoOutcomesMB3_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesMB4_invertedJetVeto_"+year;
    h_MB1VetoOutcomesMB4_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuon_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassSegMuon_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_MB1VetoOutcomesPassOneSegMuon_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassOneSegMuon_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);
    
    name = "h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto_"+year;
    h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomes_invertedJetVeto_"+year;
    h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto_"+year;
    h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto_"+year;
    h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto_"+year;
    h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_matchedSegmentTimes_invertedJetVeto_"+year;
    h_matchedSegmentTimes_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMean_invertedJetVeto_"+year;
    h_matchedSegmentTimeMean_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB2_invertedJetVeto_"+year;
    h_matchedSegmentTimesClusterStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB2_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanClusterStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB3_invertedJetVeto_"+year;
    h_matchedSegmentTimesClusterStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB3_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanClusterStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesClusterStationMB4_invertedJetVeto_"+year;
    h_matchedSegmentTimesClusterStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanClusterStationMB4_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanClusterStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB2_invertedJetVeto_"+year;
    h_matchedSegmentTimesInnerStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB2_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanInnerStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB3_invertedJetVeto_"+year;
    h_matchedSegmentTimesInnerStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB3_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanInnerStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesInnerStationMB4_invertedJetVeto_"+year;
    h_matchedSegmentTimesInnerStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanInnerStationMB4_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanInnerStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB2_invertedJetVeto_"+year;
    h_matchedSegmentTimesOuterStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB2_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanOuterStationMB2_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB3_invertedJetVeto_"+year;
    h_matchedSegmentTimesOuterStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB3_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanOuterStationMB3_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);

    name = "h_matchedSegmentTimesOuterStationMB4_invertedJetVeto_"+year;
    h_matchedSegmentTimesOuterStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",100,-200,200);

    name = "h_matchedSegmentTimeMeanOuterStationMB4_invertedJetVeto_"+year;
    h_matchedSegmentTimeMeanOuterStationMB4_invertedJetVeto[itr_year] = new TH1D(name,"",50,-100,100);



    name = "h_dPhiClusterMET_lowClusterSize_fullSelection_rpcCR_"+year;
    h_dPhiClusterMET_lowClusterSize_fullSelection_rpcCR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dPhiClusterMET_highClusterSize_fullSelection_rpcCR_"+year;
    h_dPhiClusterMET_highClusterSize_fullSelection_rpcCR[itr_year] = new TH1D(name,"",70,0,3.5);

    
    name = "h_dtRechitClusterJetVetoPt_"+year;
    h_dtRechitClusterJetVetoPt[itr_year] = new TH1D(name,"",20,0,200);

    name = "h_dtRechitClusterMuonVetoPt_"+year;
    h_dtRechitClusterMuonVetoPt[itr_year] = new TH1D(name,"",20,0,200);

    name = "h_dtRechitClusterMB1Veto_"+year;
    h_dtRechitClusterMB1Veto[itr_year] = new TH1D(name,"",60,0,60);
    

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1HitsCR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1HitsCR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB1HitsCR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB1HitsCR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB1HitsCR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB1HitsCR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB1HitsCR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB1HitsCR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB1HitsCR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB1HitsCR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB1HitsCR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB1HitsCR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1HitsCR_"+year;
    h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1HitsCR_"+year;
    h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1HitsCR_"+year;
    h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR_"+year;
    h_dPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR_"+year;
    h_dPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",70,0,3.5);


    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB3CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB3CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB4CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB4CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1Hits30CR_"+year;
    h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1Hits30CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1Hits30CR_"+year;
    h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1Hits30CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1Hits30CR_"+year;
    h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1Hits30CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dPhiClusterMET_lowClusterSize_fullSelection_MB1Hits30CR_"+year;
    h_dPhiClusterMET_lowClusterSize_fullSelection_MB1Hits30CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dPhiClusterMET_highClusterSize_fullSelection_MB1Hits30CR_"+year;
    h_dPhiClusterMET_highClusterSize_fullSelection_MB1Hits30CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR_"+year;
    h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR_"+year;
    h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR_"+year;
    h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dPhiClusterMET_lowClusterSize_fullSelection_MB1HitsNoMB1ClusterCR_"+year;
    h_dPhiClusterMET_lowClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dPhiClusterMET_highClusterSize_fullSelection_MB1HitsNoMB1ClusterCR_"+year;
    h_dPhiClusterMET_highClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[itr_year] = new TH1D(name,"",70,0,3.5);


    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB1CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB1CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1CR_"+year;
    h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1CR_"+year;
    h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1CR_"+year;
    h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dPhiClusterMET_lowClusterSize_fullSelection_MB1CR_"+year;
    h_dPhiClusterMET_lowClusterSize_fullSelection_MB1CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dPhiClusterMET_highClusterSize_fullSelection_MB1CR_"+year;
    h_dPhiClusterMET_highClusterSize_fullSelection_MB1CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacent0p8MB1Cut_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacent0p8MB1Cut_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCutNoRPC_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCutNoRPC_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCutNoRPC_MB2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCutNoRPC_MB2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB2CR_"+year;
    h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB2CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB2CR_"+year;
    h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB2CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dPhiClusterMET_lowClusterSize_fullSelection_MB2CR_"+year;
    h_dPhiClusterMET_lowClusterSize_fullSelection_MB2CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dPhiClusterMET_highClusterSize_fullSelection_MB2CR_"+year;
    h_dPhiClusterMET_highClusterSize_fullSelection_MB2CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB2withMB1CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB2withMB1CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB2withMB1CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB2withMB1CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB2withMB1CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB2withMB1CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB2withMB1CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB2withMB1CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB2withMB1CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB2withMB1CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB2withMB1CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB2withMB1CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB2withMB1CR_"+year;
    h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB2withMB1CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB2withMB1CR_"+year;
    h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB2withMB1CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dPhiClusterMET_lowClusterSize_fullSelection_MB2withMB1CR_"+year;
    h_dPhiClusterMET_lowClusterSize_fullSelection_MB2withMB1CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dPhiClusterMET_highClusterSize_fullSelection_MB2withMB1CR_"+year;
    h_dPhiClusterMET_highClusterSize_fullSelection_MB2withMB1CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1or2CR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1or2CR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1or2CR_"+year;
    h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1or2CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1or2CR_"+year;
    h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1or2CR[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dPhiClusterMET_lowClusterSize_fullSelection_MB1or2CR_"+year;
    h_dPhiClusterMET_lowClusterSize_fullSelection_MB1or2CR[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dPhiClusterMET_highClusterSize_fullSelection_MB1or2CR_"+year;
    h_dPhiClusterMET_highClusterSize_fullSelection_MB1or2CR[itr_year] = new TH1D(name,"",70,0,3.5);


    name = "h_nMatchedHitsMB2_fullSelection_SRMB3_"+year;
    h_nMatchedHitsMB2_fullSelection_SRMB3[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB4_fullSelection_SRMB3_"+year;
    h_nMatchedHitsMB4_fullSelection_SRMB3[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2and4_lowClusterSize_fullSelection_SRMB3_"+year;
    h_nMatchedHitsMB2and4_lowClusterSize_fullSelection_SRMB3[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2and4_highClusterSize_fullSelection_SRMB3_"+year;
    h_nMatchedHitsMB2and4_highClusterSize_fullSelection_SRMB3[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2_fullSelection_SRMB4_"+year;
    h_nMatchedHitsMB2_fullSelection_SRMB4[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB3_fullSelection_SRMB4_"+year;
    h_nMatchedHitsMB3_fullSelection_SRMB4[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2and3_lowClusterSize_fullSelection_SRMB4_"+year;
    h_nMatchedHitsMB2and3_lowClusterSize_fullSelection_SRMB4[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2and3_highClusterSize_fullSelection_SRMB4_"+year;
    h_nMatchedHitsMB2and3_highClusterSize_fullSelection_SRMB4[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB3_fullSelection_MB2CR_"+year;
    h_nMatchedHitsMB3_fullSelection_MB2CR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB4_fullSelection_MB2CR_"+year;
    h_nMatchedHitsMB4_fullSelection_MB2CR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB3and4_lowClusterSize_fullSelection_MB2CR_"+year;
    h_nMatchedHitsMB3and4_lowClusterSize_fullSelection_MB2CR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB3and4_highClusterSize_fullSelection_MB2CR_"+year;
    h_nMatchedHitsMB3and4_highClusterSize_fullSelection_MB2CR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2_fullSelectionWithAdjacentMB1Cut_SRMB3_"+year;
    h_nMatchedHitsMB2_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB4_fullSelectionWithAdjacentMB1Cut_SRMB3_"+year;
    h_nMatchedHitsMB4_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3_"+year;
    h_nMatchedHitsMB2and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3_"+year;
    h_nMatchedHitsMB2and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2_fullSelectionWithAdjacentMB1Cut_SRMB4_"+year;
    h_nMatchedHitsMB2_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB3_fullSelectionWithAdjacentMB1Cut_SRMB4_"+year;
    h_nMatchedHitsMB3_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2and3_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4_"+year;
    h_nMatchedHitsMB2and3_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB2and3_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4_"+year;
    h_nMatchedHitsMB2and3_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB3_fullSelectionWithAdjacentMB1Cut_MB2CR_"+year;
    h_nMatchedHitsMB3_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB4_fullSelectionWithAdjacentMB1Cut_MB2CR_"+year;
    h_nMatchedHitsMB4_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB3and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR_"+year;
    h_nMatchedHitsMB3and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year] = new TH1D(name,"",50,0,50);

    name = "h_nMatchedHitsMB3and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR_"+year;
    h_nMatchedHitsMB3and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year] = new TH1D(name,"",50,0,50);


    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);
    
    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SRMB3_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SRMB3[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);
    
    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SRMB4_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SRMB4[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacent0p8MB1Cut_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacent0p8MB1Cut_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);
    
    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SRMB3_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);
    
    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SRMB4_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);
    
    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SRMB3_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SRMB3[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);
    
    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SRMB4_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SRMB4[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);
    
    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCutNoRPC_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCutNoRPC_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCutNoRPC_SR_"+year;
    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCutNoRPC_SR[itr_year] = new TH2D(name,"",50,0,500,70,0,3.5);

    name = "h_nMB1Match_veryLowdPhiClusterMET_lowClusterSize_fullSelection_MB1CR_"+year;
    h_nMB1Match_veryLowdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_lowdPhiClusterMET_lowClusterSize_fullSelection_MB1CR_"+year;
    h_nMB1Match_lowdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_highdPhiClusterMET_lowClusterSize_fullSelection_MB1CR_"+year;
    h_nMB1Match_highdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_veryLowdPhiClusterMET_highClusterSize_fullSelection_MB1CR_"+year;
    h_nMB1Match_veryLowdPhiClusterMET_highClusterSize_fullSelection_MB1CR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_lowdPhiClusterMET_highClusterSize_fullSelection_MB1CR_"+year;
    h_nMB1Match_lowdPhiClusterMET_highClusterSize_fullSelection_MB1CR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_highdPhiClusterMET_highClusterSize_fullSelection_MB1CR_"+year;
    h_nMB1Match_highdPhiClusterMET_highClusterSize_fullSelection_MB1CR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_veryLowdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR_"+year;
    h_nMB1Match_veryLowdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_lowdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR_"+year;
    h_nMB1Match_lowdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_highdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR_"+year;
    h_nMB1Match_highdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_veryLowdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR_"+year;
    h_nMB1Match_veryLowdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_lowdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR_"+year;
    h_nMB1Match_lowdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_nMB1Match_highdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR_"+year;
    h_nMB1Match_highdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH1D(name,"",100,0,200);

    name = "h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1HitsCR_"+year;
    h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH2D(name,"",70,0,3.5,100,0,200);

    name = "h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1Hits30CR_"+year;
    h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1Hits30CR[itr_year] = new TH2D(name,"",70,0,3.5,100,0,200);

    name = "h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1HitsNoMB1ClusterCR_"+year;
    h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[itr_year] = new TH2D(name,"",70,0,3.5,100,0,200);

    name = "h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1HitsCR_"+year;
    h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1HitsCR[itr_year] = new TH2D(name,"",70,0,3.5,100,0,200);

    name = "h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1Hits30CR_"+year;
    h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1Hits30CR[itr_year] = new TH2D(name,"",70,0,3.5,100,0,200);

    name = "h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1HitsNoMB1ClusterCR_"+year;
    h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[itr_year] = new TH2D(name,"",70,0,3.5,100,0,200);

    
    name = "h_nRB1Match_dPhiJetMET_"+year;
    h_nRB1Match_dPhiJetMET[itr_year] = new TH1D(name,"",30,0,30);
    
    name = "h_nRB1Match_MB1Veto_dPhiJetMET_"+year;
    h_nRB1Match_MB1Veto_dPhiJetMET[itr_year] = new TH1D(name,"",30,0,30);
    
    name = "h_nMB1MatchAdjacent_dPhiJetMET_"+year;
    h_nMB1MatchAdjacent_dPhiJetMET[itr_year] = new TH1D(name,"",50,0,50);
    
    name = "h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET_"+year;
    h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_year] = new TH1D(name,"",30,0,30);


    name = "h_nMB1MatchPi2_dPhiClusterMET_"+year;
    h_nMB1MatchPi2_dPhiClusterMET[itr_year] = new TH1D(name,"",3,0,3);

    name = "h_nMB1MatchAdjacent_dPhiClusterMET_"+year;
    h_nMB1MatchAdjacent_dPhiClusterMET[itr_year] = new TH1D(name,"",50,0,50);
    
    name = "h_nMB1MatchAdjacent_MB1Veto_dPhiClusterMET_"+year;
    h_nMB1MatchAdjacent_MB1Veto_dPhiClusterMET[itr_year] = new TH1D(name,"",30,0,30);

    name = "h_nMB1MatchAdjacentPi2_dPhiClusterMET_"+year;
    h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_year] = new TH1D(name,"",50,0,50);
    
    name = "h_nMB1MatchAdjacentPi2_MB1Veto_dPhiClusterMET_"+year;
    h_nMB1MatchAdjacentPi2_MB1Veto_dPhiClusterMET[itr_year] = new TH1D(name,"",30,0,30);
    
    name = "h_nMB1MatchAdjacent0p8_dPhiClusterMET_"+year;
    h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_year] = new TH1D(name,"",50,0,50);
    
    name = "h_nMB1SegMatchAdjacent0p8_dPhiClusterMET_"+year;
    h_nMB1SegMatchAdjacent0p8_dPhiClusterMET[itr_year] = new TH1D(name,"",50,0,50);
    
    name = "h_nMB1MatchAdjacent0p8_MB1Veto_dPhiClusterMET_"+year;
    h_nMB1MatchAdjacent0p8_MB1Veto_dPhiClusterMET[itr_year] = new TH1D(name,"",30,0,30);

    name = "h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET_"+year;
    h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_year] = new TH1D(name,"",50,0,50);
    
    name = "h_nMB1MatchAdjacent0p8Pi2_MB1Veto_dPhiClusterMET_"+year;
    h_nMB1MatchAdjacent0p8Pi2_MB1Veto_dPhiClusterMET[itr_year] = new TH1D(name,"",30,0,30);
    

    name = "h_nMB1MatchJet_invertedJetVeto_"+year;
    h_nMB1MatchJet_invertedJetVeto[itr_year] = new TH1D(name,"",50,0,100);

    name = "h_nMB1MatchJet_MB2_invertedJetVeto_"+year;
    h_nMB1MatchJet_MB2_invertedJetVeto[itr_year] = new TH1D(name,"",50,0,100);

    name = "h_nMB1MatchJet_MB3_invertedJetVeto_"+year;
    h_nMB1MatchJet_MB3_invertedJetVeto[itr_year] = new TH1D(name,"",50,0,100);

    name = "h_nMB1MatchJet_MB4_invertedJetVeto_"+year;
    h_nMB1MatchJet_MB4_invertedJetVeto[itr_year] = new TH1D(name,"",50,0,100);

    name = "h_clusterPhi_passMB1Veto_invertedJetVeto_"+year;
    h_clusterPhi_passMB1Veto_invertedJetVeto[itr_year] = new TH1D(name,"",70,-3.5,3.5);
    
    name = "h_clusterEta_passMB1Veto_invertedJetVeto_"+year;
    h_clusterEta_passMB1Veto_invertedJetVeto[itr_year] = new TH1D(name,"",80,-2,2);

    name = "h_clusterEtaPhi_passMB1Veto_invertedJetVeto_"+year;
    h_clusterEtaPhi_passMB1Veto_invertedJetVeto[itr_year] = new TH2D(name,"",80,-2,2,70,-3.5,3.5);

    name = "h_clusterEtaPhi_MB2_passMB1Veto_invertedJetVeto_"+year;
    h_clusterEtaPhi_MB2_passMB1Veto_invertedJetVeto[itr_year] = new TH2D(name,"",80,-2,2,70,-3.5,3.5);

    name = "h_clusterEtaPhi_MB3_passMB1Veto_invertedJetVeto_"+year;
    h_clusterEtaPhi_MB3_passMB1Veto_invertedJetVeto[itr_year] = new TH2D(name,"",80,-2,2,70,-3.5,3.5);

    name = "h_clusterEtaPhi_MB4_passMB1Veto_invertedJetVeto_"+year;
    h_clusterEtaPhi_MB4_passMB1Veto_invertedJetVeto[itr_year] = new TH2D(name,"",80,-2,2,70,-3.5,3.5);

    name = "h_jetEtaPhi_passMB1Veto_invertedJetVeto_"+year;
    h_jetEtaPhi_passMB1Veto_invertedJetVeto[itr_year] = new TH2D(name,"",80,-2,2,70,-3.5,3.5);

    name = "h_jetEtaPhi_MB2_passMB1Veto_invertedJetVeto_"+year;
    h_jetEtaPhi_MB2_passMB1Veto_invertedJetVeto[itr_year] = new TH2D(name,"",80,-2,2,70,-3.5,3.5);

    name = "h_jetEtaPhi_MB3_passMB1Veto_invertedJetVeto_"+year;
    h_jetEtaPhi_MB3_passMB1Veto_invertedJetVeto[itr_year] = new TH2D(name,"",80,-2,2,70,-3.5,3.5);

    name = "h_jetEtaPhi_MB4_passMB1Veto_invertedJetVeto_"+year;
    h_jetEtaPhi_MB4_passMB1Veto_invertedJetVeto[itr_year] = new TH2D(name,"",80,-2,2,70,-3.5,3.5);

    name = "h_clusterPhi_MB2_passMB1Veto_invertedJetVeto_"+year;
    h_clusterPhi_MB2_passMB1Veto_invertedJetVeto[itr_year] = new TH1D(name,"",70,-3.5,3.5);
    
    name = "h_clusterEta_MB2_passMB1Veto_invertedJetVeto_"+year;
    h_clusterEta_MB2_passMB1Veto_invertedJetVeto[itr_year] = new TH1D(name,"",80,-2,2);

    name = "h_clusterPhi_MB3_passMB1Veto_invertedJetVeto_"+year;
    h_clusterPhi_MB3_passMB1Veto_invertedJetVeto[itr_year] = new TH1D(name,"",70,-3.5,3.5);
    
    name = "h_clusterEta_MB3_passMB1Veto_invertedJetVeto_"+year;
    h_clusterEta_MB3_passMB1Veto_invertedJetVeto[itr_year] = new TH1D(name,"",80,-2,2);

    name = "h_clusterPhi_MB4_passMB1Veto_invertedJetVeto_"+year;
    h_clusterPhi_MB4_passMB1Veto_invertedJetVeto[itr_year] = new TH1D(name,"",70,-3.5,3.5);
    
    name = "h_clusterEta_MB4_passMB1Veto_invertedJetVeto_"+year;
    h_clusterEta_MB4_passMB1Veto_invertedJetVeto[itr_year] = new TH1D(name,"",80,-2,2);

    name = "h_matchedJetPt_invertedJetVeto_"+year;
    h_matchedJetPt_invertedJetVeto[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_matchedJetPt_nMB1Match_invertedJetVeto_"+year;
    h_matchedJetPt_nMB1Match_invertedJetVeto[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_invertedJetVetoNoJetMET_"+year;
    h_matchedJetPt_nMB1Match_invertedJetVetoNoJetMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB2_invertedJetVetoNoJetMET_"+year;
    h_matchedJetPt_nMB1Match_MB2_invertedJetVetoNoJetMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB3_invertedJetVetoNoJetMET_"+year;
    h_matchedJetPt_nMB1Match_MB3_invertedJetVetoNoJetMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB4_invertedJetVetoNoJetMET_"+year;
    h_matchedJetPt_nMB1Match_MB4_invertedJetVetoNoJetMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_invertedJetVetoNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_invertedJetVetoNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_invertedJetVetoAllMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_invertedJetVetoAllMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB2_invertedJetVetoNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB2_invertedJetVetoNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB2_invertedJetVetoAllMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB2_invertedJetVetoAllMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);
    
    name = "h_matchedJetPt_nMB1Match_MB2_invertedJetVetoInvertedMB1LooseMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB2_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB3_invertedJetVetoInvertedMB1LooseMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB3_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB4_invertedJetVetoInvertedMB1LooseMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB4_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonLowClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonLowClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonLowSizeNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonLowSizeNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_matchedJetPt_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_matchedJetPt_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMETNoSize_"+year;
    h_matchedJetPt_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMETNoSize[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_MET_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_MET_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonLowClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonLowClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonLowSizeNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonLowSizeNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_matchedJetPt_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_matchedJetPt_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMETNoSize_"+year;
    h_matchedJetPt_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMETNoSize[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_MET_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_MET_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonLowClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonLowClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonLowSizeNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonLowSizeNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_matchedJetPt_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_matchedJetPt_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMETNoSize_"+year;
    h_matchedJetPt_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMETNoSize[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_MET_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_MET_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,500);

    name = "h_matchedJetPt_nMB1Match_MB3_invertedJetVetoNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB3_invertedJetVetoNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB3_invertedJetVetoAllMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB3_invertedJetVetoAllMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB4_invertedJetVetoNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB4_invertedJetVetoNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB4_invertedJetVetoAllMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB4_invertedJetVetoAllMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonNoClusterMET_"+year;
    h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonNoClusterMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_invertedJetVetoNoClusterMETNoJetMET_"+year;
    h_matchedJetPt_nMB1Match_invertedJetVetoNoClusterMETNoJetMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_invertedJetVetoAllMuonNoClusterMETNoJetMET_"+year;
    h_matchedJetPt_nMB1Match_invertedJetVetoAllMuonNoClusterMETNoJetMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_invertedJetVetoLooseMuonNoClusterMETNoJetMET_"+year;
    h_matchedJetPt_nMB1Match_invertedJetVetoLooseMuonNoClusterMETNoJetMET[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonNoClusterMETCluster80_"+year;
    h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonNoClusterMETCluster80[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonNoClusterMETCluster80_"+year;
    h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonNoClusterMETCluster80[itr_year] = new TH2D(name,"",100,0,500,100,0,100);

    name = "h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonNoClusterMETCluster80_"+year;
    h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonNoClusterMETCluster80[itr_year] = new TH2D(name,"",100,0,500,100,0,100);


    name = "h_jetNeutralHadronicEnergyFraction_passJetMET_"+year;
    h_jetNeutralHadronicEnergyFraction_passJetMET[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetNeutralHadronicEnergyFraction_passStationsWheels_"+year;
    h_jetNeutralHadronicEnergyFraction_passStationsWheels[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetNeutralHadronicEnergyFraction_SR_"+year;
    h_jetNeutralHadronicEnergyFraction_SR[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedHadronicEnergyFraction_passJetMET_"+year;
    h_jetChargedHadronicEnergyFraction_passJetMET[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedHadronicEnergyFraction_passStationsWheels_"+year;
    h_jetChargedHadronicEnergyFraction_passStationsWheels[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedHadronicEnergyFraction_SR_"+year;
    h_jetChargedHadronicEnergyFraction_SR[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetNeutralEMEnergyFraction_passJetMET_"+year;
    h_jetNeutralEMEnergyFraction_passJetMET[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetNeutralEMEnergyFraction_passStationsWheels_"+year;
    h_jetNeutralEMEnergyFraction_passStationsWheels[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetNeutralEMEnergyFraction_SR_"+year;
    h_jetNeutralEMEnergyFraction_SR[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedEMEnergyFraction_passJetMET_"+year;
    h_jetChargedEMEnergyFraction_passJetMET[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedEMEnergyFraction_passStationsWheels_"+year;
    h_jetChargedEMEnergyFraction_passStationsWheels[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedEMEnergyFraction_SR_"+year;
    h_jetChargedEMEnergyFraction_SR[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_jetChargedEMEnergyFraction_MB2CR_"+year;
    h_jetChargedEMEnergyFraction_MB2CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetChargedEMEnergyFraction_MB2withMB1CR_"+year;
    h_jetChargedEMEnergyFraction_MB2withMB1CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetChargedEMEnergyFraction_MB1CR_"+year;
    h_jetChargedEMEnergyFraction_MB1CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetChargedEMEnergyFraction_MB1HitsCR_"+year;
    h_jetChargedEMEnergyFraction_MB1HitsCR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetNeutralEMEnergyFraction_MB2CR_"+year;
    h_jetNeutralEMEnergyFraction_MB2CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetNeutralEMEnergyFraction_MB2withMB1CR_"+year;
    h_jetNeutralEMEnergyFraction_MB2withMB1CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetNeutralEMEnergyFraction_MB1CR_"+year;
    h_jetNeutralEMEnergyFraction_MB1CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetNeutralEMEnergyFraction_MB1HitsCR_"+year;
    h_jetNeutralEMEnergyFraction_MB1HitsCR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetChargedHadronicEnergyFraction_MB2CR_"+year;
    h_jetChargedHadronicEnergyFraction_MB2CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetChargedHadronicEnergyFraction_MB2withMB1CR_"+year;
    h_jetChargedHadronicEnergyFraction_MB2withMB1CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetChargedHadronicEnergyFraction_MB1CR_"+year;
    h_jetChargedHadronicEnergyFraction_MB1CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetChargedHadronicEnergyFraction_MB1HitsCR_"+year;
    h_jetChargedHadronicEnergyFraction_MB1HitsCR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetNeutralHadronicEnergyFraction_MB2CR_"+year;
    h_jetNeutralHadronicEnergyFraction_MB2CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetNeutralHadronicEnergyFraction_MB2withMB1CR_"+year;
    h_jetNeutralHadronicEnergyFraction_MB2withMB1CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetNeutralHadronicEnergyFraction_MB1CR_"+year;
    h_jetNeutralHadronicEnergyFraction_MB1CR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_jetNeutralHadronicEnergyFraction_MB1HitsCR_"+year;
    h_jetNeutralHadronicEnergyFraction_MB1HitsCR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_leadingJetChargedEMEnergyFraction_passStationsWheels_"+year;
    h_leadingJetChargedEMEnergyFraction_passStationsWheels[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetChargedEMEnergyFraction_SR_"+year;
    h_leadingJetChargedEMEnergyFraction_SR[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetNeutralHadronicEnergyFraction_passStationsWheels_"+year;
    h_leadingJetNeutralHadronicEnergyFraction_passStationsWheels[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_leadingJetNeutralHadronicEnergyFraction_SR_"+year;
    h_leadingJetNeutralHadronicEnergyFraction_SR[itr_year] = new TH1D(name,"",50,0,1);
    
    name = "h_leadingJetNeutralEMEnergyFraction_passStationsWheels_"+year;
    h_leadingJetNeutralEMEnergyFraction_passStationsWheels[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetNeutralEMEnergyFraction_SR_"+year;
    h_leadingJetNeutralEMEnergyFraction_SR[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetChargedHadronicEnergyFraction_passStationsWheels_"+year;
    h_leadingJetChargedHadronicEnergyFraction_passStationsWheels[itr_year] = new TH1D(name,"",50,0,1);

    name = "h_leadingJetChargedHadronicEnergyFraction_SR_"+year;
    h_leadingJetChargedHadronicEnergyFraction_SR[itr_year] = new TH1D(name,"",50,0,1);

    
    name = "h_dtRechitClusterSize_SRhighNHF_"+year;
    h_dtRechitClusterSize_SRhighNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_SRhighNHF_"+year;
    h_MET_SRhighNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_SRhighNHF_"+year;
    h_dPhiClusterMET_SRhighNHF[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_MB2CRhighNHF_"+year;
    h_dtRechitClusterSize_MB2CRhighNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_MB2CRhighNHF_"+year;
    h_MET_MB2CRhighNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_MB2CRhighNHF_"+year;
    h_dPhiClusterMET_MB2CRhighNHF[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_MB1CRhighNHF_"+year;
    h_dtRechitClusterSize_MB1CRhighNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_MB1CRhighNHF_"+year;
    h_MET_MB1CRhighNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_MB1CRhighNHF_"+year;
    h_dPhiClusterMET_MB1CRhighNHF[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_MB2withMB1CRhighNHF_"+year;
    h_dtRechitClusterSize_MB2withMB1CRhighNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_MB2withMB1CRhighNHF_"+year;
    h_MET_MB2withMB1CRhighNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_MB2withMB1CRhighNHF_"+year;
    h_dPhiClusterMET_MB2withMB1CRhighNHF[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_MB1HitsCRhighNHF_"+year;
    h_dtRechitClusterSize_MB1HitsCRhighNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_MB1HitsCRhighNHF_"+year;
    h_MET_MB1HitsCRhighNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_MB1HitsCRhighNHF_"+year;
    h_dPhiClusterMET_MB1HitsCRhighNHF[itr_year] = new TH1D(name,"",70,0,3.5);
    for (auto selsIt = selsStrings.begin(); selsIt != selsStrings.end(); selsIt++){
	    TString tname = *selsIt+"_"+year;
	    runNumHists[tname] = new TH1D("h_runNum_"+tname,";;",1000,271000,326000);
	    etaClusterHists[tname] = new TH1D("h_eta_"+tname,";;",100,-2,2);
	    phiClusterHists[tname] = new TH1D("h_phi_"+tname,";;",100,-3.141,3.141);
	    etaPhiClusterHists[tname] = new TH2D("h_etaPhiCluster_"+tname,";;",100,-3.141,3.141,100,-2,2);
    }

    name = "h_MET_lowNHF_"+year;
    h_MET_lowNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_MET_oneCluster_lowNHF_"+year;
    h_MET_oneCluster_lowNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_MET_oneVetoCluster_lowNHF_"+year;
    h_MET_oneVetoCluster_lowNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_MET_highNHF_"+year;
    h_MET_highNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_MET_oneCluster_highNHF_"+year;
    h_MET_oneCluster_highNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_MET_oneVetoCluster_highNHF_"+year;
    h_MET_oneVetoCluster_highNHF[itr_year] = new TH1D(name,"",100,0,1000);

    
    name = "h_leadJetPt_highNHF_"+year;
    h_leadJetPt_highNHF[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_lowNHF_"+year;
    h_leadJetPt_lowNHF[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_oneCluster_highNHF_"+year;
    h_leadJetPt_oneCluster_highNHF[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_oneCluster_lowNHF_"+year;
    h_leadJetPt_oneCluster_lowNHF[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_oneVetoCluster_highNHF_"+year;
    h_leadJetPt_oneVetoCluster_highNHF[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_oneVetoCluster_lowNHF_"+year;
    h_leadJetPt_oneVetoCluster_lowNHF[itr_year] = new TH1D(name,"",200,0,2000);
    
    name = "h_leadJetPt_MB1CR_"+year;
    h_leadJetPt_MB1CR[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_MB2withMB1CR_"+year;
    h_leadJetPt_MB2withMB1CR[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_MB1HitsCR_"+year;
    h_leadJetPt_MB1HitsCR[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_MB2CR_"+year;
    h_leadJetPt_MB2CR[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPt_SR_"+year;
    h_leadJetPt_SR[itr_year] = new TH1D(name,"",200,0,2000);

    name = "h_leadJetPtMET_highNHF_"+year;
    h_leadJetPtMET_highNHF[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_lowNHF_"+year;
    h_leadJetPtMET_lowNHF[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_oneCluster_highNHF_"+year;
    h_leadJetPtMET_oneCluster_highNHF[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_oneCluster_lowNHF_"+year;
    h_leadJetPtMET_oneCluster_lowNHF[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_oneVetoCluster_highNHF_"+year;
    h_leadJetPtMET_oneVetoCluster_highNHF[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_oneVetoCluster_lowNHF_"+year;
    h_leadJetPtMET_oneVetoCluster_lowNHF[itr_year] = new TH1D(name,"",500,0,10);
    
    name = "h_leadJetPtMET_MB1CR_"+year;
    h_leadJetPtMET_MB1CR[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_MB2withMB1CR_"+year;
    h_leadJetPtMET_MB2withMB1CR[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_MB1HitsCR_"+year;
    h_leadJetPtMET_MB1HitsCR[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_MB2CR_"+year;
    h_leadJetPtMET_MB2CR[itr_year] = new TH1D(name,"",500,0,10);

    name = "h_leadJetPtMET_SR_"+year;
    h_leadJetPtMET_SR[itr_year] = new TH1D(name,"",500,0,10);


    name = "h_dtRechitClusterSize_SRlowNHF_"+year;
    h_dtRechitClusterSize_SRlowNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_SRlowNHF_"+year;
    h_MET_SRlowNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_SRlowNHF_"+year;
    h_dPhiClusterMET_SRlowNHF[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_MB2CRlowNHF_"+year;
    h_dtRechitClusterSize_MB2CRlowNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_MB2CRlowNHF_"+year;
    h_MET_MB2CRlowNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_MB2CRlowNHF_"+year;
    h_dPhiClusterMET_MB2CRlowNHF[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_MB1CRlowNHF_"+year;
    h_dtRechitClusterSize_MB1CRlowNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_MB1CRlowNHF_"+year;
    h_MET_MB1CRlowNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_MB1CRlowNHF_"+year;
    h_dPhiClusterMET_MB1CRlowNHF[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_MB2withMB1CRlowNHF_"+year;
    h_dtRechitClusterSize_MB2withMB1CRlowNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_MB2withMB1CRlowNHF_"+year;
    h_MET_MB2withMB1CRlowNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_MB2withMB1CRlowNHF_"+year;
    h_dPhiClusterMET_MB2withMB1CRlowNHF[itr_year] = new TH1D(name,"",70,0,3.5);

    name = "h_dtRechitClusterSize_MB1HitsCRlowNHF_"+year;
    h_dtRechitClusterSize_MB1HitsCRlowNHF[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_MET_MB1HitsCRlowNHF_"+year;
    h_MET_MB1HitsCRlowNHF[itr_year] = new TH1D(name,"",100,0,1000);

    name = "h_dPhiClusterMET_MB1HitsCRlowNHF_"+year;
    h_dPhiClusterMET_MB1HitsCRlowNHF[itr_year] = new TH1D(name,"",70,0,3.5);

   
    name = "h_rpcBxMedian_MB1CR_"+year;
    h_rpcBxMedian_MB1CR[itr_year] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_rpcBxMedian_MB2CR_"+year;
    h_rpcBxMedian_MB2CR[itr_year] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_rpcBxMedian_MB2withMB1CR_"+year;
    h_rpcBxMedian_MB2withMB1CR[itr_year] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_rpcBxMedian_SR_"+year;
    h_rpcBxMedian_SR[itr_year] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_rpcBxMedian_MB1HitsCR_"+year;
    h_rpcBxMedian_MB1HitsCR[itr_year] = new TH1D(name,"",10,-4.5,5.5);

    name = "h_invertedLeadJetVeto_failedMETFilters_SR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_SR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_MB2CR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_MB2CR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR[itr_year] = new TH1D(name,"",6,0,6);

    name = "h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR_"+year;
    h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR[itr_year] = new TH1D(name,"",6,0,6);

    
    name = "h_nCAClusters_"+year;
    h_nCAClusters[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_sizeCACluster_"+year;
    h_sizeCACluster[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_radiusCACluster_"+year;
    h_radiusCACluster[itr_year] = new TH1D(name,"",50,300,800);

    name = "h_zCACluster_"+year;
    h_zCACluster[itr_year] = new TH1D(name,"",120,-600,600);

    name = "h_phiCACluster_"+year;
    h_phiCACluster[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_nDtRechitClusters_"+year;
    h_nDtRechitClusters[itr_year] = new TH1D(name,"",5,0,5);

    name = "h_dtRechitClusterSize_"+year;
    h_dtRechitClusterSize[itr_year] = new TH1D(name,"",50,0,500);

    name = "h_dtRechitClusterR_"+year;
    h_dtRechitClusterR[itr_year] = new TH1D(name,"",50,300,800);

    name = "h_dtRechitClusterZ_"+year;
    h_dtRechitClusterZ[itr_year] = new TH1D(name,"",120,-600,600);

    name = "h_dtRechitClusterPhi_"+year;
    h_dtRechitClusterPhi[itr_year] = new TH1D(name,"",70,-3.5,3.5);


    name = "h_matchedJetNHF_"+year;
    h_matchedJetNHF[itr_year] = new TH1D(name,"",100,0,1);

    name = "h_matchedJetNEMF_"+year;
    h_matchedJetNEMF[itr_year] = new TH1D(name,"",100,0,1);

    name = "h_matchedJetCHF_"+year;
    h_matchedJetCHF[itr_year] = new TH1D(name,"",100,0,1);

    name = "h_matchedJetCEMF_"+year;
    h_matchedJetCEMF[itr_year] = new TH1D(name,"",100,0,1);

    name = "h_matchedJetNConst_"+year;
    h_matchedJetNConst[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_matchedJetChargedMult_"+year;
    h_matchedJetChargedMult[itr_year] = new TH1D(name,"",20,0,20);
  


    name = "h_clusterSize_looseMuonVetoResult_ABCD_"+year;
    h_clusterSize_looseMuonVetoResult_ABCD[itr_year] = new TH2D(name,"",100,0,500,2,0,2);

    name = "h_clusterSize_recoMuonVetoResult_ABCD_"+year;
    h_clusterSize_recoMuonVetoResult_ABCD[itr_year] = new TH2D(name,"",100,0,500,2,0,2);

    name = "h_clusterSize_looseMuonVeto_D_"+year;
    h_clusterSize_looseMuonVeto_D[itr_year] = new TH1D(name,"",100,100,500);

    name = "h_clusterSize_looseMuonVeto_C_"+year;
    h_clusterSize_looseMuonVeto_C[itr_year] = new TH1D(name,"",25,51,101);

    name = "h_clusterSizeMuonMB1_looseMuonVeto_C_"+year;
    h_clusterSizeMuonMB1_looseMuonVeto_C[itr_year] = new TH1D(name,"",50,51,101);

    name = "h_clusterSize_looseMuonVeto_B_"+year;
    h_clusterSize_looseMuonVeto_B[itr_year] = new TH1D(name,"",100,100,500);

    name = "h_clusterSize_looseMuonVeto_A_"+year;
    h_clusterSize_looseMuonVeto_A[itr_year] = new TH1D(name,"",25,51,101);

    name = "h_matchedMB1Hits_looseMuonVeto_D_"+year;
    h_matchedMB1Hits_looseMuonVeto_D[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_looseMuonVeto_C_"+year;
    h_matchedMB1Hits_looseMuonVeto_C[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_looseMuonVeto_B_"+year;
    h_matchedMB1Hits_looseMuonVeto_B[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_looseMuonVeto_A_"+year;
    h_matchedMB1Hits_looseMuonVeto_A[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_invertedLooseMuonVeto_D_"+year;
    h_matchedMB1Hits_invertedLooseMuonVeto_D[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_invertedLooseMuonVeto_C_"+year;
    h_matchedMB1Hits_invertedLooseMuonVeto_C[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_invertedLooseMuonVeto_B_"+year;
    h_matchedMB1Hits_invertedLooseMuonVeto_B[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_invertedLooseMuonVeto_A_"+year;
    h_matchedMB1Hits_invertedLooseMuonVeto_A[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_invertedRecoMuonVeto_D_"+year;
    h_matchedMB1Hits_invertedRecoMuonVeto_D[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_invertedRecoMuonVeto_C_"+year;
    h_matchedMB1Hits_invertedRecoMuonVeto_C[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_invertedRecoMuonVeto_B_"+year;
    h_matchedMB1Hits_invertedRecoMuonVeto_B[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_invertedRecoMuonVeto_A_"+year;
    h_matchedMB1Hits_invertedRecoMuonVeto_A[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Segments_looseMuonVeto_D_"+year;
    h_matchedMB1Segments_looseMuonVeto_D[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_matchedMB1Segments_looseMuonVeto_C_"+year;
    h_matchedMB1Segments_looseMuonVeto_C[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_matchedMB1Segments_looseMuonVeto_B_"+year;
    h_matchedMB1Segments_looseMuonVeto_B[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_matchedMB1Segments_looseMuonVeto_A_"+year;
    h_matchedMB1Segments_looseMuonVeto_A[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_minMB1SegmentDR_looseMuonVeto_D_"+year;
    h_minMB1SegmentDR_looseMuonVeto_D[itr_year] = new TH1D(name,"",100,0,5);

    name = "h_minMB1SegmentDR_looseMuonVeto_C_"+year;
    h_minMB1SegmentDR_looseMuonVeto_C[itr_year] = new TH1D(name,"",100,0,5);

    name = "h_minMB1SegmentDR_looseMuonVeto_B_"+year;
    h_minMB1SegmentDR_looseMuonVeto_B[itr_year] = new TH1D(name,"",100,0,5);

    name = "h_minMB1SegmentDR_looseMuonVeto_A_"+year;
    h_minMB1SegmentDR_looseMuonVeto_A[itr_year] = new TH1D(name,"",100,0,5);

    name = "h_matchedMB1HitsRatio_looseMuonVeto_D_"+year;
    h_matchedMB1HitsRatio_looseMuonVeto_D[itr_year] = new TH1D(name,"",20,0,1);
    
    name = "h_matchedMB1HitsRatio_looseMuonVeto_C_"+year;
    h_matchedMB1HitsRatio_looseMuonVeto_C[itr_year] = new TH1D(name,"",20,0,1);
    
    name = "h_matchedMB1HitsRatio_looseMuonVeto_B_"+year;
    h_matchedMB1HitsRatio_looseMuonVeto_B[itr_year] = new TH1D(name,"",20,0,1);

    name = "h_matchedMB1HitsRatio_looseMuonVeto_A_"+year;
    h_matchedMB1HitsRatio_looseMuonVeto_A[itr_year] = new TH1D(name,"",20,0,1);

    name = "h_matchedMB1Hits_looseMuonVetoDoubleInverted_D_"+year;
    h_matchedMB1Hits_looseMuonVetoDoubleInverted_D[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_looseMuonVetoDoubleInverted_C_"+year;
    h_matchedMB1Hits_looseMuonVetoDoubleInverted_C[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_looseMuonVetoDoubleInverted_B_"+year;
    h_matchedMB1Hits_looseMuonVetoDoubleInverted_B[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Hits_looseMuonVetoDoubleInverted_A_"+year;
    h_matchedMB1Hits_looseMuonVetoDoubleInverted_A[itr_year] = new TH1D(name,"",20,0,40);

    name = "h_matchedMB1Segments_looseMuonVetoDoubleInverted_D_"+year;
    h_matchedMB1Segments_looseMuonVetoDoubleInverted_D[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_matchedMB1Segments_looseMuonVetoDoubleInverted_C_"+year;
    h_matchedMB1Segments_looseMuonVetoDoubleInverted_C[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_matchedMB1Segments_looseMuonVetoDoubleInverted_B_"+year;
    h_matchedMB1Segments_looseMuonVetoDoubleInverted_B[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_matchedMB1Segments_looseMuonVetoDoubleInverted_A_"+year;
    h_matchedMB1Segments_looseMuonVetoDoubleInverted_A[itr_year] = new TH1D(name,"",20,0,20);

    name = "h_minMB1SegmentDR_looseMuonVetoDoubleInverted_D_"+year;
    h_minMB1SegmentDR_looseMuonVetoDoubleInverted_D[itr_year] = new TH1D(name,"",100,0,5);

    name = "h_minMB1SegmentDR_looseMuonVetoDoubleInverted_C_"+year;
    h_minMB1SegmentDR_looseMuonVetoDoubleInverted_C[itr_year] = new TH1D(name,"",100,0,5);

    name = "h_minMB1SegmentDR_looseMuonVetoDoubleInverted_B_"+year;
    h_minMB1SegmentDR_looseMuonVetoDoubleInverted_B[itr_year] = new TH1D(name,"",100,0,5);

    name = "h_minMB1SegmentDR_looseMuonVetoDoubleInverted_A_"+year;
    h_minMB1SegmentDR_looseMuonVetoDoubleInverted_A[itr_year] = new TH1D(name,"",100,0,5);

    name = "h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_D_"+year;
    h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_D[itr_year] = new TH1D(name,"",20,0,1);
    
    name = "h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_C_"+year;
    h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_C[itr_year] = new TH1D(name,"",20,0,1);
    
    name = "h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_B_"+year;
    h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_B[itr_year] = new TH1D(name,"",20,0,1);

    name = "h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_A_"+year;
    h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_A[itr_year] = new TH1D(name,"",20,0,1);

    
    name = "h_dtRechitClusterPhi_invertedMB1AC_10GeVJet_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_10GeVJet[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1AC_jetVeto20_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_jetVeto20[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1AC_jetVeto10_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_jetVeto10[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1AC_MB3_10GeVJet_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_MB3_10GeVJet[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1AC_MB3_jetVeto20_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_MB3_jetVeto20[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1AC_MB3_jetVeto10_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_MB3_jetVeto10[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1AC_MB4_10GeVJet_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_MB4_10GeVJet[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1AC_MB4_jetVeto20_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_MB4_jetVeto20[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1AC_MB4_jetVeto10_"+year;
    h_dtRechitClusterPhi_invertedMB1AC_MB4_jetVeto10[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1A_10GeVJet_"+year;
    h_dtRechitClusterPhi_invertedMB1A_10GeVJet[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1A_jetVeto20_"+year;
    h_dtRechitClusterPhi_invertedMB1A_jetVeto20[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1A_jetVeto10_"+year;
    h_dtRechitClusterPhi_invertedMB1A_jetVeto10[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1C_10GeVJet_"+year;
    h_dtRechitClusterPhi_invertedMB1C_10GeVJet[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1C_jetVeto20_"+year;
    h_dtRechitClusterPhi_invertedMB1C_jetVeto20[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_invertedMB1C_jetVeto10_"+year;
    h_dtRechitClusterPhi_invertedMB1C_jetVeto10[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_AC_10GeVJet_"+year;
    h_dtRechitClusterPhi_AC_10GeVJet[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_AC_jetVeto20_"+year;
    h_dtRechitClusterPhi_AC_jetVeto20[itr_year] = new TH1D(name,"",70,-3.5,3.5);

    name = "h_dtRechitClusterPhi_AC_jetVeto10_"+year;
    h_dtRechitClusterPhi_AC_jetVeto10[itr_year] = new TH1D(name,"",70,-3.5,3.5);


    name = "h_jetOverlap10GeV_"+year;
    h_jetOverlap10GeV[itr_year] = new TH1D(name,"",2,0,2);

    name = "h_jetOverlap20GeV_"+year;
    h_jetOverlap20GeV[itr_year] = new TH1D(name,"",2,0,2);

    name = "h_muonOverlap10GeVLoose_"+year;
    h_muonOverlap10GeVLoose[itr_year] = new TH1D(name,"",2,0,2);

    name = "h_muonOverlap10GeVTight_"+year;
    h_muonOverlap10GeVTight[itr_year] = new TH1D(name,"",2,0,2);


    name = "h_efficiency_"+year;
    h_efficiency[itr_year] = new TH1D(name,"",20,0,20);

    TFile *_file;
    if(year == "All"){
      _file = TFile::Open(dir+year+"/v4/v4/normalized/Run2_displacedJetMuonNtupler_V1p15_Data2016_Data2017_Data2018-HighMET_goodLumi.root");
    }
    else{
      //_file = TFile::Open(dir+years[itr_year]+"/v1/v3/normalized/Run2_displacedJetMuonNtupler_V1p17_Data"+years[itr_year]+"_"+runNames[itr_year]+"-HighMET-"+dates[itr_year]+"_goodLumi.root");
      _file = TFile::Open(dir+"Run2_displacedJetMuonNtupler_V1p17_Data"+year+"_"+runNames[itr_year]+"-HighMET-"+dates[itr_year]+"_goodLumi.root");
    }

    TTreeReader treeReader("MuonSystem",_file);
    
    TTreeReaderValue<bool> Flag2_all(treeReader,"Flag2_all");
    TTreeReaderValue<bool> Flag2_HBHENoiseFilter(treeReader,"Flag2_HBHENoiseFilter");
    TTreeReaderValue<bool> Flag2_HBHEIsoNoiseFilter(treeReader,"Flag2_HBHEIsoNoiseFilter");
    TTreeReaderValue<bool> Flag2_BadPFMuonFilter(treeReader,"Flag2_BadPFMuonFilter");
    TTreeReaderValue<bool> Flag2_globalSuperTightHalo2016Filter(treeReader,"Flag2_globalSuperTightHalo2016Filter");
    TTreeReaderValue<bool> Flag2_EcalDeadCellTriggerPrimitiveFilter(treeReader,"Flag2_EcalDeadCellTriggerPrimitiveFilter");
    TTreeReaderValue<bool> Flag2_ecalBadCalibFilter(treeReader,"Flag2_ecalBadCalibFilter");
    TTreeReaderArray<bool> HLTDecision(treeReader,"HLTDecision");

    TTreeReaderValue<unsigned int> runNum(treeReader,"runNum");
    TTreeReaderValue<unsigned int> lumiSec(treeReader,"lumiSec");
    TTreeReaderValue<unsigned int> eventNum(treeReader,"evtNum");
    TTreeReaderValue<float> MET(treeReader,"met");
    TTreeReaderValue<float> METphi(treeReader,"metPhi");

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

    TTreeReaderValue<int> nDtRechits(treeReader,"nDtRechits");
    TTreeReaderArray<float> dtRechitX(treeReader,"dtRechitX");
    TTreeReaderArray<float> dtRechitY(treeReader,"dtRechitY");
    TTreeReaderArray<float> dtRechitZ(treeReader,"dtRechitZ");
    TTreeReaderArray<float> dtRechitEta(treeReader,"dtRechitEta");
    TTreeReaderArray<float> dtRechitPhi(treeReader,"dtRechitPhi");
    TTreeReaderArray<int> dtRechitStation(treeReader,"dtRechitStation");
    TTreeReaderArray<int> dtRechitWheel(treeReader,"dtRechitWheel");

    TTreeReaderValue<int> nDtSeg(treeReader,"nDtSeg");
    TTreeReaderArray<float> dtSegX(treeReader,"dtSegX");
    TTreeReaderArray<float> dtSegY(treeReader,"dtSegY");
    TTreeReaderArray<float> dtSegZ(treeReader,"dtSegZ");
    TTreeReaderArray<float> dtSegEta(treeReader,"dtSegEta");
    TTreeReaderArray<float> dtSegPhi(treeReader,"dtSegPhi");
    TTreeReaderArray<int> dtSegStation(treeReader,"dtSegStation");
    TTreeReaderArray<int> dtSegWheel(treeReader,"dtSegWheel");
    TTreeReaderArray<float> dtSegTime(treeReader,"dtSegTime");
    TTreeReaderArray<float> dtSegTimeError(treeReader,"dtSegTimeError");

    TTreeReaderValue<int> nJets(treeReader,"nJets");
    TTreeReaderArray<float> jetPt(treeReader,"jetPt");
    TTreeReaderArray<float> jetEta(treeReader,"jetEta");
    TTreeReaderArray<float> jetPhi(treeReader,"jetPhi");
    TTreeReaderArray<float> jetRechitT(treeReader,"jetRechitT");
    TTreeReaderArray<float> jetRechitT_rms(treeReader,"jetRechitT_rms");
    TTreeReaderArray<float> jetElectronEnergyFraction(treeReader,"jetElectronEnergyFraction");
    TTreeReaderArray<float> jetPhotonEnergyFraction(treeReader,"jetPhotonEnergyFraction");
    TTreeReaderArray<float> jetNeutralHadronEnergyFraction(treeReader,"jetNeutralHadronEnergyFraction");
    TTreeReaderArray<float> jetChargedHadronEnergyFraction(treeReader,"jetChargedHadronEnergyFraction");
    TTreeReaderArray<float> jetMuonEnergyFraction(treeReader,"jetMuonEnergyFraction");
    TTreeReaderArray<bool>  jetTightPassId(treeReader,"jetTightPassId");

    TTreeReaderValue<int> nRPCRechits(treeReader,"nRpc");
    TTreeReaderArray<float> RPCRechitX(treeReader,"rpcX");
    TTreeReaderArray<float> RPCRechitY(treeReader,"rpcY");
    TTreeReaderArray<float> RPCRechitZ(treeReader,"rpcZ");
    TTreeReaderArray<float> RPCRechitPhi(treeReader,"rpcPhi");
    TTreeReaderArray<float> RPCRechitEta(treeReader,"rpcEta");
    TTreeReaderArray<int> RPCRechitBx(treeReader,"rpcBx");

    TTreeReaderValue<int> nLeptons(treeReader,"nLeptons");
    TTreeReaderValue<int> nMuons(treeReader,"nMuons");
    TTreeReaderArray<int> lepPdgId(treeReader,"lepPdgId");
    TTreeReaderArray<bool> lepPassId(treeReader,"lepPassId");
    TTreeReaderArray<float> lepPhi(treeReader,"lepPhi");
    TTreeReaderArray<float> lepEta(treeReader,"lepEta");
    TTreeReaderArray<float> lepPt(treeReader,"lepPt");

    _ofile->cd();
    evtNum = 0;
    event100_rpc = 0;
    event150_rpc = 0;
    event100_clusterMET = 0;
    event150_clusterMET = 0;
    
    cout << "Data" << year << endl;
    while(treeReader.Next()){
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
      passJetMET_rpcCR = false;
      passMaxStation_rpcCR = false;
      passLepton_rpcCR = false;
      pass50Hits_rpcCR = false;
      pass25Hits_rpcCR = false;
      
      passMET = false;
      passOneJet = false;
      passJetMET = false;
      passNHFJet50 = true;
      passNHFJetLead = false;
      passStations25 = false;
      passWheels25 = false;
      passJetVeto = false;
      passMuonVeto = false;
      passMuonVetoLoose = false;
      passRpcMatch = false;
      passRpcMatchMB1CRwithNHFLead = false;
      passRpcMatchMB1HitsCR = false;
      passRpcMatchMB1HitsCRwithAdjacent = false;
      passRpcMatchMB1HitsCRwithNHF50 = false;
      passRpcMatchMB1HitsCRinvertedNHF50 = false;
      passRpcMatchMB1HitsCRwithNHFLead = false;
      passRpcMatchMB1HitsCRinvertedNHFLead = false;
      passRpcMatchMB1Hits30CR = false;
      passRpcMatchMB1Hits30MaxMB2CR = false;
      passRpcMatchMB1Hits30MaxMB3CR = false;
      passRpcMatchMB1Hits30MaxMB4CR = false;
      passRpcMatchMB1HitsNoMB1ClusterCR = false;
      passRpcMatchMB1or2CR = false;
      passRpcMatchMB2CR = false;
      passRpcMatchMB2CRwithAdjacent = false;
      passRpcMatchMB2CRwithAdjacent0p8 = false;
      passRpcMatchMB2CRwithOther = false;
      passRpcMatchMB2CRwithNHF50 = false;
      passRpcMatchMB2CRinvertedNHF50 = false;
      passRpcMatchMB2CRwithNHFLead = false;
      noRpcMatchMB2CRwithNHFLead = false;
      passRpcMatchMB2CRinvertedNHFLead = false;
      passRpcMatchMB2withMB1CR = false;
      passRpcMatchMB2withMB1CRwithAdjacent = false;
      passRpcMatchMB2withMB1CRwithNHF50 = false;
      passRpcMatchMB2withMB1CRinvertedNHF50 = false;
      passRpcMatchMB2withMB1CRwithNHFLead = false;
      passRpcMatchMB2withMB1CRinvertedNHFLead = false;
      passMB2CRwithInvertedNHFLeadNoRPC = false;
      passMB1Veto = false;
      passMaxStation = false;
      passClusterMET = false;
      passNoVetoCluster = false;
      passClusterSize = false;
      passMB1CR = false;
      clusterSizeMB1CR = 0;
      dPhiClusterMETMB1CR = 999.;      
      clusterSizeMB1CRwithNHFLead = 0;
      dPhiClusterMETMB1CRwithNHFLead = 999.;      
      clusterSizeMB1HitsCR = 0;
      dPhiClusterMETMB1HitsCR = 999.;      
      clusterSizeMB1HitsCRwithAdjacent = 0;
      dPhiClusterMETMB1HitsCRwithAdjacent = 999.;      
      clusterSizeMB1HitsCRwithNHF50 = 0;
      dPhiClusterMETMB1HitsCRwithNHF50 = 999.;      
      clusterSizeMB1HitsCRinvertedNHF50 = 0;
      dPhiClusterMETMB1HitsCRinvertedNHF50 = 999.;      
      clusterSizeMB1HitsCRwithNHFLead = 0;
      dPhiClusterMETMB1HitsCRwithNHFLead = 999.;      
      clusterSizeMB1HitsCRinvertedNHFLead = 0;
      dPhiClusterMETMB1HitsCRinvertedNHFLead = 999.;      
      clusterSizeMB1Hits30CR = 0;
      dPhiClusterMETMB1Hits30CR = 999.;      
      clusterSizeMB1Hits30MaxMB2CR = 0;
      dPhiClusterMETMB1Hits30MaxMB2CR = 999.;      
      clusterSizeMB1Hits30MaxMB3CR = 0;
      dPhiClusterMETMB1Hits30MaxMB3CR = 999.;      
      clusterSizeMB1Hits30MaxMB4CR = 0;
      dPhiClusterMETMB1Hits30MaxMB4CR = 999.;      
      clusterSizeMB1HitsNoMB1ClusterCR = 0;
      dPhiClusterMETMB1HitsNoMB1ClusterCR = 999.;      
      clusterSizeMB1or2CR = 0;
      dPhiClusterMETMB1or2CR = 999.;      
      clusterSizeMB2CR = 0;
      dPhiClusterMETMB2CR = 999.;      
      clusterSizeMB2CRwithAdjacent = 0;
      dPhiClusterMETMB2CRwithAdjacent = 999.;      
      clusterSizeMB2CRwithAdjacent0p8 = 0;
      dPhiClusterMETMB2CRwithAdjacent0p8 = 999.;      
      clusterSizeMB2CRwithOther = 0;
      dPhiClusterMETMB2CRwithOther = 999.;      
      clusterSizeMB2CRwithNHF50 = 0;
      dPhiClusterMETMB2CRwithNHF50 = 999.;      
      clusterSizeMB2CRinvertedNHF50 = 0;
      dPhiClusterMETMB2CRinvertedNHF50 = 999.;      
      clusterSizeMB2CRwithNHFLead = 0;
      dPhiClusterMETMB2CRwithNHFLead = 999.;      
      clusterSizeMB2CRwithNHFLeadNoRPC = 0;
      dPhiClusterMETMB2CRwithNHFLeadNoRPC = 999.;      
      clusterSizeMB2CRwithInvertedNHFLeadNoRPC = 0;
      dPhiClusterMETMB2CRwithInvertedNHFLeadNoRPC = 999.;      
      clusterSizeMB2CRinvertedNHFLead = 0;
      dPhiClusterMETMB2CRinvertedNHFLead = 999.;      
      clusterSizeMB2withMB1CR = 0;
      dPhiClusterMETMB2withMB1CR = 999.;      
      clusterSizeMB2withMB1CRwithAdjacent = 0;
      dPhiClusterMETMB2withMB1CRwithAdjacent = 999.;      
      clusterSizeMB2withMB1CRwithNHF50 = 0;
      dPhiClusterMETMB2withMB1CRwithNHF50 = 999.;      
      clusterSizeMB2withMB1CRinvertedNHF50 = 0;
      dPhiClusterMETMB2withMB1CRinvertedNHF50 = 999.;      
      clusterSizeMB2withMB1CRwithNHFLead = 0;
      dPhiClusterMETMB2withMB1CRwithNHFLead = 999.;      
      clusterSizeMB2withMB1CRinvertedNHFLead = 0;
      dPhiClusterMETMB2withMB1CRinvertedNHFLead = 999.;      
      nMB1MatchClusterMB1CR = 0;
      nMB1MatchClusterMB1HitsCR = 0;
      nMB1MatchClusterMB1Hits30CR = 0;
      nMB1MatchClusterMB1HitsNoMB1ClusterCR = 0;
      nMB1MatchClusterMB1or2CR = 0;
      nMB1MatchClusterMB2CR = 0;
      nMB1MatchClusterMB2withMB1CR = 0;
      nMB1MatchClusterMB2withMB1CRwithNHFLead = 0;
      nClustersVeto_dPhiJetMET = 0;
      clusterEta.clear();
      clusterPhi.clear();
      clusterSize.clear();
      clusterSizeTotal = 0;

      clusterSizeSR = 0;
      dPhiClusterMETSR = 999.;
      maxStationSR = 0;
      clusterSizeSRwithAdjacent = 0;
      dPhiClusterMETSRwithAdjacent = 999.;
      clusterSizeSRwithAdjacent0p8 = 0;
      dPhiClusterMETSRwithAdjacent0p8 = 999.;
      clusterSizeSRwithOther = 0;
      dPhiClusterMETSRwithOther = 999.;
      maxStationSRwithAdjacent = 0;
      maxStationSRwithOther = 0;
      clusterSizeSRwithNHF50 = 0;
      dPhiClusterMETSRwithNHF50 = 999.;
      clusterSizeSRinvertedNHF50 = 0;
      dPhiClusterMETSRinvertedNHF50 = 999.;
      clusterSizeSRwithNHFLead = 0;
      dPhiClusterMETSRwithNHFLead = 999.;
      clusterSizeSRwithNHFLeadNoRPC = 0;
      dPhiClusterMETSRwithNHFLeadNoRPC = 999.;
      clusterSizeSRwithInvertedNHFLeadNoRPC = 0;
      dPhiClusterMETSRwithInvertedNHFLeadNoRPC = 999.;
      clusterSizeSRinvertedNHFLead = 0;
      dPhiClusterMETSRinvertedNHFLead = 999.;
      passABCD = false;
      passABCDwithAdjacent = false;
      passABCDwithAdjacent0p8 = false;
      passABCDwithOther = false;
      passABCDwithNHF50 = false;
      passABCDinvertedNHF50 = false;
      passABCDwithNHFLead = false;
      passABCDwithNHFLeadNoRPC = false;
      passABCDwithInvertedNHFLeadNoRPC = false;
      passABCDinvertedNHFLead = false;

      passLowClusterMET_rpcCR = false;
      passHighClusterMET_rpcCR = false;
      passLowClusterMET_LowClusterSize = false;
      passHighClusterMET_LowClusterSize = false;
      passHighClusterMET_HighClusterSize = false;

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
      maxGoodClusterSize=0;

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
      
      HLT = false;
      if(year=="2016"){
	if(HLTDecision[310] || HLTDecision[467]){
	  HLT = true;
	}
      }
      else{
	if(HLTDecision[310] || HLTDecision[467] || HLTDecision[703] || HLTDecision[717] || HLTDecision[710] || HLTDecision[709]){
	  HLT = true;
	}
      }
      //if(*nLeptons>0){ HLT = false; }
      if(*MET > 200 && *Flag2_all && HLT){
      //if(*MET > 200 && HLT){
	passMET = true;
	dPhi_min = 999.;
	dPhi_min_invertedJet = 999.;
	dPhiClusterMET = 0.0;
	dPhiClusterMET_max = 0.0;
	chargedHadFraction_mindPhi = -1.0;
	chargedEMFraction_mindPhi = -1.0;
	neutralHadFraction_mindPhi = -1.0;
	neutralEMFraction_mindPhi = -1.0;
	if(*nJets>0){
	  if(jetTightPassId[0] && fabs(jetEta[0])<=2.4){ 
	    passNHFJetLead = true; 
	    h_MET_lowNHF[itr_year]->Fill(*MET);
	    h_leadJetPt_lowNHF[itr_year]->Fill(jetPt[0]);
	    h_leadJetPtMET_lowNHF[itr_year]->Fill(jetPt[0]/(*MET));
	  }
	  else{ 
	    passNHFJetLead = false; 
	    h_MET_highNHF[itr_year]->Fill(*MET);
	    h_leadJetPt_highNHF[itr_year]->Fill(jetPt[0]);
	    h_leadJetPtMET_highNHF[itr_year]->Fill(jetPt[0]/(*MET));
	  }
	}

	for(Int_t itr_dt = 0; itr_dt<*nDtRechits; itr_dt++){
	  //dtPoints.push_back( fastjet::PseudoJet( dtRechitX[itr_dt], dtRechitY[itr_dt], dtRechitZ[itr_dt], 0) );
	  dtR = sqrt(pow(dtRechitX[itr_dt],2)+pow(dtRechitY[itr_dt],2));
	  if(dtR>400. && dtR<480.){ hitStation1+=1; }
	  else if(dtR>485. && dtR<560.){ hitStation2+=1; }
	  else if(dtR>590. && dtR<650.){ hitStation3+=1; }
	  else if(dtR>690. && dtR<800.){ hitStation4+=1; }
	  if(dtRechitZ[itr_dt]>-661. && dtRechitZ[itr_dt]<-395.){ hitWheelm2+=1; }
	  else if(dtRechitZ[itr_dt]>-395. && dtRechitZ[itr_dt]<-127.){ hitWheelm1+=1; }
	  else if(fabs(dtRechitZ[itr_dt]<127.)){ hitWheel0+=1; }
	  else if(dtRechitZ[itr_dt]<395.){ hitWheel1+=1; }
	  else if(dtRechitZ[itr_dt]<661.){ hitWheel2+=1; }
	}
	for(Int_t itr_dt = 0; itr_dt<*nDtSeg; itr_dt++){
	  //dtPoints.push_back( fastjet::PseudoJet( dtRechitX[itr_dt], dtRechitY[itr_dt], dtRechitZ[itr_dt], 0) );
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
	}

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

	if(*nDtRechitClusters>0){
	  nPassNoVeto+=1;
	  for(Int_t itr_clust = 0; itr_clust<*nDtRechitClusters; itr_clust++){
	    dPhiClusterMET = dtRechitClusterPhi[itr_clust] - *METphi;
	    if(dPhiClusterMET > TMath::Pi()){ dPhiClusterMET -= 2*TMath::Pi(); }
	    if(dPhiClusterMET < -1.0*TMath::Pi()){ dPhiClusterMET += 2*TMath::Pi(); }
	    if(fabs(dPhiClusterMET)>dPhiClusterMET_max){ dPhiClusterMET_max=fabs(dPhiClusterMET); }
	    if(dtRechitClusterSize[itr_clust]>maxClusterSize){ maxClusterSize = dtRechitClusterSize[itr_clust]; }
	    if(nStations25<3 && nWheels25<3 && dtRechitClusterSize[itr_clust]>50){
	      h_dtRechitClusterSize[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
	      h_dtRechitClusterR[itr_year]->Fill(sqrt(pow(dtRechitClusterX[itr_clust],2)+pow(dtRechitClusterY[itr_clust],2)));
	      h_dtRechitClusterZ[itr_year]->Fill(dtRechitClusterZ[itr_clust]);
	      h_dtRechitClusterPhi[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
	    }
	  }
	}
	for(Int_t itr_jet = 0; itr_jet<*nJets; itr_jet++){
	  if(fabs(jetEta[itr_jet])<3.0 && jetPt[itr_jet]>30.0){
	    passOneJet = true;
	    goodInvertedJet = true;
	    if(jetNeutralHadronEnergyFraction[itr_jet]>=0.9 && jetPt[itr_jet]>50.0 && fabs(jetEta[itr_jet])<2.4){ passNHFJet50 = false; }
	    for(Int_t itr_clust=0; itr_clust<*nDtRechitClusters; itr_clust++){
	      dPhi_tmp = jetPhi[itr_jet] - dtRechitClusterPhi[itr_clust];
	      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	      if(sqrt(pow(dPhi_tmp,2)+pow(jetEta[itr_jet] - dtRechitClusterEta[itr_clust],2))<0.4){
		goodInvertedJet = false;
	      }
	    }
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
	    if(goodInvertedJet && (fabs(dPhi_tmp) < dPhi_min_invertedJet)){
	      dPhi_min_invertedJet = fabs(dPhi_tmp);
	    }
	  }
	}
	
	if(fabs(dPhiClusterMET_max)>1.0){ nPassClusterCR+=1; }
	if(fabs(dPhi_min)>0.6 && passOneJet){ 
	  //if(passOneJet){ 
	  passJetMET = true;
	  if(*nDtRechitClusters>0){
	    h_jetChargedHadronicEnergyFraction_passJetMET[itr_year]->Fill(chargedHadFraction_mindPhi);
	    h_jetChargedEMEnergyFraction_passJetMET[itr_year]->Fill(chargedEMFraction_mindPhi);
	    h_jetNeutralHadronicEnergyFraction_passJetMET[itr_year]->Fill(neutralHadFraction_mindPhi);
	    h_jetNeutralEMEnergyFraction_passJetMET[itr_year]->Fill(neutralEMFraction_mindPhi);
	  }
	  if(nStations25<3){
	    passStations25 = true;
	    if(nWheels25<3){
	      passWheels25 = true;
	      if(*nDtRechitClusters>0){
		h_jetChargedHadronicEnergyFraction_passStationsWheels[itr_year]->Fill(chargedHadFraction_mindPhi);
		h_jetChargedEMEnergyFraction_passStationsWheels[itr_year]->Fill(chargedEMFraction_mindPhi);
		h_jetNeutralHadronicEnergyFraction_passStationsWheels[itr_year]->Fill(neutralHadFraction_mindPhi);
		h_jetNeutralEMEnergyFraction_passStationsWheels[itr_year]->Fill(neutralEMFraction_mindPhi);
		h_leadingJetChargedHadronicEnergyFraction_passStationsWheels[itr_year]->Fill(jetChargedHadronEnergyFraction[0]);
		h_leadingJetChargedEMEnergyFraction_passStationsWheels[itr_year]->Fill(jetElectronEnergyFraction[0]);
		h_leadingJetNeutralHadronicEnergyFraction_passStationsWheels[itr_year]->Fill(jetNeutralHadronEnergyFraction[0]);
		h_leadingJetNeutralEMEnergyFraction_passStationsWheels[itr_year]->Fill(jetPhotonEnergyFraction[0]);
	      }
	    }
	  }
	}
	else{ h_nDtRechitClusters_dPhiJetMET[itr_year]->Fill(*nDtRechitClusters); }

	/*if(nStations25<3 && nWheels25<3){
	  fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);
	  fastjet::ClusterSequence cs(dtPoints, jet_def);
	  clustersCA = cs.inclusive_jets();
	  for(int i=0; i<clustersCA.size(); i++){
	    constituents.clear();
	    constituents = clustersCA[i].constituents();
	    if(constituents.size() > 50){ 
	      nCAClusters+=1;
	      h_sizeCACluster[itr_year]->Fill(constituents.size());
	      h_radiusCACluster[itr_year]->Fill(sqrt(pow(clustersCA[i].px(),2)+pow(clustersCA[i].py(),2)));
	      cout << "px: " << clustersCA[i].px() << " py: " << clustersCA[i].py() << " R: " << sqrt(pow(clustersCA[i].px(),2)+pow(clustersCA[i].py(),2)) << endl;
	      h_zCACluster[itr_year]->Fill(clustersCA[i].pz());
	      h_phiCACluster[itr_year]->Fill(clustersCA[i].phi_std());
	      
	      fill(CAclusterStationHits, CAclusterStationHits+4, 0);
	      fill(CAclusterWheelHits, CAclusterWheelHits+5, 0);
	      for(int j=0; j<constituents.size(); j++){
		cout << "       px: " << constituents[j].px() << " py: " << constituents[j].py() << " R: " << sqrt(pow(constituents[j].px(),2)+pow(constituents[j].py(),2)) << endl;
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
	  h_nCAClusters[itr_year]->Fill(nCAClusters);
	  h_nDtRechitClusters[itr_year]->Fill(*nDtRechitClusters);

	  if(passJetMET && passNHFJetLead){
	    if(CAclusterSize>50){
	      if(CAclusterdPhiMET>=1.0){
		if(CAclusterSize<100){
		  if(CAclusterStation==2){ SRyieldMB2A_CA+=1; }
		  else if(CAclusterStation>2){ SRyieldA_CA+=1; }
		}
		else{
		  if(CAclusterStation==2){ SRyieldMB2B_CA+=1; }
		  else if(CAclusterStation>2){ SRyieldB_CA+=1; }
		}
	      }
	      else{
		if(CAclusterSize<100){
		  if(CAclusterStation==2){ SRyieldMB2C_CA+=1; }
		  else if(CAclusterStation>2){ SRyieldC_CA+=1; }
		}
	      }
	    }
	  } 
	  }*/

	for(Int_t itr_clust=0; itr_clust<*nDtRechitClusters; itr_clust++){
	    if(nStations25<3 && nWheels25<3 && *nDtRechitClusters>0 && passOneJet){
		if(passNHFJetLead){ 
		    h_MET_oneCluster_lowNHF[itr_year]->Fill(*MET); 
		    h_leadJetPt_oneCluster_lowNHF[itr_year]->Fill(jetPt[0]);
		    h_leadJetPtMET_oneCluster_lowNHF[itr_year]->Fill(jetPt[0]/(*MET));
		    runNumHists["oneCluster_lowNHF_"+year]->Fill(*runNum);
			if(dtRechitClusterSize[itr_clust]>50){
			    etaPhiClusterHists["oneCluster_lowNHF_"+year]->Fill(dtRechitClusterPhi[itr_clust],dtRechitClusterEta[itr_clust]);
			    etaClusterHists["oneCluster_lowNHF_"+year]->Fill(dtRechitClusterEta[itr_clust]);
			    phiClusterHists["oneCluster_lowNHF_"+year]->Fill(dtRechitClusterPhi[itr_clust]);
		    }
		}
		else{ h_MET_oneCluster_highNHF[itr_year]->Fill(*MET); 
		    runNumHists["oneCluster_highNHF_"+year]->Fill(*runNum);
		    h_leadJetPt_oneCluster_highNHF[itr_year]->Fill(jetPt[0]);
		    h_leadJetPtMET_oneCluster_highNHF[itr_year]->Fill(jetPt[0]/(*MET));
			if(dtRechitClusterSize[itr_clust]>50){
			    etaPhiClusterHists["oneCluster_highNHF_"+year]->Fill(dtRechitClusterPhi[itr_clust],dtRechitClusterEta[itr_clust]);
			    etaClusterHists["oneCluster_highNHF_"+year]->Fill(dtRechitClusterEta[itr_clust]);
			    phiClusterHists["oneCluster_highNHF_"+year]->Fill(dtRechitClusterPhi[itr_clust]);
			}
		}
	    }
	  
	  if(dtRechitClusterSize[itr_clust]>50 && ((dtRechitClusterPhi[itr_clust]<0.4 || dtRechitClusterPhi[itr_clust]>0.6) || (*runNum<275.75e3 || *runNum>275.95e3))){
	    passMuon=false;
	    passMuonLoose=false;
	    passMuon_alt=false;
	    passJet=false;
	    rpcBx.clear();
	    rpcMatchStation.clear();
	    rpcSpread = 99;
	    rpcMedian = 99;
	    dPhiClusterRPC = -0.1;
	    dZClusterRPC = -1.;
	    hitsMB1 = 0;
	    invertedJetMB1 = 0;
	    nMB1MatchCluster = 0;
	    nMB2MatchCluster = 0;
	    nMB3MatchCluster = 0;
	    nMB4MatchCluster = 0;
	    nMB1MatchClusterAdjacentPlus = 0;
	    nMB1MatchClusterAdjacentMinus = 0;
	    nMB1MatchClusterAdjacent0p8Plus = 0;
	    nMB1MatchClusterAdjacent0p8Minus = 0;
	    nMB1MatchPi2AdjacentPlus = 0;
	    nMB1MatchPi2AdjacentMinus = 0;
	    nMB1MatchPi2Adjacent0p8Plus = 0;
	    nMB1MatchPi2Adjacent0p8Minus = 0;
	    nRB1MatchCluster = 0;

	    //cout << "doing rpc" << endl;
	    for(Int_t itr_rpc=0; itr_rpc<*nRPCRechits; itr_rpc++){
	      rpcStation=getRPCLayer(RPCRechitX[itr_rpc],RPCRechitY[itr_rpc]);
	      dPhi_tmp = RPCRechitPhi[itr_rpc] - dtRechitClusterPhi[itr_clust];
	      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	      if(fabs(dPhi_tmp)<dPhiClusterRPC || dPhiClusterRPC==-0.1){ dPhiClusterRPC=fabs(dPhi_tmp); }
	      if(fabs(RPCRechitZ[itr_rpc] - dtRechitClusterZ[itr_clust])<dZClusterRPC || dZClusterRPC==-1.){ dZClusterRPC=fabs(RPCRechitZ[itr_rpc]-dtRechitClusterZ[itr_rpc]); }
	      if(fabs(RPCRechitZ[itr_rpc] - dtRechitClusterZ[itr_clust])<5. && fabs(dPhi_tmp)<0.4){
		rpcBx.push_back(RPCRechitBx[itr_rpc]);
		rpcMatchStation.push_back(rpcStation);
	      }
	      /*if(itr_clust==0 && *nRPCRechits<500){
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
	    if(rpcBx.empty() || rpcSpread>0){ passRPCCR=true; }
	    
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
	    if(overlapMuon10GeVLoose){ h_muonOverlap10GeVLoose[itr_year]->Fill(1); }
	    else{ h_muonOverlap10GeVLoose[itr_year]->Fill(0.0); }
	    if(overlapMuon10GeVTight){ h_muonOverlap10GeVTight[itr_year]->Fill(1); }
	    else{ h_muonOverlap10GeVTight[itr_year]->Fill(0.0); }
	    
	    dPhiClusterMET = dtRechitClusterPhi[itr_clust] - *METphi;
	    if(dPhiClusterMET > TMath::Pi()){ dPhiClusterMET -= 2*TMath::Pi(); }
	    if(dPhiClusterMET < -1.0*TMath::Pi()){ dPhiClusterMET += 2*TMath::Pi(); }
	    
	    h_dtRechitClusterJetVetoPt[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust]);
	    h_dtRechitClusterMuonVetoPt[itr_year]->Fill(dtRechitClusterMuonVetoPt[itr_clust]);
	    if(dtRechitClusterJetVetoPt[itr_clust]<20.){ passJet = true; }
	    //if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ passMuon = true; }
	    if(passMuonLoose){ passMuon = true; }
	    if(*nLeptons==0){ passMuon_alt = true; }
	    
	    //passJet = true;
	    overlapJet10GeV = false;
	    overlapJet20GeV = false;
	    matchedJetPhi = -999.;
	    matchedJetEta = -999.;
	    matchedJetNHF = -999.;
	    matchedJetNEF = -999.;
	    matchedJetCHF = -999.;
	    matchedJetCEF = -999.;
	    for(Int_t itr_jet = 0; itr_jet<*nJets; itr_jet++){
	      //if((jetChargedHadronEnergyFraction[itr_jet]+jetElectronEnergyFraction[itr_jet])>0.1){
	      dPhi_tmp = dtRechitClusterPhi[itr_clust] - jetPhi[itr_jet] + pmRand*TMath::Pi()/2.0;
	      //dPhi_tmp = jetPhi[itr_jet] - dtRechitClusterPhi[itr_clust];
	      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	      if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust]-jetEta[itr_jet],2))<0.4){
		if(jetPt[itr_jet]>10. && jetTightPassId[itr_jet]){
		  overlapJet10GeV = true;
		  if(jetPt[itr_jet]>20.){
		    overlapJet20GeV = true;
		  }
		}
	      }
	      dPhi_tmp = jetPhi[itr_jet] - dtRechitClusterPhi[itr_clust];
	      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	      if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust]-jetEta[itr_jet],2))<0.4){
		if(jetPt[itr_jet]>20.){
		  //passJet = false;
		  //if(jetTightPassId[itr_jet]){ passJet = false; }
		  //if(jetChargedHadronEnergyFraction[itr_jet] > 0 || jetElectronEnergyFraction[itr_jet] >0){ passJet = false; }
		  if(matchedJetPhi==-999.){
		    matchedJetPhi = jetPhi[itr_jet];
		    matchedJetEta = jetEta[itr_jet];
		    matchedJetNHF = jetNeutralHadronEnergyFraction[itr_jet];
		    matchedJetCHF = jetChargedHadronEnergyFraction[itr_jet];
		    matchedJetNEF = jetPhotonEnergyFraction[itr_jet];
		    matchedJetCEF = jetElectronEnergyFraction[itr_jet];
		  }
		}
	      }
	    }
	    if(overlapJet10GeV){ h_jetOverlap10GeV[itr_year]->Fill(1); }
	    else{ h_jetOverlap10GeV[itr_year]->Fill(0.0); }
	    if(overlapJet20GeV){ h_jetOverlap20GeV[itr_year]->Fill(1); }
	    else{ h_jetOverlap20GeV[itr_year]->Fill(0.0); }
	    
	    passMB1 = true;
	    /*for(Int_t itr_dt = 0; itr_dt<*nDtRechits; itr_dt++){
	      if(sqrt(pow(dtRechitX[itr_dt],2)+pow(dtRechitY[itr_dt],2))>400. && sqrt(pow(dtRechitX[itr_dt],2)+pow(dtRechitY[itr_dt],2))<480.){
		dPhi_tmp = dtRechitPhi[itr_dt] - dtRechitClusterPhi[itr_clust];
		if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
		if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
		if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitEta[itr_dt]-dtRechitClusterEta[itr_clust],2))<0.4){
		  hitsMB1+=1;
		  //break;
		}
	      }
	      if(itr_clust==0 && *nDtRechits<750){
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
		}
		}*/
	    	    
	    nMB1MatchPi2=0;
	    nMB1MatchCluster=0;
	    nMB2MatchCluster=0;
	    nMB3MatchCluster=0;
	    nMB4MatchCluster=0;
	    for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
	      dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
	      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	      if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtRechitEta[itr_dt],2))<0.4){
		if(dtRechitStation[itr_dt]==1){	nMB1MatchCluster+=1; }
		if(dtRechitStation[itr_dt]==2){	nMB2MatchCluster+=1; }
		if(dtRechitStation[itr_dt]==3){	nMB3MatchCluster+=1; }
		if(dtRechitStation[itr_dt]==4){	nMB4MatchCluster+=1; }
	      }
	      dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt] + pmRand*TMath::Pi()/2.0;
	      if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	      if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	      if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtRechitEta[itr_dt],2))<0.4){
		nMB1MatchPi2+=1;
	      }
	    }
	    hitsMB1=nMB1MatchCluster;

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
	    for(Int_t itr_dt=0; itr_dt<*nDtSeg; itr_dt++){
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
	    }
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

	    //if(nMB1MatchCluster>1 || dtRechitClusterNSegmentStation1[itr_clust]>0){ passMB1 = false; }
	    if(nMB1MatchCluster>1 || nMB1SegMatchCluster>0){ passMB1 = false; }
	    h_dtRechitClusterMB1Veto[itr_year]->Fill(hitsMB1);

	    /*passMuon = true;
	    if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster>0){ passMuon = false; }
	    if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster>0){ passMuon = false; }
	    if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster>0){ passMuon = false; }
	    */
	    if(nStations25<3 && nWheels25<3 && passJet && passMuon && passMB1){
	      if(passNHFJetLead){
		  h_MET_oneVetoCluster_lowNHF[itr_year]->Fill(*MET); 
		  h_leadJetPt_oneVetoCluster_lowNHF[itr_year]->Fill(jetPt[0]);
		  h_leadJetPtMET_oneVetoCluster_lowNHF[itr_year]->Fill(jetPt[0]/(*MET));
		  runNumHists["oneVetoCluster_lowNHF_"+year]->Fill(*runNum);
		  etaPhiClusterHists["oneVetoCluster_lowNHF_"+year]->Fill(dtRechitClusterPhi[itr_clust],dtRechitClusterEta[itr_clust]);
		  etaClusterHists["oneVetoCluster_lowNHF_"+year]->Fill(dtRechitClusterEta[itr_clust]);
		  phiClusterHists["oneVetoCluster_lowNHF_"+year]->Fill(dtRechitClusterPhi[itr_clust]);
	      
	      }
	      else{ 
		  h_MET_oneVetoCluster_highNHF[itr_year]->Fill(*MET); 
		  h_leadJetPt_oneVetoCluster_highNHF[itr_year]->Fill(jetPt[0]);
		  h_leadJetPtMET_oneVetoCluster_highNHF[itr_year]->Fill(jetPt[0]/(*MET));
		  runNumHists["oneVetoCluster_highNHF_"+year]->Fill(*runNum);
		  etaPhiClusterHists["oneVetoCluster_highNHF_"+year]->Fill(dtRechitClusterPhi[itr_clust],dtRechitClusterEta[itr_clust]);
		  etaClusterHists["oneVetoCluster_highNHF_"+year]->Fill(dtRechitClusterEta[itr_clust]);
		  phiClusterHists["oneVetoCluster_highNHF_"+year]->Fill(dtRechitClusterPhi[itr_clust]);
	      }
	    }

	    if(!(passStations25 && passWheels25) && passJetMET && passJet && !rpcBx.empty() && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJetLead){
	      if(passMB1){
		h_MB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(0.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_year]->Fill(0.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_year]->Fill(0.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(0.0); } 
	      }
	      else if(nMB1SegMatchCluster>0 && nMB1MatchCluster>1){ 
		h_MB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(1.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_year]->Fill(1.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_year]->Fill(1.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(1.0); } 
	      }
	      else if(nMB1MatchCluster>1 && nMB1SegMatchCluster==0){ 
		h_MB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(2.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_year]->Fill(2.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_year]->Fill(2.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(2.0); } 
	      }
	      else if(nMB1MatchCluster<=1 && nMB1SegMatchCluster>0){
		h_MB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(3.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_year]->Fill(3.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_year]->Fill(3.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(3.0); } 
	      }
	      else{
		h_MB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(4.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedShowerVetoes[itr_year]->Fill(4.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedShowerVetoes[itr_year]->Fill(4.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(4.0); } 
	      }
	    }
	    if(!(passStations25 && passWheels25) && passJetMET && passJet && passMB1 && !rpcBx.empty() && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJetLead){
	      h_minSegmentDR_invertedShowerVetoes[itr_year]->Fill(minSegmentDR);
	      h_nMatchedSegments_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster);
	      h_matchedSegmentTimeMean_invertedShowerVetoes[itr_year]->Fill(meanSegTime);
	      nMB1SegMatchClusterAdjacent0p8Plus=0;
	      nMB1SegMatchClusterAdjacent0p8Minus=0;
	      for(Int_t itr_dt=0; itr_dt<*nDtSeg; itr_dt++){
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
		      h_segmentAlignmentDeltaPhi_invertedShowerVetoes[itr_year]->Fill(fabs(dPhi_tmp));
		      h_segmentAlignmentDeltaEta_invertedShowerVetoes[itr_year]->Fill(fabs(dtSegEta[itr_dt]-dtSegEta[j]));
		      if(fabs(dPhi_tmp) < 0.28 && fabs(dtSegEta[itr_dt] - dtSegEta[j]) < 0.28){
			nAlignedSegments+=1;
		      }
		    }
		  }
		  h_nAlignedSegments_invertedShowerVetoes[itr_year]->Fill(nAlignedSegments);
		}
	      }
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
		h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(0.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(0.0); } 
	      }
	      else if((nMB1MatchClusterAdjacent0p8Plus>=8 && nMB1SegMatchClusterAdjacent0p8Plus>0) || (nMB1MatchClusterAdjacent0p8Minus>=8 && nMB1SegMatchClusterAdjacent0p8Minus>0)){
		h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(1.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(1.0); } 
	      }
	      else if((nMB1MatchClusterAdjacent0p8Plus>=8 && nMB1SegMatchClusterAdjacent0p8Plus==0) || (nMB1MatchClusterAdjacent0p8Minus>=8 && nMB1SegMatchClusterAdjacent0p8Minus==0)){
		h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(2.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(2.0); } 
	      }
	      else if((nMB1MatchClusterAdjacent0p8Plus<8 && nMB1SegMatchClusterAdjacent0p8Plus>0) || (nMB1MatchClusterAdjacent0p8Minus<8 && nMB1SegMatchClusterAdjacent0p8Minus>0)){
		h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(3.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(3.0); } 
	      }
	      else{
		h_AdjacentMB1VetoOutcomes_invertedShowerVetoes[itr_year]->Fill(4.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(4.0); } 
	      }
	      if(dtRechitClusterMaxStation[itr_clust]==2){ 
		h_nMatchedSegmentsOtherStations_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster); 
		h_nMatchedSegmentsClusterStationMB2_invertedShowerVetoes[itr_year]->Fill(nMB2SegMatchCluster);
		h_nMatchedSegmentsOuterStation_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsOuterStationMB2_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsInnerStationMB2_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster);
		h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		h_matchedSegmentTimeMeanClusterStationMB2_invertedShowerVetoes[itr_year]->Fill(meanMB2SegTime);
		h_matchedSegmentTimeMeanInnerStationMB2_invertedShowerVetoes[itr_year]->Fill(meanMB1SegTime);
		h_matchedSegmentTimeMeanOuterStationMB2_invertedShowerVetoes[itr_year]->Fill(meanMB3SegTime);
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		  h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster); 
		  h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(0.0);
		}
		if(passMuonLoose){ 
		  h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster); 
		  h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(1.0);
		}
		if(nMB1SegMatchCluster!=1){ 
		  h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster); 
		  h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(2.0);
		}
		if(nMB1SegMatchCluster==0){
		  h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(3.0);
		}
		h_muonVetoOutcomesMB2_invertedShowerVetoes[itr_year]->Fill(4.0);
	      }
	      if(dtRechitClusterMaxStation[itr_clust]==3){
		h_nMatchedSegmentsClusterStationMB3_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsOtherStations_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB4SegMatchCluster); 
		h_nMatchedSegmentsInnerStation_invertedShowerVetoes[itr_year]->Fill(nMB2SegMatchCluster);
		h_nMatchedSegmentsOuterStation_invertedShowerVetoes[itr_year]->Fill(nMB4SegMatchCluster);
		h_nMatchedSegmentsInnerStationMB3_invertedShowerVetoes[itr_year]->Fill(nMB2SegMatchCluster);
		h_nMatchedSegmentsOuterStationMB3_invertedShowerVetoes[itr_year]->Fill(nMB4SegMatchCluster);
		h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedShowerVetoes[itr_year]->Fill(nMB2SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		h_matchedSegmentTimeMeanClusterStationMB3_invertedShowerVetoes[itr_year]->Fill(meanMB3SegTime);
		h_matchedSegmentTimeMeanInnerStationMB3_invertedShowerVetoes[itr_year]->Fill(meanMB2SegTime);
		h_matchedSegmentTimeMeanOuterStationMB3_invertedShowerVetoes[itr_year]->Fill(meanMB4SegTime);
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		  h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedShowerVetoes[itr_year]->Fill(nMB2SegMatchCluster); 
		  h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(0.0);
		}
		if(passMuonLoose){ 
		  h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedShowerVetoes[itr_year]->Fill(nMB2SegMatchCluster); 
		  h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(1.0);
		}
		if(nMB2SegMatchCluster!=1){ 
		  h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedShowerVetoes[itr_year]->Fill(nMB2SegMatchCluster); 
		  h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(2.0);
		}
		if(nMB2SegMatchCluster==0){
		  h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(3.0);
		}
		h_muonVetoOutcomesMB3_invertedShowerVetoes[itr_year]->Fill(4.0);
	      }
	      if(dtRechitClusterMaxStation[itr_clust]==4){
		h_nMatchedSegmentsClusterStationMB4_invertedShowerVetoes[itr_year]->Fill(nMB4SegMatchCluster);
		h_nMatchedSegmentsOtherStations_invertedShowerVetoes[itr_year]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster); 
		h_nMatchedSegmentsInnerStation_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsInnerStationMB4_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		if(nMB3SegMatchCluster==0 && nMB2SegMatchCluster==0){ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_year]->Fill(0.0); }
		else if(nMB3SegMatchCluster>0 && nMB2SegMatchCluster==0){ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_year]->Fill(1.0); }
		else if(nMB3SegMatchCluster==0 && nMB2SegMatchCluster>0){ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_year]->Fill(2.0); }
		else if(nMB2SegMatchCluster>0 && nMB3SegMatchCluster>0){ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_year]->Fill(3.0); }
		else{ h_nMatchedSegmentsInnerStationsMB4_invertedShowerVetoes[itr_year]->Fill(4.0); }
		h_matchedSegmentTimeMeanClusterStationMB4_invertedShowerVetoes[itr_year]->Fill(meanMB4SegTime);
		h_matchedSegmentTimeMeanInnerStationMB4_invertedShowerVetoes[itr_year]->Fill(meanMB3SegTime);
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		  h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster); 
		  h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(0.0);
		}
		if(passMuonLoose){ 
		  h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster); 
		  h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(1.0);
		}
		if(nMB3SegMatchCluster!=1){ 
		  h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedShowerVetoes[itr_year]->Fill(nMB3SegMatchCluster); 
		  h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(2.0);
		}
		if(nMB3SegMatchCluster==0){
		  h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(3.0);
		}
		h_muonVetoOutcomesMB4_invertedShowerVetoes[itr_year]->Fill(4.0);
	      }
	    }
	    
	    
	    if(passJetMET && passStations25 && passWheels25){
	      passNoVetoCluster = true;
	      if(nMB1MatchCluster>1){
		passMB1CR = true;
		if(passJet){
		  passJetVeto = true;
		  if(passMuonLoose){ passMuonVetoLoose = true; }
		  if(passMuon){
		    passMuonVeto = true;
		    if(!rpcBx.empty()){
		      passRpcMatchMB1HitsNoMB1ClusterCR = true;
		      clusterEta.push_back(dtRechitClusterEta[itr_clust]);
		      clusterPhi.push_back(dtRechitClusterPhi[itr_clust]);
		      clusterSize.push_back(dtRechitClusterSize[itr_clust]);
		      if(dtRechitClusterMaxStation[itr_clust]>2){
			passRpcMatchMB1HitsCR = true;
			h_jetChargedHadronicEnergyFraction_MB1HitsCR[itr_year]->Fill(chargedHadFraction_mindPhi);
			h_jetChargedEMEnergyFraction_MB1HitsCR[itr_year]->Fill(chargedEMFraction_mindPhi);
			h_jetNeutralHadronicEnergyFraction_MB1HitsCR[itr_year]->Fill(neutralHadFraction_mindPhi);
			h_jetNeutralEMEnergyFraction_MB1HitsCR[itr_year]->Fill(neutralEMFraction_mindPhi);
			if(dtRechitClusterSize[itr_clust]>clusterSizeMB1HitsCR){
			  clusterSizeMB1HitsCR = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB1HitsCR = fabs(dPhiClusterMET);
			  nMB1MatchClusterMB1HitsCR = nMB1MatchCluster;
			}
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
			  passRpcMatchMB1HitsCRwithAdjacent = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB1HitsCRwithAdjacent){
			    clusterSizeMB1HitsCRwithAdjacent = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB1HitsCRwithAdjacent = fabs(dPhiClusterMET);
			  }
			  if(passNHFJet50){
			    passRpcMatchMB1HitsCRwithNHF50 = true;
			    if(dtRechitClusterSize[itr_clust]>clusterSizeMB1HitsCRwithNHF50){
			      clusterSizeMB1HitsCRwithNHF50 = dtRechitClusterSize[itr_clust];
			      dPhiClusterMETMB1HitsCRwithNHF50 = fabs(dPhiClusterMET);
			    }
			    h_rpcBxMedian_MB1HitsCR[itr_year]->Fill(rpcMedian);
			    h_dtRechitClusterSize_MB1HitsCRlowNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			    h_MET_MB1HitsCRlowNHF[itr_year]->Fill(*MET);
			    h_dPhiClusterMET_MB1HitsCRlowNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			  }
			  else{
			    passRpcMatchMB1HitsCRinvertedNHF50 = true;
			    if(dtRechitClusterSize[itr_clust]>clusterSizeMB1HitsCRinvertedNHF50){
			      clusterSizeMB1HitsCRinvertedNHF50 = dtRechitClusterSize[itr_clust];
			      dPhiClusterMETMB1HitsCRinvertedNHF50 = fabs(dPhiClusterMET);
			    }
			    h_dtRechitClusterSize_MB1HitsCRhighNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			    h_MET_MB1HitsCRhighNHF[itr_year]->Fill(*MET);
			    h_dPhiClusterMET_MB1HitsCRhighNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			  }			    
			  if(passNHFJetLead){
			    passRpcMatchMB1HitsCRwithNHFLead = true;
			    if(dtRechitClusterSize[itr_clust]>clusterSizeMB1HitsCRwithNHFLead){
			      clusterSizeMB1HitsCRwithNHFLead = dtRechitClusterSize[itr_clust];
			      dPhiClusterMETMB1HitsCRwithNHFLead = fabs(dPhiClusterMET);
			    }
			    h_leadJetPt_MB1HitsCR[itr_year]->Fill(jetPt[0]);
			    h_leadJetPtMET_MB1HitsCR[itr_year]->Fill(jetPt[0]/(*MET));
			    if(dtRechitClusterSize[itr_clust]<100){
			      h_dtRechitClusterPhi_invertedMB1AC_jetVeto20[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
			      if(dtRechitClusterMaxStation[itr_clust]==3){ h_dtRechitClusterPhi_invertedMB1AC_MB3_jetVeto20[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
			      if(dtRechitClusterMaxStation[itr_clust]==4){ h_dtRechitClusterPhi_invertedMB1AC_MB4_jetVeto20[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
			      if(fabs(dPhiClusterMET)<1.0){ h_dtRechitClusterPhi_invertedMB1C_jetVeto20[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
			      else{ h_dtRechitClusterPhi_invertedMB1A_jetVeto20[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
			      if(dtRechitClusterJetVetoPt[itr_clust]>=10){ 
				h_dtRechitClusterPhi_invertedMB1AC_10GeVJet[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
				if(fabs(dPhiClusterMET)<1.0){ eventListInvertedMB1_10GeVJetsC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
				else{ eventListInvertedMB1_10GeVJetsA << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
				if(dtRechitClusterMaxStation[itr_clust]==3){ h_dtRechitClusterPhi_invertedMB1AC_MB3_10GeVJet[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
				if(dtRechitClusterMaxStation[itr_clust]==4){ h_dtRechitClusterPhi_invertedMB1AC_MB4_10GeVJet[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
				if(fabs(dPhiClusterMET)<1.0){ h_dtRechitClusterPhi_invertedMB1C_10GeVJet[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
				else{ h_dtRechitClusterPhi_invertedMB1A_10GeVJet[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
			      }
			      else{
				h_dtRechitClusterPhi_invertedMB1AC_jetVeto10[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
				if(dtRechitClusterMaxStation[itr_clust]==3){ h_dtRechitClusterPhi_invertedMB1AC_MB3_jetVeto10[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
				if(dtRechitClusterMaxStation[itr_clust]==4){ h_dtRechitClusterPhi_invertedMB1AC_MB4_jetVeto10[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
				if(fabs(dPhiClusterMET)<1.0){ h_dtRechitClusterPhi_invertedMB1C_jetVeto10[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
				else{ h_dtRechitClusterPhi_invertedMB1A_jetVeto10[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
				if(dtRechitClusterPhi[itr_clust]>1.0 && dtRechitClusterPhi[itr_clust]<1.2){ eventListInvertedMB1_phiSpike << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
			      }
			    }
			  }
			  else{
			    passRpcMatchMB1HitsCRinvertedNHFLead = true;
			    if(dtRechitClusterSize[itr_clust]>clusterSizeMB1HitsCRinvertedNHFLead){
			      clusterSizeMB1HitsCRinvertedNHFLead = dtRechitClusterSize[itr_clust];
			      dPhiClusterMETMB1HitsCRinvertedNHFLead = fabs(dPhiClusterMET);
			    }
			  }
			}
		      }
		      if(nMB1MatchCluster>30){
			passRpcMatchMB1Hits30CR = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeMB1Hits30CR){
			  clusterSizeMB1Hits30CR = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB1Hits30CR = fabs(dPhiClusterMET);
			  nMB1MatchClusterMB1Hits30CR = nMB1MatchCluster;
			}
			if(dtRechitClusterMaxStation[itr_clust]==2 && dtRechitClusterSize[itr_clust]>clusterSizeMB1Hits30MaxMB2CR){
			  passRpcMatchMB1Hits30MaxMB2CR = true;
			  clusterSizeMB1Hits30MaxMB2CR = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB1Hits30MaxMB2CR = fabs(dPhiClusterMET);
			}
			if(dtRechitClusterMaxStation[itr_clust]==3 && dtRechitClusterSize[itr_clust]>clusterSizeMB1Hits30MaxMB3CR){
			  passRpcMatchMB1Hits30MaxMB3CR = true;
			  clusterSizeMB1Hits30MaxMB3CR = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB1Hits30MaxMB3CR = fabs(dPhiClusterMET);
			}
			if(dtRechitClusterMaxStation[itr_clust]==4 && dtRechitClusterSize[itr_clust]>clusterSizeMB1Hits30MaxMB4CR){
			  passRpcMatchMB1Hits30MaxMB4CR = true;
			  clusterSizeMB1Hits30MaxMB4CR = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB1Hits30MaxMB4CR = fabs(dPhiClusterMET);
			}
		      }
		      if((dtRechitClusterSize[itr_clust]-dtRechitClusterNSegmentStation1[itr_clust])>clusterSizeMB1HitsNoMB1ClusterCR){
			clusterSizeMB1HitsNoMB1ClusterCR = dtRechitClusterSize[itr_clust]-dtRechitClusterNSegmentStation1[itr_clust];
			dPhiClusterMETMB1HitsNoMB1ClusterCR = fabs(dPhiClusterMET);
			nMB1MatchClusterMB1HitsNoMB1ClusterCR = nMB1MatchCluster;
		      }
		    }
		  }
		}
	      }
	      if(nMB1MatchCluster>1 && dtRechitClusterMaxStation[itr_clust]<=2){
		if(passJet){
		  if(passMuon){
		    if(!rpcBx.empty()){
		      passRpcMatch = true;
		      if(dtRechitClusterSize[itr_clust]>clusterSizeMB1CR){
			clusterSizeMB1CR = dtRechitClusterSize[itr_clust];
			dPhiClusterMETMB1CR = fabs(dPhiClusterMET);
			nMB1MatchClusterMB1CR = nMB1MatchCluster;
		      }
		      h_jetChargedHadronicEnergyFraction_MB1CR[itr_year]->Fill(chargedHadFraction_mindPhi);
		      h_jetChargedEMEnergyFraction_MB1CR[itr_year]->Fill(chargedEMFraction_mindPhi);
		      h_jetNeutralHadronicEnergyFraction_MB1CR[itr_year]->Fill(neutralHadFraction_mindPhi);
		      h_jetNeutralEMEnergyFraction_MB1CR[itr_year]->Fill(neutralEMFraction_mindPhi);
		      if(dtRechitClusterMaxStation[itr_clust]==1){
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
			  if(passNHFJet50){
			    h_dtRechitClusterSize_MB1CRlowNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			    h_MET_MB1CRlowNHF[itr_year]->Fill(*MET);
			    h_dPhiClusterMET_MB1CRlowNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			    h_rpcBxMedian_MB1CR[itr_year]->Fill(rpcMedian);
			  }
			  else{
			    h_dtRechitClusterSize_MB1CRhighNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			    h_MET_MB1CRhighNHF[itr_year]->Fill(*MET);
			    h_dPhiClusterMET_MB1CRhighNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			  }
			  if(passNHFJetLead){
			    passRpcMatchMB1CRwithNHFLead = true;
			    if(dtRechitClusterSize[itr_clust]>clusterSizeMB1CRwithNHFLead){
			      clusterSizeMB1CRwithNHFLead = dtRechitClusterSize[itr_clust];
			      dPhiClusterMETMB1CRwithNHFLead = fabs(dPhiClusterMET);
			    }
			    h_leadJetPt_MB1CR[itr_year]->Fill(jetPt[0]);
			    h_leadJetPtMET_MB1CR[itr_year]->Fill(jetPt[0]/(*MET));
			  }
			}
		      }
		    }
		  }
		}
	      }
	      if(dtRechitClusterMaxStation[itr_clust]<=2){
		if(passJet){
		  if(passMuon){
		    if(!rpcBx.empty()){
		      passRpcMatchMB1or2CR = true;
		      if(dtRechitClusterSize[itr_clust]>clusterSizeMB1or2CR){
			clusterSizeMB1or2CR = dtRechitClusterSize[itr_clust];
			dPhiClusterMETMB1or2CR = fabs(dPhiClusterMET);
			nMB1MatchClusterMB1or2CR = nMB1MatchCluster;
		      }		    
		    }
		  }
		}
	      }
	      if(nMB1MatchCluster<=1 && dtRechitClusterMaxStation[itr_clust]==2){
		if(passJet){
		  if(passMuon){
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
		    if(!rpcBx.empty()){
		      passRpcMatchMB2CR = true;
		      h_jetChargedHadronicEnergyFraction_MB2CR[itr_year]->Fill(chargedHadFraction_mindPhi);
		      h_jetChargedEMEnergyFraction_MB2CR[itr_year]->Fill(chargedEMFraction_mindPhi);
		      h_jetNeutralHadronicEnergyFraction_MB2CR[itr_year]->Fill(neutralHadFraction_mindPhi);
		      h_jetNeutralEMEnergyFraction_MB2CR[itr_year]->Fill(neutralEMFraction_mindPhi);
		      h_nMatchedHitsMB3_fullSelection_MB2CR[itr_year]->Fill(nMB3MatchCluster);
		      h_nMatchedHitsMB4_fullSelection_MB2CR[itr_year]->Fill(nMB4MatchCluster);
		      if(dtRechitClusterSize[itr_clust]>=100){ 
			h_nMatchedHitsMB3and4_highClusterSize_fullSelection_MB2CR[itr_year]->Fill(nMB3MatchCluster);
			h_nMatchedHitsMB3and4_highClusterSize_fullSelection_MB2CR[itr_year]->Fill(nMB4MatchCluster);
			if(fabs(dPhiClusterMET)>=1.0){
			  //eventListMB2CR << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			}
		      }
		      else{
			h_nMatchedHitsMB3and4_lowClusterSize_fullSelection_MB2CR[itr_year]->Fill(nMB3MatchCluster);
			h_nMatchedHitsMB3and4_lowClusterSize_fullSelection_MB2CR[itr_year]->Fill(nMB4MatchCluster);
		      }
		      if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CR){
			clusterSizeMB2CR = dtRechitClusterSize[itr_clust];
			dPhiClusterMETMB2CR = fabs(dPhiClusterMET);
			nMB1MatchClusterMB2CR = nMB1MatchCluster;
		      }
		      if(nMB3MatchCluster<5 && nMB4MatchCluster<5){
			passRpcMatchMB2CRwithOther = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRwithOther){
			  clusterSizeMB2CRwithOther = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB2CRwithOther = fabs(dPhiClusterMET);
			}
		      }
		      if(nMB1MatchClusterAdjacentPlus<5 && nMB1MatchClusterAdjacentMinus<5){
			passRpcMatchMB2CRwithAdjacent = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRwithAdjacent){
			  clusterSizeMB2CRwithAdjacent = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB2CRwithAdjacent = fabs(dPhiClusterMET);
			}
			h_nMatchedHitsMB3_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year]->Fill(nMB3MatchCluster);
			h_nMatchedHitsMB4_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year]->Fill(nMB4MatchCluster);
			if(dtRechitClusterSize[itr_clust]>=100){ 
			  h_nMatchedHitsMB3and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year]->Fill(nMB3MatchCluster);
			  h_nMatchedHitsMB3and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year]->Fill(nMB4MatchCluster);
			}
			else{
			  h_nMatchedHitsMB3and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year]->Fill(nMB3MatchCluster);
			  h_nMatchedHitsMB3and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year]->Fill(nMB4MatchCluster);
			}
		      }
		      if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
			passRpcMatchMB2CRwithAdjacent0p8 = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRwithAdjacent0p8){
			  clusterSizeMB2CRwithAdjacent0p8 = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB2CRwithAdjacent0p8 = fabs(dPhiClusterMET);
			}
			if(passNHFJet50){
			  passRpcMatchMB2CRwithNHF50 = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRwithNHF50){
			    clusterSizeMB2CRwithNHF50 = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB2CRwithNHF50 = fabs(dPhiClusterMET);
			  }
			  if(dtRechitClusterSize[itr_clust]<100 || fabs(dPhiClusterMET)>=1.0){
			    h_rpcBxMedian_MB2CR[itr_year]->Fill(rpcMedian);
			    h_dtRechitClusterSize_MB2CRlowNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			    h_MET_MB2CRlowNHF[itr_year]->Fill(*MET);
			    h_dPhiClusterMET_MB2CRlowNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			  }
			}	
			else{
			  passRpcMatchMB2CRinvertedNHF50 = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRinvertedNHF50){
			    clusterSizeMB2CRinvertedNHF50 = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB2CRinvertedNHF50 = fabs(dPhiClusterMET);
			  }
			  h_dtRechitClusterSize_MB2CRhighNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			  h_MET_MB2CRhighNHF[itr_year]->Fill(*MET);
			  h_dPhiClusterMET_MB2CRhighNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			}
			if(passNHFJetLead){
			  passRpcMatchMB2CRwithNHFLead = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRwithNHFLead){
			    clusterSizeMB2CRwithNHFLead = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB2CRwithNHFLead = fabs(dPhiClusterMET);
			  }
			  if(dtRechitClusterSize[itr_clust]<100 || fabs(dPhiClusterMET)>=1.0){
			    //eventListABC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			  }
			  h_leadJetPt_MB2CR[itr_year]->Fill(jetPt[0]);
			  h_leadJetPtMET_MB2CR[itr_year]->Fill(jetPt[0]/(*MET));
			}
			else{
			  passRpcMatchMB2CRinvertedNHFLead = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRinvertedNHFLead){
			    clusterSizeMB2CRinvertedNHFLead = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB2CRinvertedNHFLead = fabs(dPhiClusterMET);
			  }
			}
		      }
		    }
		    if(rpcBx.empty() && nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
		      if(passNHFJetLead){
			noRpcMatchMB2CRwithNHFLead = true;
			//eventListNoRPC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRwithNHFLeadNoRPC){
			  clusterSizeMB2CRwithNHFLeadNoRPC = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB2CRwithNHFLeadNoRPC = fabs(dPhiClusterMET);
			}
		      }
		    }
		    if(nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
		      if(!passNHFJetLead){
			passMB2CRwithInvertedNHFLeadNoRPC = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeMB2CRwithInvertedNHFLeadNoRPC){
			  clusterSizeMB2CRwithInvertedNHFLeadNoRPC = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB2CRwithInvertedNHFLeadNoRPC = fabs(dPhiClusterMET);
			}
		      }
		    }
		  }
		}
	      }
	      if(dtRechitClusterMaxStation[itr_clust]==2){
		if(!rpcBx.empty() && passJet){
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
		    if(passNHFJetLead){
		      if(passMB1){ 
			h_clusterSize_looseMuonVetoResult_ABCD[itr_year]->Fill(dtRechitClusterSize[itr_clust],passMuonLoose); 
			h_clusterSize_recoMuonVetoResult_ABCD[itr_year]->Fill(dtRechitClusterSize[itr_clust],dtRechitClusterMuonVetoPt[itr_clust]<10.); 
		      }
		      if(dtRechitClusterSize[itr_clust]>=100){
			if(fabs(dPhiClusterMET)<1.0){ 
			  if(dtRechitClusterMuonVetoPt[itr_clust]>=10){ h_matchedMB1Hits_invertedRecoMuonVeto_D[itr_year]->Fill(nMB1MatchCluster); }
			  if(!passMuonLoose){ h_matchedMB1Hits_invertedLooseMuonVeto_D[itr_year]->Fill(nMB1MatchCluster); }
			}
			else{ 
			  if(dtRechitClusterMuonVetoPt[itr_clust]>=10){ h_matchedMB1Hits_invertedRecoMuonVeto_B[itr_year]->Fill(nMB1MatchCluster); }
			  if(!passMuonLoose){ h_matchedMB1Hits_invertedLooseMuonVeto_B[itr_year]->Fill(nMB1MatchCluster); }
			}
		      }
		      else{
			if(fabs(dPhiClusterMET)<1.0){ 
			  if(dtRechitClusterMuonVetoPt[itr_clust]>=10){ h_matchedMB1Hits_invertedRecoMuonVeto_C[itr_year]->Fill(nMB1MatchCluster); }
			  if(!passMuonLoose){ h_matchedMB1Hits_invertedLooseMuonVeto_C[itr_year]->Fill(nMB1MatchCluster); }
			}
			else{ 
			  if(dtRechitClusterMuonVetoPt[itr_clust]>=10){ h_matchedMB1Hits_invertedRecoMuonVeto_A[itr_year]->Fill(nMB1MatchCluster); }
			  if(!passMuonLoose){ h_matchedMB1Hits_invertedLooseMuonVeto_A[itr_year]->Fill(nMB1MatchCluster); }
			}
		      }
		    }
		  }
		}
	      }
	      if(nMB1MatchCluster>1 && dtRechitClusterMaxStation[itr_clust]==2){
		if(passJet){
		  if(passMuon){
		    if(!rpcBx.empty()){
		      passRpcMatchMB2withMB1CR = true;
		      if(dtRechitClusterSize[itr_clust]>clusterSizeMB2withMB1CR){
			clusterSizeMB2withMB1CR = dtRechitClusterSize[itr_clust];
			dPhiClusterMETMB2withMB1CR = fabs(dPhiClusterMET);
			nMB1MatchClusterMB2withMB1CR = nMB1MatchCluster;
		      }
		      h_jetChargedHadronicEnergyFraction_MB2withMB1CR[itr_year]->Fill(chargedHadFraction_mindPhi);
		      h_jetChargedEMEnergyFraction_MB2withMB1CR[itr_year]->Fill(chargedEMFraction_mindPhi);
		      h_jetNeutralHadronicEnergyFraction_MB2withMB1CR[itr_year]->Fill(neutralHadFraction_mindPhi);
		      h_jetNeutralEMEnergyFraction_MB2withMB1CR[itr_year]->Fill(neutralEMFraction_mindPhi);
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
			passRpcMatchMB2withMB1CRwithAdjacent = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeMB2withMB1CRwithAdjacent){
			  clusterSizeMB2withMB1CRwithAdjacent = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETMB2withMB1CRwithAdjacent = fabs( dPhiClusterMET);
			}
			if(passNHFJet50){
			  passRpcMatchMB2withMB1CRwithNHF50 = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB2withMB1CRwithNHF50){
			    clusterSizeMB2withMB1CRwithNHF50 = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB2withMB1CRwithNHF50 = fabs(dPhiClusterMET);
			  }
			  h_rpcBxMedian_MB2withMB1CR[itr_year]->Fill(rpcMedian);
			  h_dtRechitClusterSize_MB2withMB1CRlowNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			  h_MET_MB2withMB1CRlowNHF[itr_year]->Fill(*MET);
			  h_dPhiClusterMET_MB2withMB1CRlowNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			}
			else{
			  passRpcMatchMB2withMB1CRinvertedNHF50 = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB2withMB1CRinvertedNHF50){
			    clusterSizeMB2withMB1CRinvertedNHF50 = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB2withMB1CRinvertedNHF50 = fabs(dPhiClusterMET);
			  }
			  h_dtRechitClusterSize_MB2withMB1CRhighNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			  h_MET_MB2withMB1CRhighNHF[itr_year]->Fill(*MET);
			  h_dPhiClusterMET_MB2withMB1CRhighNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			}
			if(passNHFJetLead){
			  passRpcMatchMB2withMB1CRwithNHFLead = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB2withMB1CRwithNHFLead){
			    clusterSizeMB2withMB1CRwithNHFLead = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB2withMB1CRwithNHFLead = fabs(dPhiClusterMET);
			    nMB1MatchClusterMB2withMB1CRwithNHFLead = nMB1MatchCluster;
			    if(dtRechitClusterSize[itr_clust]>=100 && fabs(dPhiClusterMET)<1.0){ cout << "SR! " << nMB1MatchClusterMB2withMB1CRwithNHFLead << " " << invertedMB1_MB2Cluster_5MB1 << endl; }
			    //cout << nMB1MatchCluster;
			    //
			    //else{ cout << endl; }
			  }
			  h_leadJetPt_MB2withMB1CR[itr_year]->Fill(jetPt[0]);
			  h_leadJetPtMET_MB2withMB1CR[itr_year]->Fill(jetPt[0]/(*MET));
			  if(dtRechitClusterSize[itr_clust] >= 100){
			    if(fabs(dPhiClusterMET) < 1.0){
			      //eventListInvertedMB1LooseMuon << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			      h_matchedMB1Hits_looseMuonVeto_D[itr_year]->Fill(nMB1MatchCluster);
			      h_clusterSize_looseMuonVeto_D[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			      h_matchedMB1Segments_looseMuonVeto_D[itr_year]->Fill(nMB1SegMatchCluster);
			      h_minMB1SegmentDR_looseMuonVeto_D[itr_year]->Fill(minMB1SegmentDR);
			      h_matchedMB1HitsRatio_looseMuonVeto_D[itr_year]->Fill(nMB1MatchCluster/float(dtRechitClusterSize[itr_clust]));
			    }
			    else{
			      h_matchedMB1Hits_looseMuonVeto_B[itr_year]->Fill(nMB1MatchCluster);
			      h_clusterSize_looseMuonVeto_B[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			      h_matchedMB1Segments_looseMuonVeto_B[itr_year]->Fill(nMB1SegMatchCluster);
			      h_minMB1SegmentDR_looseMuonVeto_B[itr_year]->Fill(minMB1SegmentDR);
			      h_matchedMB1HitsRatio_looseMuonVeto_B[itr_year]->Fill(nMB1MatchCluster/float(dtRechitClusterSize[itr_clust]));
			    }
			  }
			  else{
			    if(fabs(dPhiClusterMET) < 1.0){
			      h_matchedMB1Hits_looseMuonVeto_C[itr_year]->Fill(nMB1MatchCluster);
			      h_clusterSize_looseMuonVeto_C[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			      if(nMB1MatchCluster>=10 && nMB1MatchCluster<=20){ h_clusterSizeMuonMB1_looseMuonVeto_C[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
			      h_matchedMB1Segments_looseMuonVeto_C[itr_year]->Fill(nMB1SegMatchCluster);
			      h_minMB1SegmentDR_looseMuonVeto_C[itr_year]->Fill(minMB1SegmentDR);
			      h_matchedMB1HitsRatio_looseMuonVeto_C[itr_year]->Fill(nMB1MatchCluster/float(dtRechitClusterSize[itr_clust]));
			    }
			    else{
			      h_matchedMB1Hits_looseMuonVeto_A[itr_year]->Fill(nMB1MatchCluster);
			      h_clusterSize_looseMuonVeto_A[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			      h_matchedMB1Segments_looseMuonVeto_A[itr_year]->Fill(nMB1SegMatchCluster);
			      h_minMB1SegmentDR_looseMuonVeto_A[itr_year]->Fill(minMB1SegmentDR);
			      h_matchedMB1HitsRatio_looseMuonVeto_A[itr_year]->Fill(nMB1MatchCluster/float(dtRechitClusterSize[itr_clust]));
			    }
			  }
			}
			else{
			  passRpcMatchMB2withMB1CRinvertedNHFLead = true;
			  if(dtRechitClusterSize[itr_clust]>clusterSizeMB2withMB1CRinvertedNHFLead){
			    clusterSizeMB2withMB1CRinvertedNHFLead = dtRechitClusterSize[itr_clust];
			    dPhiClusterMETMB2withMB1CRinvertedNHFLead = fabs(dPhiClusterMET);
			  }
			  if(dtRechitClusterSize[itr_clust] >= 100){
			    if(fabs(dPhiClusterMET) < 1.0){
			      h_matchedMB1Hits_looseMuonVetoDoubleInverted_D[itr_year]->Fill(nMB1MatchCluster);
			      h_matchedMB1Segments_looseMuonVetoDoubleInverted_D[itr_year]->Fill(nMB1SegMatchCluster);
			      h_minMB1SegmentDR_looseMuonVetoDoubleInverted_D[itr_year]->Fill(minMB1SegmentDR);
			      h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_D[itr_year]->Fill(nMB1MatchCluster/float(dtRechitClusterSize[itr_clust]));
			    }
			    else{
			      h_matchedMB1Hits_looseMuonVetoDoubleInverted_B[itr_year]->Fill(nMB1MatchCluster);
			      h_matchedMB1Segments_looseMuonVetoDoubleInverted_B[itr_year]->Fill(nMB1SegMatchCluster);
			      h_minMB1SegmentDR_looseMuonVetoDoubleInverted_B[itr_year]->Fill(minMB1SegmentDR);
			      h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_B[itr_year]->Fill(nMB1MatchCluster/float(dtRechitClusterSize[itr_clust]));
			    }
			  }
			  else{
			    if(fabs(dPhiClusterMET) < 1.0){
			      h_matchedMB1Hits_looseMuonVetoDoubleInverted_C[itr_year]->Fill(nMB1MatchCluster);
			      h_matchedMB1Segments_looseMuonVetoDoubleInverted_C[itr_year]->Fill(nMB1SegMatchCluster);
			      h_minMB1SegmentDR_looseMuonVetoDoubleInverted_C[itr_year]->Fill(minMB1SegmentDR);
			      h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_C[itr_year]->Fill(nMB1MatchCluster/float(dtRechitClusterSize[itr_clust]));
			    }
			    else{
			      h_matchedMB1Hits_looseMuonVetoDoubleInverted_A[itr_year]->Fill(nMB1MatchCluster);
			      h_matchedMB1Segments_looseMuonVetoDoubleInverted_A[itr_year]->Fill(nMB1SegMatchCluster);
			      h_minMB1SegmentDR_looseMuonVetoDoubleInverted_A[itr_year]->Fill(minMB1SegmentDR);
			      h_matchedMB1HitsRatio_looseMuonVetoDoubleInverted_A[itr_year]->Fill(nMB1MatchCluster/float(dtRechitClusterSize[itr_clust]));
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	      if(nMB1MatchCluster<=1 && dtRechitClusterMaxStation[itr_clust]>2){
		if(passJet){
		  if(passMuon){
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
		    if(!rpcBx.empty()){
		      passABCD = true;
		      if(nMB1MatchClusterAdjacentPlus<5 && nMB1MatchClusterAdjacentMinus<5){
			passABCDwithAdjacent = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeSRwithAdjacent){
			  clusterSizeSRwithAdjacent = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETSRwithAdjacent = fabs(dPhiClusterMET);
			  maxStationSRwithAdjacent = dtRechitClusterMaxStation[itr_clust];
			}
		      }
		      if(nMB1MatchClusterAdjacent0p8Plus<6 && nMB1MatchClusterAdjacent0p8Minus<6){
			passABCDwithAdjacent0p8 = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeSRwithAdjacent0p8){
			  clusterSizeSRwithAdjacent0p8 = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETSRwithAdjacent0p8 = fabs(dPhiClusterMET);
			}
		      }
		      if(dtRechitClusterSize[itr_clust]<100 || fabs(dPhiClusterMET)>=1.0){
			if(dtRechitClusterSize[itr_clust]>=100 && fabs(dPhiClusterMET)>=1.5){
			  //eventListSRB << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			}
			if(dtRechitClusterMaxStation[itr_clust]==3){
			  h_nMatchedHitsMB2_fullSelection_SRMB3[itr_year]->Fill(nMB2MatchCluster);
			  h_nMatchedHitsMB4_fullSelection_SRMB3[itr_year]->Fill(nMB4MatchCluster);
			  if(nMB1MatchClusterAdjacent0p8Plus<6 && nMB1MatchClusterAdjacent0p8Minus<6){
			    h_nMatchedHitsMB2_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year]->Fill(nMB2MatchCluster);
			    h_nMatchedHitsMB4_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year]->Fill(nMB4MatchCluster);
			  }
			  if(dtRechitClusterSize[itr_clust]>=100){ 
			    h_nMatchedHitsMB2and4_highClusterSize_fullSelection_SRMB3[itr_year]->Fill(nMB2MatchCluster);
			    h_nMatchedHitsMB2and4_highClusterSize_fullSelection_SRMB3[itr_year]->Fill(nMB4MatchCluster);
			    if(nMB1MatchClusterAdjacent0p8Plus<6 && nMB1MatchClusterAdjacent0p8Minus<6){
			      h_nMatchedHitsMB2and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year]->Fill(nMB2MatchCluster);
			      h_nMatchedHitsMB2and4_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year]->Fill(nMB4MatchCluster);
			    }
			  }
			  else{
			    h_nMatchedHitsMB2and4_lowClusterSize_fullSelection_SRMB3[itr_year]->Fill(nMB2MatchCluster);
			    h_nMatchedHitsMB2and4_lowClusterSize_fullSelection_SRMB3[itr_year]->Fill(nMB4MatchCluster);
			    if(nMB1MatchClusterAdjacent0p8Plus<6 && nMB1MatchClusterAdjacent0p8Minus<6){
			      h_nMatchedHitsMB2and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year]->Fill(nMB2MatchCluster);
			      h_nMatchedHitsMB2and4_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year]->Fill(nMB4MatchCluster);
			    }
			  }
			}
			else if(dtRechitClusterMaxStation[itr_clust]==4){
			  h_nMatchedHitsMB2_fullSelection_SRMB4[itr_year]->Fill(nMB2MatchCluster);
			  h_nMatchedHitsMB3_fullSelection_SRMB4[itr_year]->Fill(nMB3MatchCluster);
			  if(nMB1MatchClusterAdjacent0p8Plus<6 && nMB1MatchClusterAdjacent0p8Minus<6){
			    h_nMatchedHitsMB2_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year]->Fill(nMB2MatchCluster);
			    h_nMatchedHitsMB3_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year]->Fill(nMB3MatchCluster);
			  }
			  if(dtRechitClusterSize[itr_clust]>=100){ 
			    h_nMatchedHitsMB2and3_highClusterSize_fullSelection_SRMB4[itr_year]->Fill(nMB2MatchCluster);
			    h_nMatchedHitsMB2and3_highClusterSize_fullSelection_SRMB4[itr_year]->Fill(nMB3MatchCluster);
			    if(nMB1MatchClusterAdjacent0p8Plus<6 && nMB1MatchClusterAdjacent0p8Minus<6){
			      h_nMatchedHitsMB2and3_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year]->Fill(nMB2MatchCluster);
			      h_nMatchedHitsMB2and3_highClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year]->Fill(nMB3MatchCluster);
			    }
			  }
			  else{
			    h_nMatchedHitsMB2and3_lowClusterSize_fullSelection_SRMB4[itr_year]->Fill(nMB2MatchCluster);
			    h_nMatchedHitsMB2and3_lowClusterSize_fullSelection_SRMB4[itr_year]->Fill(nMB3MatchCluster);
			    if(nMB1MatchClusterAdjacent0p8Plus<6 && nMB1MatchClusterAdjacent0p8Minus<6){
			      h_nMatchedHitsMB2and3_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year]->Fill(nMB2MatchCluster);
			      h_nMatchedHitsMB2and3_lowClusterSize_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year]->Fill(nMB3MatchCluster);
			    }
			  }
			}
			h_jetChargedHadronicEnergyFraction_SR[itr_year]->Fill(chargedHadFraction_mindPhi);
			h_jetChargedEMEnergyFraction_SR[itr_year]->Fill(chargedEMFraction_mindPhi);
			h_jetNeutralHadronicEnergyFraction_SR[itr_year]->Fill(neutralHadFraction_mindPhi);
			h_jetNeutralEMEnergyFraction_SR[itr_year]->Fill(neutralEMFraction_mindPhi);
			h_leadingJetChargedHadronicEnergyFraction_SR[itr_year]->Fill(jetChargedHadronEnergyFraction[0]);
			h_leadingJetChargedEMEnergyFraction_SR[itr_year]->Fill(jetElectronEnergyFraction[0]);
			h_leadingJetNeutralHadronicEnergyFraction_SR[itr_year]->Fill(jetNeutralHadronEnergyFraction[0]);
			h_leadingJetNeutralEMEnergyFraction_SR[itr_year]->Fill(jetPhotonEnergyFraction[0]);
			if(neutralHadFraction_mindPhi>=0.9){
		
			}
		      }
		      if(dtRechitClusterSize[itr_clust]>clusterSizeSR){
			clusterSizeSR = dtRechitClusterSize[itr_clust];
			dPhiClusterMETSR = fabs(dPhiClusterMET);
		      }
		      if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2MatchCluster<5 && nMB4MatchCluster<5){
			passABCDwithOther = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeSRwithOther){
			  clusterSizeSRwithOther = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETSRwithOther = fabs(dPhiClusterMET);
			  maxStationSRwithOther = dtRechitClusterMaxStation[itr_clust];
			}
		      }
		      if(dtRechitClusterMaxStation[itr_clust]==4 && nMB2MatchCluster<5 && nMB3MatchCluster<5){
			passABCDwithOther = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeSRwithOther){
			  clusterSizeSRwithOther = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETSRwithOther = fabs(dPhiClusterMET);
			  maxStationSRwithOther = dtRechitClusterMaxStation[itr_clust];
			}
		      }
		      if(passNHFJet50 && nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
			passABCDwithNHF50 = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeSRwithNHF50){
			  clusterSizeSRwithNHF50 = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETSRwithNHF50 = fabs(dPhiClusterMET);
			}
			if(dtRechitClusterSize[itr_clust]<100 || fabs(dPhiClusterMET)>=1.0){
			  h_rpcBxMedian_SR[itr_year]->Fill(rpcMedian);
			  h_dtRechitClusterSize_SRlowNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			  h_MET_SRlowNHF[itr_year]->Fill(*MET);
			  h_dPhiClusterMET_SRlowNHF[itr_year]->Fill(fabs(dPhiClusterMET));
			}
		      }
		      if(!passNHFJet50 && nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus< 8){
			passABCDinvertedNHF50 = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeSRinvertedNHF50){
			  clusterSizeSRinvertedNHF50 = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETSRinvertedNHF50 = fabs(dPhiClusterMET);
			}
			h_dtRechitClusterSize_SRhighNHF[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			h_MET_SRhighNHF[itr_year]->Fill(*MET);
			h_dPhiClusterMET_SRhighNHF[itr_year]->Fill(fabs(dPhiClusterMET));
		      }		      
		      if(passNHFJetLead && nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
			passABCDwithNHFLead = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeSRwithNHFLead){
			  clusterSizeSRwithNHFLead = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETSRwithNHFLead = fabs(dPhiClusterMET);
			  maxStationSR = dtRechitClusterMaxStation[itr_clust];
			}
			if(dtRechitClusterSize[itr_clust]<100 || fabs(dPhiClusterMET)>=1.0){
			  //eventListABC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			}
			h_leadJetPt_SR[itr_year]->Fill(jetPt[0]);
			h_leadJetPtMET_SR[itr_year]->Fill(jetPt[0]/(*MET));
			if(dtRechitClusterSize[itr_clust]<100){
			  h_dtRechitClusterPhi_AC_jetVeto20[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
			  if(dtRechitClusterJetVetoPt[itr_clust]>=10){
			    h_dtRechitClusterPhi_AC_10GeVJet[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
			    if(fabs(dPhiClusterMET)<1.0){ eventListABCD_10GeVJetsC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
			    else{ eventListABCD_10GeVJetsA << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
			  }
			  else{ h_dtRechitClusterPhi_AC_jetVeto10[itr_year]->Fill(dtRechitClusterPhi[itr_clust]); }
			}
		      }
		      if(!passNHFJetLead && nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
			passABCDinvertedNHFLead = true;
			if(dtRechitClusterSize[itr_clust]>clusterSizeSRinvertedNHFLead){
			  clusterSizeSRinvertedNHFLead = dtRechitClusterSize[itr_clust];
			  dPhiClusterMETSRinvertedNHFLead = fabs(dPhiClusterMET);
			}
		      }
		    }
		    if(rpcBx.empty() && passNHFJetLead && nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
		      passABCDwithNHFLeadNoRPC = true;
		      //eventListNoRPC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
		      if(dtRechitClusterSize[itr_clust]>clusterSizeSRwithNHFLeadNoRPC){
			clusterSizeSRwithNHFLeadNoRPC = dtRechitClusterSize[itr_clust];
			dPhiClusterMETSRwithNHFLeadNoRPC = fabs(dPhiClusterMET);
		      }
		    }
		    if(!passNHFJetLead && nMB1MatchClusterAdjacent0p8Plus<8 && nMB1MatchClusterAdjacent0p8Minus<8){
		      passABCDwithInvertedNHFLeadNoRPC = true;
		      if(dtRechitClusterSize[itr_clust]>clusterSizeSRwithInvertedNHFLeadNoRPC){
			clusterSizeSRwithInvertedNHFLeadNoRPC = dtRechitClusterSize[itr_clust];
			dPhiClusterMETSRwithInvertedNHFLeadNoRPC = fabs(dPhiClusterMET);
		      }
		    }
		  }
		}
	      }
	    }

	    if(passStations25 && passWheels25 && passJetMET && !passJet && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJetLead){
	      if(nMB1MatchPi2<2){ h_nMB1MatchPi2_dPhiClusterMET[itr_year]->Fill(0.0); }
	      else{ h_nMB1MatchPi2_dPhiClusterMET[itr_year]->Fill(1.0); }
	    }

	    if(passStations25 && passWheels25 && passJetMET && !passJet && !rpcBx.empty() && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJetLead){
	      if(passMB1){
		h_MB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(0.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_year]->Fill(0.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_year]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_year]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_year]->Fill(0.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_year]->Fill(0.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_year]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_year]->Fill(0.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_year]->Fill(0.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_year]->Fill(0.0); } 
	      }
	      else if(nMB1SegMatchCluster>0 && nMB1MatchCluster>1){ 
		h_MB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(1.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_year]->Fill(1.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_year]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_year]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_year]->Fill(1.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_year]->Fill(1.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_year]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_year]->Fill(1.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_year]->Fill(1.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_year]->Fill(1.0); } 
	      }
	      else if(nMB1SegMatchCluster==0 && nMB1MatchCluster>1){ 
		h_MB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(2.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_year]->Fill(2.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_year]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_year]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_year]->Fill(2.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_year]->Fill(2.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_year]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_year]->Fill(2.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_year]->Fill(2.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_year]->Fill(2.0); } 
	      }
	      else if(nMB1SegMatchCluster>0 && nMB1MatchCluster<=1){
		h_MB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(3.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_year]->Fill(3.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_year]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_year]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_year]->Fill(3.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_year]->Fill(3.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_year]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_year]->Fill(3.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_year]->Fill(3.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_year]->Fill(3.0); } 
	      }
	      else{
		h_MB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(4.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){
		  h_MB1VetoOutcomesPassRecoMuon_invertedJetVeto[itr_year]->Fill(4.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassRecoMuonMB2_invertedJetVeto[itr_year]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassRecoMuonMB3_invertedJetVeto[itr_year]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassRecoMuonMB4_invertedJetVeto[itr_year]->Fill(4.0); } 
		}
		if(passMuonLoose){
		  h_MB1VetoOutcomesPassLooseMuon_invertedJetVeto[itr_year]->Fill(4.0); 
		  if(dtRechitClusterMaxStation[itr_clust]==2){ h_MB1VetoOutcomesPassLooseMuonMB2_invertedJetVeto[itr_year]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==3){ h_MB1VetoOutcomesPassLooseMuonMB3_invertedJetVeto[itr_year]->Fill(4.0); } 
		  if(dtRechitClusterMaxStation[itr_clust]==4){ h_MB1VetoOutcomesPassLooseMuonMB4_invertedJetVeto[itr_year]->Fill(4.0); } 
		}
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB2_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB3_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster!=1){ h_MB1VetoOutcomesPassOneSegMuonMB4_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==2 && nMB1SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB2_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3 && nMB2SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB3_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4 && nMB3SegMatchCluster==0){ h_MB1VetoOutcomesPassSegMuonMB4_invertedJetVeto[itr_year]->Fill(4.0); } 
	      }
	      if(!passMB1){
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_nMatchedSegmentsInnerStationMB2_invertedMB1_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster); }
		else if(dtRechitClusterMaxStation[itr_clust]==3){ h_nMatchedSegmentsInnerStationMB3_invertedMB1_invertedJetVeto[itr_year]->Fill(nMB2SegMatchCluster); }
		else if(dtRechitClusterMaxStation[itr_clust]==4){ h_nMatchedSegmentsInnerStationMB4_invertedMB1_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster); }
	      }
	    }
	    if(passStations25 && passWheels25 && passJetMET && !passJet && passMB1 && !rpcBx.empty() && dtRechitClusterMaxStation[itr_clust]>=2 && passNHFJetLead){
	      h_minSegmentDR_invertedJetVeto[itr_year]->Fill(minSegmentDR);
	      h_nMatchedSegments_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster);
	      h_matchedSegmentTimeMean_invertedJetVeto[itr_year]->Fill(meanSegTime);
	      nMB1SegMatchClusterAdjacent0p8Plus=0;
	      nMB1SegMatchClusterAdjacent0p8Minus=0;
	      for(Int_t itr_dt=0; itr_dt<*nDtSeg; itr_dt++){
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
		      h_segmentAlignmentDeltaPhi_invertedJetVeto[itr_year]->Fill(fabs(dPhi_tmp));
		      h_segmentAlignmentDeltaEta_invertedJetVeto[itr_year]->Fill(fabs(dtSegEta[itr_dt]-dtSegEta[j]));
		      if(fabs(dPhi_tmp) < 0.28 && fabs(dtSegEta[itr_dt] - dtSegEta[j]) < 0.28){
			nAlignedSegments+=1;
		      }
		    }
		  }
		  h_nAlignedSegments_invertedJetVeto[itr_year]->Fill(nAlignedSegments);
		}
	      }
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
		h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(0.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(0.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(0.0); } 
	      }
	      else if((nMB1MatchClusterAdjacent0p8Plus>=8 && nMB1SegMatchClusterAdjacent0p8Plus>0) || (nMB1MatchClusterAdjacent0p8Minus>=8 && nMB1SegMatchClusterAdjacent0p8Minus>0)){
		h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(1.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(1.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(1.0); } 
	      }
	      else if((nMB1MatchClusterAdjacent0p8Plus>=8 && nMB1SegMatchClusterAdjacent0p8Plus==0) || (nMB1MatchClusterAdjacent0p8Minus>=8 && nMB1SegMatchClusterAdjacent0p8Minus==0)){
		h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(2.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(2.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(2.0); } 
	      }
	      else if((nMB1MatchClusterAdjacent0p8Plus<8 && nMB1SegMatchClusterAdjacent0p8Plus>0) || (nMB1MatchClusterAdjacent0p8Minus<8 && nMB1SegMatchClusterAdjacent0p8Minus>0)){
		h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(3.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(3.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(3.0); } 
	      }
	      else{
		h_AdjacentMB1VetoOutcomes_invertedJetVeto[itr_year]->Fill(4.0); 
		if(dtRechitClusterMaxStation[itr_clust]==2){ h_AdjacentMB1VetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==3){ h_AdjacentMB1VetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(4.0); } 
		if(dtRechitClusterMaxStation[itr_clust]==4){ h_AdjacentMB1VetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(4.0); } 
	      }
	      if(dtRechitClusterMaxStation[itr_clust]==2){ 
		h_nMatchedSegmentsClusterStationMB2_invertedJetVeto[itr_year]->Fill(nMB2SegMatchCluster);
		h_nMatchedSegmentsOtherStations_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster+nMB3SegMatchCluster+nMB4SegMatchCluster); 
		h_nMatchedSegmentsOuterStation_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsOuterStationMB2_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsInnerStationMB2_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster);
		h_nMatchedSegmentsInnerAndAdjacentStationMB2_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		h_matchedSegmentTimeMeanClusterStationMB2_invertedJetVeto[itr_year]->Fill(meanMB2SegTime);
		h_matchedSegmentTimeMeanInnerStationMB2_invertedJetVeto[itr_year]->Fill(meanMB1SegTime);
		h_matchedSegmentTimeMeanOuterStationMB2_invertedJetVeto[itr_year]->Fill(meanMB3SegTime);
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		  h_nMatchedSegmentsInnerStationPassRecoMuonMB2_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster); 
		  h_muonVetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(0.0);
		}
		if(passMuonLoose){ 
		  h_nMatchedSegmentsInnerStationPassLooseMuonMB2_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster); 
		  h_muonVetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(1.0);
		}
		if(nMB1SegMatchCluster!=1){ 
		  h_nMatchedSegmentsInnerStationPassOneSegMuonMB2_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster); 
		  h_muonVetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(2.0);
		}
		if(nMB1SegMatchCluster==0){
		  h_muonVetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(3.0);
		}
		h_muonVetoOutcomesMB2_invertedJetVeto[itr_year]->Fill(4.0);
	      }
	      if(dtRechitClusterMaxStation[itr_clust]==3){ 
		h_nMatchedSegmentsClusterStationMB3_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsOtherStations_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB4SegMatchCluster); 
		h_nMatchedSegmentsInnerStation_invertedJetVeto[itr_year]->Fill(nMB2SegMatchCluster);
		h_nMatchedSegmentsOuterStation_invertedJetVeto[itr_year]->Fill(nMB4SegMatchCluster);
		h_nMatchedSegmentsInnerStationMB3_invertedJetVeto[itr_year]->Fill(nMB2SegMatchCluster);
		h_nMatchedSegmentsOuterStationMB3_invertedJetVeto[itr_year]->Fill(nMB4SegMatchCluster);
		h_nMatchedSegmentsInnerAndAdjacentStationMB3_invertedJetVeto[itr_year]->Fill(nMB2SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		h_matchedSegmentTimeMeanClusterStationMB3_invertedJetVeto[itr_year]->Fill(meanMB3SegTime);
		h_matchedSegmentTimeMeanInnerStationMB3_invertedJetVeto[itr_year]->Fill(meanMB2SegTime);
		h_matchedSegmentTimeMeanOuterStationMB3_invertedJetVeto[itr_year]->Fill(meanMB4SegTime);
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		  h_nMatchedSegmentsInnerStationPassRecoMuonMB3_invertedJetVeto[itr_year]->Fill(nMB2SegMatchCluster); 
		  h_muonVetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(0.0);
		}
		if(passMuonLoose){ 
		  h_nMatchedSegmentsInnerStationPassLooseMuonMB3_invertedJetVeto[itr_year]->Fill(nMB2SegMatchCluster); 
		  h_muonVetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(1.0);
		}
		if(nMB2SegMatchCluster!=1){ 
		  h_nMatchedSegmentsInnerStationPassOneSegMuonMB3_invertedJetVeto[itr_year]->Fill(nMB2SegMatchCluster); 
		  h_muonVetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(2.0);
		}
		if(nMB2SegMatchCluster==0){
		  h_muonVetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(3.0);
		}
		h_muonVetoOutcomesMB3_invertedJetVeto[itr_year]->Fill(4.0);
	      }
	      if(dtRechitClusterMaxStation[itr_clust]==4){
		h_nMatchedSegmentsClusterStationMB4_invertedJetVeto[itr_year]->Fill(nMB4SegMatchCluster);
		h_nMatchedSegmentsOtherStations_invertedJetVeto[itr_year]->Fill(nMB1SegMatchCluster+nMB2SegMatchCluster+nMB3SegMatchCluster); 
		h_nMatchedSegmentsInnerStation_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsInnerStationMB4_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster);
		h_nMatchedSegmentsInnerAndAdjacentStationMB4_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster+nMB1SegMatchClusterAdjacent0p8Plus+nMB1SegMatchClusterAdjacent0p8Minus);
		if(nMB3SegMatchCluster==0 && nMB2SegMatchCluster==0){ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_year]->Fill(0.0); }
		else if(nMB3SegMatchCluster>0 && nMB2SegMatchCluster==0){ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_year]->Fill(1.0); }
		else if(nMB3SegMatchCluster==0 && nMB2SegMatchCluster>0){ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_year]->Fill(2.0); }
		else if(nMB2SegMatchCluster>0 && nMB3SegMatchCluster>0){ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_year]->Fill(3.0); }
		else{ h_nMatchedSegmentsInnerStationsMB4_invertedJetVeto[itr_year]->Fill(4.0); }
		h_matchedSegmentTimeMeanClusterStationMB4_invertedJetVeto[itr_year]->Fill(meanMB4SegTime);
		h_matchedSegmentTimeMeanInnerStationMB4_invertedJetVeto[itr_year]->Fill(meanMB3SegTime);
		if(dtRechitClusterMuonVetoPt[itr_clust]<10.){ 
		  h_nMatchedSegmentsInnerStationPassRecoMuonMB4_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster); 
		  h_muonVetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(0.0);
		}
		if(passMuonLoose){ 
		  h_nMatchedSegmentsInnerStationPassLooseMuonMB4_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster); 
		  h_muonVetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(1.0);
		}
		if(nMB3SegMatchCluster!=1){ 
		  h_nMatchedSegmentsInnerStationPassOneSegMuonMB4_invertedJetVeto[itr_year]->Fill(nMB3SegMatchCluster); 
		  h_muonVetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(2.0);
		}
		if(nMB3SegMatchCluster==0){
		  h_muonVetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(3.0);
		}
		h_muonVetoOutcomesMB4_invertedJetVeto[itr_year]->Fill(4.0);
	      }
	    }
	    
	    
	    if(passJetMET && passStations25 && passWheels25){
	      nPassCluster_JetMET_StationsWheels_InvertedJet += 1;
	      if(dtRechitClusterMaxStation[itr_clust]>=2){
		nPassMaxStation_InvertedJet += 1;
		if(!passJet){
		  nPassInvertedJetVeto_InvertedJet += 1;
		  //if(dtRechitClusterMuonVetoPt[itr_clust]<10.0){
		  if(1==1){
		    if(!rpcBx.empty()){
		      nPassRpcMatch_InvertedJet += 1;
		      invertedJetMB1=0;
		      for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
			if(dtRechitStation[itr_dt]==1){
			  dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
			  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			  if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtRechitEta[itr_dt],2))<0.4){
			    invertedJetMB1+=1;
			  }
			}
		      }
		      if(fabs(dPhiClusterMET)<1.0){
			nPassdPhiClusterMET_InvertedJet += 1;
			if(dtRechitClusterSize[itr_clust]>=100){
			  nPassClusterSize_InvertedJet += 1;
			  h_matchedJetPt_invertedJetVeto[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust]);
			  if(invertedJetMB1>1){ eventListInvertedJet << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
			  h_nMB1MatchJet_invertedJetVeto[itr_year]->Fill(invertedJetMB1);
			  if(invertedJetMB1<=1){
			    h_clusterPhi_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
			    h_clusterEta_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterEta[itr_clust]);
			    h_clusterEtaPhi_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterEta[itr_clust],dtRechitClusterPhi[itr_clust]);
			    h_jetEtaPhi_passMB1Veto_invertedJetVeto[itr_year]->Fill(matchedJetEta,matchedJetPhi);
			    //eventListInvertedJetPassMB1 << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			  }
			  if(dtRechitClusterMaxStation[itr_clust]==2){ 
			    h_nMB1MatchJet_MB2_invertedJetVeto[itr_year]->Fill(invertedJetMB1); 
			    if(invertedJetMB1<=1){
			      h_clusterPhi_MB2_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
			      h_clusterEta_MB2_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterEta[itr_clust]);
			      h_clusterEtaPhi_MB2_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterEta[itr_clust],dtRechitClusterPhi[itr_clust]);
			    }
			  }
			  if(dtRechitClusterMaxStation[itr_clust]==3){ 
			    h_nMB1MatchJet_MB3_invertedJetVeto[itr_year]->Fill(invertedJetMB1); 
			    if(invertedJetMB1<=1){
			      h_clusterPhi_MB3_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
			      h_clusterEta_MB3_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterEta[itr_clust]);
			      h_clusterEtaPhi_MB3_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterEta[itr_clust],dtRechitClusterPhi[itr_clust]);
			    }
			  }
			  if(dtRechitClusterMaxStation[itr_clust]==4){ 
			    h_nMB1MatchJet_MB4_invertedJetVeto[itr_year]->Fill(invertedJetMB1); 
			    if(invertedJetMB1<=1){
			      h_clusterPhi_MB4_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
			      h_clusterEta_MB4_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterEta[itr_clust]);
			      h_clusterEtaPhi_MB4_passMB1Veto_invertedJetVeto[itr_year]->Fill(dtRechitClusterEta[itr_clust],dtRechitClusterPhi[itr_clust]);
			    }
			  }
			  h_matchedJetPt_nMB1Match_invertedJetVeto[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1);
			}
		      }
		      if(dtRechitClusterSize[itr_clust]>=100){
			h_matchedJetPt_nMB1Match_invertedJetVetoNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1);
			if(dtRechitClusterMaxStation[itr_clust]==2){ h_matchedJetPt_nMB1Match_MB2_invertedJetVetoNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			if(dtRechitClusterMaxStation[itr_clust]==3){ h_matchedJetPt_nMB1Match_MB3_invertedJetVetoNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			if(dtRechitClusterMaxStation[itr_clust]==4){ h_matchedJetPt_nMB1Match_MB4_invertedJetVetoNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
		      }
		    }
		  }
		  if(dtRechitClusterMuonVetoPt[itr_clust]<10.0){
		    nPassMuonVeto_InvertedJet += 1;
		    if(!rpcBx.empty()){
		      if(dtRechitClusterSize[itr_clust]>=100){
			invertedJetMB1=0;
			for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
			  if(dtRechitStation[itr_dt]==1){
			    dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
			    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			    if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtRechitEta[itr_dt],2))<0.4){
			      invertedJetMB1+=1;
			    }
			  }
			}
			h_matchedJetPt_nMB1Match_invertedJetVetoAllMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1);
			if(dtRechitClusterMaxStation[itr_clust]==2){ h_matchedJetPt_nMB1Match_MB2_invertedJetVetoAllMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			if(dtRechitClusterMaxStation[itr_clust]==3){ h_matchedJetPt_nMB1Match_MB3_invertedJetVetoAllMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			if(dtRechitClusterMaxStation[itr_clust]==4){ h_matchedJetPt_nMB1Match_MB4_invertedJetVetoAllMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
		      }
		    }
		  }
		  if(passMuonLoose){
		    nPassLooseMuonVeto_InvertedJet += 1;
		    if(!rpcBx.empty()){
		      nPassRpcMatch_InvertedJetLooseMuon += 1;
		      if(fabs(dPhiClusterMET)<1.0){
			nPassdPhiClusterMET_InvertedJetLooseMuon += 1;
			if(dtRechitClusterSize[itr_clust]>=100){
			  nPassClusterSize_InvertedJetLooseMuon += 1;
			}
		      }
		      if(dtRechitClusterSize[itr_clust]>50){
			invertedJetMB1=0;
			for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
			  if(dtRechitStation[itr_dt]==1){
			    dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
			    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			    if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtRechitEta[itr_dt],2))<0.4){
			      invertedJetMB1+=1;
			    }
			  }
			}
			if(dtRechitClusterSize[itr_clust]<100){
			  if(dtRechitClusterMaxStation[itr_clust]==2){ h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonLowSizeNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			  if(dtRechitClusterMaxStation[itr_clust]==3){ h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonLowSizeNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			  if(dtRechitClusterMaxStation[itr_clust]==4){ h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonLowSizeNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
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
			if(nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8){
			  if(dtRechitClusterMaxStation[itr_clust]==2){ h_matchedJetPt_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMETNoSize[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],dtRechitClusterSize[itr_clust]); }
			  if(dtRechitClusterMaxStation[itr_clust]==3){ h_matchedJetPt_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMETNoSize[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],dtRechitClusterSize[itr_clust]); }
			  if(dtRechitClusterMaxStation[itr_clust]==4){ h_matchedJetPt_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMETNoSize[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],dtRechitClusterSize[itr_clust]); }
			  if(dtRechitClusterSize[itr_clust]>=100){
			    h_matchedJetPt_nMB1Match_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1);
			    h_matchedJetNHF[itr_year]->Fill(matchedJetNHF);
			    h_matchedJetNEMF[itr_year]->Fill(matchedJetNEF);
			    h_matchedJetCHF[itr_year]->Fill(matchedJetCHF);
			    h_matchedJetCEMF[itr_year]->Fill(matchedJetCEF);
			    if(dtRechitClusterMaxStation[itr_clust]==2){ 
			      h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); 
			      h_matchedJetPt_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],dtRechitClusterSize[itr_clust]); 
			      h_MET_clusterSize_MB2_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(*MET,dtRechitClusterSize[itr_clust]); 
			      if(invertedJetMB1>1){ h_matchedJetPt_nMB1Match_MB2_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			    }
			    if(dtRechitClusterMaxStation[itr_clust]==3){ 
			      h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); 
			      h_matchedJetPt_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],dtRechitClusterSize[itr_clust]); 
			      h_MET_clusterSize_MB3_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(*MET,dtRechitClusterSize[itr_clust]); 
			      if(invertedJetMB1>1){ h_matchedJetPt_nMB1Match_MB4_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			    }
			    if(dtRechitClusterMaxStation[itr_clust]==4){ 
			      h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); 
			      h_matchedJetPt_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],dtRechitClusterSize[itr_clust]); 
			      h_MET_clusterSize_MB4_invertedJetVetoLooseMuonNoClusterMET[itr_year]->Fill(*MET,dtRechitClusterSize[itr_clust]); 
			      if(invertedJetMB1>1){ h_matchedJetPt_nMB1Match_MB4_invertedJetVetoInvertedMB1LooseMuonNoClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			    }
			    if(fabs(dPhiClusterMET)<1.0){
			      if(dtRechitClusterMaxStation[itr_clust]==2){ h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonLowClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			      if(dtRechitClusterMaxStation[itr_clust]==3){ h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonLowClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			      if(dtRechitClusterMaxStation[itr_clust]==4){ h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonLowClusterMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			    }
			  }
			  if(dtRechitClusterSize[itr_clust]>=80){
			    if(dtRechitClusterMaxStation[itr_clust]==2){ h_matchedJetPt_nMB1Match_MB2_invertedJetVetoLooseMuonNoClusterMETCluster80[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			    if(dtRechitClusterMaxStation[itr_clust]==3){ h_matchedJetPt_nMB1Match_MB3_invertedJetVetoLooseMuonNoClusterMETCluster80[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			    if(dtRechitClusterMaxStation[itr_clust]==4){ h_matchedJetPt_nMB1Match_MB4_invertedJetVetoLooseMuonNoClusterMETCluster80[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1); }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	    if(passStations25 && passWheels25){
	      if(dtRechitClusterMaxStation[itr_clust]>=2){
		if(!passJet){
		  if(!rpcBx.empty()){	    
		    if(fabs(dPhiClusterMET)<1.0){
		      if(dtRechitClusterSize[itr_clust]>=100){
			invertedJetMB1=0;
			for(Int_t itr_dt=0; itr_dt<*nDtRechits; itr_dt++){
			  if(dtRechitStation[itr_dt]==1){
			    dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtRechitPhi[itr_dt];
			    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
			    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
			    if(sqrt(pow(dPhi_tmp,2)+pow(dtRechitClusterEta[itr_clust] - dtRechitEta[itr_dt],2))<0.4){
			      invertedJetMB1+=1;
			    }
			  }
			}
			h_matchedJetPt_nMB1Match_invertedJetVetoNoJetMET[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust],invertedJetMB1);
		      }
		    }
		  }
		}
	      }
	    }

	    if(fabs(dPhi_min)<0.6 && passJet && passMuon){
	      nClustersVeto_dPhiJetMET+=1;
	      h_dtRechitClusterNSegmentStation2_dPhiJetMET[itr_year]->Fill(dtRechitClusterNSegmentStation2[itr_clust]);
	      h_dtRechitClusterNSegmentStation3_dPhiJetMET[itr_year]->Fill(dtRechitClusterNSegmentStation3[itr_clust]);
	      h_dtRechitClusterNSegmentStation4_dPhiJetMET[itr_year]->Fill(dtRechitClusterNSegmentStation4[itr_clust]);
	      h_dtRechitClusterMaxStation_dPhiJetMET[itr_year]->Fill(dtRechitClusterMaxStation[itr_clust]);
	      h_dtRechitClusterNStation_dPhiJetMET[itr_year]->Fill(dtRechitClusterNStation[itr_clust]);
	      h_dtRechitClusterMaxStationRatio_dPhiJetMET[itr_year]->Fill(dtRechitClusterMaxStationRatio[itr_clust]);
	      h_dtRechitClusterMaxChamberRatio_dPhiJetMET[itr_year]->Fill(dtRechitClusterMaxChamberRatio[itr_clust]);
	      h_dtRechitClusterNChamber_dPhiJetMET[itr_year]->Fill(dtRechitClusterNChamber[itr_clust]);
	      h_dtRechitClusterMaxChamber_dPhiJetMET[itr_year]->Fill(dtRechitClusterMaxChamber[itr_clust]);
	      h_dtRechitClusterX_dPhiJetMET[itr_year]->Fill(dtRechitClusterX[itr_clust]);
	      h_dtRechitClusterY_dPhiJetMET[itr_year]->Fill(dtRechitClusterY[itr_clust]);
	      h_dtRechitClusterZ_dPhiJetMET[itr_year]->Fill(dtRechitClusterZ[itr_clust]);
	      h_dtRechitClusterEta_dPhiJetMET[itr_year]->Fill(dtRechitClusterEta[itr_clust]);
	      h_dtRechitClusterPhi_dPhiJetMET[itr_year]->Fill(dtRechitClusterPhi[itr_clust]);
	      h_dtRechitClusterTime_dPhiJetMET[itr_year]->Fill(dtRechitClusterTime[itr_clust]);
	      h_dtRechitClusterXSpread_dPhiJetMET[itr_year]->Fill(dtRechitClusterXSpread[itr_clust]);
	      h_dtRechitClusterYSpread_dPhiJetMET[itr_year]->Fill(dtRechitClusterYSpread[itr_clust]);
	      h_dtRechitClusterZSpread_dPhiJetMET[itr_year]->Fill(dtRechitClusterZSpread[itr_clust]);
	      h_dtRechitClusterEtaSpread_dPhiJetMET[itr_year]->Fill(dtRechitClusterEtaSpread[itr_clust]);
	      h_dtRechitClusterPhiSpread_dPhiJetMET[itr_year]->Fill(dtRechitClusterPhiSpread[itr_clust]);
	      h_dtRechitClusterTimeSpread_dPhiJetMET[itr_year]->Fill(dtRechitClusterTimeSpread[itr_clust]);
	      h_dtRechitClusterMajorAxis_dPhiJetMET[itr_year]->Fill(dtRechitClusterMajorAxis[itr_clust]);
	      h_dtRechitClusterMinorAxis_dPhiJetMET[itr_year]->Fill(dtRechitClusterMinorAxis[itr_clust]);
	      if(nStations25<3 && nWheels25<3 && dtRechitClusterMaxStation[itr_clust]>2){
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
		h_nRB1Match_dPhiJetMET[itr_year]->Fill(nRB1MatchCluster);
		if(dtRechitClusterMaxChamber[itr_clust]==-2){ h_nMB1MatchAdjacent_dPhiJetMET[itr_year]->Fill(nMB1MatchClusterAdjacentPlus); }
		else if(dtRechitClusterMaxChamber[itr_clust]==2){ h_nMB1MatchAdjacent_dPhiJetMET[itr_year]->Fill(nMB1MatchClusterAdjacentMinus); }
		else{
		  h_nMB1MatchAdjacent_dPhiJetMET[itr_year]->Fill(nMB1MatchClusterAdjacentPlus);
		  h_nMB1MatchAdjacent_dPhiJetMET[itr_year]->Fill(nMB1MatchClusterAdjacentMinus); 
		}
		if(nMB1MatchCluster<2){
		  h_nRB1Match_MB1Veto_dPhiJetMET[itr_year]->Fill(nRB1MatchCluster);
		  if(dtRechitClusterMaxChamber[itr_clust]==-2){ h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_year]->Fill(nMB1MatchClusterAdjacentPlus); }
		  else if(dtRechitClusterMaxChamber[itr_clust]==2){ h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_year]->Fill(nMB1MatchClusterAdjacentMinus); }
		  else{
		    h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_year]->Fill(nMB1MatchClusterAdjacentPlus);
		    h_nMB1MatchAdjacent_MB1Veto_dPhiJetMET[itr_year]->Fill(nMB1MatchClusterAdjacentMinus); 
		  }
		} 
	      }
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

	    if(passMB1 && dtRechitClusterMaxStation[itr_clust]>2 && fabs(dPhi_min)<0.6){
	      if(!passJet){
		h_rpcSpread_invertedJetVeto[itr_year]->Fill(rpcSpread);
		if(passMuon){
		  h_rpcSpread_invertedJetVeto_muonVeto[itr_year]->Fill(rpcSpread);
		}
	      }
	      else{
		if(passMuon){
		  h_rpcSpread_fullVeto[itr_year]->Fill(rpcSpread);
		  if(rpcMedian<0){
		    h_rpcSpread_fullVeto_negBx[itr_year]->Fill(rpcSpread);
		  }
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
	    if(passJet && passMuon && nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8 && !rpcBx.empty() && fabs(dPhiClusterMET)>1.0){ h_nMB1Matched_Nminus1_clusterMETCR[itr_year]->Fill(nMB1MatchCluster); }
	    if(passMB1 && passMuon && nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8 && !rpcBx.empty() && fabs(dPhiClusterMET)>1.0){ h_jetVetoPt_Nminus1_clusterMETCR[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust]); }
	    if(!passMB1 && passMuon && nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8 && !rpcBx.empty() && fabs(dPhiClusterMET)<1.0){ h_jetVetoPt_Nminus1_MB1CR[itr_year]->Fill(dtRechitClusterJetVetoPt[itr_clust]); }
	    if(passJet && passMB1 && nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8 && !rpcBx.empty() && fabs(dPhiClusterMET)>1.0){ 
	      h_muonVetoPt_Nminus1_clusterMETCR[itr_year]->Fill(dtRechitClusterMuonVetoPt[itr_clust]); 
	      h_muonLooseIDVetoPt_Nminus1_clusterMETCR[itr_year]->Fill(matchedLooseIDMuonPt);
	    }
	    if(passJet && !passMB1 && nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8 && !rpcBx.empty() && fabs(dPhiClusterMET)<1.0){ 
	      h_muonVetoPt_Nminus1_MB1CR[itr_year]->Fill(dtRechitClusterMuonVetoPt[itr_clust]); 
	      h_muonLooseIDVetoPt_Nminus1_MB1CR[itr_year]->Fill(matchedLooseIDMuonPt);
	    }
	    if(passJet && passMuon && nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8 && passMB1 && fabs(dPhiClusterMET)>1.0){ h_nRPCMatched_Nminus1_clusterMETCR[itr_year]->Fill(rpcBx.size()); }
	    if(passJet && passMuon && nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8 && !passMB1 && fabs(dPhiClusterMET)<1.0){ h_nRPCMatched_Nminus1_MB1CR[itr_year]->Fill(rpcBx.size()); }
	    if(passJet && passMuon && nMB1MatchClusterAdjacentPlus<8 && nMB1MatchClusterAdjacentMinus<8 && !passMB1 && !rpcBx.empty()){ h_dPhiClusterMET_Nminus1_MB1CR[itr_year]->Fill(fabs(dPhiClusterMET)); }
	    if(passJet && passMuon && passMB1 && fabs(dPhiClusterMET)>1.0){ 
	      if(dtRechitClusterMaxChamber[itr_clust]==-2){ h_nMB1MatchedAdjacent_Nminus1_clusterMETCR[itr_year]->Fill(nMB1MatchClusterAdjacentPlus); }
	      if(dtRechitClusterMaxChamber[itr_clust]==2){ h_nMB1MatchedAdjacent_Nminus1_clusterMETCR[itr_year]->Fill(nMB1MatchClusterAdjacentMinus); }
	      else{
		h_nMB1MatchedAdjacent_Nminus1_clusterMETCR[itr_year]->Fill(nMB1MatchClusterAdjacentPlus);
		h_nMB1MatchedAdjacent_Nminus1_clusterMETCR[itr_year]->Fill(nMB1MatchClusterAdjacentMinus);
	      }
	    }
	    if(passJet && passMuon && !passMB1 && fabs(dPhiClusterMET)<1.0){ 
	      if(dtRechitClusterMaxChamber[itr_clust]==-2){ h_nMB1MatchedAdjacent_Nminus1_MB1CR[itr_year]->Fill(nMB1MatchClusterAdjacentPlus); }
	      if(dtRechitClusterMaxChamber[itr_clust]==2){ h_nMB1MatchedAdjacent_Nminus1_MB1CR[itr_year]->Fill(nMB1MatchClusterAdjacentMinus); }
	      else{
		h_nMB1MatchedAdjacent_Nminus1_MB1CR[itr_year]->Fill(nMB1MatchClusterAdjacentPlus);
		h_nMB1MatchedAdjacent_Nminus1_MB1CR[itr_year]->Fill(nMB1MatchClusterAdjacentMinus);
	      }
	    }
	    
	    if(passJet && passMB1){

	      if(passMuon){ 
		h_dPhiClusterRPC_fullVeto[itr_year]->Fill(dPhiClusterRPC); 
		h_dZClusterRPC_fullVeto[itr_year]->Fill(dZClusterRPC); 
	      }

	      //cout << "doing CRs" << endl;
	      if(fabs(dPhi_min)<0.6 && dtRechitClusterMaxStation[itr_clust]>2){
		if(!rpcBx.empty() && rpcSpread==0){
		  if(fabs(dPhiClusterMET)<1.0){
		    if(passMuon){ 
		      h_dtRechitClusterSize_dPhiClusterMETLow_jetMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		      h_dtRechitClusterSize_rpcMatch_jetMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		    }
		    if(passMuon_alt){ 
		      h_dtRechitClusterSize_dPhiClusterMETLow_jetMETCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		      h_dtRechitClusterSize_rpcMatch_jetMETCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		    }
		  }
		  else{
		    if(passMuon){ h_dtRechitClusterSize_dPhiClusterMETHigh_jetMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		    if(passMuon_alt){ h_dtRechitClusterSize_dPhiClusterMETHigh_jetMETCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  }
		}
		else{
		  if(fabs(dPhiClusterMET)<1.0){
		    if(passMuon){ h_dtRechitClusterSize_rpcNoMatch_jetMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		    if(passMuon_alt){ h_dtRechitClusterSize_rpcNoMatch_jetMETCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  }
		}
	      }


	      if(fabs(dPhiClusterMET)>1.0 && dtRechitClusterMaxStation[itr_clust]>2){
		if(!rpcBx.empty() && rpcSpread==0){
		  if(fabs(dPhi_min)>0.6){
		    if(passMuon){ 
		      h_dtRechitClusterSize_dPhiJetMETHigh_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		      h_dtRechitClusterSize_rpcGood_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		      if(event100_clusterMET<4 && dtRechitClusterSize[itr_clust]>=100 && dtRechitClusterSize[itr_clust]<150){
			//eventListClusterMET << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			event100_clusterMET+=1;
		      }
		      else if(event150_clusterMET<4 && dtRechitClusterSize[itr_clust]>150){
			//eventListClusterMET << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			event150_clusterMET+=1;
		      }
		    }
		    if(passMuon_alt){ 
		      h_dtRechitClusterSize_dPhiJetMETHigh_clusterMETCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		      h_dtRechitClusterSize_rpcMatch_clusterMETCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		    }
		  }
		  else{
		    if(passMuon){ h_dtRechitClusterSize_dPhiJetMETLow_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		    if(passMuon_alt){ h_dtRechitClusterSize_dPhiJetMETLow_clusterMETCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  }
		}
		else{
		  if(fabs(dPhi_min)>0.6){
		    if(passMuon){ h_dtRechitClusterSize_rpcBad_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		    if(passMuon_alt){ h_dtRechitClusterSize_rpcNoMatch_clusterMETCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  }
		}
	      }


	      if(fabs(dPhiClusterMET)>1.0 && dtRechitClusterMaxStation[itr_clust]>2){
		if(fabs(dPhi_min)>0.6){
		  if(passMuon){
		    if(!rpcBx.empty()){ 
		      h_dtRechitClusterSize_rpcMatch_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		      if(rpcSpread==0){
			h_dtRechitClusterSize_rpcNoSpread_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			if(rpcMedian<0.0){
			  h_dtRechitClusterSize_rpcNegativeNoSpread_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			}
		      }
		      else{
			h_dtRechitClusterSize_rpcSpread_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			if(rpcMedian<0.0){
			  h_dtRechitClusterSize_rpcNegativeSpread_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			}
		      }
		    }
		    else{ h_dtRechitClusterSize_rpcNoMatch_clusterMETCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  }
		}
	      }


	      if(fabs(dPhiClusterMET)>1.0 && passMuon){
		passFullVeto_clusterCR=true;
		h_nRPCMatched_fullVeto_clusterMETCR[itr_year]->Fill(rpcBx.size());
		h_rpcSpread_fullVeto_clusterMETCR[itr_year]->Fill(rpcSpread);
		h_rpcBx_fullVeto_clusterMETCR[itr_year]->Fill(rpcMedian);
		h_dPhiJetMET_fullVeto_clusterMETCR[itr_year]->Fill(fabs(dPhi_min));
		h_dtRechitClusterMaxStation_fullVeto_clusterMETCR[itr_year]->Fill(dtRechitClusterMaxStation[itr_clust]);

		if(dtRechitClusterMaxStation[itr_clust]>2){
		  passMaxStation_clusterCR=true;
		  if(!rpcBx.empty()){
		    passRPCMatch_clusterCR=true;
		    if(0==0){
		      passRPCSpread_clusterCR=true;
		      if(0>=0.){
			passRPCBx_clusterCR=true;
			if(0==0){
			  passLepton_clusterCR=true;
			  if(nStations50<3 && nWheels50<3){
			    pass50Hits_clusterCR=true;
			    if(nStations25<3 && nWheels25<3){
			      pass25Hits_clusterCR=true;
			      h_dtRechitClusterSize_fullSelection_clusterMETCR[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			      nMB1MatchClusterAdjacentPlus = 0;
			      nMB1MatchClusterAdjacentMinus = 0;
			      nMB1MatchPi2 = 0;
			      nMB1MatchPi2AdjacentPlus = 0;
			      nMB1MatchPi2AdjacentMinus = 0;
			      nMB1MatchClusterAdjacent0p8Plus = 0;
			      nMB1MatchClusterAdjacent0p8Minus = 0;
			      nMB1SegMatchClusterAdjacent0p8Plus = 0;
			      nMB1SegMatchClusterAdjacent0p8Minus = 0;
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
			      for(Int_t itr_dt=0; itr_dt<*nDtSeg; itr_dt++){
				if(dtSegStation[itr_dt]==1){
				  dPhi_tmp = dtRechitClusterPhi[itr_clust] - dtSegPhi[itr_dt];
				  if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
				  if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
				  if(fabs(dPhi_tmp)<TMath::Pi()/4.0){
				    if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]+1){ nMB1SegMatchClusterAdjacent0p8Plus+=1; }
				    if(dtRechitWheel[itr_dt]==dtRechitClusterMaxChamber[itr_clust]-1){ nMB1SegMatchClusterAdjacent0p8Minus+=1; }
				  }
				}
			      }
			      if(dtRechitClusterMaxChamber[itr_clust]==-2){ 
				h_nMB1MatchAdjacent_dPhiClusterMET[itr_year]->Fill(nMB1MatchClusterAdjacentPlus); 
				h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_year]->Fill(nMB1MatchPi2AdjacentPlus); 
				h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_year]->Fill(nMB1MatchClusterAdjacent0p8Plus); 
				h_nMB1SegMatchAdjacent0p8_dPhiClusterMET[itr_year]->Fill(nMB1SegMatchClusterAdjacent0p8Plus); 
				h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_year]->Fill(nMB1MatchPi2Adjacent0p8Plus); 
			      }
			      else if(dtRechitClusterMaxChamber[itr_clust]==2){ 
				h_nMB1MatchAdjacent_dPhiClusterMET[itr_year]->Fill(nMB1MatchClusterAdjacentMinus); 
				h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_year]->Fill(nMB1MatchPi2AdjacentMinus); 
				h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_year]->Fill(nMB1MatchClusterAdjacent0p8Minus); 
				h_nMB1SegMatchAdjacent0p8_dPhiClusterMET[itr_year]->Fill(nMB1SegMatchClusterAdjacent0p8Minus); 
				h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_year]->Fill(nMB1MatchPi2Adjacent0p8Minus); 
			      }
			      else{
				h_nMB1MatchAdjacent_dPhiClusterMET[itr_year]->Fill(nMB1MatchClusterAdjacentPlus);
				h_nMB1MatchAdjacent_dPhiClusterMET[itr_year]->Fill(nMB1MatchClusterAdjacentMinus); 
				h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_year]->Fill(nMB1MatchPi2AdjacentPlus);
				h_nMB1MatchAdjacentPi2_dPhiClusterMET[itr_year]->Fill(nMB1MatchPi2AdjacentMinus); 
				h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_year]->Fill(nMB1MatchClusterAdjacent0p8Plus);
				h_nMB1MatchAdjacent0p8_dPhiClusterMET[itr_year]->Fill(nMB1MatchClusterAdjacent0p8Minus); 
				h_nMB1SegMatchAdjacent0p8_dPhiClusterMET[itr_year]->Fill(nMB1SegMatchClusterAdjacent0p8Plus);
				h_nMB1SegMatchAdjacent0p8_dPhiClusterMET[itr_year]->Fill(nMB1SegMatchClusterAdjacent0p8Minus); 
				h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_year]->Fill(nMB1MatchPi2Adjacent0p8Plus);
				h_nMB1MatchAdjacent0p8Pi2_dPhiClusterMET[itr_year]->Fill(nMB1MatchPi2Adjacent0p8Minus); 
			      }
			      if(nMB1MatchClusterAdjacentPlus>8 || nMB1MatchClusterAdjacentMinus>8){
				//eventListAdjacentMB1 << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			      }
			      //cout << "looping match station" << endl;
			      for(int i=0; i<rpcMatchStation.size(); i++){
				if(dtRechitClusterMaxStation[itr_clust]==3){ h_matchedRPCStation_MB3Cluster_clusterMETCR[itr_year]->Fill(rpcMatchStation[i]); }
				else{ h_matchedRPCStation_MB4Cluster_clusterMETCR[itr_year]->Fill(rpcMatchStation[i]); }
			      }
			      //cout << "done" << endl;
			    }
			  }
			}
		      }
		    }
		    else{
		      if(*nLeptons+*nMuons==0 && nStations25<3 && nWheels25<3){
			//cout << "looping match station" << endl;
			for(int i=0; i<rpcMatchStation.size(); i++){
			  if(dtRechitClusterMaxStation[itr_clust]==3){ h_matchedRPCStation_MB3Cluster_spreadRPC_clusterMETCR[itr_year]->Fill(rpcMatchStation[i]); }
			  else{ h_matchedRPCStation_MB4Cluster_spreadRPC_clusterMETCR[itr_year]->Fill(rpcMatchStation[i]); }
			}
			//cout << "done" << endl;
		      }
		    }
		  }
		}

		if(!rpcBx.empty() && rpcSpread==0 && rpcMedian>=0.){ h_dtRechitClusterMaxStation_Nminus1_clusterMETCR[itr_year]->Fill(dtRechitClusterMaxStation[itr_clust]); }
		if(!rpcBx.empty() && rpcSpread==0 && rpcMedian>=0.){ h_dPhiJetMET_Nminus1_clusterMETCR[itr_year]->Fill(fabs(dPhi_min)); }
		if(!rpcBx.empty() && rpcSpread==0 && dtRechitClusterMaxStation[itr_clust]>2){ h_rpcBx_Nminus1_clusterMETCR[itr_year]->Fill(rpcMedian); }
		if(!rpcBx.empty() && rpcMedian>=0. && dtRechitClusterMaxStation[itr_clust]>2){ h_rpcSpread_Nminus1_clusterMETCR[itr_year]->Fill(rpcSpread); }
		if(dtRechitClusterMaxStation[itr_clust]>2){ h_nRPCMatched_Nminus1_clusterMETCR[itr_year]->Fill(rpcBx.size()); }
	      }
  
	      //cout << "doing rpc CRs" << endl;
	      if((rpcBx.empty() || rpcSpread>0) && dtRechitClusterMaxStation[itr_clust]>2){
		if(fabs(dPhi_min)>0.6){
		  if(passMuon){ 
		    h_dtRechitClusterSize_dPhiJetMETHigh_rpcCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		    if(fabs(dPhiClusterMET<1.0)){ 
		      h_dtRechitClusterSize_dPhiClusterMETLow_rpcCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); 
		      if(dtRechitClusterSize[itr_clust]>=100 && dtRechitClusterSize[itr_clust]<150 && event100_rpc<5){
			//eventListRPC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			event100_rpc+=1;
		      }
		      else if(dtRechitClusterSize[itr_clust]>150 && event150_rpc<5){
			//eventListRPC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl;
			event150_rpc+=1;
		      }
		    }
		    else{ h_dtRechitClusterSize_dPhiClusterMETHigh_rpcCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  }
		  if(passMuon_alt){ 
		    h_dtRechitClusterSize_dPhiJetMETHigh_rpcCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
		    if(fabs(dPhiClusterMET<1.0)){ h_dtRechitClusterSize_dPhiClusterMETLow_rpcCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		    else{ h_dtRechitClusterSize_dPhiClusterMETHigh_rpcCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  }
		}
		else{
		  if(passMuon){ h_dtRechitClusterSize_dPhiJetMETLow_rpcCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  if(passMuon_alt){ h_dtRechitClusterSize_dPhiJetMETLow_rpcCRnoLepton[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		}
	      }

	      //cout << "rpc Spread CR" << endl;
	      if(rpcSpread>0 && rpcSpread<99 && dtRechitClusterMaxStation[itr_clust]>2){
		if(fabs(dPhi_min)>0.6){
		  if(passMuon){
		    if(fabs(dPhiClusterMET<1.0)){ h_dtRechitClusterSize_dPhiClusterMETLow_rpcSpreadCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		    else{ h_dtRechitClusterSize_dPhiClusterMETHigh_rpcSpreadCRmuonVeto[itr_year]->Fill(dtRechitClusterSize[itr_clust]); }
		  }
		}
	      }


	      if((rpcBx.empty() || rpcSpread>0 || rpcMedian<0)){
		if(passMuon){
		  if(dtRechitClusterMaxStation[itr_clust]>2){
		    if(fabs(dPhi_min)>0.6){
		      if(*nLeptons==0){
			if(nStations25<3 && nWheels25<3){
			  if(dtRechitClusterSize[itr_clust]<100){ h_dPhiClusterMET_lowClusterSize_fullSelection_rpcCR[itr_year]->Fill(fabs(dPhiClusterMET)); }
			  else{ h_dPhiClusterMET_highClusterSize_fullSelection_rpcCR[itr_year]->Fill(fabs(dPhiClusterMET)); }
			  h_dtRechitClusterSize_dPhiClusterMET_fullSelection_rpcCR[itr_year]->Fill(dtRechitClusterSize[itr_clust],fabs(dPhiClusterMET));
			  if(fabs(dPhiClusterMET)<1.0){
			    passLowClusterMET_rpcCR = true;
			    if(dtRechitClusterSize[itr_clust]>maxGoodClusterSize){
			      maxGoodClusterSize=dtRechitClusterSize[itr_clust];
			    }
			  }
			  else{
			    passHighClusterMET_rpcCR = true;
			    if(dtRechitClusterSize[itr_clust]>maxGoodClusterSize){
			      maxGoodClusterSize=dtRechitClusterSize[itr_clust];
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }

	      if(!rpcBx.empty()){
		if(passMuon){
		  if(dtRechitClusterMaxStation[itr_clust]>2){
		    if(fabs(dPhi_min)>0.6){
		      if(nStations25<3 && nWheels25<3){
			if(fabs(dPhiClusterMET)<1.0){
			  h_dtRechitClusterSize_dPhiClusterMET_fullSelection_rpcCR[itr_year]->Fill(dtRechitClusterSize[itr_clust],fabs(dPhiClusterMET));
			  if(dtRechitClusterSize[itr_clust]<100){
			    passLowClusterMET_LowClusterSize = true;
			    if(dtRechitClusterSize[itr_clust]>maxGoodClusterSize){
			      maxGoodClusterSize=dtRechitClusterSize[itr_clust];
			    }
			  }
			}
			else{
			  h_dtRechitClusterSize_dPhiClusterMET_fullSelection_rpcCR[itr_year]->Fill(dtRechitClusterSize[itr_clust],fabs(dPhiClusterMET));
			  if(dtRechitClusterSize[itr_clust]>maxGoodClusterSize){
			    maxGoodClusterSize=dtRechitClusterSize[itr_clust];
			  }	
			  if(dtRechitClusterSize[itr_clust]<100){
			    passHighClusterMET_LowClusterSize = true;
			  }
			  else{
			    passHighClusterMET_HighClusterSize = true;
			  }
			}
		      }
		    }
		  }
		}
	      }
	      

	      if(rpcBx.empty() && passMuon){
		passFullVeto_rpcCR=true;
		h_dPhiClusterMET_fullVeto_rpcCR[itr_year]->Fill(fabs(dPhiClusterMET));
		h_dPhiJetMET_fullVeto_rpcCR[itr_year]->Fill(fabs(dPhi_min));
		h_dtRechitClusterMaxStation_fullVeto_rpcCR[itr_year]->Fill(dtRechitClusterMaxStation[itr_clust]);
		if(fabs(dPhi_min)>0.6 && passJet && passMB1){
		  if(dtRechitClusterMaxStation[itr_clust]==2){
		    h_dtRechitClusterPhi_fullSelectionMB2_rpcCR[itr_year]->Fill(fabs(dtRechitClusterPhi[itr_clust]));
		    h_metPhi_fullSelectionMB2_rpcCR[itr_year]->Fill(fabs(*METphi));
		  }
		  else{
		    h_dtRechitClusterPhi_fullSelectionMB34_rpcCR[itr_year]->Fill(fabs(dtRechitClusterPhi[itr_clust]));
		    h_metPhi_fullSelectionMB34_rpcCR[itr_year]->Fill(fabs(*METphi));
		  }
		}
		if(dtRechitClusterMaxStation[itr_clust]>2){
		  passMaxStation_rpcCR=true;
		  h_dPhiClusterMET_Nminus1_rpcCR[itr_year]->Fill(fabs(dPhiClusterMET));
		  h_dPhiJetMET_Nminus1_rpcCR[itr_year]->Fill(fabs(dPhi_min));
		  if(fabs(dPhiClusterMET)<1.0){ 
		    passClusterMET_rpcCR=true;
		    if(fabs(dPhi_min)>0.6){
		      passJetMET_rpcCR=true;
		      if(*nLeptons+*nMuons==0){
			passLepton_rpcCR=true;
			if(nStations50<3 && nWheels50<3){
			  pass50Hits_rpcCR=true;
			  if(nStations25<3 && nWheels25<3){
			    pass25Hits_rpcCR=true;
			    h_dtRechitClusterSize_fullSelection_rpcCR[itr_year]->Fill(dtRechitClusterSize[itr_clust]);
			  }
			}
		      }
		    }
		  }
		}
		if(fabs(dPhiClusterMET)<1.0){ 
		  h_dtRechitClusterMaxStation_Nminus1_rpcCR[itr_year]->Fill(dtRechitClusterMaxStation[itr_clust]); 
		}
	      }


	    }
	  }
	}
      }
      if(passFullVeto_clusterCR){ nPassFullVeto_clusterCR+=1; }
      if(passRPCMatch_clusterCR){ nPassRPCMatch_clusterCR+=1; }
      if(passRPCSpread_clusterCR){ nPassRPCSpread_clusterCR+=1; }
      if(passRPCBx_clusterCR){ nPassRPCBx_clusterCR+=1; }
      if(passMaxStation_clusterCR){ 
	nPassMaxStation_clusterCR+=1; 
	if(maxClusterSize>150){
	  h_nDtSegsStation_150hits_clusterMETCR[itr_year]->Fill(segStation1);
	  h_nDtSegsStation_150hits_clusterMETCR[itr_year]->Fill(segStation2);
	  h_nDtSegsStation_150hits_clusterMETCR[itr_year]->Fill(segStation3);
	  h_nDtSegsStation_150hits_clusterMETCR[itr_year]->Fill(segStation4);

	  h_nDtSegsWheel_150hits_clusterMETCR[itr_year]->Fill(segWheel2);
	  h_nDtSegsWheel_150hits_clusterMETCR[itr_year]->Fill(segWheel1);
	  h_nDtSegsWheel_150hits_clusterMETCR[itr_year]->Fill(segWheel0);
	  h_nDtSegsWheel_150hits_clusterMETCR[itr_year]->Fill(segWheelm1);
	  h_nDtSegsWheel_150hits_clusterMETCR[itr_year]->Fill(segWheelm2);

	  for(int i=0; i<4; i++){
	    for(int j=0; j<5; j++){
	      h_nDtSegsChamber_150hits_clusterMETCR[itr_year]->Fill(segChambers[i][j]);
	    }
	  }

	  h_nDtSegs_150hits_clusterMETCR[itr_year]->Fill(*nDtSeg);
	  
	  h_nStations1_150hits_clusterMETCR[itr_year]->Fill(nStations1);
	  h_nStations25_150hits_clusterMETCR[itr_year]->Fill(nStations25);
	  h_nStations50_150hits_clusterMETCR[itr_year]->Fill(nStations50);
	  h_nWheels1_150hits_clusterMETCR[itr_year]->Fill(nWheels1);
	  h_nWheels25_150hits_clusterMETCR[itr_year]->Fill(nWheels25);
	  h_nWheels50_150hits_clusterMETCR[itr_year]->Fill(nWheels50);

	  h_nStations1Seg_150hits_clusterMETCR[itr_year]->Fill(nStations1Seg);
	  h_nStations5Seg_150hits_clusterMETCR[itr_year]->Fill(nStations5Seg);
	  h_nStations10Seg_150hits_clusterMETCR[itr_year]->Fill(nStations10Seg);
	  h_nWheels1Seg_150hits_clusterMETCR[itr_year]->Fill(nWheels1Seg);
	  h_nWheels5Seg_150hits_clusterMETCR[itr_year]->Fill(nWheels5Seg);
	  h_nWheels10Seg_150hits_clusterMETCR[itr_year]->Fill(nWheels10Seg);

	  h_nRPCStations1_150hits_clusterMETCR[itr_year]->Fill(nRPCStations1);
	  h_nRPCStations5_150hits_clusterMETCR[itr_year]->Fill(nRPCStations5);
	  h_nRPCStations10_150hits_clusterMETCR[itr_year]->Fill(nRPCStations10);
	  h_nRPCWheels1_150hits_clusterMETCR[itr_year]->Fill(nRPCWheels1);
	  h_nRPCWheels5_150hits_clusterMETCR[itr_year]->Fill(nRPCWheels5);
	  h_nRPCWheels10_150hits_clusterMETCR[itr_year]->Fill(nRPCWheels10);
	}
	else if(maxClusterSize>=100){
	  h_nDtSegsStation_100hits_clusterMETCR[itr_year]->Fill(segStation1);
	  h_nDtSegsStation_100hits_clusterMETCR[itr_year]->Fill(segStation2);
	  h_nDtSegsStation_100hits_clusterMETCR[itr_year]->Fill(segStation3);
	  h_nDtSegsStation_100hits_clusterMETCR[itr_year]->Fill(segStation4);

	  h_nDtSegsWheel_100hits_clusterMETCR[itr_year]->Fill(segWheel2);
	  h_nDtSegsWheel_100hits_clusterMETCR[itr_year]->Fill(segWheel1);
	  h_nDtSegsWheel_100hits_clusterMETCR[itr_year]->Fill(segWheel0);
	  h_nDtSegsWheel_100hits_clusterMETCR[itr_year]->Fill(segWheelm1);
	  h_nDtSegsWheel_100hits_clusterMETCR[itr_year]->Fill(segWheelm2);

	  for(int i=0; i<4; i++){
	    for(int j=0; j<5; j++){
	      h_nDtSegsChamber_100hits_clusterMETCR[itr_year]->Fill(segChambers[i][j]);
	    }
	  }

	  h_nDtSegs_100hits_clusterMETCR[itr_year]->Fill(*nDtSeg);
	  
	  h_nStations1_100hits_clusterMETCR[itr_year]->Fill(nStations1);
	  h_nStations25_100hits_clusterMETCR[itr_year]->Fill(nStations25);
	  h_nStations50_100hits_clusterMETCR[itr_year]->Fill(nStations50);
	  h_nWheels1_100hits_clusterMETCR[itr_year]->Fill(nWheels1);
	  h_nWheels25_100hits_clusterMETCR[itr_year]->Fill(nWheels25);
	  h_nWheels50_100hits_clusterMETCR[itr_year]->Fill(nWheels50);

	  h_nStations1Seg_100hits_clusterMETCR[itr_year]->Fill(nStations1Seg);
	  h_nStations5Seg_100hits_clusterMETCR[itr_year]->Fill(nStations5Seg);
	  h_nStations10Seg_100hits_clusterMETCR[itr_year]->Fill(nStations10Seg);
	  h_nWheels1Seg_100hits_clusterMETCR[itr_year]->Fill(nWheels1Seg);
	  h_nWheels5Seg_100hits_clusterMETCR[itr_year]->Fill(nWheels5Seg);
	  h_nWheels10Seg_100hits_clusterMETCR[itr_year]->Fill(nWheels10Seg);

	  h_nRPCStations1_100hits_clusterMETCR[itr_year]->Fill(nRPCStations1);
	  h_nRPCStations5_100hits_clusterMETCR[itr_year]->Fill(nRPCStations5);
	  h_nRPCStations10_100hits_clusterMETCR[itr_year]->Fill(nRPCStations10);
	  h_nRPCWheels1_100hits_clusterMETCR[itr_year]->Fill(nRPCWheels1);
	  h_nRPCWheels5_100hits_clusterMETCR[itr_year]->Fill(nRPCWheels5);
	  h_nRPCWheels10_100hits_clusterMETCR[itr_year]->Fill(nRPCWheels10);
	}
	else{

	  h_nDtSegsStation_50hits_clusterMETCR[itr_year]->Fill(segStation1);
	  h_nDtSegsStation_50hits_clusterMETCR[itr_year]->Fill(segStation2);
	  h_nDtSegsStation_50hits_clusterMETCR[itr_year]->Fill(segStation3);
	  h_nDtSegsStation_50hits_clusterMETCR[itr_year]->Fill(segStation4);

	  h_nDtSegsWheel_50hits_clusterMETCR[itr_year]->Fill(segWheel2);
	  h_nDtSegsWheel_50hits_clusterMETCR[itr_year]->Fill(segWheel1);
	  h_nDtSegsWheel_50hits_clusterMETCR[itr_year]->Fill(segWheel0);
	  h_nDtSegsWheel_50hits_clusterMETCR[itr_year]->Fill(segWheelm1);
	  h_nDtSegsWheel_50hits_clusterMETCR[itr_year]->Fill(segWheelm2);

	  for(int i=0; i<4; i++){
	    for(int j=0; j<5; j++){
	      h_nDtSegsChamber_50hits_clusterMETCR[itr_year]->Fill(segChambers[i][j]);
	    }
	  }

	  h_nDtSegs_50hits_clusterMETCR[itr_year]->Fill(*nDtSeg);
	  

	  h_nStations1_50hits_clusterMETCR[itr_year]->Fill(nStations1);
	  h_nStations25_50hits_clusterMETCR[itr_year]->Fill(nStations25);
	  h_nStations50_50hits_clusterMETCR[itr_year]->Fill(nStations50);
	  h_nWheels1_50hits_clusterMETCR[itr_year]->Fill(nWheels1);
	  h_nWheels25_50hits_clusterMETCR[itr_year]->Fill(nWheels25);
	  h_nWheels50_50hits_clusterMETCR[itr_year]->Fill(nWheels50);

	  h_nStations1Seg_50hits_clusterMETCR[itr_year]->Fill(nStations1Seg);
	  h_nStations5Seg_50hits_clusterMETCR[itr_year]->Fill(nStations5Seg);
	  h_nStations10Seg_50hits_clusterMETCR[itr_year]->Fill(nStations10Seg);
	  h_nWheels1Seg_50hits_clusterMETCR[itr_year]->Fill(nWheels1Seg);
	  h_nWheels5Seg_50hits_clusterMETCR[itr_year]->Fill(nWheels5Seg);
	  h_nWheels10Seg_50hits_clusterMETCR[itr_year]->Fill(nWheels10Seg);

	  h_nRPCStations1_50hits_clusterMETCR[itr_year]->Fill(nRPCStations1);
	  h_nRPCStations5_50hits_clusterMETCR[itr_year]->Fill(nRPCStations5);
	  h_nRPCStations10_50hits_clusterMETCR[itr_year]->Fill(nRPCStations10);
	  h_nRPCWheels1_50hits_clusterMETCR[itr_year]->Fill(nRPCWheels1);
	  h_nRPCWheels5_50hits_clusterMETCR[itr_year]->Fill(nRPCWheels5);
	  h_nRPCWheels10_50hits_clusterMETCR[itr_year]->Fill(nRPCWheels10);
	}
      }
      if(passLepton_clusterCR){ nPassLepton_clusterCR+=1; }
      if(pass50Hits_clusterCR){ nPass50Hits_clusterCR+=1; }
      if(pass25Hits_clusterCR){ nPass25Hits_clusterCR+=1; }

      if(passRPCCR){ nPassRPCCR+=1; }
      if(passFullVeto_rpcCR){ nPassFullVeto_rpcCR+=1; }
      if(passClusterMET_rpcCR){ nPassClusterMET_rpcCR+=1; }
      if(passMaxStation_rpcCR){ nPassMaxStation_rpcCR+=1; }
      if(passJetMET_rpcCR){ nPassJetMET_rpcCR+=1; }
      if(passLepton_rpcCR){ nPassLepton_rpcCR+=1; }
      if(pass50Hits_rpcCR){ nPass50Hits_rpcCR+=1; }
      if(pass25Hits_rpcCR){ nPass25Hits_rpcCR+=1; }

      
      if(passLowClusterMET_rpcCR){ h_dtRechitClusterSize_lowClusterMET_fullSelection_rpcCR[itr_year]->Fill(maxGoodClusterSize); }
      if(passHighClusterMET_rpcCR){ h_dtRechitClusterSize_highClusterMET_fullSelection_rpcCR[itr_year]->Fill(maxGoodClusterSize); }
      if(passLowClusterMET_LowClusterSize){ h_dtRechitClusterSize_lowClusterMET_lowClusterSize_SR[itr_year]->Fill(maxGoodClusterSize); }
      if(passHighClusterMET_LowClusterSize){ h_dtRechitClusterSize_highClusterMET_lowClusterSize_SR[itr_year]->Fill(maxGoodClusterSize); }
      if(passHighClusterMET_HighClusterSize){ h_dtRechitClusterSize_highClusterMET_highClusterSize_SR[itr_year]->Fill(maxGoodClusterSize); }


      h_nDtRechitClustersVeto_dPhiJetMET[itr_year]->Fill(nClustersVeto_dPhiJetMET);
      if(clusterPhi.size()>1){
	for(Int_t i=0; i<clusterPhi.size()-1; i++){
	  for(Int_t j=i+1; j<clusterPhi.size(); j++){
	    dPhi_tmp = clusterPhi[i] - clusterPhi[j];
	    if(dPhi_tmp > TMath::Pi()){ dPhi_tmp -= 2*TMath::Pi(); }
	    if(dPhi_tmp < -1.0*TMath::Pi()){ dPhi_tmp += 2*TMath::Pi(); }
	    h_dtRechitClustersVetoDR_dPhiJetMET[itr_year]->Fill(sqrt(pow(dPhi_tmp,2)+pow(clusterEta[i]-clusterEta[j],2)));
	    if(sqrt(pow(dPhi_tmp,2)+pow(clusterEta[i]-clusterEta[j],2))<0.6){
	      clusterSizeTotal+=clusterSize[i];
	      clusterSizeTotal+=clusterSize[j];
	      break;
	    }
	  }
	}
      }

      if(passMET){ nPassMET+=1; }
      if(passOneJet){ nPassOneJet+=1; }
      if(passJetMET){ nPassJetMET+=1; }
      if(passStations25){ nPassStations25+=1; }
      if(passWheels25){ nPassWheels25+=1; }
      if(passNoVetoCluster){ nPassNoVetoCluster+=1; }
      if(passMB1CR){ nPassMB1CR+=1; }
      if(passJetVeto){ nPassJetVeto+=1; }
      if(passMuonVeto){ nPassMuonVeto+=1; }
      if(passMuonVetoLoose){ nPassMuonVetoLoose+=1; }
      if(passRpcMatchMB1HitsCR){
	nPassRpcMatch+=1;
	h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1HitsCR[itr_year]->Fill(clusterSizeMB1HitsCR,fabs(dPhiClusterMETMB1HitsCR));
	if(fabs(dPhiClusterMETMB1HitsCR)<1.0){
	  h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1HitsCR[itr_year]->Fill(clusterSizeMB1HitsCR);
	  h_dtRechitClusterSizeTotal_lowdPhiClusterMET_fullSelection_MB1HitsCR[itr_year]->Fill(max(clusterSizeMB1HitsCR,clusterSizeTotal));
	  nPassLowClusterMETMB1HitsCR+=1;
	  if(clusterSizeMB1HitsCR>=100){ 
	    nPassLowClusterMETHighClusterSizeMB1HitsCR+=1; 
	    if(fabs(dPhiClusterMETMB1HitsCR)<0.5){ h_nMB1Match_veryLowdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(nMB1MatchClusterMB1HitsCR); }
	    else{ h_nMB1Match_lowdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(nMB1MatchClusterMB1HitsCR); }
	  }
	  else{ 
	    nPassLowClusterMETLowClusterSizeMB1HitsCR+=1; 
	    if(fabs(dPhiClusterMETMB1HitsCR)<0.5){ h_nMB1Match_veryLowdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(nMB1MatchClusterMB1HitsCR); }
	    else{ h_nMB1Match_lowdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(nMB1MatchClusterMB1HitsCR); }
	  }
	}
	else{
	  nPassHighClusterMETMB1HitsCR+=1;
	  h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1HitsCR[itr_year]->Fill(clusterSizeMB1HitsCR);
	  if(clusterSizeMB1HitsCR>=100){ 
	    nPassHighClusterMETHighClusterSizeMB1HitsCR+=1; 
	    h_nMB1Match_highdPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(nMB1MatchClusterMB1HitsCR);
	  }
	  else{ 
	    nPassHighClusterMETLowClusterSizeMB1HitsCR+=1; 
	    h_nMB1Match_highdPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(nMB1MatchClusterMB1HitsCR);
	  }
	}
	if(clusterSizeMB1HitsCR>=100){ 
	  h_dPhiClusterMET_highClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(fabs(dPhiClusterMETMB1HitsCR)); 
	  h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(fabs(dPhiClusterMETMB1HitsCR),nMB1MatchClusterMB1HitsCR);
	}
	else{ 
	  h_dPhiClusterMET_lowClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(fabs(dPhiClusterMETMB1HitsCR)); 
	  h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1HitsCR[itr_year]->Fill(fabs(dPhiClusterMETMB1HitsCR),nMB1MatchClusterMB1HitsCR);
	}
	//if(clusterSizeMB1HitsCR>=250){ eventListMB1 << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	if(passRpcMatchMB1HitsCRwithAdjacent){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB1HitsCR[itr_year]->Fill(clusterSizeMB1HitsCRwithAdjacent,fabs(dPhiClusterMETMB1HitsCRwithAdjacent));
	}
	if(passRpcMatchMB1HitsCRwithNHF50){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB1HitsCR[itr_year]->Fill(clusterSizeMB1HitsCRwithNHF50,fabs(dPhiClusterMETMB1HitsCRwithNHF50));
	}
	if(passRpcMatchMB1HitsCRinvertedNHF50){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB1HitsCR[itr_year]->Fill(clusterSizeMB1HitsCRinvertedNHF50,fabs(dPhiClusterMETMB1HitsCRinvertedNHF50));
	}
	if(passRpcMatchMB1HitsCRwithNHFLead){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB1HitsCR[itr_year]->Fill(clusterSizeMB1HitsCRwithNHFLead,fabs(dPhiClusterMETMB1HitsCRwithNHFLead));
	}
	if(passRpcMatchMB1HitsCRinvertedNHFLead){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB1HitsCR[itr_year]->Fill(clusterSizeMB1HitsCRinvertedNHFLead,fabs(dPhiClusterMETMB1HitsCRinvertedNHFLead));
	  //if(clusterSizeMB1HitsCRinvertedNHFLead>=100 && fabs(dPhiClusterMETMB1HitsCRinvertedNHFLead)<1.0){ eventListInvertedJetIDSR << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	  //else{ eventListInvertedJetIDABC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	}
      }
      if(passRpcMatch){
	h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1CR[itr_year]->Fill(clusterSizeMB1CR,fabs(dPhiClusterMETMB1CR));
	if(fabs(dPhiClusterMETMB1CR)<1.0){ 
	  h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1CR[itr_year]->Fill(clusterSizeMB1CR); 
	  if(fabs(dPhiClusterMETMB1HitsCR)<0.5){ 
	    if(clusterSizeMB1CR>=100){ h_nMB1Match_veryLowdPhiClusterMET_highClusterSize_fullSelection_MB1CR[itr_year]->Fill(nMB1MatchClusterMB1CR); }
	    else{ h_nMB1Match_veryLowdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[itr_year]->Fill(nMB1MatchClusterMB1CR); }
	  }
	  else{ 
	    if(clusterSizeMB1CR>=100){ h_nMB1Match_lowdPhiClusterMET_highClusterSize_fullSelection_MB1CR[itr_year]->Fill(nMB1MatchClusterMB1CR); }
	    else{ h_nMB1Match_lowdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[itr_year]->Fill(nMB1MatchClusterMB1CR); }
	  }
	}
	else{ 
	  h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1CR[itr_year]->Fill(clusterSizeMB1CR); 
	  if(clusterSizeMB1CR>=100){ h_nMB1Match_highdPhiClusterMET_highClusterSize_fullSelection_MB1CR[itr_year]->Fill(nMB1MatchClusterMB1CR); }
	  else{ h_nMB1Match_highdPhiClusterMET_lowClusterSize_fullSelection_MB1CR[itr_year]->Fill(nMB1MatchClusterMB1CR); }
	}
      	if(clusterSizeMB1CR>=100){ h_dPhiClusterMET_highClusterSize_fullSelection_MB1CR[itr_year]->Fill(fabs(dPhiClusterMETMB1CR)); }
	else{ h_dPhiClusterMET_lowClusterSize_fullSelection_MB1CR[itr_year]->Fill(fabs(dPhiClusterMETMB1CR)); }
	if(passRpcMatchMB1CRwithNHFLead){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB1CR[itr_year]->Fill(clusterSizeMB1CRwithNHFLead,fabs(dPhiClusterMETMB1CRwithNHFLead));
	}
      }
      if(passRpcMatchMB1HitsNoMB1ClusterCR){
	if(clusterSizeMB1HitsNoMB1ClusterCR>50){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[itr_year]->Fill(clusterSizeMB1HitsNoMB1ClusterCR,fabs(dPhiClusterMETMB1HitsNoMB1ClusterCR));
	  if(fabs(dPhiClusterMETMB1HitsNoMB1ClusterCR)<1.0){ h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[itr_year]->Fill(clusterSizeMB1HitsNoMB1ClusterCR); }
	  else{ h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1HitsNoMB1ClusterCR[itr_year]->Fill(clusterSizeMB1HitsNoMB1ClusterCR); }
	  if(clusterSizeMB1HitsNoMB1ClusterCR>=100){ 
	    h_dPhiClusterMET_highClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[itr_year]->Fill(fabs(dPhiClusterMETMB1HitsNoMB1ClusterCR)); 
	    h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[itr_year]->Fill(fabs(dPhiClusterMETMB1HitsNoMB1ClusterCR),nMB1MatchClusterMB1HitsNoMB1ClusterCR);
	  }
	  else{ 
	    h_dPhiClusterMET_lowClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[itr_year]->Fill(fabs(dPhiClusterMETMB1HitsNoMB1ClusterCR)); 
	    h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1HitsNoMB1ClusterCR[itr_year]->Fill(fabs(dPhiClusterMETMB1HitsNoMB1ClusterCR),nMB1MatchClusterMB1HitsNoMB1ClusterCR);
	  }
	}
      }
      if(passRpcMatchMB1Hits30CR){
	h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30CR[itr_year]->Fill(clusterSizeMB1Hits30CR,fabs(dPhiClusterMETMB1Hits30CR));
	if(fabs(dPhiClusterMETMB1Hits30CR)<1.0){ h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1Hits30CR[itr_year]->Fill(clusterSizeMB1Hits30CR); }
	else{ h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1Hits30CR[itr_year]->Fill(clusterSizeMB1Hits30CR); }
	if(clusterSizeMB1Hits30CR>=100){ 
	  h_dPhiClusterMET_highClusterSize_fullSelection_MB1Hits30CR[itr_year]->Fill(fabs(dPhiClusterMETMB1Hits30CR)); 
	  h_dPhiClusterMET_nMB1Match_highClusterSize_fullSelection_MB1Hits30CR[itr_year]->Fill(fabs(dPhiClusterMETMB1Hits30CR),nMB1MatchClusterMB1Hits30CR);
	}
	else{ 
	  h_dPhiClusterMET_lowClusterSize_fullSelection_MB1Hits30CR[itr_year]->Fill(fabs(dPhiClusterMETMB1Hits30CR)); 
	  h_dPhiClusterMET_nMB1Match_lowClusterSize_fullSelection_MB1Hits30CR[itr_year]->Fill(fabs(dPhiClusterMETMB1Hits30CR),nMB1MatchClusterMB1Hits30CR);
	}
	if(passRpcMatchMB1Hits30MaxMB2CR){ 
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB2CR[itr_year]->Fill(clusterSizeMB1Hits30MaxMB2CR,fabs(dPhiClusterMETMB1Hits30MaxMB2CR));
	}
	if(passRpcMatchMB1Hits30MaxMB3CR){ 
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB3CR[itr_year]->Fill(clusterSizeMB1Hits30MaxMB3CR,fabs(dPhiClusterMETMB1Hits30MaxMB3CR));
	}
	if(passRpcMatchMB1Hits30MaxMB4CR){ 
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1Hits30MaxMB4CR[itr_year]->Fill(clusterSizeMB1Hits30MaxMB4CR,fabs(dPhiClusterMETMB1Hits30MaxMB4CR));
	}
      }
      if(passRpcMatchMB1or2CR){
	h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB1or2CR[itr_year]->Fill(clusterSizeMB1or2CR,fabs(dPhiClusterMETMB1or2CR));
	if(fabs(dPhiClusterMETMB1or2CR)<1.0){ h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB1or2CR[itr_year]->Fill(clusterSizeMB1or2CR); }
	else{ h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB1or2CR[itr_year]->Fill(clusterSizeMB1or2CR); }
	if(clusterSizeMB1or2CR>=100){ h_dPhiClusterMET_highClusterSize_fullSelection_MB1or2CR[itr_year]->Fill(fabs(dPhiClusterMETMB1or2CR)); }
	else{ h_dPhiClusterMET_lowClusterSize_fullSelection_MB1or2CR[itr_year]->Fill(fabs(dPhiClusterMETMB1or2CR)); }
      }
      if(passRpcMatchMB2CR){
	h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB2CR[itr_year]->Fill(clusterSizeMB2CR,fabs(dPhiClusterMETMB2CR));
	if(fabs(dPhiClusterMETMB2CR)<1.0){ h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB2CR[itr_year]->Fill(clusterSizeMB2CR); }
	else{ h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB2CR[itr_year]->Fill(clusterSizeMB2CR); }
	if(clusterSizeMB2CR>=100){ h_dPhiClusterMET_highClusterSize_fullSelection_MB2CR[itr_year]->Fill(fabs(dPhiClusterMETMB2CR)); }
	else{ h_dPhiClusterMET_lowClusterSize_fullSelection_MB2CR[itr_year]->Fill(fabs(dPhiClusterMETMB2CR)); }
	if(passRpcMatchMB2CRwithAdjacent){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB2CR[itr_year]->Fill(clusterSizeMB2CRwithAdjacent,fabs(dPhiClusterMETMB2CRwithAdjacent));
	}
	if(passRpcMatchMB2CRwithAdjacent0p8){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacent0p8MB1Cut_MB2CR[itr_year]->Fill(clusterSizeMB2CRwithAdjacent0p8,fabs(dPhiClusterMETMB2CRwithAdjacent0p8));
	}
	if(passRpcMatchMB2CRwithOther){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_MB2CR[itr_year]->Fill(clusterSizeMB2CRwithOther,fabs(dPhiClusterMETMB2CRwithOther));
	}
	if(passRpcMatchMB2CRwithNHF50){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB2CR[itr_year]->Fill(clusterSizeMB2CRwithNHF50,fabs(dPhiClusterMETMB2CRwithNHF50));
	}
	if(passRpcMatchMB2CRinvertedNHF50){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB2CR[itr_year]->Fill(clusterSizeMB2CRinvertedNHF50,fabs(dPhiClusterMETMB2CRinvertedNHF50));
	}
	if(passRpcMatchMB2CRwithNHFLead){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB2CR[itr_year]->Fill(clusterSizeMB2CRwithNHFLead,fabs(dPhiClusterMETMB2CRwithNHFLead));
	}
	if(passRpcMatchMB2CRinvertedNHFLead){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB2CR[itr_year]->Fill(clusterSizeMB2CRinvertedNHFLead,fabs(dPhiClusterMETMB2CRinvertedNHFLead));
	  //if(clusterSizeMB2CRinvertedNHFLead>=100 && fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ eventListInvertedJetIDSR << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	  //else{ eventListInvertedJetIDABC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	  if(!(*Flag2_HBHENoiseFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_MB2CR[itr_year]->Fill(0.0);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR[itr_year]->Fill(0.0); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR[itr_year]->Fill(0.0); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR[itr_year]->Fill(0.0); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR[itr_year]->Fill(0.0); }
	  }
	  if(!(*Flag2_HBHEIsoNoiseFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_MB2CR[itr_year]->Fill(1);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR[itr_year]->Fill(1); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR[itr_year]->Fill(1); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR[itr_year]->Fill(1); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR[itr_year]->Fill(1); }
	  }
	  if(!(*Flag2_BadPFMuonFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_MB2CR[itr_year]->Fill(2);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR[itr_year]->Fill(2); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR[itr_year]->Fill(2); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR[itr_year]->Fill(2); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR[itr_year]->Fill(2); }
	  }
	  if(!(*Flag2_globalSuperTightHalo2016Filter)){
	    h_invertedLeadJetVeto_failedMETFilters_MB2CR[itr_year]->Fill(3);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR[itr_year]->Fill(3); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR[itr_year]->Fill(3); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR[itr_year]->Fill(3); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR[itr_year]->Fill(3); }
	  }
	  if(!(*Flag2_EcalDeadCellTriggerPrimitiveFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_MB2CR[itr_year]->Fill(4);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR[itr_year]->Fill(4); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR[itr_year]->Fill(4); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR[itr_year]->Fill(4); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR[itr_year]->Fill(4); }
	  }
	  if(year=="2016" && !(*Flag2_ecalBadCalibFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_MB2CR[itr_year]->Fill(5);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeMB2CR[itr_year]->Fill(5); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeMB2CR[itr_year]->Fill(5); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETMB2CR[itr_year]->Fill(5); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETMB2CR[itr_year]->Fill(5); }
	  }
	}
      }
      if(noRpcMatchMB2CRwithNHFLead){
	h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCutNoRPC_MB2CR[itr_year]->Fill(clusterSizeMB2CRwithNHFLeadNoRPC,fabs(dPhiClusterMETMB2CRwithNHFLeadNoRPC));
      }
      if(passMB2CRwithInvertedNHFLeadNoRPC){
	h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCutNoRPC_MB2CR[itr_year]->Fill(clusterSizeMB2CRwithInvertedNHFLeadNoRPC,fabs(dPhiClusterMETMB2CRwithInvertedNHFLeadNoRPC));
      }
      if(passRpcMatchMB2withMB1CR){
	h_dtRechitClusterSize_dPhiClusterMET_fullSelection_MB2withMB1CR[itr_year]->Fill(clusterSizeMB2withMB1CR,fabs(dPhiClusterMETMB2withMB1CR));
	if(fabs(dPhiClusterMETMB2withMB1CR)<1.0){ h_dtRechitClusterSize_lowdPhiClusterMET_fullSelection_MB2withMB1CR[itr_year]->Fill(clusterSizeMB2withMB1CR); }
	else{ h_dtRechitClusterSize_highdPhiClusterMET_fullSelection_MB2withMB1CR[itr_year]->Fill(clusterSizeMB2withMB1CR); }
	if(clusterSizeMB2withMB1CR>=100){ h_dPhiClusterMET_highClusterSize_fullSelection_MB2withMB1CR[itr_year]->Fill(fabs(dPhiClusterMETMB2withMB1CR)); }
	else{ h_dPhiClusterMET_lowClusterSize_fullSelection_MB2withMB1CR[itr_year]->Fill(fabs(dPhiClusterMETMB2withMB1CR)); }
	if(passRpcMatchMB2withMB1CRwithAdjacent){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_MB2withMB1CR[itr_year]->Fill(clusterSizeMB2withMB1CRwithAdjacent,fabs(dPhiClusterMETMB2withMB1CRwithAdjacent));
	}
	if(passRpcMatchMB2withMB1CRwithNHF50){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_MB2withMB1CR[itr_year]->Fill(clusterSizeMB2withMB1CRwithNHF50,fabs(dPhiClusterMETMB2withMB1CRwithNHF50));
	}
	if(passRpcMatchMB2withMB1CRinvertedNHF50){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_MB2withMB1CR[itr_year]->Fill(clusterSizeMB2withMB1CRinvertedNHF50,fabs(dPhiClusterMETMB2withMB1CRinvertedNHF50));
	}
	if(passRpcMatchMB2withMB1CRwithNHFLead){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_MB2withMB1CR[itr_year]->Fill(clusterSizeMB2withMB1CRwithNHFLead,fabs(dPhiClusterMETMB2withMB1CRwithNHFLead));
	  if(clusterSizeMB2withMB1CRwithNHFLead>=80 && fabs(dPhiClusterMETMB2withMB1CRwithNHFLead)<1.0){ 
	    if(nMB1MatchClusterMB2withMB1CRwithNHFLead<=5 && nMB1MatchClusterMB2withMB1CRwithNHFLead>=2){
	      invertedMB1_MB2Cluster_5MB1 += 1;
	    }
	  }
	}
	if(passRpcMatchMB2withMB1CRinvertedNHFLead){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_MB2withMB1CR[itr_year]->Fill(clusterSizeMB2withMB1CRinvertedNHFLead,fabs(dPhiClusterMETMB2withMB1CRinvertedNHFLead));
	  //if(clusterSizeMB2withMB1CRinvertedNHFLead>=100 && fabs(dPhiClusterMETMB2withMB1CRinvertedNHFLead)<1.0){ eventListInvertedJetIDSR << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	  //else{ eventListInvertedJetIDABC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	}
      }

      if(passABCD){
	if(clusterSizeSR<100 || fabs(dPhiClusterMETSR)>=1.0){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SR[itr_year]->Fill(clusterSizeSR,fabs(dPhiClusterMETSR));
	  if(fabs(dPhiClusterMETSR)>=1.0 && clusterSizeSR<100){ SRyieldA+=1; }
	  if(fabs(dPhiClusterMETSR)>=1.0 && clusterSizeSR>=100){ SRyieldB+=1; }
	  if(fabs(dPhiClusterMETSR)<1.0 && clusterSizeSR<100){ SRyieldC+=1; }
	}
	if(passABCDwithOther){
	  if(clusterSizeSRwithOther<100 || fabs(dPhiClusterMETSRwithOther)>=1.0){
	    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SR[itr_year]->Fill(clusterSizeSRwithOther,fabs(dPhiClusterMETSRwithOther));
	    if(maxStationSRwithOther==3){ h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SRMB3[itr_year]->Fill(clusterSizeSRwithOther,fabs(dPhiClusterMETSRwithOther)); }
	    else if(maxStationSRwithOther==4){ h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithOtherStationsCut_SRMB4[itr_year]->Fill(clusterSizeSRwithOther,fabs(dPhiClusterMETSRwithOther)); }	  
	  }
	}
	if(passABCDwithAdjacent){
	  if(clusterSizeSRwithAdjacent<100 || fabs(dPhiClusterMETSRwithAdjacent)>=1.0){
	    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SR[itr_year]->Fill(clusterSizeSRwithAdjacent,fabs(dPhiClusterMETSRwithAdjacent));
	    if(maxStationSRwithAdjacent==3){ h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SRMB3[itr_year]->Fill(clusterSizeSRwithAdjacent,fabs(dPhiClusterMETSRwithAdjacent)); }
	    else if(maxStationSRwithAdjacent==4){ h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacentMB1Cut_SRMB4[itr_year]->Fill(clusterSizeSRwithAdjacent,fabs(dPhiClusterMETSRwithAdjacent)); }
	  }
	}
	if(passABCDwithAdjacent0p8){
	  if(clusterSizeSRwithAdjacent0p8<100 || fabs(dPhiClusterMETSRwithAdjacent0p8)>=1.0){
	    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithAdjacent0p8MB1Cut_SR[itr_year]->Fill(clusterSizeSRwithAdjacent0p8,fabs(dPhiClusterMETSRwithAdjacent0p8));
	  }
	}
	if(passABCDwithNHF50){
	  if(clusterSizeSRwithNHF50<100 || fabs(dPhiClusterMETSRwithNHF50)>=1.0){
	    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHF50Cut_SR[itr_year]->Fill(clusterSizeSRwithNHF50,fabs(dPhiClusterMETSRwithNHF50));
	  }
	}
	if(passABCDinvertedNHF50){
	  if(clusterSizeSRinvertedNHF50<100 || fabs(dPhiClusterMETSRinvertedNHF50)>=1.0){
	    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionInvertedNHF50Cut_SR[itr_year]->Fill(clusterSizeSRinvertedNHF50,fabs(dPhiClusterMETSRinvertedNHF50));
	  }
	}
	if(passABCDwithNHFLead){
	  //if(clusterSizeSRwithNHFLead<100 || fabs(dPhiClusterMETSRwithNHFLead)>=1.0){
	    h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut_SR[itr_year]->Fill(clusterSizeSRwithNHFLead,fabs(dPhiClusterMETSRwithNHFLead));
	    if(maxStationSR==3){ h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SRMB3[itr_year]->Fill(clusterSizeSRwithNHFLead,fabs(dPhiClusterMETSRwithNHFLead)); }
	    else if(maxStationSR==4){ h_dtRechitClusterSize_dPhiClusterMET_fullSelection_SRMB4[itr_year]->Fill(clusterSizeSRwithNHFLead,fabs(dPhiClusterMETSRwithNHFLead)); }
	    //}
	}
	if(passABCDinvertedNHFLead){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCut_SR[itr_year]->Fill(clusterSizeSRinvertedNHFLead,fabs(dPhiClusterMETSRinvertedNHFLead));
	  //if(clusterSizeSRinvertedNHFLead>=100 && fabs(dPhiClusterMETSRinvertedNHFLead)<1.0){ eventListInvertedJetIDSR << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	  //else{ eventListInvertedJetIDABC << *runNum << ":" << *lumiSec << ":" << *eventNum << endl; }
	  if(!(*Flag2_HBHENoiseFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_SR[itr_year]->Fill(0.0);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR[itr_year]->Fill(0.0); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR[itr_year]->Fill(0.0); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR[itr_year]->Fill(0.0); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR[itr_year]->Fill(0.0); }
	  }
	  if(!(*Flag2_HBHEIsoNoiseFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_SR[itr_year]->Fill(1);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR[itr_year]->Fill(1); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR[itr_year]->Fill(1); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR[itr_year]->Fill(1); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR[itr_year]->Fill(1); }
	  }
	  if(!(*Flag2_BadPFMuonFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_SR[itr_year]->Fill(2);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR[itr_year]->Fill(2); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR[itr_year]->Fill(2); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR[itr_year]->Fill(2); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR[itr_year]->Fill(2); }
	  }
	  if(!(*Flag2_globalSuperTightHalo2016Filter)){
	    h_invertedLeadJetVeto_failedMETFilters_SR[itr_year]->Fill(3);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR[itr_year]->Fill(3); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR[itr_year]->Fill(3); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR[itr_year]->Fill(3); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR[itr_year]->Fill(3); }
	  }
	  if(!(*Flag2_EcalDeadCellTriggerPrimitiveFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_SR[itr_year]->Fill(4);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR[itr_year]->Fill(4); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR[itr_year]->Fill(4); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR[itr_year]->Fill(4); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR[itr_year]->Fill(4); }
	  }
	  if(year=="2016" && !(*Flag2_ecalBadCalibFilter)){
	    h_invertedLeadJetVeto_failedMETFilters_SR[itr_year]->Fill(5);
	    if(clusterSizeMB2CRinvertedNHFLead>=100){ h_invertedLeadJetVeto_failedMETFilters_highClusterSizeSR[itr_year]->Fill(5); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_lowClusterSizeSR[itr_year]->Fill(5); }
	    if(fabs(dPhiClusterMETMB2CRinvertedNHFLead)<1.0){ h_invertedLeadJetVeto_failedMETFilters_lowdPhiClusterMETSR[itr_year]->Fill(5); }
	    else{ h_invertedLeadJetVeto_failedMETFilters_highdPhiClusterMETSR[itr_year]->Fill(5); }
	  }
	}
      }
      if(passABCDwithNHFLeadNoRPC){
	if(clusterSizeSRwithNHFLeadNoRPC<100 || fabs(dPhiClusterMETSRwithNHFLeadNoRPC)>=1.0){
	  h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCutNoRPC_SR[itr_year]->Fill(clusterSizeSRwithNHFLeadNoRPC,fabs(dPhiClusterMETSRwithNHFLeadNoRPC));
	}
      }
      if(passABCDwithInvertedNHFLeadNoRPC){
	h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithInvertedNHFLeadCutNoRPC_SR[itr_year]->Fill(clusterSizeSRwithInvertedNHFLeadNoRPC,fabs(dPhiClusterMETSRwithInvertedNHFLeadNoRPC));
      }
      
      evtNum+=1;
    }
  }
  
  cout << "nPassMET: " << nPassMET << endl;
  cout << "nPassOneJet: " << nPassOneJet << endl;
  cout << "nPassJetMET: " << nPassJetMET << endl;
  cout << "nPassStations25: " << nPassStations25 << endl;
  cout << "nPassWheels25: " << nPassWheels25 << endl;
  cout << "nPassNoVetoCluster: " << nPassNoVetoCluster << endl;
  cout << "nPassNoVeto: " << nPassNoVeto << endl;
  cout << " " << endl;
  cout << "nPassMB1CR: " << nPassMB1CR << endl;
  cout << "nPassJetVetoMB1CR: " << nPassJetVeto << endl;
  cout << "nPassMuonVetoMB1CR: " << nPassMuonVeto << endl;
  cout << "   nPassMuonVetoLooseMB1CR: " << nPassMuonVetoLoose << endl;
  cout << "nPassRpcMatchMB1HitsCR: " << nPassRpcMatch << endl;
  cout << "nPassLowClusterMETMB1HitsCR: " << nPassLowClusterMETMB1HitsCR << endl;
  cout << "   nPassLowClusterMETHighClusterSizeMB1HitsCR: " << nPassLowClusterMETHighClusterSizeMB1HitsCR << endl;
  cout << "   nPassLowClusterMETLowClusterSizeMB1HitsCR: " << nPassLowClusterMETLowClusterSizeMB1HitsCR << endl;
  cout << "nPassHighClusterMETMB1HitsCR: " << nPassHighClusterMETMB1HitsCR << endl;
  cout << "   nPassHighClusterMETHighClusterSizeMBHits1CR: " << nPassHighClusterMETHighClusterSizeMB1HitsCR << endl;
  cout << "   nPassHighClusterMETLowClusterSizeMB1HitsCR: " << nPassHighClusterMETLowClusterSizeMB1HitsCR << endl;
  cout << " " << endl;
  cout << "nPassCluster_JetMET_StationsWheels_InvertedJet: " << nPassCluster_JetMET_StationsWheels_InvertedJet << endl;
  cout << "nPassMaxStation_InvertedJet: " << nPassMaxStation_InvertedJet << endl;
  cout << "nPassInvertedJetVeto_InvertedJet: " << nPassInvertedJetVeto_InvertedJet << endl;
  cout << "nPassRpcMatch_InvertedJet: " << nPassRpcMatch_InvertedJet << endl;
  cout << "nPassdPhiClusterMET_InvertedJet: " << nPassdPhiClusterMET_InvertedJet << endl;
  cout << "nPassClusterSize_InvertedJet: " << nPassClusterSize_InvertedJet << endl;
  cout << " " << endl;
  cout << "nPassMuonVeto_InvertedJet: " << nPassMuonVeto_InvertedJet << endl;
  cout << "nPassLooseMuonVeto_InvertedJet: " << nPassLooseMuonVeto_InvertedJet << endl;
  cout << "nPassRpcMatch_InvertedJetLooseMuon: " << nPassRpcMatch_InvertedJetLooseMuon << endl;
  cout << "nPassdPhiClusterMET_InvertedJetLooseMuon: " << nPassdPhiClusterMET_InvertedJetLooseMuon << endl;
  cout << "nPassClusterSize_InvertedJetLooseMuon: " << nPassClusterSize_InvertedJetLooseMuon << endl;
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
  cout << " " << endl;
  cout << "A: " << SRyieldA << endl;
  cout << "B: " << SRyieldB << endl;
  cout << "C: " << SRyieldC << endl;
  cout << "D (pred): " << float(SRyieldC)*SRyieldB/SRyieldA << " +/- " << float(SRyieldC)*SRyieldB/SRyieldA*sqrt(1.0/SRyieldC+1.0/SRyieldB+1.0/SRyieldA) << endl;
  cout << " " << endl;
  cout << "CA clusters: " << endl;
  cout << "   MB2: " << SRyieldMB2A_CA << " " << SRyieldMB2B_CA << " " << SRyieldMB2C_CA << " " << float(SRyieldMB2C_CA)*SRyieldMB2B_CA/SRyieldMB2A_CA << " +/- " << float(SRyieldMB2C_CA)*SRyieldMB2B_CA/SRyieldMB2A_CA*sqrt(1.0/SRyieldMB2C_CA+1.0/SRyieldMB2B_CA+1.0/SRyieldMB2A_CA) << endl;
  cout << "   MB3-4: " << SRyieldA_CA << " " << SRyieldB_CA << " " << SRyieldC_CA << " " << float(SRyieldC_CA)*SRyieldB_CA/SRyieldA_CA << " +/- " << float(SRyieldC_CA)*SRyieldB_CA/SRyieldA_CA*sqrt(1.0/SRyieldC_CA+1.0/SRyieldB_CA+1.0/SRyieldA_CA) << endl;
  cout << " " << endl;
  cout << "inverted MB1 2-5 hits: " << invertedMB1_MB2Cluster_5MB1 << endl;

  _ofile->Write();
  _ofile->Close();
  /*eventListRPC.close();
  eventListNoRPC.close();
  eventListMB1.close();
  eventListClusterMET.close();
  */
  eventListInvertedJet.close();
  eventListInvertedMB1_phiSpike.close();
  eventListInvertedMB1_10GeVJetsA.close();
  eventListInvertedMB1_10GeVJetsC.close();
  eventListABCD_10GeVJetsA.close();
  eventListABCD_10GeVJetsC.close();
  /*eventListInvertedJetPassMB1.close();
  eventListSRB.close();
  eventListABC.close();
  eventListMB2CR.close();
  eventListAdjacentMB1.close();
  eventListInvertedJetIDSR.close();
  eventListInvertedJetIDABC.close();
  eventListInvertedMB1LooseMuon.close();
  */
}
