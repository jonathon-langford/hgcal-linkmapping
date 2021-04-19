#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <map>
#include <sys/stat.h>
#include "json.hpp" // Author: Niels Lohmann https://github.com/nlohmann/json
using nlohmann::json;

unsigned tpgRawDataPacketWords(unsigned ntc) {
  unsigned nWords;
  assert(ntc<=48);
  
  if(ntc==0) { // No hits
    nWords=1;    
  } else if(ntc<8) { // Low occupancy
    nWords=(16 + 3 + 13*ntc + 15)/16;
  } else { // High occupancy
    nWords=(16 + 48 + 7*ntc + 15)/16;
  }
  
  return nWords;
}
TH3D * convert_tcs_to_words( TH3D * hist ){
  
  TH3D * words = (TH3D*) hist->Clone( TString(hist->GetName()) + "_words" );
  
  for ( int i = 1; i < words->GetNbinsX()+1; i++){
    for ( int j = 1; j < words->GetNbinsY()+1; j++){
      for ( int k = 1; k < words->GetNbinsZ()+1; k++){
	
	unsigned ntcs = words->GetBinContent(i,j,k);
	unsigned nwords = tpgRawDataPacketWords(ntcs);

	if (nwords < 0.5){
	  std::cout << "nwords = " << nwords << std::endl;
	}
	words->SetBinContent(i,j,k,nwords);
	
      }
    }
  }
  
  return words;
  
}

unsigned uvMapping(unsigned layer, std::pair<int,int> &uv) {
  unsigned sector(0);
  int offset;
  
  if(layer<=28) { // CE-E    
    if(uv.first>0 && uv.second>=0) return sector;
    
    offset=0;
    if(uv.first>=uv.second && uv.second<0) sector=2;
    else sector=1;
    
  } else if((layer%2)==1) { // CE-H Odd
    if(uv.first>=0 && uv.second>=0) return sector;
    
    offset=-1;    
    if(uv.first>uv.second && uv.second<0) sector=2;
    else sector=1;
    
  } else { // CE-H Even
    if(uv.first>=1 && uv.second>=1) return sector;
    
    offset=1;
    if(uv.first>=uv.second && uv.second<1) sector=2;
    else sector=1;
  }
  
  int up,vp;
  
  if(sector==1) {
    up=uv.second-uv.first;
    vp=-uv.first+offset;    
    
  } else {
    up=-uv.second+offset;
    vp=uv.first-uv.second+offset;
  }
  
  uv.first=up;
  uv.second=vp;
  return sector;
}

//Rotate and convert from cell to tile numbering
unsigned etaphiMapping(unsigned layer, std::pair<int,int> &etaphi, std::string configFileVersion = "TpgV7") {
  unsigned sector(0);
  
  if (etaphi.second > 24 && etaphi.second <= 72){
    sector = 0;
  }
  else if (etaphi.second > 72 && etaphi.second <= 120){
    sector = 2;
  }
  else {
    sector = 1;
  }
  
  int ep;
  int pp;

  if (sector==0) {
    pp=etaphi.second-24;
  }
  else if (sector==2) {
    pp=etaphi.second-72;
  }
  else if (sector==1) {
    if (etaphi.second<=24){
      etaphi.second+=144;
    }
    pp=etaphi.second-120;
  }
  
  pp = (pp-1)/4; //Phi index 1-12

  if ( configFileVersion == "V7" ){
    if ( etaphi.first <= 3 ){
      ep = 0;
    }
    else if ( etaphi.first <= 9 ){
      ep = 1;
    }
    else if ( etaphi.first <= 13 ){
      ep = 2;
    }
    else if ( etaphi.first <= 17 ){
      ep = 3;
    }
    else{
      ep = 4;
    }
  }
  else if ( configFileVersion == "TpgV7" ){

    int split = 12;
    if ( layer > 40 ){
      split = 8;
    }
    if ( etaphi.first <= split ){
      ep = 0;
    }
    else{
      ep = 1;
    }   
  }
  else{
    std::cout << "Expected config file version to be either V7 or TpgV7" << std::endl;
  }
  
  etaphi.first=ep;
  etaphi.second=pp;
  
  
  return sector;
}

std::pair<float,float> getEtaPhi(float x, float y, float z){

  float r = std::sqrt( x*x + y*y  );
  float phi = std::atan2(y,x);
  float eta = -std::log(std::tan(std::atan2( r,z)/2));
  return std::make_pair(eta,phi);

}
std::pair<float,float> getROverZPhi(float x, float y, float z, unsigned sector = 0){

  if ( z > 0 ){
    x *= -1;
  }
  
  float r = std::sqrt( x*x + y*y  );
  float phi = std::atan2(y,x);
  
  if (sector == 1)
    if ( phi < M_PI && phi > 0)
      phi = phi-(2*M_PI/3);
    else
      phi = phi+(4*M_PI/3);
  else if (sector == 2)
    phi = phi+(2*M_PI/3);
        
  return std::make_pair(r/z,phi);

}

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

unsigned get45degSector(float x, float y, float z){

  std::pair<float,float> etaphi = getEtaPhi(x, y, z);

  unsigned sector = 0;

  if ( etaphi.second > -M_PI && etaphi.second <=  -M_PI/2 ){
    sector = 0;
  }
  else if  ( etaphi.second > -M_PI/2 && etaphi.second <= 0){
    sector = 1;
  }
  else if  ( etaphi.second > 0 && etaphi.second <= M_PI/2){
    sector = 2;
  }
  else if  ( etaphi.second > M_PI/2 && etaphi.second <= M_PI){
    sector = 3;
  }

  if ( z > 0 ){
    sector += 4;
  }
  
  return sector;
  
}

int main(int argc, char **argv){  
  TH1::SetDefaultSumw2();
  TH1::AddDirectory(kFALSE);

  //Load json config
  std::string configfile = "config/extract_data.json";
  if (argc > 1){
    configfile = argv[1];
  }
  if ( !file_exists(configfile) ){
    std::cout << "Config file " << configfile << " does not exist, exiting" << std::endl;
    exit(1);
  }
  json config;
  std::ifstream inputjson(configfile);
  inputjson >> config;
  
  std::string input_file = config["inputfile"];
  std::string flat_file_silicon = config["flat_file_silicon"];
  std::string flat_file_scintillator = config["flat_file_scintillator"];
  std::string file_ROverZHistograms = config["file_ROverZHistograms"];
  std::string file_nTCsPerEvent = config["file_nTCsPerEvent"];
  std::string average_tcs_sil = config["average_tcs_sil"];
  std::string average_tcs_scin = config["average_tcs_scin"];
  std::string configFileVersion = config["configFileVersion"];
  bool createFlatFile = config["createFlatFile"];
  
  TFile * file = new TFile(TString(input_file),"READ");
  TTree * tree = (TTree*)file->Get("HGCalTriggerNtuple");

  int nPhiBins = 32*3;//5/3 degree bins
  double phiMax = 7.*M_PI/9.;//140 degrees
  double phiMin = -1.*M_PI/9.;//-20 degrees
  
  // Declaration of leaf types
  std::vector<int>     *tc_layer = 0;
  std::vector<int>     *tc_waferu = 0;
  std::vector<int>     *tc_waferv = 0;
  std::vector<int>     *tc_cellu = 0;
  std::vector<int>     *tc_cellv = 0;
  std::vector<int>     *tc_panel_number = 0;
  std::vector<int>     *tc_panel_sector = 0;
  std::vector<int>     *tc_zside = 0;
  std::vector<float>     *tc_x = 0;
  std::vector<float>     *tc_y = 0;
  std::vector<float>     *tc_z = 0;
  // List of branches
  TBranch        *b_tc_layer = 0;   
  TBranch        *b_tc_waferu = 0;   
  TBranch        *b_tc_waferv = 0;
  TBranch        *b_tc_cellu = 0;   
  TBranch        *b_tc_cellv = 0;
  TBranch        *b_tc_panel_number = 0;   
  TBranch        *b_tc_panel_layer = 0;
  TBranch        *b_tc_zside = 0;
  TBranch        *b_tc_x = 0;
  TBranch        *b_tc_y = 0;
  TBranch        *b_tc_z = 0;   
    
  tree->SetBranchAddress("tc_layer" , &tc_layer , &b_tc_layer);
  tree->SetBranchAddress("tc_waferu", &tc_waferu, &b_tc_waferu);
  tree->SetBranchAddress("tc_waferv", &tc_waferv, &b_tc_waferv);
  tree->SetBranchAddress("tc_cellu", &tc_cellu, &b_tc_cellu);
  tree->SetBranchAddress("tc_cellv", &tc_cellv, &b_tc_cellv);
  tree->SetBranchAddress("tc_zside", &tc_zside, &b_tc_zside);

  tree->SetBranchAddress("tc_x", &tc_x, &b_tc_x);
  tree->SetBranchAddress("tc_y", &tc_y, &b_tc_y);
  tree->SetBranchAddress("tc_z", &tc_z, &b_tc_z);
  
  std::vector<TH3D*> per_event_plus(3);
  std::vector<TH3D*> per_event_minus(3);
  
  std::vector<TH3D*> per_event_plus_scin(3);
  std::vector<TH3D*> per_event_minus_scin(3);

  std::vector<int> ntcs_per_event(8);//divide into 45 degree sectors
  
  for (int i = 0;i<per_event_plus.size();i++ ){
    per_event_plus.at(i) = new TH3D(TString("per_event_plus_hist" + std::to_string(i)),"",15,-0.5,14.5,15,-0.5,14.5,52,0.5,52.5) ;
    per_event_minus.at(i) = new TH3D(TString("per_event_minus_hist" + std::to_string(i)),"",15,-0.5,14.5,15,-0.5,14.5,52,0.5,52.5) ;
    per_event_plus_scin.at(i) = new TH3D(TString("per_event_plus_hist_scin" + std::to_string(i)),"",5,-0.5,4.5,12,-0.5,11.5,52,0.5,52.5) ;
    per_event_minus_scin.at(i) = new TH3D(TString("per_event_minus_hist_scin" + std::to_string(i)),"",5,-0.5,4.5,12,-0.5,11.5,52,0.5,52.5) ;
  }
  
    TH3D * out_words = new TH3D("out_words_hist","",15,-0.5,14.5,15,-0.5,14.5,52,0.5,52.5);
    TH3D * out_tcs = new TH3D("out_tcs_hist","",15,-0.5,14.5,15,-0.5,14.5,52,0.5,52.5);
    TH1D * out_ntcs_per_event = new TH1D("out_ntcs_per_event_hist","",10000,-0.5,9999.5);

    TH3D * out_words_scin = new TH3D("out_words_hist_scin","",5,-0.5,4.5,12,-0.5,11.5,52,0.5,52.5);
    TH3D * out_tcs_scin = new TH3D("out_tcs_hist_scin","",5,-0.5,4.5,12,-0.5,11.5,52,0.5,52.5);

    //R/Z Histograms

    TH2D * ROverZ_Inclusive = new TH2D("ROverZ_Inclusive","",42,0.076,0.58,nPhiBins,phiMin,phiMax);
    std::map<std::tuple<int,int,int,int>,TH2D*> ROverZ_per_module;

    //Create one for each module (silicon at first)
    for ( int i = 0; i < 15; i++){//u
      for ( int j = 0; j < 15; j++){//v
	for ( int k = 1; k < 53; k++){//layer

	  if ( k < 28 && k%2 == 0 ) continue;
	  ROverZ_per_module[std::make_tuple(0,i,j,k)] = new TH2D( ("ROverZ_silicon_" + std::to_string(i) + "_" +  std::to_string(j) +"_"+ std::to_string(k)).c_str(),"",42,0.076,0.58,nPhiBins,phiMin,phiMax);

	}
      }
    }
	  
    for ( int i = 0; i < 5; i++){
      for ( int j = 0; j < 12; j++){
	for ( int k = 37; k < 53; k++){

	  ROverZ_per_module[std::make_tuple(1,i,j,k)] = new TH2D( ("ROverZ_scintillator_" + std::to_string(i) + "_" +  std::to_string(j) +"_"+ std::to_string(k)).c_str(),"",42,0.076,0.58,nPhiBins,phiMin,phiMax);

	}
      }
    }


    std::ofstream f_flattree_silicon;
    std::ofstream f_flattree_scintillator;
    if (createFlatFile){
      f_flattree_silicon.open (flat_file_silicon);
      f_flattree_scintillator.open (flat_file_scintillator);
    }


    
    Int_t nentries = (Int_t)tree->GetEntries();
    Int_t nentries_looped = 0;
    Long64_t nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb = tree->GetEntry(jentry);   
      if (jentry % 100 == 0) std::cout << jentry << " / " << nentries << std::endl;;
      //if (jentry > 100 )break;

      //Reset nTCs counter
      for (auto& x: ntcs_per_event) {
	x = 0;
      }

      
      for (int j = 0;j<tc_waferu->size();j++){
	int u = tc_waferu->at(j);
	int v = tc_waferv->at(j);

	//For counting the number of TCs in a 45 degree sector
	unsigned phisector = get45degSector(tc_x->at(j), tc_y->at(j), tc_z->at(j));
	ntcs_per_event.at(phisector)++;

       	//u,v for silicon and eta,phi for scintillator
	std::pair<int,int> coordinates = std::make_pair(u,v);
	unsigned sector = 0;
	
	if ( u > -990 ){//Silicon

	  sector = uvMapping(tc_layer->at(j),coordinates);
	  
	  if ( tc_zside->at(j) > 0 ){
	    per_event_plus.at(sector)->Fill(coordinates.first , coordinates.second, tc_layer->at(j) );
	  }
	  else if ( tc_zside->at(j) < 0 ){
	    per_event_minus.at(sector)->Fill(coordinates.first , coordinates.second, tc_layer->at(j) );
	  }

	}
	else{
	  int eta = tc_cellu->at(j);
	  int phi = tc_cellv->at(j);

	  coordinates = std::make_pair(eta,phi);
	  
	  sector = etaphiMapping(tc_layer->at(j),coordinates,configFileVersion);

	  if ( tc_zside->at(j) > 0 ){
	    per_event_plus_scin.at(sector)->Fill(coordinates.first , coordinates.second, tc_layer->at(j) );
	  }
	  else if ( tc_zside->at(j) < 0 ){
	    per_event_minus_scin.at(sector)->Fill(coordinates.first , coordinates.second, tc_layer->at(j) );
	  }

	}

	//Fill R/Z-Phi Histograms
	std::pair<float,float> roverz_phi = getROverZPhi(tc_x->at(j),tc_y->at(j),tc_z->at(j),sector);
	ROverZ_Inclusive->Fill(std::abs(roverz_phi.first),roverz_phi.second);

        if ( u > -990 ){//Silicon
	  ROverZ_per_module[std::make_tuple(0,coordinates.first,coordinates.second,tc_layer->at(j))]->Fill(std::abs(roverz_phi.first),roverz_phi.second);
	}
	else{
	  ROverZ_per_module[std::make_tuple(1,coordinates.first,coordinates.second,tc_layer->at(j))]->Fill(std::abs(roverz_phi.first),roverz_phi.second);
	}
      	
      }


      
      std::vector<TH3D*> words_plus;
      std::vector<TH3D*> words_minus;
      std::vector<TH3D*> words_plus_scin;
      std::vector<TH3D*> words_minus_scin;
      //Convert number of trigger cells, to words
      for (int i = 0;i<per_event_plus.size();i++ ){
	words_plus.push_back(convert_tcs_to_words(per_event_plus.at(i)));
	words_minus.push_back(convert_tcs_to_words(per_event_minus.at(i)));
	words_plus_scin.push_back(convert_tcs_to_words(per_event_plus_scin.at(i)));
	words_minus_scin.push_back(convert_tcs_to_words(per_event_minus_scin.at(i)));	
      }




      //If creating flat text file (for Paul Dauncey's studies)
      if (createFlatFile){
	for (int z = 0;z<words_plus_scin.size();z++ ){
	  
	  for ( int i = 0; i < 5; i++){
	    for ( int j = 0; j < 12; j++){
	      for ( int k = 37; k < 53; k++){
		f_flattree_scintillator << words_plus_scin.at(z)->GetBinContent(words_plus_scin.at(z)->FindBin(i,j,k)) << " ";
	      }
	    }
	  }
	  
	  f_flattree_scintillator <<  std::endl;
	  
	  for ( int i = 0; i < 5; i++){
	    for ( int j = 0; j < 12; j++){
	      for ( int k = 37; k < 53; k++){
		f_flattree_scintillator << words_minus_scin.at(z)->GetBinContent(words_minus_scin.at(z)->FindBin(i,j,k)) << " ";
	      }
	    }
	  }
	  
	  f_flattree_scintillator <<  std::endl;
	  
	}
	
	
	for (int z = 0;z<words_plus.size();z++ ){
	  
	  for ( int i = 0; i < 15; i++){
	    for ( int j = 0; j < 15; j++){
	      for ( int k = 1; k < 53; k++){
		f_flattree_silicon << words_plus.at(z)->GetBinContent(words_plus.at(z)->FindBin(i,j,k)) << " ";
	      }
	    }
	  }
	  
	  f_flattree_silicon <<  std::endl;
	  
	  for ( int i = 0; i < 15; i++){
	    for ( int j = 0; j < 15; j++){
	      for ( int k = 1; k < 53; k++){
		f_flattree_silicon << words_minus.at(z)->GetBinContent(words_minus.at(z)->FindBin(i,j,k)) << " ";
	      }
	    }
	  }
	  
	  f_flattree_silicon <<  std::endl;
	  
	}
	
      }

            
      //Add plus and minus sides and all rotated histograms together
      for (int i = 0;i<per_event_plus.size();i++ ){
	out_tcs->Add(per_event_plus.at(i));
	out_tcs->Add(per_event_minus.at(i));
	out_words->Add(words_plus.at(i));
	out_words->Add(words_minus.at(i));
	out_tcs_scin->Add(per_event_plus_scin.at(i));
	out_tcs_scin->Add(per_event_minus_scin.at(i));
	out_words_scin->Add(words_plus_scin.at(i));
	out_words_scin->Add(words_minus_scin.at(i));
      }
      for (int i = 0;i<per_event_plus.size();i++ ){
	per_event_plus.at(i)->Reset();
	per_event_minus.at(i)->Reset();
	words_plus.at(i)->Delete();
	words_minus.at(i)->Delete();
	per_event_plus_scin.at(i)->Reset();
	per_event_minus_scin.at(i)->Reset();
	words_plus_scin.at(i)->Delete();
	words_minus_scin.at(i)->Delete();
      }
      for (auto& x: ntcs_per_event) {
	out_ntcs_per_event->Fill(x);
      }

      nentries_looped++;
    }

    if (createFlatFile){
      f_flattree_silicon.close();
      f_flattree_scintillator.close();
    }

    out_tcs->Scale(1./double(nentries_looped*6.));
    out_words->Scale(1./double(nentries_looped*6.));
    out_tcs_scin->Scale(1./double(nentries_looped*6.));
    out_words_scin->Scale(1./double(nentries_looped*6.));


    //Save Eta Histograms	
    TFile * file_out = new TFile(TString(file_ROverZHistograms),"RECREATE");
    file_out->cd();
    ROverZ_Inclusive->Write();

    for (auto& x: ROverZ_per_module) {
      x.second->Write();
    }
    out_ntcs_per_event->Write();
    file_out->Close();

    //Create output csv
    std::ofstream fout;
    fout.open (average_tcs_sil);

    for ( int i = 0; i < 15; i++){
      for ( int j = 0; j < 15; j++){
	for ( int k = 1; k < 53; k++){

	  double ntcs = out_tcs->GetBinContent(out_tcs->FindBin(i,j,k));
	  double words = out_words->GetBinContent(out_words->FindBin(i,j,k));
	  fout << i << "," << j << "," << k << "," << ntcs << "," << words << "\n";
	}
      }
    }

    fout.close();

    fout.open (average_tcs_scin);
    for ( int i = 0; i < 5; i++){
      for ( int j = 0; j < 12; j++){
	for ( int k = 37; k < 53; k++){

	  double ntcs = out_tcs_scin->GetBinContent(out_tcs_scin->FindBin(i,j,k));
	  double words = out_words_scin->GetBinContent(out_words_scin->FindBin(i,j,k));
	  fout << i << "," << j << "," << k << "," << ntcs << "," << words << "\n";
	}
      }
    }

    fout.close();
    
    out_tcs->Delete();
    out_words->Delete();
    file->Close();
    
    return 0;
}
