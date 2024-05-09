#include <iostream>
#include "CaloAna.h"
using namespace std;

  bool emcal    = 1;
  bool clusters = 0;
  bool ihcal    = 0;
  bool ohcal    = 0;
  bool zdc      = 0;
  bool mbd      = 0;
  bool sepd     = 0;

  CaloAna::CaloAna(const std::string& name, const char* outname, bool _isMC, bool _isHI)
: SubsysReco(name)
  , detector("HCALIN")
{
  outfilename = Form("%s",outname);
//  calibfile = new TFile(Form("/sphenix/user/dlis/Projects/centrality/calib/calib_mbd_%d.root",runN),"read");
//  tcalibfile = new TFile(Form("/sphenix/user/dlis/Projects/centrality/calib/t0_calib_mbd_%d.root",runN),"read");
  isMC = _isMC;
  isHI = _isHI;
}

CaloAna::~CaloAna()
{
    //delete hm;
    delete g4hitntuple;
    delete g4cellntuple;
    delete towerntuple;
    delete clusterntuple;
}

int CaloAna::Init(PHCompositeNode*)
{ 
  try {
    InitOutputFile();
    InitTree();
  }  
  catch (const std::exception& e){
    std::cerr << "CaloAna::Init - Exception during init! For " << e.what() << " return -1" << std::endl;
    return -1;
  }

  return 0;
}

void CaloAna::InitOutputFile(){
  std::cout << "output filename : " << outfilename.c_str() << std::endl;
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  if(!outfile || outfile->IsZombie()){
    throw std::runtime_error("Failed open file");
  }
}

void CaloAna::InitTree(){
    towerntuple = new TTree("towerntup", "Towers");
    towerntuple->Branch("vz",&vz);

    if (emcal) {  
        towerntuple->Branch("emcal_totenergy",&m_emcal_totenergy);
        towerntuple->Branch("emcal_energy",&m_emcal_energy);
        towerntuple->Branch("emcal_etabin",&m_emcal_etabin);
        towerntuple->Branch("emcal_phibin",&m_emcal_phibin);
        towerntuple->Branch("emcal_time",&m_emcal_time);
        towerntuple->Branch("emcal_etacorr",&m_emcal_etacorr);
        towerntuple->Branch("emcal_chi2",&m_emcal_chi2);
//        towerntuple->Branch("emcal_pedestal",&m_emcal_pedestal);
    }
    
    if (clusters) {
      towerntuple->Branch("emcal_cluster_E",&m_emcal_cluster_E);
      towerntuple->Branch("emcal_cluster_Ecore",&m_emcal_cluster_Ecore);
      towerntuple->Branch("emcal_cluster_eta",&m_emcal_cluster_eta);
      towerntuple->Branch("emcal_cluster_phi",&m_emcal_cluster_phi);
      towerntuple->Branch("emcal_cluster_ntower",&m_emcal_cluster_ntower);
      towerntuple->Branch("emcal_cluster_chi2",&m_emcal_cluster_chi2);
      towerntuple->Branch("emcal_cluster_prob",&m_emcal_cluster_prob);
      towerntuple->Branch("emcal_cluster_iso_unsubR01",&m_emcal_cluster_iso_unsubR01);
      towerntuple->Branch("emcal_cluster_iso_unsubR02",&m_emcal_cluster_iso_unsubR02);
      towerntuple->Branch("emcal_cluster_iso_unsubR03",&m_emcal_cluster_iso_unsubR03);
      towerntuple->Branch("emcal_cluster_iso_unsubR04",&m_emcal_cluster_iso_unsubR04);
      towerntuple->Branch("emcal_cluster_iso_subR01",&m_emcal_cluster_iso_subR01);
      towerntuple->Branch("emcal_cluster_iso_subR02",&m_emcal_cluster_iso_subR02);
      towerntuple->Branch("emcal_cluster_iso_subR03",&m_emcal_cluster_iso_subR03);
      towerntuple->Branch("emcal_cluster_iso_subR04",&m_emcal_cluster_iso_subR04);
    }

    if (ihcal) {
        towerntuple->Branch("hcalin_totenergy",&m_hcalin_totenergy);
        towerntuple->Branch("hcalin_energy",&m_hcalin_energy);
        towerntuple->Branch("hcalin_etabin",&m_hcalin_etabin);
        towerntuple->Branch("hcalin_phibin",&m_hcalin_phibin);
        towerntuple->Branch("hcalin_time",&m_hcalin_time);
        towerntuple->Branch("hcalin_etacorr",&m_hcalin_etacorr);
        towerntuple->Branch("hcalin_chi2",&m_hcalin_chi2);
    }

    if (ohcal) {
        towerntuple->Branch("hcalout_totenergy",&m_hcalout_totenergy);
        towerntuple->Branch("hcalout_energy",&m_hcalout_energy);
        towerntuple->Branch("hcalout_etabin",&m_hcalout_etabin);
        towerntuple->Branch("hcalout_phibin",&m_hcalout_phibin);
        towerntuple->Branch("hcalout_time",&m_hcalout_time);
        towerntuple->Branch("hcalout_etacorr",&m_hcalout_etacorr);
        towerntuple->Branch("hcalout_chi2",&m_hcalout_chi2);
    }

    if (zdc) {
        towerntuple->Branch("zdc_energy",&m_zdc_energy);
        towerntuple->Branch("zdc_index",&m_zdc_index);
        towerntuple->Branch("zdc_side",&m_zdc_side);
        towerntuple->Branch("zdc_type",&m_zdc_type);
    }

    if (mbd) {
        towerntuple->Branch("mbd_totcharge",&m_mbd_totcharge);
        towerntuple->Branch("mbd_totcharge_south",&m_mbd_totcharge_south);
        towerntuple->Branch("mbd_totcharge_north",&m_mbd_totcharge_north);
        towerntuple->Branch("mbd_charge_ch",&m_mbd_ch);
        towerntuple->Branch("nMBDHitNorth",&nMBDHitNorth);
        towerntuple->Branch("nMBDHitSouth",&nMBDHitSouth);
    }

    if (sepd) {
        towerntuple->Branch("sepd_energy",&m_sepd_energy);
        towerntuple->Branch("sepd_sector",&m_sepd_sector);
        towerntuple->Branch("sepd_channel",&m_sepd_channel);
    }


    if (isMC){
      towerntuple->Branch("truthpar_n",&truthpar_n);
      towerntuple->Branch("truth_energy",&truth_energy);
      towerntuple->Branch("truth_pt",&truth_pt);
      towerntuple->Branch("truth_eta",&truth_eta);
      towerntuple->Branch("truth_phi",&truth_phi);
      towerntuple->Branch("truth_id",&truth_id);
      if(isHI){
        towerntuple->Branch("npart",&npart);
        towerntuple->Branch("ncoll",&ncoll);
        towerntuple->Branch("bimp",&bimp);
        towerntuple->Branch("Cent_impactparam",&Cent_impactparam);
        towerntuple->Branch("Centrality",&Centrality);
      }
      towerntuple->Branch("signalId",&processId);
    }
    else if(!isMC){
      if(isHI){
        towerntuple->Branch("isMinBias",&isMinBias);
        towerntuple->Branch("centbin",&centbin);
      }
    }

}

int CaloAna::process_event(PHCompositeNode* topNode)
{
    if(!topNode){
      std::cerr << "CaloAna::Init - topnode PHCompositeNode not valid! return -1" << std::endl;
      return -1; 
    }
    process_towers(topNode);
    return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_towers(PHCompositeNode* topNode)
{
  ProcessGlobalEventInfo(topNode);
  if(isMC) ProcessFillTruthParticle(topNode);

  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if (!geomEM) throw std::runtime_error("failed to find EMCAL TOWERGEOM node in RawClusterDeadHotMask::CreateNodeTree");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if (!geomIH) throw std::runtime_error("failed to find IHCAL TOWERGEOM node in RawClusterDeadHotMask::CreateNodeTree");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if (!geomOH) throw std::runtime_error("failed to find OHCAL TOWERGEOM node in RawClusterDeadHotMask::CreateNodeTree");

  if(emcal){
    TowerInfoData emcalData(m_emcal_energy, m_emcal_time, m_emcal_etabin, m_emcal_phibin, m_emcal_etacorr, m_emcal_totenergy);
    TowerInfoContainer* offlinetowersEM = (isMC) ? static_cast<TowerInfoContainer*>(findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC")) : static_cast<TowerInfoContainer*>(findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC"));
    ProcessTowerInfoFill(geomEM, offlinetowersEM, emcalData, "CEMC");
  }
  if(ihcal){
    TowerInfoData hcalinData(m_hcalin_energy, m_hcalin_time, m_hcalin_etabin, m_hcalin_phibin, m_hcalin_etacorr, m_hcalin_totenergy);
    TowerInfoContainer* offlinetowersIH = (isMC) ? static_cast<TowerInfoContainer*>(findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALIN")) : static_cast<TowerInfoContainer*>(findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALIN"));
    ProcessTowerInfoFill(geomIH, offlinetowersIH, hcalinData, "HCALIN");
  }
  if(ohcal){
    TowerInfoData hcaloutData(m_hcalout_energy, m_hcalout_time, m_hcalout_etabin, m_hcalout_phibin, m_hcalout_etacorr, m_hcalout_totenergy);
    TowerInfoContainer* offlinetowersOH = (isMC) ? static_cast<TowerInfoContainer*>(findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALOUT")) : static_cast<TowerInfoContainer*>(findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALOUT"));
    ProcessTowerInfoFill(geomOH, offlinetowersOH, hcaloutData, "HCALOUT");
  }

  if (clusters) {
    RawClusterContainer *clustersEM = (isMC) ? static_cast<RawClusterContainer*>(findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC")) : static_cast<RawClusterContainer*>(findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC"));
    RawClusterContainer::ConstIterator hiter;
    RawClusterContainer::ConstRange begin_end = clustersEM->getClusters(); 

    int nCluster=0;
    for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      RawCluster* cluster = hiter->second;
      CLHEP::Hep3Vector vertex(vx, vy, vz);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      //float clPt = E_vec_cluster.perp();
      float clEta = E_vec_cluster.pseudoRapidity();
      float clPhi = E_vec_cluster.phi();
      float clE = E_vec_cluster.mag();
      float clChi2 = cluster->get_chi2();

      if(clE < _emcal_cluster_emincut) continue;

      m_emcal_cluster_E.push_back(cluster->get_energy());
      m_emcal_cluster_Ecore.push_back(clE);
      m_emcal_cluster_eta.push_back(clEta);
      m_emcal_cluster_phi.push_back(clPhi);
      m_emcal_cluster_ntower.push_back(cluster->getNTowers());
      m_emcal_cluster_chi2.push_back(clChi2);
      m_emcal_cluster_prob.push_back(cluster->get_prob());
      m_emcal_cluster_iso_unsubR01.push_back(cluster->get_et_iso(1,0,1));
      m_emcal_cluster_iso_unsubR02.push_back(cluster->get_et_iso(2,0,1));
      m_emcal_cluster_iso_unsubR03.push_back(cluster->get_et_iso(3,0,1));
      m_emcal_cluster_iso_unsubR04.push_back(cluster->get_et_iso(4,0,1));
      m_emcal_cluster_iso_subR01.push_back(cluster->get_et_iso(1,1,1));
      m_emcal_cluster_iso_subR02.push_back(cluster->get_et_iso(2,1,1));
      m_emcal_cluster_iso_subR03.push_back(cluster->get_et_iso(3,1,1));
      m_emcal_cluster_iso_subR04.push_back(cluster->get_et_iso(4,1,1));
      nCluster++;
    }


    m_emcal_cluster_E.shrink_to_fit();
    m_emcal_cluster_Ecore.shrink_to_fit();
    m_emcal_cluster_eta.shrink_to_fit();
    m_emcal_cluster_phi.shrink_to_fit();
    m_emcal_cluster_ntower.shrink_to_fit();
    m_emcal_cluster_chi2.shrink_to_fit();
    m_emcal_cluster_prob.shrink_to_fit();
    m_emcal_cluster_iso_unsubR01.shrink_to_fit();
    m_emcal_cluster_iso_unsubR02.shrink_to_fit();
    m_emcal_cluster_iso_unsubR03.shrink_to_fit();
    m_emcal_cluster_iso_unsubR04.shrink_to_fit();
    m_emcal_cluster_iso_subR01.shrink_to_fit();
    m_emcal_cluster_iso_subR02.shrink_to_fit();
    m_emcal_cluster_iso_subR03.shrink_to_fit();
    m_emcal_cluster_iso_subR04.shrink_to_fit();
  } 

  if(zdc) {
    TowerInfoContainer* offlinetowers = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_ZDC");
    int size = offlinetowers->size(); //online towers should be the same!

    for (int channel = 0; channel < size;channel++)
    {
      TowerInfo* offlinetower = offlinetowers->get_tower_at_channel(channel);
      if(offlinetower->get_isHot()) continue;
      if(!offlinetower->get_isGood()) continue;
      float offlineenergy = offlinetower->get_energy();

      m_zdc_energy.push_back(offlineenergy);
      m_zdc_side.push_back(channel / 8);
      m_zdc_index.push_back(channel);
      if (channel == 6 || channel == 7 || channel == 14 || channel == 15) m_zdc_type.push_back(0);
      else m_zdc_type.push_back(1);


      if(channel == 0 || channel == 2 || channel == 4 || channel == 8 || channel == 10 || channel == 12){
      }
    }
  }

  if (sepd) {
    TowerInfoContainer* offlinetowers = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_EPD");
    int size = offlinetowers->size(); //online towers should be the same!

    int overflow=0;
    for (int channel = 0; channel < size;channel++)
    {
      TowerInfo* offlinetower = offlinetowers->get_tower_at_channel(channel);
      if(offlinetower->get_isHot()) continue;
      if(offlinetower->get_isBadChi2()) break;
      float offlineenergy = offlinetower->get_energy();

      if(channel<64) continue;
      tile = -1;//sepd_mapping_function(channel);

      m_sepd_energy.push_back(offlineenergy);
      m_sepd_channel.push_back(tile);

      if(channel<96) m_sepd_sector.push_back(5);
      else if(channel>=96) m_sepd_sector.push_back(0);
      else m_sepd_sector.push_back(-1);

      if(offlineenergy>sepdHighECut){
        overflow++;
      }

    }
  }
  if (mbd) {
    m_mbd_totcharge=0;
    m_mbd_totcharge_south=0;
    m_mbd_totcharge_north=0;
    PHNodeIterator iter(topNode);
    PHCompositeNode *mbdNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "MBD"));
    if(!mbdNode){std::cout << "no MBD node" << std::endl; return Fun4AllReturnCodes::ABORTEVENT;}
    MbdPmtContainer *mbdpmts = findNode::getClass<MbdPmtContainerV1>(mbdNode, "MbdPmtContainer");
    if(!mbdpmts){std::cout << "no MbdPmtContainer" << std::endl; return Fun4AllReturnCodes::ABORTEVENT;}

    nMBDHitNorth=0;
    nMBDHitSouth=0;
    m_mbd_ch=0;
    for(int ipmt=0;ipmt<128; ++ipmt){ //mbdpmts->get_npmts() function not working..
      MbdPmtHit *mbdhit = mbdpmts->get_pmt(ipmt);
      m_mbd_totcharge += mbdhit->get_q();
      if(ipmt==2) m_mbd_ch = mbdhit->get_q();
      if(ipmt<64 && mbdhit->get_q()>0.4){ nMBDHitNorth++; m_mbd_totcharge_north+=mbdhit->get_q();}
      if(ipmt>=64 && mbdhit->get_q()>0.4){ nMBDHitSouth++; m_mbd_totcharge_south+=mbdhit->get_q();}
    }
  }

  towerntuple->Fill();
  ProcessClearBranchVar();

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::ProcessGlobalEventInfo(PHCompositeNode* topNode){

  //vertex
  vz=-999;
  vx=-999;
  vy=-999;

  if(isMC){
    MbdVertexMap *mbdvtxmap = findNode::getClass<MbdVertexMap>(topNode,"MbdVertexMap");
    bool isglbvtx=true;
    if(!mbdvtxmap){std::cout << "no mbdvertex map..."<<std::endl;}
    if(mbdvtxmap->empty()){ 
      std::cout << "Empty mbdmap or mbdvtx node" << std::endl;
      isglbvtx=false;
    }
    if(isglbvtx){
      MbdVertex *bvertex= nullptr;
      if (mbdvtxmap && m_UseAltZVertex == 1)
      {
        for (MbdVertexMap::ConstIter mbditer= mbdvtxmap->begin(); mbditer != mbdvtxmap->end(); ++mbditer)
        {
          bvertex = mbditer->second;
        }
        if(!bvertex){std::cout << "could not find globalvtxmap iter :: set vtx to (-999,-999,-999)" << std::endl;}
        else if(bvertex){
          vz = bvertex->get_z();
          vy = bvertex->get_y();
          vx = bvertex->get_x();
        }
      }
    }
  }
  else if(!isMC){
    GlobalVertexMapv1 *globalvtxmap = findNode::getClass<GlobalVertexMapv1>(topNode,"GlobalVertexMap");
    bool isglbvtx=true;
    if(!globalvtxmap || globalvtxmap->empty()){ 
      std::cout << "Empty mbdmap " << std::endl;
      isglbvtx=false;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if(isglbvtx){
      GlobalVertex *bvertex= nullptr;
      if (globalvtxmap && m_UseAltZVertex == 1)
      {
        for (GlobalVertexMap::ConstIter globaliter= globalvtxmap->begin(); globaliter != globalvtxmap->end(); ++globaliter)
        {
          bvertex = globaliter->second;
        }
        if(!bvertex){std::cout << "could not find globalvtxmap iter :: set vtx to (-999,-999,-999)" << std::endl;}
        else if(bvertex){
          vz = bvertex->get_z();
          vy = bvertex->get_y();
          vx = bvertex->get_x();
        }
      }
    }
  }
    
  //Event info
  if(isMC){
    PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    if(!genevtmap){std::cout << "no PHHepMCGenEventMap" << std::endl; return Fun4AllReturnCodes::ABORTEVENT;}
    PHHepMCGenEvent *genevt = nullptr;

    for (PHHepMCGenEventMap::ReverseIter iter = genevtmap->rbegin(); iter != genevtmap->rend(); ++iter)
    {   
      genevt = iter->second;
      if(!genevt) {cout<<"ERROR: no PHHepMCGenEvent!" << endl; return Fun4AllReturnCodes::ABORTEVENT;}
    }   

    HepMC::GenEvent *event = genevt->getEvent();
    if (!event) { cout << PHWHERE << "ERROR: no HepMC::GenEvent!" << endl; return Fun4AllReturnCodes::ABORTEVENT;}
    processId = event->signal_process_id();

    if(isHI){
      EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
      if ( event_header ) {
        npart = event_header->get_intval("npart");
        ncoll = event_header->get_intval("ncoll");
        bimp = event_header->get_floatval("bimp");
      } 
      CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
      if (!cent_node)
      {
        std::cout
          << "Error can not find centrality node "
          << std::endl;
        exit(-1);
      }   
      Centrality =  cent_node->get_centile(CentralityInfo::PROP::mbd_NS);
      Cent_impactparam =  cent_node->get_quantity(CentralityInfo::PROP::bimp);
    }
  }
  else if(!isMC){
    if(isHI){
      CentralityInfo *centrality = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
      if (!centrality)
      {
        std::cout << "no centrality node " << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      MinimumBiasInfo *minimumbiasinfo = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
      if (!minimumbiasinfo)
      {
        std::cout << "no minimumbias node " << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      float centile = (centrality->has_centile(CentralityInfo::PROP::mbd_NS) ? centrality->get_centile(CentralityInfo::PROP::mbd_NS) : -999.99);
      centbin = centile*100;
      isMinBias = minimumbiasinfo->isAuAuMinimumBias(); 
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloAna::ProcessFillTruthParticle(PHCompositeNode* topNode){
    PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");
    if(!truthinfo){std::cout << "no truth info node... just skip the whole part.." << std::endl; return;}
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {

      PHG4Particle* g4particle = iter->second;
      if(isHI){
        if(truthinfo->isEmbeded(g4particle->get_track_id())<=0) continue;
      }

      TLorentzVector t;
      t.SetPxPyPzE (g4particle->get_px (), g4particle->get_py (), g4particle->get_pz (), g4particle->get_e ());

      // take only particles from primary Pythia event                                                            
      truth_energy.push_back(t.E());
      truth_pt.push_back(t.Pt());
      truth_eta.push_back(t.Eta());
      truth_phi.push_back(t.Phi());
      truth_id.push_back(g4particle->get_pid());
      truthpar_n++;
    }
}

void CaloAna::ProcessTowerInfoFill(RawTowerGeomContainer* geom, TowerInfoContainer* offlinetowers, TowerInfoData& data, const std::string& sectionName){
  int size = offlinetowers->size();

  data.totalEnergy = 0;
  for (int channel = 0; channel < size; channel++) {
    TowerInfo* offlinetower = offlinetowers->get_tower_at_channel(channel);

    if (!isMC) {
      if (offlinetower->get_isHot()) continue;
      if (!offlinetower->get_isGood()) continue;
    }

    float offlineenergy = offlinetower->get_energy();
    data.totalEnergy += offlineenergy;

    unsigned int towerkey = offlinetowers->encode_key(channel);
    int ieta = offlinetowers->getTowerEtaBin(towerkey);
    int iphi = offlinetowers->getTowerPhiBin(towerkey);
    float _time = offlinetowers->get_tower_at_channel(channel)->get_time_float();

    RawTowerDefs::CalorimeterId caloId;
    if (sectionName == "HCALOUT") {
          caloId = RawTowerDefs::CalorimeterId::HCALOUT;
    } else if (sectionName == "HCALIN") {
          caloId = RawTowerDefs::CalorimeterId::HCALIN;
    } else if (sectionName == "CEMC") {
          caloId = RawTowerDefs::CalorimeterId::CEMC;
    }
    const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(caloId, ieta, iphi);
    RawTowerGeom* tower_geom = geom->get_tower_geometry(geomkey);

    float posx = tower_geom->get_center_x();
    float posy = tower_geom->get_center_y();
    float posz = tower_geom->get_center_z();

    float newx = posx - vx;
    float newy = posy - vy;
    float newz = posz - vz;

    data.etabinVec.push_back(ieta);
    data.phibinVec.push_back(iphi);
    data.etacorrVec.push_back(asinh(newz / sqrt(newx * newx + newy * newy)));
    data.energyVec.push_back(offlineenergy);
    data.timeVec.push_back(_time);
  }
  data.energyVec.shrink_to_fit();
  data.etabinVec.shrink_to_fit();
  data.phibinVec.shrink_to_fit();
  data.etacorrVec.shrink_to_fit();
}

void CaloAna::ProcessClearBranchVar(){
  if (emcal) {
    m_emcal_etabin.clear();
    m_emcal_phibin.clear();
    m_emcal_energy.clear();
    m_emcal_time.clear();
    m_emcal_etacorr.clear();
    m_emcal_chi2.clear();
    m_emcal_pedestal.clear();
  }
  if(clusters){
    m_emcal_cluster_E.clear();
    m_emcal_cluster_Ecore.clear();
    m_emcal_cluster_eta.clear();
    m_emcal_cluster_phi.clear();
    m_emcal_cluster_ntower.clear();
    m_emcal_cluster_chi2.clear();
    m_emcal_cluster_prob.clear();
    m_emcal_cluster_iso_unsubR01.clear();
    m_emcal_cluster_iso_unsubR02.clear();
    m_emcal_cluster_iso_unsubR03.clear();
    m_emcal_cluster_iso_unsubR04.clear();
    m_emcal_cluster_iso_subR01.clear();
    m_emcal_cluster_iso_subR02.clear();
    m_emcal_cluster_iso_subR03.clear();
    m_emcal_cluster_iso_subR04.clear();
  }
  
  if (ihcal) {
    m_hcalin_etabin.clear();
    m_hcalin_phibin.clear();
    m_hcalin_energy.clear();
    m_hcalin_time.clear();
    m_hcalin_etacorr.clear();
    m_hcalin_chi2.clear();
  }
  if (ohcal) {
    m_hcalout_etabin.clear();
    m_hcalout_phibin.clear();
    m_hcalout_energy.clear();
    m_hcalout_time.clear();
    m_hcalout_etacorr.clear();
    m_hcalout_chi2.clear();
  }
  if (zdc) {
    m_zdc_energy.clear();
    m_zdc_index.clear();
    m_zdc_side.clear();
    m_zdc_type.clear();
  }
  if (sepd) {
    m_sepd_energy.clear();
    m_sepd_channel.clear();
    m_sepd_sector.clear();
  }
  if (mbd) {
    m_mbd_energy.clear();
    m_mbd_ch=0; 
    m_mbd_totcharge=0;
    m_mbd_totcharge_south=0;
    m_mbd_totcharge_north=0;
  }

  if(isMC){
    truthpar_n = 0;
    truth_energy.clear();
    truth_pt.clear();
    truth_eta.clear();
    truth_phi.clear();
    truth_id.clear();
  }
}

int CaloAna::End(PHCompositeNode* /*topNode*/)
{
  std::cout << "end..! " << std::endl;
  if(outfile){
    outfile->cd();
    towerntuple->Write();
    outfile->Write();
    outfile->Close();
    delete outfile;
    outfile=nullptr;
  }
  return 0;
}
