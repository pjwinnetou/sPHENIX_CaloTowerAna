#ifndef CALOANA_H__
#define CALOANA_H__

// Utility
#include <vector>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH2.h>
#include <cassert>
#include <sstream>
#include <string>
#include <TLorentzVector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

// Fun4All
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

// Event
#include <Event/Event.h>
#include <Event/packet.h>
#include <ffaobjects/EventHeaderv1.h>

// Jet base
#include <jetbase/Jet.h>
#include <jetbase/Jetv1.h>
#include <jetbase/JetAlgo.h>
#include <jetbase/JetMap.h>
#include <jetbackground/TowerBackgroundv1.h>

// Global vertex
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>

// G4
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfov2.h>
#include <calobase/TowerInfoDefs.h>

// MBD
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdGeom.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>

// phool
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

// HepMC
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>              // for GenEvent::particle_const_ite...
#include <HepMC/GenVertex.h>             // for GenVertex, GenVertex::partic...
#include <HepMC/GenParticle.h>           // for GenParticle
#include <HepMC/IteratorRange.h>         // for ancestors, children, descend...
#include <HepMC/SimpleVector.h>   // for FourVector
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

// Centrality MB
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <centrality/CentralityInfov1.h>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH2F;
class TH1;

class CaloAna : public SubsysReco
{
 public:
  //! constructor
  CaloAna(const std::string &name = "CaloAna",  const char* outname="DST-00021615-0000.root", bool _isMC=false, bool _iHI=false);

  //! destructor
  virtual ~CaloAna();

  //! full initialization
  int Init(PHCompositeNode *);
  void InitOutputFile();
  void InitTree();
  
  //! event processing method
  int process_event(PHCompositeNode *);
  int ProcessGlobalEventInfo(PHCompositeNode *);
  void ProcessFillTruthParticle(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);
  int Getpeaktime(TH1 *h);

  void Detector(const std::string &name) { detector = name; }
  void ProcessClearBranchVar();

  struct TowerInfoData{
    std::vector<float>& energyVec;
    std::vector<float>& timeVec;
    std::vector<int>& etabinVec;
    std::vector<int>& phibinVec;
    std::vector<float>& etacorrVec;
    float& totalEnergy;

    TowerInfoData(std::vector<float>& energy, std::vector<float>& time, std::vector<int>& etabin, std::vector<int>& phibin, std::vector<float>& etacorr, float& totEnergy)
      : energyVec(energy), timeVec(time), etabinVec(etabin), phibinVec(phibin), etacorrVec(etacorr), totalEnergy(totEnergy) {}
  };
  void ProcessTowerInfoFill(RawTowerGeomContainer* geom, TowerInfoContainer* offlinetowers, TowerInfoData& data, const std::string& sectionName); 

//  int sepd_mapping_function(int input);


  
 protected:
  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm = nullptr;
  TFile *outfile = nullptr;
  TFile *calibfile= nullptr;
  TFile *tcalibfile= nullptr;
  TNtuple *calibntuple= nullptr;
  TNtuple *tcalibntuple= nullptr;
  TNtuple *g4hitntuple = nullptr;
  TNtuple *g4cellntuple = nullptr;
  TTree *towerntuple = nullptr;
  TNtuple *clusterntuple = nullptr;

  int eventcount=0;
  double z_offset=28.3;

  float vx, vy, vz;
  int centbin;
  bool isMinBias;
  
  float m_emcal_totenergy;
  std::vector<float> m_emcal_energy;
  std::vector<int> m_emcal_etabin;
  std::vector<int> m_emcal_phibin;
  std::vector<float> m_emcal_time;
  std::vector<float> m_emcal_etacorr;
  std::vector<float> m_emcal_chi2;
  std::vector<float> m_emcal_pedestal;

  std::vector<float> m_emcal_cluster_E;
  std::vector<float> m_emcal_cluster_Ecore;
  std::vector<float> m_emcal_cluster_eta;
  std::vector<float> m_emcal_cluster_phi;
  std::vector<float> m_emcal_cluster_ntower;
  std::vector<float> m_emcal_cluster_chi2;
  std::vector<float> m_emcal_cluster_prob;
  std::vector<float> m_emcal_cluster_iso_unsubR01;
  std::vector<float> m_emcal_cluster_iso_unsubR02;
  std::vector<float> m_emcal_cluster_iso_unsubR03;
  std::vector<float> m_emcal_cluster_iso_unsubR04;
  std::vector<float> m_emcal_cluster_iso_subR01;
  std::vector<float> m_emcal_cluster_iso_subR02;
  std::vector<float> m_emcal_cluster_iso_subR03;
  std::vector<float> m_emcal_cluster_iso_subR04;

  float m_hcalin_totenergy;
  std::vector<float> m_hcalin_energy;
  std::vector<int> m_hcalin_etabin;
  std::vector<int> m_hcalin_phibin;
  std::vector<float> m_hcalin_time;
  std::vector<float> m_hcalin_etacorr;
  std::vector<float> m_hcalin_chi2;

  float m_hcalout_totenergy;
  std::vector<float> m_hcalout_energy;
  std::vector<int> m_hcalout_etabin;
  std::vector<int> m_hcalout_phibin;
  std::vector<float> m_hcalout_time;
  std::vector<float> m_hcalout_etacorr;
  std::vector<float> m_hcalout_chi2;

  std::vector<float> m_zdc_energy;
  std::vector<int> m_zdc_index;
  std::vector<int> m_zdc_side;
  std::vector<float> m_zdc_type;

  std::vector<float> m_mbd_energy;
  float m_mbd_totcharge;
  float m_mbd_totcharge_south;
  float m_mbd_totcharge_north;
  bool minbias=false;

  std::vector<float> m_sepd_energy;
  std::vector<int> m_sepd_sector;
  std::vector<int> m_sepd_channel;

  int truthpar_n = 0;
  std::vector<float> truth_energy;
  std::vector<float> truth_pt;
  std::vector<float> truth_eta;
  std::vector<float> truth_phi;
  std::vector<int> truth_id;

  float _emcal_cluster_emincut = 0.5;

  int npart;
  int ncoll;
  float bimp;

  int Centrality;
  float Cent_impactparam;

  bool isMC;
  bool isHI;

  int processId;

  float gaincorr[128] {};
  float time_shift_corr[128] {};
  float toa_shift_corr[128] {};
  float time_scale_corr[128] {};
  float _mbd_charge_threshold=0.4;
  int tile=-1;
  float sepdHighECut=14000;

  double tthresh = 16;
  double cthresh = 0.4;
  int central_cut = 4;
  float sigma_cut = 1.5;

  std::vector<float> time_sum_n;
  std::vector<float> time_sum_s;
  float sum_n = 0.;
  float sum_s = 0.;
  float sum_n2 = 0.;
  float sum_s2 = 0.;
  int hits_s_t=0;
  int hits_n_t=0;
  int nHitNorth=0;
  int nHitSouth=0;

  int nEvent[3]= {0,0,0};
  int nEventPass[3]= {0,0,0};
  
  int m_UseAltZVertex = 1;

};

#endif
