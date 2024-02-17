#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>

#include <G4_Global.C>
//#include <calotowerbuilder/CaloTowerBuilder.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloWaveformProcessing.h>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <clusteriso/ClusterIso.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>

// #include <runtowerinfo/RunTowerInfo.h>
#include <caloana/CaloAna.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/GlobalVertexMap.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4centrality/PHG4CentralityReco.h>

#include <centrality/CentralityReco.h>
#include <calotrigger/MinimumBiasClassifier.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libcaloana.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4vertex.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libclusteriso.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libFROG.so)


#endif

void Fun4All_macro(const char* infile="DST_CALO_run1auau_ana402_2023p009-00023727-0006.root", bool isMC=false, bool isHI=false)
{
    string outdir = "outputfiles";
    void * dirf = gSystem->OpenDirectory(outdir.c_str());
    if(dirf) gSystem->FreeDirectory(dirf);
    else {gSystem->mkdir(outdir.c_str(), kTRUE);}
    const char *outfile = Form("%s/outtree_%s",outdir.c_str(),infile);
    Fun4AllServer *se = Fun4AllServer::instance();
    int verbosity = 0;

    se->Verbosity(verbosity);
    recoConsts *rc = recoConsts::instance();

    //===============
    // conditions DB flags
    //===============

    // global tag
    rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2023"); 
    rc->set_uint64Flag("TIMESTAMP",23727);

    if(isMC && isHI){
      Global_Reco();
      PHG4CentralityReco *cent = new PHG4CentralityReco();
      cent->Verbosity(verbosity);
      cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
      se->registerSubsystem( cent );
    }

    RetowerCEMC *rcemc = new RetowerCEMC();
    rcemc->Verbosity(verbosity);
    rcemc->set_towerinfo(true);
    se->registerSubsystem(rcemc);

    JetReco *towerjetreco = new JetReco();
    towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER));
    towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO));
    towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO));
    towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.2), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
    towerjetreco->set_algo_node("ANTIKT");
    towerjetreco->set_input_node("TOWER");
    towerjetreco->Verbosity(verbosity);
    se->registerSubsystem(towerjetreco);

    DetermineTowerBackground *dtb = new DetermineTowerBackground();
    dtb->SetBackgroundOutputName("TowerInfoBackground_Sub1");
    dtb->SetFlow(false);
    dtb->SetSeedType(0);
    dtb->SetSeedJetD(3);
    dtb->set_towerinfo(true);
    dtb->Verbosity(verbosity);
    se->registerSubsystem(dtb);

    CopyAndSubtractJets *casj = new CopyAndSubtractJets();
    casj->SetFlowModulation(false);
    casj->Verbosity(verbosity);
    casj->set_towerinfo(true);
    se->registerSubsystem(casj);

    DetermineTowerBackground *dtb2 = new DetermineTowerBackground();
    dtb2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
    dtb2->SetFlow(false);
    dtb2->SetSeedType(1);
    dtb2->SetSeedJetPt(7);
    dtb2->Verbosity(verbosity);
    dtb2->set_towerinfo(true);
    se->registerSubsystem(dtb2);

    SubtractTowers *st = new SubtractTowers();
    st->SetFlowModulation(false);
    st->Verbosity(verbosity);
    st->set_towerinfo(true);
    se->registerSubsystem(st);

    ClusterIso *cliso1 = new ClusterIso("ClusterIso1",0.3,1,true,true);
    cliso1->Verbosity(verbosity);
    se->registerSubsystem( cliso1 );
    ClusterIso *cliso2 = new ClusterIso("ClusterIso2",0.3,2,true,true);
    cliso2->Verbosity(verbosity);
    se->registerSubsystem( cliso2 );
    ClusterIso *cliso3 = new ClusterIso("ClusterIso3",0.3,3,true,true);
    cliso3->Verbosity(verbosity);
    se->registerSubsystem( cliso3 );
    ClusterIso *cliso4 = new ClusterIso("ClusterIso4",0.3,4,true,true);
    cliso4->Verbosity(verbosity);
    se->registerSubsystem( cliso4 );

    std::cout << "infile : " << infile <<std::endl;
    Fun4AllInputManager *intrue = new Fun4AllDstInputManager("DST_TOWERS");
    intrue->AddFile(infile);
    se->registerInputManager(intrue);

    CaloAna *ca = new CaloAna("caloana",outfile,isMC,isHI);
    se->registerSubsystem(ca);

    std::cout << "now run..." << std::endl;
    se->run();
    se->End();
    std::cout << "ok done.. " << std::endl;

    delete se;
}
