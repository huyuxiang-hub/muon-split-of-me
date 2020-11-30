//--------------------------------------------------------------------------
//                            junoSD_PMT_v2
//
// PMTs are difined as sensitive detector. They collect hits on them.
// The data members of hits are set up here using the information of G4Step.
// -------------------------------------------------------------------------
// Author: Liang Zhan, 2006/01/27
// Modified by: Weili Zhong, 2006/03/01
// -------------------------------------------------------------------------

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"

#include "BufferMemMgr/IDataMemMgr.h"
#include "G4RunManager.hh"
#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavigator.h"

#include "junoSD_PMT_muon.hh"
//#include "junoHit_PMT.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include <cassert>
#include "NormalTrackInfo.hh"
#include "G4OpticalPhoton.hh"
#include "Randomize.hh"
#include "G4DataInterpolation.hh"
#include "G4VProcess.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"

//#include <TTimeStamp.h>

/*#ifdef WITH_G4OPTICKS
#include <iomanip>
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Opticks.hh"

#include "PMTEfficiency.hh"
#include "PMTEfficiencyTable.hh"

#define WITH_G4OPTICKS_EFFICIENCY_CHECK 1

#endif
*/

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////////////

junoSD_PMT_muon::junoSD_PMT_muon(const std::string& name)
:G4VSensitiveDetector(name)
{
/*    G4String HCname;
    collectionName.insert(HCname="hitCollection");
    collectionName.insert(HCname="hitCollectionMuon");
#ifdef WITH_G4OPTICKS
    if(m_opticksMode > 0){
        collectionName.insert(HCname="hitCollectionOpticks");
    }
    m_PMTEfficiency = new PMTEfficiency ; 
    m_PMTEfficiencyTable = new PMTEfficiencyTable(m_PMTEfficiency) ; 
#endif
  */
   minmax_initialized = false;
   max_CDLPMT_hittime = 0;
   min_CDLPMT_hittime = 0; 
  
   m_timewindow=0;
   n_hit=0;

    m_debug = true;
  //  m_time_window = 1; // 1ns

  //  m_pmthitmerger = 0;
  //  m_pmthitmerger_opticks = 0;
   // m_hit_type = 1; // 1 for normal, 2 for muon

    m_qescale = 1.0/0.8 ;  // scale back for 0.8 QE-scale in LsExpDetectorConstruction
    m_angle_response = 1.0;


    m_ce_flat_value = 0.9;

    MCP20inch_m_ce_flat_value = 0.85;
    MCP8inch_m_ce_flat_value = 0.85;
    Ham20inch_m_ce_flat_value = 0.95;
    Ham8inch_m_ce_flat_value = 0.7;
    HZC9inch_m_ce_flat_value = 0.67;

    MCP20inch_m_EAR_value = 1.;
    MCP8inch_m_EAR_value = 1.;
    Ham20inch_m_EAR_value = 0.93;
    Ham8inch_m_EAR_value = 0.88;
    HZC9inch_m_EAR_value = 0.92;

    m_disable = false;
    // 20inchfunc (1D)
    // FIXME: following are not used in current code
    m_ce_func_str = "0.9*[0]/(1+[1]*exp(-[2]*x))";
    //OLD
    //m_ce_func_params.push_back(1.006); // p0
    //m_ce_func_params.push_back(0.9023);// p1
    //m_ce_func_params.push_back(0.1273);// p2

    //NEW from the PMT test benches at Zhongshan
    m_ce_func_params.push_back(0.9194); // p0
    m_ce_func_params.push_back(0.504);  // p1
    m_ce_func_params.push_back(0.08076);// p2

    //These fit the Hamamatsu CE better..
    //m_ce_func_params.push_back(1.02557);// p0
    //m_ce_func_params.push_back(6.77639);// p1
    //m_ce_func_params.push_back(0.16419);// p2

    m_ce_func = 0;
    m_merge_count = 0 ; 
    //split out put
    iotaskname="muon_io";
    iotask=0;
    maxhit=100;
    hit_count=0;
  
    m_split=true;

}

junoSD_PMT_muon::~junoSD_PMT_muon()
{;}

void junoSD_PMT_muon::Initialize(G4HCofThisEvent *HCE)
{
      LogInfo<<"begin of event action"<<std::endl;
      sim_header =new JM::SimHeader;
      sim_event = new JM::SimEvent;

     // m_PMThit.clear();

   minmax_initialized = false;
   max_CDLPMT_hittime = 0;
   min_CDLPMT_hittime = 0;
   
   Task* m_scope = getScope();
   std::cout<<"m_split="<<m_split<<std::endl;
   if(m_split){
    iotask = dynamic_cast<Task*>(m_scope->find(iotaskname));
    if (iotask == 0) {
        LogError << "Can't find the task for I/O." << std::endl;
        throw SniperException("Make sure the IO task is created");
    }
    // check the BufferMgr in iotask
    SniperPtr<IDataMemMgr> mMgr(*iotask, "BufferMemMgr");
     if ( mMgr.invalid() ) {
           LogError << "Failed to get BufferMemMgr!" << std::endl;
           throw SniperException("Make sure you have load the BufferMemMgr.");
         }
    m_bufmgr=mMgr.data();
    }
   
    else {
          
            SniperPtr<IDataMemMgr> mMgr(*getScope(), "BufferMemMgr");
            if ( mMgr.invalid() ) {
           LogError << "Failed to get BufferMemMgr!" << std::endl;
           throw SniperException("Make sure you have load the BufferMemMgr.");
           }
         m_bufmgr=mMgr.data();
 
      std::cout<<"Task buffer"<<std::endl;  
     }
   
    // m_PMTParamsvc = 0;
    
   /* LogInfo << "Retrieving PMTParamSvc." << std::endl;
    SniperPtr<PMTParamSvc> svc(*getScope(), "PMTParamSvc");
      if (svc.invalid()) { 
        LogError << "Can't get PMTParamSvc. We can't initialize PMT." << std::endl;
        assert(0);
      } else {
        LogInfo << "Retrieve PMTParamSvc successfully." << std::endl;
        m_PMTParamsvc = svc.data();
      } 
   */

}

G4bool junoSD_PMT_muon::ProcessHits(G4Step * step,G4TouchableHistory*)
{
   //std::cout<<"G4UniformRand"<<G4UniformRand()<<std::endl; 
  // return true; 
  
   if (m_disable) {
        return false;
    }
    // TODO: now it only support the single PE.
    // = only accept the optical photon
    G4Track* track = step->GetTrack();
    if (track->GetDefinition() != G4OpticalPhoton::Definition()) {
        return false;
    }
    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
    double edep = step->GetTotalEnergyDeposit();
    // = only when the photon is detected by the surface, the edep is non-zero.
    // = the QE is already applied in the OpBoundaryProcess::DoAbsorption
    if (edep<=0.0) {
        return false;
    }

    // LT
    // = Due to update of Geant4, now OpAbsorption will also cause non-zero edep.
    // = Hence we need to check the OP boundary.
    G4bool isOnBoundary = (postStepPoint->GetStepStatus() == fGeomBoundary);
    if (not isOnBoundary) {
        return false;
    }

    static G4ThreadLocal G4OpBoundaryProcess* boundary_proc=NULL;
    if (!boundary_proc) {
        G4ProcessManager* OpManager =
            G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
        if (OpManager) {
            G4int MAXofPostStepLoops =
                OpManager->GetPostStepProcessVector()->entries();
            G4ProcessVector* fPostStepDoItVector =
                OpManager->GetPostStepProcessVector(typeDoIt);
            for ( G4int i=0; i<MAXofPostStepLoops; i++) {
                G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
                G4OpBoundaryProcess* op =  dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
                if (op) {
                    boundary_proc = op;
                    break;
                }
            }
        }
     
    }

    if (!boundary_proc) {
        G4cout << "Can't locate OpBoundaryProcess." << G4endl;
        return false;
    }

    G4OpBoundaryProcessStatus theStatus = Undefined;
    theStatus = boundary_proc->GetStatus();

    if (theStatus != Detection) {
        return false;
    }

 

    // TODO: get CE and angle response from data.
    // = decide the CE (collection efficiency)
    // = the CE can be different at different position
    // == volume name
    std::string volname = track->GetVolume()->GetName(); // physical volume
    // == position
    const G4AffineTransform& trans = track->GetTouchable()->GetHistory()->GetTopTransform();
    const G4ThreeVector& global_pos = postStepPoint->GetPosition();
    G4ThreeVector local_pos = trans.TransformPoint(global_pos);
  
    double qe = 1;
    bool pmt_type = true;
    bool qe_type  = true;


    // == get the copy number -> pmt id
    int pmtID = get_pmtid(track);
    if(pmtID<18000){
      qe =  m_PMTParamsvc->getPMTQE(pmtID);
      pmt_type =  m_PMTParamsvc->isHamamatsu(pmtID);
      qe_type  =  m_PMTParamsvc->isHighQE(pmtID);
    }else if(pmtID<300000){
      qe = 0.3;
    }else if(pmtID>=300000){
      qe = m_PMTParamsvc->getPMTQE(pmtID);
    }

    double ce = get_ce(volname, local_pos, pmt_type, qe_type);  

    double f_angle_response = 1.0;
    // = final DE = QE * CE, but QE is already applied, so only CE is important.
    // = DE: Detection Efficiency
    double de = qe*ce*f_angle_response*m_qescale;
    //double de = qe*ce*f_angle_response;
    //std::cout << "test de: " << qe*ce << std::endl;

   if (G4UniformRand() > de) {
        return false;
    }
 
    // ========================================================================
    // create the transient PMT Hit Object
    // ========================================================================
    // == get the copy number -> pmt id
    int pmtid = get_pmtid(track);
    

     double hittime = postStepPoint->GetGlobalTime();
    int copyno =pmtid;
   

    if ((copyno < 30000) or (copyno >= 300000)) 
    {
        n_hit++;
    }
    if(copyno < 30000)
       {
         if (!minmax_initialized)
            {
              max_CDLPMT_hittime =hittime;
              min_CDLPMT_hittime =hittime;
              minmax_initialized = true;
            }
         else
            {
              if (hittime < min_CDLPMT_hittime)
                 {
                   min_CDLPMT_hittime = hittime;
                 }
               if (hittime > max_CDLPMT_hittime)
                 {
                   max_CDLPMT_hittime = hittime;
                 }
            }
        }
    
  
    /* if(m_merge_flag)
       {
         bool ok= domerge(pmtid,hittime);
         if (ok)
            {return true;}  
       
       }*/

    
    hit_count++;   
     
    sim_hit=0;
     if ((copyno < 30000) or (copyno >= 300000)) 
        {
           sim_hit = sim_event->addCDHit();
        } 
     else if (copyno >= 30000) 
        {
           sim_hit = sim_event->addWPHit();
        }

        sim_hit->setPMTID(pmtid);
        sim_hit->setHitTime(hittime);
        sim_hit->setLocalTheta(local_pos.theta());
        sim_hit->setLocalPhi(local_pos.phi());
        sim_hit->setNPE(1);
       // m_PMThit[pmtid].push_back(sim_hit); 
      
     if (hit_count==maxhit and m_split ){
      std::cout<<"muon_io is start,fire!!"<<std::endl;
      JM::EvtNavigator* nav = new JM::EvtNavigator();
      TTimeStamp ts;
      nav->setTimeStamp(ts);
      
      m_bufmgr->adopt(nav, "/Event");
      sim_header->setEvent(sim_event);
  
      sim_header->setCDLPMTtimeWindow(0);
      sim_header->setCDLPMTtotalHits(0);
      nav->addHeader("/Event/Sim", sim_header);
      Incident::fire(*getScope(), iotaskname);
       
      sim_header =new JM::SimHeader;
      sim_event = new JM::SimEvent;
      
      hit_count=0;
      }





      

    return true;  
}




void junoSD_PMT_muon::EndOfEvent(G4HCofThisEvent* HCE){
    
      LogInfo << "Begin Fill Tracks" << std::endl;
  const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent() ;

   G4int nVertex = event -> GetNumberOfPrimaryVertex();
   
   for (G4int index=0; index < nVertex; ++index) 
      {
         G4PrimaryVertex* vtx = event->GetPrimaryVertex( index );
    
        double x = vtx->GetX0();
        double y = vtx->GetY0();
        double z = vtx->GetZ0();
        double t = vtx->GetT0();
      
         G4PrimaryParticle* pp = vtx -> GetPrimary();
         while (pp) 
           {

            int trkid = pp -> GetTrackID();
            int pdgid = pp -> GetPDGcode();
            double px = pp -> GetPx();
            double py = pp -> GetPy();
            double pz = pp -> GetPz();
            double mass = pp -> GetMass();


            JM::SimTrack* jm_trk = sim_event->addTrack();
            jm_trk->setPDGID(pdgid);
            jm_trk->setTrackID(trkid);
            jm_trk->setInitPx(px);
            jm_trk->setInitPy(py);
            jm_trk->setInitPz(pz);
            jm_trk->setInitMass(mass);
            jm_trk->setInitX(x);
            jm_trk->setInitY(y);
            jm_trk->setInitZ(z);
            jm_trk->setInitT(t);

            pp = pp->GetNext();

           }
      }
    

   SniperDataPtr<JM::NavBuffer>  navBuf(*m_scope, "/Event");
   if (navBuf.invalid()) {
        LogError << "Can't find the NavBuffer." << std::endl;
        return;
    }
   
  JM::EvtNavigator*  evt_nav = navBuf->curEvt();
   if (not evt_nav) {
        LogError << "Can't find the event navigator." << std::endl;
        return;
    }
  
   JM::EvtNavigator* nav = new JM::EvtNavigator();
   TTimeStamp ts = evt_nav->TimeStamp();
  // TTimeStamp ts;
   nav->setTimeStamp(ts);
   LogInfo << "current Timestamp: '"
            << ts
            << "'." << std::endl;

 // SniperPtr<IDataMemMgr> mMgr(*iotask, "BufferMemMgr");
   m_bufmgr->adopt(nav, "/Event");

   LogInfo<< "End Fill Tracks" << std::endl;

   m_timewindow=max_CDLPMT_hittime - min_CDLPMT_hittime;
    
     if (m_timewindow>0) 
       {
        sim_header->setCDLPMTtimeWindow(m_timewindow);
       } 
     if (n_hit>0) 
       {
        std::cout<<"n_hit"<<n_hit<<std::endl;
        sim_header->setCDLPMTtotalHits(n_hit);
       } 

     sim_header->setEvent(sim_event);
     nav->addHeader("/Event/Sim", sim_header); 
  
    if(m_split)
      { Incident::fire(*getScope(), iotaskname);}  

    LogInfo<<"tian  chong"<<std::endl;   
    
  
}

  bool junoSD_PMT_muon::domerge(int pmtid,double hittime)
{
   if(m_PMThit.empty())   
      {
         return false;
      }
  else
      {
         std::map<int, std::vector<JM::SimPMTHit*>>::iterator pmt =m_PMThit.find(pmtid);
         if (pmt != m_PMThit.end()) 
            {
              int time1 = static_cast<int>(hittime/m_time_window);
              std::vector<JM::SimPMTHit*>::iterator it = pmt->second.begin();
              for ( ; it != pmt->second.end(); ++it) 
                   {
                     if (time1 == static_cast<int>((*it)->getHitTime()/m_time_window)) 
                        {
                           if (hittime < (*it)->getHitTime()) (*it)->setHitTime(hittime);
                           (*it)->setNPE(1 + (*it)->getNPE());
                            
                           return true ;

                        }
                   }
               return false;  
            }  
            else 
            {
              return false; 
            }              
      }

}










void junoSD_PMT_muon::clear(){}

void junoSD_PMT_muon::DrawAll(){} 

void junoSD_PMT_muon::PrintAll(){} 

void junoSD_PMT_muon::SimpleHit(const ParamsForSD_PMT&){}

int junoSD_PMT_muon::get_pmtid(G4Track* track) {
    int ipmt= -1;
    // find which pmt we are in
    // The following doesn't work anymore (due to new geometry optimization?)
    //  ipmt=fastTrack.GetEnvelopePhysicalVolume()->GetMother()->GetCopyNo();
    // so we do this:
    {
        const G4VTouchable* touch= track->GetTouchable();
        int nd= touch->GetHistoryDepth();
        int id=0;
        for (id=0; id<nd; id++) {
            if (touch->GetVolume(id)==track->GetVolume()) {
                int idid=1;
                G4VPhysicalVolume* tmp_pv=NULL;
                for (idid=1; idid < (nd-id); ++idid) {
                    tmp_pv = touch->GetVolume(id+idid);

                    G4LogicalVolume* mother_vol = tmp_pv->GetLogicalVolume();
                    G4LogicalVolume* daughter_vol = touch->GetVolume(id+idid-1)->
                        GetLogicalVolume();
                    int no_daugh = mother_vol -> GetNoDaughters();
                    if (no_daugh > 1) {
                        int count = 0;
                        for (int i=0; (count<2) &&(i < no_daugh); ++i) {
                            if (daughter_vol->GetName()
                                    ==mother_vol->GetDaughter(i)->GetLogicalVolume()->GetName()) {
                                ++count;
                            }
                        }
                        if (count > 1) {
                            break;
                        }
                    }
                    // continue to find
                }
                ipmt= touch->GetReplicaNumber(id+idid-1);
                break;
            }
        }
        if (ipmt < 0) {
            G4Exception("junoPMTOpticalModel: could not find envelope -- where am I !?!", // issue
                    "", //Error Code
                    FatalException, // severity
                    "");
        }
    }

    return ipmt;
}

// ============================================================================
// = Collection Efficiency Related
// ============================================================================
// == change the Collection Efficiency Mode
void junoSD_PMT_muon::setCEMode(const std::string& mode) {
    m_ce_mode = mode;
}

// == get the Collection Efficiency 
double junoSD_PMT_muon::get_ce(const std::string& volname, const G4ThreeVector& localpos, bool pmt_type, bool qe_type) {
    // volname:
    // * PMT_20inch_body_phys
    // * PMT_3inch_body_phys
    if (m_ce_mode == "None") {
        return 1.0;
    } else if (m_ce_mode == "20inch") {
        // only 20inch PMT will be affected
         //G4cout << volname << G4endl;
        if (volname == "PMT_20inch_body_phys") {
            // calculate the angle theta
            double theta = localpos.theta();
            // do interpolate
            static double s_theta_NNVT[] = {
                0.*deg, 14.*deg, 30.*deg, 42.5*deg, 55.*deg, 67.*deg,
                77.5*deg, 85.*deg, 90.*deg,
            };
            static double s_ce_NNVT[] =    {
                0.9,    0.9,   0.845,     0.801,    0.775,    0.802,
                0.802,   0.771,    0.66,
            };
            static double s_theta_hamamatsu[] = {
                0.*deg, 13.*deg, 28.*deg, 41.*deg, 55.*deg, 66.*deg,
                79.*deg, 85.*deg, 90.*deg,
            };
            static double s_ce_hamamatsu[] =    {
                0.873,    0.873,   0.888,     0.896,    0.881,    0.9,
                0.881,     0.627,    0.262,
            };
            static G4DataInterpolation s_di(s_theta_NNVT, s_ce_NNVT, 9, 0., 0.);
            if(pmt_type){
            static G4DataInterpolation s_di(s_theta_hamamatsu, s_ce_hamamatsu, 9, 0., 0.);
            }

             return s_di.CubicSplineInterpolation(theta);
        }
/*
     
*/
        else if (volname == "HamamatsuR12860_PMT_20inch_body_phys") {
            double theta = localpos.theta();

            static double s_theta_hamamatsu[] = {
                0.*deg, 13.*deg, 28.*deg, 41.*deg, 55.*deg, 66.*deg,
                79.*deg, 85.*deg, 90.*deg,
            };
            static double s_ce_hamamatsu[] =    {
                0.911,    0.911,    0.9222,     0.9294,     0.9235,     0.93,
                0.9095, 0.6261, 0.2733, 
            };
            static G4DataInterpolation s_di(s_theta_hamamatsu, s_ce_hamamatsu, 9, 0., 0.);  
            return s_di.CubicSplineInterpolation(theta);
        }


        else if (volname == "NNVTMCPPMT_PMT_20inch_body_phys") {
            // calculate the angle theta
            double theta = localpos.theta();
            // do interpolate
            static double s_theta_NNVT[] = {
                0.*deg, 14.*deg, 30.*deg, 42.5*deg, 55.*deg, 67.*deg,
                77.5*deg, 85.*deg, 90.*deg,
            };
            static double s_ce_NNVT[] =    {
                1.0,    1.0,    0.9453,     0.9105,     0.8931,     0.9255, 
                0.9274,     0.8841,     0.734,  
            };
            static double s_ce_NNVT_highQE[] = {
               1.0,     1.0,    0.9772,     0.9723,     0.9699,     0.9697, 
               0.9452,  0.9103,     0.734,   
            };

            if(!pmt_type && !qe_type){
                static G4DataInterpolation s_di(s_theta_NNVT, s_ce_NNVT, 9, 0., 0.);

                return s_di.CubicSplineInterpolation(theta);
            }
            else if(!pmt_type && qe_type) {
                static G4DataInterpolation s_di(s_theta_NNVT, s_ce_NNVT_highQE, 9, 0., 0.);
                return s_di.CubicSplineInterpolation(theta);
            }
        }

    } else if (m_ce_mode == "20inchflat"){
        // This is a flat mode which means no matter where the photon
        // hits, use the same CE.
        if (volname == "PMT_20inch_body_phys") {
            // FIXME It's a fixed number here, we can make it a variable
            // if it is needed to be modified.
            // -- 2015.10.10 Tao Lin <lintao@ihep.ac.cn>
            static double mean_val = m_ce_flat_value;
            return mean_val;
        }
    }else if (m_ce_mode == "flat"){
        // This is a flat mode which means no matter where the photon
        // hits, use the same CE.
        // G4cout << "PMT volume name : "<<volname << G4endl;
        if (volname == "R12860TorusPMTManager_body_phys") {
            // FIXME It's a fixed number here, we can make it a variable
            // if it is needed to be modified.
            static double Ham20inch_R12860_mean_val = Ham20inch_m_ce_flat_value*Ham20inch_m_EAR_value;
            return Ham20inch_R12860_mean_val;
        }
        else if (volname == "MCP20inchPMTManager_body_phys") {
            // FIXME It's a fixed number here, we can make it a variable
            // if it is needed to be modified.
            static double MCP20inch_mean_val = MCP20inch_m_ce_flat_value*MCP20inch_m_EAR_value;
            return MCP20inch_mean_val;
        }
        else if (volname == "Ham8inchPMTManager_body_phys") {
            // FIXME It's a fixed number here, we can make it a variable
            // if it is needed to be modified.
            static double Ham8inch_mean_val = Ham8inch_m_ce_flat_value*Ham8inch_m_EAR_value;
            return Ham8inch_mean_val;
        }
        else if (volname == "MCP8inchPMTManager_body_phys") {
            // FIXME It's a fixed number here, we can make it a variable
            // if it is needed to be modified.
            static double MCP8inch_mean_val = MCP8inch_m_ce_flat_value*MCP8inch_m_EAR_value;
            return MCP8inch_mean_val;
        }
        else if (volname == "HZC9inchPMTManager_body_phys") {
            // FIXME It's a fixed number here, we can make it a variable
            // if it is needed to be modified.
            static double HZC9inch_mean_val = HZC9inch_m_ce_flat_value*HZC9inch_m_EAR_value;
            return HZC9inch_mean_val;
        }
    }else if (m_ce_mode == "20inchfunc") {
        // In this mode, user needs to input:
        // 1. a function, which can be interpret by ROOT TF1.
        // 2. a list of parameters
        if (!m_ce_func) {
            G4cout << "WARNING: The CE Function is not defined." << G4endl;
            assert(m_ce_func);
        }
        // calculate the angle theta
        double theta = localpos.theta(); // unit: radians
        if (theta>CLHEP::halfpi) { theta = CLHEP::halfpi; }
        // convert angle
        // NOTE: the angle needs to be converted
        // 1. pi/2-theta
        // 2. radian -> degree
        theta = (CLHEP::halfpi-theta)/degree;
        return m_ce_func->Eval(theta);
    } else {
        G4cout << "WARNING: unknown CE mode " << m_ce_mode << G4endl;
    }

    return 1.0;
}

void
junoSD_PMT_muon::setCEFunc(const std::string& func, const std::vector<double>& param)
{
    // detele origial function
    if (m_ce_func) {
        delete m_ce_func;
        m_ce_func = 0;
    }

    // Info:
    std::cout << "Following is the CE Function detail:" << std::endl;
    std::cout << "CE Function: " << func << std::endl;
    // angle from 0 to 90 deg.
    // angle is from equator
    m_ce_func = new TF1("ce", func.c_str(), 0, 90);
    std::cout << "CE Params: ";
    for (size_t i = 0; i < param.size(); ++i) {
        std::cout << param[i] << " "; 
        m_ce_func->SetParameter(i, param[i]);
    }
    std::cout << std::endl;
}


