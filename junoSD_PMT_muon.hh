//--------------------------------------------------------------------------
//                            junoSD_PMT_v2
//
// PMTs are difined as sensitive detector. They collect hits on them.
// The data members of hits are set up here using the information of G4Step.
// -------------------------------------------------------------------------
// Author: Liang Zhan, 2006/01/27
// Modified by: Weili Zhong, 2006/03/01
// -------------------------------------------------------------------------

#ifndef junoSD_PMT_muon_h
#define junoSD_PMT_muon_h 1

#include "Event/SimHeader.h"
#include "Event/SimEvent.h"
#include "EvtNavigator/NavBuffer.h"
#include "G4Event.hh"
#include "DataRegistritionSvc/DataRegistritionSvc.h"

#include "BufferMemMgr/IDataMemMgr.h"


#include "globals.hh"

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
//#include "junoHit_PMT.hh"
//#include "junoHit_PMT_muon.hh"
#include "IToolForSD_PMT.h"
//#include "PMTHitMerger.hh"
#include <map>
#include <vector>
#include <TF1.h>
#include "Geometry/PMTParamSvc.h"
//////////////////////////////////////////////////////////////////////////

class G4Step;
class G4Track;
class Task;
//class G4HCofThisEvent;
/*
#ifdef WITH_G4OPTICKS
class PMTEfficiency ; 
class PMTEfficiencyTable ; 
#endif
*/
class G4Event;
/*namespace JM {
    class SimEvent;
}
*/


class junoSD_PMT_muon : public G4VSensitiveDetector, public IToolForSD_PMT
{
    public:
        junoSD_PMT_muon(const std::string& name);
        ~junoSD_PMT_muon();

        void Initialize(G4HCofThisEvent*HCE);
        G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
        void EndOfEvent(G4HCofThisEvent*HCE);
     //   void EndOfEvent_Opticks(G4HCofThisEvent*HCE);

        void clear();
        void DrawAll();
        void PrintAll();
        void SimpleHit( const ParamsForSD_PMT& );

        void setCEMode(const std::string& mode);
        void setCEFlatValue(double v) {m_ce_flat_value = v;}
        void setMergeFlag(bool f) { m_merge_flag = f; }
        void setMergeWindows(double t) { m_time_window = t; }
       // void setMerger(PMTHitMerger* phm) { m_pmthitmerger=phm; }
      //  void setMergerOpticks(PMTHitMerger* phm) { m_pmthitmerger_opticks=phm; }
        void setPMTParamSvc(PMTParamSvc* para){ m_PMTParamsvc=para; }
        PMTParamSvc* getPMTParamSvc() const { return m_PMTParamsvc ; }

       // void setHitType(int i) { m_hit_type = i; }
       // int getHitType() { return m_hit_type; }

        void disableSD() { m_disable = true; }
        void enableSD() { m_disable = false; }

        void setCEFunc(const std::string& func, const std::vector<double>& param);
        //split
        void setmaxhit(int m) { maxhit=m;}
        void setsplit(bool f) { m_split=f;}
  
  public:
     void setScope(Task* scope) {m_scope = scope;}
     Task* getScope() {return m_scope;}
   private:
       Task* m_scope;

       // double getEfficiencyScale() const ; 
    private:
        int get_pmtid(G4Track*);
        double get_ce(const std::string& volname, const G4ThreeVector& localpos, bool pmt_type, bool qe_type);
   /* private:
        junoHit_PMT_Collection* hitCollection;
        junoHit_PMT_muon_Collection* hitCollection_muon;
        junoHit_PMT_Collection* hitCollection_opticks ;
  */ 

   private:
        bool m_debug;
        std::string m_ce_mode;
        double m_qescale ; 
        double m_angle_response ; 

        // if flat mode enabled, this is used to set the fixed number
        double m_ce_flat_value;
        double MCP20inch_m_ce_flat_value;
        double MCP8inch_m_ce_flat_value;
        double Ham20inch_m_ce_flat_value;
        double Ham8inch_m_ce_flat_value;
        double HZC9inch_m_ce_flat_value;

        double MCP20inch_m_EAR_value;
        double MCP8inch_m_EAR_value;
        double Ham20inch_m_EAR_value;
        double Ham8inch_m_EAR_value;
        double HZC9inch_m_EAR_value;

        // 20inchfunc mode, function mode
        std::string m_ce_func_str;
        std::vector<double> m_ce_func_params;
        TF1* m_ce_func;

        // flag to enable/disable the sensitive detector
        // disable is true, means disable the SD.
        bool m_disable;

    private:
        // ========================================================================
        // merge related
        // ========================================================================
        // = merge flag 
        bool m_merge_flag;
        int  m_merge_count ; 
        // = the time window is used when merge is enabled.
        double m_time_window;

       // typedef std::multimap<int, int> PMTID2COLIDS;
       // typedef std::pair< PMTID2COLIDS::iterator, PMTID2COLIDS::iterator > PMTITER;
       // PMTID2COLIDS m_pmtid2idincol;

        // new merger
       // PMTHitMerger* m_pmthitmerger;
       // PMTHitMerger* m_pmthitmerger_opticks;
       // PMTParamSvc* m_PMTParamsvc;
    private:

       JM::SimEvent*   sim_event;
       JM::SimHeader*  sim_header;
       JM::SimPMTHit*  sim_hit;
       PMTParamSvc*  m_PMTParamsvc;   
       IDataMemMgr* m_bufmgr; 
      // double m_time_window;
      // int merge_flag;
        double m_timewindow;
        int n_hit;

        bool minmax_initialized;
        double max_CDLPMT_hittime ;
        double min_CDLPMT_hittime;

        std::map<int, std::vector<JM::SimPMTHit*>> m_PMThit;

      //  int get_pmtid(G4Track*);
        bool domerge(int pmtid ,double hittime);
      //split-out-put
        int maxhit;
        int hit_count;
        bool m_split;
        std::string iotaskname;
        Task * iotask;     


};

#endif

