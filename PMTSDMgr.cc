

#include "PMTSDMgr.hh"
#include "junoSD_PMT.hh"
#include "junoSD_PMT_v2.hh"
#include "junoSD_PMT_muon.hh"



#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"

using namespace CLHEP;
DECLARE_TOOL(PMTSDMgr);

G4ThreadLocal PMTHitMerger* PMTSDMgr::m_pmthitmerger = NULL;
G4ThreadLocal PMTHitMerger* PMTSDMgr::m_pmthitmerger_opticks = NULL;


PMTSDMgr::PMTSDMgr(const std::string& name)
    : ToolBase(name)
{
    m_merge_flag = false;
    m_time_window = 1000*ns;

    declProp("EnableMergeHit", m_merge_flag);
    declProp("MergeTimeWindow", m_time_window);
    // XXX need to switch.
    declProp("PMTSD", m_pmt_sd="junoSD_PMT_v2");
    // declProp("PMTSD", m_pmt_sd="junoSD_PMT");
    declProp("CollEffiMode", m_ce_mode="None");

    declProp("HitType", m_hit_type=1);

    declProp("CEFlatValue", m_ce_flat_value=0.9);
    // Ref to JUNO-DOC-1245
    declProp("CEFunction", m_ce_func="0.9*[0]/(1+[1]*exp(-[2]*x))");
    //NEW from the NNVT PMT test benches at Zhongshan -- GS, 2018/09/11
    m_ce_func_params.push_back(0.9194); // p0
    m_ce_func_params.push_back(0.504);  // p1
    m_ce_func_params.push_back(0.08076);// p2
    //Hamamatsu have different ones, the following: 
    //m_ce_func_params.push_back(1.02557);// p0
    //m_ce_func_params.push_back(6.77639);// p1
    //m_ce_func_params.push_back(0.16419);// p2
    
    declProp("CEFuncParams", m_ce_func_params);

    declProp("DisableSD", m_disableSD=false);
    declProp("OpticksMode", m_opticksMode=0);

    declProp("UsePMTOpticalModel", m_enable_optical_model=false);
    declProp("maxhit",m_maxhit=100);
    declProp("split", m_split=false);
}

PMTSDMgr::~PMTSDMgr() 
{
}

G4VSensitiveDetector* 
PMTSDMgr::getSD()
{
    m_pmt_param_svc = 0;
    LogInfo << "Retrieving PMTParamSvc." << std::endl;
    SniperPtr<PMTParamSvc> svc(*getParent(), "PMTParamSvc");
      if (svc.invalid()) { 
        LogError << "Can't get PMTParamSvc. We can't initialize PMT." << std::endl;
        assert(0);
      } else {
        LogInfo << "Retrieve PMTParamSvc successfully." << std::endl;
        m_pmt_param_svc = svc.data();
      } 

    m_pmt_sim_param_svc = 0;
    LogInfo << "Retrieving PMTSimParamSvc." << std::endl;
    SniperPtr<IPMTSimParamSvc> simsvc(*getParent(), "PMTSimParamSvc");
      if (svc.invalid()) { 
        LogError << "Can't get PMTSimParamSvc. We can't initialize PMT." << std::endl;
        assert(0);
      } else {
        LogInfo << "Retrieve PMTSimParamSvc successfully." << std::endl;
        m_pmt_sim_param_svc = simsvc.data();
      } 


    G4VSensitiveDetector* ifsd = 0;
    if (m_pmt_sd == "junoSD_PMT") {
        junoSD_PMT* sd = new junoSD_PMT(objName());

        sd->setMergeFlag(m_merge_flag);
        sd->setMergeWindows(m_time_window);
        ifsd = sd;
    } else if(m_pmt_sd == "junoSD_PMT_muon"){
     junoSD_PMT_muon* sd = new junoSD_PMT_muon(objName());
     sd->setScope(getParent());
     sd->setmaxhit(m_maxhit);
     sd->setsplit(m_split);

     sd->setCEMode(m_ce_mode);
     sd->setCEFlatValue(m_ce_flat_value);
     sd->setCEFunc(m_ce_func, m_ce_func_params);
     sd->setMergeFlag(m_merge_flag);
     sd->setMergeWindows(m_time_window);
     sd->setPMTParamSvc(m_pmt_param_svc);
     ifsd=sd;
     if (m_disableSD) {
            LogInfo << "junoSD_PMT_muon::ProcessHits is disabled now. " << std::endl;
            sd->disableSD();
        }   
        
} 


     else if (m_pmt_sd == "junoSD_PMT_v2") {
        junoSD_PMT_v2* sd = new junoSD_PMT_v2(objName(), m_opticksMode);
        // As a merger is attached to a specific SD, so also create new merger for the new SD.
        PMTHitMerger* pmthitmerger = new PMTHitMerger();
#ifdef WITH_G4OPTICKS
        PMTHitMerger* pmthitmerger_opticks = new PMTHitMerger();
#else
        PMTHitMerger* pmthitmerger_opticks = NULL ; 
#endif

        if (m_pmthitmerger) {
            G4cout << "WARNING: PMTSDMgr::m_pmthitmerger already exists." << G4endl;
        }
        if (m_pmthitmerger_opticks) {
            G4cout << "WARNING: PMTSDMgr::m_pmthitmerger_opticks already exists." << G4endl;
        }

        m_pmthitmerger = pmthitmerger;
        m_pmthitmerger_opticks = pmthitmerger_opticks;

        sd->setCEMode(m_ce_mode);
        // if flat mode
        sd->setCEFlatValue(m_ce_flat_value);
        // func mode
        sd->setCEFunc(m_ce_func, m_ce_func_params);
        sd->setMergeFlag(m_merge_flag);
        sd->setMergeWindows(m_time_window);
        sd->setMerger(pmthitmerger);
        sd->setPMTParamSvc(m_pmt_param_svc);
        sd->setPMTSimParamSvc(m_pmt_sim_param_svc);
        sd->setHitType(m_hit_type);
        // configure the merger
        pmthitmerger->setMergeFlag(m_merge_flag);
        pmthitmerger->setTimeWindow(m_time_window);

#ifdef WITH_G4OPTICKS
        pmthitmerger_opticks->setMergeFlag(m_merge_flag);
        pmthitmerger_opticks->setTimeWindow(m_time_window);
        sd->setMergerOpticks(pmthitmerger_opticks);
#endif

        ifsd = sd;
        if (m_disableSD) {
            LogInfo << "junoSD_PMT_v2::ProcessHits is disabled now. " << std::endl;
            sd->disableSD();
        }

        if(m_enable_optical_model){
            LogInfo << "junoSD_PMT_v2::The new PMT optical model is enabled now." << std::endl;
            sd->enableOpticalModel();
        }

    } else {
        LogError << "Unsupport PMTSD type: " << m_pmt_sd << std::endl;
        LogError << "The available types are "
                 << "junoSD_PMT, junoSD_PMT_v2"
                 << std::endl;
    }

    return ifsd;
}
