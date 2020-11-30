#ifndef PMTSDMgr_hh
#define PMTSDMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/ISensitiveDetectorElement.h"
#include "PMTHitMerger.hh"
#include <vector>
#include "Geometry/PMTParamSvc.h"
#include "IPMTSimParamSvc/IPMTSimParamSvc.h"

class PMTSDMgr: public ISensitiveDetectorElement,
                public ToolBase {
public:
    PMTSDMgr(const std::string& name);
    ~PMTSDMgr();

    // can only call once
    virtual G4VSensitiveDetector* getSD();
    // return the Hit Merger;
    PMTHitMerger* getPMTMerger() { return m_pmthitmerger; }
    PMTHitMerger* getPMTMergerOpticks() { return m_pmthitmerger_opticks ; }
    PMTParamSvc* gerPMTParSvc(){ return m_pmt_param_svc;}
    IPMTSimParamSvc* getPMTSimParamSvc() {return m_pmt_sim_param_svc; }
private:
    bool m_merge_flag;
    // the time window is used when merge is enabled.
    double m_time_window;
    std::string m_pmt_sd;
    std::string m_ce_mode;

    static G4ThreadLocal PMTHitMerger* m_pmthitmerger;
    static G4ThreadLocal PMTHitMerger* m_pmthitmerger_opticks ;
    PMTParamSvc* m_pmt_param_svc;
    IPMTSimParamSvc* m_pmt_sim_param_svc;
    G4int m_hit_type;

    // if flat mode enabled, this is used to set the fixed number
    G4double m_ce_flat_value;
    // 20inchfunc, function mode
    std::string m_ce_func;
    std::vector<double> m_ce_func_params;

    // disable the sd
    G4bool m_disableSD;
    G4int m_opticksMode ;

    // enable the new pmt optical model
    G4bool m_enable_optical_model;
    
    int m_maxhit;
    bool m_split;

};

#endif
