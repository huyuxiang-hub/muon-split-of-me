#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao

from __future__ import print_function

import sys
import os
import Sniper
import logging
log = logging.getLogger(__name__)

import textwrap
tools_ = lambda s:filter(None,textwrap.dedent(s).split("\n"))


# print the current machine info
import platform
print("*** NODE Info:", platform.uname(), "***")
print("*** CURRENT PID: ", os.getpid(), "***")

# This function is used to supress the help message
import argparse
helpmore = False
import sys
if '--help-more' in sys.argv:
    helpmore = True
def mh(helpstr):
    if helpmore:
        return helpstr
    return argparse.SUPPRESS

def get_parser():

    import argparse

    class MakeTVAction(argparse.Action):
        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            #print "__init__ begin"
            #print option_strings 
            #print dest 
            #print nargs
            #print kwargs
            #print "__init__ end"
            super(MakeTVAction, self).__init__(option_strings, dest, nargs, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            #print "__call__ begin"
            #print parser
            #print namespace
            #print values
            # == convert a list into 3-tuple ==
            if len(values) % 3:
                print("please set %s like x1 y1 z1 x2 y2 z2 ..." %option_string)
                sys.exit(-1)
            it = iter(values)
            values = list(zip(*([it]*3)))
            setattr(namespace, self.dest, values)
            #print option_string
            #print "__call__ end"

    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()
    argparse.ArgumentParser.convert_arg_line_to_args = convert_arg_line_to_args

    parser = argparse.ArgumentParser(description='Run JUNO Detector Simulation.',
                                     fromfile_prefix_chars='@')
    parser.add_argument("--help-more", action='store_true', help="print more options")


    parser.add_argument("--loglevel", default="Info", 
                            choices=["Test", "Debug", "Info", "Warn", "Error", "Fatal"],
                            help="Set the Log Level")
    parser.add_argument("--evtmax", type=int, default=10, help='events to be processed')
    parser.add_argument("--seed", type=int, default=42, help='seed')
    parser.add_argument("--start-evtid", type=int, default=0, help='start event number.')
    parser.add_argument("--restore-seed-status", default=None, 
                                         help=mh('restore the random engine, both '
                                         'a list of integer or a file contains '
                                         'the list of integer is supported. '
                                         'such as:'
                                         '   --restore-seed-status "1,2,3..."'))
    parser.add_argument("--output", default="sample_detsim.root", help="output file name")
    parser.add_argument("--user-output", default="sample_detsim_user.root", help="output user data file name")
   
       
    parser.add_argument("--muon-maxhit",type=int,default=100,help="TBD")


    parser.add_argument("--dbtype", default="File", help="PMTSimParamSvc db type", choices=["File", "MySQL"])
    # = Split Mode =
    grp_split = parser.add_argument_group(mh("splitmode"), mh("Split mode related"))
    grp_split.add_argument("--output-split", dest="splitoutput", action="store_true", help=mh("enable split output"))
    grp_split.add_argument("--no-output-split", dest="splitoutput", action="store_false", help=mh("disable split output"))
    grp_split.add_argument("--split-maxhits", type=int, default=100, help=mh("Max hits in sub event."))
    grp_split.set_defaults(splitoutput=False)
    # split primary track step by step. For muon simulation
    grp_split.add_argument("--track-split", dest="splittrack", action="store_true", help=mh("enable split track"))
    grp_split.add_argument("--no-track-split", dest="splittrack", action="store_false", help=mh("disable split track"))
    grp_split.set_defaults(splittrack=False)
    grp_split.add_argument("--track-split-mode", default="PrimaryTrack", 
                                   choices=["PrimaryTrack",
                                            "EveryTrack",
                                            "Time"
                                            ],
                                   help=mh("Choose differet mode for track split."))
    grp_split.add_argument("--track-split-time", default=3000., type=float,
                            help=mh("Time cut for track split mode."))

    parser.add_argument("--mac", default="run.mac", help="mac file")
    parser.add_argument("--vis", default=False, action="store_true", help="start vis.")
    parser.add_argument("--vis-mac", default="vis.mac", help="visualization macro file.")

    parser.add_argument("--detoption", default="Acrylic", 
                                       choices=["Acrylic", "Balloon"],
                                       help=mh("Det Option"))
    parser.add_argument("--qescale", default=1.0, type=float, 
                                     help=mh("QE scale for ElecSim."))
    parser.add_argument("--light-yield", default=-1, type=float, 
                                     help=mh("Light Yield. If -1, default value will be used."))
    parser.add_argument("--gdml", dest="gdml", action="store_true", help="Save GDML.")
    parser.add_argument("--no-gdml", dest="gdml", action="store_false",
                                                  help="Don't Save GDML. (Default option since version J20)")
    parser.set_defaults(gdml=False)
    parser.add_argument("--dae", dest="dae", action="store_true", help=mh("Save DAE."))
    parser.add_argument("--no-dae", dest="dae", action="store_false",
                                                  help=mh("Don't Save DAE."))
    parser.set_defaults(dae=False)

    parser.add_argument("--muon-split", dest="muon_split", action="store_true", help=mh("TBD"))
    parser.add_argument("--no-muon-split", dest="muon_split", action="store_false", help=mh("TBD"))
    parser.set_defaults(muon_split=False)

    # = Calibration =
    grp_calib_unit = parser.add_argument_group(mh("calibunits"), mh("Calibration Units."))
    grp_calib_unit.add_argument("--pelletron", action="store_true",
                                       help=mh("enable pelletron in Central Detector."))

    grp_calib_unit.add_argument("--source", action="store_true",
                                    help=mh("enable source enclosure Central Detector."))
    grp_calib_unit.add_argument("--source_weights", action="store_true",
                                            help=mh("enable source enclosure and two weights in Central Detector."))
    grp_calib_unit.add_argument("--source_weight_QC", action="store_true",
                                              help=mh("enable source enclosure, bottom weight, and quick connector in Central Detector."))
    grp_calib_unit.add_argument("--submarine", action="store_true",
                                       help=mh("enable ROV source in Central Detector."))
    grp_calib_unit.add_argument("--OffsetInZ",type=float, default=0,
                                      help=mh("source assembly position in Z direction."))
    grp_calib_unit.add_argument("--OffsetInX",type=float, default=0,
                                      help=mh("source assembly position in X direction."))
    grp_calib_unit.add_argument("--OffsetInY",type=float, default=0,
                                      help=mh("source assembly position in Y direction."))
    grp_calib_unit.add_argument("--GT_source_theta",type=float, default=0,
                                      help=mh("GT source assembly position in theta direction."))
    grp_calib_unit.add_argument("--guide_tube", dest="guide_tube", action="store_true", help=mh("Add Guide Tube Into Detector"))
    grp_calib_unit.add_argument("--no-guide_tube", dest="guide_tube", action="store_false",
                                                  help=mh("Don't Add Guide Tube"))
    grp_calib_unit.set_defaults(guide_tube=True)

    # = Detector Components =
    grp_det_comp = parser.add_argument_group(mh("detcomp"), mh("Detector Components."))
    # == CD ==
    grp_det_comp.add_argument("--cd", dest="cd_enabled", action="store_true",
                                      help=mh("Enable CD."))
    grp_det_comp.add_argument("--no-cd", dest="cd_enabled", action="store_false",
                                      help=mh("Disable CD."))
    grp_det_comp.set_defaults(cd_enabled=True)
    # == WP ==
    grp_det_comp.add_argument("--wp", dest="wp_enabled", action="store_true",
                                      help=mh("Enable WP."))
    grp_det_comp.add_argument("--no-wp", dest="wp_enabled", action="store_false",
                                      help=mh("Disable WP."))
    grp_det_comp.set_defaults(wp_enabled=True)
    grp_det_comp.add_argument("--wp-pmt", dest="wp_pmt_enabled", action="store_true",
                                      help=mh("Enable PMTs in WP."))
    grp_det_comp.add_argument("--no-wp-pmt", dest="wp_pmt_enabled", action="store_false",
                                      help=mh("Disable PMTs in WP."))
    grp_det_comp.set_defaults(wp_pmt_enabled=True)
    # == TT ==
    grp_det_comp.add_argument("--tt", dest="tt_enabled", action="store_true",
                                      help=mh("Enable TT."))
    grp_det_comp.add_argument("--no-tt", dest="tt_enabled", action="store_false",
                                      help=mh("Disable TT."))
    grp_det_comp.set_defaults(tt_enabled=True)
    # == Chimney ==
    grp_det_comp.add_argument("--shutter", dest="shutter", action="store_true",
                              help=mh("Enable Shutter"))
    grp_det_comp.add_argument("--no-shutter", dest="shutter", action="store_false",
                              help=mh("Disable Shutter"))
    grp_det_comp.set_defaults(shutter=False) # disable shutter by default


    # == Commissioning ==
    grp_det_comp.add_argument("--commissioning", dest="commissioning_enabled", action="store_true",
                              help=mh("Enable commissioning"))
    grp_det_comp.add_argument("--no-commissioning", dest="commissioning_enabled", action="store_false",
                              help=mh("Disable commissioning"))
    grp_det_comp.set_defaults(commissioning_enabled=False) # disable commissioning by default
    grp_det_comp.add_argument("--below-z-is-water", type=float, default=0.0,
                              help=mh("During commissioning, below z is water. If z=0, the surface between LS and water is equator. (unit is mm)"))

    # = Material Property =
    grp_mat_prop = parser.add_argument_group(mh("matprop"), mh("Material Properties."))
    grp_mat_prop.add_argument("--mat-LS-abslen", type=float, default=77.,
                                help=mh("scale to LS.AbsLen (m) at 430nm"))
    grp_mat_prop.add_argument("--mat-LS-raylen", type=float, default=27.,
                                help=mh("scale to LS.RayleighLen (m) at 430nm"))
    # = PMT and Optical Progress related =
    grp_pmt_op = parser.add_argument_group(mh("pmtop"), mh("PMT and Optical Progress"))
    grp_pmt_op.add_argument("--pmt20inch", dest="pmt20inch", action="store_true",
                                      help=mh("Enable 20inch PMTs."))
    grp_pmt_op.add_argument("--no-pmt20inch", dest="pmt20inch", action="store_false",
                                      help=mh("Disable 20inch PMTs."))
    grp_pmt_op.set_defaults(pmt20inch=True)
    grp_pmt_op.add_argument("--pmt20inch-name", default="PMTMask",
                                      choices = ["R12860", "OnlyPMT", "20inchPMT",
                                               "R3600", "PMTMask",
                                               "HamamatsuMask", "NNVTMask"
                                               ],
                                      help=mh("20inch PMT name."))
    grp_pmt_op.add_argument("--pmt20inch-polycone-neck", default=False, action="store_true", 
                                      help=mh("Use economical polycone 20inch PMT neck shape replacing cylinder-torus."))
    grp_pmt_op.add_argument("--pmt20inch-simplify-csg", default=False, action="store_true", 
                                      help=mh("Simplify CSG modelling of 20inch PMTs, avoiding Inner_Separator anti-pattern, see HamamatsuR12860PMTManager + NNVTMCPPMTManager "))
    grp_pmt_op.add_argument("--pmt20inch-extra", default="TWO-mask",
                                      choices = ["ONE", "TWO", "TWO-mask"],
                                      help=mh("ONE category or TWO categories of LPMT. TWO: pmts without mask. TWO-mask: pmts with mask"))
    grp_pmt_op.add_argument("--pmtmask-top-thick", default=10., type=float,
                                      help=mh("PMT Mask's thickness at top"))
    grp_pmt_op.add_argument("--pmt3inch", dest="pmt3inch", action="store_true",
                                      help=mh("Enable 3inch PMTs."))
    grp_pmt_op.add_argument("--no-pmt3inch", dest="pmt3inch", action="store_false",
                                      help=mh("Disable 3inch PMTs."))
    grp_pmt_op.set_defaults(pmt3inch=True)
    grp_pmt_op.add_argument("--pmt3inch-name", default="Tub3inchV3",
                                      choices = ["Tub3inch", 
                                                 "Tub3inchV2",
                                                 "Tub3inchV3",
                                                 "Hello3inch"],
                                      help=mh("Enable 3inch PMTs."))
    grp_pmt_op.add_argument("--pmt3inch-offset", type=float, default=-50.0,
                        help=mh("The offset of the 3inch PMT (mm)."))

    # add new optical model
    
    grp_pmt_op.add_argument("--new-optical-model", dest="new_optical_model", action="store_true",
                  help=mh("Use the new optical model."))
    grp_pmt_op.add_argument("--old-optical-model", dest="new_optical_model", action="store_false",
                  help=mh("Use the old optical model"))
    grp_pmt_op.set_defaults(new_optical_model=False)
   
    



    grp_pmt_op.add_argument("--pmtsd-v2", dest="pmtsd_v2", action="store_true",
                  help=mh("Use the new PMT SD v2. (without old PMT Optical Model)"))
    grp_pmt_op.add_argument("--no-pmtsd-v2", dest="pmtsd_v2", action="store_false",
                  help=mh("Don't use the new PMT SD v2."))
    grp_pmt_op.set_defaults(pmtsd_v2=True)
    grp_pmt_op.add_argument("--ce-mode", default="20inch",
                                     choices=["None",
                                              "20inch",
                                              "20inchflat",
                                              "20inchfunc"],
                 help=mh("Different CE mode for PMTs. Only available in PMTSD-v2"))
    grp_pmt_op.add_argument("--ce-flat-value", default=0.9, type=float,
                        help=mh("Set the CE using a fixed number when 20inchflat is enabled."))
    grp_pmt_op.add_argument("--ce-func", default=None,
                        help=mh("a TF1 style string to specify CE function"))
    grp_pmt_op.add_argument("--ce-func-par", default=None, 
                        action="append", metavar="par", type=float,
                        help=mh("parameters for CE function. The first one is [0], second is [1]..."))
    grp_pmt_op.add_argument("--pmtsd-merge-twindow", type=float, default=0.0,
                        help=mh("Merge the hits in a PMT within the time window (ns)"))
    grp_pmt_op.add_argument("--optical", dest="useoptical", action="store_true",
                        help=mh("Enable Optical Progress"))
    grp_pmt_op.add_argument("--no-optical", dest="useoptical", action="store_false",
                        help=mh("Disable Optical Progress"))
    grp_pmt_op.set_defaults(useoptical=True)
    grp_pmt_op.add_argument("--cerenkov-only", dest="cerenkov_only", action="store_true",
                        help=mh("Only enable Cerenkov generation. Note: Reemission is also enabled."))
    grp_pmt_op.set_defaults(cerenkov_only=False)
    # == enable/disable cerenkov ==
    # note: cerenkov_only means disable the scintillation.
    grp_pmt_op.add_argument("--cerenkov", dest="cerenkov", action="store_true",
                        help=mh("Enable Cerenkov (default is enable)"))
    grp_pmt_op.add_argument("--no-cerenkov", dest="cerenkov", action="store_false",
                        help=mh("Disable Cerenkov (default is enable)"))
    grp_pmt_op.set_defaults(cerenkov=True)
    # == enable/disable new optical model, absreemit ==
    # note: cerenkov_only means disable the scintillation.
    '''
    grp_pmt_op.add_argument("--absreemit", dest="absreemit", action="store_true",
                        help=mh("Enable Absorption and Reemission by PPO and bis_MSB (default is disabled)"))
    grp_pmt_op.add_argument("--no-absreemit", dest="absreemit", action="store_false",
                        help=mh("Disable Absorption and Reemission by PPO and bis_MSB (default is disabled)"))
    grp_pmt_op.set_defaults(absreemit=False)
    # == use simple scint or not. Default scint also impl re-emission. Simple scint only do scint. ==
    grp_pmt_op.add_argument("--scintsimple", dest="scintsimple", action="store_true",
                        help=mh("Scint process without re-emission in it."))
    grp_pmt_op.add_argument("--no-scintsimple", dest="scintsimple", action="store_false",
                        help=mh("Disable scint process without re-emission in it. Use default one."))
    grp_pmt_op.set_defaults(scintsimple=False)
    '''
    # == track optical photons first or not ==
    grp_pmt_op.add_argument("--track-op-first", dest="track_op_first", action="store_true",
                        help=mh("Track optical photon first."))
    grp_pmt_op.add_argument("--no-track-op-first", dest="track_op_first", action="store_false",
                        help=mh("Do not track optical photon first."))
    grp_pmt_op.set_defaults(track_op_first=True)

    # == OP simulator: deferred OP simulation ==
    grp_pmt_op.add_argument("--deferred-op", dest="deferred_op", action="store_true",
                            help=mh("Enable deferred Optical simulation"))
    grp_pmt_op.add_argument("--no-deferred-op", dest="deferred_op", action="store_false",
                            help=mh("Disable deferred Optical simulation"))
    grp_pmt_op.set_defaults(deferred_op=False)

    # == opticks ==
    grp_pmt_op.add_argument("--opticks-mode", type=int, dest="opticks_mode", default=0,
                            help=mh("Control Opticks GPU Optical simulation"))
    grp_pmt_op.add_argument("--opticks-anamgr", action="store_true", dest="opticks_anamgr", default=False,
                            help=mh("Enable G4OpticksAnaMgr for extra optional Opticks analysis, requires non-standard G4OpticksBridge package."))


    # == Load optical parameters from service. ==
    grp_pmt_op.add_argument("--use-param-svc", action="store_true", default=False,
                        help=mh("Load optical parameters from param service."))
    # == enable/disable quenching ==
    grp_pmt_op.add_argument("--quenching", dest="quenching", action="store_true",
                        help=mh("Enable Quenching (default is enable)"))
    grp_pmt_op.add_argument("--no-quenching", dest="quenching", action="store_false",
                        help=mh("Disable Quenching (default is enable)"))
    grp_pmt_op.set_defaults(quenching=True)

    grp_pmt_op.add_argument("--pmt-hit-type", type=int, default=1, choices=[1,2],
                        help=mh("1 for normal hit, 2 for muon"))
    grp_pmt_op.add_argument("--pmt-disable-process", action="store_true",
                            help=mh("disable SD::ProcessHits"))
    grp_pmt_op.set_defaults(pmt_disable_process=False)

    # = physics lists tuning =
    grp_physics_list = parser.add_argument_group(mh("physicslist"), mh("physics list"))
    grp_physics_list.add_argument("--physics-list", default="JUNO", help=mh("physics list selector"))
    grp_physics_list.add_argument("--positronium", dest="positronium", action="store_true",
                                  help=mh("Enable positronium (default is enabled.)"))
    grp_physics_list.add_argument("--no-positronium", dest="positronium", action="store_false",
                                  help=mh("Disable positronium"))
    grp_physics_list.set_defaults(positronium=True)

    # = struts and fastener =
    grp_struts_fastener = parser.add_argument_group(mh("strutsfastener"), mh("struts and fastener"))
    grp_struts_fastener.add_argument("--disable-struts-fastens", 
                            default="none",
                            dest="flag_struts_fasteners",
                            choices=["all", "strut", "fastener", "none"],
                            help=mh("disable struts and fasteners"))

    # = analysis manager =
    grp_anamgr = parser.add_argument_group(mh("anamgr"), mh("analysis manager"))

    # == radioactivity related ==
    grp_anamgr.add_argument("--anamgr-grdm", action="store_true", dest="anamgr_grdm", help=mh("enable Geant4 Radioactivity Decay Module related anamgr (Default is enabled)"))
    grp_anamgr.add_argument("--no-anamgr-grdm", action="store_false", dest="anamgr_grdm", help=mh("disable Geant4 Radioactivity Decay Module related anamgr (Default is enabled)"))
    grp_anamgr.set_defaults(anamgr_grdm=True)
    grp_anamgr.add_argument("--stop-at-Pa234m", action="store_true", dest="stopAtPa234m", help=mh('During Geant4 Radioactivity Decay simulation, force stop at Pa234m and kill all secondaries.'))
    grp_anamgr.add_argument("--no-stop-at-Pa234m", action="store_false", dest="stopAtPa234m", help=mh('During Geant4 Radioactivity Decay simulation, force stop at Pa234m and kill all secondaries.'))
    grp_anamgr.set_defaults(stopAtPa234m=True)

    # == event data model ==
    grp_anamgr.add_argument("--anamgr-edm", action="store_true", dest="anamgr_edm", help=mh("enable event data model writer, including the writer with split"))
    grp_anamgr.add_argument("--no-anamgr-edm", action="store_false", dest="anamgr_edm", help=mh("disable event data model writer, including the writer with split"))
    grp_anamgr.set_defaults(anamgr_edm=True)

    # == tt ==
    grp_anamgr.add_argument("--anamgr-tt", action="store_true", dest="anamgr_tt", help=mh("enable TT output"))
    grp_anamgr.add_argument("--no-anamgr-tt", action="store_false", dest="anamgr_tt", help=mh("disable TT output"))
    grp_anamgr.set_defaults(anamgr_tt=False)

    # == normal ==
    grp_anamgr.add_argument("--anamgr-normal", action="store_true", dest="anamgr_normal", help=mh("TBD"))
    grp_anamgr.add_argument("--no-anamgr-normal", action="store_false", dest="anamgr_normal", help=mh("TBD"))
    grp_anamgr.set_defaults(anamgr_normal=True)
   
    grp_anamgr.add_argument("--anamgr-normal-hit",action="store_true",dest="anamgr_normal_hit",help=mh("TBD"))
    grp_anamgr.add_argument("--no-anamgr-normal-hit",action="store_false",dest="anamgr_normal_hit",help=mh("TBD"))
    grp_anamgr.set_defaults(anamgr_normal_hit=False) 
   

    # == genevt ==
    grp_anamgr.add_argument("--anamgr-genevt", action="store_true", dest="anamgr_genevt", help=mh("TBD"))
    grp_anamgr.add_argument("--no-anamgr-genevt", action="store_false", dest="anamgr_genevt", help=mh("TBD"))
    grp_anamgr.set_defaults(anamgr_genevt=True)
    # == deposit ==
    grp_anamgr.add_argument("--anamgr-deposit", action="store_true", dest="anamgr_deposit", help=mh("TBD"))
    grp_anamgr.add_argument("--no-anamgr-deposit", action="store_false", dest="anamgr_deposit", help=mh("TBD"))
    grp_anamgr.set_defaults(anamgr_deposit=True)
    grp_anamgr.add_argument("--anamgr-deposit-tt", action="store_true", dest="anamgr_deposit_tt", help=mh("Enable TT output"))
    grp_anamgr.add_argument("--no-anamgr-deposit-tt", action="store_false", dest="anamgr_deposit_tt", help=mh("Disable TT output"))
    grp_anamgr.set_defaults(anamgr_deposit_tt=True)
    # == interesting process ==
    grp_anamgr.add_argument("--anamgr-interesting-process", action="store_true", dest="anamgr_interesting_process", help=mh("TBD"))
    grp_anamgr.add_argument("--no-anamgr-interesting-process", action="store_false", dest="anamgr_interesting_process", help=mh("TBD"))
    grp_anamgr.set_defaults(anamgr_interesting_process=True)
    # == optical parameter ==
    grp_anamgr.add_argument("--anamgr-optical-parameter", action="store_true", dest="anamgr_optical_parameter", help=mh("TBD"))
    grp_anamgr.add_argument("--no-anamgr-optical-parameter", action="store_false", dest="anamgr_optical_parameter", help=mh("TBD"))
    grp_anamgr.set_defaults(anamgr_optical_parameter=True)
    # == timer ==
    grp_anamgr.add_argument("--anamgr-timer", action="store_true", dest="anamgr_timer", help=mh("TBD"))
    grp_anamgr.add_argument("--no-anamgr-timer", action="store_false", dest="anamgr_timer", help=mh("TBD"))
    grp_anamgr.set_defaults(anamgr_timer=True)
    # == photon tracking ==
    grp_anamgr.add_argument("--anamgr-photon-tracking", action="store_true", dest="anamgr_photon_tracking", help=mh("TBD"))
    grp_anamgr.add_argument("--no-anamgr-photon-tracking", action="store_false", dest="anamgr_photon_tracking", help=mh("TBD"))
    grp_anamgr.set_defaults(anamgr_photon_tracking=False)
    # == extend the anamgr ==
    grp_anamgr.add_argument("--anamgr-list", action="append", metavar="anamgr", default=[],
            help=mh("append anamgr to the anamgr list. You can specify anamgrs multiple times. "
                "such as: "
                " \"--anamgr-list anamgr1 --anamgr-list anamgr2\", so that "
                "both anamgr1 and anamgr2 are added to the list. "
                ))
    grp_anamgr.add_argument("--anamgr-config-file", 
            help=mh("configure the anamgr from file."))
    # = Voxel Method =
    # == enable/disable voxel method ==
    grp_voxel = parser.add_argument_group(mh("voxel"), mh("Voxel Method Related"))
    grp_voxel.add_argument("--voxel-fast-sim", dest="voxel_fast_sim", action="store_true", help=mh("TBD"))
    grp_voxel.add_argument("--no-voxel-fast-sim", dest="voxel_fast_sim", action="store_false", help=mh("TBD"))
    grp_voxel.set_defaults(voxel_fast_sim=False)
    # == enable/disable merge mode and time window ==
    grp_voxel.add_argument("--voxel-merge-flag", action="store_true", dest="voxel_merge_flag", help=mh("TBD"))
    grp_voxel.add_argument("--voxel-no-merge-flag", action="store_false", dest="voxel_merge_flag", help=mh("TBD"))
    grp_voxel.set_defaults(voxel_merge_flag=False)
    grp_voxel.add_argument("--voxel-merge-twin", default=1, type=float, help=mh("TBD"))
    # == debug mode: fill ntuple ==
    grp_voxel.add_argument("--voxel-fill-ntuple", action="store_true", dest="voxel_fill_ntuple", help=mh("TBD"))
    grp_voxel.add_argument("--voxel-no-fill-ntuple", action="store_false", dest="voxel_fill_ntuple", help=mh("TBD"))
    grp_voxel.set_defaults(voxel_fill_ntuple=False)
    grp_voxel.add_argument("--voxel-fast-dir", help=mh("Stored data for fast simulation."))
    # === gen npe ===
    grp_voxel.add_argument("--voxel-gen-npe-on", action="store_true", dest="voxel_gen_npe", help=mh("TBD"))
    grp_voxel.add_argument("--voxel-gen-npe-off", action="store_false", dest="voxel_gen_npe", help=mh("TBD"))
    grp_voxel.set_defaults(voxel_gen_npe=True)
    # === gen time ===
    grp_voxel.add_argument("--voxel-gen-time-on", action="store_true", dest="voxel_gen_time", help=mh("TBD"))
    grp_voxel.add_argument("--voxel-gen-time-off", action="store_false", dest="voxel_gen_time", help=mh("TBD"))
    grp_voxel.set_defaults(voxel_gen_time=True)
    # === save hits ===
    grp_voxel.add_argument("--voxel-save-hits-on", action="store_true", dest="voxel_save_hits", help=mh("TBD"))
    grp_voxel.add_argument("--voxel-save-hits-off", action="store_false", dest="voxel_save_hits", help=mh("TBD"))
    grp_voxel.set_defaults(voxel_save_hits=True)
    # == no PMTs and Structs ==
    grp_voxel.add_argument("--voxel-pmts-structs", action="store_true", dest="voxel_pmts_structs", help=mh("TBD"))
    grp_voxel.add_argument("--voxel-no-pmts-structs", action="store_false", dest="voxel_pmts_structs", help=mh("TBD"))
    grp_voxel.set_defaults(voxel_pmts_structs=True)
    grp_voxel.add_argument("--voxel-quenching-scale", type=float, default=0.93,
                           help=mh("Quenching factor, Qedep->edep. gamma 0.93, e- 0.98"))

    # = global time =
    grp_globaltime = parser.add_argument_group(mh("globaltime"), mh("Global time related"))
    grp_globaltime.add_argument("--global-time-begin", 
        default="1970-01-01 00:00:01",
        help=mh("Global time begin"))
    grp_globaltime.add_argument("--global-time-end",
        default="2038-01-19 03:14:07",
        help=mh("Global time end"))
    grp_globaltime.add_argument("--global-event-rate",
        default=0.0, type=float,
        help=mh("Event rate. if greater than 0, global time mode is enabled."))

    # = For different generator, we try to reuse some common part
    # == used by GtPositionerTool
    base_parser_positioner = argparse.ArgumentParser('positioner', add_help=False)
    base_parser_positioner.add_argument("--material", default="None", help="material")
    base_parser_positioner.add_argument("--volume", default="None", 
                                     choices=list(DATA_MATERIALS.keys()),
                                     help="Volume name")
    base_parser_positioner.add_argument("--volume-radius-min", default=0.0, type=float,
                                     help="min of the radius")
    base_parser_positioner.add_argument("--volume-radius-max", default=0.0, type=float,
                                     help="max of the radius")
    # z cut
    base_parser_positioner.add_argument("--volume-z-min", default=None, type=float,
                                     required='--volume-z-max' in sys.argv,
                                     help="min Z")
    base_parser_positioner.add_argument("--volume-z-max", default=None, type=float,
                                     required='--volume-z-min' in sys.argv,
                                     help="max Z")
    # x cut
    base_parser_positioner.add_argument("--volume-x-min", default=None, type=float,
                                     required='--volume-x-max' in sys.argv,
                                     help="min X")
    base_parser_positioner.add_argument("--volume-x-max", default=None, type=float,
                                     required='--volume-x-min' in sys.argv,
                                     help="max X")
    # y cut
    base_parser_positioner.add_argument("--volume-y-min", default=None, type=float,
                                     required='--volume-y-max' in sys.argv,
                                     help="min Y")
    base_parser_positioner.add_argument("--volume-y-max", default=None, type=float,
                                     required='--volume-y-min' in sys.argv,
                                     help="max Y")
    base_parser_positioner.add_argument("--global-position", default=None,
                                     nargs='+', type=float, action=MakeTVAction,
                                     help="Global Postion. It will omit the volume and material")

    # = gentool =
    subparsers = parser.add_subparsers(help='Please select the generator mode', 
                                       dest='gentool_mode')
    # == gun mode ==
    parser_gun = subparsers.add_parser("gun", help="gun mode", parents=[base_parser_positioner])
    parser_gun.add_argument("--particles",default="gamma", nargs='+',
                            help="Particles to do the simulation.")
    parser_gun.add_argument("--momentums",default=1.0, nargs='+',
                            type=float, 
                            help="Momentums(MeV) p1 p2 ....")
    parser_gun.add_argument("--momentums-mode", default="Fix",
                                choices=["Fix", "Uniform", "Range", "Gaus"],
                                help="different momentum modes")
    parser_gun.add_argument("--momentums-extra-params", nargs='+',
                            type=float, 
                            help="Extra Momentums Parameters(MeV) p1 p2 .... when mode is different, meaning is different."
                                 " Uniform: [mom-param, mom+param];"
                                 " Range: [mom, param];"
                                 " Gaus: Gaus(mom, param);"
                            )
    parser_gun.add_argument("--momentums-interp", default="Momentum",
                                choices=["Momentum", "KineticEnergy", "TotalEnergy"],
                                help="Interpret momentum.")
    parser_gun.add_argument("--positions",default=[(0,0,0)], nargs='+',
                            type=float, action=MakeTVAction,
                            help="Positions(mm) x1 y1 z1 x2 y2 z2 ....")
    parser_gun.add_argument("--directions",default=None, nargs='+',
                            type=float, action=MakeTVAction,
                            help="If you don't set, the directions are randoms. "
                                 "Directions dx1 dy1 dz1 dx2 dy2 dz2 ....")
    # == optical photon mode ==
    parser_photon = subparsers.add_parser("photon", help="optical photon mode", parents=[base_parser_positioner])
    parser_photon.add_argument("--totalphotons", type=int, default=11522, help="total generated photons")
    parser_photon.add_argument("--cos-theta-lower", type=float, default=-1.)
    parser_photon.add_argument("--cos-theta-upper", type=float, default=+1.)
    parser_photon.add_argument("--fixed-energy", type=float, help="If specified, the energy will be fixed.")
    # == gendecay mode ==
    parser_gendecay = subparsers.add_parser("gendecay", help="GenDecay mode", parents=[base_parser_positioner])
    parser_gendecay.add_argument("--nuclear", default="U-238", help="mother nuclide name")
    parser_gendecay.add_argument("--stop-nuclear", default="", help="Stop decay when reach the nuclide.")
    parser_gendecay.add_argument("-t", "--correlation-time", default=None, type=float, help="correlation time (ns).")
    parser_gendecay.add_argument("-d", "--decay-depth", default=-1, type=int, help="decay depth")
    # == hepevt mode ==
    parser_hepevt = subparsers.add_parser("hepevt", help="HepEvt mode", parents=[base_parser_positioner])
    parser_hepevt.add_argument("--exe", default="IBD", 
                                         choices=list(GENERATOR_EXEC.keys()),
                                         help="select the Generator to run")
    parser_hepevt.add_argument("--file", default=None,
                                         help="specify the HepEvt filename.")
    # == pelletron beam ==
    parser_beam = subparsers.add_parser("beam", help="Pelletron Beam mode")
    parser_beam.add_argument("--particle", default="e+", help="Particle Name")
    # === position of plane===
    parser_beam.add_argument("--plane-r", default=10., type=float,
                                          help="Plane Radius (mm)")
    parser_beam.add_argument("--plane-x", default=0, type=float,
                                          help="Plane position X (mm)")
    parser_beam.add_argument("--plane-y", default=0, type=float,
                                          help="Plane position Y (mm)")
    parser_beam.add_argument("--plane-z", default=1e3, type=float,
                                          help="Plane position Z (mm)")
    # === direction of plane ===
    parser_beam.add_argument("--plane-dirx", default=0, type=float,
                                          help="Plane direction X (global coord)")
    parser_beam.add_argument("--plane-diry", default=0, type=float,
                                          help="Plane direction Y (global coord)")
    parser_beam.add_argument("--plane-dirz", default=-1, type=float,
                                          help="Plane direction Z (global coord)")
    # === beam momentum ===
    parser_beam.add_argument("--momentum", default=1., type=float,
                                          help="Momentum (MeV)")
    parser_beam.add_argument("--momentum-spread", default=1.e-2, type=float,
                                          help="Momentum Spread (MeV)")
    parser_beam.add_argument("--divergence", default=0.10, type=float,
                                          help="Beam divergence (deg)")
    # == supernova mode ==
    parser_sn = subparsers.add_parser("sn", help="supernova mode", parents=[base_parser_positioner])
    parser_sn.add_argument("--input", help="supernova input file")
    parser_sn.add_argument("--index", type=int, default=0, help='supernova start index')
    parser_sn.add_argument("--relative-hittime", help="Use relative hit time in each event, while save the event time in event navigator. (default)",
                           dest="relative_hittime", action="store_true")
    parser_sn.add_argument("--absolute-hittime", help="Use absolute hit time in each event",
                           dest="relative_hittime", action="store_false")
    parser_sn.set_defaults(relative_hittime=True)

    # Solar Neutrino
    parser_nusol = subparsers.add_parser("nusol", help="solar neutrino mode", parents=[base_parser_positioner])
    parser_nusol.add_argument("--type", default="B8",help="neutrino type", choices = ["pp", "Be7", "Be7_862", "Be7_384", "B8", "N13", "O15", "F17", "pep", "hep"]) 

    # == atmospheric mode ==
    parser_atmo = subparsers.add_parser("atmo", help="atmospheric mode", parents=[base_parser_positioner])
    parser_atmo.add_argument("--input", help="atmospheric input file, you can find the default data: $JUNOTOP/data/Generator/NuAtm/data/tree_100000100.root")
    parser_atmo.add_argument("--index", type=int, default=0, help='atmospheric start index')

    # == spallation neutron mode ==
    parser_neutron = subparsers.add_parser("neutron", help="spallation neutron mode")
    parser_neutron.add_argument("--input", help="spallation neutron input file")
    parser_neutron.add_argument("--index", type=int, default=0, help='start index')
    parser_neutron.add_argument("--energy", type=float, default=None, help='Ek (MeV)')

    return parser

# helper for GtPositionerTool.
# 
def helper_positioner(gt, default="pTarget"):
    if args.volume == "None" and default is None:
        return
    elif args.volume == "None":
        args.volume = default

    print("SETUP POSITIONER")

    # = enable the gen in volume mode =
    # == positioner related ==
    gun_pos = gt.createTool("GtPositionerTool")
    gun_pos.property("GenInVolume").set(args.volume)
    if args.material == "None":
        gun_pos.property("Material").set(DATA_MATERIALS[args.volume])
    else:
        gun_pos.property("Material").set(args.material)
    # === volume cut ===
    radius_vec = []
    if args.volume_radius_min != 0.0:
        radius_vec.append(args.volume_radius_min)
    if args.volume_radius_max != 0.0:
        radius_vec.append(args.volume_radius_max)
    gun_pos.property("RadiusCut").set(radius_vec)
    # === Z cut ===
    z_vec = []
    if args.volume_z_min:
        z_vec.append(args.volume_z_min)
    if args.volume_z_max:
        z_vec.append(args.volume_z_max)
    gun_pos.property("ZCut").set(z_vec)
    # === X cut ===
    x_vec = []
    if args.volume_x_min:
        x_vec.append(args.volume_x_min)
    if args.volume_x_max:
        x_vec.append(args.volume_x_max)
    gun_pos.property("XCut").set(x_vec)
    # === Y cut ===
    y_vec = []
    if args.volume_y_min:
        y_vec.append(args.volume_y_min)
    if args.volume_y_max:
        y_vec.append(args.volume_y_max)
    gun_pos.property("YCut").set(y_vec)


    # == global positions ==
    if args.global_position:
        if len(args.global_position) != 1:
            assert(len(args.global_position) != 1)
        gun_pos.property("PositionMode").set("GenInGlobal")
        gun_pos.property("Positions").set(args.global_position[0])
    # = append it =
    gt.property("GenToolNames").append(gun_pos.objName())

def setup_generator(task):
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")

    gun = gt.createTool("GtGunGenTool/gun")
    gun.property("particleNames").set(args.particles)
    gun.property("particleMomentums").set(args.momentums)
    gun.property("particleMomentumMode").set(args.momentums_mode)
    if args.momentums_extra_params:
        gun.property("particleMomentumParams").set(args.momentums_extra_params)
    gun.property("particleMomentumInterp").set(args.momentums_interp)
    if args.directions:
        gun.property("DirectionMode").set("Fix")
        gun.property("Directions").set([makeTV(px,py,pz) for px,py,pz in args.directions])
    print(args.positions)
    if len(args.positions) == 1:
        gun.property("PositionMode").set("FixOne")
    else:
        gun.property("PositionMode").set("FixMany")
    gun.property("Positions").set([makeTV(x,y,z) for x,y,z in args.positions])

    gt.property("GenToolNames").set([gun.objName()])
    
    # positioner
    # Note: by default, we choose fixed position instead of using positioner
    helper_positioner(gt, None)


def setup_generator_photon(task):
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")
    # optical photon gun (using LS emission spectrum)
    gun = gt.createTool("GtOpScintTool/gun")
    gun.property("PhotonsPerEvent").set(args.totalphotons)
    gun.property("cosThetaLower").set(args.cos_theta_lower)
    gun.property("cosThetaUpper").set(args.cos_theta_upper)

    # Some user need to fix the wavelength
    if args.fixed_energy:
        gun.property("EnergyMode").set("Fixed")
        gun.property("FixedEnergy").set(args.fixed_energy)
        gun.property("TimeMode").set("Fixed")

    gt.property("GenToolNames").set([gun.objName()])
    # positioner
    helper_positioner(gt)

def setup_generator_gendecay(task):
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")
    # == gendecay related ==
    Sniper.loadDll("libGenDecay.so")
    era = gt.createTool("GtDecayerator")
    era.property("ParentNuclide").set(args.nuclear)
    era.property("StopNuclide").set(args.stop_nuclear)
    correlation_time = args.correlation_time
    if correlation_time is None and args.nuclear in DECAY_DATA:
        correlation_time = DECAY_DATA[args.nuclear]
    # if correlation_time is still None, assign default value 1s
    if correlation_time is None:
        correlation_time = 1e9
    era.property("CorrelationTime").set(correlation_time)
    era.property("ParentAbundance").set(5e16)
    era.property("DecayDepth").set(args.decay_depth)
    gt.property("GenToolNames").set([era.objName()])
    # == positioner related ==
    helper_positioner(gt)
    # == GtTimeOffsetTool ==
    toffset = gt.createTool("GtTimeOffsetTool")
    gt.property("GenToolNames").append(toffset.objName())

def setup_generator_hepevt(task):
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")
    # == HepEvt to HepMC ==
    gun = gt.createTool("GtHepEvtGenTool/gun")
    #gun.property("Source").set("K40.exe -seed 42 -n 100|")
    source = GENERATOR_EXEC[args.exe].format(SEED=args.seed,
                                            EVENT=args.evtmax)
    if args.file:
        source = args.file
    gun.property("Source").set(source)
    gt.property("GenToolNames").set([gun.objName()])
    # == positioner related ==
    # === if muon event, use the hepevt file's position ===
    if args.exe == "Muon":
        pass
    else:
        helper_positioner(gt)
    # == GtTimeOffsetTool ==
    toffset = gt.createTool("GtTimeOffsetTool")
    gt.property("GenToolNames").append([toffset.objName()])

def setup_generator_beam(task):
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")

    from GenTools import makeTV
    gun = gt.createTool("GtPelletronBeamerTool/gun")
    gun.property("particleName").set(args.particle)
    gun.property("planeCentrePos").set(makeTV(args.plane_x,
                                              args.plane_y,
                                              args.plane_z)) # (0,0,1m)
    gun.property("planeDirection").set(makeTV(args.plane_dirx,
                                              args.plane_diry,
                                              args.plane_dirz)) # down
    gun.property("planeRadius").set(args.plane_r) # 20mm
    import math
    gun.property("beamThetaMax").set(math.radians(args.divergence)) # 10deg -> rad
    gun.property("beamMomentum").set(args.momentum) # 1MeV
    gun.property("beamMomentumSpread").set(args.momentum_spread) # 0.1MeV

    gt.property("GenToolNames").set([gun.objName()])

# == additional Calib Unit ==
def setup_calib_unit(acrylic_conf, enable):

    # old one: Calib_GuideTube
    guide_tube_option = "Calib_GuideTube_V1"

    detsim0 = acrylic_conf.detsimfactory()
    detsimalg = acrylic_conf.detsimalg()
    Sniper.loadDll("libCalibUnit.so")
    detsim0.property("CalibUnitEnable").set(enable)
    if not enable:
        print("setup_calib_unit exit as not enabled") 
        return
    pass 

    detsim0.property("CalibUnitName").set(guide_tube_option) # this is default one
    detsim0.property("CalibUnitExtras").set([]) # Enable more calibration units here

    if guide_tube_option == "Calib_GuideTube_V1":
        guide_tube_v1_0 = detsimalg.createTool("Calib_GuideTube_Construction/Calib_GuideTube_Construction_V1_0")
        guide_tube_v1_0.property("Option").set("V1_0")

        if args.GT_source_theta > 0:
            guide_tube_v1_0.property("Theta").set(args.GT_source_theta)
            guide_tube_v1_0.property("UseSource").set(True)

        guide_tube_v1_1 = detsimalg.createTool("Calib_GuideTube_Construction/Calib_GuideTube_Construction_V1_1")
        guide_tube_v1_1.property("Option").set("V1_1")

        if args.GT_source_theta < 0:
            guide_tube_v1_1.property("Theta").set(args.GT_source_theta)
            guide_tube_v1_1.property("UseSource").set(True)

        import math

        # the default guide tube is constructed at -x-z plane. 
        # V1_0: 123.40 deg. Rotate from 180 to 123.40, hence 123.40-180
        # V1_0: 267.40 deg. Rotate from 180 to 267.40, hence 267.40-180

        guide_tube_place_v1_0 = detsimalg.createTool("Calib_GuideTube_Placement/Calib_GuideTube_Placement_V1_0")
        # guide_tube_place_v1_0.property("Phi").set(-(90.+2.6)*math.pi/180.)
        guide_tube_place_v1_0.property("Phi").set((123.40-180.0)*math.pi/180.)

        guide_tube_place_v1_1 = detsimalg.createTool("Calib_GuideTube_Placement/Calib_GuideTube_Placement_V1_1")
        # guide_tube_place_v1_1.property("Phi").set((303.4-180.)*math.pi/180.)
        guide_tube_place_v1_1.property("Phi").set((267.40-180.)*math.pi/180.)

    elif guide_tube_option == "Calib_GuideTube" and args.GT_source_theta:
        guide_tube_source_place = detsimalg.createTool("Calib_GuideTube_Construction")
        guide_tube_source_place.property("Theta").set(args.GT_source_theta);
        guide_tube_source_place.property("UseSource").set(True);

    # == Add calib source geometry into detector if necessary
    if args.source or args.source_weight_QC or args.source_weights or args.submarine: 
     detsim0.property("CalibUnitExtras").set(["lSourceWorld"]) #
    
     Calib_Source = detsimalg.createTool("GDMLDetElemConstruction/lSourceWorld")
    
     path = os.environ.get("CALIBUNITROOT")
     if args.source: 
      Calib_Source.property("GdmlFilename").set(path+"/share/source.gdml") # for only source enclosure
     if args.source_weight_QC: 
      Calib_Source.property("GdmlFilename").set(path+"/share/source_weight_QuickConnector.gdml")#include source enclosure, bottom weight and top quick connector
     if args.source_weights: 
      Calib_Source.property("GdmlFilename").set(path+"/share/source_weights.gdml") # source enclosure and two weights
     if args.submarine: 
      Calib_Source.property("GdmlFilename").set(path+"/share/submarine.gdml") # this is for ROV calibration system
    
     calibsourceplace = detsimalg.createTool("CalibSourcePlacement")
     calibsourceplace.property("OffsetInZ").set(args.OffsetInZ)
     calibsourceplace.property("OffsetInY").set(args.OffsetInY)
     calibsourceplace.property("OffsetInX").set(args.OffsetInX)


# == additional Calib Unit ==
def setup_calib_pelletron(acrylic_conf, enable):
    detsim0 = acrylic_conf.detsimfactory()
    detsimalg = acrylic_conf.detsimalg()
    Sniper.loadDll("libCalibUnit.so")
    detsim0.property("CalibUnitEnable").set(enable)
    if not enable:
        return
    pass

    detsim0.property("CalibUnitName").set("CalibTube")
    # Calib Unit Related
    calibtube = detsimalg.createTool("CalibTubeConstruction")
    print(calibtube)
    calibtubeplace = detsimalg.createTool("CalibTubePlacement")
    # FIXME a more general geometry service is needed.
    calibTubeLength1 = 17.3e3; # 17.3m
    calibTubeLength2 = 0.3e3   #  0.3m
    offset_z_in_cd = (calibTubeLength1+calibTubeLength2)/2.
    calibtubeplace.property("OffsetInZ").set(offset_z_in_cd)

    acrylic_conf.add_anamgr("DepositEnergyCalibAnaMgr")
    calib_anamgr = detsimalg.createTool("DepositEnergyCalibAnaMgr")
    calib_anamgr.property("EnableNtuple").set(True)

# = setup supernova =
def setup_generator_sn(task):
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")
    # supernova
    gun = gt.createTool("GtSNTool/gun")
    # check the input file
    import os.path
    if not args.input or not os.path.exists(args.input):
        print("can't find the supernova input file '%s'"%args.input)
        sys.exit(-1)
    gun.property("inputSNFile").set(args.input)
    gun.property("StartIndex").set(args.index)

    gt.property("GenToolNames").set([gun.objName()])
    # positioner
    helper_positioner(gt)
    # time
    gt.property("EnableSNTime").set(args.relative_hittime)

# = setup solar neutrinos=
def setup_generator_nusol(task):
    import GenTools
    import NuSolGen
    from GenTools import makeTV
    gt = task.createAlg("GenTools")
    # nusol
    gun = gt.createTool("GtNuSolTool/gun")
    gun.property("neutrinoType").set(args.type);

    gt.property("GenToolNames").set([gun.objName()])
    # positioner
    helper_positioner(gt)

# = setup atmospheric neutrinos=
def setup_generator_atmo(task):
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")
    #gst is genie format of atmospheric neutrino generator
    gun_gst = gt.createTool("GtGstTool")
    # check the input file
    import os.path
    if not args.input or not os.path.exists(args.input):
        print("can't find atmospheric input file '%s'"%args.input)
        sys.exit(-1)

    gun_gst.property("inputGstFile").set(args.input)
    gun_gst.property("GstStartIndex").set(args.index)
    gt.property("GenToolNames").set([gun_gst.objName()])
    # positioner
    helper_positioner(gt)

def setup_generator_neutron(task):
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")
    gun = gt.createTool("GtNeutronTool")
    # check the input file
    import os.path
    if not args.input or not os.path.exists(args.input):
        print("can't find neutron input file '%s'"%args.input)
        sys.exit(-1)

    gun.property("inputFile").set(args.input)
    gun.property("startIndex").set(args.index)
    if args.energy:
        gun.property("neutronEnergy").set(args.energy)
    gt.property("GenToolNames").set([gun.objName()])

##############################################################################
DEFAULT_GDML_OUTPUT = {"Acrylic": "geometry_acrylic.gdml", 
                       "Balloon": "geometry_balloon.gdml"}
DEFAULT_DAE_OUTPUT = {"Acrylic": "geometry_acrylic.dae", 
                       "Balloon": "geometry_balloon.dae"}

# Note: please register the Physical Volume name. 
#       You could get the list of PV names when running simulation with volumes mode.
DATA_MATERIALS = {"PMT_20inch_body_phys": "Pyrex",
                  "PMT_3inch_body_phys": "Pyrex",
                  "pCentralDetector": "Steel",
                  "pTarget": "LS",
                  "pSource": "Analcime",
                  "pTopRock": "Rock",
                  "pBtmRock": "Rock",
                  "pBar": "Scintillator",
                  # for TWO-mask mode. note: PMT_20inch_body_phys is for WP.
                  "pLPMT_Hamamatsu_R12860": "Pyrex",
                  "pLPMT_NNVT_MCPPMT": "Pyrex",
                  }
DECAY_DATA = {"U-238": 1.5e5, "U238": 1.5e5, # alias
              "Th-232": 280, "Th232": 280,
              "K-40": 1e9, "K40": 1e9,
              "Co-60": 1e9, "Co60": 1e9,
              } # unit: ns

GENERATOR_EXEC = {"IBD": "IBD.exe -n {EVENT} -seed {SEED}|",
                  "IBD-NH": "IBD.exe -n {EVENT} -seed {SEED} -NH|",
                  "IBD-IH": "IBD.exe -n {EVENT} -seed {SEED} -IH|",
                  "IBD-eplus": "IBD.exe -n {EVENT} -seed {SEED} -eplus_only|",
                  "IBD-neutron": "IBD.exe -n {EVENT} -seed {SEED} -neutron_only|",
                  "AmC": "AmC.exe -n {EVENT} -seed {SEED}|",
                  "Muon": "Muon.exe -n {EVENT} -seed {SEED} -s juno|",
                  "Co60": "Co60.exe -n {EVENT} -seed {SEED}|",
                  "Cs137": "Cs137.exe -n {EVENT} -seed {SEED}|",
                  "Ge68": "Ge68.exe -n {EVENT} -seed {SEED}|",
                  "Ge68-geom": "Ge68.exe -n {EVENT} -seed {SEED} -geom 1|",
                  "DSNB-NC": "DSNB-NC.exe -n {EVENT} -seed {SEED} |",
                  }
DATA_LOG_MAP = {
        "Test":0, "Debug":2, "Info":3, "Warn":4, "Error":5, "Fatal":6
        }

if __name__ == "__main__":
    parser = get_parser()
    import sys
    if '--help-more' in sys.argv:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    pylogfmt = '[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s'
    pyloglevel = "NOTSET" if args.loglevel == 'Test' else args.loglevel.upper()
    logging.basicConfig(level=getattr(logging,pyloglevel), format=pylogfmt)

    log.info(args)
    #import sys; sys.exit(0)
    gdml_filename = None
    dae_filename = None
    if args.gdml:
        gdml_filename = DEFAULT_GDML_OUTPUT[args.detoption]
    if args.dae:
        dae_filename = DEFAULT_DAE_OUTPUT[args.detoption]

    # initial seed check
    if args.seed == 0:
        args.seed = 0x19900418 # MAGIC NUMBER
    elif args.seed < 0:
        args.seed = -args.seed

    print("= INITIALIZATION =")
    task = Sniper.TopTask("detsimtask")
    # task.asTop()
    task.setEvtMax(args.evtmax)
    task.setLogLevel(DATA_LOG_MAP[args.loglevel])
    # = I/O Related =
    print("== Data Registrition Svc ==")
    import DataRegistritionSvc
    log.info("task.createSvc DataRegistritionSvc")
    task.createSvc("DataRegistritionSvc")
    
    # if split output, we need to create another iotask for output
    iotask = task
    
    if args.muon_split:
        iotask=task.createTask("Task/muon_io")
        import  DataRegistritionSvc
        dr = iotask.createSvc("DataRegistritionSvc")
        import BufferMemMgr
        bufMgr = iotask.createSvc("BufferMemMgr")
   
    if args.splitoutput:
        iotask = task.createTask("Task/detsimiotask")
        import DataRegistritionSvc
        dr = iotask.createSvc("DataRegistritionSvc")
        #dr.property("EventToPath").set({"JM::SimEvent": "/Event/Sim"})
        import BufferMemMgr
        bufMgr = iotask.createSvc("BufferMemMgr")

    print("== ROOT IO Svc ==")
    import RootIOSvc
    ro = iotask.createSvc("RootOutputSvc/OutputSvc")
    output_streams = {}
    if args.anamgr_edm:
        output_streams["/Event/Sim"] = args.output
    ro.property("OutputStreams").set(output_streams)
    # = Data Buffer =
    print("== Buffer Memory Management ==")
    import BufferMemMgr
    bufMgr = task.createSvc("BufferMemMgr")

    # = random svc =
    print("== Random Svc ==")
    import RandomSvc
    task.property("svcs").append("RandomSvc")
    rndm = task.find("RandomSvc")
    rndm.property("Seed").set(args.seed)
    if args.restore_seed_status:
        # == maybe this is a file? ==
        import os.path
        if os.path.exists(args.restore_seed_status):
            filename = args.restore_seed_status
            with open(filename) as f:
                for line in f:
                    print(line)
                    l = line.strip()
                    break
        else:
            l = args.restore_seed_status
        import re
        l = re.split(',\s*|\s+', l)
        seedstatus = [int(i) for i in l if i.isdigit()]
        print("loaded seed status: ", seedstatus)
        rndm.property("SeedStatusInputVector").set(seedstatus)

    
    # = root writer =
    print("== Root Writer ==")
    import RootWriter
    rootwriter = task.createSvc("RootWriter")

    rootwriter.property("Output").set({"SIMEVT":args.user_output})

    # = global time =
    global_time_enabled = args.global_event_rate > 0.
    if global_time_enabled:
        import MCGlobalTimeSvc
        globaltime = task.createSvc("MCGlobalTimeSvc")
        globaltime.property("BeginTime").set(args.global_time_begin)
        globaltime.property("EndTime").set(args.global_time_end)
        globaltime.property("EventRate").set(args.global_event_rate) # Hz

    # = timer svc =
    try:
        import JunoTimer
        task.createSvc("JunoTimerSvc")
    except:
        pass

    # = geometry service for PMT =
    import Geometry
    pmt_param_svc = task.createSvc("PMTParamSvc")
    tt_geom_svc = task.createSvc("TTGeomSvc")

    # = pmt parameter service =
    import PMTSimParamSvc
    print(" == PMTSimParamSvc == ")
    pmt_sim_param_svc = task.createSvc("PMTSimParamSvc")
    pmt_sim_param_svc.property("DBType").set(args.dbtype)

    # = generator related =
    print("GENTOOL MODE: ", args.gentool_mode)
    if args.gentool_mode == "gun":
        setup_generator(task)
    elif args.gentool_mode == "photon":
        # using optical photon
        setup_generator_photon(task)
        # disable several anamgrs
        args.anamgr_genevt = False
        args.anamgr_deposit = False
    elif args.gentool_mode == "gendecay":
        setup_generator_gendecay(task)
    elif args.gentool_mode == "hepevt":
        setup_generator_hepevt(task)
    elif args.gentool_mode == "beam":
        setup_generator_beam(task)
    elif args.gentool_mode == "sn":
        setup_generator_sn(task)
    elif args.gentool_mode == "nusol":
        setup_generator_nusol(task)
    elif args.gentool_mode == "atmo":
        setup_generator_atmo(task)
    elif args.gentool_mode == "neutron":
        setup_generator_neutron(task)
    gt = task.find("GenTools")
    gt.property("EnableGlobalTime").set(global_time_enabled)
    #gt.setLogLevel(2)
    # = setup start index =
    #   Note: even though we could set event index in DetSimAlg, it would be
    #         more consistent when EvtNavigator is created with event id.
    gt.property("EvtID").set(args.start_evtid)
    
    ##########################################################################
    # = MC parameters =
    ##########################################################################
    Sniper.loadDll("libMCParamsSvc.so")
    mcsvc = task.createSvc("MCParamsFileSvc/MCParamsSvc")

    ##########################################################################
    # = OP Simulator =
    ##########################################################################
    Sniper.loadDll("libOPSimulator.so")
    opsim = task.createSvc("OPSimSvc")

    ##########################################################################
    # = opticks =
    ##########################################################################

    if args.opticks_anamgr:
        g4okbr_root = os.environ.get("G4OPTICKSBRIDGEROOT",None) 
        if g4okbr_root is None:
            log.fatal(" --opticks-anamgr can only be used when non-standard G4OpticksBridge package is built and setup")
            assert 0
        else:
            log.info("[loadDll libG4OpticksBridge.so --opticks-anamgr ")   
            Sniper.loadDll("libG4OpticksBridge.so")   
            log.info("]loadDll libG4OpticksBridge.so")   
        pass
    else:
        log.info(" not loading libG4OpticksBridge as --opticks-anamgr not requested" )
    pass 

    ##########################################################################
    # = geant4 related =
    ##########################################################################
    import DetSimOptions
    sim_conf = None
    if args.detoption == "Acrylic":
        from DetSimOptions.ConfAcrylic import ConfAcrylic
        acrylic_conf = ConfAcrylic(task)
        acrylic_conf.configure()
        sim_conf = acrylic_conf
        ## Chimney is enabled by default
        #enable the blow 2 lines only when you want to define new geomery parameters.
        #acrylic_conf.set_top_chimney(3.5, 0.1) #(upper_chimney_height,inner_reflectivity) 
        #acrylic_conf.set_lower_chimney(0.3, 0.1) # (blocker_Z_position, inner_reflectivity)

        #acrylic_conf.disable_chimney() # enable this line to simulate without chimney

    elif args.detoption == "Balloon":
        from DetSimOptions.ConfBalloon import ConfBalloon
        balloon_conf = ConfBalloon(task)
        balloon_conf.configure()
        sim_conf = balloon_conf

    if sim_conf:
        cd = sim_conf.cd()
        cd.property("CheckOverlap").set(False)
        # = commissioning =
        if args.commissioning_enabled:
            cd.property("IsCommissioning").set(True)
            cd.property("ZbelowWater").set(args.below_z_is_water)
        
        # cd.setLogLevel(0) # enable when debug contruction of cd.
        detsimfactory = sim_conf.detsimfactory()

        detsimfactory.property("PhysicsList").set(args.physics_list)
        # = detector components =
        detsimfactory.property("CDEnabled").set(args.cd_enabled)
        detsimfactory.property("WPEnabled").set(args.wp_enabled)
        detsimfactory.property("TTEnabled").set(args.tt_enabled)
        if args.shutter:
            acrylic_conf.enable_shutter()
        # = analysis manager control =
        # == reset the anamgr list to data model writer ==
        detsimfactory.property("AnaMgrList").set([])
        # == op simulator interface ==
        if args.deferred_op:
            detsimfactory.property("AnaMgrList").append("OPSimAnaMgr")
        # == opticks
        if args.opticks_mode > 0 and args.opticks_anamgr:
            if g4okbr_root == None:
                print("WARNING : Opticks needs non-standard G4OpticksBridge package to be installed and setup")    
            else:
                detsimfactory.property("AnaMgrList").append("G4OpticksAnaMgr")
                g4ok_anamgr = sim_conf.tool("G4OpticksAnaMgr")
                g4ok_anamgr.setLogLevel(4)
            pass
        pass
        # == grdm (Geant4 Radioactivity Decay Module) ==
        if args.anamgr_grdm:
            detsimfactory.property("AnaMgrList").append("RadioAnaMgr")
            ram = sim_conf.tool("RadioAnaMgr")
            ram.property("StopAtPa234m").set(args.stopAtPa234m)

        # == edm (event data model) ==
        if args.anamgr_edm:
            detsimfactory.property("AnaMgrList").append("DataModelWriter")
        # == TT ==
        if args.anamgr_tt:
            sim_conf.set_tt_edep_output()
        # == if split mode enable, disable others ==
        if args.anamgr_edm and args.splitoutput:
            detsimfactory = sim_conf.detsimfactory()
            detsimfactory.property("AnaMgrList").append(
                                        "DataModelWriterWithSplit",
                                        )
            dmwws = sim_conf.tool("DataModelWriterWithSplit")
            dmwws.property("HitsMax").set(args.split_maxhits)
        # == normal anamgr ==
        if args.anamgr_normal:
            detsimfactory.property("AnaMgrList").append("NormalAnaMgr")
            if args.anamgr_normal_hit:
               normal_anamgr=sim_conf.tool("NormalAnaMgr")
               normal_anamgr.property("EnableHitInfo").set(args.anamgr_normal_hit)
        # == genevt anamgr ==
        if args.anamgr_genevt:
            detsimfactory.property("AnaMgrList").append("GenEvtInfoAnaMgr")
        # == deposit anamgr ==
        if args.anamgr_deposit:
            detsimfactory.property("AnaMgrList").append("DepositEnergyAnaMgr")
        if args.anamgr_deposit_tt:
            detsimfactory.property("AnaMgrList").append("DepositEnergyTTAnaMgr")
            tt_anamgr = sim_conf.tool("DepositEnergyTTAnaMgr")
            tt_anamgr.property("EnableNtuple").set(True)
        # == interesting process ==
        if args.anamgr_interesting_process:
            detsimfactory.property("AnaMgrList").append("InteresingProcessAnaMgr")
        # == optical parameter ==
        if args.anamgr_optical_parameter:
            detsimfactory.property("AnaMgrList").append("OpticalParameterAnaMgr")
        # == timer anamgr ==
        if args.anamgr_timer:
            detsimfactory.property("AnaMgrList").append("TimerAnaMgr")
            timer = acrylic_conf.tool("TimerAnaMgr")
            timer.setLogLevel(3)
        # == photon tracking ==
        if args.anamgr_photon_tracking:
            detsimfactory.property("AnaMgrList").append("PhotonTrackingAnaMgr")
        # == append other anamgr into the list ==
        for tmp_anamgr in args.anamgr_list:
            detsimfactory.property("AnaMgrList").append(tmp_anamgr)
        # == for extension: load config file =
        import os.path
        if args.anamgr_config_file and os.path.exists(args.anamgr_config_file):
            # 
            print("Loading config file: '%s'"%args.anamgr_config_file)
            execfile(args.anamgr_config_file)
            # exec(compile(open(args.anamgr_config_file, "rb").read(), args.anamgr_config_file, 'exec'))


        # global geom info
        geom_info = acrylic_conf.tool("GlobalGeomInfo")
        geom_info.property("LS.AbsLen").set(args.mat_LS_abslen)
        geom_info.property("LS.RayleighLen").set(args.mat_LS_raylen)
        geom_info.property("LS.LightYield").set(args.light_yield)
        geom_info.property("UseParamSvc").set(args.use_param_svc)

        # voxel fast simulation. need to disable the optical progress
        if args.voxel_fast_sim:
            print("voxel method enabled")
            # disable pmts and struts
            if not args.voxel_pmts_structs:
                print("disable pmts and structs")
                acrylic_conf.disable_pmts_and_struts_in_cd()
            import os
            if os.environ.get("DETSIMOPTIONSROOT", None) is None:
                print("Missing DETSIMOPTIONSROOT")
                sys.exit(-1)
            dp = lambda f: os.path.join(os.environ.get("DETSIMOPTIONSROOT"),
                                        "share", "examples", "voxelmuon", f)
            if (os.environ.get("VOXELFASTDIR")):
                dp = lambda f: os.path.join(os.environ.get("VOXELFASTDIR"), f)
            if args.voxel_fast_dir:
                dp = lambda f: os.path.join(args.voxel_fast_dir, f)
            acrylic_conf.add_anamgr("MuonFastSimVoxel")
            mfsv = acrylic_conf.tool("MuonFastSimVoxel")
            mfsv.property("GeomFile").set(dp("geom-geom-20pmt.root"))
            mfsv.property("NPEFile").set(dp("npehist3d_single.root"))
            mfsv.property("HitTimeMean").set(dp("hist3d.root"))
            mfsv.property("HitTimeRes").set(dp("dist_tres_single.root"))
            mfsv.property("MergeFlag").set(args.voxel_merge_flag)
            mfsv.property("MergeTimeWindow").set(args.voxel_merge_twin)
            mfsv.property("EnableNtuple").set(args.voxel_fill_ntuple)
            mfsv.property("QuenchingFactor").set(args.voxel_quenching_scale)
            mfsv.property("SampleNPE").set(args.voxel_gen_npe)
            mfsv.property("SampleTime").set(args.voxel_gen_time)
            mfsv.property("SaveHits").set(args.voxel_save_hits)

        # split the primary track. step -> sub event
        if args.splittrack:
            acrylic_conf.add_anamgr("PostponeTrackAnaMgr")
            pta = acrylic_conf.tool("PostponeTrackAnaMgr")
            pta.property("SplitMode").set(args.track_split_mode)
            pta.property("TimeCut").set(args.track_split_time)
            #pta.setLogLevel(2)

        # physics list
        # = em =
        em_process = acrylic_conf.em_process()
        em_process.property("UsePositronium").set(args.positronium)

        # disable the optical progress
        op_process = acrylic_conf.optical_process()
        op_process.property("OpticksMode").set(args.opticks_mode)  # see DsPhysConsOptical
        op_process.property("UseCerenkov").set(args.cerenkov)
        if not args.useoptical or args.voxel_fast_sim:
            print("Disable Optical Process")
            op_process.property("UseScintillation").set(False)
            op_process.property("UseCerenkov").set(False)
        if args.cerenkov_only:
            print("Enable Cerenkov. (note: Scintillation is used to do reemission only)")
            op_process.property("UseScintillation").set(True)
            op_process.property("UseCerenkov").set(True)
            op_process.property("ScintDoReemissionOnly").set(True)

        # For testing of OPSimulator only
        if args.deferred_op:
            op_process.property("UseScintillation").set(True)
            op_process.property("UseCerenkov").set(False)
            op_process.property("ScintDoReemissionOnly").set(True)

        op_process.property("UseQuenching").set(args.quenching)
        # new optical model
        if args.new_optical_model:
            op_process.property("UseAbsReemit").set(True)
            op_process.property("UseScintSimple").set(True)
        else:
            op_process.property("UseAbsReemit").set(False)
            op_process.property("UseScintSimple").set(False)

       # op_process.property("UseAbsReemit").set(args.absreemit)
       # op_process.property("UseScintSimple").set(args.scintsimple)
        # other flags:
        op_process.property("doTrackSecondariesFirst").set(args.track_op_first)
            
        # == beam mode ==
        setup_calib_pelletron(acrylic_conf, args.pelletron)
        setup_calib_unit(acrylic_conf, args.guide_tube)
        # == geant4 run mac ==
        detsimalg = sim_conf.detsimalg()
        detsimalg.property("RunCmds").set([
                     #"/run/initialize",
                     #"/tracking/verbose 2",
                     #"/process/inactivate Scintillation",
                     #"/process/inactivate Cerenkov",
                 ])
        if args.mac:
            if os.path.exists(args.mac):
                detsimalg.property("RunMac").set(args.mac)
            else:
                print("WARNING: mac file %s does not exist."%args.mac)
                detsimalg.property("RunMac").set("")
        if args.vis:
            if os.path.exists(args.vis_mac):
                detsimalg.property("VisMac").set(args.vis_mac)
            else:
                print("WARNING: vis mac file %s does not exist."%args.vis_mac)
                detsimalg.property("VisMac").set("")
        

        # == QE scale ==
        sim_conf.set_qe_scale(args.qescale)
        # == 20inch PMT ==
        log.info("PMTName %s --pmt20inch-name " % args.pmt20inch_name)
        log.info("LPMTExtra %s --pmt20inch-extra " % args.pmt20inch_extra)

        detsimfactory.property("OpticksMode").set(args.opticks_mode)

        # NB: these are internal envvars, externally set envvars with these keys are ignored
        # use command line arguments --pmt20inch-polycone-neck or --pmt20inch-simplify-csg to set them 
        for key in "JUNO_PMT20INCH_POLYCONE_NECK JUNO_PMT20INCH_SIMPLIFY_CSG".split(" "):
            attn = key.replace("JUNO_","").lower()
            att = getattr( args, attn, None )
            # assert att, key 
            if att == True:
                log.info("setting key %s from args.%s  " % (key, attn))
                os.environ[key] = "ENABLED"  
            else:
                if key in os.environ:
                    log.info("un-setting key %s from args.%s  " % (key, attn))
                    os.environ.pop(key)
                pass
            pass
        pass    

        # add LS AbsLength Mode #
        if args.new_optical_model:
            detsimfactory.property("GdLSAbsLengthMode").set(1) # 1 for LAB AbsLength
        else:
            detsimfactory.property("GdLSAbsLengthMode").set(0) # 0 for old LS AbsLength
        #-----------------------------------------#     
       
        detsimfactory.property("PMTName").set(args.pmt20inch_name)
        detsimfactory.property("LPMTExtra").set(args.pmt20inch_extra)
        if args.pmt20inch_name == "R12860":
            r12860 = sim_conf.tool("R12860PMTManager/PMT_20inch")
            r12860.property("FastCover").set(True)
            r12860.property("FastCoverMaterial").set("Water")
        elif args.pmt20inch_name == "PMTMask":
            mask = sim_conf.tool("R12860MaskManager")
            mask.property("TopThickness").set(args.pmtmask_top_thick)
        pass

        ## debug: add cover for NNVT MCP-PMT

        if args.pmt20inch_extra == "TWO":
            log.info("TWO . args.pmt20inch_extra %s " % args.pmt20inch_extra)  
            nnvt_mcp_pmt = sim_conf.tool("NNVTMCPPMTManager/NNVTMCPPMT")
            nnvt_mcp_pmt.property("FastCover").set(True)
            nnvt_mcp_pmt.property("FastCoverMaterial").set("Water")

            hamamatsu_pmt = sim_conf.tool("HamamatsuR12860PMTManager/HamamatsuR12860")
            hamamatsu_pmt.property("FastCover").set(True)
            hamamatsu_pmt.property("FastCoverMaterial").set("Water")
        else:
            log.info("args.pmt20inch_extra %s " % args.pmt20inch_extra)  
        pass

        # == enable or disable 20inch PMTs ==
        if not args.pmt20inch:
            detsimfactory.property("PMTPosFile").set("")
        # == enable or disable 3inch PMTs ==
        if not args.pmt3inch:
            sim_conf.disable_3inch_PMT()
        else:
            print("3inch PMT type: ", args.pmt3inch_name)
            sim_conf.set_3inch_pmt_name(args.pmt3inch_name)
            sim_conf.set_3inch_pmt_offset(args.pmt3inch_offset)

        # == enable or disable PMTs in WP ==
        if not args.wp_pmt_enabled:
            detsimfactory.property("VetoPMTPosMode").set("")
        if args.pmtsd_v2:
            sim_conf.enable_PMTSD_v2()
            pmtsdmgr = sim_conf.pmtsd_mgr()
            pmtsdmgr.property("CollEffiMode").set(args.ce_mode)
            pmtsdmgr.property("CEFlatValue").set(args.ce_flat_value)
            pmtsdmgr.property("OpticksMode").set(args.opticks_mode)
        if args.pmtsd_merge_twindow>0:
            pmtsdmgr = sim_conf.pmtsd_mgr()
            pmtsdmgr.property("EnableMergeHit").set(True)
            pmtsdmgr.property("MergeTimeWindow").set(args.pmtsd_merge_twindow)
        # pmt hit type
        pmtsdmgr = sim_conf.pmtsd_mgr()
        pmtsdmgr.property("HitType").set(args.pmt_hit_type)
        pmtsdmgr.property("DisableSD").set(args.pmt_disable_process)
        pmtsdmgr.property("PMTSD").set("junoSD_PMT_muon")
        if args.muon_split:
            pmtsdmgr.property("split").set(True)
            pmtsdmgr.property("maxhit").set(args.muon_maxhit)
        else:
            pmtsdmgr.property("split").set(False)
        if args.ce_func:
            pmtsdmgr.property("CEFunction").set(args.ce_func)
        if args.ce_func_par:
            pmtsdmgr.property("CEFuncParams").set(args.ce_func_par)
        # disable struts and fasteners
        if args.flag_struts_fasteners != "none":
            sim_conf.disable_struts_in_cd(args.flag_struts_fasteners)
        # gdml output
        if gdml_filename:
            sim_conf.set_gdml_output(gdml_filename)
        if dae_filename:
            sim_conf.set_dae_output(dae_filename)
        pass
    pass

    ##########################################################################
    # = begin run =
    ##########################################################################
    log.info("task.show")
    task.show()
    log.info("task.run")
    task.run()
    log.info("DONE")
