// $Id: EventAnlProc.h 755 2011-05-20 08:04:11Z linev $
// Adapted for DESPEC by A.K.Mistry 2020
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum f�r Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef EVENTANLPROCESSOR_H
#define EVENTANLPROCESSOR_H

#include "TGo4EventProcessor.h"

#include "EventUnpackStore.h"
#include "EventAnlStore.h"
#include "CalibParameter.h"
#include "CorrelParameter.h"
#include "AIDA_Headers.h"
#include "AIDA_Event.h"
#include "AIDA_Data_Types.h"
#include "TCutG.h"
#include "TGraph.h"
#include "Go4ConditionsBase/TGo4WinCond.h"
#include "Go4ConditionsBase/TGo4PolyCond.h"
#include "Configuration_Files/DESPEC_General_Setup/DESPEC_Setup_File.h"
#include "TFRSParameter.h"

//These are for TAMEX
#define CYCLE_TIME  (Double_t) 5000
#define COARSE_CT_RANGE  0x800  // 11 bits

//Uncomment this to align the AIDA ASICs with a pulser
// Only needed if the ASICs didn't align properly
//#define AIDA_PULSER_ALIGN

class EventAnlStore;


class EventAnlProc : public TGo4EventProcessor {
public:
  EventAnlProc();
  EventAnlProc(const char * name);
  virtual ~EventAnlProc();
  CalibParameter *fCal;
  CorrelParameter *fCorrel;
  //AidaAnlData pAida;

  virtual Bool_t BuildEvent(TGo4EventElement* dest);
  virtual void UserPostLoop();
  TIDParameter* frs_id;


  void get_used_systems();
	long long lastGeWR = 0;
  long long lastbPlastWR = 0;
  long long lastFRSWR = 0;

  //double lead_lead_bplas[48][100], trail_trail_bplas[48][100];

  double ToT_bplas[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;
  double lead_bplas[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;
  double trail_bplas[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;
  //double lead_lead_bplas_gated[48][100];
  double lead_lead_bplas_Ref1[bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS];
  double lead_lead_bplas_Ref2[bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS];
  double lead_lead_bplas_Ref3[bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS];

  double SC41L_ANA_lead_bPlas[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS], SC41R_ANA_lead_bPlas[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS];
  double SC41L_DIG_lead_bPlas[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS], SC41R_DIG_lead_bPlas[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS];
  double Ge_Trig_lead_bPlas[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS];
  int hits_bplas_lead = 0, hits_bplas_trail=0;
  int hits_bplas_lead_fast = 0, hits_bplas_trail_fast=0;
  int hits_bplas_lead_slow = 0, hits_bplas_trail_slow=0;
        
  double ToT_bplas_Fast[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;
  double ToT_bplas_Slow[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;
  double lead_bplas_fast[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;
  double lead_bplas_slow[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;
  double trail_bplas_fast[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;
  double trail_bplas_slow[4][bPLASTIC_CHAN_PER_DET][bPLASTIC_TAMEX_HITS] ;

  double bPlas_TAM_SC41L_ANA;
  double bPlas_TAM_SC41R_ANA;
  double bPlas_TAM_SC41L_DIG;
  double bPlas_TAM_SC41R_DIG;
  double bPlas_TAM_Ge_TRIG;

  double bPlas_RefCh0_Det1[bPLASTIC_TAMEX_HITS];
  double bPlas_RefCh0_Det2[bPLASTIC_TAMEX_HITS];
  double bPlas_RefCh0_Det3[bPLASTIC_TAMEX_HITS];

  //FRS Histograms and conditions setup
  TH1I* MakeH1I(const char* foldername,
    const char* histoname,
    Int_t nbins,
    Float_t xmin, Float_t xmax,
    const char* xtitle = "channels",
    Color_t linecolor = 2,
    Color_t fillcolor = 6,
    const char* ytitle = 0);

  TH2I* MakeH2I(const char* foldername,
    const char* histoname,
    Int_t nbinsx, Float_t xmin, Float_t xmax,
    Int_t nbinsy, Float_t ymin, Float_t ymax,
    const char* xtitle, const char* ytitle,
    Color_t marker);

  TGo4WinCond* MakeWindowCond(const char* foldername,
    const char* condname,
    float left = 0.,
    float right = 4096.,
    const char* HistoName = 0);

  TGo4PolyCond* MakePolyCond(const char* foldername,
    const char* condname,
    Int_t size,
    Float_t (*points)[2],
    const char* HistoName = 0);

  Bool_t Check_PolyCond_Multi_X_Y(Float_t X, Float_t Y, Float_t*** V, int n, int cond_num);
  Bool_t Check_PolyCond_X_Y(Float_t X, Float_t Y, Float_t** V, int n );
  void FRS_Gates();

  void ProcessAida(EventUnpackStore* pInput, EventAnlStore* pOutput);

  int PrcID_Conv[7];
  int Used_Systems[7];
  int event_number;
  bool VMEorTAMEX_bPlas;
  void get_used_Systems();
  void load_GermaniumMap_File();
  //void checkTAMEXorVME();

  Float_t  FRS_dE[2],  FRS_dE_cor[2];
  Float_t  FRS_sci_l[12],  FRS_sci_r[12],  FRS_sci_e[12],  FRS_sci_tx[12],  FRS_sci_x[12];
  Float_t FRS_tpc_x[7],FRS_tpc_y[7];
  Float_t  FRS_sci_tofll2,  FRS_sci_tofll3, FRS_sci_tof2, FRS_sci_tofrr2, FRS_sci_tofrr3, FRS_sci_tof3;
  Float_t  FRS_ID_x2, FRS_ID_y2, FRS_ID_a2, FRS_ID_b2;
  Float_t  FRS_ID_x4, FRS_ID_y4, FRS_ID_a4, FRS_ID_b4;
  Int_t    FRS_sci_dt_21l_21r, FRS_sci_dt_41l_41r, FRS_sci_dt_42l_42r, FRS_sci_dt_43l_43r;
  Int_t    FRS_sci_dt_21l_41l, FRS_sci_dt_21r_41r, FRS_sci_dt_21l_42l, FRS_sci_dt_21r_42r;
  Float_t  FRS_ID_brho[2], FRS_ID_rho[2];
  Float_t  FRS_beta, FRS_beta3, FRS_gamma;
  Float_t  FRS_AoQ, FRS_AoQ_corr;
  Float_t  FRS_z, FRS_z2, FRS_z3, z1_corr, z2_corr;
  Float_t  FRS_dEdeg, FRS_dEdegoQ;
  //Float_t  FRS_AoQ_mhtdc, FRS_AoQ_corr_mhtdc;
  //Float_t  FRS_z_mhtdc, FRS_z2_mhtdc;
  Float_t  FRS_dEdeg_mhtdc[10], FRS_dEdegoQ_mhtdc[10];
  Float_t  FRS_beta_mhtdc[10];
  Float_t  FRS_timestamp, FRS_ts, FRS_ts2;
  bool FRS_spill;
  Long64_t FRS_time_mins;
  Float_t  FRS_AoQ_mhtdc[10], FRS_AoQ_corr_mhtdc[10] , FRS_z_mhtdc[10], FRS_z2_mhtdc[10];
  int num_ID_Z_AoQ = 0;
  int num_ID_x2AoQ = 0;
  int num_ID_x4AoQ = 0;
  int num_ID_Z_Z2 = 0;
  int num_ID_x4AoQ_mhtdc = 0;
  int num_ID_x2AoQ_mhtdc = 0;
  int num_ID_Z_Z2_mhtdc = 0;
  int num_ID_Z_AoQ_mhtdc = 0;
  int num_ID_dEvsBRho = 0;
  int num_ID_dEvsZ = 0;
  int num_ID_dEvsZ_x4AoQ_p1 = 0;
  int num_ID_dEvsZ_x4AoQ_m1 = 0;

  //Float_t  FRS_z_shifted, FRS_z2_shifted;
  int Z_Shift_array;
  Float_t FRS_WR_a[200];
  Float_t FRS_WR_b[200] ;
  Float_t Z1_shift_value[200] ;
  Float_t Z2_shift_value[200] ;
  //FRS_AoQ_shifted;
  int AoQ_Shift_array;
  Float_t FRS_WR_i[200];
  Float_t FRS_WR_j[200] ;
  Float_t AoQ_shift_value[200] ;
  Float_t AoQ_shift_TPC_value[200] ;
  Float_t AoQ_shift_Sci21_value[200] ;
  Float_t AoQ_shift_Sci22_value[200] ;

  Int_t ZAoQgnum,Z1Z2gnum,X2AoQgnum,X4AoQgnum,dEdeggnum;
  Float_t X_ZAoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints],Y_ZAoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t X_ZZ2[MAX_FRS_GATE][MAX_FRS_PolyPoints],Y_ZZ2[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t X_ZZ2_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints],Y_ZZ2_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t XX4_AoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints], YX4_AoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t XX4_AoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints], YX4_AoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t XX2_AoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints], YX2_AoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t XX2_AoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints], YX2_AoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t X_dEvsBRho[MAX_FRS_GATE][MAX_FRS_PolyPoints], Y_dEvsBRho[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t X_dEvsZ[MAX_FRS_GATE][MAX_FRS_PolyPoints], Y_dEvsZ[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t X_dEvsZ_x4AoQ_p1[MAX_FRS_GATE][MAX_FRS_PolyPoints], Y_dEvsZ_x4AoQ_p1[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t X_dEvsZ_x4AoQ_m1[MAX_FRS_GATE][MAX_FRS_PolyPoints], Y_dEvsZ_x4AoQ_m1[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  Float_t X_ZAoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints],Y_ZAoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints];
  


  int Aida_Fired;

  //int    bPlasQDCFired;
  //int    bPlasQDCID[32];
  //double bPlasQDC[32];
  //double bPlasQDCGain[32];

  //int    bPlasTDCFired;
  //int    bPlasTDCID[50];
  //double bPlasTDC_TS[50][32];
  //double bPlasTDC_ref;
  //int    bPlas_TDC_Multiplicity[32];
  //double bPlas_TDC_diff[50];
  //double bPlas_TDC_diff_sum;
  
  int    ScalarFired;
  int    ScalarID;
  //int    Fatmult;
  int    SC40mult;
  int    SC41mult;
  int    stdcmult;
  int    sqdcmult;


  // From unpacker
  int    GeFired;
  int    GeDet[Germanium_MAX_CHANNELS];
  int    GeCrys[Germanium_MAX_CHANNELS];
  double GeE[Germanium_MAX_CHANNELS];
  double GeE_Cal[Germanium_MAX_CHANNELS];
  ULong64_t GeT[Germanium_MAX_CHANNELS];
  ULong64_t GeCF_T[Germanium_MAX_CHANNELS];
  ULong64_t RefTGe, RefCFDGe;
  Long64_t Ge_time_mins;
  double GeEventT[Germanium_MAX_CHANNELS];
  bool GePileUp[Germanium_MAX_CHANNELS];
  bool GeOverFlow[Germanium_MAX_CHANNELS];


  int Gam_mult;
  Long64_t dT_Ge,dT_Ge_cfd;
  Long64_t dT_Align,dT_CFD_Align;
  Long64_t Ge_Talign[Germanium_MAX_HITS],Ge_cfd_Talign[Germanium_MAX_HITS];
  Long64_t dT_Addback;

  Long64_t Ge_WR;

  int    Fing_firedTamex;
  int    Fing_leadHits[4];
  int    Fing_trailHits[4];
  int    Fing_iterator[4];
  double Fing_trig[4];
  int    Fing_tamex_ch[4][32];
  int    Fing_leadChan[4][32];
  int    Fing_trailChan[4][32];
  double Fing_lead_coarse[4][32];
  double Fing_lead_fine[4][32];
  double Fing_trail_coarse[4][32];
  double Fing_trail_fine[4][32];

  double Fing_leadT[4][32];
  double Fing_trailT[4][32];
  double Fing_TOT[4][32];
  double Fing_TOT_added[4][32];
  int    Fing_chID[4][32];
  double dataSetPerEvent[50];
  double pmtSetPerEvent[50];
  double totaltimeEvent[50];
  double maxToT;
  int    maxToTChan;
  double maxToT_added;
  int    maxToT_added_Chan;

  void Make_FRS_Histos();
  void Make_Aida_Histos();
  //void Make_Plastic_VME_Histos();
  void Make_Plastic_Twinpeaks_Histos();
  void Make_Plastic_Tamex_Histos();
  void Make_Germanium_Histos();
  void Make_Finger_Histos();
  void Make_WR_Histos();

  void Make_Fing_Plas_Histos();
  void FRS_GainMatching();

  void Process_FRS_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput);
  //void Process_Plastic_VME_Histos(EventAnlStore* pOutput);
  void Process_Plastic_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput);
  void Process_Plastic_Twinpeaks_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput);
  void Process_Germanium_Histos(EventAnlStore* pOutput);
  void Process_Finger_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput);
  void Process_Fing_Plas_Histos(EventAnlStore* pOutput);
  void Process_WR_Histos(EventUnpackStore* pInput);
  // TH1 *GermaniumCal;



  long lastTime;
  int ID;
  AidaEvent old;
  AidaEvent evt;
  long startTime;
  long stopTime;

  //GeEvent gal;

  /* Multiplexer correction */

  int totalEvents;
  int implantEvents;
  int goodImplantEvents;
  int stoppedEvents;
  int decayEvents;
  int pulserEvents;
  int nonsenseEvents;

  std::vector<TH2*> implants_strip_xy;
	std::vector<TH2*> implants_strip_xy_stopped;
  std::vector<TH2*> implants_pos_xy;
	std::vector<TH2*> implants_pos_xy_stopped;
  std::vector<TH1*> implants_e;
  std::vector<TH2*> implants_e_xy;
  std::vector<TH1*> implants_time_delta;
  std::vector<TH1*> implants_strip_1d;
  std::vector<TH1*> implants_per_event;
  std::vector<TH1*> implants_channels;
  std::vector<TH2*> implants_x_ex;
  std::vector<TH2*> implants_y_ey;
  #ifdef AIDA_PULSER_ALIGN
  TH2* aida_pulser_time;
  #endif

  std::vector<TH2*> decays_strip_xy;
  std::vector<TH2*> decays_pos_xy;
  std::vector<TH1*> decays_e;
  std::vector<TH2*> decays_e_xy;
  std::vector<TH1*> decays_time_delta;
  std::vector<TH1*> decays_strip_1d;
  std::vector<TH1*> decays_per_event;
  std::vector<TH1*> decays_channels;

  std::vector<AidaCluster> EventsToClusters(std::vector<AidaEvent> const&);
  AidaHit ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const&);

  int      IsData(std::ifstream &f);


  TH1 *hGe_deadtime;
  TH1 *hbPlast_deadtime;
  TH1 *hFRS_deadtime;
  TH1 *hAida_Ge_WRdT;
  TH1 *hAida_bPlas_WRdT;
  TH1 *hbPlas_Ge_WRdT;
  TH1 *hFRS_Ge_WRdT;
  TH1 *hFRS_bPlas_WRdT;

  //FRS Histograms
  TH1  *hID_Z1_corr;
  TH1  *hID_Z2_corr;
  TH2  *hID_Z1_vs_T;
  TH2  *hID_Z2_vs_T;
  TH2  *hID_AoQ_vs_T;
  TH2  *hID_Z_AoQ;
  TH2  *hID_Z_AoQ_corr;
  TH2  *hID_Z_AoQ_zsame; 
  TH2  *hID_x4AoQ_zsame;
  TH2  *hID_x2AoQ_zsame;
  TH2  *hID_x2AoQ;
  TH2  *hID_x4AoQ;
  TH2  *hID_x2AoQ_corr;
  TH2  *hID_x4AoQ_corr;
  TH2  *hID_a2AoQ;
  TH2  *hID_a4AoQ;
  TH2  *hID_dEdegoQ_Z;
  TH2  *hID_dEdeg_Z;
  TH2  *hID_Z_Z2;
  TH2  *hID_Z_dE2;
  TH2  *hSCI_dE24;
  TH2  *hID_Z_Sc21E; 
  TH2  *hID_Z_Sc21E_mhtdc;
  TH2  *hID_x2z1;
  TH2  *hID_x4z1;
  TH2  *hID_x4z2;
  TH2  *hID_E_x4;
  TH2  *hID_E_x2;
  TH2  *hID_x2a2;
  TH2  *hID_y2b2;
  TH2  *hID_x4a4;
  TH2  *hID_y4b4;
  TH2  *hID_x2x4;
  TH2  *hID_SC41dE_AoQ;
  TH2  *hID_dEToF;
  TH2  *hID_dEBRho; //ADDED BY GEE 07/12/22
  TH2  *hID_dBrho_dE4;
  TH2  *hID_BRho1v2;
  TH2  *hID_gamgam;
	    
	TH2 *hID_Z_E_SC41;
	TH2 *hID_Z_E_SC42;
                 
  TH2  *hID_AoQ_mhtdc_T;
  TH2  *hID_Z_mhtdc_T;
  TH2  *hID_Z_AoQ_mhtdc;
  TH2  *hID_Z_AoQ_corr_mhtdc;
  TH2  *hID_Z_Z2_mhtdc;
  TH2  *hID_Z_AoQ_zsame_mhtdc;
  TH2  *hID_x4AoQ_zsame_mhtdc;
  TH2  *hID_x2AoQ_mhtdc;
  TH2  *hID_x2AoQ_zsame_mhtdc;
  TH2  *hID_x4AoQ_mhtdc;
  TH2  *hID_a2AoQ_mhtdc;
  TH2  *hID_a4AoQ_mhtdc;
  TH2  *hID_dEdegoQ_Z_mhtdc;
  TH2  *hID_dEdeg_Z_mhtdc;
  TH2  *hID_Z_dE2_mhtdc;
  TH2  *hID_x2z_mhtdc;
  TH2  *hID_x4z_mhtdc;
  // TH2 *hID_Z_AoQ_mhtdc_first_hit;
  // TH2 *hID_Z_AoQ_mhtdc_excluding_first_hit;
  // TH2 *hID_Z_AoQ_corr_mhtdc_first_hit;
  // TH2 *hID_Z_AoQ_corr_mhtdc_excluding_first_hit;
            


  TGo4PolyCond  *cID_Z_AoQ[MAX_FRS_GATE];
  TGo4PolyCond  *cID_Z_Z2gate[MAX_FRS_GATE];
  TGo4PolyCond  *cID_x2AoQ[MAX_FRS_GATE];
  TGo4PolyCond  *cID_x4AoQ[MAX_FRS_GATE];
  TGo4PolyCond  *cID_dEvsBRho[MAX_FRS_GATE];
  TGo4PolyCond  *cID_dEvsZ[MAX_FRS_GATE];
  TGo4PolyCond  *cID_dEvsZ_x4AoQ_p1[MAX_FRS_GATE];
  TGo4PolyCond  *cID_dEvsZ_x4AoQ_m1[MAX_FRS_GATE];
  
  TGo4PolyCond  *cID_Z_AoQ_mhtdc[MAX_FRS_GATE];
  TGo4PolyCond  *cID_Z_Z2gate_mhtdc[MAX_FRS_GATE];
  TGo4PolyCond  *cID_x2AoQ_mhtdc[MAX_FRS_GATE];
  TGo4PolyCond  *cID_x4AoQ_mhtdc[MAX_FRS_GATE];
                 
  ///Z AoQ gated
  TH2 *hID_ZAoQ_ZAoQgate[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_ZAoQgate[MAX_FRS_GATE];
  TH2 *hID_x2AoQ_Z1AoQgate[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_Z1AoQgate[MAX_FRS_GATE];
  TH1 *hID_a2_Z1AoQgate[MAX_FRS_GATE];
  TH1 *hID_a4_Z1AoQgate[MAX_FRS_GATE];
                   
  ///Z1 Z2 gated
  TH2 *hID_Z1_Z2gate[MAX_FRS_GATE];
  TH2 *hID_x2AoQ_Z1Z2gate[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_Z1Z2gate[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_Z1Z2gate[MAX_FRS_GATE];
  TH1 *hID_a2_Z1Z2gate[MAX_FRS_GATE];
  TH1 *hID_a4_Z1Z2gate[MAX_FRS_GATE];

                   
  ///X2 AoQ gated
  TH2 *hID_x2AoQ_x2AoQgate[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_x2AoQgate[MAX_FRS_GATE];
  ///Z1 Z2 + X2 AoQ gated
  TH2 *hID_x2AoQ_Z1Z2x2AoQgate[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_Z1Z2x2AoQgate [MAX_FRS_GATE];
  TH2 *hID_ZAoQ_Z1Z2x2AoQgate[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_Z1Z2x2AoQgate[MAX_FRS_GATE];
  TH1 *hID_a2_Z1Z2x2AoQgate[MAX_FRS_GATE];
  TH1 *hID_a4_Z1Z2x2AoQgate[MAX_FRS_GATE];
             
  ///X4 AoQ gated
  TH2 *hID_x4AoQ_x4AoQgate[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_x4AoQgate[MAX_FRS_GATE];
  TH2 *hID_Z1AoQ_x4AoQgate[MAX_FRS_GATE];

  ///Z1 Z2 + X2 AoQ gated
  TH2 *hID_x2AoQ_Z1Z2x4AoQgate[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_Z1Z2x4AoQgate[MAX_FRS_GATE]; 
  TH2 *hID_ZAoQ_Z1Z2x4AoQgate[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_Z1Z2x4AoQgate[MAX_FRS_GATE];
  TH1 *hID_a2_Z1Z2x4AoQgate[MAX_FRS_GATE];
  TH1 *hID_a4_Z1Z2x4AoQgate[MAX_FRS_GATE];
             

  //dEvsBRho gated
  TH2 *hID_ZAoQ_dEvsBRhogate[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_dEvsBRhogate[MAX_FRS_GATE];
  TH2 *hID_dEvsBRho_dEvsBRhogate[MAX_FRS_GATE];

  // dEvsZ gated
  TH2 *hID_ZAoQ_dEvsZgate[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_dEvsZgate[MAX_FRS_GATE];
  TH2 *hID_dEvsZ_dEvsZgate[MAX_FRS_GATE];

  //dEvsZ_x4AoQ deltaQ=0
  TH2 *hID_x4AoQ_dEvsZ_x4AoQgate[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_dEvsZ_x4AoQgate[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_dEvsZ_x4AoQgate[MAX_FRS_GATE];

  //dEvsZ_x4AoQ Z1Z2 deltaQ=0
  TH2 *hID_x4AoQ_dEvsZ_Z1Z2_x4AoQgate[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_dEvsZ_Z1Z2_x4AoQgate[MAX_FRS_GATE];

  //dEvsZ_x4AoQ deltaQ=+1
  TH2 *hID_x4AoQ_dEvsZ_x4AoQ_p1gate[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_dEvsZ_x4AoQ_p1gate[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_dEvsZ_x4AoQ_p1gate[MAX_FRS_GATE];

  //dEvsZ_x4AoQ Z1Z2 deltaQ=+1
  TH2 *hID_x4AoQ_dEvsZ_Z1Z2_p1gate[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_dEvsZ_Z1Z2_p1gate[MAX_FRS_GATE];

  //dEvsZ_x4AoQ deltaQ=-1
  TH2 *hID_x4AoQ_dEvsZ_x4AoQ_m1gate[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_dEvsZ_x4AoQ_m1gate[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_dEvsZ_x4AoQ_m1gate[MAX_FRS_GATE];

  //dEvsZ_x4AoQ Z1Z2 deltaQ=-1
  TH2 *hID_x4AoQ_dEvsZ_Z1Z2_m1gate[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_dEvsZ_Z1Z2_m1gate[MAX_FRS_GATE];


  ///MHTDC
  ///Z vs AoQ Gated
  TH2 *hID_ZAoQ_ZAoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_ZAoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_x2AoQ_Z1AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_Z1AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_dEdegoQ_Z1_Z1AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_dEdegZ1_Z1AoQgate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a2_Z1AoQgate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a4_Z1AoQgate_mhtdc[MAX_FRS_GATE];
             
  ///Z1 vs Z2 Gated
  TH2 *hID_Z1_Z2gate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_x2AoQ_Z1Z2gate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_Z1Z2gate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_Z1Z2gate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a2_Z1Z2gate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a4_Z1Z2gate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_dEdegZ1_Z1Z2gate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_dEdegoQ_Z1_Z1Z2gate_mhtdc[MAX_FRS_GATE];
             
  ///X2 vs AoQ Gated
  TH2 *hID_x2AoQ_x2AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_x2AoQgate_mhtdc[MAX_FRS_GATE];
  ///Z1 Z2 + X2 vs AoQ Gated
  TH2 *hID_x2AoQ_Z1Z2x2AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_Z1Z2x2AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_Z1Z2x2AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_x2AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_dEdegZ1_Z1Z2x2AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_dEdegoQ_Z1_Z1Z2x2AoQgate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a2_Z1Z2x2AoQgate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a4_Z1Z2x2AoQgate_mhtdc[MAX_FRS_GATE];
              
  ///X4 vs AoQ Gated
  TH2 *hID_x4AoQ_x4AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_x4AoQgate_mhtdc[MAX_FRS_GATE];
  ///Z1 Z2 + X4 vs AoQ Gated
  TH2 *hID_x2AoQ_Z1Z2x4AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_Z1Z2x4AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_Z1Z2x4AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_ZAoQ_x4AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_dEdegZ1_Z1Z2x4AoQgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_dEdegoQ_Z1_Z1Z2x4AoQgate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a2_Z1Z2x4AoQgate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a4_Z1Z2x4AoQgate_mhtdc[MAX_FRS_GATE];
              
  ///dE S2 deg vs Z1 Gated
  TH2 *hID_dEdegZ1_dEdegZ1Gated_mhtdc[MAX_FRS_GATE];
  TH2 *hID_Z1AoQ_dEdegZgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_Z1AoQ_zsame_dEdegZgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_Z1Z2_dEdegZgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_x2AoQ_dEdegZgate_mhtdc[MAX_FRS_GATE];
  TH2 *hID_x4AoQ_dEdegZgate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a2_dEdegZgate_mhtdc[MAX_FRS_GATE];
  TH1 *hID_a4_dEdegZgate_mhtdc[MAX_FRS_GATE];
             
  ///bPlast Histograms
  TH1 *hbPlas_ToT_det[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_Multiplicity_Chan[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_Lead_T[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_Trail_T[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_Lead_dT_coinc[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_Multiplicity_Det1;
  TH1 *hbPlas_Multiplicity_Det2;
  TH1 *hbPlas_Multiplicity_Det3;
  TH1 *hbPlas_lead_lead_ref_det[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_lead_lead_gated[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_ToT[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_trail_trail[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_Energy_Calib[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_SC41L_Digi_lead[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_SC41R_Digi_lead[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_Ge_Trig_lead[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_ToT_Sum[4];
  //TH1 *hbPlas_ToT_Sum_FibreCorr[4];

  TH1 *hbPlas_hit_pattern_det[4];
            
  ///Additional for twinpeaks
  TH1 *hbPlas_ToT_Sum_Slow[4];
  TH1 *hbPlas_ToT_Sum_Fast[4];
  TH1 *hbPlas_Lead_T_Slow[4][bPLASTIC_CHAN_PER_DET];
  TH1 *hbPlas_Lead_T_Fast[4][bPLASTIC_CHAN_PER_DET];         
  TH1 *hbPlas_ToT_det_Slow[4][bPLASTIC_CHAN_PER_DET];  
  TH1 *hbPlas_ToT_det_Fast[4][bPLASTIC_CHAN_PER_DET];  
  TH2 *hbPlas_ToT_Slow_vs_Fast_Det1;
  TH2 *hbPlas_ToT_Slow_vs_Fast_Det2;
	//NEW H.M.A.
  //TH2 *hFIMP_ToT_Correlation_Comb1;
  //TH2 *hFIMP_ToT_Correlation_Comb2;




  ///Germanium Histograms
  TH1 *hGe_Chan_E[Germanium_MAX_DETS][Germanium_CRYSTALS];
  TH1 *hGe_ERaw[Germanium_MAX_DETS][Germanium_CRYSTALS];
  TH1 *hGe_Chan_E_halfkev[Germanium_MAX_DETS][Germanium_CRYSTALS];
  //TH1 *hGe_Chan_E2;
  TH1 *hGe_Chan_Egate;
  TH1 *hGe_AddbackSum;
  TH1 *hGe_AddbackSum_halfkev;
  TH1 *hGe_dTaddback;
  TH1 *hGe_dTgammagamma;
  TH1 *hGe_CFdT_gammagamma;
  TH1 *hGe_Energy_GainMonitor;
  TH1 *hGe_Chan_E_vsTime;
  TH1 *hGe_SC41L;
  TH1 *hGe_SC41R;
  TH1 *hGe_SC41L_digi;
  TH1 *hGe_SC41R_digi;
  TH1 *hGe_Chan_Time_Diff[Germanium_MAX_DETS][Germanium_CRYSTALS];
  TH1 *hGe_Chan_Time_Diff_CF[Germanium_MAX_DETS][Germanium_CRYSTALS];
  //TH1 *hGe_Time_Diff_vs_Energy[32];
  TH1 *hGe_ESum;
  TH1 *hGe_Mult;
  TH1 *hGe_ESum_halfkev;
  TH1 *hGe_ESum_largerange_OF;
  TH1 *hGe_ESum_largerange_PU;
  TH1 *hGe_Hit_Pat;
  TH2 *hGe_Chan_E_Mat;
  TH2 *hGe_Chan_E_vsCrys;
  TH2 *hGe_Time_Diff_vs_Energy_sum;
  TH2 *hGe_Chan_E_vsDet;
  TH2 *hGe_MultvsdT;
  //TH1 *hGe_Chan_E_gate[32];
  TH1 *hGe_ge_time_difference_gg;


  //Finger Histograms
  TH1 *hFING_Hit_Pat;
  TH1 *hFING_Multiplicity;
  TH1 *hFING_ToT_Strip[52];
  TH1 *hFING_MaxToT_Strip[52];
  TH2 *hFING_ToT_v_Strip;
  TH2 *hFING_MaxToT_v_Strip;

  TH1 *hFING_ToT_PMT[52];
  TH2 *hFING_ToT_v_PMT;
  TH1 *hFING_lead_lead[52];
  TH1 *hFING_trail_trail[52];
  TH1 *hFING_Sc41lead_leadmaxtot[52];
  TH2* hFing_ToTRatio_v_Strip;
  TH2* hFing_MaxToTRatio_v_Strip;

  TH1 *hFING_SC41_lead_lead;
  TH1 *hFING_SC41_trail_trail;
  TH1 *hFING_SC41_tot;
  //TH1 *hSC41FatPlas;

  private :
  ///TEMP
  bool fired_det1, fired_det2;
  int bPlas_tot_hits[4][bPLASTIC_CHAN_PER_DET] ;

  bool same_ring_exclusion; // Read from General Setup File
  bool output_position_matrix; // Read from General Setup File

  ClassDef(EventAnlProc, 1)
};
#endif //EVENTANLPROCESSOR_H
