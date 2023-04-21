// $Id: EventAnlProc.cxx 754 2011-05-18 11:04:52Z adamczew $
// Adapted for DESPEC by A.K.Mistry 2022
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum fuer Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
    //-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

//Uncomment this to align the AIDA ASICs with a pulser
//Only needed if the ASICs didn't align properly
//#define AIDA_PULSER_ALIGN
#include "EventAnlProc.h"

#include <cstdlib>
#include <math.h>
#include <unordered_map>

#include "TH1.h"
#include "TH2.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"

#include "TGo4WinCond.h"
#include "TGo4Analysis.h"
#include "TGo4Picture.h"

#include "EventAnlStore.h"
#include "EventUnpackStore.h"

#include "TAidaConfiguration.h"
#include "TFRSParameter.h"
#include "DESPECAnalysis.h"
#include "FRS_Detector_System.h"

#define ABS(x)  ((x)>=0 ? (x):-(x))  // absolute_value(x)
#define ZERO_ARRAY(x) memset(x, 0, sizeof(x)) //reset arrays to 0
//-----------------------------------------------------------
EventAnlProc::EventAnlProc() :
  TGo4EventProcessor(){
  }
//-----------------------------------------------------------
EventAnlProc::EventAnlProc(const char* name) :
TGo4EventProcessor(name){
  //Clear up for AIDA
  implantEvents = 0;
  goodImplantEvents = 0;
  stoppedEvents = 0;
  decayEvents = 0;
  pulserEvents = 0;
  nonsenseEvents = 0;

  FRS_spill=0;



  cout << "**** EventAnlProc: Create" << endl;

  //checkTAMEXorVME();

  fCal = (CalibParameter*) GetParameter("CalibPar");
  fCorrel = (CorrelParameter*) GetParameter("CorrelPar");

  DESPECAnalysis* an = dynamic_cast<DESPECAnalysis*> (TGo4Analysis::Instance());
  frs_id = dynamic_cast<TIDParameter*> (an->GetParameter("IDPar"));
     

  // read_setup_parameters();
  get_used_systems();
  FRS_Gates();
  FRS_GainMatching();
}
//-----------------------------------------------------------
EventAnlProc::~EventAnlProc(){
  cout << "**** EventAnlProc: Delete" << endl;
}

void EventAnlProc::UserPostLoop(){
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
  if (!conf->ShowStats()) return;
  std::cout << "AIDA Analysis Statistics" << std::endl;
  std::cout << "Incoming Implant Events: " << implantEvents << std::endl;
  std::cout << "Good Implant Events    : " << goodImplantEvents << " (" << (100. * goodImplantEvents / implantEvents) << "%)" << std::endl;
  std::cout << "Stopped Implant Events : " << stoppedEvents << " (" << (100. * stoppedEvents / implantEvents) << "%)" << std::endl;
}

/**----------------------------------------------------------------------------------------------**/
/**-------------------------------------EVENT BUILDER--------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/

Bool_t EventAnlProc::BuildEvent(TGo4EventElement* dest){
  for(int i=0; i<7; i++){
    PrcID_Conv[i]=-1;
  }

  Bool_t isValid=kFALSE; // validity of output event

  EventUnpackStore* pInput  = (EventUnpackStore*) GetInputEvent();
  EventAnlStore* pOutput = (EventAnlStore*) dest;

  //pAida.Implants.clear();
  //pAida.Decays.clear();

  if((pInput==0) || !pInput->IsValid()){ // input invalid
    pOutput->SetValid(isValid); // invalid
    return isValid; // must be same is for SetValid
  }
  isValid=kTRUE;

  //general inputs from the unpacker
  event_number = pInput->fevent_number;
  pOutput->pOnSpill = FRS_spill;
  pOutput->pTrigger = pInput->fTrigger;
  pOutput->pEvent_Number = event_number;

  for (int i = 0; i<7; i++){
    if(pInput->fProcID[i]>-1){
      PrcID_Conv[i] = pInput->fProcID[i];
      pOutput->pPrcID_Conv[i] = pInput->fProcID[i];
    }
  }

  static bool create =false;
  //Create histograms
  if (!create){
    Make_WR_Histos();
    if(Used_Systems[0])  Make_FRS_Histos();
    if(Used_Systems[1])  Make_Aida_Histos();
    if(Used_Systems[2] && bPLASTIC_TWINPEAKS==1)  Make_Plastic_Twinpeaks_Histos();
    if(Used_Systems[5]) Make_Germanium_Histos();
  }

  create = true;
  Process_WR_Histos(pInput);
  //Fat_TimeCorrection(pInput);
  /**Now extract the data from the stored Unpacker array (root tree)**/

///--------------------------------------/**FRS Input**/------------------------------------------///

  if (PrcID_Conv[0]==0 && Used_Systems[0]==1){
    for(int i=0; i<10; i++){
      FRS_AoQ_corr_mhtdc[i]=0;
      FRS_z_mhtdc[i]=0;
      FRS_z2_mhtdc[i]=0;
      FRS_AoQ_mhtdc[i]=0;
    }
    pOutput->pFRS_WR = pInput->fFRS_WR;

    ///MUSIC
    for(int i =0; i<2; ++i){
      FRS_dE[i] = pInput->fFRS_Music_dE[i];
      FRS_dE_cor[i] = pInput->fFRS_Music_dE_corr[i];
    }

    ///SCI
    for(int l=0;l<12;++l){
      FRS_sci_l[l] = pInput->fFRS_sci_l[l];
      //if(pInput->fFRS_sci_l[l]!=0)cout<<"pInput->fFRS_sci_l[l] " << pInput->fFRS_sci_l[l] <<" l " << l << endl;
      FRS_sci_r[l] = pInput->fFRS_sci_r[l];
      FRS_sci_e[l] = pInput->fFRS_sci_e[l];
      FRS_sci_tx[l] = pInput->fFRS_sci_tx[l];
      FRS_sci_x[l] = pInput->fFRS_sci_x[l];
    }

    ///VFTX
    for(int i = 0; i < VFTX_MAX_HITS; i++){
      pOutput->pTRaw_vftx_21l[i] = pInput->fTRaw_vftx_21l[i];
      pOutput->pTRaw_vftx_21r[i] = pInput->fTRaw_vftx_21r[i];
      pOutput->pTRaw_vftx_22l[i] = pInput->fTRaw_vftx_22l[i];
      pOutput->pTRaw_vftx_22r[i] = pInput->fTRaw_vftx_22r[i];
      pOutput->pTRaw_vftx_41l[i] = pInput->fTRaw_vftx_41l[i];
      pOutput->pTRaw_vftx_41r[i] = pInput->fTRaw_vftx_41r[i];
      pOutput->pTRaw_vftx_42l[i] = pInput->fTRaw_vftx_42l[i];
      pOutput->pTRaw_vftx_42r[i] = pInput->fTRaw_vftx_42r[i];
    }
        
    ///TPC
    for(int i = 0; i < 7; i++){
      pOutput->pFRS_tpc_x[i] = pInput->fFRS_TPC_x[i];
      pOutput->pFRS_tpc_y[i] = pInput->fFRS_TPC_y[i];
    }
        
    ///Scaler
    for (int i = 0; i < 64; i++){
      pOutput->pFRS_scaler[i] = pInput->fFRS_scaler[i];
      pOutput->pFRS_scaler_delta[i] = pInput->fFRS_scaler_delta[i];
    }
    if (pOutput->pFRS_scaler_delta[8] > 0){
      FRS_spill = true;
      pOutput->pOnSpill = true;
	    //cout<<"OnSpill " <<pOutput->pOnSpill << endl;
    }
    if (pOutput->pFRS_scaler_delta[9] > 0){
      FRS_spill = false;
      pOutput->pOnSpill = false;
	    // cout<<"Offspill " <<pOutput->pOnSpill << endl;
    }

    ///MUSIC 
    for(int i =0; i<2; ++i){
      pOutput->pFRS_Music_dE[i] = pInput->fFRS_Music_dE[i];
    }
                
    ///SCI TOF
    //FRS_sci_tofll2 = pInput->fFRS_sci_tofll2;
    //FRS_sci_tofll3 = pInput->fFRS_sci_tofll3;
    FRS_sci_tof2 = pInput->fFRS_sci_tof2;
    //FRS_sci_tofrr2 = pInput->fFRS_sci_tofrr2;
    //FRS_sci_tofrr3 = pInput->fFRS_sci_tofrr3;
    //FRS_sci_tof3 = pInput->fFRS_sci_tof3;

    ///ID 2 4
    FRS_ID_x2 = pInput->fFRS_ID_x2;
    FRS_ID_y2 = pInput->fFRS_ID_y2;
    FRS_ID_a2 = pInput->fFRS_ID_a2;
    FRS_ID_b2 = pInput->fFRS_ID_b2;
    FRS_ID_x4 = pInput->fFRS_ID_x4;
    FRS_ID_y4 = pInput->fFRS_ID_y4;
    FRS_ID_a4 = pInput->fFRS_ID_a4;
    FRS_ID_b4 = pInput->fFRS_ID_b4;

    ///SCI dT
    //FRS_sci_dt_21l_21r = pInput->fFRS_sci_dt_21l_21r;
    //FRS_sci_dt_41l_41r = pInput->fFRS_sci_dt_41l_41r;
    //FRS_sci_dt_42l_42r = pInput->fFRS_sci_dt_42l_42r;
    //FRS_sci_dt_43l_43r = pInput->fFRS_sci_dt_43l_43r;
    //FRS_sci_dt_21l_41l = pInput->fFRS_sci_dt_21l_41l;
    //FRS_sci_dt_21r_41r = pInput->fFRS_sci_dt_21r_41r;
    //FRS_sci_dt_21l_42l = pInput->fFRS_sci_dt_21l_42l;
    //FRS_sci_dt_21r_42r = pInput->fFRS_sci_dt_21r_42r;

    ///ID Beta Rho
    for(int i =0; i<2; ++i){
      FRS_ID_brho[i] = pInput->fFRS_ID_brho[i];
      FRS_ID_rho[i] = pInput->fFRS_ID_rho[i];
    }
        
    //FRS_beta3 = pInput->fFRS_beta3;
    //FRS_gamma  = pInput->fFRS_gamma;

    ///ID Z AoQ
    FRS_AoQ = pInput->fFRS_AoQ;
    FRS_AoQ_corr = pInput->fFRS_AoQ_corr;
    FRS_z = pInput->fFRS_z;
    FRS_z2 = pInput->fFRS_z2;
    cout << "Initial Z1: " << FRS_z << endl;
    cout << "Initial Z2: " << FRS_z2 << endl;

    FRS_time_mins = 0;
    if(pOutput->pFRS_WR>0) FRS_time_mins =(pOutput->pFRS_WR/60E9);//-26900000;

    //Gainmatch Z
    for(int i=0; i<Z_Shift_array; i++){
      if(FRS_time_mins >=FRS_WR_a[i] && FRS_time_mins < FRS_WR_b[i]){
        cout << "Before Z1: " << FRS_z << endl;
        FRS_z = FRS_z + Z1_shift_value[i];
        cout << "After Z1: " << FRS_z << endl;
        cout << "Before Z2: " << FRS_z2 << endl;
        FRS_z2 = FRS_z2 + Z2_shift_value[i];
        cout << "After Z2: " << FRS_z2 << endl;
      }
    }
    


    //x4size=(sizeof(FRS_ID_x4)/sizeof(*FRS_ID_x4));
		//x4size = ( ( sizeof FRS_ID_x4 ) / ( sizeof FRS_ID_x4[0] ) );
    
		//GEEBART Z CORRECTIONS
		int z1 = 74;
		double a0 =0.000000000203148044;
		double a1 =-0.00000000217563685;
		double a2=-0.00000133322001;
		double a3=-0.00003582538447;
		double a4=0.000724027053;
		double a5=74.406645;
		
		FRS_z= z1*FRS_z/(a0*pow(FRS_ID_x4,5)+a1*pow(FRS_ID_x4,4)+a2*pow(FRS_ID_x4,3)+a3*pow(FRS_ID_x4,2)+a4*FRS_ID_x4+a5);
		
		int z2 = 74;
		double b0 =0.00000000012958873;
		double b1 =-0.0000000058809294;
		double b2=-0.00000072300288;
		double b3=-0.000029715951;
		double b4=0.000096556729;
		double b5=74.481789;
		
		FRS_z2= z2*FRS_z2/(b0*pow(FRS_ID_x4,5)+b1*pow(FRS_ID_x4,4)+b2*pow(FRS_ID_x4,3)+b3*pow(FRS_ID_x4,2)+b4*FRS_ID_x4+b5);


    FRS_dEdeg = pInput->fFRS_dEdeg;
    FRS_dEdegoQ = pInput->fFRS_dEdegoQ;
    FRS_beta = pInput->fFRS_beta;

    for (int i=0; i<10; i++){
      FRS_AoQ_mhtdc[i] = pInput->fFRS_AoQ_mhtdc[i];  
      //cout<<"pInput->fFRS_AoQ_mhtdc[i]" << pInput->fFRS_AoQ_mhtdc[i]   << endl;
      FRS_AoQ_corr_mhtdc[i] = pInput->fFRS_AoQ_corr_mhtdc[i];
      FRS_z_mhtdc[i] = pInput->fFRS_z_mhtdc[i];
      FRS_z2_mhtdc[i] = pInput->fFRS_z2_mhtdc[i];
      FRS_dEdeg_mhtdc[i] = pInput->fFRS_dEdeg_mhtdc[i];
      FRS_dEdegoQ_mhtdc[i] = pInput->fFRS_dEdegoQ_mhtdc[i];
      FRS_beta_mhtdc[i] = pInput->fFRS_beta_mhtdc[i];        
      //fFRS_tof4121_mhtdc[i]
    }
     
    Process_FRS_Histos(pInput, pOutput);
  }


///--------------------------------------/**AIDA Input**/------------------------------------------///
  if (Used_Systems[1]&&  PrcID_Conv[1]==1){
    ProcessAida(pInput, pOutput);
    Aida_Fired = 0;
    Aida_Fired = pInput->fAIDAHits;
    pOutput->pAidaScalers = pInput->fAidaScalers;
  }
   
///--------------------------------------/**bPlastic TAMEX Input**/------------------------------------------///
// if (PrcID_Conv[2] ==2 && Used_Systems[2]==1 && bPLASTIC_TWINPEAKS==0){
// 
//  for(int i=0; i<bPLASTIC_TAMEX_HITS; i++){
//    bPlas_RefCh0_Det1[i] =0;
//    bPlas_RefCh0_Det2[i] =0;
//    //bPlas_RefCh0_Det3[i] =0;
//  }
//
//  for (int j = 0; j < bPLASTIC_CHAN_PER_DET; j++){
//    for(int k=0; k<bPLASTIC_TAMEX_HITS;k++){
//      lead_lead_bplas_Ref1[j][k]=0;
//      lead_lead_bplas_Ref2[j][k]=0;
//      lead_lead_fat_Ref0[j][k]=0;
//    }
//  }
// 
//  for (int i=1; i<4; i++){
//    for (int j = 0; j < bPLASTIC_CHAN_PER_DET; j++){
//      for(int k=0; k<bPLASTIC_TAMEX_HITS;k++){
//        SC41L_ANA_lead_bPlas[i][j][k] = 0;
//        SC41R_ANA_lead_bPlas[i][j][k] = 0;
//        SC41L_DIG_lead_bPlas[i][j][k] = 0;
//        SC41R_DIG_lead_bPlas[i][j][k] = 0;
//      }
//    }
//  }
// 
//  pOutput->pbPLAS_WR = pInput->fbPlas_WR;
// 
//  bPlas_TAM_FATTAM = pInput->fbPlas_Lead_PMT[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_FATTAMEX][0];
//  bPlas_TAM_FATVME = pInput->fbPlas_Lead_PMT[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_FATVME][0];
//  bPlas_TAM_SC41L_DIG = pInput->fbPlas_Lead_PMT[bPLASTIC_ADDITIONAL_CH_MOD][SC41L_bPLASTIC][0];
//  bPlas_TAM_SC41R_DIG = pInput->fbPlas_Lead_PMT[bPLASTIC_ADDITIONAL_CH_MOD][SC41R_bPLASTIC][0];
// 
//  //bPlas_AND_Coinc[j] = pInput->fFat_Lead_PMT[9][j];
// 
//  for(int i=1; i<4; i++){ ///Detector number
//    for (int j = 0; j <bPLASTIC_CHAN_PER_DET ; j++){  ///Channel number
//      for(int k=0; k<bPLASTIC_TAMEX_HITS; k++){
//        //Fat_RefCh[j] = pInput->fFat_Lead_PMT[1][j];
//        bPlas_RefCh0_Det1[k] = pInput->fbPlas_Lead_PMT[1][bPlastRefCh_Det1][k];
//        bPlas_RefCh0_Det2[k] = pInput->fbPlas_Lead_PMT[2][bPlastRefCh_Det2][k];
//        //bPlas_RefCh0_Det3[k] = pInput->fbPlas_Lead_PMT[3][bPlastRefCh_Det3][k];
//      }
//    }
//   }
//   Process_Plastic_Tamex_Histos(pInput,pOutput);
// 
//}

///--------------------------------------/**bPlastic TwinPeaks TAMEX Input**/------------------------------------------///
  if (PrcID_Conv[2] ==2 && Used_Systems[2]==1 && bPLASTIC_TWINPEAKS==1){

      for(int i=0; i<bPLASTIC_TAMEX_HITS; i++){
         bPlas_RefCh0_Det1[i] =0;
         bPlas_RefCh0_Det2[i] =0;
        // bPlas_RefCh0_Det3[i] =0;
      }
      for (int j = 0; j < bPLASTIC_CHAN_PER_DET; j++){   
        for(int k=0; k<bPLASTIC_TAMEX_HITS;k++){    
            lead_lead_bplas_Ref1[j][k]=0;  
            lead_lead_bplas_Ref2[j][k]=0;  
            //lead_lead_fat_Ref0[j][k]=0;   
        }
      }
                 
      for (int i=1; i<4; i++){
        for (int j = 0; j < bPLASTIC_CHAN_PER_DET; j++){   
          for(int k=0; k<bPLASTIC_TAMEX_HITS;k++){    
            SC41L_ANA_lead_bPlas[i][j][k] = 0;    
            SC41R_ANA_lead_bPlas[i][j][k] = 0;    
            SC41L_DIG_lead_bPlas[i][j][k] = 0;    
            SC41R_DIG_lead_bPlas[i][j][k] = 0;    
          }
        }   
      }
  
      pOutput->pbPLAS_WR = pInput->fbPlas_WR;
                
      //bPlas_TAM_FATTAM = pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_FATTAMEX][0];
      //bPlas_TAM_FATVME = pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_FATVME][0];
      bPlas_TAM_SC41L_DIG = pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][SC41L_bPLASTIC][0];
      bPlas_TAM_SC41R_DIG = pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][SC41R_bPLASTIC][0];
      bPlas_TAM_Ge_TRIG = pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_Ge_TRIGGER][0];

      //bPlas_AND_Coinc[j] = pInput->fFat_Lead_PMT[9][j];
      
      for(int i=1; i<pInput->fbPlasDetNum_Fast; i++){ ///Detector number
        for (int j = 0; j < pInput->fbPlas_FastChan[i]; j++){  ///Channel number      
          for(int k=0; k<= pInput->fbPlast_Fast_Lead_N[i][j]; k++){ 
            //Fat_RefCh[j] = pInput->fFat_Lead_PMT[1][j]; 
            bPlas_RefCh0_Det1[k] = pInput->fbPlast_Fast_Lead[1][bPlastRefCh_Det1][k];
            bPlas_RefCh0_Det2[k] = pInput->fbPlast_Fast_Lead[2][bPlastRefCh_Det2][k];
            //bPlas_RefCh0_Det3[k] = pInput->fbPlas_Lead_PMT[3][bPlastRefCh_Det3][k];
          }      
        }
      Process_Plastic_Twinpeaks_Histos(pInput,pOutput); 
      }
  }
///--------------------------------------/**Germanium Input**/------------------------------------------///
  GeFired = -1;
  RefTGe=0;
  //Ge_WR = 0;
  for(int g = 0; g<Germanium_MAX_CHANNELS; g++){
    GeDet[g] = -1;
    GeCrys[g] = -1;
    GeE[g] = -1;
    GeE_Cal[g] = -1;
    GeT[g] =-1;
    GePileUp[g] = false;
    GeOverFlow[g] = false;
    Ge_Talign[g]=0;
  }

  if ( PrcID_Conv[5]==5 && Used_Systems[5]==1){
    GeFired =  pInput->fGe_fired;
    //GePileup = pInput->fGe_Pileup;
    Ge_WR = pInput->fGe_WR;

    for (int i = 0; i<GeFired; i++){
      GeDet[i] = pInput->fGe_Detector[i];
      GeCrys[i] = pInput->fGe_Crystal[i];
      GePileUp[i] = pInput->fGe_Pileup[i];
      GeOverFlow[i] = pInput->fGe_Overflow[i];
      GeE[i] = pInput->fGe_E[i];
      GeEventT[i]=pInput->fGe_Event_T[i];
      GeT[i] = pInput->fGe_T[i];
      GeCF_T[i] = pInput->fGe_cfT[i];

      // Maybe keep this for auto calibration program purposes
      int id = GeDet[i] * Germanium_CRYSTALS + GeCrys[i];
	    if(id<Germanium_MAX_CHANNELS){
        GeE_Cal[i] = (fCal->AGe[id]* pow( GeE[i],2) + fCal->BGe[id]*  GeE[i] + fCal->CGe[id]);
        //  printf("detIDGe %02d  %lE %lE %lE   calculated Energy %lf\n", id, fCal->AGe[id], fCal->BGe[id], fCal->CGe[id], GeE_Cal[i]  );
      }
      // int id_Ge = det * Germanium_CRYSTALS +  crys;

      if(id==1){
        RefTGe=GeT[i];
        RefCFDGe=GeCF_T[i];
      }

      Ge_Talign[i] = (GeT[i]-fCal->Ge_T_align_par[id]);

      pOutput->pGe_T_Aligned[GeDet[i]][GeCrys[i]]=Ge_Talign[i];
      //Do the CFD time alignment (detector vs all)
      Ge_cfd_Talign[i] = GeCF_T[i]+fCal->Ge_cfd_align_par[id];
      pOutput->pGe_CF_T_Aligned[GeDet[i]][GeCrys[i]] = Ge_cfd_Talign[i];
    }

    Process_Germanium_Histos(pOutput);
  }

///--------------------------------------/**FINGER Input**/------------------------------------------///

  pOutput->SetValid(isValid);
  return isValid;

}  
//End of BuildEvent




/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------WHITE RABBIT--------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_WR_Histos(){
  //Added on 01/03 from MP to check the subsystems' deadtime

  //hFat_deadtime = MakeTH1('I',"WR/DeadTime/Fat_deadtime","Dead Time Fatima VME", 500, 0, 500,"WR dT(Fatima VME)[us]", "Counts");
  //hFatTAM_deadtime = MakeTH1('I',"WR/DeadTime/FatTAM_deadtime","Dead Time Fatima TAMEX", 500, 0, 500,"WR dT(Fatima TAMEX)[us]", "Counts");

  hGe_deadtime = MakeTH1('I',"WR/DeadTime/HPGe_deadtime","Dead Time Germanium", 500, 0, 500,"WR dT(Germanium)[us]", "Counts");

  hbPlast_deadtime = MakeTH1('I',"WR/DeadTime/bPlast_deadtime","Dead Time bPlastic", 500, 0, 500,"WR dT(bPlastic)[us]", "Counts");

  hFRS_deadtime = MakeTH1('I',"WR/DeadTime/FRS_deadtime","Dead Time FRS", 500, 0, 500,"WR dT(FRS)[us]", "Counts");

  hbPlas_Ge_WRdT = MakeTH1('I',"WR/bPlast-Germanium_dT","White Rabbit bPlast-Germanium",10000,-50000,50000,"WR dT(bPlast-Germanium)[ns]", "Counts");

  hFRS_Ge_WRdT =  MakeTH1('I',"WR/Germanium-FRS_dT","White Rabbit FRS WR -Germanium WR ",10000,-50000,50000,"WR dT(Germanium-FRS)[ns]", "Counts");

  hFRS_bPlas_WRdT = MakeTH1('I',"WR/bPlast-FRS_dT","White Rabbit FRS_bPlas",10000,-50000,50000,"WR dT(bPlast-FRS)[ns]", "Counts");
}
 
void EventAnlProc::Process_WR_Histos(EventUnpackStore* pInput){
  /// Germanium DEAD TIME
  if (pInput->fGe_WR > 0){
    if (lastGeWR == 0){
      lastGeWR = pInput->fGe_WR;
    } 
    else{
      hGe_deadtime->Fill((long long)(pInput->fGe_WR - lastGeWR)/1000);
      lastGeWR = pInput->fGe_WR;
    }
  }

  /// bPlastic DEAD TIME
  if (pInput->fbPlas_WR > 0){
    if (lastbPlastWR == 0){
      lastbPlastWR = pInput->fbPlas_WR;
    } 
    else{
      hbPlast_deadtime->Fill((long long)(pInput->fbPlas_WR - lastbPlastWR)/1000);
      lastbPlastWR = pInput->fbPlas_WR;
    }
  }

  /// FRS DEAD TIME
    if (pInput->fFRS_WR > 0){
      if (lastFRSWR == 0){
        lastFRSWR = pInput->fFRS_WR;
      }
      else{
        hFRS_deadtime->Fill((long long)(pInput->fFRS_WR - lastFRSWR)/1000);
        lastFRSWR = pInput->fFRS_WR;
      }
    }

  ///Aida deadtime in Correl proc

  ///WR dT(bPlas-Germanium)
  if(pInput->fGe_WR>0 && pInput->fbPlas_WR>0)hbPlas_Ge_WRdT->Fill(pInput->fbPlas_WR - pInput->fGe_WR);

  ///WR dT(Germanium - FRS)
  if(pInput->fGe_WR>0 && pInput->fFRS_WR>0)hFRS_Ge_WRdT->Fill(pInput->fGe_WR - pInput->fFRS_WR );

  ///WR dT(bPlast- FRS)
  if(pInput->fbPlas_WR>0 && pInput->fFRS_WR>0)hFRS_bPlas_WRdT->Fill(pInput->fbPlas_WR-pInput->fFRS_WR);

}

/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------    FRS    ---------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_FRS_Histos(){

///--------------------------------------/**TAC HISTOGRAMS**/------------------------------------------///

  hID_Z1_corr = MakeTH1('D',"FRS/ID/ID_Z1_corr","ID_Z1_corr",1000,70,90,"Z1_corr s2-s4");

  hID_Z2_corr = MakeTH1('D',"FRS/ID/ID_Z2_corr","ID_Z2_corr",1000,70,90,"Z2_corr s2-s4");
  
  hID_Z1_vs_T = MakeTH2('D',"FRS/ID/ID_Z1_Time", "Z1 vs Time",2000,27538600,27546200, 1500,frs_id->min_z_plot,frs_id->max_z_plot,"Time (/10 mins)", "Z1 s2-s4");

  hID_Z2_vs_T = MakeTH2('D',"FRS/ID/ID_Z2_Time", "Z2 vs Time",2000,27538600,27546200, 1500,frs_id->min_z_plot,frs_id->max_z_plot,"Time (/10 mins)", "Z2 s2-s4");

  hID_AoQ_vs_T = MakeTH2('D',"FRS/ID/ID_AoQ_Time", "AoQ vs Time",2000,27450000,27560000,1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot,"Time (/10 mins)", "AoQ mhtdc s2-s4");
  
  hID_Z_AoQ = MakeTH2('D',"FRS/ID/ID_Z1_AoQ", "Z1 vs A/Q",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z1 s2-s4");
  
  hID_Z_AoQ_corr = MakeTH2('D',"FRS/ID/ID_Z1_AoQ_corr", "",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot, "A/Q s2-s4", "Z s2-s4");
  
  hID_Z_AoQ_zsame = MakeTH2('D',"FRS/ID/ID_Z1_AoQ_zsame","Z1 vs Z2: mod(Z1-Z2)<0.4", 1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"Z1==Z2 A/Q s2-s4", "Z s2-s4");
  
  hID_Z_Z2 = MakeTH2('D',"FRS/ID/ID_Z1_Z2","Z1 vs. Z2", 2000,frs_id->min_z_plot,frs_id->max_z_plot, 400,frs_id->min_z_plot,frs_id->max_z_plot,"Z1", "Z2");

  hID_x2AoQ = MakeTH2('D',"FRS/ID/ID_x2AoQ", "X2 vs A/Q",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 300,-100.,100.,"A/Q s2-s4", "X at S2 [mm]");

  hID_x2AoQ_corr = MakeTH2('D',"FRS/ID/ID_x2AoQ_corr", "X2 vs A/Q_corr",1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q_corr s2-s4", "X at S2 [mm]");
  
  hID_x2AoQ_zsame = MakeTH2('D',"FRS/ID/ID_x2AoQ_zsame","X2 vs A/Q Zsame", 1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 300,-100.,100.,"A/Q s2-s4", "X at S2 [mm]");

  hID_x4AoQ = MakeTH2('D',"FRS/ID/ID_x4AoQ", "X4 vs A/Q",1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 300,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

  hID_x4AoQ_corr = MakeTH2('D',"FRS/ID/ID_x4AoQ_corr", "X4 vs A/Q_corr",1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q_corr s2-s4", "X at S4 [mm]");

  hID_x4AoQ_zsame = MakeTH2('D',"FRS/ID/ID_x4AoQ_zsame","X4 vs A/Q Zsame", 1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 300,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");
  
  hID_a2AoQ = MakeTH2('D',"FRS/ID/ID_Angle_s2_AoQ", "A/Q vs Angle",500,frs_id->min_aoq_plot,frs_id->max_aoq_plot,500,-25,25,"AoQ s2-s4", "Angle s2 (mrad)");
  
  hID_a4AoQ = MakeTH2('D',"FRS/ID/ID_Angle_s4_AoQ", "A/Q vs Angle",500,frs_id->min_aoq_plot,frs_id->max_aoq_plot,500,-25,25,"AoQ s2-s4", "Angle s4 (mrad)");
  
  hID_dEdegoQ_Z = MakeTH2('D',"FRS/ID/ID_dEdegoQZ1","dE in s2 degrader/Q vs. Z1", 1000,frs_id->min_z_plot,frs_id->max_z_plot, 1000, 0,1.2, "Z from MUSIC41", "dE(S2deg)/Q [a.u.]");

  hID_dEdeg_Z   = MakeTH2('D',"FRS/ID/ID_dEdegZ1" ,"dE in s2 degrader vs. Z1" , 1000,frs_id->min_z_plot,frs_id->max_z_plot, 1000, 10.,100., "Z from MUSIC41", "dE(S2deg) [a.u.]");
  
  hID_Z_dE2 = MakeTH2('D',"FRS/ID/ID_Z1_dE2","ID_Z1_dE2", 400,frs_id->min_z_plot,frs_id->max_z_plot, 250,0.,4000.,"Z1", "MUSIC2_dE");

  hID_Z_E_SC41 = MakeTH2('D',"FRS/ID/ID_Z1_E_SC41","E SC41 vs. Z1", 400,0,4000, 2000,frs_id->min_z_plot,frs_id->max_z_plot,"E SC41", "Z1");
  
  hID_Z_E_SC42 = MakeTH2('D',"FRS/ID/ID_Z1_E_SC42","E SC42 vs. Z1", 400,0,4000, 2000,frs_id->min_z_plot,frs_id->max_z_plot,"E SC42", "Z1");
  
  hID_x2z1 = MakeTH2('D',"FRS/ID/ID_x2Z1","ID_x2Z1", 400,frs_id->min_z_plot,frs_id->max_z_plot, 200,-100.,100., "Z1 s2-s4", "X at S2 [mm]");

  hID_x4z1 = MakeTH2('D',"FRS/ID/ID_x4Z1","ID_x4Z1", 400,frs_id->min_z_plot,frs_id->max_z_plot, 200,-100.,100., "Z1 s2-s4", "X at S4 [mm]");

  hID_x4z2 = MakeTH2('D',"FRS/ID/ID_x4Z2","ID_x4Z2", 400,frs_id->min_z_plot,frs_id->max_z_plot, 200,-100.,100., "Z2 s2-s4", "X at S4 [mm]"); 

  hID_E_x2 = MakeTH2('D',"FRS/ID/ID_dE_x2","ID_dE_x2", 200,-100.,100., 400,0.,4000., "X s2 [mm]", "Delta E");

  hID_E_x4 = MakeTH2('D',"FRS/ID/ID_dE_x4","ID_dE_x4", 200,-100.,100., 400,0.,4000., "X s4 [mm]", "Delta E");

  hID_x2a2 = MakeTH2('D',"FRS/ID/ID_x2_a2", "ID_x2_a2", 200, -100., 100., 200, -100., 100., "X s2 [mm]", "AngleX s2 [mrad]");

  hID_y2b2 = MakeTH2('D',"FRS/ID/ID_y2_b2", "ID_y2_b2", 200, -100., 100., 200, -100., 100., "Y s2 [mm]", "AngleY s2 [mrad]");

  hID_x4a4 = MakeTH2('D',"FRS/ID/ID_x4_a4", "ID_x4_a4", 200, -100., 100., 200, -100., 100., "X s4 [mm]", "AngleX s4 [mrad]");

  hID_y4b4 = MakeTH2('D',"FRS/ID/ID_y4_b4", "ID_y4_b4", 200, -100., 100., 200, -100., 100., "Y s4 [mm]", "AngleY s4 [mrad]");

  hID_x2x4 = MakeTH2('D',"FRS/ID/ID_x2_x4","ID_x2_x4",200,-200,100,200,-200,100,"x2 mm","x4 mm");

  hID_SC41dE_AoQ = MakeTH2('D',"FRS/ID/ID_SC41dE_AoQ","ID_SC41dE_AoQ", 1000,1.2,3.0, 1000,0.,4000.,"A/Q s2-s4", "SC41 dE");

  hID_dEToF = MakeTH2('D',"FRS/ID/ID_dEToF","ID_dEToF", 2000, 0.,70000.,400,0,4000, "tof S2-S4 Sci.Tof(2)", "Music_dE(1)");

  hID_Z_Sc21E = MakeTH2('D',"FRS/ID/ID_Z1_Sc21E","ID_Z1_Sc21E", 300,0,80.,400,0,1000.,"Z s2-s4", "sqrt(Sc21_L*sC21_R)");

  hID_dEBRho = MakeTH2('D',"FRS/ID/ID_dEBRho","ID_dEBRho", 400, 3.0 ,5.0, 400,frs_id->min_z_plot,frs_id->max_z_plot , "BRho", "dE"); //GEEBART

  hID_dBrho_dE4 = MakeTH2('D', "FRS/ID/ID_dBrho_dE4", "ID_dBrho_dE4", 400, 3.0, 5.0, 5000, 0., 4000., "dBrho", "MUSIC4_dE"); //GEEBART

  hID_BRho1v2 = MakeTH2('D',"FRS/ID/ID_BRho1v2","ID_BRho1v2", 1000, 13, 14, 1000, 9, 11, "BRho1", "BRho2"); //GEEBART
  
  //hTPC_X_AX_S4=MakeTH2('D',"FRS/TPC/S4_focus/S4focus_X_angleX","S4focus_X_angleX", 400,-100.,100., 250,-50.0,50.0,"X at S4 [mm] ","x angle [mrad] ");

  //hTPC_Y_AY_S4=MakeTH2('D',"FRS/TPC/S4_focus/S4focus_Y_angleY","S4_Y_angleY", 400,-100.,100., 250,-50.0,50.0,"Y at S4 [mm] ","y angle [mrad] ");

   
///--------------------------------------/**MHTDC HISTOGRAMS**/------------------------------------------///

  hID_Z_mhtdc_T = MakeTH2('D',"FRS/MHTDC/ID/ID_Z_Time_mhtdc", "Z mhtdc vs Time",1200,17000,29000,2000,frs_id->min_z_plot,frs_id->max_z_plot,"Time (/10 mins)", "Z mhtdc");

  hID_AoQ_mhtdc_T = MakeTH2('D',"FRS/MHTDC/ID/ID_AoQ_Corr_Time_mhtdc", "AoQ mhtdc vs Time",1200,17000,29000,1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot,"Time (/10 mins)", "AoQ mhtdc s2-s4");
  
  hID_Z_AoQ_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_Z1_AoQ_mhtdc", "Z1 vs A/Q (mhtdc)",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q mhtdc s2-s4", "Z1 mhtdc s2-s4");
  
  hID_Z_AoQ_corr_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_Z1_AoQ_corr_mhtdc", "Z1 vs A/Q corr (mhtdc)",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q corr mhtdc s2-s4", "Z1 mhtdc s2-s4");
  
  hID_Z_AoQ_zsame_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_Z1_AoQ_zsame_mhtdc","Z1 vs Z2: mod(Z1-Z2)<0.4 MHTDC", 1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"Z1==Z2 A/Q s2-s4", "Z s2-s4");

  hID_Z_Z2_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_Z1_Z2_mhtdc","Z1 vs. Z2 MHTDC", 2000,frs_id->min_z_plot,frs_id->max_z_plot, 1500,frs_id->min_z_plot,frs_id->max_z_plot,"Z1", "Z2");
   
  hID_x2AoQ_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_x2AoQ_mhtdc", "X2 vs A/Q MHTDC",1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S2 [mm]");

  hID_x2AoQ_zsame_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_x2AoQ_zsame_mhtdc","X2 vs A/Q Zsame MHTDC", 1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 300,-150.,100.,"A/Q s2-s4", "X at S2 [mm]");

  hID_x4AoQ_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_x4AoQ_mhtdc", "X4 vs A/Q MHTDC",1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");
  
  hID_x4AoQ_zsame_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_x4AoQ_zsame_mhtdc","X4 vs A/Q Zsame MHTDC", 1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 300,-150.,100.,"A/Q s2-s4", "X at S4 [mm]");

  hID_a2AoQ_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_Angle_s2_AoQ_mhtdc", "A/Q vs Angle",500,frs_id->min_aoq_plot,frs_id->max_aoq_plot,500,-25,25,"AoQ s2-s4", "Angle s2 (mrad)");
  
  hID_a4AoQ_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_Angle_s4_AoQ_mhtdc", "A/Q vs Angle",500,frs_id->min_aoq_plot,frs_id->max_aoq_plot,500,-25,25,"AoQ s2-s4", "Angle s2 (mrad)");

  hID_dEdegoQ_Z_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_dEdegoQZ1_mhtdc","dE in s2 degrader/Q vs. Z1 MHTDC", 1000,frs_id->min_z_plot,frs_id->max_z_plot, 1000, 0.1,0.8, "Z from MUSIC41 MHTDC", "dE(S2deg)/Q [a.u.]");
  
  hID_dEdeg_Z_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_dEdegoQZ1_mhtdc","dE in s2 degradervs. Z1 MHTDC", 1000,frs_id->min_z_plot,frs_id->max_z_plot, 1000, 10.,100., "Z from MUSIC41 MHTDC", "dE(S2deg) [a.u.]");
  
  hID_Z_dE2_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_Z1_dE2_mhtdc","ID_Z1_dE2_mhtdc", 400,frs_id->min_z_plot,frs_id->max_z_plot, 250,0.,4000.,"Z1", "MUSIC2_dE");

  hID_Z_Sc21E_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_Z1_Sc21E_mhtdc","ID_Z1_Sc21E_mhtdc", 300,0,25.,400,0,4000.,"Z s2-s4", "sqrt(Sc21_L*sC21_R)");
  
  hID_x2z_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_x2Z1_mhtdc","ID_x2Z1_mhtdc", 400,frs_id->min_z_plot,frs_id->max_z_plot, 200,-100.,100., "Z1 s2-s4", "X at S2 [mm]");

  hID_x4z_mhtdc = MakeTH2('D',"FRS/MHTDC/ID/ID_x4Z1_mhtdc","ID_x4Z1_mhtdc", 400,frs_id->min_z_plot,frs_id->max_z_plot, 200,-100.,100., "Z1 s2-s4", "X at S4 [mm]");
  
  //hID_Z_AoQ_mhtdc_first_hit = MakeTH2('D',"FRS/MHTDC/ID_Z1_AoQ_mhtdc_first_hit", "Z1 vs A/Q First hit (mhtdc)",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q mhtdc first hit s2-s4", "Z1 mhtdc first hit s2-s4");

  //hID_Z_AoQ_mhtdc_excluding_first_hit = MakeTH2('D',"FRS/MHTDC/ID_Z1_AoQ_mhtdc_exc_first_hit_elif", "Z1 vs A/Q Excl First hit (mhtdc)",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q mhtdc non-first hit s2-s4", "Z1 mhtdc non-first hit s2-s4");

  //hID_Z_AoQ_corr_mhtdc_first_hit = MakeTH2('D',"FRS/MHTDC//ID_Z1_AoQ_corr_mhtdc_first_hit", "Z1 vs A/Q corr First hit (mhtdc)",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q corr mhtdc first hit s2-s4", "Z1 mhtdc first hit s2-s4");

  //hID_Z1_AoQ_corr_mhtdc_exc_first_hit", "Z1 vs A/Q corr Excl First hit (mhtdc)",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q mhtdc corr non-first hit s2-s4", "Z1 mhtdc  non-first hit s2-s4");

  //hID_A2_AoQ_Corr_mhtdc= MakeTH2('D',"FRS/MHTDC/ID_AoQ_corr_a2_mhtdc", "A/Q corr vs angle",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,-100,100,"A/Q mhtdc corr ", "Angle S2 (mrad)");

  //hID_A4_AoQ_Corr_mhtdc= MakeTH2('D',"FRS/MHTDC/ID_AoQ_corr_a4_mhtdc", "A/Q corr vs angle",1500,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 1000,-100,100,"A/Q mhtdc corr ", "Angle S4 (mrad)");
  
  ///Now define the 2D Polygon PID Gated

  Float_t init_ID_x2AoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_x4AoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_x2AoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_x4AoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_Z_Z2[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_Z_Z2_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_Z_AoQ[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_Z_AoQ_mhtdc[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_dEvsBRho[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_dEvsZ[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];
  Float_t init_ID_dEvsZ_x4AoQ_p1[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];  
  Float_t init_ID_dEvsZ_x4AoQ_m1[MAX_FRS_GATE][MAX_FRS_PolyPoints][2];




  ///FRS gates initialisation
  for(int i=0; i<MAX_FRS_GATE; i++){
    for(int j=0; j<MAX_FRS_PolyPoints; j++){
      init_ID_x2AoQ[i][j][0] = XX2_AoQ[i][j];
      init_ID_x2AoQ[i][j][1] = YX2_AoQ[i][j];
      init_ID_x4AoQ[i][j][0] = XX4_AoQ[i][j];
      init_ID_x4AoQ[i][j][1] = YX4_AoQ[i][j];
      init_ID_Z_Z2[i][j][0] = X_ZZ2[i][j];
      init_ID_Z_Z2[i][j][1] = Y_ZZ2[i][j];
      init_ID_Z_AoQ[i][j][0] =X_ZAoQ[i][j];
      init_ID_Z_AoQ[i][j][1] =Y_ZAoQ[i][j];
      init_ID_dEvsBRho[i][j][0] =X_dEvsBRho[i][j];
      init_ID_dEvsBRho[i][j][1] =Y_dEvsBRho[i][j];
      init_ID_dEvsZ[i][j][0] =X_dEvsZ[i][j];
      init_ID_dEvsZ[i][j][1] =Y_dEvsZ[i][j];
      init_ID_dEvsZ_x4AoQ_p1[i][j][0] =X_dEvsZ_x4AoQ_p1[i][j];
      init_ID_dEvsZ_x4AoQ_p1[i][j][1] =Y_dEvsZ_x4AoQ_p1[i][j];
      init_ID_dEvsZ_x4AoQ_m1[i][j][0] =X_dEvsZ_x4AoQ_m1[i][j];
      init_ID_dEvsZ_x4AoQ_m1[i][j][1] =Y_dEvsZ_x4AoQ_m1[i][j];
    


      init_ID_Z_AoQ_mhtdc[i][j][0] =X_ZAoQ_mhtdc[i][j];   
      init_ID_Z_AoQ_mhtdc[i][j][1] =Y_ZAoQ_mhtdc[i][j];
      init_ID_Z_Z2_mhtdc[i][j][0] =X_ZZ2_mhtdc[i][j];   
      init_ID_Z_Z2_mhtdc[i][j][1] =Y_ZZ2_mhtdc[i][j];
      init_ID_x2AoQ_mhtdc[i][j][0] = XX2_AoQ_mhtdc[i][j]; 
      init_ID_x2AoQ_mhtdc[i][j][1] = YX2_AoQ_mhtdc[i][j];
      init_ID_x4AoQ_mhtdc[i][j][0] = XX4_AoQ_mhtdc[i][j];
      init_ID_x4AoQ_mhtdc[i][j][1] = YX4_AoQ_mhtdc[i][j];
    }
  }
  char name[50], title[100];
    
///--------------------------------------/**TAC GATES**/------------------------------------------///
  for(int i=0; i<MAX_FRS_GATE; i++){
    //Z vs AoQ
    sprintf(name,"cID_Z1_AoQ%d",i);
    cID_Z_AoQ[i] = MakePolyCond("FRS_Z1_AoQ_Gates", name, num_ID_Z_AoQ, init_ID_Z_AoQ[i], hID_Z_AoQ->GetName());

    hID_ZAoQ_ZAoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/Z1AoQ/Z1AoQ_Z1AoQGated/ID_ZAoQ_ZAoQgate%d",i), Form("ID_ZAoQ_ZAoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");
      
    hID_Z1Z2_ZAoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/Z1Z2/Z1Z2_Z1AoQGated/ID_Z1Z2_ZAoQgate%d",i), Form("ID_Z1Z2_ZAoQgate%d", i), FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"MUSIC Z1", "gate on Z1 vs AoQ: MUSIC Z2");
      
    hID_x2AoQ_Z1AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x2AoQ/x2AoQ_Z1AoQGated/ID_x2AoQ_Z1AoQgate%d",i), Form("ID_x2AoQ_Z1AoQgate%d",i),1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q_corr s2-s4", "X at S2 [mm]");
       
    hID_x4AoQ_Z1AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x4AoQ/x4AoQ_Z1AoQGated/ID_x4AoQ_Z1AoQgate%d",i), Form("ID_x4AoQ_Z1AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on Z1 vs AoQ: X at S4 [mm]");
            
    hID_a2_Z1AoQgate[i] = MakeTH1('F', Form("FRS/ID_Gated/a2/a2_Z1AoQGated/ID_a2_Z1AoQGated%d",i), Form("ID_a2_Z1AoQGate%d",i), 100, -1000, 1000, "Angle S2 (mrad)");
            
    hID_a4_Z1AoQgate[i] = MakeTH1('F', Form("FRS/ID_Gated/a4/a4_Z1AoQGated/ID_a4_Z1AoQGated%d",i), Form("ID_a4_Z1AoQGate%d",i), 100, -1000, 1000, "Angle S4 (mrad)");
      
    ///Z vs Z2
    sprintf(name,"cID_Z1_Z2_Gate%d",i);
    cID_Z_Z2gate[i] = MakePolyCond("FRS_Z1_Z2_Gates",name,num_ID_Z_Z2,init_ID_Z_Z2[i], hID_Z_Z2 ->GetName());
      
    hID_Z1_Z2gate[i] = MakeTH2('I', Form("FRS/ID_Gated/Z1Z2/Z1Z2Gated/ID_Z1_Z2gate%d",i), Form("ID_Z1_Z2gate%d",i),FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"Z1 s2-s4", "Z2 s2-s4");

    hID_ZAoQ_Z1Z2gate[i] = MakeTH2('I', Form("FRS/ID_Gated/Z1AoQ/Z1AoQ_Z1Z2Gated/ID_Z1AoQ_Z1Z2gate%d",i), Form("ID_Z1AoQ_Z1Z2gate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", " Z music2");
      
    hID_x2AoQ_Z1Z2gate[i] = MakeTH2('I', Form("FRS/ID_Gated/x2AoQ/x2AoQ_Z1Z2Gated/ID_x2AoQ_Z1Z2gate%d",i), Form("ID_x2AoQ_Z1Z2gate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q_corr s2-s4", "X at S2 [mm]");

    hID_x4AoQ_Z1Z2gate[i] = MakeTH2('I', Form("FRS/ID_Gated/x4AoQ/x4AoQ_Z1Z2Gated/ID_x4AoQ_Z1Z2gate%d",i), Form("ID_x4AoQ_Z1Z2gate%d",i),1000,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q_corr s2-s4", "X at S4 [mm]");
                   
    hID_a2_Z1Z2gate[i] = MakeTH1('F', Form("FRS/ID_Gated/a2/a2_Z1Z2gate/ID_a2_Z1Z2gate%d",i), Form("ID_a2_Z1Z2gate%d",i), 100, -1000, 1000, "Angle S2 (mrad)");
       
    hID_a4_Z1Z2gate[i] = MakeTH1('F', Form("FRS/ID_Gated/a4/a4_Z1Z2gate/ID_a4_Z1Z2gate%d",i), Form("ID_a4_Z1Z2gate%d",i), 100, -1000, 1000, "Angle S4 (mrad)");
            
    ///X2 vs AoQ
    sprintf(name,"cID_x2AoQ%d",i);
    cID_x2AoQ[i] = MakePolyCond("FRS_X2_Gates",name,num_ID_x2AoQ,init_ID_x2AoQ[i], hID_x2AoQ->GetName());
      
    hID_x2AoQ_x2AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x2AoQ/x2AoQ_x2AoQGated/ID_x2AoQ_x2AoQgate%d",i), Form("ID_x2AoQ_x2AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on X2 AoQ, ID X2: X at S2 [mm]");
      
    hID_Z1Z2_x2AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/Z1Z2/Z1Z2_x2AoQGated/ID_Z1Z2_x2AoQgate%d",i), Form("ID_Z1Z2_x2AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"MUSIC Z1", "gate on x2 vs AoQ: MUSIC Z2");

    hID_ZAoQ_Z1Z2x2AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/Z1AoQ/Z1AoQ_x2AoQGated/ID_Z1AoQ_Z1Z2x2AoQgate%d",i), Form("ID_Z1AoQ_Z1Z2x2AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", " gate on Z1 Z2 X2 AoQ Z music");
    
    hID_x2AoQ_Z1Z2x2AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x2AoQ/x2AoQ_Z1Z2x2AoQGated/ID_x2AoQ_Z1Z2x2AoQgate%d",i), Form("ID_x2AoQ_Z1Z2x2AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on Z1 Z2 X2 AoQ, ID X2: X at S2 [mm]");

    hID_x4AoQ_Z1Z2x2AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x4AoQ/x4AoQ_Z1Z2x2AoQGated/ID_x4AoQ_Z1Z2x2AoQgate%d",i), Form("ID_x4AoQ_Z1Z2x2AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on  Z1 Z2 X2 AoQ, ID X2: X at S4 [mm]");
      
    hID_a2_Z1Z2x2AoQgate[i] = MakeTH1('F', Form("FRS/ID_Gated/a2/a2_Z1Z2x2AoQGated/ID_a2_Z1Z2x2AoQgate%d",i), Form("ID_a2_Z1Z2x2AoQgate%d",i), 100, -1000, 1000, "gate on Z1 Z2 X2 AoQ Angle S2 (mrad)");

    hID_a4_Z1Z2x2AoQgate[i] = MakeTH1('F', Form("FRS/ID_Gated/a4/a4_Z1Z2x2AoQGated/ID_a4_Z1Z2x2AoQgate%d",i), Form("ID_a4_Z1Z2x2AoQgate%d",i), 100, -1000, 1000, "gate on Z1 Z2 X2 AoQ Angle S4 (mrad)");
      
    ///X4 vs AoQ
    sprintf(name,"cID_x4AoQ%d",i);
    cID_x4AoQ[i] = MakePolyCond("FRS_X4_Gates",name,num_ID_x4AoQ,init_ID_x4AoQ[i], hID_x4AoQ->GetName());
      
    hID_x4AoQ_x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x4AoQ/x4AoQ_x4AoQGated/ID_x4AoQ_x4AoQgate%d",i), Form("ID_x4AoQ_42AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on X4 AoQ, ID X4: X at S4 [mm]");
      
    hID_Z1Z2_x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/Z1Z2/Z1Z2_x4AoQGated/ID_Z1Z2_x4AoQgate%d",i), Form("ID_Z1Z2_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"MUSIC Z1", "gate on x4 vs AoQ: MUSIC Z2");
      
    hID_x4AoQ_Z1Z2x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x4AoQ/x4AoQ_Z1Z2x4AoQGated/ID_x4AoQ_Z1Z2x4AoQgate%d",i), Form("ID_x4AoQ_Z1Z2x4AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X4: X at S4 [mm]");
      
    hID_x2AoQ_Z1Z2x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x2AoQ/x2AoQ_Z1Z2x4AoQGated/ID_x2AoQ_Z1Z2x4AoQgate%d",i), Form("ID_x2AoQ_Z1Z2x4AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X4: X at S2 [mm]");
      
    hID_Z1Z2_Z1Z2x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/x4AoQ/Z1Z2_Z1Z2x2AoQGated/ID_Z1Z2_Z1Z2x4AoQgate%d",i), Form("ID_Z1Z2_Z1Z2x4AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"Z music41", "gate on FRS AoQ, ID X4: Z music42");
       
    hID_ZAoQ_Z1Z2x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/Z1AoQ/Z1AoQ_Z1Z2x4AoQGated/ID_Z1AoQ_Z1Z2x4AoQgate%d",i), Form("ID_Z1AoQ_Z1Z2x4AoQgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", " Z music41");
                   
    hID_a2_Z1Z2x4AoQgate[i] = MakeTH1('F', Form("FRS/ID_Gated/a2/a2_Z1Z2x4AoQGated/ID_a2_Z1Z2x4AoQgate%d",i), Form("ID_a2_Z1Z2x4AoQgate%d",i), 100, -1000, 1000, "gate on Z1 Z2 X4 AoQ Angle S2 (mrad)");
      
    hID_a4_Z1Z2x4AoQgate[i] = MakeTH1('F', Form("FRS/ID_Gated/a4/a4_Z1Z2x4AoQGated/ID_a4_Z1Z2x4AoQgate%d",i), Form("ID_a4_Z1Z2x4AoQgate%d",i), 100, -1000, 1000, "gate on Z1 Z2 X4 AoQ Angle S4 (mrad)");

    //dEvsBRho
    sprintf(name,"cID_dEvsBRho%d",i);
    cID_dEvsBRho[i] = MakePolyCond("FRS_dEvsBRho",name, num_ID_dEvsBRho,init_ID_dEvsBRho[i], hID_dEBRho->GetName());
      
    hID_ZAoQ_dEvsBRhogate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsBRho/Z1AoQ/ID_ZAoQ_dEvsBRhogate%d",i), Form("ZAoQ_dEvsBRhogate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");

    hID_x4AoQ_dEvsBRhogate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsBRho/x4AoQ/ID_x4AoQ_dEvsBRhogate%d",i), Form("x4AoQ_dEvsBRhogate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

    hID_dEvsBRho_dEvsBRhogate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsBRho/dEvsBRho/ID_dEvsBRho_dEvsBRhogate%d",i), Form("dEvsBRho_dEvsBRhogate%d", i), 5000, 3.0 ,5.0, 5000,frs_id->min_z_plot,frs_id->max_z_plot , "BRho", "dE");

    //dEvsZ
    sprintf(name,"cID_dEvsZ%d",i);
    cID_dEvsZ[i] = MakePolyCond("FRS_dEvsZ",name, num_ID_dEvsZ,init_ID_dEvsZ[i], hID_dEdeg_Z->GetName());
      
    hID_ZAoQ_dEvsZgate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ/Z1AoQ/ID_ZAoQ_dEvsZgate%d",i), Form("ZAoQ_dEvsZgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");

    hID_x4AoQ_dEvsZgate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ/x4AoQ/ID_x4AoQ_dEvsZgate%d",i), Form("x4AoQ_dEvsZgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

    hID_dEvsZ_dEvsZgate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ/dEvsZ/ID_dEvsZ_dEvsZgate%d",i), Form("dEvsZ_dEvsZgate%d", i), FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, 1000, 10.,100., "Z from MUSIC41", "dE(S2deg) [a.u.]");

    //Delta Q=0
    hID_x4AoQ_dEvsZ_x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/zero/x4AoQ/ID_x4AoQ_dEvsZ_x4AoQgate%d",i), Form("x4AoQ_dEvsZ_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

    hID_ZAoQ_dEvsZ_x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/zero/Z1AoQ/ID_ZAoQ_dEvsZ_x4AoQgate%d",i), Form("ZAoQ_dEvsZ_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");

    hID_Z1Z2_dEvsZ_x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/zero/Z1Z2/ID_Z1Z2_dEvsZ_x4AoQgate%d",i), Form("Z1Z2_dEvsZ_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, 1000, 10.,100., "Z from MUSIC41", "dE(S2deg) [a.u.]");

    //Delta Q=0 & Z1Z2 gate
    hID_x4AoQ_dEvsZ_Z1Z2_x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/zero/Z1Z2_gated/x4AoQ/ID_x4AoQ_dEvsZ_Z1Z2_x4AoQgate%d",i), Form("x4AoQ_dEvsZ_Z1Z2_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

    hID_ZAoQ_dEvsZ_Z1Z2_x4AoQgate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/zero/Z1Z2_gated/Z1AoQ/ID_ZAoQ_dEvsZ_Z1Z2_x4AoQgate%d",i), Form("ZAoQ_dEvsZ_Z1Z2_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");


    //dEvsZ & x4AoQ Q+1
    sprintf(name,"cID_dEvsZ_x4AoQ_p1%d",i);
    cID_dEvsZ_x4AoQ_p1[i] = MakePolyCond("FRS_dEvsZ_x4AoQ_p1",name, num_ID_dEvsZ_x4AoQ_p1,init_ID_dEvsZ_x4AoQ_p1[i], hID_x4AoQ->GetName());
    hID_x4AoQ_dEvsZ_x4AoQ_p1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/p1/x4AoQ/ID_x4AoQ_dEvsZ_x4AoQ_p1gate%d",i), Form("x4AoQ_dEvsZ_x4AoQ_p1gate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

    hID_ZAoQ_dEvsZ_x4AoQ_p1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/p1/Z1AoQ/ID_ZAoQ_dEvsZ_x4AoQ_p1gate%d",i), Form("ZAoQ_dEvsZ_x4AoQ_p1gate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");

    hID_Z1Z2_dEvsZ_x4AoQ_p1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/p1/Z1Z2/ID_Z1Z2_dEvsZ_x4AoQ_p1gate%d",i), Form("Z1Z2_dEvsZ_x4AoQ_p1gate%d", i), FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, 1000, 10.,100., "Z from MUSIC41", "dE(S2deg) [a.u.]");

    //Delta Q=+1 & Z1Z2 gate
    hID_x4AoQ_dEvsZ_Z1Z2_p1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/p1/Z1Z2_gated/x4AoQ/ID_x4AoQ_dEvsZ_Z1Z2_p1gate%d",i), Form("x4AoQ_dEvsZ_Z1Z2_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

    hID_ZAoQ_dEvsZ_Z1Z2_p1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/p1/Z1Z2_gated/Z1AoQ/ID_ZAoQ_dEvsZ_Z1Z2_p1gate%d",i), Form("ZAoQ_dEvsZ_Z1Z2_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");


    //dEvsZ & x4AoQ Q-1
    sprintf(name,"cID_dEvsZ_x4AoQ_m1%d",i);
    cID_dEvsZ_x4AoQ_m1[i] = MakePolyCond("FRS_dEvsZ_x4AoQ_m1",name, num_ID_dEvsZ_x4AoQ_m1,init_ID_dEvsZ_x4AoQ_m1[i], hID_x4AoQ->GetName());
    hID_x4AoQ_dEvsZ_x4AoQ_m1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/m1/x4AoQ/ID_x4AoQ_dEvsZ_x4AoQ_m1gate%d",i), Form("x4AoQ_dEvsZ_x4AoQ_m1gate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

    hID_ZAoQ_dEvsZ_x4AoQ_m1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/m1/Z1AoQ/ID_ZAoQ_dEvsZ_x4AoQ_m1gate%d",i), Form("ZAoQ_dEvsZ_x4AoQ_m1gate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");

    hID_Z1Z2_dEvsZ_x4AoQ_m1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/m1/Z1Z2/ID_Z1Z2_dEvsZ_x4AoQ_m1gate%d",i), Form("Z1Z2_dEvsZ_x4AoQ_m1gate%d", i), FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, 1000, 10.,100., "Z from MUSIC41", "dE(S2deg) [a.u.]");

    //Delta Q=-1 & Z1Z2 gate
    hID_x4AoQ_dEvsZ_Z1Z2_m1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/m1/Z1Z2_gated/x4AoQ/ID_x4AoQ_dEvsZ_Z1Z2_m1gate%d",i), Form("x4AoQ_dEvsZ_Z1Z2_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "X at S4 [mm]");

    hID_ZAoQ_dEvsZ_Z1Z2_m1gate[i] = MakeTH2('I', Form("FRS/ID_Gated/dEvsZ_x4AoQ/m1/Z1Z2_gated/Z1AoQ/ID_ZAoQ_dEvsZ_Z1Z2_m1gate%d",i), Form("ZAoQ_dEvsZ_Z1Z2_x4AoQgate%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");


///--------------------------------------/**MHTDC GATES**/------------------------------------------///
      
    //Z vs AoQ
    sprintf(name,"cID_Z_AoQ_mhtdc%d",i);
    cID_Z_AoQ_mhtdc[i] = MakePolyCond("FRS_Z1_AoQ_Gates_mhtdc", name, num_ID_Z_AoQ, init_ID_Z_AoQ_mhtdc[i], hID_Z_AoQ_mhtdc->GetName());

    hID_ZAoQ_ZAoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/Z1AoQ/Z1AoQGated_mhtdc/ID_Z1_AoQ_mhtdc_gate%d",i), Form("ID_ZAoQ_ZAoQgate_mhtdc%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", "Z s2-s4");
        
    hID_Z1Z2_ZAoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/Z1Z2/Z1Z2_Z1AoQGated/ID_Z1Z2_ZAoQgate_mhtdc%d",i), Form("ID_Z1Z2_ZAoQgate_mhtdc%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"MUSIC Z1", "gate on Z1 vs AoQ: MUSIC Z2");
        
    hID_x2AoQ_Z1AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x2AoQ/x2AoQ_Z1AoQGated_mhtdc/ID_x2AoQ_Z1AoQ_mhtdcgate%d",i), Form("ID_x2AoQ_Z1AoQ_mhtdcgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q mhtdc s2-s4", " X at S2 [mm]");
         
    hID_x4AoQ_Z1AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x4AoQ/x4AoQ_Z1AoQGated_mhtdc/ID_x4AoQ_Z1AoQ_mhtdcgate%d",i), Form("ID_x4AoQ_Z1AoQ_mhtdcgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q mhtdc s2-s4", " gate on Z X at S4 [mm]");
                
    hID_a2_Z1AoQgate_mhtdc[i] = MakeTH1('F', Form("FRS/MHTDC/ID_Gated_MHTDC/a2/a2_Z1AoQGated/ID_a2_Z1AoQGated_mhtdc%d",i), Form("ID_a2_Z1AoQGate_mhtdc%d",i), 100, -1000, 1000, "Angle S2 (mrad)");
               
    hID_a4_Z1AoQgate_mhtdc[i] = MakeTH1('F', Form("FRS/MHTDC/ID_Gated_MHTDC/a4/a4_Z1AoQGated/ID_a4_Z1AoQGated_mhtdc%d",i), Form("ID_a4_Z1AoQGate_mhtdc%d",i), 100, -1000, 1000, "Angle S4 (mrad)");


    /// Z1 vs Z2 
    sprintf(name,"cID_Z_Z2gate_mhtdc%d",i);
    cID_Z_Z2gate_mhtdc[i] = MakePolyCond("FRS_Z1_Z2_mhtdc_Gates",name,num_ID_Z_Z2,init_ID_Z_Z2_mhtdc[i], hID_Z_Z2_mhtdc ->GetName());
        
    hID_Z1_Z2gate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/Z1Z2/Z1Z2Gated/ID_Z1_Z2_mhtdcgate%d",i), Form("ID_Z1_Z2_mhtdcgate%d",i),2000,frs_id->min_z_plot,frs_id->max_z_plot, 2000,frs_id->min_z_plot,frs_id->max_z_plot,"Z1 s2-s4", "Z2 s2-s4");
        
    hID_x2AoQ_Z1Z2gate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x2AoQ/x2AoQ_Z1Z2Gated/ID_x2AoQ_Z1Z2_mhtdcgate%d",i), Form("ID_x2AoQ_Z1Z2_mhtdcgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on Z X at S2 [mm]");

    hID_x4AoQ_Z1Z2gate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x4AoQ/x4AoQ_Z1Z2Gated/ID_x4AoQ_Z1Z2_mhtdcgate%d",i), Form("ID_x4AoQ_Z1Z2_mhtdcgate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on Z X at S4 [mm]");

    hID_ZAoQ_Z1Z2gate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/Z1AoQ/Z1AoQ_Z1Z2Gated/ID_Z1AoQ_Z1Z2gate_mhtdc%d",i), Form("ID_Z1AoQ_Z1Z2gate%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", " Z music2");
        
    hID_a2_Z1Z2gate_mhtdc[i] = MakeTH1('F', Form("FRS/MHTDC/ID_Gated_MHTDC/a2/a2_Z1Z2Gated/ID_a2_Z1Z2gate_mhtdc%d",i), Form("ID_a2_Z1Z2gate_mhtdc%d",i), 100, -1000, 1000, "Angle S2 (mrad)");
        
    hID_a4_Z1Z2gate_mhtdc[i] = MakeTH1('F', Form("FRS/MHTDC/ID_Gated_MHTDC/a4/a4_Z1Z2Gated/ID_a4_Z1Z2gate_mhtdc%d",i), Form("ID_a4_Z1Z2gate_mhtdc%d",i), 100, -1000, 1000, "Angle S4 (mrad)");

    
    ///X2 AoQ
    sprintf(name,"cID_x2AoQ_mhtdc%d",i);
    cID_x2AoQ_mhtdc[i] = MakePolyCond("FRS_X2_Gates_mhtdc",name,num_ID_x2AoQ,init_ID_x2AoQ_mhtdc[i], hID_x2AoQ_mhtdc->GetName());
      
    hID_x2AoQ_x2AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x2AoQ/x2AoQ_x2AoQGated/ID_x2AoQ_x2AoQgate_mhtdc%d",i), Form("ID_x2AoQ_x2AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on X2 AoQ, ID X2: X at S2 [mm]");
      
    hID_Z1Z2_x2AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/Z1Z2/Z1Z2_x2AoQGated/ID_Z1Z2_x2AoQgate_mhtdc%d",i), Form("ID_Z1Z2_x2AoQgate_mhtdc%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"MUSIC Z1", "gate on x2 vs AoQ: MUSIC Z2");
      
    hID_x2AoQ_Z1Z2x2AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x2AoQ/x2AoQ_Z1Z2x2AoQGated/ID_x2AoQ_Z1Z2x2AoQgate_mhtdc%d",i), Form("ID_x2AoQ_Z1Z2x2AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on Z1Z2 X2AoQ, ID X2: X at S4 [mm]");
      
    hID_x4AoQ_Z1Z2x4AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x4AoQ/x4AoQ_Z1Z2x4AoQGated/ID_x4AoQ_Z1Z2x2AoQgate_mhtdc%d",i), Form("ID_x4AoQ_Z1Z2x4AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on Z1Z2 X2AoQ, ID X2: X at S4 [mm]");

    hID_Z1Z2_Z1Z2x2AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x2AoQ/Z1Z2_Z1Z2x2AoQGated/ID_Z1Z2_Z1Z2x2AoQgate_mhtdc%d",i), Form("ID_Z1Z2_Z1Z2x2AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"Z music41", "gate on Z1Z2 X2AoQ, ID X2: Z music41");
       
    hID_ZAoQ_x2AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/Z1AoQ/Z1AoQ_Z1Z2x2AoQGated/ID_Z1AoQ_Z1Z2x2AoQgate_mhtdc%d",i), Form("ID_Z1AoQ_Z1Z2x2AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", " gate on Z1Z2 X2AoQ Z music41");
      
    hID_a2_Z1Z2x2AoQgate_mhtdc[i] = MakeTH1('F', Form("FRS/MHTDC/ID_Gated_MHTDC/a2/a2_Z1Z2x2AoQGated/ID_a2_Z1Z2x2AoQgate_mhtdc%d",i), Form("ID_a2_Z1Z2x2AoQgate_mhtdc%d",i), 100, -1000, 1000, "Angle S2 (mrad)");
       
    hID_a4_Z1Z2x2AoQgate_mhtdc[i] = MakeTH1('F', Form("FRS/MHTDC/ID_Gated_MHTDC/a4/a4_Z1Z2x2AoQGated/ID_a4_Z1Z2x2AoQgate_mhtdc%d",i), Form("ID_a4_Z1Z2x2AoQgate_mhtdc%d",i), 100, -1000, 1000, "Angle S4 (mrad)");


    ///X4 vs AoQ
    sprintf(name,"cID_x4AoQ_mhtdc%d",i);
    cID_x4AoQ_mhtdc[i] = MakePolyCond("FRS_X4_Gates",name,num_ID_x4AoQ,init_ID_x4AoQ[i], hID_x4AoQ_mhtdc->GetName());
      
    hID_x4AoQ_x4AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x4AoQ/x4AoQ_x4AoQGated/ID_x4AoQ_x4AoQgate_mhtdc%d",i), Form("ID_x2AoQ_x2AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on X4 AoQ, ID X4: X at S4 [mm]");
      
    hID_Z1Z2_x4AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/Z1Z2/Z1Z2_x4AoQGated/ID_Z1Z2_x4AoQgate_mhtdc%d",i), Form("ID_Z1Z2_x4AoQgate_mhtdc%d", i), FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"MUSIC Z1", "gate on x4 vs AoQ: MUSIC Z2");

    hID_x4AoQ_Z1Z2x4AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x4AoQ/x4AoQ_Z1Z2x4AoQGated/ID_x4AoQ_Z1Z2x4AoQgate_mhtdc%d",i), Form("ID_x4AoQ_Z1Z2x4AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on Z1Z2 X4AoQ, ID X4: X at S4 [mm]");
      
    hID_x2AoQ_Z1Z2x4AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x2AoQ/x2AoQ_Z1Z2x4AoQGated/ID_x2AoQ_Z1Z2x4AoQgate_mhtdc%d",i), Form("ID_x2AoQ_Z1Z2x4AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, 200,-100.,100.,"A/Q s2-s4", "gate on Z1Z2 X4AoQ, ID X4: X at S2 [mm]");

    hID_Z1Z2_Z1Z2x4AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/x4AoQ/Z1Z2_Z1Z2x4AoQGated/ID_Z1Z2_Z1Z2x4AoQgate_mhtdc%d",i), Form("ID_Z1Z2_Z1Z2x4AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"Z music41", "gate on Z1Z2 X4AoQ, ID X4: Z music41");
       
    hID_ZAoQ_x4AoQgate_mhtdc[i] = MakeTH2('I', Form("FRS/MHTDC/ID_Gated_MHTDC/Z1AoQ/Z1AoQ_Z1Z2x4AoQGated/ID_Z1AoQ_Z1Z2x4AoQgate_mhtdc%d",i), Form("ID_Z1AoQ_Z1Z2x4AoQgate_mhtdc%d",i),FRS_HISTO_BIN,frs_id->min_aoq_plot,frs_id->max_aoq_plot, FRS_HISTO_BIN,frs_id->min_z_plot,frs_id->max_z_plot,"A/Q s2-s4", " gate on Z1Z2 X4AoQ Z music41");
            
    hID_a2_Z1Z2x4AoQgate_mhtdc[i] = MakeTH1('F', Form("FRS/MHTDC/ID_Gated_MHTDC/a2/a2_Z1Z2x4AoQGated/ID_a2_Z1Z2x4AoQgate_mhtdc%d",i), Form("ID_a2_Z1Z2x4AoQgate_mhtdc%d",i), 100, -1000, 1000, "Angle S2 (mrad)");
      
    hID_a4_Z1Z2x4AoQgate_mhtdc[i] = MakeTH1('F', Form("FRS/MHTDC/ID_Gated_MHTDC/a4/a4_Z1Z2x4AoQGated/ID_a4_Z1Z2x4AoQgate_mhtdc%d",i), Form("ID_a4_Z1Z2x4AoQgate_mhtdc%d",i), 100, -1000, 1000, "Angle S4 (mrad)");

  }
}

/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------FRS HISTOS ---------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Process_FRS_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput){


    
  ///Defined in the setup file.
  if(FRS_CORR==true){
    FRS_AoQ = FRS_AoQ_corr;
    for (int i=0; i<10;i++) FRS_AoQ_mhtdc[i] = FRS_AoQ_corr_mhtdc[i];
  }

  pOutput->pFRS_AoQ = FRS_AoQ;
  pOutput->pFRS_z = FRS_z;
  pOutput->pFRS_z2 = FRS_z2;
  pOutput->pFRS_beta = FRS_beta;
  pOutput->pFRS_dEdeg = FRS_dEdeg;
  pOutput->pFRS_dEdegoQ = FRS_dEdegoQ;

  for (int i=0; i<10;i++){ 
    pOutput->pFRS_AoQ_mhtdc[i] = FRS_AoQ_mhtdc[i];
    pOutput->pFRS_z_mhtdc[i] = FRS_z_mhtdc[i];
    pOutput->pFRS_z2_mhtdc[i] = FRS_z2_mhtdc[i]; 
    pOutput->pFRS_beta_mhtdc[i] = FRS_beta_mhtdc[i];
    pOutput->pFRS_dEdeg_mhtdc[i] = FRS_dEdeg_mhtdc[i];
    pOutput->pFRS_dEdegoQ_mhtdc[i] = FRS_dEdegoQ_mhtdc[i];
  }

  pOutput->pFRS_ID_x2 = FRS_ID_x2;
  pOutput->pFRS_ID_y2 = FRS_ID_y2;
  pOutput->pFRS_ID_a2 = FRS_ID_a2;
  pOutput->pFRS_ID_b2 = FRS_ID_b2;
  pOutput->pFRS_ID_x4 = FRS_ID_x4;
  pOutput->pFRS_ID_y4 = FRS_ID_y4;
  pOutput->pFRS_ID_a4 = FRS_ID_a4;
  pOutput->pFRS_ID_b4 = FRS_ID_b4;


  for(int l=0; l<12;l++){
    pOutput->pFRS_sci_l[l] = FRS_sci_l[l];
    pOutput->pFRS_sci_r[l] = FRS_sci_r[l];
    pOutput->pFRS_sci_e[l] = FRS_sci_e[l];
    pOutput->pFRS_sci_tx[l] = FRS_sci_tx[l];
    pOutput->pFRS_sci_x[l] = FRS_sci_x[l];
  }

  if(pOutput->pFRS_WR > 0) pOutput->pt_lastSC41 = pOutput->pFRS_WR;

///--------------------------------------/**FRS TAC Histograms and PID Gates**/------------------------------------------///

  /****  S4  (MUSIC 1)   */
  if(FRS_z!=0)hID_Z1_corr->Fill(FRS_z);

  /****  S4  (MUSIC 2)   */
  if(FRS_z2!=0) hID_Z2_corr->Fill(FRS_z2);

  ///Z and A/Q Time in mins
  if(FRS_z>0 && FRS_time_mins>0)hID_Z1_vs_T->Fill(FRS_time_mins,FRS_z);
  if(FRS_z2>0 && FRS_time_mins>0)hID_Z2_vs_T->Fill(FRS_time_mins,FRS_z2);
     
  if(FRS_z>0 && FRS_time_mins>0)hID_AoQ_vs_T->Fill(FRS_time_mins,FRS_AoQ);
   
  ///AoQ vs Z
   
  if(pInput->fFRS_AoQ>0 && FRS_z>0) hID_Z_AoQ->Fill(pInput->fFRS_AoQ, FRS_z);
  if(pInput->fFRS_AoQ_corr>0 && FRS_z>0)   hID_Z_AoQ_corr->Fill(pInput->fFRS_AoQ_corr, FRS_z);  //S2-S4 correlated
 
  ///Z1 Z2
  if(FRS_z>0 && FRS_z2>0) hID_Z_Z2->Fill(FRS_z,FRS_z2);
    
  if(TMath::Abs(FRS_z-FRS_z2)<0.4){
    hID_Z_AoQ_zsame->Fill(FRS_AoQ, FRS_z);
    hID_x4AoQ_zsame->Fill(FRS_AoQ, FRS_ID_x4);
    hID_x2AoQ_zsame->Fill(FRS_AoQ, FRS_ID_x2);
  }

  //sci    
  if(FRS_z>0 && FRS_sci_e[5]) hID_Z_E_SC41->Fill(FRS_sci_e[5],FRS_z);
  if(FRS_z>0 && FRS_sci_e[6]) hID_Z_E_SC42->Fill(FRS_sci_e[6],FRS_z);
             
  ///AoQ vs X2
  if(FRS_AoQ>0 && FRS_ID_x2>-100 && FRS_ID_x2<100)  hID_x2AoQ->Fill(pInput->fFRS_AoQ, FRS_ID_x2);
  if(FRS_AoQ>0 && FRS_ID_x4>-100 && FRS_ID_x4<100)  hID_x4AoQ->Fill(pInput->fFRS_AoQ, FRS_ID_x4);
  if(pInput->fFRS_AoQ_corr>0 && FRS_ID_x2>-100 && FRS_ID_x2<100)  hID_x2AoQ_corr->Fill(pInput->fFRS_AoQ_corr, FRS_ID_x2);
  if(pInput->fFRS_AoQ_corr>0 && FRS_ID_x4>-100 && FRS_ID_x4<100)  hID_x4AoQ_corr->Fill(pInput->fFRS_AoQ_corr, FRS_ID_x4);
 
  ///Charge states
  if(FRS_z>0 && FRS_dEdegoQ!=0)   hID_dEdegoQ_Z->Fill(FRS_z, FRS_dEdegoQ);
  if(FRS_z>0 && FRS_dEdeg!=0)  hID_dEdeg_Z->Fill(FRS_z, FRS_dEdeg);
   
  //Angles vs AoQ 
  if(FRS_AoQ!=0 && FRS_ID_a2!=0)hID_a2AoQ->Fill(FRS_AoQ,FRS_ID_a2);
  if(FRS_AoQ!=0 && FRS_ID_a4!=0)hID_a4AoQ->Fill(FRS_AoQ,FRS_ID_a4);
   
  //Z vs Energy loss Music 2    
  if(FRS_z!=0 && pInput->fFRS_Music_dE[1]!=0)hID_Z_dE2->Fill(FRS_z,pInput->fFRS_Music_dE[1]);
 
  //X2 vs X4
  if(FRS_ID_x2!=0&&FRS_ID_x4!=0 ) hID_x2x4->Fill(FRS_ID_x2, FRS_ID_x4);
  if(FRS_AoQ!=0 && FRS_sci_e[5]!=0) hID_SC41dE_AoQ->Fill(FRS_AoQ, FRS_sci_e[5]);

  if(FRS_sci_tof2!=0 && pInput->fFRS_Music_dE[0]!=0) hID_dEToF->Fill(FRS_sci_tof2, pInput->fFRS_Music_dE[0]);

  if(FRS_z!=0 && FRS_ID_x2!=0) hID_x2z1->Fill(FRS_z, FRS_ID_x2);
  if(FRS_z!=0 && FRS_ID_x4!=0) hID_x4z1->Fill(FRS_z, FRS_ID_x4);
  if(FRS_z2!=0 && FRS_ID_x4!=0) hID_x4z2->Fill(FRS_z2, FRS_ID_x4);//GEEBART
     
  if (FRS_ID_brho[0]!= 0 && FRS_ID_brho[1]!=0 && FRS_z!=0) hID_dEBRho->Fill(FRS_ID_brho[0]-FRS_ID_brho[1], z1_corr);//GEEBART
  if (FRS_ID_brho[0]!= 0 && FRS_ID_brho[1]!=0 && pInput->fFRS_Music_dE[0]!= 0) hID_dBrho_dE4->Fill(FRS_ID_brho[0]-FRS_ID_brho[1], FRS_dE[0]);
  if (FRS_ID_brho[0]!= 0 && FRS_ID_brho[1]!=0) hID_BRho1v2->Fill(FRS_ID_brho[0], FRS_ID_brho[1]);

  if(FRS_ID_x2!=0 && pInput->fFRS_Music_dE[0]!=0)hID_E_x2->Fill(FRS_ID_x2,pInput->fFRS_Music_dE[0]);
  if(FRS_ID_x4!=0 && pInput->fFRS_Music_dE[0]!=0)hID_E_x4->Fill(FRS_ID_x4,pInput->fFRS_Music_dE[0]);
  
  if(FRS_ID_x2!=0 && FRS_ID_a2!=0)hID_x2a2->Fill(FRS_ID_x2,FRS_ID_a2);
  if(FRS_ID_y2!=0 && FRS_ID_b2!=0)hID_y2b2->Fill(FRS_ID_y2,FRS_ID_b2);
  if(FRS_ID_x4!=0 && FRS_ID_a4!=0)hID_x4a4->Fill(FRS_ID_x4,FRS_ID_a4);
  if(FRS_ID_x4!=0 && FRS_ID_b4!=0) hID_y4b4->Fill(FRS_ID_y4,FRS_ID_b4);
     
  if(FRS_z!=0 && FRS_sci_l[2]!=0 && FRS_sci_r[2]!=0)hID_Z_Sc21E->Fill(FRS_z, sqrt(FRS_sci_l[2]*FRS_sci_r[2]));
    
  ///Define PID Gates 
  for(int g=0; g<MAX_FRS_GATE; g++){
          
    ///GATE: Z1 vs A/Q
    if(cID_Z_AoQ[g]->Test(FRS_AoQ, FRS_z)==true){
      pOutput->pFRS_ZAoQ_pass[g] =true;

      hID_ZAoQ_ZAoQgate[g]->Fill(FRS_AoQ, FRS_z);
      hID_Z1Z2_ZAoQgate[g]->Fill(FRS_z, FRS_z2);
      hID_x2AoQ_Z1AoQgate[g]->Fill(FRS_AoQ, FRS_ID_x2);
      hID_x4AoQ_Z1AoQgate[g]->Fill(FRS_AoQ, FRS_ID_x4);
      hID_a2_Z1AoQgate[g]->Fill(FRS_ID_a2);
      hID_a4_Z1AoQgate[g]->Fill(FRS_ID_a4);
    }
           
    ///GATE: Z1 vs Z2
    if(cID_Z_Z2gate[g]->Test(FRS_z, FRS_z2)==true){
      pOutput->pFRS_Z_Z2_pass[g] = true;

      hID_Z1_Z2gate[g]->Fill(FRS_z,FRS_z2);
      hID_a2_Z1Z2gate[g]->Fill(FRS_ID_a2);
      hID_a4_Z1Z2gate[g]->Fill(FRS_ID_a4);
      
      if(FRS_AoQ>1 && FRS_AoQ<3 &&   FRS_ID_x2 > -100 && FRS_ID_x2<100){
        hID_x2AoQ_Z1Z2gate[g]->Fill(FRS_AoQ, FRS_ID_x2);
      }

      if(FRS_AoQ>1 && FRS_AoQ<3   && FRS_ID_x4 > -100 && FRS_ID_x4<100){
        hID_x4AoQ_Z1Z2gate[g]->Fill(FRS_AoQ, FRS_ID_x4);
      }
         
      if(FRS_AoQ>1 && FRS_AoQ<3){
        hID_ZAoQ_Z1Z2gate[g] ->Fill(FRS_AoQ, FRS_z);
      }
    }

       
    ///GATE: ID vs x2AoQ
    if(cID_x2AoQ[g]->Test(FRS_AoQ, FRS_ID_x2)==true){
      pOutput->pFRS_x2AoQ_pass[g] = true;

      hID_x2AoQ_x2AoQgate[g]->Fill(FRS_AoQ, FRS_ID_x2);
      hID_Z1Z2_x2AoQgate[g]->Fill(FRS_z, FRS_z2);
      ///Z1 Z2 + x2AoQ
      ///The selected Z1 Z2 gate for this part can be found in the Correlations_config.dat file 

      if(cID_Z_Z2gate[g]->Test(FRS_z, FRS_z2)==true){
        if(cID_x2AoQ[g]->Test(FRS_AoQ, FRS_ID_x2)==true){
          hID_x2AoQ_Z1Z2x2AoQgate[g]->Fill(FRS_AoQ, FRS_ID_x2);
          hID_x4AoQ_Z1Z2x2AoQgate[g]->Fill(FRS_AoQ, FRS_ID_x4);
          hID_ZAoQ_Z1Z2x2AoQgate[g]->Fill(FRS_AoQ, FRS_z);
          hID_a2_Z1Z2x2AoQgate[g]->Fill(FRS_ID_a2);
          hID_a4_Z1Z2x2AoQgate[g]->Fill(FRS_ID_a4); 
        }   
      }
    }

    ///GATE: ID vs x4AoQ
    if(cID_x4AoQ[g]->Test(FRS_AoQ, FRS_ID_x4)==true){ 
      pOutput->pFRS_x4AoQ_pass[g] = true;

      hID_x4AoQ_x4AoQgate[g]->Fill(FRS_AoQ, FRS_ID_x4);
      hID_Z1Z2_x4AoQgate[g]->Fill(FRS_z, FRS_z2);
      ///Z1 Z2 + x4AoQ
      ///The selected Z1 Z2 gate for this part can be found in the Correlations_config.dat file 
      if(cID_Z_Z2gate[g]->Test(FRS_z, FRS_z2)==true){        
        hID_x2AoQ_Z1Z2x4AoQgate[g]->Fill(FRS_AoQ, FRS_ID_x2);
        hID_x4AoQ_Z1Z2x4AoQgate[g]->Fill(FRS_AoQ, FRS_ID_x4);
        hID_ZAoQ_Z1Z2x4AoQgate[g]->Fill(FRS_AoQ, FRS_z);
        hID_a2_Z1Z2x4AoQgate[g]->Fill(FRS_ID_a2);
        hID_a4_Z1Z2x4AoQgate[g]->Fill(FRS_ID_a4);
      }
    }
        


    //dEvsZ Gates
    if(cID_dEvsZ[g]->Test(FRS_z, FRS_dEdeg)==true){
      pOutput->pFRS_dEvsZ_pass[g] =true;
      hID_ZAoQ_dEvsZgate[g]->Fill(FRS_AoQ_corr, FRS_z);
      hID_x4AoQ_dEvsZgate[g]->Fill(FRS_AoQ_corr, FRS_ID_x4);
      hID_dEvsZ_dEvsZgate[g]->Fill(FRS_z, FRS_dEdeg);
    
    }

    //devsZ deltaQ=0
    if (cID_dEvsZ[1]->Test(FRS_z, FRS_dEdeg)==true){
      if(cID_x4AoQ[g]->Test(FRS_AoQ_corr,FRS_ID_x4)==true){
        pOutput->pFRS_dEvsZ_x4AoQ_pass[g] =true;
        hID_x4AoQ_dEvsZ_x4AoQgate[g]->Fill(FRS_AoQ_corr, FRS_ID_x4);
        hID_ZAoQ_dEvsZ_x4AoQgate[g]->Fill(FRS_AoQ_corr, FRS_z);
        hID_Z1Z2_dEvsZ_x4AoQgate[g]->Fill(FRS_z, FRS_z2);

        if(cID_Z_Z2gate[g]->Test(FRS_z, FRS_z2)==true){
          hID_x4AoQ_dEvsZ_Z1Z2_x4AoQgate[g]->Fill(FRS_AoQ_corr, FRS_ID_x4);
          hID_ZAoQ_dEvsZ_Z1Z2_x4AoQgate[g]->Fill(FRS_AoQ_corr, FRS_z);
        }
      }
    }  

    //dEvsZ Q+1 Gates
    if(cID_dEvsZ[3]->Test(FRS_z, FRS_dEdeg)==true){
      if(cID_dEvsZ_x4AoQ_p1[g]->Test(FRS_AoQ_corr,FRS_ID_x4)==true){
        pOutput->pFRS_dEvsZ_x4AoQ_p1_pass[g] =true;
        hID_x4AoQ_dEvsZ_x4AoQ_p1gate[g]->Fill(FRS_AoQ_corr, FRS_ID_x4);
        hID_ZAoQ_dEvsZ_x4AoQ_p1gate[g]->Fill(FRS_AoQ_corr, FRS_z);
        hID_Z1Z2_dEvsZ_x4AoQ_p1gate[g]->Fill(FRS_z, FRS_z2);

        if(cID_Z_Z2gate[g]->Test(FRS_z, FRS_z2)==true){
          hID_x4AoQ_dEvsZ_Z1Z2_p1gate[g]->Fill(FRS_AoQ_corr, FRS_ID_x4);
          hID_ZAoQ_dEvsZ_Z1Z2_p1gate[g]->Fill(FRS_AoQ_corr, FRS_z);
        }
      }
    }

    //dEvsZ Q-1 Gates
    if(cID_dEvsZ[2]->Test(FRS_z, FRS_dEdeg)==true){
      if(cID_dEvsZ_x4AoQ_m1[g]->Test(FRS_AoQ_corr,FRS_ID_x4)==true){
        pOutput->pFRS_dEvsZ_x4AoQ_m1_pass[g] =true;
        hID_x4AoQ_dEvsZ_x4AoQ_m1gate[g]->Fill(FRS_AoQ_corr, FRS_ID_x4);
        hID_ZAoQ_dEvsZ_x4AoQ_m1gate[g]->Fill(FRS_AoQ_corr, FRS_z);
        hID_Z1Z2_dEvsZ_x4AoQ_m1gate[g]->Fill(FRS_z, FRS_z2);

        if(cID_Z_Z2gate[g]->Test(FRS_z, FRS_z2)==true){
          hID_x4AoQ_dEvsZ_Z1Z2_m1gate[g]->Fill(FRS_AoQ_corr, FRS_ID_x4);
          hID_ZAoQ_dEvsZ_Z1Z2_m1gate[g]->Fill(FRS_AoQ_corr, FRS_z);
        }
      }
    }


    //dEvsBRho Gates
    if(cID_dEvsBRho[g]->Test(FRS_ID_brho[0]-FRS_ID_brho[1], FRS_dE[0])==true){
      pOutput->pFRS_dEvsBRho_pass[g] =true;

      hID_ZAoQ_dEvsBRhogate[g]->Fill(FRS_AoQ_corr, FRS_z);
      hID_x4AoQ_dEvsBRhogate[g]->Fill(FRS_AoQ_corr, FRS_ID_x4);
      hID_dEvsBRho_dEvsBRhogate[g]->Fill(FRS_ID_brho[0]-FRS_ID_brho[1], FRS_dE[0]);
    }
  }

    
        
///--------------------------------------/**FRS MHTDC Histograms and PID Gates**/------------------------------------------///
/// Loop on MHTDC arrays
  for (int i=0; i<10 ; i++){

    ///Z1 vs Time 
    if(FRS_z_mhtdc[i]>0 && FRS_time_mins>0) hID_Z_mhtdc_T->Fill(FRS_time_mins, FRS_z_mhtdc[i]);
    
    ///AoQ vs Time  
    if(FRS_AoQ_mhtdc[i]>0 && FRS_time_mins>0) hID_AoQ_mhtdc_T->Fill(FRS_time_mins, FRS_AoQ_mhtdc[i]);
    
    ///AoQ vs Z 
    if(pInput->fFRS_AoQ_mhtdc[i]>0 && FRS_z_mhtdc[i]>0) hID_Z_AoQ_mhtdc->Fill(pInput->fFRS_AoQ_mhtdc[i], FRS_z_mhtdc[i]);
    if(pInput->fFRS_AoQ_corr_mhtdc[i]>0 && FRS_z_mhtdc[i]>0) hID_Z_AoQ_corr_mhtdc->Fill(pInput->fFRS_AoQ_corr_mhtdc[i], FRS_z_mhtdc[i]);
    
    //Z1 Z2 
    if(FRS_z_mhtdc[i]>0 && FRS_z2_mhtdc[i]>0) hID_Z_Z2_mhtdc->Fill(FRS_z_mhtdc[i],FRS_z2_mhtdc[i]);
    
    if(TMath::Abs(FRS_z_mhtdc[i]-FRS_z2_mhtdc[i])<0.4){
      hID_Z_AoQ_zsame_mhtdc->Fill(FRS_AoQ_mhtdc[i], FRS_z_mhtdc[i]);
      hID_x4AoQ_zsame_mhtdc->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x4);
      hID_x2AoQ_zsame_mhtdc->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x2);
    }
    
    ///AoQ vs X2
    if(FRS_AoQ_mhtdc[i]>0 && FRS_ID_x2>-100 && FRS_ID_x2<100)  hID_x2AoQ_mhtdc->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x2);
    
    ///AoQ vs X4
    if(FRS_AoQ_mhtdc[i]>0 && FRS_ID_x4>-100 && FRS_ID_x4<100)  hID_x4AoQ_mhtdc->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x4);

    ///Charge states
    if(FRS_z_mhtdc[i]>0 && FRS_dEdegoQ!=0)   hID_dEdegoQ_Z_mhtdc->Fill(FRS_z_mhtdc[i], FRS_dEdegoQ_mhtdc[i]);
    if(FRS_z_mhtdc[i]>0 && FRS_dEdeg!=0)  hID_dEdeg_Z_mhtdc->Fill(FRS_z_mhtdc[i], FRS_dEdegoQ_mhtdc[i]);
    
    //Angles vs AoQ 
    if(FRS_AoQ_mhtdc[i]!=0 && FRS_ID_a2!=0) hID_a2AoQ_mhtdc->Fill(FRS_AoQ_mhtdc[i],FRS_ID_a2);
    if(FRS_AoQ_mhtdc[i]!=0 && FRS_ID_a4!=0) hID_a4AoQ_mhtdc->Fill(FRS_AoQ_mhtdc[i],FRS_ID_a4);
    
    if(FRS_z_mhtdc[i]!=0 && pInput->fFRS_Music_dE[1]!=0) hID_Z_dE2_mhtdc->Fill(FRS_z_mhtdc[i],pInput->fFRS_Music_dE[1]);
    
    if(FRS_z_mhtdc[i]!=0 && FRS_sci_l[2]!=0 && FRS_sci_r[2]!=0)hID_Z_Sc21E_mhtdc->Fill(FRS_z_mhtdc[i], sqrt(FRS_sci_l[2]*FRS_sci_r[2]));
     
    hID_x2z_mhtdc->Fill(FRS_z_mhtdc[i],FRS_ID_x2);
     
    hID_x4z_mhtdc->Fill(FRS_z_mhtdc[i],FRS_ID_x4);

    //if(FRS_AoQ_mhtdc[i]>0 && FRS_z_mhtdc[i]>0 && i==0) hID_Z_AoQ_mhtdc_first_hit->Fill(FRS_AoQ_mhtdc[i], FRS_z_mhtdc[i]);

    //if(FRS_AoQ_mhtdc[i]>0 && FRS_z_mhtdc[i]>0 && i>0) hID_Z_AoQ_mhtdc_excluding_first_hit->Fill(FRS_AoQ_mhtdc[i], FRS_z_mhtdc[i]);

    //if(FRS_AoQ_mhtdc[i]>0 && FRS_z_mhtdc[i]>0 && i==0) hID_Z_AoQ_corr_mhtdc_first_hit->Fill(FRS_AoQ_mhtdc[i], FRS_z_mhtdc[i]);

    //if(FRS_AoQ_mhtdc[i]>0 && FRS_z_mhtdc[i]>0 && i>0) hID_Z_AoQ_corr_mhtdc_excluding_first_hit->Fill(FRS_AoQ_mhtdc[i], FRS_z_mhtdc[i]);
  }

      
  ///MHTDC PID gates
  for(int g=0; g<MAX_FRS_GATE; g++){

    for (int i=0; i<10; i++){

      ///GATE: AoQ vs Z
      if(cID_Z_AoQ_mhtdc[g]->Test(FRS_AoQ_mhtdc[i], FRS_z_mhtdc[i])==true){
        pOutput->pFRS_ZAoQ_pass_mhtdc[g] =true;

        hID_x2AoQ_Z1AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x2);
        hID_x4AoQ_Z1AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x4);
        hID_ZAoQ_ZAoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_z_mhtdc[i]);
        hID_Z1Z2_ZAoQgate_mhtdc[g]->Fill(FRS_z_mhtdc[i], FRS_z2_mhtdc[i]);
        hID_a2_Z1AoQgate_mhtdc[g]->Fill(FRS_ID_a2);
        hID_a4_Z1AoQgate_mhtdc[g]->Fill(FRS_ID_a4);
      }
        
      ///GATE: Z1 vs Z2
      if(cID_Z_Z2gate_mhtdc[g]->Test(FRS_z_mhtdc[i], FRS_z2_mhtdc[i])==true){
        pOutput->pFRS_Z_Z2_pass_mhtdc[g] = true;
             
        hID_Z1_Z2gate_mhtdc[g]->Fill(FRS_z_mhtdc[i],FRS_z2_mhtdc[i]);

        if(FRS_ID_a2!=0) hID_a2_Z1Z2gate_mhtdc[g] ->Fill(FRS_ID_a2);
        if(FRS_ID_a4!=0) hID_a4_Z1Z2gate_mhtdc[g] ->Fill(FRS_ID_a4);
               
        //GATE: X2 AoQ on Z1 Z2
        if(FRS_ID_x2 > -100 && FRS_ID_x2<100){
          hID_x2AoQ_Z1Z2gate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x2);
        }
        if(FRS_ID_x4 > -100 && FRS_ID_x4<100){
          hID_x4AoQ_Z1Z2gate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x4);
          //Z1 AoQ gated on Z1 Z2
          if(FRS_AoQ_mhtdc[i]!=0 &&FRS_z_mhtdc[i]!=0){
            hID_ZAoQ_Z1Z2gate_mhtdc[g] ->Fill(FRS_AoQ_mhtdc[i], FRS_z_mhtdc[i]);
          }
        }
      }    
                
      /// GATE: x2 vs AoQ
      if(cID_x2AoQ_mhtdc[g]->Test(FRS_AoQ_mhtdc[i], FRS_ID_x2)==true){
        pOutput->pFRS_x2AoQ_pass_mhtdc[g] = true;
                     
        hID_x2AoQ_x2AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x2);
        hID_Z1Z2_x2AoQgate[g]->Fill(FRS_z_mhtdc[i], FRS_z2_mhtdc[i]);
        ///The selected Z1 Z2 gate for this part can be found in the Correlations_config.dat file
        /// Z1 Z2 +X2 AoQ
        if(cID_Z_Z2gate_mhtdc[fCorrel->GZ1Z2_Gate]->Test(FRS_z_mhtdc[i], FRS_z2_mhtdc[i])==true){
          hID_x2AoQ_Z1Z2x2AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x2);
          hID_x4AoQ_Z1Z2x2AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x4);
          hID_Z1Z2_Z1Z2x2AoQgate_mhtdc[g]->Fill(FRS_z_mhtdc[i], FRS_z2_mhtdc[i]);
          hID_a2_Z1Z2x2AoQgate_mhtdc[g]->Fill(FRS_ID_a2);
          hID_a4_Z1Z2x2AoQgate_mhtdc[g]->Fill(FRS_ID_a4);
          hID_Z1Z2_x4AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[g], FRS_z_mhtdc[i]);
        }
      }

      /// GATE: x4 vs AoQ
      if(cID_x4AoQ_mhtdc[g]->Test(FRS_AoQ_mhtdc[i], FRS_ID_x4)==true){
        pOutput->pFRS_x4AoQ_pass_mhtdc[g] = true;
                     
        hID_x4AoQ_x4AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x4);
        hID_Z1Z2_x4AoQgate[g]->Fill(FRS_z_mhtdc[i], FRS_z2_mhtdc[i]);
        ///The selected Z1 Z2 gate for this part can be found in the Correlations_config.dat file
        /// Z1 Z2 +X4 AoQ
        if(cID_Z_Z2gate_mhtdc[fCorrel->GZ1Z2_Gate]->Test(FRS_z_mhtdc[i], FRS_z2_mhtdc[i])==true){
          hID_x2AoQ_Z1Z2x4AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x2);
          hID_x4AoQ_Z1Z2x4AoQgate_mhtdc[g]->Fill(FRS_AoQ_mhtdc[i], FRS_ID_x4);
          hID_Z1Z2_Z1Z2x4AoQgate_mhtdc[g]->Fill(FRS_z_mhtdc[i], FRS_z2_mhtdc[i]);
          hID_a2_Z1Z2x4AoQgate_mhtdc[g]->Fill(FRS_ID_a2);
          hID_a4_Z1Z2x4AoQgate_mhtdc[g]->Fill(FRS_ID_a4);   
        }
      }
             
             
    }
  }
}


/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------    AIDA   ---------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/

void EventAnlProc::Make_Aida_Histos(){
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
  implants_strip_xy.resize(conf->DSSDs());
  implants_pos_xy.resize(conf->DSSDs());
  implants_strip_xy_stopped.resize(conf->DSSDs());
  implants_pos_xy_stopped.resize(conf->DSSDs());
  implants_e.resize(conf->DSSDs());
  implants_e_xy.resize(conf->DSSDs());
  implants_time_delta.resize(conf->DSSDs());
  implants_strip_1d.resize(conf->DSSDs());
  implants_per_event.resize(conf->DSSDs());
  decays_strip_xy.resize(conf->DSSDs());
  decays_pos_xy.resize(conf->DSSDs());
  decays_e.resize(conf->DSSDs());
  decays_e_xy.resize(conf->DSSDs());
  decays_time_delta.resize(conf->DSSDs());
  decays_strip_1d.resize(conf->DSSDs());
  decays_per_event.resize(conf->DSSDs());
  //implants_channels.resize(conf->DSSDs());
  //decays_channels.resize(conf->DSSDs());
  implants_x_ex.resize(conf->DSSDs());
  implants_y_ey.resize(conf->DSSDs());

  #ifdef AIDA_PULSER_ALIGN
    aida_pulser_time = MakeTH2('I', "AIDA/Pulser_Time", "AIDA Pulser Time Comparison", 768, 0, 768, 2000, -4000, 4000);
  #endif

  int xstrips = 128;
  if (conf->Wide()) xstrips = 386;
  double xmax = 37.8;
  if (conf->Wide()) xmax = 113.4;

  for (int i = 0; i < conf->DSSDs(); ++i){
    implants_strip_xy[i] = MakeTH2('I', Form("AIDA/Implants/DSSD%d_implants_strip_XY", i+1), Form("DSSD %d implant hit pattern", i+1), xstrips, 0, xstrips, 128, 0, 128, "X strip", "Y strip");
    implants_strip_xy_stopped[i]= MakeTH2('I', Form("AIDA/Implants_Stopped/DSSD%d_implants_stopped_strip_XY", i+1), Form("DSSD %d implant stopped hit pattern", i+1), xstrips, 0, xstrips, 128, 0, 128, "X strip", "Y strip");
    implants_pos_xy[i] = MakeTH2('D', Form("AIDA/Implants/DSSD%d_implants_pos_XY", i+1), Form("DSSD %d implant position", i+1), xstrips, -xmax, xmax, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    implants_pos_xy_stopped[i]= MakeTH2('I', Form("AIDA/Implants_Stopped/DSSD%d_implants_stopped_pos_XY", i+1), Form("DSSD %d implant stopped position hit pattern", i+1), xstrips, -xmax, xmax, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    implants_e[i] = MakeTH1('F', Form("AIDA/Implants/DSSD%d_implants_energy", i+1), Form("DSSD %d implant energy", i+1), 2000, 0, 20000, "Implant Energy/MeV");
    //implants_e_xy[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_energy_XY", i+1), Form("DSSD %d implant front energy vs back energy", i+1), 1000, 0, 10000, 1000, 0, 10000, "X Energy", "Y Energy");
    implants_time_delta[i] = MakeTH1('F', Form("AIDA/Implants/DSSD%d_implants_time_delta", i+1), Form("DSSD %d implant front vs back time", i+1), 1000, -10000, 10000, "Time Difference/ns");
    implants_strip_1d[i] = MakeTH1('I', Form("AIDA/Implants/DSSD%d_implants_strip_1d", i+1), Form("DSSD %d implant 1D hit pattern", i+1), 128 + xstrips, 0, 128 + xstrips, "Strip number");
    implants_per_event[i] = MakeTH1('I', Form("AIDA/Implants/DSSD%d_implants_per_event", i+1), Form("DSSD %d implants per event", i+1), 100, 0, 100, "Number of implants");
    implants_x_ex[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_x_ex", i+1), Form("DSSD %d Ex vs X position", i+1), 128, 0, 128, 2000, 0, 20000, "X Strip", "X Energy");
    implants_y_ey[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_y_ey", i+1), Form("DSSD %d Ey vs Y position", i+1), 128, 0, 128, 2000, 0, 20000, "Y Strip", "Y Energy");

    decays_strip_xy[i] = MakeTH2('I', Form("AIDA/Decays/DSSD%d_decays_strip_XY", i+1), Form("DSSD %d decay hit pattern", i+1), xstrips, 0, xstrips, 128, 0, 128, "X strip", "Y strip");
    decays_pos_xy[i] = MakeTH2('D', Form("AIDA/Decays/DSSD%d_decays_pos_XY", i+1), Form("DSSD %d decay position", i+1), xstrips, -xmax, xmax, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    decays_e[i] = MakeTH1('F', Form("AIDA/Decays/DSSD%d_decays_energy", i+1), Form("DSSD %d decay energy", i+1), 1000, 0, 20000, "Decay Energy/keV");
    decays_e_xy[i] = MakeTH2('F', Form("AIDA/Decays/DSSD%d_decays_energy_XY", i+1), Form("DSSD %d decay front energy vs back energy", i+1), 1000, 0, 10000, 1000, 0, 20000, "X Energy", "Y Energy");
    decays_time_delta[i] = MakeTH1('F', Form("AIDA/Decays/DSSD%d_decays_time_delta", i+1), Form("DSSD %d decay front vs back time", i+1), 1000, -10000, 10000, "Time Difference/ns");
    decays_strip_1d[i] = MakeTH1('I', Form("AIDA/Decays/DSSD%d_decays_strip_1d", i+1), Form("DSSD %d decay 1D hit pattern", i+1), 128 + xstrips, 0, 128 + xstrips, "Strip number");
    decays_per_event[i] = MakeTH1('I', Form("AIDA/Decays/DSSD%d_decays_per_event", i+1), Form("DSSD %d decays per event", i+1), 100, 0, 100, "Number of decays");

    //implants_channels[i] = MakeTH1('I', Form("AIDA/DSSD%d_implants_channels", i+1), Form("DSSD %d number of implant channels", i+1), 769, 0, 769);
    //decays_channels[i] = MakeTH1('I', Form("AIDA/DSSD%d_decays_channels", i+1), Form("DSSD %d number of decay channels", i+1), 769, 0, 769);
  }
}

///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///

using aida_coord_t = std::tuple<int, int, int>;
// https://stackoverflow.com/a/7222201/916549
template<typename T>

inline void hash_combine(std::size_t& seed, const T& val){
  std::hash<T> hasher;
  seed ^= hasher(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct aida_coord_hash{
  inline size_t operator()(const aida_coord_t& val) const{
    size_t seed = 0;
    hash_combine(seed, std::get<0>(val));
    hash_combine(seed, std::get<1>(val));
    hash_combine(seed, std::get<2>(val));
    return seed;
  }
};

void EventAnlProc::ProcessAida(EventUnpackStore* pInputMain, EventAnlStore* pOutput){
  // int Aida_hits =0;
  //double bPlasQDCGainMatch_AIDA[32] ={0};
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();

  int ev = 0;
  for(AidaUnpackData& pInputD : pInputMain->Aida){
    ev++;
    AidaUnpackData* pInput = &pInputD;
    //pOutput->pAida = pAida;
    pOutput->pAIDA_WR = pInputMain->fAIDA_WR;

    AidaAnlData aida;
    if (pInput->ImplantEvents.size() > 1){

      //cout << " pInput->ImplantEvents.size() " << pInput->ImplantEvents.size() <<  endl;

      implantEvents++;

      // Cluster events on adjecent strips into one
      std::vector<AidaCluster> clusters = EventsToClusters(pInput->ImplantEvents);
      // Match front-back clusters which define physical hits on the detector
      std::vector<std::pair<AidaCluster, AidaCluster>> hits;
      // Track the number of implants in every DSSD
      std::vector<int> counts(conf->DSSDs(), 0);
      int max_dssd = 0;

      for (auto& i : clusters){
        if (i.DSSD == -1) continue;

        counts[i.DSSD - 1]++;
        if (i.DSSD > max_dssd) max_dssd = i.DSSD;

        if(i.Side != conf->DSSD(i.DSSD -1).XSide) continue;

        for (auto& j : clusters){
          if(j.DSSD != i.DSSD || j.Side != conf->DSSD(j.DSSD -1).YSide) continue;
          /// Gates (set in TAidaConfiguration)
          if (abs(i.Energy - j.Energy) < conf->FrontBackEnergyH() && i.IsGoodTime(j, conf->FrontBackWindow())){
            hits.push_back({i, j});
          }
        }
      }

      std::vector<int> channels(conf->FEEs() * 64);
      for (auto& i : pInput->ImplantEvents){
        channels[i.Module * 64 + i.Channel]++;
      }
      int channelM = 0;
      for (int i = 0; i < 64 * conf->FEEs(); ++i){
        if (channels[i]) ++channelM;
      }
      //implants_channels[0]->Fill(channelM);

      // Generate stored data for hits and plot the histograms

      for (auto& i : hits){
        AidaHit hit = ClusterPairToHit(i);
        hit.Event = ev;
        hit.Stopped = (hit.DSSD == max_dssd);

        // Check that every DSSD before has at least one implant event
        for(int j = hit.DSSD - 1; j > 0; j--){
          if (counts[j - 1] == 0) hit.Stopped = false;
        }

        if (hit.Stopped) stoppedEvents++;
        //  pInput->Implants.push_back(hit);
        //  pOutput->pAida.Implants.push_back(hit);

        aida.Implants.push_back(hit);
        implants_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
        if (hit.Stopped){
          implants_strip_xy_stopped[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
          implants_pos_xy_stopped[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        }

        // Helena
        if(hit.Stopped){
          if(hit.DSSD==1){
            pOutput->pt_lastIMP_DSSD1 = hit.TimeFront;
            pOutput->plastIMP_DSSD1_StripX = hit.StripX;
          }
          if(hit.DSSD==2){
            pOutput->pt_lastIMP_DSSD2 = hit.TimeFront;
            pOutput->plastIMP_DSSD2_StripX = hit.StripX;
          }
          if(hit.DSSD==3){
            pOutput->pt_lastIMP_DSSD3 = hit.TimeFront;
            pOutput->plastIMP_DSSD3_StripX = hit.StripX;
          }
        }

        //cout<<"ANALYSIS AIDA " <<pOutput->pEvent_Number<< " Energy " <<  hit.Energy << " hit.DSSD - 1 " <<hit.DSSD - 1 << endl;
        implants_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        implants_e[hit.DSSD - 1]->Fill(hit.Energy);

        //implants_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
        //implants_time_delta[hit.DSSD - 1]->Fill(hit.TimeFront - hit.TimeBack);
        implants_time_delta[hit.DSSD - 1]->Fill(hit.FastTimeFront - hit.FastTimeBack);

        implants_x_ex[hit.DSSD - 1]->Fill(hit.StripX, hit.EnergyFront);
        implants_y_ey[hit.DSSD - 1]->Fill(hit.StripY, hit.EnergyBack);

        int channel = i.first.Strip;
        implants_strip_1d[hit.DSSD - 1]->Fill(channel);
        channel = i.second.Strip + (conf->Wide() ? 386 : 128);
        implants_strip_1d[hit.DSSD - 1]->Fill(channel);
      }

      if (hits.size() > 0){
        implants_per_event[0]->Fill(hits.size());
        goodImplantEvents++;
      }
    }
    else if (pInput->DecayEvents.size() > 1){
      decayEvents++;

      std::vector<int> channels(conf->FEEs() * 64);

      #ifdef AIDA_PULSER_ALIGN
        int64_t wr_base = 0;
      #endif
      for (auto& i : pInput->DecayEvents){
        channels[i.Module * 64 + i.Channel]++;

        #ifdef AIDA_PULSER_ALIGN
          if(i.Module == 0 && i.Channel == 0) wr_base = i.Time;
        #endif
      }

      int channelM = 0;
      for (int i = 0; i < 768; ++i)
      if (channels[i]) ++channelM;

      //decays_channels[0]->Fill(channelM);

      if (channelM > 400){
        decayEvents--;
        pulserEvents++;
        #ifdef AIDA_PULSER_ALIGN
          if(channelM > 700){
          //if(pInputMain->fTrigger == 3)
            std::cout << "Identified a pulser event!" << std::endl;
            std::vector<int> offsets;
            offsets.resize(conf->FEEs());
            for (int i = 0; i < offsets.size(); i++) offsets[i] = 0;

            for (auto& i : pInput->DecayEvents){
              if (i.Energy < 1000) continue;
              std::cout << i.Module << ", " << i.Channel << ", " << i.Energy << ", " << std::hex << i.Time << "," << i.FastTime << std::dec << std::endl;
              int offset = (i.Time - wr_base) % 2000;
              if (offsets[i.Module] == 0){
                offsets[i.Module] = offset;
              }
              else if (offset > offsets[i.Module]){
                // confirm the offset is 2us out
                if (abs(offset - offsets[i.Module]) % 2000 != 0)
                std::cout << "LOGICAL MISTAKE IN ALIGNMENT" << std::endl;
              }
              else if (offset < offsets[i.Module]){
                if (abs(offset - offsets[i.Module]) % 2000 != 0)
                std::cout << "LOGICAL MISTAKE IN ALIGNMENT" << std::endl;
                offsets[i.Module] = offset;
              }

              aida_pulser_time->Fill(i.Module * 64 + i.Channel, offsets[i.Module]);
            }

            std::cout << std::endl;

            std::cout << "Put this into AIDA_Times.txt" << std::endl;
            for (int i  = 0; i < offsets.size(); i++){
              std::cout << i << " " << offsets[i] << std::endl;
            }

            throw TGo4UserException(3,"");
          }
        #endif
        return;
      }

      // Clean up huge event buffers - for now we just destroy them
      if (pInput->DecayEvents.size() > 400){
        decayEvents--;
        nonsenseEvents++;
        //pInput->SetValid(kFALSE); ///NEEDED!?
        return;
      }

      std::vector<AidaCluster> clusters = EventsToClusters(pInput->DecayEvents);

      std::vector<std::pair<AidaCluster, AidaCluster>> hits;
      std::unordered_map<aida_coord_t, int, aida_coord_hash> mults;

      for (auto& i : clusters){

        if(i.DSSD == -1 || i.Side != conf->DSSD(i.DSSD -1).XSide) continue;

        //if(i.Energy < 100) continue;
        for (auto& j : clusters){
          if(j.DSSD != i.DSSD || j.Side != conf->DSSD(j.DSSD -1).YSide) continue;
          //if(j.Energy < 100) continue;
          // Gates
          if (abs(i.Energy - j.Energy) < conf->FrontBackEnergyL() && i.IsGoodTime(j, conf->FrontBackWindow())){
            aida_coord_t x{i.DSSD, i.Side, i.Strip};
            aida_coord_t y{j.DSSD, j.Side, j.Strip};
            mults[x]++;
            mults[y]++;
            hits.push_back({i, j});
          }
        }
      }

      for (auto& i : hits){
        AidaHit hit = ClusterPairToHit(i);
        hit.Event = ev;

        aida_coord_t x{i.first.DSSD, i.first.Side, i.first.Strip};
        aida_coord_t y{i.second.DSSD, i.second.Side, i.second.Strip};
        if (mults[x] > 1 || mults[y] > 1)
        continue;

        //pInput->Decays.push_back(hit);
        //pOutput->pAida.push_back(hit); ///TEST
        //if (hit.DSSD != 3)
        //pOutput->pAida.Decays.push_back(hit);

        if (hit.DSSD != 3 || (hit.EnergyBack > 150 && hit.EnergyFront > 150 ))
        aida.Decays.push_back(hit);
        decays_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
        decays_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        decays_e[hit.DSSD - 1]->Fill(hit.Energy);
        decays_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
        decays_time_delta[hit.DSSD - 1]->Fill(hit.TimeFront - hit.TimeBack);
        //decays_time_delta[hit.DSSD - 1]->Fill(hit.FastTimeFront - hit.FastTimeBack);

        int channel = i.first.Strip;
        decays_strip_1d[hit.DSSD - 1]->Fill(channel);
        channel = i.second.Strip + (conf->Wide() ? 386 : 128);
        decays_strip_1d[hit.DSSD - 1]->Fill(channel);
      }

      if (clusters.size() > 0){
        decays_per_event[0]->Fill(clusters.size());
      }
      //cout <<" decayEvents " << decayEvents << endl;
    }
    else{
      nonsenseEvents++;
      //pInput->SetValid(kFALSE);
    }
    pOutput->pAida.push_back(aida);
  }
}


std::vector<AidaCluster> EventAnlProc::EventsToClusters(std::vector<AidaEvent> const& events){
  // track strip multiplicity and reject bad strips
  std::unordered_map<aida_coord_t, int, aida_coord_hash> stripm;

  std::vector<AidaCluster> clusters;
  for (auto& i : events){
    // Don't cluster invalid events
    if (i.DSSD == -1) continue;

    aida_coord_t coord{i.DSSD, i.Side, i.Strip};
    if(++stripm[coord] > 1)
      continue;

    bool added = false;
    AidaCluster* cluster;
    // Try to add the event to an existing cluster
    // Weird loop because we can remove clusters too
    auto it = std::begin(clusters);
    while (it != std::end(clusters)){
      auto& j = *it;
      if(j.IsAdjacent(i) && j.IsGoodTime(i)){
        // Add the event to the cluster
        if (!added){
          j.AddEvent(i);
          cluster = &j;
          added = true;
        }
        // If we match again the two clusters are merged into one
        // The old cluster is then removed
        else{
          cluster->AddCluster(j);
          it = clusters.erase(it);
          continue;
        }
      }
      ++it;
    }
    // Otherwise make a new cluster for the event
    if (!added){
      AidaCluster c_test;
      c_test.AddEvent(i);
      clusters.push_back(c_test);
    }
  }

  // remove bad strips again
  auto it = std::begin(clusters);
  while (it != std::end(clusters)){
    auto&j = *it;
    bool destroy = false;
    for(int s = j.StripMin; s <= j.StripMax; s++){
      aida_coord_t coord{j.DSSD, j.Side, s};
      if (stripm[coord] > 1){
        destroy = true;
      }
    }
    if (destroy){
      it = clusters.erase(it);
      continue;
    }
    ++it;
  }

  return clusters;
}

AidaHit EventAnlProc::ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const& i){
  AidaHit hit;

  hit.DSSD = i.first.DSSD;
  hit.StripX = i.first.Strip;
  hit.StripY = i.second.Strip;
  if (TAidaConfiguration::GetInstance()->Wide()){
    hit.PosX = 226.8 * i.first.Strip / 386. - 113.45;
  }
  else{
    hit.PosX = 75.6 * i.first.Strip / 128. - 37.75;
  }

  hit.PosY = 75.6 * i.second.Strip / 128. - 37.75;

  hit.StripXMin = i.first.StripMin;
  hit.StripXMax = i.first.StripMax;
  hit.StripYMin = i.second.StripMin;
  hit.StripYMax = i.second.StripMax;
  hit.ClusterSizeX  = i.first.N;
  hit.ClusterSizeY = i.second.N;

  //hit.Energy = (i.first.Energy + i.second.Energy) / 2;
  hit.EnergyFront = i.first.Energy;
  hit.EnergyBack = i.second.Energy;
  hit.Energy = hit.EnergyFront; // AIDA Triple n+n (back) has terrible energy resolution!

  hit.Time = std::min(i.first.Time, i.second.Time);
  hit.TimeFront = i.first.Time;
  hit.TimeBack = i.second.Time;
  hit.FastTime = std::min(i.first.FastTime, i.second.FastTime);
  hit.FastTimeFront = i.first.FastTime;
  hit.FastTimeBack = i.second.FastTime;

  return hit;
}


/**----------------------------------------------------------------------------------------------**/
/**-----------------------------------bplastic TWINPEAKS-----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
        
void EventAnlProc::Make_Plastic_Twinpeaks_Histos(){ 
        
  for (int i =1; i<4; i++){ 
    hbPlas_ToT_Sum_Slow[i] = MakeTH1('D', Form("bPlastic/Slow/ToT_Sum_Slow_Det.%2d",i), Form("bPlastic Sum ToT Slow Det. %2d",i), 4000, 0, 4000000);
                
    hbPlas_ToT_Sum_Fast[i] = MakeTH1('D', Form("bPlastic/Fast/ToT_Sum_Fast_Det.%2d",i), Form("bPlastic Sum ToT Fast Det. %2d",i), 4000, 0, 4000000);
  
    hbPlas_hit_pattern_det[i]= MakeTH1('D', Form("bPlastic/Stats/HitPattern_Det.%2d",i), Form("bPlastic Hit pattern Det. %2d",i), bPLASTIC_CHAN_PER_DET, 0, bPLASTIC_CHAN_PER_DET);
             
    for(int j=0; j<bPLASTIC_CHAN_PER_DET; j++){  
              
      hbPlas_Lead_T_Slow[i][j] = MakeTH1('D', Form("bPlastic/Slow/Lead/Lead T Slow Plas Det. %2d Ch.%2d",  i,j), Form("Lead - Time Det %2d Ch. %2d", i,j),2500, 0, 2000);
             
      hbPlas_Lead_T_Fast[i][j] = MakeTH1('D', Form("bPlastic/Fast/Lead/Lead T Fast Plas Det. %2d Ch.%2d",  i,j), Form("Lead - Time Det %2d Ch. %2d", i,j),2500, 0, 2000);
               
      if(i<3){//Take only cards with bPlast channels
        hbPlas_Lead_dT_coinc[i][j] = MakeTH1('D', Form("bPlastic/Lead-LeadCoincidenceChan/Lead dT Plas Det. %2d Ch %d Ref Coincidence",i,j), Form("Lead dT Plas Det. %2d Ch %d Ref Coincidence",i,j),500,-200,200);
      }
               
      hbPlas_ToT_det_Slow[i][j] = MakeTH1('D', Form("bPlastic/Slow/ToT/ToT Slow Plas Det. %2d Ch. %2d",  i,j), Form("ToT Slow Det. %2d Ch. %2d", i,j),40000, 0., 4000000.); 
            
      hbPlas_ToT_det_Fast[i][j] = MakeTH1('D', Form("bPlastic/Fast/ToT/ToT Fast Plas Det. %2d Ch. %2d",  i,j), Form("ToT Fast Det. %2d Ch. %2d", i,j),40000, 0., 4000000.); 
            
      hbPlas_Multiplicity_Chan[i][j] = MakeTH1('D', Form("bPlastic/Stats/Mulitplicity_per_Chan/bPlast Multiplicity Det. %2d Ch. %2d",  i,j), Form("ToT Det. %2d Ch. %2d", i,j),50, 0., 50.);   

      hbPlas_lead_lead_ref_det[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Ref/Lead-Lead Plas Det. %2d RefCh. - Ch. %d", i,j), Form("Lead Ref Ch.0 - Lead Det.%2d Ch. %2d", i,j),2500, -50000., 50000.);
        
      //hbPlas_lead_lead_ref_det2[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Ref/Lead-Lead Plas Det. %2d RefCh. %2d", i,j), Form("Lead Ref Ch.0 - Lead Det.%2d Ch. %2d", i,j),2500, -50000., 50000.);

      //hbPlas_lead_lead_ref_det3[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Ref/Lead-Lead Plas Det. %2d RefCh. %2d", i,j), Form("Lead Ref Ch.0 - Lead Det.%2d Ch. %2d", i,j),2500, -50000., 50000.);
          
      //hbPlas_lead_lead_gated[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Egated/Lead-Lead Egated Plas Det. %2d Ch. %2d",  i,j), Form("Lead - Lead Energy gated Det. %2d Ch.  %2d", i,j),2500, -50000., 50000.); 
            
	    //hbPlas_SC41L_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41-Lead_Plas/SC41_Lead Plas Det. %2d Ch.%02d", i,j), Form("SC41 Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);    
	    
      //hbPlas_fatimatamex_dT[i][j] = MakeTH1('D', Form("bPlastic/AdditionalChannels/FatTamChan-Lead_bPlas/(bPlast)FatTamChan_Lead bPlas Det.%2d Ch.%02d", i,j), Form("FATIMA tamex channel in bPlast: Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4000, -100000., 100000.); 
	    
	    //hbPlas_fatimavme_dT[i][j] = MakeTH1('D', Form("bPlastic/AdditionalChannels/FatVMEChan-Lead_bPlas/(bPlast)FatVMEChan_Lead bPlas Det.%2d Ch.%02d", i,j), Form("FATIMA VME channel in bPlast: Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4000, -100000., 100000.); 
		
	    //hbPlas_fatimavme_dT[i][j] = MakeTH1('D', Form("bPlastic/AdditionalChannels/SC41R_Anal-Lead_bPlas/SC41R_Ana_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41R Analogue Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4000, -100000., 100000.); 
		
	    hbPlas_SC41L_Digi_lead[i][j] = MakeTH1('D', Form("bPlastic/AdditionalChannels/SC41L_Digi-Lead_bPlas/SC41L_Digi_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41L Digital Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4000, -100000., 100000.); 
		
	    hbPlas_SC41R_Digi_lead[i][j] = MakeTH1('D', Form("bPlastic/AdditionalChannels/SC41R_Digi-Lead_bPlas/SC41R_Digi_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41R Digital Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4000, -100000., 100000.);    
        
      hbPlas_Ge_Trig_lead[i][j] = MakeTH1('D', Form("bPlastic/AdditionalChannels/Ge_Trig-Lead_bPlas/Ge_Trig_Lead bPlas Det. %2d Ch.%02d", i,j), Form("Ge Trigger - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4000, -100000., 100000.);    
        
    }
  }
        
  //New histogram looking at correlation between Fibre implanted channels 1 and 9 (2 and 10 if A=1->4 etc). H.M.A (don't blame me though...)
        
  //hFIMP_ToT_Correlation_Comb1 = MakeTH2('D', "bPlastic/FIMP_ToT_Correlation_Comb1", "ToT vs ToT for 2 FIMP channels, combination 1",500,0,100000,500,0,100000);
  
  //hFIMP_ToT_Correlation_Comb2 = MakeTH2('D', "bPlastic/FIMP_ToT_Correlation_Comb2", "ToT vs ToT for 2 FIMP channels, combination 2",500,0,100000,500,0,100000); 
                       
  //hSC41_Analogue_Tamex = MakeTH1('D',"bPlastic/SC41/Analogue L-R","SC41 Analogue L - R",4002, -100000., 100000.); 
  
  //hSC41_Digital_Tamex = MakeTH1('D',"bPlastic/SC41/Digital L-R","SC41 Analogue L - R",4002, -100000., 100000.);   
                   
  hbPlas_Multiplicity_Det1 = MakeTH1('D',"bPlastic/Stats/Multiplicity_Det1","bPlastic Multiplicity Det 1",32,0,32);

  hbPlas_Multiplicity_Det2 = MakeTH1('D',"bPlastic/Stats/Multiplicity_Det2","bPlastic Multiplicity Det 2",32,0,32); 

  hbPlas_ToT_Slow_vs_Fast_Det1=MakeTH2('D',"bPlastic/Fast_vs_Slow_bPlast_Det1Ch24","Fast ToT vs Slow ToT bPlast det1 Ch.24",3000,0,3000000,1000,0,1000000);
  
  hbPlas_ToT_Slow_vs_Fast_Det2=MakeTH2('D',"bPlastic/Fast_vs_Slow_bPlast_Det2Ch24","Fast ToT vs Slow ToT bPlast det2 Ch.24",3000,0,3000000,1000,0,1000000);
  
  // hbPlas_Multiplicity_Det3 = MakeTH1('D',"bPlastic/Stats/Multiplicity_Fibre","bPlastic Multiplicity Fibre",32,0,32);

}   


///--------------------------------------/**Process bplastic histos**/------------------------------------------///
void EventAnlProc::Process_Plastic_Twinpeaks_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput){   
    
  fired_det1=false, fired_det2=false;    
  ZERO_ARRAY(bPlas_tot_hits);

  for(int a=1; a<4; a++){
    for (int b = 0; b < bPLASTIC_CHAN_PER_DET; b++){  
      for(int k=0; k<bPLASTIC_TAMEX_HITS; k++){
        lead_bplas_fast[a][b][k]=0; 
        ToT_bplas_Slow[a][b][k] = 0;   
        ToT_bplas_Fast[a][b][k] = 0;   
      }
    }
  }
         
        
///--------------------------------------/**LEAD (Fast)**/------------------------------------------///
     
  ///Loop over channels 
  pOutput->pbPlasDetNum_Fast= pInput->fbPlasDetNum_Fast;
           
  for(int a=1; a<4; a++){ ///Detector number (this is crappy coding... A.K.M)
    pOutput->pbPlas_FastChan[a]= pInput->fbPlas_FastChan[a];
    pOutput->pbPlas_SlowChan[a]= pInput->fbPlas_SlowChan[a];

    for (int b = 0; b < bPLASTIC_CHAN_PER_DET; b++){  ///Channel number 

///--------------------------------------/**Plastic fast lead time**/------------------------------------------///

      if(pInput->fbPlast_Fast_Lead_N[a][b]<bPLASTIC_TAMEX_HITS){
        pOutput->pbPlas_Fast_Lead_N[a][b] = pInput->fbPlast_Fast_Lead_N[a][b]; 

        for(int j=0; j< bPLASTIC_TAMEX_HITS; j++){ ///Hits 
          // if(j<20){
            lead_bplas_fast[a][b][j] = pInput->fbPlast_Fast_Lead[a][b][j];  

            if(lead_bplas_fast[a][b][j]!=0){
            
              hits_bplas_lead_fast++;  
              pOutput->pbPlas_FastLeadT[a][b][j] = lead_bplas_fast[a][b][j]; 
              //cout<<"pOutput->pbPlas_FastChan[a] " <<pOutput->pbPlas_FastChan[a] << " a " << a <<" b " << b << "  pOutput->pbPlas_FastLeadT[a][b][j] " << pOutput->pbPlas_FastLeadT[a][b][j] <<  endl;    
              hbPlas_Lead_T_Fast[a][b]->Fill(lead_bplas_fast[a][b][j]);

              pOutput->pbPlas_FastLeadHits = hits_bplas_lead_fast; 
              //pOutput->pbPlas_FastLeadT_Avg = lead_bplas_fast[a][b][j]/hits_bplas_lead;

///--------------------------------------/**bplast channels and signal**/------------------------------------------///

              if(a==bPLASTIC_DOWNSTREAM_DET && pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_DOWN_COIN][0] !=0 && lead_bplas_fast[bPLASTIC_DOWNSTREAM_DET][b][j]>0){ // s452 a==1 is the downstream detector. Coincidence is in a=3, b=0.
                hbPlas_Lead_dT_coinc[bPLASTIC_DOWNSTREAM_DET][b]->Fill((lead_bplas_fast[bPLASTIC_DOWNSTREAM_DET][b][j] - pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_DOWN_COIN][0])*5); ///In ps
              }

              if(a==bPLASTIC_UPSTREAM_DET && pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_UP_COIN][0] != 0&&  lead_bplas_fast[bPLASTIC_UPSTREAM_DET][b][j]>0){ // s452 a==2 is the upstream detector. Coincidence is in a=3, b=1.
                hbPlas_Lead_dT_coinc[bPLASTIC_UPSTREAM_DET][b]->Fill((lead_bplas_fast[bPLASTIC_UPSTREAM_DET][b][j] - pInput->fbPlast_Fast_Lead[bPLASTIC_ADDITIONAL_CH_MOD][bPLASTIC_UP_COIN][0])*5); ///In ns    
              }      
            }
///--------------------------------------/**Plastic fast kead ref dT**/------------------------------------------///

            //if(i>15 && pInput->fbPlas_Lead_PMT[16][j]>0 && pInput->fbPlas_Lead_PMT[a][b][j]>0){   

            if(bPlas_RefCh0_Det1[j]>0 && lead_bplas_fast[1][b][j]>0){
              lead_lead_bplas_Ref1[b][j] = (bPlas_RefCh0_Det1[j] -  lead_bplas_fast[1][b][j])*CYCLE_TIME;
            }
            if(bPlas_RefCh0_Det2[j]>0 && lead_bplas_fast[2][b][j]>0){
              lead_lead_bplas_Ref2[b][j] = (bPlas_RefCh0_Det2[j] -  lead_bplas_fast[2][b][j])*CYCLE_TIME;
            }
            if(lead_lead_bplas_Ref1[b][j]!=0 && a==1) hbPlas_lead_lead_ref_det[1][b] ->Fill(lead_lead_bplas_Ref1[b][j]);

            //cout<<"ANL EVENT " << pInput->fevent_number << " a " << a << " b " << b << " j " << j << endl;
            if(lead_lead_bplas_Ref2[b][j]!=0 && a==2) hbPlas_lead_lead_ref_det[2][b] ->Fill(lead_lead_bplas_Ref2[b][j]);


            ///Reference signals
            ///Fatima tamex OR in bplast tamex
            //if(bPlas_TAM_FATTAM>0 && lead_bplas_fast[a][b][j]>0){
              //bPlas_fatimatamex_dT[a][b][j] = (bPlas_TAM_FATTAM -  lead_bplas_fast[a][b][j])*CYCLE_TIME;
              //hbPlas_fatimatamex_dT[a][b]->Fill(bPlas_fatimatamex_dT[a][b][j] );
            //}
            /////Fatima VME OR in bplast tamex
            //if(bPlas_TAM_FATVME>0 && lead_bplas_fast[a][b][j]>0){
              //bPlas_fatimavme_dT[a][b][j] = (bPlas_TAM_FATVME -  lead_bplas_fast[a][b][j])*CYCLE_TIME;
              //hbPlas_fatimavme_dT[a][b]->Fill(bPlas_fatimavme_dT[a][b][j] );
            //}
            ///SC41L in bplast tamex  
            if(bPlas_TAM_SC41L_DIG>0 && lead_bplas_fast[a][b][j]>0){

              SC41L_DIG_lead_bPlas[a][b][j] = (bPlas_TAM_SC41L_DIG -  lead_bplas_fast[a][b][j])*CYCLE_TIME;
              hbPlas_SC41L_Digi_lead[a][b]->Fill(SC41L_DIG_lead_bPlas[a][b][j] );
            }
            ///SC41R in bplast tamex
            if(bPlas_TAM_SC41R_DIG>0 && lead_bplas_fast[a][b][j]>0){
              SC41R_DIG_lead_bPlas[a][b][j] = (bPlas_TAM_SC41R_DIG -  lead_bplas_fast[a][b][j])*CYCLE_TIME;
              hbPlas_SC41R_Digi_lead[a][b]->Fill(SC41R_DIG_lead_bPlas[a][b][j] );
            }
            /// Ge trigger in bplastic
            if(bPlas_TAM_Ge_TRIG>0 && lead_bplas_fast[a][b][j]>0){
              Ge_Trig_lead_bPlas[a][b][j] = (bPlas_TAM_Ge_TRIG -  lead_bplas_fast[a][b][j])*CYCLE_TIME;
              hbPlas_Ge_Trig_lead[a][b]->Fill(Ge_Trig_lead_bPlas[a][b][j] );
            }
          //}
        }

///--------------------------------------/**Plastic Fast/Slow Trail**/------------------------------------------///  

        pOutput->pbPlast_Fast_Trail_N[a][b] = pInput->fbPlast_Fast_Trail_N[a][b]; 
        pOutput->pbPlast_Slow_Trail_N[a][b] = pInput->fbPlast_Slow_Trail_N[a][b]; 
        if(a<4 && b<bPLASTIC_CHAN_PER_DET){
          for(int j=0; j< bPLASTIC_TAMEX_HITS; j++){ ///Hits 
            // if(j<20){          
            //trail_bplas[a][b][j] = ;  
            //hbPlas_Trail_T[a][b]->Fill(trail_bplas[a][b][j]);
            if(pInput->fbPlast_Fast_Trail[a][b][j]!=0)  hits_bplas_trail_fast++;  
            pOutput->pbPlas_Fast_TrailT[a][b][j] = pInput->fbPlast_Fast_Trail[a][b][j];  
            if(pInput->fbPlast_Slow_Trail[a][b][j]!=0)  hits_bplas_trail_slow++;  
            pOutput->pbPlas_Slow_TrailT[a][b][j] = pInput->fbPlast_Slow_Trail[a][b][j];  
          }
        }

///--------------------------------------/**Plastic Fast ToT**/------------------------------------------///  

        for(int j=0; j< bPLASTIC_TAMEX_HITS; j++){ 

          if(pInput->fbPlast_Fast_Trail[a][b][j] >0 && pInput->fbPlast_Fast_Lead[a][b][j]>0){        
            ToT_bplas_Fast[a][b][j] = (pInput->fbPlast_Fast_Trail[a][b][j] - pInput->fbPlast_Fast_Lead[a][b][j]);   

            ///Correction for overflows 
            if(ABS(ToT_bplas_Fast[a][b][j]) >(double)(COARSE_CT_RANGE>>1)) {   
              ToT_bplas_Fast[a][b][j] = CYCLE_TIME*(ToT_bplas_Fast[a][b][j] + COARSE_CT_RANGE);    
            } 
            else{  
              ToT_bplas_Fast[a][b][j]= CYCLE_TIME*ToT_bplas_Fast[a][b][j];                         
            }    

            ///Gain matching  
            pOutput-> pbPlas_Fast_ToTCalib[a][b][j] =ToT_bplas_Fast[a][b][j];

            if(ToT_bplas_Fast[a][b][j]>0) {
              hbPlas_ToT_det_Fast[a][b] ->Fill(ToT_bplas_Fast[a][b][j]);   
              hbPlas_ToT_Sum_Fast[a]->Fill(ToT_bplas_Fast[a][b][j]);   
              hbPlas_hit_pattern_det[a]->Fill(b); 

              //bPlas_tot_hits++; 
              bPlas_tot_hits[a][b]++;      
              hbPlas_Multiplicity_Chan[a][b] ->Fill(bPlas_tot_hits[a][b]);  

              if(a==1) hbPlas_Multiplicity_Det1->Fill(bPlas_tot_hits[1][b]);
              if(a==2) hbPlas_Multiplicity_Det2->Fill(bPlas_tot_hits[2][b]);
              //if(a==3) hbPlas_Multiplicity_Det3->Fill(bPlas_tot_hits);                      
            }//ToT>0
          }//Lead+Trail>0         
        }//hits
      }//Limit hits loop

///--------------------------------------/**Plastic SLOW BRANCH**/------------------------------------------///  

      pOutput->pbPlas_Slow_Lead_N[a][b] = pInput->fbPlast_Slow_Lead_N[a][b]; 
      if(pInput->fbPlast_Slow_Lead_N[a][b]<bPLASTIC_TAMEX_HITS){
        ///Slow Lead
        for(int j=0; j< bPLASTIC_TAMEX_HITS; j++){ 
          lead_bplas_slow[a][b][j] = pInput->fbPlast_Slow_Lead[a][b][j]; 
          if(lead_bplas_slow[a][b][j]!=0){
            hbPlas_Lead_T_Slow[a][b]->Fill(lead_bplas_slow[a][b][j]);
            hits_bplas_lead_slow++;  
            pOutput->pbPlas_SlowLeadT[a][b][j] = lead_bplas_slow[a][b][j];    
            pOutput->pbPlas_SlowLeadHits = hits_bplas_lead_slow; 
          }

///--------------------------------------/**Plastic Slow ToT**/------------------------------------------///  

          if(pInput->fbPlast_Slow_Trail[a][b][j] >0 && pInput->fbPlast_Slow_Lead[a][b][j]>0){      
            ToT_bplas_Slow[a][b][j] = (pInput->fbPlast_Slow_Trail[a][b][j] - pInput->fbPlast_Slow_Lead[a][b][j]);   

            ///Correction for overflows 
            if(ABS(ToT_bplas_Slow[a][b][j]) >(double)(COARSE_CT_RANGE>>1)) {   
              ToT_bplas_Slow[a][b][j] = CYCLE_TIME*(ToT_bplas_Slow[a][b][j] + COARSE_CT_RANGE);    
            } 
            else{  
              ToT_bplas_Slow[a][b][j]= CYCLE_TIME*ToT_bplas_Slow[a][b][j];                         
            }    
            ///Gain matching  
            // pOutput-> pbPlas_ToTCalib[a][b][j] = fCal->Abplas_TAMEX_ZAoQ[i]* ToT_bplas[a][b][j] + fCal->Bbplas_TAMEX_ZAoQ[i];
            pOutput-> pbPlas_Slow_ToTCalib[a][b][j] =ToT_bplas_Slow[a][b][j];

            if(ToT_bplas_Slow[a][b][j]>0) {
              hbPlas_ToT_det_Slow[a][b] ->Fill(ToT_bplas_Slow[a][b][j]);   
              //cout<<"ToT_bplas_Slow[a][b][j] " <<ToT_bplas_Slow[a][b][j] << endl;          
              hbPlas_ToT_Sum_Slow[a]->Fill(ToT_bplas_Slow[a][b][j]);   
              if(pOutput-> pbPlas_Fast_ToTCalib[1][24][j]>0 && a==1 && b==24)hbPlas_ToT_Slow_vs_Fast_Det1->Fill(ToT_bplas_Slow[1][24][j],pOutput-> pbPlas_Fast_ToTCalib[1][24][j]);
              if(pOutput-> pbPlas_Fast_ToTCalib[2][24][j]>0 && a==2 && b==24)hbPlas_ToT_Slow_vs_Fast_Det2->Fill(ToT_bplas_Slow[2][24][j],pOutput-> pbPlas_Fast_ToTCalib[2][24][j]);
            }//ToT>0
          }//Lead+Trail>0               
        }///hits loop
      }///hit limit 
    }  ///Channel
  }  ///Detector 
} ///Function bPlast Twin Peaks



/**----------------------------------------------------------------------------------------------**/
/**-------------------------------------   Germanium   ------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_Germanium_Histos(){

  hGe_ESum = MakeTH1('I',"Germanium/Sum/Germanium_ESum_1keV","Germanium Energy Sum",20000,0,20000);

  hGe_ESum_halfkev = MakeTH1('I',"Germanium/Sum/Germanium_ESum_0_5keV","Germanium Energy Sum 0_5keV",10000,0,5000);

  hGe_Mult = MakeTH1('I',"Germanium/Stats/Germanium_Multiplicity","Germanium Multiplicity",30,0,30);

  hGe_SC41L = MakeTH1('I',"Germanium/SCI41/Germanium_SC41L_Ana","Germanium SC41 Left",5000,0,5000);

  hGe_SC41R = MakeTH1('I',"Germanium/SCI41/Germanium_SC41R_Ana","Germanium SC41 Right",5000,0,5000);

  hGe_SC41L_digi = MakeTH1('I',"Germanium/SCI41/Germanium_SC41L_Digi","Germanium SC41 Left",5000,0,5000);

  hGe_SC41R_digi = MakeTH1('I',"Germanium/SCI41/Germanium_SC41R_Digi","Germanium SC41 Right",5000,0,5000);

  hGe_ESum_largerange_OF = MakeTH1('I',"Germanium/Sum/Germanium_largerange_OF","Germanium Energy Sum (Overflow)",5000,0,5000);

  hGe_ESum_largerange_PU = MakeTH1('I',"Germanium/Sum/Germanium_largerange_PU","Germanium Energy Sum (Pileup)",5000,0,5000);

  hGe_Hit_Pat = MakeTH1('I',"Germanium/Stats/Germanium_Hit_Pat","Germanium Hit Pattern",Germanium_MAX_HITS,0,Germanium_MAX_HITS);

  hGe_Chan_E_Mat = MakeTH2('D',"Germanium/Sum/Germanium_E_Mat","Germanium Energy-Energy Matrix",2000,0,2000,2000,0,2000);

  hGe_Chan_E_vsDet = MakeTH2('D',"Germanium/Sum/Germanium_E_DetectorID","Germanium Energy vs Dectector Matrix",Germanium_MAX_HITS,0,Germanium_MAX_HITS,5000,0,5000);

  hGe_MultvsdT = MakeTH2('D',"Germanium/Stats/Germanium_Mult_vs_GamGamdT","Germanium Multiplicity vs Gamma Gamma dT",Germanium_MAX_HITS,0,Germanium_MAX_HITS,200,-1000,1000);

  hGe_AddbackSum = MakeTH1('I',"Germanium/Sum/Germanium_Addback_1keV","Germanium Addback Energy Sum 1keV",5000,0,5000);

  hGe_AddbackSum_halfkev = MakeTH1('I',"Germanium/Sum/Germanium_Addback_0_5keV","Germanium Addback Energy Sum 0.5keV",10000,0,5000);

  hGe_dTaddback = MakeTH1('I',"Germanium/Sum/Germanium_Addback_dT","Germanium Addback dT",100,-1000,1000);

  //hGe_dTgammagamma = MakeTH1('D',"Germanium/Sum/Germanium_GammaGamma_dT","Germanium Gamma-Gamma dT",200,-1000,1000);

  //hGe_CFdT_gammagamma = MakeTH1('D',"Germanium/Sum/Germanium_GammaGamma_CF_dT","Germanium Gamma-Gamma dT",2000,-1000,1000);

  hGe_Energy_GainMonitor= MakeTH2('D',"Germanium/Stats/Germanium_GainMonitor","Germanium Energy-Time Matrix",1240,16600,29000,2500,0,5000,"Time (minutes)","Ge Energy (keV)");

  for (int i=0; i<Germanium_MAX_DETS; i++){
    for (int j = 0; j < Germanium_CRYSTALS; j++){
          
      hGe_Chan_E[i][j] = MakeTH1('D',Form("Germanium/Energy_Ch_1keV/Germanium_E_Det_%1d_%1d",i, j), Form("Germanium Channel Energy Detector %1d Crystal %1d",i, j),5000,0,5000);

      ///I moved this back to unpack proc
      //hGe_ERaw[i][j]= MakeTH1('D',Form("Germanium/Raw/Germanium_ERaw_Det_%1d_%1d",i, j), Form("Germanium Channel Energy Detector %1d Crystal %1d",i, j),20000,0,10000);

      hGe_Chan_E_halfkev[i][j] = MakeTH1('D',Form("Germanium/Energy_Ch_0_5keV/Germanium_E_0_5keV_Det_%1d_%1d",i, j), Form("Germanium Channel 0.5keV Energy Detector %1d Crystal %1d",i, j),10000,0,5000);

      hGe_Chan_Time_Diff[i][j] = MakeTH1('D',Form("Germanium/Time_diff/Germanium_Chan_Time_Diff_Det_%1d_Chan_%1d",i,j), Form("Germanium Channel Time Difference for Detector %1d Channel %1d",i,j),200,-1000,1000);

      hGe_Chan_Time_Diff_CF[i][j] = MakeTH1('D',Form("Germanium/Time_diff_CF/Germanium_Chan_Time_CF_Diff_Det_%1d_Chan_%1d",i,j), Form("Germanium Channel Time Difference with Const Frac for Detector %1d Channel %1d",i,j),2000,-1000,1000);
    }
  }

  for (int i = 0; i < GeFired; i++){
      for (int j = 0; j < GeFired; j++){
        if (i!=j){

          hID_gamgam[i][j]= MakeTH2('D',Form("Germanium/Ge_Gamma_Gamma%1d_%1d",i, j),Form("Ge_Gamma_Gamma_all%1d_%1d",i, j),1240,16600,29000,2500,0,5000,Form("Ge%1d Energy (keV)",i),Form("Ge%1d Energy (keV)",j));

        }
      }
    }


    // hGe_Chan_E2[j] = MakeTH1('D',Form("Germanium/Germanium_Energy2/Germanium_E2%2d",j), Form("Germanium Channel Energy Channel %2d",j),5000,0,5000);

    //hGe_Chan_Egate[j] = MakeTH1('D',Form("Germanium/gated energy/Germanium_Egate%2d",j), Form("Germanium Channel Energy Channel %2d",j),5000,0,5000);

  //for (int k=0; k<32; k++){
    //hGe_Chan_Timedifference_new[k] = MakeTH1('D',Form("Germanium/Timediff_new/Germanium_Chan_T_Diff%2d",k), Form("Germanium Channel T Difference for %2d",k),2000,-1000,1000);
    // Ge_Time_Diff_vs_Energy[k] = MakeTH2('D',Form("Germanium/Germanium_dT_vs_Energy_Spectra/Germanium_dT_vs_E%2d",k), Form("Germanium Time Difference Vs Channel Energy Channel %2d",k),5000,0,5000,100,-1000,1000);
  //}

}

///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///
void EventAnlProc::Process_Germanium_Histos(EventAnlStore* pOutput){
  Ge_time_mins=0;
  if(Ge_WR>0) Ge_time_mins =(Ge_WR/60E9)-26900000;

  pOutput->pGe_WR=Ge_WR;
  Gam_mult=0;
  dT_Ge=0;
  //timeFIRST_Ge=0;

  dT_Addback=0;
  dT_Ge_cfd=0;
  //Process hits once
	
	RefTGe=0;
	RefCFDGe=0;


  for (int i = 0; i < GeFired; i++){
    int det = GeDet[i];
    int crys = GeCrys[i];

    //GAM GAM
    for (int j = 0; j < GeFired; j++){
      if (i!=j && GeE_Cal[i]>0 && GeE_Cal[j]>0){
        hID_gamgam[g]->Fill(GeE_Cal[i],GeE_Cal[j]); //GEEBART GAMGAM
        /*if (det[i]>det[j] && crys[i]>crys[j]){
          llT_array[0]= GeE_Cal[i];
          llT_array[1]= GeE_Cal[j];
          llT_array[2]= (GeT[i]-GeT[j]);
        }*/
      }
    }

    pOutput->pGePileUp[i]=GePileUp[i];
    pOutput->pGeOverFlow[i]=GeOverFlow[i];
    //Skip pileup/overflow events
    if (GePileUp[i] && (det!=Germanium_SC41_Det  && det!=Germanium_SC41_Det_Digi && det!=Germanium_TimeMachine_Det)){
      hGe_ESum_largerange_PU->Fill(GeE_Cal[i]);
      continue;
    }

    if (GeOverFlow[i]&& (det!=Germanium_SC41_Det  && det!=Germanium_SC41_Det_Digi && det!=Germanium_TimeMachine_Det)){
      hGe_ESum_largerange_OF->Fill(GeE_Cal[i]);
      continue;
    }

    if(det>-1){
      pOutput->pGe_Event_T[det][crys] = GeEventT[i];
      pOutput->pGe_T[det][crys] = Ge_Talign[i];
      pOutput->pGe_CF_T[det][crys] = Ge_cfd_Talign[i];
      pOutput->pGe_E[det][crys] = GeE_Cal[i];
      pOutput->pGe_E_Raw[det][crys] = GeE[i];

      //hGe_ERaw[det][crys]->Fill(GeE[i]);


      if(det!=Germanium_SC41_Det  && det!=Germanium_SC41_Det_Digi && det!=Germanium_TimeMachine_Det){
        hGe_ESum->Fill(GeE_Cal[i]);

        hGe_Hit_Pat->Fill(det * Germanium_CRYSTALS + crys);

        hGe_ESum_halfkev->Fill(GeE_Cal[i]);
        //cout<<"det " << det <<endl;
        hGe_Chan_E_vsDet->Fill(det,GeE_Cal[i]);

        if(Ge_time_mins>0 && GeE_Cal[i]>0)  hGe_Energy_GainMonitor->Fill(Ge_time_mins,GeE_Cal[i]);

        // printf("   it is det%02d %lf\n", det, GeE_Cal[i]);
        hGe_Chan_E[det][crys]->Fill(GeE_Cal[i]);
        hGe_Chan_E_halfkev[det][crys]->Fill(GeE_Cal[i]);

        //cout<<"2GeE_Cal[i] " << GeE_Cal[i] <<" det " << det << " crys " << crys <<endl;

      }
      // Cross correlations
      if ( det!=Germanium_TimeMachine_Det&& det!=Germanium_SC41_Det  && det!=Germanium_SC41_Det_Digi){
        Gam_mult++;

        //Time alignment (detector vs all)

        //Aligned dT
        //if(det!=0 && crys!=1){
        //if(RefTGe!=0){
          //dT_Align= Ge_Talign[i]-RefTGe;
          //hGe_dTgammagamma->Fill(dT_Align);
          //hGe_Chan_Time_Diff[det][crys]->Fill(dT_Align);
          //cout<<"2event " << event_number <<" det " << det << " crys " << crys << " dT_Align " << dT_Align <<" Ge_Talign[i] " << Ge_Talign[i] << " RefTGe " << R/*efTGe << " i " << i << endl;
        //}
        /*///Constant fraction time testing 03.02.21
        if(RefCFDGe!=0){
          //CFD Time alignment

          //Aligned CFD dT
          dT_CFD_Align= Ge_cfd_Talign[i]-RefCFDGe;
          hGe_CFdT_gammagamma->Fill(dT_CFD_Align);
          hGe_Chan_Time_Diff_CF[det][crys]->Fill(dT_CFD_Align);
        }*/

        for (int j = 0; j < GeFired; j++){
          if (i == j) continue;

          if(GeE_Cal[i]>0 && GeE_Cal[j]>0){
            //dT Germanium raw
            //dT_Ge= GeT[i]-GeT[j];
                 
            hGe_Mult->Fill(Gam_mult);
            hGe_MultvsdT->Fill(Gam_mult,dT_CFD_Align);

            ///Gamma-Gamma Time gate
            if((dT_CFD_Align)>fCorrel->GGe1_Ge2_Low && (dT_CFD_Align)< fCorrel->GGe1_Ge2_High && det!=Germanium_SC41_Det  && det!=Germanium_SC41_Det_Digi ){
              hGe_Chan_E_Mat->Fill(GeE_Cal[i], GeE_Cal[j]);
            }
          }
        }
      }
    }
  }


	for (int i = 0; i < GeFired; i++){
	
    if(GeDet[i]==1 && GeCrys[i]==0){ 
      if(GeE_Cal[i]>700){
      RefTGe= Ge_Talign[i];
      RefCFDGe= Ge_cfd_Talign[i];

        for (int j = 0; j < GeFired; j++){
        	if(j==i) continue;
        	if(GeE_Cal[j]>700){
          	dT_Align= Ge_Talign[j]-RefTGe;
            hGe_Chan_Time_Diff[GeDet[j]][GeCrys[j]]->Fill(dT_Align);
            dT_CFD_Align= Ge_cfd_Talign[j]-RefCFDGe;
            hGe_Chan_Time_Diff_CF[GeDet[j]][GeCrys[j]]->Fill(dT_CFD_Align);
      	  }
        }
      }
    }
	else continue;
  }

  ////ADDBACK
  static const long dT_addback = 100;
  static const double energy_thres_addback = 100;

  // Detector addback
  for (int i = 0; i < Germanium_MAX_DETS; i++){
    double E[Germanium_CRYSTALS] = { 0 };
    long T[Germanium_CRYSTALS] = { 0 };
    int n[Germanium_CRYSTALS] = { 0 };
    int crys[Germanium_CRYSTALS] = { 0 };
    int v = 0;
    for (int j = 0; j < Germanium_CRYSTALS; j++){
      if (pOutput->pGe_E[i][j] == 0) continue;
      bool added = false;

      // Try to addback to an existing hit
      for (int k = 0; k < v; k++){
        dT_Addback=T[k] - pOutput->pGe_T[i][j];
        hGe_dTaddback->Fill(dT_Addback);
        if (T[k] - pOutput->pGe_T[i][j] < dT_addback){
          hGe_dTaddback->Fill(T[k] - pOutput->pGe_T[i][j]);
          if(pOutput->pGe_E[i][j]>energy_thres_addback){
            E[k] += pOutput->pGe_E[i][j];
            T[k] = (T[k] + pOutput->pGe_T[i][j]) / (n[k] + 1);
            n[k] += 1;
            added = true;
          }
        }
      }

      // Add to a new hit
      if (!added){
        T[v] = pOutput->pGe_T[i][j];
        E[v] = pOutput->pGe_E[i][j];
        crys[v] = j;
        n[v] = 1;
        v++;
      }
    }
    // Fill and write to Tree the addback energies
    for (int j = 0; j < v; j++){

      pOutput->pGe_EAddback[i][crys[j]] = E[j];
      pOutput->pGe_T[i][crys[j]] = T[j];

      ///  SCI41 Signals
      if(pOutput->pGe_EAddback[Germanium_SC41_Det][Germanium_SC41L_Crystal]>0) hGe_SC41L->Fill(pOutput->pGe_EAddback[Germanium_SC41_Det][Germanium_SC41L_Crystal]);

      if(pOutput->pGe_EAddback[Germanium_SC41_Det][Germanium_SC41R_Crystal]>0) hGe_SC41R->Fill(pOutput->pGe_EAddback[Germanium_SC41_Det][Germanium_SC41R_Crystal]);

      if(pOutput->pGe_EAddback[Germanium_SC41_Det_Digi][Germanium_SC41R_Crystal_Digi]>0) hGe_SC41R_digi->Fill(pOutput->pGe_EAddback[Germanium_SC41_Det_Digi][Germanium_SC41R_Crystal_Digi]);

      if(pOutput->pGe_EAddback[Germanium_SC41_Det_Digi][Germanium_SC41L_Crystal_Digi]>0) hGe_SC41L_digi->Fill(pOutput->pGe_EAddback[Germanium_SC41_Det_Digi][Germanium_SC41L_Crystal_Digi]);

      if(GePileUp[i]==0 && GeOverFlow[i]==0 ){

        if(i!=Germanium_SC41_Det && i!=Germanium_SC41_Det_Digi && i!=Germanium_TimeMachine_Det){
          hGe_AddbackSum->Fill(E[j]);
          hGe_AddbackSum_halfkev->Fill(E[j]);
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------------------------//
// TH1I* EventAnlProc::MakeH1I(const char* fname,
//                             const char* hname,
//                             Int_t nbinsx,
//                             Float_t xmin, Float_t xmax,
//                             const char* xtitle,
//                             Color_t linecolor,
//                             Color_t fillcolor,
//                             const char* ytitle) {
// //    TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetHistogram, fname, hname);
// //    if (res!=0) return dynamic_cast<TH1I*>(res);
//
//    TH1I* histo = new TH1I(hname, hname, nbinsx, xmin, xmax);
//    histo->SetXTitle(xtitle);
//    if (ytitle) histo->SetYTitle(ytitle);
//    histo->SetLineColor(linecolor);
//    histo->SetFillColor(fillcolor);
//    AddHistogram(histo, fname);
//    return histo;
// }
//-----------------------------------------------------------------------------------------------------------------------------//

// TH2I* EventAnlProc::MakeH2I(const char* fname,
//                              const char* hname,
//                              Int_t nbinsx, Float_t xmin, Float_t xmax,
//                              Int_t nbinsy, Float_t ymin, Float_t ymax,
//                              const char* xtitle, const char* ytitle,
//                              Color_t markercolor) {
// //    TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetHistogram, fname, hname);
// //    if (res!=0) return dynamic_cast<TH2I*>(res);
//
//    TH2I* histo = new TH2I(hname, hname, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
//    histo->SetMarkerColor(markercolor);
//    histo->SetXTitle(xtitle);
//    histo->SetYTitle(ytitle);
//    AddHistogram(histo, fname);
//    return histo;
// }
//-----------------------------------------------------------------------------------------------------------------------------//

TGo4WinCond* EventAnlProc::MakeWindowCond(const char* fname,
                                           const char* cname,
                                           float left,
                                           float right,
                                           const char* HistoName) {
  // TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetAnalysisCondition, fname, cname);
   //if (res!=0) return dynamic_cast<TGo4WinCond*>(res);

   TGo4WinCond* cond = new TGo4WinCond((Text_t*)cname);
   cond->SetValues(left, right);
   cond->Enable();
   if (HistoName!=0)
     cond->SetHistogram(HistoName);
   AddAnalysisCondition(cond, fname);
   return cond;
}
//-----------------------------------------------------------------------------------------------------------------------------//

TGo4PolyCond* EventAnlProc::MakePolyCond(const char* fname,
                                          const char* cname,
                                          Int_t size,
                                          Float_t (*points)[2],
                                          const char* HistoName) {
   //TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetAnalysisCondition, fname, cname);
   //if (res!=0) return dynamic_cast<TGo4PolyCond*>(res);

   Float_t fullx[size+1], fully[size+1];
   int numpoints = size;

   for (int i=0;i<numpoints;i++) {
     fullx[i] = points[i][0];
     fully[i] = points[i][1];
   }

   // connect first and last points
   if ((fullx[0]!=fullx[numpoints-1]) || (fully[0]!=fully[numpoints-1])) {
      fullx[numpoints] = fullx[0];
      fully[numpoints] = fully[0];
      numpoints++;
   }

   TCutG mycat("initialcut", numpoints, fullx, fully);
   TGo4PolyCond* cond = new TGo4PolyCond((Text_t*)cname);
   cond->SetValues(&mycat);
   cond->Enable();
   if (HistoName!=0)
     cond->SetHistogram(HistoName);
   AddAnalysisCondition(cond, fname);
   return cond;
}
///-------------------------------------------------------------------------------------------------------
void EventAnlProc::FRS_Gates(){
  Int_t i;
  Int_t j;
  Int_t first;
  ifstream    file;
  string line;
  stringstream ss;


//   ///--------------------------------------------------------------------------------
// 
//       file.open("Configuration_Files/FRS/AoQ_Shift.txt");
//    while(file.good()){
//     getline(file,line,'\n');
//     if(line[0] == '#') continue;
//     sscanf(line.c_str(),"%f %f %f %f %f",&frs_wr_i,&frs_wr_j,&aoq_shift_tpc_value,&aoq_shift_sci21_value,&aoq_shift_sci22_value);
// 
//     FRS_WR_i[d]=frs_wr_i;
//     FRS_WR_j[d]=frs_wr_j;
//     AoQ_shift_TPC_value[d]=aoq_shift_tpc_value;
//     AoQ_shift_Sci21_value[d]=aoq_shift_sci21_value;
//     AoQ_shift_Sci22_value[d]=aoq_shift_sci22_value;
//     AoQ_Shift_array=d;
//        d++;
//      }
// 
//   file.close();

 ///--------------------------------------------------------------------------------  
   file.open("Configuration_Files/2D_Gates/ID_x2AoQ.txt");
   
    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
     if(line.empty()) break;
     if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
     }
     if(line.rfind("#",0) != 0){
        first = 1;
        ss.clear();
        ss << line;
        if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
           ss >> XX2_AoQ[j][i] >> YX2_AoQ[j][i];
        }
        else cout << "Warning: problem with ID_x2AoQ.txt." << endl;
        i++;
     }
    }
    num_ID_x2AoQ = j;
    file.close();   


 ///--------------------------------------------------------------------------------
  file.open("Configuration_Files/2D_Gates/ID_x4AoQ.txt");
  
    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
     if(line.empty()) break;
     if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
     }
     if(line.rfind("#",0) != 0){
        first = 1;
        ss.clear();
        ss << line;
        if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
           ss >> XX4_AoQ[j][i] >> YX4_AoQ[j][i];
        }
        else cout << "Warning: problem with ID_x4AoQ.txt." << endl;
        i++;
     }
    }
    num_ID_x4AoQ = j;
    file.close();    
  

  ///--------------------------------------------------------------------------------
   file.open("Configuration_Files/2D_Gates/ID_x4AoQ_mhtdc.txt");
   
    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
     if(line.empty()) break;  
     if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
     }
     if(line.rfind("#",0) != 0){
        first = 1;
        ss.clear();
        ss << line;
        if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
           ss >> XX4_AoQ_mhtdc[j][i] >> YX4_AoQ_mhtdc[j][i];
        }
        else cout << "Warning: problem with ID_x4AoQ_mhtdc." << endl;
        i++;
     }
    }
    num_ID_x4AoQ_mhtdc = j;
    file.close();     
   

   ///--------------------------------------------------------------------------------
    file.open("Configuration_Files/2D_Gates/ID_x2AoQ_mhtdc.txt");
    
    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
      if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
     }
     if(line.rfind("#",0) != 0){
        first = 1;
        ss.clear();
        ss << line;
        if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
           ss >> XX2_AoQ_mhtdc[j][i] >> YX2_AoQ_mhtdc[j][i];
        }
        else cout << "Warning: problem with ID_x2AoQ_mhtdc.txt." << endl;
        i++;
     }
    }
    num_ID_x2AoQ_mhtdc = j;
    file.close(); 

 ///--------------------------------------------------------------------------------
  file.open("Configuration_Files/2D_Gates/ID_Z_Z2.txt");
  
    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
      if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
     }
     if(line.rfind("#",0) != 0){
        first = 1;
        ss.clear();
        ss << line;
        if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
           ss >> X_ZZ2[j][i] >> Y_ZZ2[j][i];
        }
        else cout << "Warning: problem with ID_Z_Z2.txt." << endl;
        i++;
     }
    }
    num_ID_Z_Z2 = j;
    file.close();   
  

  ///--------------------------------------------------------------------------------
  file.open("Configuration_Files/2D_Gates/ID_Z_Z2_mhtdc.txt");
  
    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
      if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
     }
     if(line.rfind("#",0) != 0){
        first = 1;
        ss.clear();
        ss << line;
        if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
           ss >> X_ZZ2_mhtdc[j][i] >> Y_ZZ2_mhtdc[j][i];
        }
        else cout << "Warning: problem with ID_Z_Z2_mhtdc.txt." << endl;
        i++;
     }
    }
    num_ID_Z_Z2_mhtdc = j;
    file.close();  


 ///--------------------------------------------------------------------------------
    file.open("Configuration_Files/2D_Gates/ID_ZvsAoQ.txt");

    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
       if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
       }
       if(line.rfind("#",0) != 0){
          first = 1;
          ss.clear();
          ss << line;
          if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
             ss >> X_ZAoQ[j][i] >> Y_ZAoQ[j][i];
          }
        else cout << "Warning: problem with ID_ZvsAoQ.txt." << endl;
        i++;
     }
    }
    num_ID_Z_AoQ = j;
    file.close();


 ///------------------------------------Elif--------------------------------------------
    file.open("Configuration_Files/2D_Gates/ID_ZvsAoQ_mhtdc.txt");
      
    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
      if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
     }
     if(line.rfind("#",0) != 0){
        first = 1;
        ss.clear();
        ss << line;
        if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
           ss >> X_ZAoQ_mhtdc[j][i] >> Y_ZAoQ_mhtdc[j][i];
        }
        else cout << "Warning: problem with ID_ZvsAoQ_mhtdc.txt." << endl;
        i++;
     }
    }
    num_ID_Z_AoQ_mhtdc = j;
    file.close();      

  
  ///--------------------------------------------------------------------------------   
    file.open("Configuration_Files/2D_Gates/ID_dEvsZ.txt");

    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
       if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
       }
       if(line.rfind("#",0) != 0){
          first = 1;
          ss.clear();
          ss << line;
          if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
             ss >> X_dEvsZ[j][i] >> Y_dEvsZ[j][i];
          }
        else cout << "Warning: problem with ID_dEvsZ.txt." << endl;
        i++;
     }
    }
    num_ID_dEvsZ = j;
    file.close();


 ///--------------------------------------------------------------------------------
    file.open("Configuration_Files/2D_Gates/ID_dEvsZ_x4AoQ_p1.txt");

    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
       if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
       }
       if(line.rfind("#",0) != 0){
          first = 1;
          ss.clear();
          ss << line;
          if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
             ss >> X_dEvsZ_x4AoQ_p1[j][i] >> Y_dEvsZ_x4AoQ_p1[j][i];
          }
        else cout << "Warning: problem with ID_dEvsZ_x4AoQ_p1.txt." << endl;
        i++;
     }
    }
    num_ID_dEvsZ_x4AoQ_p1 = j;
    file.close();
    

 ///--------------------------------------------------------------------------------
    file.open("Configuration_Files/2D_Gates/ID_dEvsZ_x4AoQ_m1.txt");

    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
       if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
       }
       if(line.rfind("#",0) != 0){
          first = 1;
          ss.clear();
          ss << line;
          if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
             ss >> X_dEvsZ_x4AoQ_m1[j][i] >> Y_dEvsZ_x4AoQ_m1[j][i];
          }
        else cout << "Warning: problem with ID_dEvsZ_x4AoQ_m1.txt." << endl;
        i++;
     }
    }
    num_ID_dEvsZ_x4AoQ_m1 = j;
    file.close();
  

 ///--------------------------------------------------------------------------------
    file.open("Configuration_Files/2D_Gates/ID_dEvsBRho.txt");

    first = 0;
    i=0;
    j=0; // j is number of gates in the file
    ss.clear();
    while (getline (file,line)){
       if(line.empty()) break;
       if(line.rfind("#",0) == 0){
         i = 0; if(first>0)j++;
       }
       if(line.rfind("#",0) != 0){
          first = 1;
          ss.clear();
          ss << line;
          if(j<MAX_FRS_GATE && i<MAX_FRS_PolyPoints){
             ss >> X_dEvsBRho[j][i] >> Y_dEvsBRho[j][i];
          }
        else cout << "Warning: problem with ID_dEvsBRho.txt." << endl;
        i++;
     }
    }
    num_ID_dEvsBRho = j;
    file.close();
}

void EventAnlProc::FRS_GainMatching(){
  //Float_t frs_wr_a;
  //Float_t frs_wr_b;
  Float_t frs_wr_i;                    
  Float_t frs_wr_j;
  //Float_t z1_shift_value;
  //Float_t z2_shift_value;
  Float_t aoq_shift_value;
  Float_t aoq_shift_tpc_value;
  Float_t aoq_shift_sci21_value;
  Float_t aoq_shift_sci22_value;


  ifstream file;
  string line;
  float frs_wr_a, frs_wr_b, z1_shift_value, z2_shift_value;
  int line_num = 0;

  file.open("Configuration_Files/FRS/Z1_Z2_Shift.txt");

  while (getline(file, line)) {
    if(line[0] == '#') continue;
    line_num++;
    if (sscanf(line.c_str(), "%f %f %f %f", &frs_wr_a, &frs_wr_b, &z1_shift_value, &z2_shift_value) == 4) {
      cout << "Line " << line_num << ": " << frs_wr_a << ", " << frs_wr_b << ", " << z1_shift_value << ", " << z2_shift_value << endl;
      FRS_WR_a[line_num]=frs_wr_a;
      FRS_WR_b[line_num]=frs_wr_b;
      Z1_shift_value[line_num]=z1_shift_value;
      Z2_shift_value[line_num]=z2_shift_value;
    } else {
      cout << "Error reading line " << line_num << endl;
    }
  }
  Z_Shift_array=line_num;

  file.close();

  ///--------------------------------------------------------------------------------

  file.open("Configuration_Files/FRS/AoQ_Shift.txt");
  int d=0;
  while(file.good()){   
    getline(file,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%f %f %f %f %f",&frs_wr_i,&frs_wr_j,&aoq_shift_tpc_value,&aoq_shift_sci21_value,&aoq_shift_sci22_value);

    FRS_WR_i[d]=frs_wr_i;
    FRS_WR_j[d]=frs_wr_j;
    AoQ_shift_TPC_value[d]=aoq_shift_tpc_value;
    AoQ_shift_Sci21_value[d]=aoq_shift_sci21_value;
    AoQ_shift_Sci22_value[d]=aoq_shift_sci22_value;
    d++;  
  }
  AoQ_Shift_array=d;

  file.close();    
}

///-------------------------------------------------------------------------------------------------------
int EventAnlProc::IsData(ifstream &f) {
  char dum;
  char dumstr[300];
  int retval = 0;

  /*'operator >>' does not read End-of-line, therefore check if read character is not EOL (10) */
  do{
    dum = f.get();
    if (dum == '#')    // comment line => read whole line and throw it away
    f.getline(dumstr,300);
  }
  while ((dum == '#') || ((int)dum == 10));

  f.unget();   // go one character back
  retval = 1;
  return retval;
}
///-------------------------------------------------------------------------------------------------------

void EventAnlProc::get_used_systems(){
  for(int i = 0;i < 6;i++) Used_Systems[i] = false;

  ifstream data("Configuration_Files/DESPEC_General_Setup/Used_Systems.txt");
  if(data.fail()){
    cerr << "Could not find Used_Systems config file!" << endl;
    exit(0);
  }
  int i = 0;
  int id = 0;
  string line;
  char s_tmp[100];
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%s %d",s_tmp,&id);
    Used_Systems[i] = (id == 1);
    i++;
  }
  string DET_NAME[7] = {"FRS","AIDA","PLASTIC","FATIMA_VME","FATIMA_TAMEX","Germanium","FINGER"};
}

//-----------------------------------------------------------------------------------------------------------------------------//
//                                                            END                                                              //
//-----------------------------------------------------------------------------------------------------------------------------//
