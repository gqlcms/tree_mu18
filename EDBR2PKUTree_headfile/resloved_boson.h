#ifndef _resloved_boson_
#define _resloved_boson_


void EDBR2PKUTree::find_boson_W() {   
  // find AK4 pair which mass is closest to W boson's mass
  double W_PDG_mass = 80.4;
  double Z_PDG_mass = 91.2;
  double mass_distance = 1000;
  resolved_AK4_id[0] = -99;
  resolved_AK4_id[1] = -99;
    W_mass_resolved_1AK8 = -99;
    W_phi_resolved_1AK8 = -99;
    W_pt_resolved_1AK8 = -99;
    W_eta_resolved_1AK8 = -99;
    W_pt_resolved_1AK8_JEC_up = -99;
    W_pt_resolved_1AK8_JEC_down = -99;
    W_pt_resolved_1AK8_JER_up = -99;
    W_pt_resolved_1AK8_JER_down = -99;
    gen_w_resolved_id = -99;
   
  for(Int_t ii=0; ii<7; ii++)  {
       
    for(Int_t jj=ii+1; jj<8; jj++){
        if (pt_AK4_ex_DR0p8[jj]>0){
            TLorentzVector ak4_1,ak4_2,ak4_1_JEC_up,ak4_2_JEC_up,ak4_1_JEC_down,ak4_2_JEC_down,ak4_1_JER_up,ak4_2_JER_up,ak4_1_JER_down,ak4_2_JER_down,resolved_W,resolved_W_JEC_up,resolved_W_JEC_down,resolved_W_JER_up,resolved_W_JER_down;
            ak4_1.SetPtEtaPhiE(pt_AK4_ex_DR0p8[ii],eta_AK4_ex_DR0p8[ii],phi_AK4_ex_DR0p8[ii],e_AK4_ex_DR0p8[ii]);
            ak4_2.SetPtEtaPhiE(pt_AK4_ex_DR0p8[jj],eta_AK4_ex_DR0p8[jj],phi_AK4_ex_DR0p8[jj],e_AK4_ex_DR0p8[jj]);
            ak4_1_JEC_up.SetPtEtaPhiE(pt_AK4_ex_DR0p8_JEC_up[ii],eta_AK4_ex_DR0p8[ii],phi_AK4_ex_DR0p8[ii],e_AK4_ex_DR0p8[ii]);
            ak4_2_JEC_up.SetPtEtaPhiE(pt_AK4_ex_DR0p8_JEC_up[jj],eta_AK4_ex_DR0p8[jj],phi_AK4_ex_DR0p8[jj],e_AK4_ex_DR0p8[jj]);
            ak4_1_JEC_down.SetPtEtaPhiE(pt_AK4_ex_DR0p8_JEC_down[ii],eta_AK4_ex_DR0p8[ii],phi_AK4_ex_DR0p8[ii],e_AK4_ex_DR0p8[ii]);
            ak4_2_JEC_down.SetPtEtaPhiE(pt_AK4_ex_DR0p8_JEC_down[jj],eta_AK4_ex_DR0p8[jj],phi_AK4_ex_DR0p8[jj],e_AK4_ex_DR0p8[jj]);
            ak4_1_JER_up.SetPtEtaPhiE(pt_AK4_ex_DR0p8_JER_up[ii],eta_AK4_ex_DR0p8[ii],phi_AK4_ex_DR0p8[ii],e_AK4_ex_DR0p8[ii]);
            ak4_2_JER_up.SetPtEtaPhiE(pt_AK4_ex_DR0p8_JER_up[jj],eta_AK4_ex_DR0p8[jj],phi_AK4_ex_DR0p8[jj],e_AK4_ex_DR0p8[jj]);
            ak4_1_JER_down.SetPtEtaPhiE(pt_AK4_ex_DR0p8_JER_down[ii],eta_AK4_ex_DR0p8[ii],phi_AK4_ex_DR0p8[ii],e_AK4_ex_DR0p8[ii]);
            ak4_2_JER_down.SetPtEtaPhiE(pt_AK4_ex_DR0p8_JER_down[jj],eta_AK4_ex_DR0p8[jj],phi_AK4_ex_DR0p8[jj],e_AK4_ex_DR0p8[jj]);
            resolved_W = ak4_1 + ak4_2;
            resolved_W_JEC_up = ak4_1_JEC_up + ak4_2_JEC_up;
            resolved_W_JEC_down = ak4_1_JEC_down + ak4_2_JEC_down;
            resolved_W_JER_up = ak4_1_JER_up + ak4_2_JER_up;
            resolved_W_JER_down = ak4_1_JER_down + ak4_2_JER_down;
            if (fabs (resolved_W.Mag() - W_PDG_mass) < mass_distance){
                mass_distance = fabs (resolved_W.Mag() - W_PDG_mass);
                resolved_AK4_id[0] = ii;
                resolved_AK4_id[1] = jj;
                W_mass_resolved_1AK8 = resolved_W.Mag();
                W_phi_resolved_1AK8 = resolved_W.Phi();
                W_pt_resolved_1AK8 = resolved_W.Pt();
                W_eta_resolved_1AK8 = resolved_W.Eta();
                W_pt_resolved_1AK8_JEC_up = resolved_W_JEC_up.Pt();
                W_pt_resolved_1AK8_JEC_down = resolved_W_JEC_down.Pt();
                W_pt_resolved_1AK8_JER_up = resolved_W_JER_up.Pt();
                W_pt_resolved_1AK8_JER_down = resolved_W_JER_down.Pt();
                
                }
                        if (fabs (resolved_W.Mag() - Z_PDG_mass) < mass_distance){
                mass_distance = fabs (resolved_W.Mag() - W_PDG_mass);
                resolved_AK4_id[0] = ii;
                resolved_AK4_id[1] = jj;
                W_mass_resolved_1AK8 = resolved_W.Mag();
                W_phi_resolved_1AK8 = resolved_W.Phi();
                W_pt_resolved_1AK8 = resolved_W.Pt();
                W_eta_resolved_1AK8 = resolved_W.Eta();
                W_pt_resolved_1AK8_JEC_up = resolved_W_JEC_up.Pt();
                W_pt_resolved_1AK8_JEC_down = resolved_W_JEC_down.Pt();
                W_pt_resolved_1AK8_JER_up = resolved_W_JER_up.Pt();
                W_pt_resolved_1AK8_JER_down = resolved_W_JER_down.Pt();
                
                }
            }
        }
  }
}


// void EDBR2PKUTree::find_boson_W() {   
//   // find AK4 pair which mass is closest to W boson's mass
//   double W_PDG_mass = 80.4;
//   double Z_PDG_mass = 91.2;
//   double mass_distance = 1000;
//   double temp_mass_distance = 1000;
   
//   resolved_AK4_id[0] = -99;
//   resolved_AK4_id[1] = -99;
//     W_mass_resolved_1AK8 = -99;
//     W_phi_resolved_1AK8 = -99;
//     W_pt_resolved_1AK8 = -99;
//     W_eta_resolved_1AK8 = -99;
//     gen_w_resolved_id = -99;
   
//   for(Int_t ii=0; ii<7; ii++)  {
       
//     for(Int_t jj=ii+1; jj<8; jj++){
//         if (pt_AK4_ex_DR0p8[jj]>0){
//             TLorentzVector ak4_1,ak4_2,resolved_W;
//             ak4_1.SetPtEtaPhiE(pt_AK4_ex_DR0p8[ii],eta_AK4_ex_DR0p8[ii],phi_AK4_ex_DR0p8[ii],e_AK4_ex_DR0p8[ii]);
//             ak4_2.SetPtEtaPhiE(pt_AK4_ex_DR0p8[jj],eta_AK4_ex_DR0p8[jj],phi_AK4_ex_DR0p8[jj],e_AK4_ex_DR0p8[jj]);
//             resolved_W = ak4_1 + ak4_2;
//             temp_mass_distance = TMath::Min(fabs (resolved_W.Mag() - W_PDG_mass),fabs (resolved_W.Mag() - Z_PDG_mass));
//             if (temp_mass_distance < mass_distance){
//                 mass_distance = temp_mass_distance;
//                 resolved_AK4_id[0] = ii;
//                 resolved_AK4_id[1] = jj;
//                 W_mass_resolved_1AK8 = resolved_W.Mag();
//                 W_phi_resolved_1AK8 = resolved_W.Phi();
//                 W_pt_resolved_1AK8 = resolved_W.Pt();
//                 W_eta_resolved_1AK8 = resolved_W.Eta();
                
//                 }
//             }
//         }
//   }
// }

void EDBR2PKUTree::resolved_w_boson_gen_w_matching() {
// this part matches gen w pt with 2AK4's vector sum pt
 min_deltaR_resolved_1AK8 = 1000;
 gen_W_pt_resolved_1AK8 = -99.;
 gen_w_resolved_id =-99;
 
for (int i=0;i<5;i++){
    if (W_pt_resolved_1AK8>0.&&ptgenwl[i]>0.){
           double temp_deltaR_resolved_1AK8;
           temp_deltaR_resolved_1AK8 = sqrt(pow(fabs(etagenwl[i]-W_eta_resolved_1AK8),2)+pow(TMath::Min(fabs(phigenwl[i]-W_phi_resolved_1AK8),2*Pi-fabs(phigenwl[i]-W_phi_resolved_1AK8)),2));
           if (temp_deltaR_resolved_1AK8<min_deltaR_resolved_1AK8){
               min_deltaR_resolved_1AK8 = temp_deltaR_resolved_1AK8;
               gen_W_pt_resolved_1AK8 = ptgenwl[i];
               gen_w_resolved_id = i;
           }
    }
        }


}

void EDBR2PKUTree::resolved_w_boson_gen_w_daughter_matching() {
// this part matches gen w pt with 2AK4's vector sum pt
 min_deltaR_first_AK4_resolved_1AK8 = 1000;
 min_deltaR_second_AK4_resolved_1AK8 = 1000;
 first_AK4_matched_gen_w_daughter_pt = -99;
 first_AK4_matched_w_daughter_pt = -99;
 second_AK4_matched_gen_w_daughter_pt = -99;
 second_AK4_matched_w_daughter_pt = -99;
 
 
 double temp_deltaR_resolved_1AK8_1,temp_deltaR_resolved_1AK8_2,temp_deltaR_resolved_1AK8_3,temp_deltaR_resolved_1AK8_4;

    if (pt_AK4_ex_DR0p8[resolved_AK4_id[0]]>0.&&genw_q1_pt[gen_w_resolved_id]>0.&&genw_q2_pt[gen_w_resolved_id]>0.){
        
        first_AK4_matched_w_daughter_pt = pt_AK4_ex_DR0p8[resolved_AK4_id[0]];
           
           temp_deltaR_resolved_1AK8_1 = sqrt(pow(fabs(eta_AK4_ex_DR0p8[resolved_AK4_id[0]]-genw_q1_eta[gen_w_resolved_id]),2)+pow(TMath::Min(fabs(phi_AK4_ex_DR0p8[resolved_AK4_id[0]]-genw_q1_phi[gen_w_resolved_id]),2*Pi-fabs(phi_AK4_ex_DR0p8[resolved_AK4_id[0]]-genw_q1_phi[gen_w_resolved_id])),2));
           temp_deltaR_resolved_1AK8_2 = sqrt(pow(fabs(eta_AK4_ex_DR0p8[resolved_AK4_id[0]]-genw_q2_eta[gen_w_resolved_id]),2)+pow(TMath::Min(fabs(phi_AK4_ex_DR0p8[resolved_AK4_id[0]]-genw_q2_phi[gen_w_resolved_id]),2*Pi-fabs(phi_AK4_ex_DR0p8[resolved_AK4_id[0]]-genw_q2_phi[gen_w_resolved_id])),2));
           temp_deltaR_resolved_1AK8_3 = sqrt(pow(fabs(eta_AK4_ex_DR0p8[resolved_AK4_id[1]]-genw_q1_eta[gen_w_resolved_id]),2)+pow(TMath::Min(fabs(phi_AK4_ex_DR0p8[resolved_AK4_id[1]]-genw_q1_phi[gen_w_resolved_id]),2*Pi-fabs(phi_AK4_ex_DR0p8[resolved_AK4_id[1]]-genw_q1_phi[gen_w_resolved_id])),2));
           temp_deltaR_resolved_1AK8_4 = sqrt(pow(fabs(eta_AK4_ex_DR0p8[resolved_AK4_id[1]]-genw_q2_eta[gen_w_resolved_id]),2)+pow(TMath::Min(fabs(phi_AK4_ex_DR0p8[resolved_AK4_id[1]]-genw_q2_phi[gen_w_resolved_id]),2*Pi-fabs(phi_AK4_ex_DR0p8[resolved_AK4_id[1]]-genw_q2_phi[gen_w_resolved_id])),2));
           
           if (temp_deltaR_resolved_1AK8_1<temp_deltaR_resolved_1AK8_2){
               min_deltaR_first_AK4_resolved_1AK8 = temp_deltaR_resolved_1AK8_1;
               min_deltaR_second_AK4_resolved_1AK8 = temp_deltaR_resolved_1AK8_4;
               first_AK4_matched_gen_w_daughter_pt = genw_q1_pt[gen_w_resolved_id];
               second_AK4_matched_gen_w_daughter_pt = genw_q2_pt[gen_w_resolved_id];
           }
           if (temp_deltaR_resolved_1AK8_1>temp_deltaR_resolved_1AK8_2){
               min_deltaR_first_AK4_resolved_1AK8 = temp_deltaR_resolved_1AK8_2;
               min_deltaR_second_AK4_resolved_1AK8 = temp_deltaR_resolved_1AK8_3;
               first_AK4_matched_gen_w_daughter_pt = genw_q2_pt[gen_w_resolved_id];
               second_AK4_matched_gen_w_daughter_pt = genw_q1_pt[gen_w_resolved_id];
               
           }
    }
    



}

void EDBR2PKUTree::find_boson_Z() {   
   // find AK4 pair which mass is closest to W boson's mass
   double W_PDG_mass = 80.4;
   double mass_distance = 1000;
   resolved_AK4_id_Z[0] = -99;
   resolved_AK4_id_Z[1] = -99;
    Z_mass_resolved_1AK8 = -99;
    Z_phi_resolved_1AK8 = -99;
    Z_pt_resolved_1AK8 = -99;
    Z_eta_resolved_1AK8 = -99;
   
   for(Int_t ii=0; ii<7; ii++)  {
       
    for(Int_t jj=ii+1; jj<8; jj++){
        if (pt_AK4_ex_DR0p8[jj]>0){
            TLorentzVector ak4_1,ak4_2,resolved_W;
            ak4_1.SetPtEtaPhiE(pt_AK4_ex_DR0p8[ii],eta_AK4_ex_DR0p8[ii],phi_AK4_ex_DR0p8[ii],e_AK4_ex_DR0p8[ii]);
            ak4_2.SetPtEtaPhiE(pt_AK4_ex_DR0p8[jj],eta_AK4_ex_DR0p8[jj],phi_AK4_ex_DR0p8[jj],e_AK4_ex_DR0p8[jj]);
            resolved_W = ak4_1 + ak4_2;
            if (fabs (resolved_W.Mag() - W_PDG_mass) < mass_distance){
                mass_distance = fabs (resolved_W.Mag() - W_PDG_mass);
                resolved_AK4_id_Z[0] = ii;
                resolved_AK4_id_Z[1] = jj;
                Z_mass_resolved_1AK8 = resolved_W.Mag();
                Z_phi_resolved_1AK8 = resolved_W.Phi();
                Z_pt_resolved_1AK8 = resolved_W.Pt();
                Z_eta_resolved_1AK8 = resolved_W.Eta();
                
                }
            }
        }
   }
}

void EDBR2PKUTree::resolved_Z_boson_gen_Z_matching() {
// this part matches gen w pt with 2AK4's vector sum pt
 min_deltaR_resolved_1AK8_Z = 1000;
 gen_Z_pt_resolved_1AK8 = -99.;
 
 
for (int i=0;i<5;i++){
    if (Z_pt_resolved_1AK8>0.&&ptgenzl[i]>0.){
           double temp_deltaR_resolved_1AK8;
           temp_deltaR_resolved_1AK8 = sqrt(pow(fabs(etagenzl[i]-Z_eta_resolved_1AK8),2)+pow(TMath::Min(fabs(phigenzl[i]-Z_phi_resolved_1AK8),2*Pi-fabs(phigenzl[i]-Z_phi_resolved_1AK8)),2));
           if (temp_deltaR_resolved_1AK8<min_deltaR_resolved_1AK8_Z){
               min_deltaR_resolved_1AK8_Z = temp_deltaR_resolved_1AK8;
               gen_Z_pt_resolved_1AK8 = ptgenwl[i];
               
           }
    }
        }


}

void EDBR2PKUTree::use_w_mass_resolved(){
    find_boson_W();
    resolved_w_boson_gen_w_matching();
    resolved_w_boson_gen_w_daughter_matching();
}

void EDBR2PKUTree::use_z_mass_resolved(){
    find_boson_Z();
    resolved_Z_boson_gen_Z_matching();
}

#endif
