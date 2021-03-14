#ifndef _resolved_boson_0AK8_
#define _resolved_boson_0AK8_

/*
void EDBR2PKUTree::find_w_boson_0AK8() {   
    double W_PDG_mass = 80.4;
    double mass_distance = 1000;
    AK4_1_id_0AK8[0] = -99;
    AK4_2_id_0AK8[0] = -99;
    AK4_1_id_0AK8[1] = -99;
    AK4_2_id_0AK8[1] = -99;
    W_mass_resolved_1_0AK8 = -99;
    W_phi_resolved_1_0AK8 = -99;
    W_pt_resolved_1_0AK8 = -99;
    W_eta_resolved_1_0AK8 = -99;
    W_mass_resolved_2_0AK8 = -99;
    W_phi_resolved_2_0AK8 = -99;
    W_pt_resolved_2_0AK8 = -99;
    W_eta_resolved_2_0AK8 = -99;
    
    
    
    for(Int_t ii=0; ii<7; ii++)  {
       
    for(Int_t jj=ii+1; jj<8; jj++){
        if (ak4jet_pt[jj]>0){
            TLorentzVector ak4_1,ak4_2,resolved_W;
            ak4_1.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
            ak4_2.SetPtEtaPhiE(ak4jet_pt[jj],ak4jet_eta[jj],ak4jet_phi[jj],ak4jet_e[jj]);
            resolved_W = ak4_1 + ak4_2;
            if (fabs (resolved_W.Mag() - W_PDG_mass) < mass_distance){
                mass_distance = fabs (resolved_W.Mag() - W_PDG_mass);
                AK4_1_id_0AK8[0] = ii;
                AK4_1_id_0AK8[1] = jj;
                W_mass_resolved_1_0AK8 = resolved_W.Mag();
                W_phi_resolved_1_0AK8 = resolved_W.Phi();
                W_pt_resolved_1_0AK8 = resolved_W.Pt();
                W_eta_resolved_1_0AK8 = resolved_W.Eta();
                
                }
            }
        }
   }
   
   mass_distance = 1000;
   
       for(Int_t ii=0; ii<7; ii++)  {
       if (ii!=AK4_1_id_0AK8[0]&&ii!=AK4_1_id_0AK8[1]){
        for(Int_t jj=ii+1; jj<8; jj++){
            if (jj!=AK4_1_id_0AK8[0]&&jj!=AK4_1_id_0AK8[1]){
        if (ak4jet_pt[jj]>0){
            TLorentzVector ak4_1,ak4_2,resolved_W;
            ak4_1.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
            ak4_2.SetPtEtaPhiE(ak4jet_pt[jj],ak4jet_eta[jj],ak4jet_phi[jj],ak4jet_e[jj]);
            resolved_W = ak4_1 + ak4_2;
            if (fabs (resolved_W.Mag() - W_PDG_mass) < mass_distance){
                mass_distance = fabs (resolved_W.Mag() - W_PDG_mass);
                AK4_2_id_0AK8[0] = ii;
                AK4_2_id_0AK8[1] = jj;
                W_mass_resolved_2_0AK8 = resolved_W.Mag();
                W_phi_resolved_2_0AK8 = resolved_W.Phi();
                W_pt_resolved_2_0AK8 = resolved_W.Pt();
                W_eta_resolved_2_0AK8 = resolved_W.Eta();
                
                }
            }
        }
       }
}
   }
   
   
}
*/

void EDBR2PKUTree::find_w_boson_0AK8() {   
    double W_PDG_mass = 80.4;
   double Z_PDG_mass = 91.2;
   double mass_distance = 1000;
   double temp_mass_distance = 1000;
    AK4_1_id_0AK8[0] = -99;
    AK4_2_id_0AK8[0] = -99;
    AK4_1_id_0AK8[1] = -99;
    AK4_2_id_0AK8[1] = -99;
    W_mass_resolved_1_0AK8 = -99;
    W_phi_resolved_1_0AK8 = -99;
    W_pt_resolved_1_0AK8 = -99;
    W_eta_resolved_1_0AK8 = -99;
    W_mass_resolved_2_0AK8 = -99;
    W_phi_resolved_2_0AK8 = -99;
    W_pt_resolved_2_0AK8 = -99;
    W_eta_resolved_2_0AK8 = -99;
    
    
    
    for(Int_t ii=0; ii<7; ii++)  {
       
    for(Int_t jj=ii+1; jj<8; jj++){
        if (ak4jet_pt[jj]>0){
            TLorentzVector ak4_1,ak4_2,resolved_W;
            ak4_1.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
            ak4_2.SetPtEtaPhiE(ak4jet_pt[jj],ak4jet_eta[jj],ak4jet_phi[jj],ak4jet_e[jj]);
            resolved_W = ak4_1 + ak4_2;
            temp_mass_distance = TMath::Min(fabs (resolved_W.Mag() - W_PDG_mass),fabs (resolved_W.Mag() - Z_PDG_mass));
            if (temp_mass_distance < mass_distance){
                mass_distance = temp_mass_distance;
                AK4_1_id_0AK8[0] = ii;
                AK4_1_id_0AK8[1] = jj;
                W_mass_resolved_1_0AK8 = resolved_W.Mag();
                W_phi_resolved_1_0AK8 = resolved_W.Phi();
                W_pt_resolved_1_0AK8 = resolved_W.Pt();
                W_eta_resolved_1_0AK8 = resolved_W.Eta();
                
                }
            }
        }
   }
   
   mass_distance = 1000;
   temp_mass_distance = 1000;
   
       for(Int_t ii=0; ii<7; ii++)  {
       if (ii!=AK4_1_id_0AK8[0]&&ii!=AK4_1_id_0AK8[1]){
        for(Int_t jj=ii+1; jj<8; jj++){
            if (jj!=AK4_1_id_0AK8[0]&&jj!=AK4_1_id_0AK8[1]){
        if (ak4jet_pt[jj]>0){
            TLorentzVector ak4_1,ak4_2,resolved_W;
            ak4_1.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
            ak4_2.SetPtEtaPhiE(ak4jet_pt[jj],ak4jet_eta[jj],ak4jet_phi[jj],ak4jet_e[jj]);
            resolved_W = ak4_1 + ak4_2;
            temp_mass_distance = TMath::Min(fabs (resolved_W.Mag() - W_PDG_mass),fabs (resolved_W.Mag() - Z_PDG_mass));
            if (temp_mass_distance < mass_distance){
                mass_distance = temp_mass_distance;
                AK4_2_id_0AK8[0] = ii;
                AK4_2_id_0AK8[1] = jj;
                W_mass_resolved_2_0AK8 = resolved_W.Mag();
                W_phi_resolved_2_0AK8 = resolved_W.Phi();
                W_pt_resolved_2_0AK8 = resolved_W.Pt();
                W_eta_resolved_2_0AK8 = resolved_W.Eta();
                
                }
            }
        }
       }
}
   }
   
   
}


void EDBR2PKUTree::resolved_w_boson_gen_w_matching_0AK8() {
// this part matches gen w pt with 2AK4's vector sum pt
 min_deltaR_resolved_w_1_0AK8 = 1000;
 min_deltaR_resolved_w_2_0AK8 = 1000;
 gen_W_1_pt_resolved_0AK8 = -99.;
 gen_W_2_pt_resolved_0AK8 = -99.;
 gen_w_1_resolved_id =-99;
 gen_w_2_resolved_id =-99;
 
for (int i=0;i<5;i++){
    if (W_pt_resolved_1_0AK8>0.&&ptgenwl[i]>0.){
           double temp_deltaR_resolved_1AK8;
           temp_deltaR_resolved_1AK8 = sqrt(pow(fabs(etagenwl[i]-W_eta_resolved_1_0AK8),2)+pow(TMath::Min(fabs(phigenwl[i]-W_phi_resolved_1_0AK8),2*Pi-fabs(phigenwl[i]-W_phi_resolved_1_0AK8)),2));
           if (temp_deltaR_resolved_1AK8<min_deltaR_resolved_w_1_0AK8){
               min_deltaR_resolved_w_1_0AK8 = temp_deltaR_resolved_1AK8;
               gen_W_1_pt_resolved_0AK8 = ptgenwl[i];
               gen_w_1_resolved_id = i;
           }
    }
        }
        
for (int i=0;i<5;i++){
    if (i!=gen_w_1_resolved_id){
    if (W_pt_resolved_2_0AK8>0.&&ptgenwl[i]>0.){
           double temp_deltaR_resolved_1AK8;
           temp_deltaR_resolved_1AK8 = sqrt(pow(fabs(etagenwl[i]-W_eta_resolved_2_0AK8),2)+pow(TMath::Min(fabs(phigenwl[i]-W_phi_resolved_2_0AK8),2*Pi-fabs(phigenwl[i]-W_phi_resolved_2_0AK8)),2));
           if (temp_deltaR_resolved_1AK8<min_deltaR_resolved_w_2_0AK8){
               min_deltaR_resolved_w_2_0AK8 = temp_deltaR_resolved_1AK8;
               gen_W_2_pt_resolved_0AK8 = ptgenwl[i];
               gen_w_2_resolved_id = i;
           }
    }
        }
}


}

void EDBR2PKUTree::find_z_boson_0AK8() {   
    double W_PDG_mass = 80.4;
    double mass_distance = 1000;
    AK4_1_id_0AK8_Z[0] = -99;
    AK4_2_id_0AK8_Z[0] = -99;
    AK4_1_id_0AK8_Z[1] = -99;
    AK4_2_id_0AK8_Z[1] = -99;
    Z_mass_resolved_1_0AK8 = -99;
    Z_phi_resolved_1_0AK8 = -99;
    Z_pt_resolved_1_0AK8 = -99;
    Z_eta_resolved_1_0AK8 = -99;
    Z_mass_resolved_2_0AK8 = -99;
    Z_phi_resolved_2_0AK8 = -99;
    Z_pt_resolved_2_0AK8 = -99;
    Z_eta_resolved_2_0AK8 = -99;
    
    
    
    for(Int_t ii=0; ii<7; ii++)  {
       
    for(Int_t jj=ii+1; jj<8; jj++){
        if (ak4jet_pt[jj]>0){
            TLorentzVector ak4_1,ak4_2,resolved_W;
            ak4_1.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
            ak4_2.SetPtEtaPhiE(ak4jet_pt[jj],ak4jet_eta[jj],ak4jet_phi[jj],ak4jet_e[jj]);
            resolved_W = ak4_1 + ak4_2;
            if (fabs (resolved_W.Mag() - W_PDG_mass) < mass_distance){
                mass_distance = fabs (resolved_W.Mag() - W_PDG_mass);
                AK4_1_id_0AK8_Z[0] = ii;
                AK4_1_id_0AK8_Z[1] = jj;
                Z_mass_resolved_1_0AK8 = resolved_W.Mag();
                Z_phi_resolved_1_0AK8 = resolved_W.Phi();
                Z_pt_resolved_1_0AK8 = resolved_W.Pt();
                Z_eta_resolved_1_0AK8 = resolved_W.Eta();
                
                }
            }
        }
   }
   
       for(Int_t ii=0; ii<7; ii++)  {
       if (ii!=AK4_1_id_0AK8_Z[0]&&ii!=AK4_1_id_0AK8_Z[0]){
        for(Int_t jj=ii+1; jj<8; jj++){
            if (jj!=AK4_1_id_0AK8_Z[0]&&jj!=AK4_1_id_0AK8_Z[0]){
        if (ak4jet_pt[jj]>0){
            TLorentzVector ak4_1,ak4_2,resolved_W;
            ak4_1.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
            ak4_2.SetPtEtaPhiE(ak4jet_pt[jj],ak4jet_eta[jj],ak4jet_phi[jj],ak4jet_e[jj]);
            resolved_W = ak4_1 + ak4_2;
            if (fabs (resolved_W.Mag() - W_PDG_mass) < mass_distance){
                mass_distance = fabs (resolved_W.Mag() - W_PDG_mass);
                AK4_2_id_0AK8_Z[0] = ii;
                AK4_2_id_0AK8_Z[1] = jj;
                Z_mass_resolved_2_0AK8 = resolved_W.Mag();
                Z_phi_resolved_2_0AK8 = resolved_W.Phi();
                Z_pt_resolved_2_0AK8 = resolved_W.Pt();
                Z_eta_resolved_2_0AK8 = resolved_W.Eta();
                
                }
            }
        }
       }
}
   }
   
   
}

void EDBR2PKUTree::resolved_z_boson_gen_z_matching_0AK8() {
// this part matches gen w pt with 2AK4's vector sum pt
 min_deltaR_resolved_z_1_0AK8 = 1000;
 min_deltaR_resolved_z_2_0AK8 = 1000;
 gen_Z_1_pt_resolved_0AK8 = -99.;
 gen_Z_2_pt_resolved_0AK8 = -99.;
 gen_Z_1_resolved_id =-99;
 gen_Z_2_resolved_id =-99;
 
for (int i=0;i<5;i++){
    if (Z_pt_resolved_1_0AK8>0.&&ptgenwl[i]>0.){
           double temp_deltaR_resolved_1AK8;
           temp_deltaR_resolved_1AK8 = sqrt(pow(fabs(etagenzl[i]-Z_eta_resolved_1_0AK8),2)+pow(TMath::Min(fabs(phigenzl[i]-Z_phi_resolved_1_0AK8),2*Pi-fabs(phigenzl[i]-Z_phi_resolved_1_0AK8)),2));
           if (temp_deltaR_resolved_1AK8<min_deltaR_resolved_z_1_0AK8){
               min_deltaR_resolved_z_1_0AK8 = temp_deltaR_resolved_1AK8;
               gen_Z_1_pt_resolved_0AK8 = ptgenzl[i];
               gen_Z_1_resolved_id = i;
           }
    }
        }
        
for (int i=0;i<5;i++){
    if (i!=gen_Z_1_resolved_id){
    if (Z_pt_resolved_2_0AK8>0.&&ptgenwl[i]>0.){
           double temp_deltaR_resolved_1AK8;
           temp_deltaR_resolved_1AK8 = sqrt(pow(fabs(etagenzl[i]-Z_eta_resolved_2_0AK8),2)+pow(TMath::Min(fabs(phigenzl[i]-Z_phi_resolved_2_0AK8),2*Pi-fabs(phigenzl[i]-Z_phi_resolved_2_0AK8)),2));
           if (temp_deltaR_resolved_1AK8<min_deltaR_resolved_z_2_0AK8){
               min_deltaR_resolved_z_2_0AK8 = temp_deltaR_resolved_1AK8;
               gen_Z_2_pt_resolved_0AK8 = ptgenzl[i];
               gen_Z_2_resolved_id = i;
           }
    }
        }
}


}

void EDBR2PKUTree::use_w_mass_resolved_0AK8(){
    find_w_boson_0AK8();
    resolved_w_boson_gen_w_matching_0AK8();
}

void EDBR2PKUTree::use_z_mass_resolved_0AK8(){
    find_z_boson_0AK8();
    resolved_z_boson_gen_z_matching_0AK8();
}

#endif