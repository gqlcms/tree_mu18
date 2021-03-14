#ifndef _resolved_MJJlv_1AK8_
#define _resolved_MJJlv_1AK8_

void EDBR2PKUTree::resolved_MJJlv_1AK8() {   
    if (Nj8==1){
    TLorentzVector ak8jet1_3,ak8jet2_3,Wlv;
    TLorentzVector ak8jet1_3_JEC_up,ak8jet2_3_JEC_up,Wlv_JEC_up;
    TLorentzVector ak8jet1_3_JEC_down,ak8jet2_3_JEC_down,Wlv_JEC_down;
    TLorentzVector ak8jet1_3_JER_up,ak8jet2_3_JER_up,Wlv_JER_up;
    TLorentzVector ak8jet1_3_JER_down,ak8jet2_3_JER_down,Wlv_JER_down;
        
        
    ak8jet1_3.SetPtEtaPhiM(jetAK8puppi_ptJEC,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi);
    ak8jet2_3.SetPtEtaPhiM(W_pt_resolved_1AK8,W_eta_resolved_1AK8,W_phi_resolved_1AK8,W_mass_resolved_1AK8);
    Wlv.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
    resolved_MJJlv = (ak8jet1_3+ak8jet2_3+Wlv).Mag();
        
    ak8jet1_3_JEC_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_JEC_up,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi);
    ak8jet2_3_JEC_up.SetPtEtaPhiM(W_pt_resolved_1AK8_JEC_up,W_eta_resolved_1AK8,W_phi_resolved_1AK8,W_mass_resolved_1AK8);
    Wlv_JEC_up.SetPtEtaPhiM(ptVlepJEC_JEC_up,yVlepJEC_JEC_up,phiVlepJEC_JEC_up,massVlepJEC_JEC_up);
    resolved_MJJlv_JEC_up = (ak8jet1_3_JEC_up+ak8jet2_3_JEC_up+Wlv_JEC_up).Mag();
        
    ak8jet1_3_JEC_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_JEC_down,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi);
    ak8jet2_3_JEC_down.SetPtEtaPhiM(W_pt_resolved_1AK8_JEC_down,W_eta_resolved_1AK8,W_phi_resolved_1AK8,W_mass_resolved_1AK8);
    Wlv_JEC_down.SetPtEtaPhiM(ptVlepJEC_JEC_down,yVlepJEC_JEC_down,phiVlepJEC_JEC_down,massVlepJEC_JEC_down);
    resolved_MJJlv_JEC_down = (ak8jet1_3_JEC_down+ak8jet2_3_JEC_down+Wlv_JEC_down).Mag();
    
        
    ak8jet1_3_JER_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_JER_up,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi);
    ak8jet2_3_JER_up.SetPtEtaPhiM(W_pt_resolved_1AK8_JER_up,W_eta_resolved_1AK8,W_phi_resolved_1AK8,W_mass_resolved_1AK8);
    Wlv_JER_up.SetPtEtaPhiM(ptVlepJEC_JER_up,yVlepJEC_JER_up,phiVlepJEC_JER_up,massVlepJEC_JER_up);
    resolved_MJJlv_JER_up = (ak8jet1_3_JER_up+ak8jet2_3_JER_up+Wlv_JER_up).Mag();
        
    ak8jet1_3_JER_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_JER_down,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi);
    ak8jet2_3_JER_down.SetPtEtaPhiM(W_pt_resolved_1AK8_JER_down,W_eta_resolved_1AK8,W_phi_resolved_1AK8,W_mass_resolved_1AK8);
    Wlv_JER_down.SetPtEtaPhiM(ptVlepJEC_JER_down,yVlepJEC_JER_down,phiVlepJEC_JER_down,massVlepJEC_JER_down);
    resolved_MJJlv_JER_down = (ak8jet1_3_JER_down+ak8jet2_3_JER_down+Wlv_JER_down).Mag();
    
    }
}

void EDBR2PKUTree::resolved_MJJlv_cal() {  
resolved_MJJlv = -99.;
resolved_MJJlv_1AK8();
}


#endif