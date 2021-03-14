#ifndef _num_of_bJet_
#define _num_of_bJet_

void EDBR2PKUTree::num_of_bJet() {   
   //number of bjet calculation
		num_bJet=0.;
		num_bJet_loose=0.;
		num_bJet_tight=0.;
        num_bJet_deep=0.;
        num_bJet_loose_deep=0.;
        num_bJet_tight_deep=0.;

		for(Int_t i=0; i<8; i++)  {
                        deltaRAK4AK8_new[i]=0.;
                        deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
            if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));

			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.8838 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet=num_bJet+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.5426 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_loose=num_bJet_loose+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.9693 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_tight=num_bJet_tight+1;}
            if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb_bb[i])>0.6321 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_deep=num_bJet_deep+1;}
            if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb_bb[i])>0.2217 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_loose_deep=num_bJet_loose_deep+1;}
            if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb_bb[i])>0.8953 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_tight_deep=num_bJet_tight_deep+1;}

		}
		nbtag=num_bJet;
        nbtag_deep=num_bJet_deep;

		//number of bjet calculation Done
}

#endif