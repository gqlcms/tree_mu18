#ifndef _lepton_
#define _lepton_

void EDBR2PKUTree::lepton() {
        pfMET             = Double_t(met);
		pfMETPhi          = Double_t(metPhi);
    
        l_pt              = Double_t(ptlep1);
        l_pt_2              = Double_t(ptlep2);

		l_eta             = Double_t(etalep1);
		l_phi             = Double_t(philep1);   
}

#endif