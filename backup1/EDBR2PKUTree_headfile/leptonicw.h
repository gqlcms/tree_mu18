#ifndef _leptonicW_
#define _leptonicW_




void EDBR2PKUTree::leptonicW() { 
		W_pt              = Double_t(ptVlepJEC);
		W_eta             = Double_t(yVlep);
		W_phi             = Double_t(phiVlep);
		
		                mtVlepnew         = Double_t(sqrt(2*ptlep1*met*(1.0-cos(philep1-metPhi))));
                MTVlep            = Double_t(sqrt(2*ptlep1*MET_et*(1.0-cos(philep1-MET_phi))));
}

#endif