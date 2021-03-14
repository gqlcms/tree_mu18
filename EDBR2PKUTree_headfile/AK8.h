#ifndef _AK8Jet_
#define _AK8Jet_
void EDBR2PKUTree::AK8Jet() {            
                jet_mass_puppi    = Double_t(jetAK8puppi_sd);
                jet_mass_puppi_un = Double_t(jetAK8puppi_sd);
                jet_pt_puppi      = Double_t(jetAK8puppi_ptJEC);
                jet_mass_puppi_corr = Double_t(jetAK8puppi_sd);
                //		jet_mass_softdrop = Double_t(sdropJEC);
                
                		jet_mass_puppi_2    = Double_t(jetAK8puppi_sd_2);
                jet_mass_puppi_un_2 = Double_t(jetAK8puppi_sd_2);
                jet_pt_puppi_2      = Double_t(jetAK8puppi_ptJEC_2);
                jet_mass_puppi_corr_2 = Double_t(jetAK8puppi_sd_2);

        
        jet_mass_puppi_3    = Double_t(jetAK8puppi_sd_3);
        jet_mass_puppi_un_3 = Double_t(jetAK8puppi_sd_3);
        jet_pt_puppi_3      = Double_t(jetAK8puppi_ptJEC_3);
        
        
        
                
}
#endif