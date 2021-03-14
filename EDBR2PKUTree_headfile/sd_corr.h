void EDBR2PKUTree::sd_corr() {
//puppi+softdrop recorrected by Thea's "JEC"
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;
        jet_mass_puppi_corr=-99;jet_mass_puppi_corr_2=-99;jet_mass_puppi_corr_3=-99;
//                gencorrect=1.0-0.321*pow(jet_pt_puppi*0.0354,-1.1);
//                recocorrect_0eta1p3=1.09-1.69e-04*jet_pt_puppi+3.34e-07*pow(jet_pt_puppi,2)-2.47e-10*pow(jet_pt_puppi,3)+7.8e-14*pow(jet_pt_puppi,4)-8.83e-18*pow(jet_pt_puppi,5);
//                recocorrect_1p3eta2p5=1.3-7.76e-04*jet_pt_puppi+1.11e-06*pow(jet_pt_puppi,2)-6.79e-10*pow(jet_pt_puppi,3)+1.87e-13*pow(jet_pt_puppi,4)-1.9e-17*pow(jet_pt_puppi,5);
                gencorrect=1.006-1.062*pow(jet_pt_puppi*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jet_pt_puppi+3.449e-07*pow(jet_pt_puppi,2)-2.681e-10*pow(jet_pt_puppi,3)+8.674e-14*pow(jet_pt_puppi,4)-1.001e-17*pow(jet_pt_puppi,5);
                recocorrect_1p3eta2p5=1.272-5.72e-04*jet_pt_puppi+8.37e-07*pow(jet_pt_puppi,2)-5.204e-10*pow(jet_pt_puppi,3)+1.454e-13*pow(jet_pt_puppi,4)-1.504e-17*pow(jet_pt_puppi,5);
                if (fabs(jetAK8puppi_eta)<=1.3){jet_mass_puppi_corr=jet_mass_puppi_un*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta)<2.5 && fabs(jetAK8puppi_eta)>1.3){jet_mass_puppi_corr=jet_mass_puppi_un*gencorrect*recocorrect_1p3eta2p5;}


Double_t gencorrect_2=1.0;
                Double_t recocorrect_0eta1p3_2=1.0;
                Double_t recocorrect_1p3eta2p5_2=1.0;
gencorrect_2=1.006-1.062*pow(jet_pt_puppi_2*0.08,-1.2);
                recocorrect_0eta1p3_2=1.093-1.501e-04*jet_pt_puppi_2+3.449e-07*pow(jet_pt_puppi_2,2)-2.681e-10*pow(jet_pt_puppi_2,3)+8.674e-14*pow(jet_pt_puppi_2,4)-1.001e-17*pow(jet_pt_puppi_2,5);
                recocorrect_1p3eta2p5_2=1.272-5.72e-04*jet_pt_puppi_2+8.37e-07*pow(jet_pt_puppi_2,2)-5.204e-10*pow(jet_pt_puppi_2,3)+1.454e-13*pow(jet_pt_puppi_2,4)-1.504e-17*pow(jet_pt_puppi_2,5);
                if (fabs(jetAK8puppi_eta_2)<=1.3){jet_mass_puppi_corr_2=jet_mass_puppi_un_2*gencorrect_2*recocorrect_0eta1p3_2;}
                else if (fabs(jetAK8puppi_eta_2)<2.5 && fabs(jetAK8puppi_eta_2)>1.3){jet_mass_puppi_corr_2=jet_mass_puppi_un_2*gencorrect_2*recocorrect_1p3eta2p5_2;}

        Double_t gencorrect_3=1.0;
        Double_t recocorrect_0eta1p3_3=1.0;
        Double_t recocorrect_1p3eta2p5_3=1.0;
        gencorrect_3=1.006-1.062*pow(jet_pt_puppi_3*0.08,-1.2);
        recocorrect_0eta1p3_3=1.093-1.501e-04*jet_pt_puppi_3+3.449e-07*pow(jet_pt_puppi_3,2)-2.681e-10*pow(jet_pt_puppi_3,3)+8.674e-14*pow(jet_pt_puppi_3,4)-1.001e-17*pow(jet_pt_puppi_3,5);
        recocorrect_1p3eta2p5_3=1.272-5.72e-04*jet_pt_puppi_3+8.37e-07*pow(jet_pt_puppi_3,2)-5.204e-10*pow(jet_pt_puppi_3,3)+1.454e-13*pow(jet_pt_puppi_3,4)-1.504e-17*pow(jet_pt_puppi_3,5);
        if (fabs(jetAK8puppi_eta_3)<=1.3){jet_mass_puppi_corr_3=jet_mass_puppi_un_3*gencorrect_3*recocorrect_0eta1p3_3;}
        else if (fabs(jetAK8puppi_eta_3)<2.5 && fabs(jetAK8puppi_eta_3)>1.3){jet_mass_puppi_corr_3=jet_mass_puppi_un_3*gencorrect_3*recocorrect_1p3eta2p5_3;}
        jet_mass_puppi_corr=jetAK8puppi_sd;
        jet_mass_puppi_corr_2=jetAK8puppi_sd_2;
        jet_mass_puppi_corr_3=jetAK8puppi_sd_3;
}