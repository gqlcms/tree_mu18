#ifndef _Mj_oder_
#define _Mj_oder_

void EDBR2PKUTree::Mj_oder() {   
           Mj_min=-99; Mj_max=-99;pt_min=-99; pt_max=-99; t21_max=-99; t31_max=-99; t41_max=-99;  t21_min=-99; t31_min=-99; t41_min=-99; PTw2divPTw1=-99; PTw3divPTw1=-99; PTw3divPTw2=-99;phi_max=-99;phi_min=-99;eta_max=-99;eta_min=-99;t32_max=-99;t42_max=-99;t43_max=-99;t32_min=-99;t42_min=-99;t43_min=-99; t1_max=-99;t2_max=-99;t3_max=-99;t4_max=-99;t1_mid=-99;t2_mid=-99;t3_mid=-99;t4_mid=-99;t1_min=-99;t2_min=-99;t3_min=-99;t4_min=-99;
        Mj_mid=-99;pt_mid=-99; t21_mid=-99; t31_mid=-99; t41_mid=-99;phi_mid=-99;eta_mid=-99;t32_mid=-99;t42_mid=-99;t43_mid=-99;
        dnnt_max=-99;dnnw_max=-99;dnnh_max=-99;dnndect_max=-99;dnndecw_max=-99;dnndech_max=-99;
        dnnt_mid=-99;dnnw_mid=-99;dnnh_mid=-99;dnndect_mid=-99;dnndecw_mid=-99;dnndech_mid=-99;
        dnnt_min=-99;dnnw_min=-99;dnnh_min=-99;dnndect_min=-99;dnndecw_min=-99;dnndech_min=-99;
        dnnz_max=-99;dnnzbb_max=-99;dnnhbb_max=-99;dnndecz_max=-99;dnndeczbb_max=-99;dnndechbb_max=-99;dnndecbb_max=-99;dnndecbbnog_max=-99;dnndeccc_max=-99;dnndecccnog_max=-99;
        dnnz_mid=-99;dnnzbb_mid=-99;dnnhbb_mid=-99;dnndecz_mid=-99;dnndeczbb_mid=-99;dnndechbb_mid=-99;dnndecbb_mid=-99;dnndecbbnog_mid=-99;dnndeccc_mid=-99;dnndecccnog_mid=-99;
        dnnz_min=-99;dnnzbb_min=-99;dnnhbb_min=-99;dnndecz_min=-99;dnndeczbb_min=-99;dnndechbb_min=-99;dnndecbb_min=-99;dnndecbbnog_min=-99;dnndeccc_min=-99;dnndecccnog_min=-99;
        dnndecrawqcd_max=-99;dnndecraww_max=-99;dnndecrawz_max=-99;dnndecrawzbb_max=-99;dnndecrawhbb_max=-99;dnndecrawt_max=-99;dnndecrawh4q_max=-99;
        dnndecrawqcd_mid=-99;dnndecraww_mid=-99;dnndecrawz_mid=-99;dnndecrawzbb_mid=-99;dnndecrawhbb_mid=-99;dnndecrawt_mid=-99;dnndecrawh4q_mid=-99;
        dnndecrawqcd_min=-99;dnndecraww_min=-99;dnndecrawz_min=-99;dnndecrawzbb_min=-99;dnndecrawhbb_min=-99;dnndecrawt_min=-99;dnndecrawh4q_min=-99;
        Double_t Mj_sort[3]={jet_mass_puppi_corr,jet_mass_puppi_corr_2,jet_mass_puppi_corr_3},
            pt_sort[3]={jet_pt_puppi,jet_pt_puppi_2,jet_pt_puppi_3},
            eta_sort[3]={jetAK8puppi_eta,jetAK8puppi_eta_2,jetAK8puppi_eta_3},
            phi_sort[3]={jetAK8puppi_phi,jetAK8puppi_phi_2,jetAK8puppi_phi_3},
        t1_sort[3]={jetAK8puppi_tau1,jetAK8puppi_tau1_2,jetAK8puppi_tau1_3},
        t2_sort[3]={jetAK8puppi_tau2,jetAK8puppi_tau2_2,jetAK8puppi_tau2_3},
        t3_sort[3]={jetAK8puppi_tau3,jetAK8puppi_tau3_2,jetAK8puppi_tau3_3},
        t4_sort[3]={jetAK8puppi_tau4,jetAK8puppi_tau4_2,jetAK8puppi_tau4_3},
            t21_sort[3]={jetAK8puppi_tau21,jetAK8puppi_tau21_2,jetAK8puppi_tau21_3},
            t31_sort[3]={jetAK8puppi_tau31,jetAK8puppi_tau31_2,jetAK8puppi_tau31_3},
            t41_sort[3]={jetAK8puppi_tau41,jetAK8puppi_tau41_2,jetAK8puppi_tau41_3},
            t32_sort[3]={jetAK8puppi_tau32,jetAK8puppi_tau32_2,jetAK8puppi_tau32_3},
            t42_sort[3]={jetAK8puppi_tau42,jetAK8puppi_tau42_2,jetAK8puppi_tau42_3},
            t43_sort[3]={jetAK8puppi_tau43,jetAK8puppi_tau43_2,jetAK8puppi_tau43_3},
            dnnt_sort[3]={jetAK8puppi_dnnTop,jetAK8puppi_dnnTop_2,jetAK8puppi_dnnTop_3},
            dnnw_sort[3]={jetAK8puppi_dnnW,jetAK8puppi_dnnW_2,jetAK8puppi_dnnW_3},
            dnnh_sort[3]={jetAK8puppi_dnnH4q,jetAK8puppi_dnnH4q_2,jetAK8puppi_dnnH4q_3},
            dnndect_sort[3]={jetAK8puppi_dnnDecorrTop,jetAK8puppi_dnnDecorrTop_2,jetAK8puppi_dnnDecorrTop_3},
            dnndecw_sort[3]={jetAK8puppi_dnnDecorrW,jetAK8puppi_dnnDecorrW_2,jetAK8puppi_dnnDecorrW_3},
            dnndech_sort[3]={jetAK8puppi_dnnDecorrH4q,jetAK8puppi_dnnDecorrH4q_2,jetAK8puppi_dnnDecorrH4q_3},
        dnnz_sort[3]={jetAK8puppi_dnnZ,jetAK8puppi_dnnZ_2,jetAK8puppi_dnnZ_3},
        dnnzbb_sort[3]={jetAK8puppi_dnnZbb,jetAK8puppi_dnnZbb_2,jetAK8puppi_dnnZbb_3},
        dnnhbb_sort[3]={jetAK8puppi_dnnHbb,jetAK8puppi_dnnHbb_2,jetAK8puppi_dnnHbb_3},
        dnndecz_sort[3]={jetAK8puppi_dnnDecorrZ,jetAK8puppi_dnnDecorrZ_2,jetAK8puppi_dnnDecorrZ_3},
        dnndeczbb_sort[3]={jetAK8puppi_dnnDecorrZbb,jetAK8puppi_dnnDecorrZbb_2,jetAK8puppi_dnnDecorrZbb_3},
        dnndechbb_sort[3]={jetAK8puppi_dnnDecorrHbb,jetAK8puppi_dnnDecorrHbb_2,jetAK8puppi_dnnDecorrHbb_3},
        dnndecbb_sort[3]={jetAK8puppi_dnnDecorrbb,jetAK8puppi_dnnDecorrbb_2,jetAK8puppi_dnnDecorrbb_3},
        dnndecbbnog_sort[3]={jetAK8puppi_dnnDecorrbbnog,jetAK8puppi_dnnDecorrbbnog_2,jetAK8puppi_dnnDecorrbbnog_3},
        dnndeccc_sort[3]={jetAK8puppi_dnnDecorrcc,jetAK8puppi_dnnDecorrcc_2,jetAK8puppi_dnnDecorrcc_3},
        dnndecccnog_sort[3]={jetAK8puppi_dnnDecorrccnog,jetAK8puppi_dnnDecorrccnog_2,jetAK8puppi_dnnDecorrccnog_3},
        dnndecrawqcd_sort[3]={jetAK8puppi_dnnDecorrqcd,jetAK8puppi_dnnDecorrqcd_2,jetAK8puppi_dnnDecorrqcd_3},
        dnndecrawt_sort[3]={jetAK8puppi_dnnDecorrtop,jetAK8puppi_dnnDecorrtop_2,jetAK8puppi_dnnDecorrtop_3},
        dnndecrawz_sort[3]={jetAK8puppi_dnnDecorrz,jetAK8puppi_dnnDecorrz_2,jetAK8puppi_dnnDecorrz_3},
        dnndecraww_sort[3]={jetAK8puppi_dnnDecorrw,jetAK8puppi_dnnDecorrw_2,jetAK8puppi_dnnDecorrw_3},
        dnndecrawzbb_sort[3]={jetAK8puppi_dnnDecorrzbb,jetAK8puppi_dnnDecorrzbb_2,jetAK8puppi_dnnDecorrzbb_3},
        dnndecrawhbb_sort[3]={jetAK8puppi_dnnDecorrhbb,jetAK8puppi_dnnDecorrhbb_2,jetAK8puppi_dnnDecorrhbb_3},
        dnndecrawh4q_sort[3]={jetAK8puppi_dnnDecorrh4q,jetAK8puppi_dnnDecorrh4q_2,jetAK8puppi_dnnDecorrh4q_3}
        ;
        
        if(Nj8==1){
            Mj_max=Mj_sort[0];pt_max=pt_sort[0];eta_max=eta_sort[0];phi_max=phi_sort[0];t21_max=t21_sort[0];t31_max=t31_sort[0];t41_max=t41_sort[0];t32_max=t32_sort[0];t42_max=t42_sort[0];t43_max=t43_sort[0];
            Mj_min=Mj_sort[0];pt_min=pt_sort[0];eta_min=eta_sort[0];phi_min=phi_sort[0];t21_min=t21_sort[0];t31_min=t31_sort[0];t41_min=t41_sort[0];t32_min=t32_sort[0];t42_min=t42_sort[0];t43_min=t43_sort[0];
            dnnt_max=dnnt_sort[0];dnnw_max=dnnw_sort[0];dnnh_max=dnnh_sort[0];dnndect_max=dnndect_sort[0];dnndecw_max=dnndecw_sort[0];dnndech_max=dnndech_sort[0];
            dnnt_min=dnnt_sort[0];dnnw_min=dnnw_sort[0];dnnh_min=dnnh_sort[0];dnndect_min=dnndect_sort[0];dnndecw_min=dnndecw_sort[0];dnndech_min=dnndech_sort[0];
            dnnz_max=dnnz_sort[0];dnnzbb_max=dnnzbb_sort[0];dnnhbb_max=dnnhbb_sort[0];dnndecz_max=dnndecz_sort[0];dnndeczbb_max=dnndeczbb_sort[0];dnndechbb_max=dnndechbb_sort[0];dnndecbb_max=dnndecbb_sort[0];dnndecbbnog_max=dnndecbbnog_sort[0];dnndeccc_max=dnndeccc_sort[0];dnndecccnog_max=dnndecccnog_sort[0];
            dnnz_min=dnnz_sort[0];dnnzbb_min=dnnzbb_sort[0];dnnhbb_min=dnnhbb_sort[0];dnndecz_min=dnndecz_sort[0];dnndeczbb_min=dnndeczbb_sort[0];dnndechbb_min=dnndechbb_sort[0];dnndecbb_min=dnndecbb_sort[0];dnndecbbnog_min=dnndecbbnog_sort[0];dnndeccc_min=dnndeccc_sort[0];dnndecccnog_min=dnndecccnog_sort[0];
            dnndecrawqcd_max=dnndecrawqcd_sort[0];dnndecrawt_max=dnndecrawt_sort[0];dnndecraww_max=dnndecraww_sort[0];dnndecrawz_max=dnndecrawz_sort[0];dnndecrawzbb_max=dnndecrawzbb_sort[0];dnndecrawhbb_max=dnndecrawhbb_sort[0];dnndecrawh4q_max=dnndecrawh4q_sort[0];
            dnndecrawqcd_min=dnndecrawqcd_sort[0];dnndecrawt_min=dnndecrawt_sort[0];dnndecraww_min=dnndecraww_sort[0];dnndecrawz_min=dnndecrawz_sort[0];dnndecrawzbb_min=dnndecrawzbb_sort[0];dnndecrawhbb_min=dnndecrawhbb_sort[0];dnndecrawh4q_min=dnndecrawh4q_sort[0];
            t1_max=t1_sort[0];t2_max=t2_sort[0];t3_max=t3_sort[0];t4_max=t4_sort[0];
            t1_min=t1_sort[0];t2_min=t2_sort[0];t3_min=t3_sort[0];t4_min=t4_sort[0];
        }
        if(Nj8==2){
            Int_t *indexmass=new Int_t[3];TMath::Sort(2,Mj_sort,indexmass,1);
            Mj_max=Mj_sort[indexmass[0]];pt_max=pt_sort[indexmass[0]];eta_max=eta_sort[indexmass[0]];phi_max=phi_sort[indexmass[0]];t21_max=t21_sort[indexmass[0]];t31_max=t31_sort[indexmass[0]];t41_max=t41_sort[indexmass[0]];t32_max=t32_sort[indexmass[0]];t42_max=t42_sort[indexmass[0]];t43_max=t43_sort[indexmass[0]];
            Mj_min=Mj_sort[indexmass[1]];pt_min=pt_sort[indexmass[1]];eta_min=eta_sort[indexmass[1]];phi_min=phi_sort[indexmass[1]];t21_min=t21_sort[indexmass[1]];t31_min=t31_sort[indexmass[1]];t41_min=t41_sort[indexmass[1]];t32_min=t32_sort[indexmass[1]];t42_min=t42_sort[indexmass[1]];t43_min=t43_sort[indexmass[1]];
            dnnt_max=dnnt_sort[indexmass[0]];dnnw_max=dnnw_sort[indexmass[0]];dnnh_max=dnnh_sort[indexmass[0]];dnndect_max=dnndect_sort[indexmass[0]];dnndecw_max=dnndecw_sort[indexmass[0]];dnndech_max=dnndech_sort[indexmass[0]];
            dnnt_min=dnnt_sort[indexmass[1]];dnnw_min=dnnw_sort[indexmass[1]];dnnh_min=dnnh_sort[indexmass[1]];dnndect_min=dnndect_sort[indexmass[1]];dnndecw_min=dnndecw_sort[indexmass[1]];dnndech_min=dnndech_sort[indexmass[1]];
            dnnz_max=dnnz_sort[indexmass[0]];dnnzbb_max=dnnzbb_sort[indexmass[0]];dnnhbb_max=dnnhbb_sort[indexmass[0]];dnndecz_max=dnndecz_sort[indexmass[0]];dnndeczbb_max=dnndeczbb_sort[indexmass[0]];dnndechbb_max=dnndechbb_sort[indexmass[0]];dnndecbb_max=dnndecbb_sort[indexmass[0]];dnndecbbnog_max=dnndecbbnog_sort[indexmass[0]];dnndeccc_max=dnndeccc_sort[indexmass[0]];dnndecccnog_max=dnndecccnog_sort[indexmass[0]];
            dnnz_min=dnnz_sort[indexmass[1]];dnnzbb_min=dnnzbb_sort[indexmass[1]];dnnhbb_min=dnnhbb_sort[indexmass[1]];dnndecz_min=dnndecz_sort[indexmass[1]];dnndeczbb_min=dnndeczbb_sort[indexmass[1]];dnndechbb_min=dnndechbb_sort[indexmass[1]];dnndecbb_min=dnndecbb_sort[indexmass[1]];dnndecbbnog_min=dnndecbbnog_sort[indexmass[1]];dnndeccc_min=dnndeccc_sort[indexmass[1]];dnndecccnog_min=dnndecccnog_sort[indexmass[1]];
            dnndecrawqcd_max=dnndecrawqcd_sort[indexmass[0]];dnndecrawt_max=dnndecrawt_sort[indexmass[0]];dnndecraww_max=dnndecraww_sort[indexmass[0]];dnndecrawz_max=dnndecrawz_sort[indexmass[0]];dnndecrawzbb_max=dnndecrawzbb_sort[indexmass[0]];dnndecrawhbb_max=dnndecrawhbb_sort[indexmass[0]];dnndecrawh4q_max=dnndecrawh4q_sort[indexmass[0]];
            dnndecrawqcd_min=dnndecrawqcd_sort[indexmass[1]];dnndecrawt_min=dnndecrawt_sort[indexmass[1]];dnndecraww_min=dnndecraww_sort[indexmass[1]];dnndecrawz_min=dnndecrawz_sort[indexmass[1]];dnndecrawzbb_min=dnndecrawzbb_sort[indexmass[1]];dnndecrawhbb_min=dnndecrawhbb_sort[indexmass[1]];dnndecrawh4q_min=dnndecrawh4q_sort[indexmass[1]];
            t1_max=t1_sort[indexmass[0]];t2_max=t2_sort[indexmass[0]];t3_max=t3_sort[indexmass[0]];t4_max=t4_sort[indexmass[0]];
            t1_min=t1_sort[indexmass[1]];t2_min=t2_sort[indexmass[1]];t3_min=t3_sort[indexmass[1]];t4_min=t4_sort[indexmass[1]];
        }
        if(Nj8>=3){
            Int_t *indexmass=new Int_t[3];TMath::Sort(3,Mj_sort,indexmass,1);
            Mj_max=Mj_sort[indexmass[0]];pt_max=pt_sort[indexmass[0]];eta_max=eta_sort[indexmass[0]];phi_max=phi_sort[indexmass[0]];t21_max=t21_sort[indexmass[0]];t31_max=t31_sort[indexmass[0]];t41_max=t41_sort[indexmass[0]];t32_max=t32_sort[indexmass[0]];t42_max=t42_sort[indexmass[0]];t43_max=t43_sort[indexmass[0]];
            Mj_mid=Mj_sort[indexmass[1]];pt_mid=pt_sort[indexmass[1]];eta_mid=eta_sort[indexmass[1]];phi_mid=phi_sort[indexmass[1]];t21_mid=t21_sort[indexmass[1]];t31_mid=t31_sort[indexmass[1]];t41_mid=t41_sort[indexmass[1]];t32_mid=t32_sort[indexmass[1]];t42_mid=t42_sort[indexmass[1]];t43_mid=t43_sort[indexmass[1]];
            Mj_min=Mj_sort[indexmass[2]];pt_min=pt_sort[indexmass[2]];eta_min=eta_sort[indexmass[2]];phi_min=phi_sort[indexmass[2]];t21_min=t21_sort[indexmass[2]];t31_min=t31_sort[indexmass[2]];t41_min=t41_sort[indexmass[2]];t32_min=t32_sort[indexmass[2]];t42_min=t42_sort[indexmass[2]];t43_min=t43_sort[indexmass[2]];
            dnnt_max=dnnt_sort[indexmass[0]];dnnw_max=dnnw_sort[indexmass[0]];dnnh_max=dnnh_sort[indexmass[0]];dnndect_max=dnndect_sort[indexmass[0]];dnndecw_max=dnndecw_sort[indexmass[0]];dnndech_max=dnndech_sort[indexmass[0]];
            dnnt_mid=dnnt_sort[indexmass[1]];dnnw_mid=dnnw_sort[indexmass[1]];dnnh_mid=dnnh_sort[indexmass[1]];dnndect_mid=dnndect_sort[indexmass[1]];dnndecw_mid=dnndecw_sort[indexmass[1]];dnndech_mid=dnndech_sort[indexmass[1]];
            dnnt_min=dnnt_sort[indexmass[2]];dnnw_min=dnnw_sort[indexmass[2]];dnnh_min=dnnh_sort[indexmass[2]];dnndect_min=dnndect_sort[indexmass[2]];dnndecw_min=dnndecw_sort[indexmass[2]];dnndech_min=dnndech_sort[indexmass[2]];
            dnnz_max=dnnz_sort[indexmass[0]];dnnzbb_max=dnnzbb_sort[indexmass[0]];dnnhbb_max=dnnhbb_sort[indexmass[0]];dnndecz_max=dnndecz_sort[indexmass[0]];dnndeczbb_max=dnndeczbb_sort[indexmass[0]];dnndechbb_max=dnndechbb_sort[indexmass[0]];dnndecbb_max=dnndecbb_sort[indexmass[0]];dnndecbbnog_max=dnndecbbnog_sort[indexmass[0]];dnndeccc_max=dnndeccc_sort[indexmass[0]];dnndecccnog_max=dnndecccnog_sort[indexmass[0]];
            dnnz_mid=dnnz_sort[indexmass[1]];dnnzbb_mid=dnnzbb_sort[indexmass[1]];dnnhbb_mid=dnnhbb_sort[indexmass[1]];dnndecz_mid=dnndecz_sort[indexmass[1]];dnndeczbb_mid=dnndeczbb_sort[indexmass[1]];dnndechbb_mid=dnndechbb_sort[indexmass[1]];dnndecbb_mid=dnndecbb_sort[indexmass[1]];dnndecbbnog_mid=dnndecbbnog_sort[indexmass[1]];dnndeccc_mid=dnndeccc_sort[indexmass[1]];dnndecccnog_mid=dnndecccnog_sort[indexmass[1]];
            dnnz_min=dnnz_sort[indexmass[2]];dnnzbb_min=dnnzbb_sort[indexmass[2]];dnnhbb_min=dnnhbb_sort[indexmass[2]];dnndecz_min=dnndecz_sort[indexmass[2]];dnndeczbb_min=dnndeczbb_sort[indexmass[2]];dnndechbb_min=dnndechbb_sort[indexmass[2]];dnndecbb_min=dnndecbb_sort[indexmass[2]];dnndecbbnog_min=dnndecbbnog_sort[indexmass[2]];dnndeccc_min=dnndeccc_sort[indexmass[2]];dnndecccnog_min=dnndecccnog_sort[indexmass[2]];
            dnndecrawqcd_max=dnndecrawqcd_sort[indexmass[0]];dnndecrawt_max=dnndecrawt_sort[indexmass[0]];dnndecraww_max=dnndecraww_sort[indexmass[0]];dnndecrawz_max=dnndecrawz_sort[indexmass[0]];dnndecrawzbb_max=dnndecrawzbb_sort[indexmass[0]];dnndecrawhbb_max=dnndecrawhbb_sort[indexmass[0]];dnndecrawh4q_max=dnndecrawh4q_sort[indexmass[0]];
            dnndecrawqcd_mid=dnndecrawqcd_sort[indexmass[1]];dnndecrawt_mid=dnndecrawt_sort[indexmass[1]];dnndecraww_mid=dnndecraww_sort[indexmass[1]];dnndecrawz_mid=dnndecrawz_sort[indexmass[1]];dnndecrawzbb_mid=dnndecrawzbb_sort[indexmass[1]];dnndecrawhbb_mid=dnndecrawhbb_sort[indexmass[1]];dnndecrawh4q_mid=dnndecrawh4q_sort[indexmass[1]];
            dnndecrawqcd_min=dnndecrawqcd_sort[indexmass[2]];dnndecrawt_min=dnndecrawt_sort[indexmass[2]];dnndecraww_min=dnndecraww_sort[indexmass[2]];dnndecrawz_min=dnndecrawz_sort[indexmass[2]];dnndecrawzbb_min=dnndecrawzbb_sort[indexmass[2]];dnndecrawhbb_min=dnndecrawhbb_sort[indexmass[2]];dnndecrawh4q_min=dnndecrawh4q_sort[indexmass[2]];
            t1_max=t1_sort[indexmass[0]];t2_max=t2_sort[indexmass[0]];t3_max=t3_sort[indexmass[0]];t4_max=t4_sort[indexmass[0]];
            t1_mid=t1_sort[indexmass[1]];t2_mid=t2_sort[indexmass[1]];t3_mid=t3_sort[indexmass[1]];t4_mid=t4_sort[indexmass[1]];
            t1_min=t1_sort[indexmass[2]];t2_min=t2_sort[indexmass[2]];t3_min=t3_sort[indexmass[2]];t4_min=t4_sort[indexmass[2]];
        }

}

#endif