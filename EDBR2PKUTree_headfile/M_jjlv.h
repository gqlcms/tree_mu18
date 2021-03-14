#ifndef _Mjjlv_
#define _Mjjlv_

void EDBR2PKUTree::Mjjlv() {
        m_JJlv             = Double_t(candMasspuppiJEC);
        m_Jlv_new=-99;
        m_JJlv_new=-99;
        m_Jlv_JEC_up=-99;
        m_Jlv_JEC_down=-99;
        m_Jlv_JER_up=-99;
        m_Jlv_JER_down=-99;
        m_JJlv_JEC_up=-99;
        m_JJlv_JEC_down=-99;
        m_JJlv_JER_up=-99;
        m_JJlv_JER_down=-99;
        
        if (candMasspuppiJEC_new>0) m_JJlv_new             = Double_t(candMasspuppiJEC_new);
        if (candMasspuppiJEC_JEC_up>0) m_JJlv_JEC_up             = Double_t(candMasspuppiJEC_JEC_up);
        if (candMasspuppiJEC_JEC_down>0) m_JJlv_JEC_down             = Double_t(candMasspuppiJEC_JEC_down);
        if (candMasspuppiJEC_JER_up>0) m_JJlv_JER_up             = Double_t(candMasspuppiJEC_JER_up);
        if (candMasspuppiJEC_JER_down>0) m_JJlv_JER_down             = Double_t(candMasspuppiJEC_JER_down);
        
        
        m_Jlv = m_jlv;
        if (m_jlv_new>0) m_Jlv_new = m_jlv_new;
        if (m_jlv_JEC_up>0) m_Jlv_JEC_up = m_jlv_JEC_up;
        if (m_jlv_JEC_down>0) m_Jlv_JEC_down = m_jlv_JEC_down;
        if (m_jlv_JER_up>0) m_Jlv_JER_up = m_jlv_JER_up;
        if (m_jlv_JER_down>0) m_Jlv_JER_down = m_jlv_JER_down;
        
        
}

void EDBR2PKUTree::inv_mass() {

        m_Jlvj=-99;m_JJlvj=-99;m_Jlvjj=-99;m_JJlvjj=-99;m_JJJlv=-99;m_JJJlvj=-99;m_JJJlvjj=-99;
        TLorentzVector ak8jet1_3,ak8jet2_3,ak8jet3_3,ak4exjet1,ak4exjet2,Wlv;
        ak8jet1_3.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
        ak8jet2_3.SetPtEtaPhiM(jet_pt_puppi_2,jetAK8puppi_eta_2,jetAK8puppi_phi_2,jet_mass_puppi_corr_2);
        ak8jet3_3.SetPtEtaPhiM(jet_pt_puppi_3,jetAK8puppi_eta_3,jetAK8puppi_phi_3,jet_mass_puppi_corr_3);
        Wlv.SetPtEtaPhiM(W_pt,W_eta,W_phi,W_Mass);
        // Wlv.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
        m_JJJlv=(ak8jet1_3+ak8jet2_3+ak8jet3_3+Wlv).Mag();
        m_JJlv=(ak8jet1_3+ak8jet2_3+Wlv).Mag();
        if(Nj4ex_dR0p8==0){
            m_Jlvj=m_jlv;m_Jlvjj=m_jlv;
            if(Nj8>=2)      {m_JJlvj=m_JJlv;m_JJlvjj=m_JJlv;}
            if(Nj8>=3)      {m_JJJlvj=m_JJJlv;m_JJJlvjj=m_JJJlv;}
        }
        if(Nj4ex_dR0p8==1){
            ak4exjet1.SetPtEtaPhiE(ak4expt1_dR0p8,ak4exeta1_dR0p8,ak4exphi1_dR0p8,ak4exe1_dR0p8);
            m_Jlvj = (ak8jet1_3+Wlv+ak4exjet1).Mag(); m_Jlvjj=m_Jlvj;
            if(Nj8==2)  {m_JJlvj=(ak8jet1_3+ak8jet2_3+Wlv+ak4exjet1).Mag();m_JJlvjj=m_JJlvj;}
            if(Nj8>=3)  {m_JJJlvj=(ak8jet1_3+ak8jet2_3+ak8jet3_3+Wlv+ak4exjet1).Mag();m_JJJlvjj=m_JJlvj;}
        }
        if(Nj4ex_dR0p8>=2){
            ak4exjet1.SetPtEtaPhiE(ak4expt1_dR0p8,ak4exeta1_dR0p8,ak4exphi1_dR0p8,ak4exe1_dR0p8);
            ak4exjet2.SetPtEtaPhiE(ak4expt2_dR0p8,ak4exeta2_dR0p8,ak4exphi2_dR0p8,ak4exe2_dR0p8);
            m_Jlvj = (ak8jet1_3+Wlv+ak4exjet1).Mag(); m_Jlvjj=(ak8jet1_3+Wlv+ak4exjet1+ak4exjet2).Mag();
            if(Nj8==2)  {m_JJlvj=(ak8jet1_3+ak8jet2_3+Wlv+ak4exjet1).Mag();m_JJlvjj=(ak8jet1_3+ak8jet2_3+Wlv+ak4exjet1+ak4exjet2).Mag();}
            if(Nj8>=3)  {m_JJJlvj=(ak8jet1_3+ak8jet2_3+ak8jet3_3+Wlv+ak4exjet1).Mag();m_JJJlvjj=(ak8jet1_3+ak8jet2_3+ak8jet3_3+Wlv+ak4exjet1+ak4exjet2).Mag();}
        }

}


#endif