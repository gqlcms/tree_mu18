    fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("nLooseEle", &nLooseEle, &b_nLooseEle);
   fChain->SetBranchAddress("nLooseMu", &nLooseMu, &b_nLooseMu);
    fChain->SetBranchAddress("L1prefiring", &L1prefiring, &b_L1prefiring);
    fChain->SetBranchAddress("L1prefiringup", &L1prefiringup, &b_L1prefiringup);
    fChain->SetBranchAddress("L1prefiringdown", &L1prefiringdown, &b_L1prefiringdown);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("ptVlep", &ptVlep, &b_ptVlep);
   fChain->SetBranchAddress("yVlep", &yVlep, &b_yVlep);
   fChain->SetBranchAddress("phiVlep", &phiVlep, &b_phiVlep);
   fChain->SetBranchAddress("massVlep", &massVlep, &b_massVlep);
   fChain->SetBranchAddress("mtVlep", &mtVlep, &b_mtVlep);


   fChain->SetBranchAddress("lep", &lep, &b_lep);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("ptlep1", &ptlep1, &b_ptlep1);
   fChain->SetBranchAddress("ptlep2", &ptlep2, &b_ptlep2);
   fChain->SetBranchAddress("etalep1", &etalep1, &b_etalep1);
   fChain->SetBranchAddress("etalep2", &etalep2, &b_etalep2);
   fChain->SetBranchAddress("philep1", &philep1, &b_philep1);
   fChain->SetBranchAddress("philep2", &philep2, &b_philep2);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("theWeight", &theWeight, &b_theWeight);
   fChain->SetBranchAddress("nump", &nump, &b_nump);
   fChain->SetBranchAddress("numm", &numm, &b_numm);
   fChain->SetBranchAddress("npT", &npT, &b_npT);
   fChain->SetBranchAddress("npIT", &npIT, &b_npIT);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);




   fChain->SetBranchAddress("vbftag", &vbftag, &b_vbftag);
   fChain->SetBranchAddress("IDLoose", &IDLoose, &b_IDLoose);
   fChain->SetBranchAddress("IDTight", &IDTight, &b_IDTight);
    fChain->SetBranchAddress("IDLoose_2", &IDLoose_2, &b_IDLoose_2);
    fChain->SetBranchAddress("IDTight_2", &IDTight_2, &b_IDTight_2);
    fChain->SetBranchAddress("IDLoose_3", &IDLoose_3, &b_IDLoose_3);
    fChain->SetBranchAddress("IDTight_3", &IDTight_3, &b_IDTight_3);

         fChain->SetBranchAddress("trackIso", &trackIso, &b_trackIso);
   fChain->SetBranchAddress("muchaiso", &muchaiso, &b_muchaiso);
   fChain->SetBranchAddress("muneuiso", &muneuiso, &b_muneuiso);
   fChain->SetBranchAddress("muphoiso", &muphoiso, &b_muphoiso);
   fChain->SetBranchAddress("muPU", &muPU, &b_muPU);
   fChain->SetBranchAddress("muisolation", &muisolation, &b_muisolation);
   fChain->SetBranchAddress("METraw_et", &METraw_et, &b_METraw_et);
   fChain->SetBranchAddress("METraw_phi", &METraw_phi, &b_METraw_phi);
   fChain->SetBranchAddress("METraw_sumEt", &METraw_sumEt, &b_METraw_sumEt);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("jetAK8puppi_pt1", jetAK8puppi_pt1, &b_jetAK8puppi_pt1);
   fChain->SetBranchAddress("jetAK8puppi_eta1", jetAK8puppi_eta1, &b_jetAK8puppi_eta1);
   fChain->SetBranchAddress("jetAK8puppi_mass1", jetAK8puppi_mass1, &b_jetAK8puppi_mass1);
   fChain->SetBranchAddress("ptVlepJEC", &ptVlepJEC, &b_ptVlepJEC);
   fChain->SetBranchAddress("yVlepJEC", &yVlepJEC, &b_yVlepJEC);
   fChain->SetBranchAddress("phiVlepJEC", &phiVlepJEC, &b_phiVlepJEC);
   fChain->SetBranchAddress("massVlepJEC", &massVlepJEC, &b_massVlepJEC);

   fChain->SetBranchAddress("jetAK8puppi_sd", &jetAK8puppi_sd, &b_jetAK8puppi_sd);
         fChain->SetBranchAddress("jetAK8puppi_ptJEC", &jetAK8puppi_ptJEC, &b_jetAK8puppi_ptJEC);
   fChain->SetBranchAddress("jetAK8puppi_eta", &jetAK8puppi_eta, &b_jetAK8puppi_eta);
   fChain->SetBranchAddress("jetAK8puppi_phi", &jetAK8puppi_phi, &b_jetAK8puppi_phi);


    fChain->SetBranchAddress("m_jlv", &m_jlv, &b_m_jlv);
   fChain->SetBranchAddress("candMasspuppiJEC", &candMasspuppiJEC, &b_candMasspuppiJEC);

   fChain->SetBranchAddress("jetAK8puppi_sd_2", &jetAK8puppi_sd_2, &b_jetAK8puppi_sd_2);
         fChain->SetBranchAddress("jetAK8puppi_ptJEC_2", &jetAK8puppi_ptJEC_2, &b_jetAK8puppi_ptJEC_2);
   fChain->SetBranchAddress("jetAK8puppi_eta_2", &jetAK8puppi_eta_2, &b_jetAK8puppi_eta_2);
   fChain->SetBranchAddress("jetAK8puppi_phi_2", &jetAK8puppi_phi_2, &b_jetAK8puppi_phi_2);

    

    fChain->SetBranchAddress("jetAK8puppi_sd_3", &jetAK8puppi_sd_3, &b_jetAK8puppi_sd_3);
          fChain->SetBranchAddress("jetAK8puppi_ptJEC_3", &jetAK8puppi_ptJEC_3, &b_jetAK8puppi_ptJEC_3);
    fChain->SetBranchAddress("jetAK8puppi_eta_3", &jetAK8puppi_eta_3, &b_jetAK8puppi_eta_3);
    fChain->SetBranchAddress("jetAK8puppi_phi_3", &jetAK8puppi_phi_3, &b_jetAK8puppi_phi_3);


fChain->SetBranchAddress("massww",&massww);


   fChain->SetBranchAddress("mtVlepJEC", &mtVlepJEC, &b_mtVlepJEC);
   fChain->SetBranchAddress("HLT_Ele1", &HLT_Ele1, &b_HLT_Ele1);
   fChain->SetBranchAddress("HLT_Ele2", &HLT_Ele2, &b_HLT_Ele2);
   fChain->SetBranchAddress("HLT_Ele3", &HLT_Ele3, &b_HLT_Ele3);
   fChain->SetBranchAddress("HLT_Ele4", &HLT_Ele4, &b_HLT_Ele4);
   fChain->SetBranchAddress("HLT_Ele5", &HLT_Ele5, &b_HLT_Ele5);
   fChain->SetBranchAddress("HLT_Ele6", &HLT_Ele6, &b_HLT_Ele6);
   fChain->SetBranchAddress("HLT_Ele7", &HLT_Ele7, &b_HLT_Ele7);
   fChain->SetBranchAddress("HLT_Ele8", &HLT_Ele8, &b_HLT_Ele8);
   fChain->SetBranchAddress("HLT_Mu1", &HLT_Mu1, &b_HLT_Mu1);
   fChain->SetBranchAddress("HLT_Mu2", &HLT_Mu2, &b_HLT_Mu2);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu4", &HLT_Mu4, &b_HLT_Mu4);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu6", &HLT_Mu6, &b_HLT_Mu6);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
   fChain->SetBranchAddress("HLT_Mu10", &HLT_Mu10, &b_HLT_Mu10);
   fChain->SetBranchAddress("HLT_Mu11", &HLT_Mu11, &b_HLT_Mu11);
   fChain->SetBranchAddress("HLT_Mu12", &HLT_Mu12, &b_HLT_Mu12);
   
   fChain->SetBranchAddress("passFilter_HBHE", &passFilter_HBHE, &b_passFilter_HBHE_);
   fChain->SetBranchAddress("passFilter_HBHEIso", &passFilter_HBHEIso, &b_passFilter_HBHEIso_);
   fChain->SetBranchAddress("passFilter_GlobalHalo", &passFilter_GlobalHalo, &b_passFilter_GlobalHalo_);
   fChain->SetBranchAddress("passFilter_ECALDeadCell", &passFilter_ECALDeadCell, &b_passFilter_ECALDeadCell_);
   fChain->SetBranchAddress("passFilter_GoodVtx", &passFilter_GoodVtx, &b_passFilter_GoodVtx_);
   fChain->SetBranchAddress("passFilter_EEBadSc", &passFilter_EEBadSc, &b_passFilter_EEBadSc_);
   fChain->SetBranchAddress("passFilter_badMuon", &passFilter_badMuon, &b_passFilter_badMuon_);
   fChain->SetBranchAddress("passFilter_badChargedHadron", &passFilter_badChargedHadron, &b_passFilter_badChargedHadron_);
    fChain->SetBranchAddress("passecalBadCalibFilterUpdate", &passecalBadCalibFilterUpdate, &b_passecalBadCalibFilterUpdate_);
fChain->SetBranchAddress("ak4jet_hf", ak4jet_hf, &b_ak4jet_hf);
   fChain->SetBranchAddress("ak4jet_pf", ak4jet_pf, &b_ak4jet_pf);
   fChain->SetBranchAddress("ak4jet_pt", ak4jet_pt, &b_ak4jet_pt);
   fChain->SetBranchAddress("ak4jet_pt_JEC_up", ak4jet_pt_JEC_up, &b_ak4jet_pt_JEC_up);
   fChain->SetBranchAddress("ak4jet_pt_JEC_down", ak4jet_pt_JEC_down, &b_ak4jet_pt_JEC_down);
   fChain->SetBranchAddress("ak4jet_pt_JER_up", ak4jet_pt_JER_up, &b_ak4jet_pt_JER_up);
   fChain->SetBranchAddress("ak4jet_pt_JER_down", ak4jet_pt_JER_down, &b_ak4jet_pt_JER_down);
   fChain->SetBranchAddress("ak4jet_pt_uncorr", ak4jet_pt_uncorr, &b_ak4jet_pt_uncorr);
   fChain->SetBranchAddress("ak4jet_eta", ak4jet_eta, &b_ak4jet_eta);
   fChain->SetBranchAddress("ak4jet_phi", ak4jet_phi, &b_ak4jet_phi);
   fChain->SetBranchAddress("ak4jet_e", ak4jet_e, &b_ak4jet_e);
   fChain->SetBranchAddress("ak4jet_dr", ak4jet_dr, &b_ak4jet_dr);
   fChain->SetBranchAddress("ak4jet_csv", ak4jet_csv, &b_ak4jet_csv);
   fChain->SetBranchAddress("ak4jet_icsv", ak4jet_icsv, &b_ak4jet_icsv);
    fChain->SetBranchAddress("ak4jet_deepcsvb", ak4jet_deepcsvb, &b_ak4jet_deepcsvb);
    fChain->SetBranchAddress("ak4jet_deepcsvbb", ak4jet_deepcsvbb, &b_ak4jet_deepcsvbb);
   fChain->SetBranchAddress("ak4jet_IDLoose", ak4jet_IDLoose, &b_ak4jet_IDLoose);
   fChain->SetBranchAddress("ak4jet_IDTight", ak4jet_IDTight, &b_ak4jet_IDTight);
    
                                            fChain->SetBranchAddress("jetAK8puppi_dnnTop", &jetAK8puppi_dnnTop, &b_jetAK8puppi_dnnTop);
    fChain->SetBranchAddress("jetAK8puppi_dnnW", &jetAK8puppi_dnnW  , &b_jetAK8puppi_dnnW  );
    fChain->SetBranchAddress("jetAK8puppi_dnnH4q", &jetAK8puppi_dnnH4q, &b_jetAK8puppi_dnnH4q);
    fChain->SetBranchAddress("jetAK8puppi_dnnTop_2", &jetAK8puppi_dnnTop_2, &b_jetAK8puppi_dnnTop_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnW_2", &jetAK8puppi_dnnW_2  , &b_jetAK8puppi_dnnW_2  );
    fChain->SetBranchAddress("jetAK8puppi_dnnH4q_2", &jetAK8puppi_dnnH4q_2, &b_jetAK8puppi_dnnH4q_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnTop_3", &jetAK8puppi_dnnTop_3, &b_jetAK8puppi_dnnTop_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnW_3", &jetAK8puppi_dnnW_3  , &b_jetAK8puppi_dnnW_3  );
    fChain->SetBranchAddress("jetAK8puppi_dnnH4q_3", &jetAK8puppi_dnnH4q_3, &b_jetAK8puppi_dnnH4q_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrTop", &jetAK8puppi_dnnDecorrTop, &b_jetAK8puppi_dnnDecorrTop);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrW", &jetAK8puppi_dnnDecorrW  , &b_jetAK8puppi_dnnDecorrW  );
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrH4q", &jetAK8puppi_dnnDecorrH4q, &b_jetAK8puppi_dnnDecorrH4q);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrTop_2", &jetAK8puppi_dnnDecorrTop_2, &b_jetAK8puppi_dnnDecorrTop_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrW_2", &jetAK8puppi_dnnDecorrW_2  , &b_jetAK8puppi_dnnDecorrW_2  );
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrH4q_2", &jetAK8puppi_dnnDecorrH4q_2, &b_jetAK8puppi_dnnDecorrH4q_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrTop_3", &jetAK8puppi_dnnDecorrTop_3, &b_jetAK8puppi_dnnDecorrTop_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrW_3", &jetAK8puppi_dnnDecorrW_3  , &b_jetAK8puppi_dnnDecorrW_3  );
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrH4q_3", &jetAK8puppi_dnnDecorrH4q_3, &b_jetAK8puppi_dnnDecorrH4q_3);
    
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbb", &jetAK8puppi_dnnDecorrbb, &b_jetAK8puppi_dnnDecorrbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbbnog", &jetAK8puppi_dnnDecorrbbnog, &b_jetAK8puppi_dnnDecorrbbnog);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrcc", &jetAK8puppi_dnnDecorrcc, &b_jetAK8puppi_dnnDecorrcc);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrccnog", &jetAK8puppi_dnnDecorrccnog, &b_jetAK8puppi_dnnDecorrccnog);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbb_2", &jetAK8puppi_dnnDecorrbb_2, &b_jetAK8puppi_dnnDecorrbb_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbbnog_2", &jetAK8puppi_dnnDecorrbbnog_2, &b_jetAK8puppi_dnnDecorrbbnog_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrcc_2", &jetAK8puppi_dnnDecorrcc_2, &b_jetAK8puppi_dnnDecorrcc_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrccnog_2", &jetAK8puppi_dnnDecorrccnog_2, &b_jetAK8puppi_dnnDecorrccnog_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbb_3", &jetAK8puppi_dnnDecorrbb_3, &b_jetAK8puppi_dnnDecorrbb_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbbnog_3", &jetAK8puppi_dnnDecorrbbnog_3, &b_jetAK8puppi_dnnDecorrbbnog_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrcc_3", &jetAK8puppi_dnnDecorrcc_3, &b_jetAK8puppi_dnnDecorrcc_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrccnog_3", &jetAK8puppi_dnnDecorrccnog_3, &b_jetAK8puppi_dnnDecorrccnog_3);

    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrqcd", &jetAK8puppi_dnnDecorrqcd, &b_jetAK8puppi_dnnDecorrqcd);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrtop", &jetAK8puppi_dnnDecorrtop, &b_jetAK8puppi_dnnDecorrtop);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrz", &jetAK8puppi_dnnDecorrz, &b_jetAK8puppi_dnnDecorrz);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrw", &jetAK8puppi_dnnDecorrw, &b_jetAK8puppi_dnnDecorrw);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrzbb", &jetAK8puppi_dnnDecorrzbb, &b_jetAK8puppi_dnnDecorrzbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrhbb", &jetAK8puppi_dnnDecorrhbb, &b_jetAK8puppi_dnnDecorrhbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrh4q", &jetAK8puppi_dnnDecorrh4q, &b_jetAK8puppi_dnnDecorrh4q);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrqcd_2", &jetAK8puppi_dnnDecorrqcd_2, &b_jetAK8puppi_dnnDecorrqcd_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrtop_2", &jetAK8puppi_dnnDecorrtop_2, &b_jetAK8puppi_dnnDecorrtop_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrz_2", &jetAK8puppi_dnnDecorrz_2, &b_jetAK8puppi_dnnDecorrz_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrw_2", &jetAK8puppi_dnnDecorrw_2, &b_jetAK8puppi_dnnDecorrw_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrzbb_2", &jetAK8puppi_dnnDecorrzbb_2, &b_jetAK8puppi_dnnDecorrzbb_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrhbb_2", &jetAK8puppi_dnnDecorrhbb_2, &b_jetAK8puppi_dnnDecorrhbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrh4q_2", &jetAK8puppi_dnnDecorrh4q_2, &b_jetAK8puppi_dnnDecorrh4q_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrqcd_3", &jetAK8puppi_dnnDecorrqcd_3, &b_jetAK8puppi_dnnDecorrqcd_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrtop_3", &jetAK8puppi_dnnDecorrtop_3, &b_jetAK8puppi_dnnDecorrtop_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrz_3", &jetAK8puppi_dnnDecorrz_3, &b_jetAK8puppi_dnnDecorrz_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrw_3", &jetAK8puppi_dnnDecorrw_3, &b_jetAK8puppi_dnnDecorrw_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrzbb_3", &jetAK8puppi_dnnDecorrzbb_3, &b_jetAK8puppi_dnnDecorrzbb_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrhbb_3", &jetAK8puppi_dnnDecorrhbb_3, &b_jetAK8puppi_dnnDecorrhbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrh4q_3", &jetAK8puppi_dnnDecorrh4q_3, &b_jetAK8puppi_dnnDecorrh4q_3);