#ifndef ForAK4_switchtoarray_h
#define ForAK4_switchtoarray_h

void EDBR2PKUTree::ForAK4_switchtoarray() {
    int number_of_jets = 8;
    int count = 0;
    std::vector<float> test = *ak4jet_pt_v;
    for(int ijet = 0;ijet<number_of_jets;ijet++){
        ak4jet_pt[ijet] = -99.;
        ak4jet_eta[ijet] = -99.;
        ak4jet_phi[ijet] = -99.;
        ak4jet_e[ijet] = -99.;
        ak4jet_deepcsvb_bb[ijet] = -99.;
        // ak4jet_IDLoose[ijet] = -99.;
    }

    for(auto pt = test.begin(); pt != test.end(); pt++)
    {
        if (count<number_of_jets){
            // cout<<count<<endl;

            ak4jet_pt[count] = *pt;
            count++;
            // cout<<count<<endl;
            // cout<<ak4jet_pt_v->size()<<endl;
        }
    }
    
    // cout<<"first"<<endl;

    count = 0;
    for(auto pt = ak4jet_eta_v->begin(); pt != ak4jet_eta_v->end(); pt++)
    {
        if (count<number_of_jets){
            ak4jet_eta[count] = *pt;
            count++;
            // cout<<count<<endl;
        }
    }

    count = 0;
    for(auto pt = ak4jet_phi_v->begin(); pt != ak4jet_phi_v->end(); pt++)
    {
        if (count<number_of_jets){
            ak4jet_phi[count] = *pt;
            count++;
        }
    }

    count = 0;
    for(auto pt = ak4jet_e_v->begin(); pt != ak4jet_e_v->end(); pt++)
    {
        if (count<number_of_jets){
            ak4jet_e[count] = *pt;
            count++;
        }
    }

     count = 0;
    for(auto pt = Jet_btagDeepB->begin(); pt != Jet_btagDeepB->end(); pt++)
    {
        if (count<number_of_jets){
            ak4jet_deepcsvb_bb[count] = *pt;
            ak4jet_IDLoose[count] = 1;
        

            count++;
            // cout<<count<<endl;
        }
    }

     count = 0;
    for(auto pt = Common_LHEWeight_mg_reweighting->begin(); pt != Common_LHEWeight_mg_reweighting->end(); pt++)
    {
        if (count<10){
            pweight[count] = *pt; // the Common_LHEWeight_mg_reweighting is (madgraph reweight) / original weight
            // cout <<*pt<<endl;
        }
count++;
    }


    jetAK8puppi_ptJEC = jetAK8puppi_ptJEC_f;
    jetAK8puppi_eta = jetAK8puppi_eta_f;
    jetAK8puppi_phi = jetAK8puppi_phi_f;
    jetAK8puppi_sd = jetAK8puppi_sd_f;
    jetAK8puppi_ptJEC_2 = jetAK8puppi_ptJEC_2_f;
    jetAK8puppi_eta_2 = jetAK8puppi_eta_2_f;
    jetAK8puppi_phi_2 = jetAK8puppi_phi_2_f;
    jetAK8puppi_sd_2 = jetAK8puppi_sd_2_f;
    jetAK8puppi_dnnDecorrW_2 = jetAK8puppi_dnnDecorrW_2_f;
    jetAK8puppi_dnnDecorrW = jetAK8puppi_dnnDecorrW_f;

    ptlep1 =         ptlep1_f;
    ptlep2 =         ptlep2_f;
    etalep1 =         etalep1_f;
    etalep2 =         etalep2_f;
    philep1 =       philep1_f;
    philep2 =         philep2_f;

    ptVlepJEC =         ptVlepJEC_f;
    yVlepJEC =        yVlepJEC_f;
    phiVlepJEC =         phiVlepJEC_f;
    massVlepJEC =        massVlepJEC_f;

   ptVlepJEC_new =         ptVlepJEC_f;
    yVlepJEC_new =        yVlepJEC_f;
    phiVlepJEC_new =         phiVlepJEC_f;
    massVlepJEC_new =        massVlepJEC_f;
    mtVlepJEC_new =        mtVlepJEC_f;

    RelIso03_NewID = RelIso03_NewID_f;


    npT = npT_f;
    MET_et = MET_et_f;
    Common_event_lepSF = Common_event_lepSF_f;



}

#endif
