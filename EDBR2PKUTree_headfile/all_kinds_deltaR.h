#ifndef _all_kinds_deltaR_
#define _all_kinds_deltaR_

void EDBR2PKUTree::all_kinds_deltaR() {   
        deltaRjet1jet2=sqrt(pow(fabs(jetAK8puppi_eta_2-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(jetAK8puppi_phi_2-jetAK8puppi_phi),2*Pi-fabs(jetAK8puppi_phi_2-jetAK8puppi_phi)),2));
        deltaetajet1jet2=jetAK8puppi_eta-jetAK8puppi_eta_2;
        deltaetajet1lep=jetAK8puppi_eta-etalep1;
        deltaetajet2lep=jetAK8puppi_eta_2-etalep1;
        delPhilepmet_m=TMath::Min(fabs(philep1-MET_phi),2*Pi-fabs(philep1-MET_phi));
        delPhiWjet1=TMath::Min(fabs(W_phi-jetAK8puppi_phi),2*Pi-fabs(W_phi-jetAK8puppi_phi));
        delPhiWjet2=TMath::Min(fabs(W_phi-jetAK8puppi_phi_2),2*Pi-fabs(W_phi-jetAK8puppi_phi_2));
        delPhiWjetclose=(fabs(delPhiWjet1)<fabs(delPhiWjet2))?delPhiWjet1:delPhiWjet2;
        delPhiWjetfar=(fabs(delPhiWjet1)>fabs(delPhiWjet2))?delPhiWjet1:delPhiWjet2;
}

void EDBR2PKUTree::all_kinds_deltaR_unknown() {   
        for(int i=0;i<8;i++){
            deltaRjet1ak4[i]=-99;deltaRjet1ak4_1[i]=-99;deltaRphiak4[i]=-99;}
        for(Int_t i=0;i<8;i++){
        deltaRjet1ak4[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
        deltaRjet1ak4_1[i]=deltaR(ak4jet_eta[i],ak4jet_phi[i],jetAK8puppi_eta,jetAK8puppi_phi);

        }
        Double_t tmp_deltaRjet1ak4=999,tmp_deltaRjet1ak4_1=999;
        Int_t tmp_mini=-99,tmp_mini_1=-99;
        for(Int_t i=0;i<8;i++){
            if((tmp_deltaRjet1ak4>=deltaRjet1ak4[i])&&deltaRjet1ak4[i]<1.6&&deltaRjet1ak4[i]>0.8&&deltaRjet1ak4[i]>0&&ak4jet_pt[i]>0)
            {tmp_mini=i;
                tmp_deltaRjet1ak4=deltaRjet1ak4[i];}
            if((tmp_deltaRjet1ak4_1>=deltaRjet1ak4_1[i])&&deltaRjet1ak4_1[i]<1.6&&deltaRjet1ak4_1[i]>0.8&&ak4jet_pt[i]>0)
            {tmp_mini_1=i;
                tmp_deltaRjet1ak4_1=deltaRjet1ak4[i];}


        }
        jet_mass_puppi_casea=-99;        jet_mass_puppi_casea_1=-99;

        //if(tmp_mini<0) cout<<tmp_mini<<"  "<<jet_mass_puppi_casea<<endl;
        if(tmp_mini>0){
            TLorentzVector ak8jet1,ak4closejet1;
            ak8jet1.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
            ak4closejet1.SetPtEtaPhiE(ak4jet_pt[tmp_mini],ak4jet_eta[tmp_mini],ak4jet_phi[tmp_mini],ak4jet_e[tmp_mini]);
            jet_mass_puppi_casea = (ak8jet1+ak4closejet1).Mag();

        }
        //cout<<tmp_mini<<"  "<<jet_mass_puppi_casea<<"  "<<jet_mass_puppi_corr<<endl;
        numak4outsideak8=0;
        Double_t tmpak4jet_pt[8],tmpak4jet_eta[8],tmpak4jet_phi[8],tmpak4jet_e[8];
        Double_t tmpdRjet1ak4[8]={999,999,999,999,999,999,999,999};
        for(Int_t i=0;i<8;i++){
            if(deltaRjet1ak4_1[i]<1.6&&deltaRjet1ak4_1[i]>0.8&&ak4jet_pt[i]>0)
            {
                tmpak4jet_pt[numak4outsideak8]=ak4jet_pt[i];
                tmpak4jet_eta[numak4outsideak8]=ak4jet_eta[i];
                tmpak4jet_phi[numak4outsideak8]=ak4jet_phi[i];
                tmpak4jet_e[numak4outsideak8]=ak4jet_e[i];
                tmpdRjet1ak4[numak4outsideak8]=deltaRjet1ak4_1[i];
            numak4outsideak8++;
            }
        }
        if(numak4outsideak8==0) jet_mass_puppi_casea_1=-99;

        if(numak4outsideak8==1)
        {   TLorentzVector ak8jet1_1,ak4closejet1_1;
            ak8jet1_1.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
            ak4closejet1_1.SetPtEtaPhiE(tmpak4jet_pt[0],tmpak4jet_eta[0],tmpak4jet_phi[0],tmpak4jet_e[0]);
            jet_mass_puppi_casea_1 = (ak8jet1_1+ak4closejet1_1).Mag();
}
        if(numak4outsideak8>=2){
            Int_t *indexx=new Int_t[8];
            TMath::Sort(8,tmpdRjet1ak4,indexx,0);
            TLorentzVector ak8jet1_1,ak4closejet1_1,ak4closejet1_2;
            ak8jet1_1.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
            ak4closejet1_1.SetPtEtaPhiE(ak4jet_pt[indexx[0]],ak4jet_eta[indexx[0]],ak4jet_phi[indexx[0]],ak4jet_e[indexx[0]]);
            ak4closejet1_2.SetPtEtaPhiE(ak4jet_pt[indexx[1]],ak4jet_eta[indexx[1]],ak4jet_phi[indexx[1]],ak4jet_e[indexx[1]]);
            
            jet_mass_puppi_casea_1 = (ak8jet1_1+ak4closejet1_1+ak4closejet1_2).Mag();
            
        }
 

 
        for(Int_t i=0;i<8;i++){
            //deltaRphiak4[i]=sqrt(pow(fabs(-ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-(-Pi+jetAK8puppi_phi)),2*Pi-fabs(Pi+ak4jet_phi[i]-jetAK8puppi_phi)),2));
            deltaRphiak4[i]=deltaR(-jetAK8puppi_eta,((jetAK8puppi_phi>0)?(-Pi+jetAK8puppi_phi):(Pi+jetAK8puppi_phi)),ak4jet_eta[i],ak4jet_phi[i]);
        }
        numak4insidephi=0;
        Double_t tmp1ak4jet_pt[8],tmp1ak4jet_eta[8],tmp1ak4jet_phi[8],tmp1ak4jet_e[8];
        Double_t tmpdRphiak4[8]={-999,-999,-999,-999,-999,-999,-999,-999};
        for(Int_t i=0;i<8;i++){
            if(deltaRphiak4[i]<1.6&&ak4jet_pt[i]>0)
            {
                tmp1ak4jet_pt[numak4outsideak8]=ak4jet_pt[i];
                tmp1ak4jet_eta[numak4outsideak8]=ak4jet_eta[i];
                tmp1ak4jet_phi[numak4outsideak8]=ak4jet_phi[i];
                tmp1ak4jet_e[numak4outsideak8]=ak4jet_e[i];
                tmpdRphiak4[numak4outsideak8]=tmpdRphiak4[i];
                numak4insidephi++;
            }
        }
        Int_t *indexa=new Int_t[8];
        TMath::Sort(8,deltaRphiak4,indexa,0);

        //if((status_1==4 &&(status_2!=4||status_3!=4) )&&(! (status_2==4&&status_3==4)))
        //cout<<jetAK8puppi_phi<<"   "<<((jetAK8puppi_phi>0)?(-Pi+jetAK8puppi_phi):(Pi+jetAK8puppi_phi))<<"  "<<deltaRphiak4[indexa[0]]<<" "<<deltaR(jetAK8puppi_eta,jetAK8puppi_phi,W_eta,W_phi)<<" "<<deltaR(jetAK8puppi_eta,jetAK8puppi_phi,ak4jet_eta[0],ak4jet_phi[0])<<" "<<numak4insidephi<<endl;
        if(numak4insidephi==0)
        {                    jet_mass_puppi_caseb = massVlepJEC;
        }
       
        if(numak4insidephi==1)
        {   TLorentzVector phireco_1,ak4closephi_1;
            phireco_1.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
            ak4closephi_1.SetPtEtaPhiE(tmp1ak4jet_pt[0],tmp1ak4jet_eta[0],tmp1ak4jet_phi[0],tmp1ak4jet_e[0]);
            jet_mass_puppi_caseb = (phireco_1+ak4closephi_1).Mag();
        }
        if(numak4insidephi>=2){
            Int_t *index=new Int_t[8];
            TMath::Sort(8,tmpdRphiak4,index,1);
            TLorentzVector phireco_1,ak4closephi_1,ak4closephi_2;
            phireco_1.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
            ak4closephi_1.SetPtEtaPhiE(tmp1ak4jet_pt[index[0]],tmp1ak4jet_eta[index[0]],tmp1ak4jet_phi[index[0]],tmp1ak4jet_e[index[0]]);
            ak4closephi_2.SetPtEtaPhiE(tmp1ak4jet_pt[index[1]],tmp1ak4jet_eta[index[1]],tmp1ak4jet_phi[index[1]],tmp1ak4jet_e[index[1]]);

            jet_mass_puppi_caseb = (phireco_1+ak4closephi_1+ak4closephi_2).Mag();
        }
}

void EDBR2PKUTree::two_w_boson_deltaR() {   
        dR_min_WW=-99;dPhi_min_WW=-99;dR_max_WW=-99;dPhi_max_WW=-99;dEta_max_WW=-99;dEta_max_WW=-99;
        dR_j1W=-99;dR_j2W=-99;dR_j3W=-99;dR_j1j2=-99;dR_j1j3=-99;dR_j2j3=-99;
        dPhi_j1W=-99;dPhi_j2W=-99;dPhi_j3W=-99;dPhi_j1j2=-99;dPhi_j1j3=-99;dPhi_j2j3=-99;
        dEta_j1W=-99;dEta_j2W=-99;dEta_j3W=-99;dEta_j1j2=-99;dEta_j1j3=-99;dEta_j2j3=-99;
        dR_min_WW_mass=-99;dPhi_min_WW_mass=-99;dR_max_WW_mass=-99;dPhi_max_WW_mass=-99;dEta_max_WW_mass=-99;dEta_max_WW_mass=-99;
        dR_jmaxW=-99;dR_jmidW=-99;dR_jminW=-99;dR_jmaxjmid=-99;dR_jmaxjmin=-99;dR_jmidjmin=-99;
        dPhi_jmaxW=-99;dPhi_jmidW=-99;dPhi_jminW=-99;dPhi_jmaxjmid=-99;dPhi_jmaxjmin=-99;dPhi_jmidjmin=-99;
        dEta_jmaxW=-99;dEta_jmidW=-99;dEta_jminW=-99;dEta_jmaxjmid=-99;dEta_jmaxjmin=-99;dEta_jmidjmin=-99;

        double dRWW[6],dPhiWW[6],dEtaWW[6];
        dRWW[0]=deltaR(W_eta,W_phi,jetAK8puppi_eta,jetAK8puppi_phi);
        dRWW[1]=deltaR(W_eta,W_phi,jetAK8puppi_eta_2,jetAK8puppi_phi_2);
        dRWW[3]=deltaR(W_eta,W_phi,jetAK8puppi_eta_3,jetAK8puppi_phi_3);
        dRWW[2]=deltaR(jetAK8puppi_eta,jetAK8puppi_phi,jetAK8puppi_eta_2,jetAK8puppi_phi_2);
        dRWW[4]=deltaR(jetAK8puppi_eta,jetAK8puppi_phi,jetAK8puppi_eta_3,jetAK8puppi_phi_3);
        dRWW[5]=deltaR(jetAK8puppi_eta_2,jetAK8puppi_phi_2,jetAK8puppi_eta_3,jetAK8puppi_phi_3);
        
        dPhiWW[0]=deltaPhi(W_phi,jetAK8puppi_phi);
        dPhiWW[1]=deltaPhi(W_phi,jetAK8puppi_phi_2);
        dPhiWW[3]=deltaPhi(W_phi,jetAK8puppi_phi_3);
        dPhiWW[2]=deltaPhi(jetAK8puppi_phi,jetAK8puppi_phi_2);
        dPhiWW[4]=deltaPhi(jetAK8puppi_phi,jetAK8puppi_phi_3);
        dPhiWW[5]=deltaPhi(jetAK8puppi_phi_2,jetAK8puppi_phi_3);

        dEtaWW[0]=(W_eta-jetAK8puppi_eta);
        dEtaWW[1]=(W_eta-jetAK8puppi_eta_2);
        dEtaWW[3]=(W_eta-jetAK8puppi_eta_3);
        dEtaWW[2]=(jetAK8puppi_eta-jetAK8puppi_eta_2);
        dEtaWW[4]=(jetAK8puppi_eta-jetAK8puppi_eta_3);
        dEtaWW[5]=(jetAK8puppi_eta_2-jetAK8puppi_eta_3);
        dR_j1W=dRWW[0];dR_j2W=dRWW[1];dR_j3W=dRWW[3];dR_j1j2=dRWW[2];dR_j1j3=dRWW[4];dR_j2j3=dRWW[5];
        dPhi_j1W=dPhiWW[0];dPhi_j2W=dPhiWW[1];dPhi_j3W=dPhiWW[3];dPhi_j1j2=dPhiWW[2];dPhi_j1j3=dPhiWW[4];dPhi_j2j3=dPhiWW[5];
        dEta_j1W=dEtaWW[0];dEta_j2W=dEtaWW[1];dEta_j3W=dEtaWW[3];dEta_j1j2=dEtaWW[2];dEta_j1j3=dEtaWW[4];dEta_j2j3=dEtaWW[5];
        if(Nj8==1) {dR_min_WW=dRWW[0];dR_max_WW=dRWW[0];
        dPhi_min_WW=dPhiWW[0];dPhi_max_WW=dPhiWW[0];
        dEta_min_WW=dEtaWW[0];dEta_max_WW=dEtaWW[0];}
        if(Nj8==2) {
            Int_t *indexxwwr=new Int_t[6];
            TMath::Sort(3,dRWW,indexxwwr,1);
            dR_min_WW=dRWW[indexxwwr[2]];dR_max_WW=dRWW[indexxwwr[0]];
            Int_t *indexxwwPhi=new Int_t[6];
            TMath::Sort(3,dPhiWW,indexxwwPhi,1);
            dPhi_min_WW=dPhiWW[indexxwwPhi[2]];dPhi_max_WW=dPhiWW[indexxwwPhi[0]];
            Int_t *indexxwwEta=new Int_t[6];
            TMath::Sort(3,dEtaWW,indexxwwEta,1);
            dEta_min_WW=dEtaWW[indexxwwEta[2]];dEta_max_WW=dEtaWW[indexxwwEta[0]];
        }
    if(Nj8>=3) {
        Int_t *indexxwwr=new Int_t[6];
        TMath::Sort(6,dRWW,indexxwwr,1);
        dR_min_WW=dRWW[indexxwwr[5]];dR_max_WW=dRWW[indexxwwr[0]];
        Int_t *indexxwwPhi=new Int_t[6];
        TMath::Sort(6,dPhiWW,indexxwwPhi,1);
        dPhi_min_WW=dPhiWW[indexxwwPhi[5]];dPhi_max_WW=dPhiWW[indexxwwPhi[0]];
        Int_t *indexxwwEta=new Int_t[6];
        TMath::Sort(6,dEtaWW,indexxwwEta,1);
        dEta_min_WW=dEtaWW[indexxwwEta[5]];dEta_max_WW=dEtaWW[indexxwwEta[0]];
    }
        //cout<<"NJ  " <<Nj8<<dR_min_WW<<"    "<<dR_max_WW<<"    "<<dR_j1W<<"    "<<dR_j2W<<"    "<<dR_j3W<<"    "<<dR_j1j2<<"    "<<dR_j1j3<<"    "<<dR_j2j3<<"    "<<endl;
    double dRWW_mass[6],dPhiWW_mass[6],dEtaWW_mass[6];
    dRWW_mass[0]=deltaR(W_eta,W_phi,eta_max,phi_max);
    dRWW_mass[3]=deltaR(W_eta,W_phi,eta_mid,phi_mid);
    dRWW_mass[1]=deltaR(W_eta,W_phi,eta_min,phi_min);
    dRWW_mass[4]=deltaR(eta_max,phi_max,eta_mid,phi_mid);
    dRWW_mass[2]=deltaR(eta_max,phi_max,eta_min,phi_min);
    dRWW_mass[5]=deltaR(eta_mid,phi_mid,eta_min,phi_min);
    
    dPhiWW_mass[0]=deltaPhi(W_phi,phi_max);
    dPhiWW_mass[3]=deltaPhi(W_phi,phi_mid);
    dPhiWW_mass[1]=deltaPhi(W_phi,phi_min);
    dPhiWW_mass[4]=deltaPhi(phi_max,phi_mid);
    dPhiWW_mass[2]=deltaPhi(phi_max,phi_min);
    dPhiWW_mass[5]=deltaPhi(phi_mid,phi_min);
    
    dEtaWW_mass[0]=(W_eta-eta_max);
    dEtaWW_mass[3]=(W_eta-eta_mid);
    dEtaWW_mass[2]=(W_eta-eta_min);
    dEtaWW_mass[4]=(eta_max-eta_mid);
    dEtaWW_mass[2]=(eta_max-eta_min);
    dEtaWW_mass[5]=(eta_mid-eta_min);
    dR_jmaxW=dRWW_mass[0];dR_jmidW=dRWW_mass[3];dR_jminW=dRWW_mass[1];dR_jmaxjmid=dRWW_mass[4];dR_jmaxjmin=dRWW_mass[2];dR_jmidjmin=dRWW_mass[5];
    dPhi_jmaxW=dPhiWW_mass[0];dPhi_jmidW=dPhiWW_mass[3];dPhi_jminW=dPhiWW_mass[2];dPhi_jmaxjmid=dPhiWW_mass[4];dPhi_jmaxjmin=dPhiWW_mass[2];dPhi_jmidjmin=dPhiWW_mass[5];
    dEta_jmaxW=dEtaWW_mass[0];dEta_jmidW=dEtaWW_mass[3];dEta_jminW=dEtaWW_mass[1];dEta_jmaxjmid=dEtaWW_mass[4];dEta_jmaxjmin=dEtaWW_mass[2];dEta_jmidjmin=dEtaWW_mass[5];
    if(Nj8==1) {dR_min_WW_mass=dRWW_mass[0];dR_max_WW_mass=dRWW_mass[0];
        dPhi_min_WW_mass=dPhiWW_mass[0];dPhi_max_WW_mass=dPhiWW_mass[0];
        dEta_min_WW_mass=dEtaWW_mass[0];dEta_max_WW_mass=dEtaWW_mass[0];}
if(Nj8==2) {
    Int_t *indexxwwr_mass=new Int_t[6];
    TMath::Sort(3,dRWW_mass,indexxwwr_mass,1);
    dR_min_WW_mass=dRWW_mass[indexxwwr_mass[2]];dR_max_WW_mass=dRWW_mass[indexxwwr_mass[0]];
    Int_t *indexxwwPhi_mass=new Int_t[6];
    TMath::Sort(3,dPhiWW_mass,indexxwwPhi_mass,1);
    dPhi_min_WW_mass=dPhiWW_mass[indexxwwPhi_mass[2]];dPhi_max_WW_mass=dPhiWW_mass[indexxwwPhi_mass[0]];
    Int_t *indexxwwEta_mass=new Int_t[6];
    TMath::Sort(3,dEtaWW_mass,indexxwwEta_mass,1);
    dEta_min_WW_mass=dEtaWW_mass[indexxwwEta_mass[2]];dEta_max_WW_mass=dEtaWW_mass[indexxwwEta_mass[0]];
}
if(Nj8>=3) {
    Int_t *indexxwwr_mass=new Int_t[6];
    TMath::Sort(6,dRWW_mass,indexxwwr_mass,1);
    dR_min_WW_mass=dRWW_mass[indexxwwr_mass[5]];dR_max_WW_mass=dRWW_mass[indexxwwr_mass[0]];
    Int_t *indexxwwPhi_mass=new Int_t[6];
    TMath::Sort(6,dPhiWW_mass,indexxwwPhi_mass,1);
    dPhi_min_WW_mass=dPhiWW_mass[indexxwwPhi_mass[5]];dPhi_max_WW_mass=dPhiWW_mass[indexxwwPhi_mass[0]];
    Int_t *indexxwwEta_mass=new Int_t[6];
    TMath::Sort(6,dEtaWW_mass,indexxwwEta_mass,1);
    dEta_min_WW_mass=dEtaWW_mass[indexxwwEta_mass[5]];dEta_max_WW_mass=dEtaWW_mass[indexxwwEta_mass[0]];
}
}

#endif