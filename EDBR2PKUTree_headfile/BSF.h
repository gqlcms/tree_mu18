#ifndef _BSF_
#define _BSF_

void EDBR2PKUTree::BSF(Int_t IsData) {   
   //--BSF----------------------------
      btagweight_center=1.0;
      btagweight_up=1.0;
      btagweight_down=1.0;
        btagweight_center_deep=1.0;
        btagweight_up_deep=1.0;
        btagweight_down_deep=1.0;
        btagweight_center_tig_deep=1.0;
        btagweight_up_tig_deep=1.0;
        btagweight_down_tig_deep=1.0;
        btagweight_center_los_deep=1.0;
        btagweight_up_los_deep=1.0;
        btagweight_down_los_deep=1.0;

        

        if(IsData>0){
      double  bweight=1.0, dweight=1.0;
      double  bweightup=1.0;
      double  bweightdown=1.0;
            double  bweight_deep=1.0, dweight_deep=1.0;
            double  bweightup_deep=1.0;
            double  bweightdown_deep=1.0;
            double  bweight_tig_deep=1.0, dweight_tig_deep=1.0;
            double  bweightup_tig_deep=1.0;
            double  bweightdown_tig_deep=1.0;
            double  bweight_los_deep=1.0, dweight_los_deep=1.0;
            double  bweightup_los_deep=1.0;
            double  bweightdown_los_deep=1.0;

      double  beff=1.0, ceff=1.0, leff=1.0;
      
        
      for(Int_t i=0; i<8; i++)  {
       deltaRAK4AK8_new[i]=0.;
       deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
          if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));

       if(ak4jet_pt[i]>30 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1||deltaRAK4AK8_new_2[i]>=0.8)) {
         if(abs(ak4jet_hf[i])==5 && ak4jet_pt[i]<3000.) 
          { 
            double jet_sf    = bsv(1,ak4jet_pt[i]); 
            double jet_sfu    = bsv(2,ak4jet_pt[i]);
            double jet_sfd    = bsv(3,ak4jet_pt[i]);
              double jet_sf_deep    = bsv_deep(1,ak4jet_pt[i]);
              double jet_sfu_deep    = bsv_deep(2,ak4jet_pt[i]);
              double jet_sfd_deep    = bsv_deep(3,ak4jet_pt[i]);
              double jet_sf_tig_deep    = bsv_tig_deep(1,ak4jet_pt[i]);
              double jet_sfu_tig_deep    = bsv_tig_deep(2,ak4jet_pt[i]);
              double jet_sfd_tig_deep    = bsv_tig_deep(3,ak4jet_pt[i]);
              double jet_sf_los_deep    = bsv_los_deep(1,ak4jet_pt[i]);
              double jet_sfu_los_deep    = bsv_los_deep(2,ak4jet_pt[i]);
              double jet_sfd_los_deep    = bsv_los_deep(3,ak4jet_pt[i]);

//            int bbin=hb->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            beff=hb->GetBinContent(bbin);
              //b jet middle wp csv_2
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.647332650923;}
              else if (ak4jet_eta[i]<-0.8){beff=0.574428968825;}
              else if (ak4jet_eta[i]>-0.8){beff=0.580733478893;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.681764680183;}
              else if (ak4jet_eta[i]<-0.8){beff=0.607928437101;}
              else if (ak4jet_eta[i]>-0.8){beff=0.615291984235;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.695375820567;}
              else if (ak4jet_eta[i]<-0.8){beff=0.622144210455;}
              else if (ak4jet_eta[i]>-0.8){beff=0.629384525268;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.696089410422;}
              else if (ak4jet_eta[i]<-0.8){beff=0.619753472781;}
              else if (ak4jet_eta[i]>-0.8){beff=0.627676830787;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.680199371866;}
              else if (ak4jet_eta[i]<-0.8){beff=0.611069273396;}
              else if (ak4jet_eta[i]>-0.8){beff=0.616045887084;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.639681418349;}
              else if (ak4jet_eta[i]<-0.8){beff=0.563957042283;}
              else if (ak4jet_eta[i]>-0.8){beff=0.568304340723;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.560957890044;}
              else if (ak4jet_eta[i]<-0.8){beff=0.523701843813;}
              else if (ak4jet_eta[i]>-0.8){beff=0.526011078231;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.46878094322;}
              else if (ak4jet_eta[i]<-0.8){beff=0.47157152604;}
              else if (ak4jet_eta[i]>-0.8){beff=0.471656125556;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.41957313899;}
              else if (ak4jet_eta[i]<-0.8){beff=0.450354609929;}
              else if (ak4jet_eta[i]>-0.8){beff=0.425266903915;}
              }
              if(ak4jet_icsv[i]>0.8838) {
                  bweight=bweight*beff*jet_sf; dweight=dweight*beff;
                  bweightup=bweightup*beff*jet_sfu;
                  bweightdown=bweightdown*beff*jet_sfd;
              }
              else {
                  bweight=bweight*(1-beff*jet_sf); dweight=dweight*(1-beff);
                  bweightup=bweightup*(1-beff*jet_sfu);
                  bweightdown=bweightdown*(1-beff*jet_sfd);
              }
              //b jet middle wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.69289054503;}
              else if (ak4jet_eta[i]<-0.8){beff=0.626876762927;}
              else if (ak4jet_eta[i]>-0.8){beff=0.637414613848;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.726955160895;}
              else if (ak4jet_eta[i]<-0.8){beff=0.663040270851;}
              else if (ak4jet_eta[i]>-0.8){beff=0.672930844849;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.739178943248;}
              else if (ak4jet_eta[i]<-0.8){beff=0.671524512238;}
              else if (ak4jet_eta[i]>-0.8){beff=0.681053968587;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.742580432054;}
              else if (ak4jet_eta[i]<-0.8){beff=0.670239757334;}
              else if (ak4jet_eta[i]>-0.8){beff=0.677429556976;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.735916249616;}
              else if (ak4jet_eta[i]<-0.8){beff=0.662461650933;}
              else if (ak4jet_eta[i]>-0.8){beff=0.667615935214;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.706052668703;}
              else if (ak4jet_eta[i]<-0.8){beff=0.628799388706;}
              else if (ak4jet_eta[i]>-0.8){beff=0.631470724493;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.652342361457;}
              else if (ak4jet_eta[i]<-0.8){beff=0.587623735846;}
              else if (ak4jet_eta[i]>-0.8){beff=0.588928672311;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.561566159168;}
              else if (ak4jet_eta[i]<-0.8){beff=0.537073003542;}
              else if (ak4jet_eta[i]>-0.8){beff=0.532968340871;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.498954521694;}
              else if (ak4jet_eta[i]<-0.8){beff=0.479573712256;}
              else if (ak4jet_eta[i]>-0.8){beff=0.467857142857;}
              }
              if(ak4jet_icsv[i]>0.6321) {
                  bweight_deep=bweight_deep*beff*jet_sf_deep; dweight_deep=dweight_deep*beff;
                  bweightup_deep=bweightup_deep*beff*jet_sfu_deep;
                  bweightdown_deep=bweightdown_deep*beff*jet_sfd_deep;
              }
              else {
                  bweight_deep=bweight_deep*(1-beff*jet_sf_deep); dweight_deep=dweight_deep*(1-beff);
                  bweightup_deep=bweightup_deep*(1-beff*jet_sfu_deep);
                  bweightdown_deep=bweightdown_deep*(1-beff*jet_sfd_deep);
              }
              
              
              //b jet loose wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.841210764847;}
              else if (ak4jet_eta[i]<-0.8){beff=0.810264618953;}
              else if (ak4jet_eta[i]>-0.8){beff=0.816790149527;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.864728995719;}
              else if (ak4jet_eta[i]<-0.8){beff=0.834669348882;}
              else if (ak4jet_eta[i]>-0.8){beff=0.839980838846;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.876020129011;}
              else if (ak4jet_eta[i]<-0.8){beff=0.844805459205;}
              else if (ak4jet_eta[i]>-0.8){beff=0.849756794514;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.882550883747;}
              else if (ak4jet_eta[i]<-0.8){beff=0.851453339837;}
              else if (ak4jet_eta[i]>-0.8){beff=0.854276310286;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.884411924905;}
              else if (ak4jet_eta[i]<-0.8){beff=0.857190397309;}
              else if (ak4jet_eta[i]>-0.8){beff=0.858229372025;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.877699871512;}
              else if (ak4jet_eta[i]<-0.8){beff=0.853839246023;}
              else if (ak4jet_eta[i]>-0.8){beff=0.854308470831;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.866537339179;}
              else if (ak4jet_eta[i]<-0.8){beff=0.855038623192;}
              else if (ak4jet_eta[i]>-0.8){beff=0.854764213243;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.851237560602;}
              else if (ak4jet_eta[i]<-0.8){beff=0.866103050383;}
              else if (ak4jet_eta[i]>-0.8){beff=0.859676267555;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.861474124412;}
              else if (ak4jet_eta[i]<-0.8){beff=0.852575488455;}
              else if (ak4jet_eta[i]>-0.8){beff=0.857142857143;}
              }
              if(ak4jet_icsv[i]>0.2217) {
                  bweight_los_deep=bweight_los_deep*beff*jet_sf_los_deep; dweight_los_deep=dweight_los_deep*beff;
                  bweightup_los_deep=bweightup_los_deep*beff*jet_sfu_los_deep;
                  bweightdown_los_deep=bweightdown_los_deep*beff*jet_sfd_los_deep;
              }
              else {
                  bweight_los_deep=bweight_los_deep*(1-beff*jet_sf_los_deep); dweight_los_deep=dweight_los_deep*(1-beff);
                  bweightup_los_deep=bweightup_los_deep*(1-beff*jet_sfu_los_deep);
                  bweightdown_los_deep=bweightdown_los_deep*(1-beff*jet_sfd_los_deep);
              }
              
              
              //b jet tight wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.520869700243;}
              else if (ak4jet_eta[i]<-0.8){beff=0.43937276358;}
              else if (ak4jet_eta[i]>-0.8){beff=0.452455595379;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.563580454628;}
              else if (ak4jet_eta[i]<-0.8){beff=0.483288908117;}
              else if (ak4jet_eta[i]>-0.8){beff=0.495775694296;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.57455874948;}
              else if (ak4jet_eta[i]<-0.8){beff=0.492246553058;}
              else if (ak4jet_eta[i]>-0.8){beff=0.503952520111;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.571409970752;}
              else if (ak4jet_eta[i]<-0.8){beff=0.486007686416;}
              else if (ak4jet_eta[i]>-0.8){beff=0.495542760314;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.554259448319;}
              else if (ak4jet_eta[i]<-0.8){beff=0.467166648355;}
              else if (ak4jet_eta[i]>-0.8){beff=0.474168752141;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.501582653719;}
              else if (ak4jet_eta[i]<-0.8){beff=0.409965311191;}
              else if (ak4jet_eta[i]>-0.8){beff=0.415229321835;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.411013411297;}
              else if (ak4jet_eta[i]<-0.8){beff=0.337036158045;}
              else if (ak4jet_eta[i]>-0.8){beff=0.339300846649;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.292875166567;}
              else if (ak4jet_eta[i]<-0.8){beff=0.261624585856;}
              else if (ak4jet_eta[i]>-0.8){beff=0.259223994287;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.2232096184;}
              else if (ak4jet_eta[i]<-0.8){beff=0.241563055062;}
              else if (ak4jet_eta[i]>-0.8){beff=0.180357142857;}
              }
              if(ak4jet_icsv[i]>0.8953) {
                  bweight_tig_deep=bweight_tig_deep*beff*jet_sf_tig_deep; dweight_tig_deep=dweight_tig_deep*beff;
                  bweightup_tig_deep=bweightup_tig_deep*beff*jet_sfu_tig_deep;
                  bweightdown_tig_deep=bweightdown_tig_deep*beff*jet_sfd_tig_deep;
              }
              else {
                  bweight_tig_deep=bweight_tig_deep*(1-beff*jet_sf_tig_deep); dweight_tig_deep=dweight_tig_deep*(1-beff);
                  bweightup_tig_deep=bweightup_tig_deep*(1-beff*jet_sfu_tig_deep);
                  bweightdown_tig_deep=bweightdown_tig_deep*(1-beff*jet_sfd_tig_deep);
              }
              
              

          }
        else if(abs(ak4jet_hf[i])==4 && ak4jet_pt[i]<3000.)
          {   
            double jet_sf    = csv(1,ak4jet_pt[i]);
            double jet_sfu    = csv(2,ak4jet_pt[i]);
            double jet_sfd    = csv(3,ak4jet_pt[i]);
              double jet_sf_deep    = csv_deep(1,ak4jet_pt[i]);
              double jet_sfu_deep    = csv_deep(2,ak4jet_pt[i]);
              double jet_sfd_deep    = csv_deep(3,ak4jet_pt[i]);
              double jet_sf_tig_deep    = csv_tig_deep(1,ak4jet_pt[i]);
              double jet_sfu_tig_deep    = csv_tig_deep(2,ak4jet_pt[i]);
              double jet_sfd_tig_deep    = csv_tig_deep(3,ak4jet_pt[i]);
              double jet_sf_los_deep    = csv_los_deep(1,ak4jet_pt[i]);
              double jet_sfu_los_deep    = csv_los_deep(2,ak4jet_pt[i]);
              double jet_sfd_los_deep    = csv_los_deep(3,ak4jet_pt[i]);

//            int cbin=hc->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            ceff=hc->GetBinContent(cbin);
              //c jet middle wp csv_2
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.134641720495;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.121166477571;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.121321670268;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.131091961461;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.117347838189;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.118499893002;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.136123314875;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.124065366633;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.125565528143;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.143276930452;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.130265186364;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.132288574989;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.141473999343;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.137933649736;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.139164307147;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.134456602805;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.129486943226;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.130881506146;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.12118518329;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.139912735705;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.136867817858;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.125290422511;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.152911249293;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.144447536877;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.120165745856;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.174418604651;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.156756756757;}
              }
              if(ak4jet_icsv[i]>0.8838) {
                  bweight=bweight*ceff*jet_sf; dweight=dweight*ceff;
                  bweightup=bweightup*ceff*jet_sfu;
                  bweightdown=bweightdown*ceff*jet_sfd;
              }
              else {
                  bweight=bweight*(1-ceff*jet_sf); dweight=dweight*(1-ceff);
                  bweightup=bweightup*(1-ceff*jet_sfu);
                  bweightdown=bweightdown*(1-ceff*jet_sfd);
              }
              //c jet middle wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.129222981653;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.12950337414;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.131247230564;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.117511658516;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.122105267909;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.123202524755;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.11615914238;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.118675833028;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.119520303458;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.121660554992;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.121737542076;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.121667035755;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.129610177752;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.129728045772;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.12830101997;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.135430225509;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.130401636845;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.128649540396;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.144650004058;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.141602206632;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.13753036073;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.152528207271;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.15451577802;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.160846851801;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.138200782269;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.154255319149;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.178378378378;}
              }
              if(ak4jet_icsv[i]>0.6321) {
                  bweight_deep=bweight_deep*ceff*jet_sf_deep; dweight_deep=dweight_deep*ceff;
                  bweightup_deep=bweightup_deep*ceff*jet_sfu_deep;
                  bweightdown_deep=bweightdown_deep*ceff*jet_sfd_deep;
              }
              else {
                  bweight_deep=bweight_deep*(1-ceff*jet_sf_deep); dweight_deep=dweight_deep*(1-ceff);
                  bweightup_deep=bweightup_deep*(1-ceff*jet_sfu_deep);
                  bweightdown_deep=bweightdown_deep*(1-ceff*jet_sfd_deep);
              }
              
              
              //c jet loose wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.407532419829;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.414756479434;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.41552030184;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.396198586996;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.406929414522;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.40737402914;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.395266574669;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.405871869675;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.406265413537;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.404117752108;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.413086131753;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.411639878628;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.417210082833;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.429569432325;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.425682017993;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.42772089264;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.442946517278;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.436942165497;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.464082567948;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.496432212029;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.486341419533;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.504972837443;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.55359085963;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.554577948859;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.552803129074;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.622340425532;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.578378378378;}
              }
              if(ak4jet_icsv[i]>0.2217) {
                  bweight_los_deep=bweight_los_deep*ceff*jet_sf_los_deep; dweight_los_deep=dweight_los_deep*ceff;
                  bweightup_los_deep=bweightup_los_deep*ceff*jet_sfu_los_deep;
                  bweightdown_los_deep=bweightdown_los_deep*ceff*jet_sfd_los_deep;
              }
              else {
                  bweight_los_deep=bweight_los_deep*(1-ceff*jet_sf_los_deep); dweight_los_deep=dweight_los_deep*(1-ceff);
                  bweightup_los_deep=bweightup_los_deep*(1-ceff*jet_sfu_los_deep);
                  bweightdown_los_deep=bweightdown_los_deep*(1-ceff*jet_sfd_los_deep);
              }
              
              
              //c jet tight wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0235693322031;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0245835119498;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0251726314782;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0221002199997;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0236972600643;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0241362178695;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0227875353133;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0234240379638;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0238317547109;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0248622843829;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0242885304979;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0248959753838;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0273376534458;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0269927593883;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0266163842966;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0290333494669;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0261314787673;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0257881920605;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0297079776972;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.028077591893;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0269123511949;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0308399498537;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0329162132753;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0299697552928;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0286831812256;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0212765957447;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0540540540541;}
              }
              if(ak4jet_icsv[i]>0.8953) {
                  bweight_tig_deep=bweight_tig_deep*ceff*jet_sf_tig_deep; dweight_tig_deep=dweight_tig_deep*ceff;
                  bweightup_tig_deep=bweightup_tig_deep*ceff*jet_sfu_tig_deep;
                  bweightdown_tig_deep=bweightdown_tig_deep*ceff*jet_sfd_tig_deep;
              }
              else {
                  bweight_tig_deep=bweight_tig_deep*(1-ceff*jet_sf_tig_deep); dweight_tig_deep=dweight_tig_deep*(1-ceff);
                  bweightup_tig_deep=bweightup_tig_deep*(1-ceff*jet_sfu_tig_deep);
                  bweightdown_tig_deep=bweightdown_tig_deep*(1-ceff*jet_sfd_tig_deep);
              }

          }
        else if (abs(ak4jet_hf[i])==0 && ak4jet_pt[i]<3000.)
          {  
            double jet_sf    = lsv(1,ak4jet_pt[i]);
            double jet_sfu    = lsv(2,ak4jet_pt[i]);
            double jet_sfd    = lsv(3,ak4jet_pt[i]);
              double jet_sf_deep    = lsv_deep(1,ak4jet_pt[i]);
              double jet_sfu_deep    = lsv_deep(2,ak4jet_pt[i]);
              double jet_sfd_deep    = lsv_deep(3,ak4jet_pt[i]);
              double jet_sf_tig_deep    = lsv_tig_deep(1,ak4jet_pt[i]);
              double jet_sfu_tig_deep    = lsv_tig_deep(2,ak4jet_pt[i]);
              double jet_sfd_tig_deep    = lsv_tig_deep(3,ak4jet_pt[i]);
              double jet_sf_los_deep    = lsv_los_deep(1,ak4jet_pt[i]);
              double jet_sfu_los_deep    = lsv_los_deep(2,ak4jet_pt[i]);
              double jet_sfd_los_deep    = lsv_los_deep(3,ak4jet_pt[i]);

//            int lbin=hl->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            leff=hl->GetBinContent(lbin);
              //l jet middle wp csv_2
              if (ak4jet_pt[i]>=20.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00883804539947;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0126624606715;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0118134413393;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<500.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0108378571471;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0209253471333;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0202266537294;}
              }
              else if (ak4jet_pt[i]>=500.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0194805814586;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0362205015431;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0329325280287;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0339058666072;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0450112528132;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0479452054795;}
              }
              if(ak4jet_icsv[i]>0.8838) {
                  bweight=bweight*leff*jet_sf; dweight=dweight*leff;
                  bweightup=bweightup*leff*jet_sfu;
                  bweightdown=bweightdown*leff*jet_sfd;
              }
              else {
                  bweight=bweight*(1-leff*jet_sf); dweight=dweight*(1-leff);
                  bweightup=bweightup*(1-leff*jet_sfu);
                  bweightdown=bweightdown*(1-leff*jet_sfd);
              }
              //l jet middle wp DeepCSV
              if (ak4jet_pt[i]>=20.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.011375461494;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0166255497284;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0163500521057;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<500.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0136983493182;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0207888628102;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0196059086432;}
              }
              else if (ak4jet_pt[i]>=500.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0297503594554;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0386133774486;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0353440504808;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0398640996602;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0518292682927;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0468509984639;}
              }
              if(ak4jet_icsv[i]>0.6321) {
                  bweight_deep=bweight_deep*leff*jet_sf_deep; dweight_deep=dweight_deep*leff;
                  bweightup_deep=bweightup_deep*leff*jet_sfu_deep;
                  bweightdown_deep=bweightdown_deep*leff*jet_sfd_deep;
              }
              else {
                  bweight_deep=bweight_deep*(1-leff*jet_sf_deep); dweight_deep=dweight_deep*(1-leff);
                  bweightup_deep=bweightup_deep*(1-leff*jet_sfu_deep);
                  bweightdown_deep=bweightdown_deep*(1-leff*jet_sfd_deep);
              }
              
              
              //l jet loose wp DeepCSV
              if (ak4jet_pt[i]>=20.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0955787532571;}
              else if (ak4jet_eta[i]<-0.8){leff=0.148789610047;}
              else if (ak4jet_eta[i]>-0.8){leff=0.145664194575;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<500.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.113664749688;}
              else if (ak4jet_eta[i]<-0.8){leff=0.180917356463;}
              else if (ak4jet_eta[i]>-0.8){leff=0.173081118715;}
              }
              else if (ak4jet_pt[i]>=500.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.216589950689;}
              else if (ak4jet_eta[i]<-0.8){leff=0.2958792345;}
              else if (ak4jet_eta[i]>-0.8){leff=0.286283052885;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.334314835787;}
              else if (ak4jet_eta[i]<-0.8){leff=0.383384146341;}
              else if (ak4jet_eta[i]>-0.8){leff=0.357142857143;}
              }
              if(ak4jet_icsv[i]>0.2217) {
                  bweight_los_deep=bweight_los_deep*leff*jet_sf_los_deep; dweight_los_deep=dweight_los_deep*leff;
                  bweightup_los_deep=bweightup_los_deep*leff*jet_sfu_los_deep;
                  bweightdown_los_deep=bweightdown_los_deep*leff*jet_sfd_los_deep;
              }
              else {
                  bweight_los_deep=bweight_los_deep*(1-leff*jet_sf_los_deep); dweight_los_deep=dweight_los_deep*(1-leff);
                  bweightup_los_deep=bweightup_los_deep*(1-leff*jet_sfu_los_deep);
                  bweightdown_los_deep=bweightdown_los_deep*(1-leff*jet_sfd_los_deep);
              }
              
              
              //l jet tight wp DeepCSV
              if (ak4jet_pt[i]>=20.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00115701280759;}
              else if (ak4jet_eta[i]<-0.8){leff=0.00153811907691;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0015196417698;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<500.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00174898076121;}
              else if (ak4jet_eta[i]<-0.8){leff=0.00228349093204;}
              else if (ak4jet_eta[i]>-0.8){leff=0.00217531498754;}
              }
              else if (ak4jet_pt[i]>=500.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00353182989318;}
              else if (ak4jet_eta[i]<-0.8){leff=0.00475617550852;}
              else if (ak4jet_eta[i]>-0.8){leff=0.00368088942308;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00724801812005;}
              else if (ak4jet_eta[i]<-0.8){leff=0.00533536585366;}
              else if (ak4jet_eta[i]>-0.8){leff=0.00614439324117;}
              }
              if(ak4jet_icsv[i]>0.8953) {
                  bweight_tig_deep=bweight_tig_deep*leff*jet_sf_tig_deep; dweight_tig_deep=dweight_tig_deep*leff;
                  bweightup_tig_deep=bweightup_tig_deep*leff*jet_sfu_tig_deep;
                  bweightdown_tig_deep=bweightdown_tig_deep*leff*jet_sfd_tig_deep;
              }
              else {
                  bweight_tig_deep=bweight_tig_deep*(1-leff*jet_sf_tig_deep); dweight_tig_deep=dweight_tig_deep*(1-leff);
                  bweightup_tig_deep=bweightup_tig_deep*(1-leff*jet_sfu_tig_deep);
                  bweightdown_tig_deep=bweightdown_tig_deep*(1-leff*jet_sfd_tig_deep);
              }

          }
      //std::cout<<" yy1 "<<abs(ak4jet_hf[i])<<" "<<dweight<<" "<<bweight<<std::endl;
        } 
      //std::cout<<" yy2 "<<abs(ak4jet_hf[i])<<" "<<ak4jet_pt[i]<<" "<<fabs(ak4jet_eta[i])<<" "<<ak4jet_IDLoose[i]<<" "<<deltaRAK4AK8[i]<<std::endl;
     }
//     input->Close();
//     std::cout<<" xx "<<bweight/dweight<<" "<<bweightup/dweight<<" "<<bweightdown/dweight<<std::endl;
        

     btagweight_center=bweight/dweight;
     btagweight_up=bweightup/dweight;
     btagweight_down=bweightdown/dweight;
     btagweight_center_deep=bweight_deep/dweight_deep;
     btagweight_up_deep=bweightup_deep/dweight_deep;
     btagweight_down_deep=bweightdown_deep/dweight_deep;

            btagweight_center_tig_deep=bweight_tig_deep/dweight_tig_deep;
            btagweight_up_tig_deep=bweightup_tig_deep/dweight_tig_deep;
            btagweight_down_tig_deep=bweightdown_tig_deep/dweight_tig_deep;
            btagweight_center_los_deep=bweight_los_deep/dweight_los_deep;
            btagweight_up_los_deep=bweightup_los_deep/dweight_los_deep;
            btagweight_down_los_deep=bweightdown_los_deep/dweight_los_deep;

//     std::cout<<btagweight_center<<std::endl;
 
     
//--BSF----------------------------
        }
}

#endif