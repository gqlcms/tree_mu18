#ifndef _in_and_ex_ak4_
#define _in_and_ex_ak4_

void EDBR2PKUTree::in_and_ex_ak4() {   
           double deltaRjet1ak4p[8]={-99,-99,-99,-99,-99,-99,-99,-99};
        double deltaRjet2ak4p[8]={-99,-99,-99,-99,-99,-99,-99,-99};
        double deltaRjet3ak4p[8]={-99,-99,-99,-99,-99,-99,-99,-99};

        TLorentzVector tempak4,ak8jet1_2;
        ak8jet1_2.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
        
        Nj4ex=0;
        Nj4in=0;
            ak4expt1=-99,ak4exeta1=-99,ak4exphi1=-99,ak4exe1=-99;
            ak4expt2=-99,ak4exeta2=-99,ak4exphi2=-99,ak4exe2=-99;
            
            Nj4ex_dR0p5=0;
            Nj4in_dR0p5=0;
            ak4expt1_dR0p5=-99,ak4exeta1_dR0p5=-99,ak4exphi1_dR0p5=-99,ak4exe1_dR0p5=-99;
            ak4expt2_dR0p5=-99,ak4exeta2_dR0p5=-99,ak4exphi2_dR0p5=-99,ak4exe2_dR0p5=-99;

        
        Nj4ex_dR0p8=0;
        Nj4in_dR0p8=0;
        ak4expt1_dR0p8=-99,ak4exeta1_dR0p8=-99,ak4exphi1_dR0p8=-99,ak4exe1_dR0p8=-99;
        ak4expt2_dR0p8=-99,ak4exeta2_dR0p8=-99,ak4exphi2_dR0p8=-99,ak4exe2_dR0p8=-99;
        
        Nj4ex_dR1p0=0;
        Nj4in_dR1p0=0;
        ak4expt1_dR1p0=-99,ak4exeta1_dR1p0=-99,ak4exphi1_dR1p0=-99,ak4exe1_dR1p0=-99;
        ak4expt2_dR1p0=-99,ak4exeta2_dR1p0=-99,ak4exphi2_dR1p0=-99,ak4exe2_dR1p0=-99;
        
        Nj4ex_dR1p2=0;
        Nj4in_dR1p2=0;
        ak4expt1_dR1p2=-99,ak4exeta1_dR1p2=-99,ak4exphi1_dR1p2=-99,ak4exe1_dR1p2=-99;
        ak4expt2_dR1p2=-99,ak4exeta2_dR1p2=-99,ak4exphi2_dR1p2=-99,ak4exe2_dR1p2=-99;

            double dR0p8=0.8;
        
        for(Int_t ii=0; ii<8; ii++)  {
        pt_AK4_ex_DR0p8[ii] =-99; eta_AK4_ex_DR0p8[ii] =-99.;phi_AK4_ex_DR0p8[ii] =-99.;e_AK4_ex_DR0p8[ii] =-99.;
        pt_AK4_ex_DR0p8_JEC_up[ii] =-99;pt_AK4_ex_DR0p8_JEC_down[ii] =-99;pt_AK4_ex_DR0p8_JER_up[ii] =-99;pt_AK4_ex_DR0p8_JER_down[ii] =-99;
        }
        
        for(Int_t ii=0; ii<8; ii++)  {
            if(ak4jet_pt[ii]>30 && fabs(ak4jet_eta[ii])<2.4 && ak4jet_IDLoose[ii]>0){
                deltaRjet1ak4p[ii]=sqrt(pow(fabs(ak4jet_eta[ii]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[ii]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[ii]-jetAK8puppi_phi)),2));
                if(Nj8>=2)
                deltaRjet2ak4p[ii]=sqrt(pow(fabs(ak4jet_eta[ii]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[ii]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[ii]-jetAK8puppi_phi_2)),2));
                if(Nj8>=3)
                deltaRjet3ak4p[ii]=sqrt(pow(fabs(ak4jet_eta[ii]-jetAK8puppi_eta_3),2)+pow(TMath::Min(fabs(ak4jet_phi[ii]-jetAK8puppi_phi_3),2*Pi-fabs(ak4jet_phi[ii]-jetAK8puppi_phi_3)),2));
                

            if((Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]<=dR0p8)||(Nj8==2&&(ak4jet_pt[ii]>0)&&(deltaRjet1ak4p[ii]<=dR0p8||deltaRjet2ak4p[ii]<=dR0p8))||(Nj8>=3&&(ak4jet_pt[ii]>0)&&(deltaRjet1ak4p[ii]<=dR0p8||deltaRjet2ak4p[ii]<=dR0p8||deltaRjet3ak4p[ii]<=dR0p8))) {
                ++Nj4in_dR0p8;
                if(Nj4in_dR0p8==1&&(ak4jet_pt[ii]>0)) {
                    //ak4expt1_dR0p8=ak4jet_pt[ii];ak4exeta1_dR0p8=ak4jet_eta[ii];ak4exphi1_dR0p8=ak4jet_phi[ii];ak4exe1_dR0p8=ak4jet_e[ii];ak4expt1_uncorr_dR0p8=ak4jet_pt_uncorr[ii];
                    ak4inpt1_dR0p8=ak4jet_pt[ii];
                    deltaRjet1inak4p_dR0p8=deltaRjet1ak4p[ii];
                    deltaphijet1inak4p_dR0p8=deltaPhi(ak4jet_phi[ii],jetAK8puppi_phi);
                    deltaRinak4mu_dR0p8=deltaR(l_eta,l_phi,ak4jet_eta[ii],ak4jet_phi[ii]);
                    deltaPhiinak4mu_dR0p8=deltaPhi(l_phi,ak4jet_phi[ii]);

            }

            }
            if((Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR0p8)||(Nj8==2&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR0p8&&deltaRjet2ak4p[ii]>dR0p8)||(Nj8>=3&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR0p8&&deltaRjet2ak4p[ii]>dR0p8&&deltaRjet3ak4p[ii]>dR0p8))
            {++Nj4ex_dR0p8;
            
            pt_AK4_ex_DR0p8[Nj4ex_dR0p8]=ak4jet_pt[ii];eta_AK4_ex_DR0p8[Nj4ex_dR0p8] =ak4jet_eta[ii];phi_AK4_ex_DR0p8[Nj4ex_dR0p8] =ak4jet_phi[ii];e_AK4_ex_DR0p8[Nj4ex_dR0p8] =ak4jet_e[ii];
            pt_AK4_ex_DR0p8_JEC_up[Nj4ex_dR0p8]=ak4jet_pt_JEC_up[ii];
            pt_AK4_ex_DR0p8_JEC_down[Nj4ex_dR0p8]=ak4jet_pt_JEC_down[ii];
            pt_AK4_ex_DR0p8_JER_up[Nj4ex_dR0p8]=ak4jet_pt_JER_up[ii];
            pt_AK4_ex_DR0p8_JER_down[Nj4ex_dR0p8]=ak4jet_pt_JER_down[ii];
            if(Nj4ex_dR0p8==1&&(ak4jet_pt[ii]>0)) {
                    ak4expt1_dR0p8=ak4jet_pt[ii];ak4exeta1_dR0p8=ak4jet_eta[ii];ak4exphi1_dR0p8=ak4jet_phi[ii];ak4exe1_dR0p8=ak4jet_e[ii];ak4expt1_uncorr_dR0p8=ak4jet_pt_uncorr[ii];
                deltaRjet1ak4p_dR0p8=deltaRjet1ak4p[ii];
                deltaphijet1ak4p_dR0p8=deltaPhi(ak4exphi1_dR0p8,jetAK8puppi_phi);
                deltaRak4mu_dR0p8=deltaR(l_eta,l_phi,ak4exeta1_dR0p8,ak4exphi1_dR0p8);
                //cout<<" deltaRjet1ak4  "<<deltaR(jetAK8puppi_eta,jetAK8puppi_phi,ak4jet_eta[ii],ak4jet_phi[ii])<<"  deltaRak4lep  "<<deltaR(l_eta,l_phi,ak4jet_eta[ii],ak4jet_phi[ii])<<"  deltaRjet1lep    "<<deltaR(jetAK8puppi_eta,jetAK8puppi_phi,l_eta,l_phi)<<endl;
                /*if(ak4expt1_uncorr_dR0p8>100){
                    cout<<"deltaRjet1ak4p_dR0p8 "<<deltaRjet1ak4p_dR0p8<<"  "<<jet_pt_puppi<<endl;
                for(int i=0;i<numsj1;i++)
                    cout<<"deltaRak4sj "<<i<<" "<<deltaRak4sj[i][ii]<<endl;
                for(int i=0;i<numsj2;i++)
                    cout<<"deltaRak4sj2 "<<i<<" "<<deltaRak4sj2[i][ii]<<endl;
                
                }*/
                }
            if(Nj4ex_dR0p8==2&&(ak4jet_pt[ii]>0)) {
                    ak4expt2_dR0p8=ak4jet_pt[ii];ak4exeta2_dR0p8=ak4jet_eta[ii];ak4exphi2_dR0p8=ak4jet_phi[ii];ak4exe2_dR0p8=ak4jet_e[ii];
                deltaRjet2ak4p_dR0p8=deltaRjet2ak4p[ii];

                }
            }
            double dR1p0=1.0;
            
            if(Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR1p0) ++Nj4ex_dR1p0;
            if(Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]<dR1p0) ++Nj4in_dR1p0;
            if(Nj8==2&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR1p0&&deltaRjet2ak4p[ii]>dR1p0) ++Nj4ex_dR1p0;
            if(Nj8==2&&(ak4jet_pt[ii]>0)&&(deltaRjet1ak4p[ii]<dR1p0||deltaRjet2ak4p[ii]<dR1p0)) ++Nj4in_dR1p0;
            
            if(Nj4ex_dR1p0==1&&(ak4jet_pt[ii]>0)) {
                ak4expt1_dR1p0=ak4jet_pt[ii];ak4exeta1_dR1p0=ak4jet_eta[ii];ak4exphi1_dR1p0=ak4jet_phi[ii];ak4exe1_dR1p0=ak4jet_e[ii];
            }
            if(Nj4ex_dR1p0==2&&(ak4jet_pt[ii]>0)) {
                ak4expt2_dR1p0=ak4jet_pt[ii];ak4exeta2_dR1p0=ak4jet_eta[ii];ak4exphi2_dR1p0=ak4jet_phi[ii];ak4exe2_dR1p0=ak4jet_e[ii];
            }

            double dR1p2=1.2;
            
            if(Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR1p2) ++Nj4ex_dR1p2;
            if(Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]<dR1p2) ++Nj4in_dR1p2;
            if(Nj8==2&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]<dR1p2&&deltaRjet2ak4p[ii]<dR1p2) ++Nj4ex_dR1p2;
            if(Nj8==2&&(ak4jet_pt[ii]>0)&&(deltaRjet1ak4p[ii]>dR1p2||deltaRjet2ak4p[ii]>dR1p2)) ++Nj4in_dR1p2;
            
            if(Nj4ex_dR1p2==1&&(ak4jet_pt[ii]>0)) {
                ak4expt1_dR1p2=ak4jet_pt[ii];ak4exeta1_dR1p2=ak4jet_eta[ii];ak4exphi1_dR1p2=ak4jet_phi[ii];ak4exe1_dR1p2=ak4jet_e[ii];
            }
            if(Nj4ex_dR1p2==2&&(ak4jet_pt[ii]>0)) {
                ak4expt2_dR1p2=ak4jet_pt[ii];ak4exeta2_dR1p2=ak4jet_eta[ii];ak4exphi2_dR1p2=ak4jet_phi[ii];ak4exe2_dR1p2=ak4jet_e[ii];
            }


        }
        //jet_mass_puppi_caseb = ak8jet1_2.Mag();
        }
}

#endif
