#ifndef _match_particle_for_each_jet_
#define _match_particle_for_each_jet_

void EDBR2PKUTree::match_particle_for_each_jet(Int_t IsData) {   
           for(int i=0;i<3;i++){
            matchedt[i]=-99;matchedw[i]=-99;unmatched[i]=-99;matchedz[i]=-99;matchedg[i]=-99;matchedr[i]=-99;
            matchedft[i]=-99;matchedfw[i]=-99;unmatchedf[i]=-99;
            matched0p4t[i]=-99;matched0p4w[i]=-99;unmatched0p4[i]=-99;matched0p4z[i]=-99;matched0p4g[i]=-99;matched0p4r[i]=-99;
            matched0p4ft[i]=-99;matched0p4fw[i]=-99;unmatched0p4f[i]=-99;
            matched0p8t[i]=-99;matched0p8w[i]=-99;unmatched0p8[i]=-99;matched0p8z[i]=-99;matched0p8g[i]=-99;matched0p8r[i]=-99;
            matched0p8ft[i]=-99;matched0p8fw[i]=-99;unmatched0p8f[i]=-99;}
        double deltaRtj[3]={99,99,99},deltaRantitj[3]={99,99,99},deltaRwj[5][3],deltaRwfj[5][3],deltaRzj[5][3],deltaRgj[10][3],deltaRzfj[5][3],deltaRgfj[10][3],deltaRrj[3]={99,99,99};
        if (gentop_pt>0){
            deltaRtj[0] = deltaR(gentop_eta,gentop_phi,eta_max,phi_max);
            deltaRtj[1] = deltaR(gentop_eta,gentop_phi,eta_mid,phi_mid);
            deltaRtj[2] = deltaR(gentop_eta,gentop_phi,eta_min,phi_min);}
        if(IsData==10&&gen_rad_pt>0){
            deltaRrj[0] = deltaR(gen_rad_eta,gen_rad_phi,eta_max,phi_max);
            deltaRrj[1] = deltaR(gen_rad_eta,gen_rad_phi,eta_mid,phi_mid);
            deltaRrj[2] = deltaR(gen_rad_eta,gen_rad_phi,eta_min,phi_min);}
        if (genantitop_pt>0){
            deltaRantitj[0] = deltaR(genantitop_eta,genantitop_phi,eta_max,phi_max);
            deltaRantitj[1] = deltaR(genantitop_eta,genantitop_phi,eta_mid,phi_mid);
            deltaRantitj[2] = deltaR(genantitop_eta,genantitop_phi,eta_min,phi_min);}
        for(int i=0;i<5;i++){
            for (int j=0;j<3;j++){
                deltaRwj[i][j]=99;
                deltaRwfj[i][j]=99;
                deltaRzj[i][j]=99;
                deltaRzfj[i][j]=99;
                }
            if(ptgenwl[i]>0){
                deltaRwj[i][0] = deltaR(etagenwl[i],phigenwl[i],eta_max,phi_max);
                deltaRwj[i][1] = deltaR(etagenwl[i],phigenwl[i],eta_mid,phi_mid);
                deltaRwj[i][2] = deltaR(etagenwl[i],phigenwl[i],eta_min,phi_min);}
            if(ptgenwf[i]>0){
                deltaRwfj[i][0] = deltaR(etagenwf[i],phigenwf[i],eta_max,phi_max);
                deltaRwfj[i][1] = deltaR(etagenwf[i],phigenwf[i],eta_mid,phi_mid);
                deltaRwfj[i][2] = deltaR(etagenwf[i],phigenwf[i],eta_min,phi_min);}
            if(ptgenzl[i]>0){
                deltaRzj[i][0] = deltaR(etagenzl[i],phigenzl[i],eta_max,phi_max);
                deltaRzj[i][1] = deltaR(etagenzl[i],phigenzl[i],eta_mid,phi_mid);
                deltaRzj[i][2] = deltaR(etagenzl[i],phigenzl[i],eta_min,phi_min);}
        }
        for(int i=0;i<10;i++){
            for (int j=0;j<3;j++){
                deltaRgj[i][j]=99;
                deltaRgfj[i][j]=99;}
            if(ptgengl[i]>0){
                deltaRgj[i][0] = deltaR(etagengl[i],phigengl[i],eta_max,phi_max);
                deltaRgj[i][1] = deltaR(etagengl[i],phigengl[i],eta_mid,phi_mid);
                deltaRgj[i][2] = deltaR(etagengl[i],phigengl[i],eta_min,phi_min);}
        }
        
        double deltaRwj1[5],deltaRwj2[5],deltaRwj3[5];
        for(int i=0;i<5;i++){
            deltaRwj1[i]=deltaRwj[i][0];
            deltaRwj2[i]=deltaRwj[i][1];
            deltaRwj3[i]=deltaRwj[i][2];
            if(Nj8==2) deltaRwj2[i]=deltaRwj[i][2];
        }
        Int_t *mindrwj1=new Int_t[5];
        TMath::Sort(5,deltaRwj1,mindrwj1,0);
        Int_t *mindrwj2=new Int_t[5];
        TMath::Sort(5,deltaRwj2,mindrwj2,0);
        Int_t *mindrwj3=new Int_t[5];
        TMath::Sort(5,deltaRwj3,mindrwj3,0);
        
        double rmatched=0.6;
        double masstorw=100;//110->100
        if(Nj8==1){
        if (deltaRtj[0]<rmatched||deltaRantitj[0]<rmatched) {matchedt[0]=1;matchedft[0]=1;}
        if (deltaRrj[0]<rmatched) {matchedr[0]=1;}
        for (int i=0;i<5;i++){
            if (deltaRwj[i][0]<rmatched) matchedw[0]=1;
            if (deltaRwfj[i][0]<rmatched) matchedfw[0]=1;
            if (deltaRzj[i][0]<rmatched) matchedz[0]=1;
        }
        for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched) matchedg[0]=1;}

        if(matchedt[0]<=0&&matchedr[0]<=0&&matchedw[0]<=0 && matchedz[0]<=0&&matchedg[0]<=0)  unmatched[0]=1;
            matchedtorg[0]=matchedt[0];matchedrorg[0]=matchedr[0];
            if(IsData>0 && IsData!=10&&jet_mass_puppi_corr<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matchedt[0]==1) matchedt[0]=0;
            if(IsData==10&&jet_mass_puppi_corr<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matchedr[0]==1) matchedr[0]=0;
        if(matchedft[0]<=0&&matchedfw[0]<=0)  unmatchedf[0]=1;
        }
        if(Nj8==2){
            if (deltaRtj[0]<rmatched||deltaRantitj[0]<rmatched) {matchedt[0]=1;matchedft[0]=1;}
            if (deltaRrj[0]<rmatched) {matchedr[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched) matchedw[0]=1;
                if (deltaRwfj[i][0]<rmatched) matchedfw[0]=1;
                if (deltaRzj[i][0]<rmatched) matchedz[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched) matchedg[0]=1;}
            if(matchedt[0]<=0&&matchedr[0]<=0&&matchedw[0]<=0 && matchedz[0]<=0&&matchedg[0]<=0)  unmatched[0]=1;
            matchedtorg[0]=matchedt[0];matchedrorg[0]=matchedr[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matchedt[0]==1) matchedt[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matchedr[0]==1) matchedr[0]=0;

            if(matchedft[0]<=0&&matchedfw[0]<=0)  unmatchedf[0]=1;
            
            if (deltaRtj[2]<rmatched||deltaRantitj[2]<rmatched) {matchedt[1]=1;matchedft[1]=1;}
            if (deltaRrj[2]<rmatched) {matchedr[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched) matchedw[1]=1;
                if (deltaRwfj[i][2]<rmatched) matchedfw[1]=1;
                if (deltaRzj[i][2]<rmatched) matchedz[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched) matchedg[2]=1;}

            if(matchedt[1]<=0&&matchedr[1]<=0&&matchedw[1]<=0 && matchedz[1]<=0&&matchedg[1]<=0)  unmatched[1]=1;
            matchedtorg[1]=matchedt[1];matchedrorg[1]=matchedr[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj2[mindrwj2[0]]&&matchedt[1]==1) matchedt[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[2]> deltaRwj2[mindrwj2[0]]&&matchedr[1]==1) matchedr[1]=0;

            if(matchedft[1]<=0&&matchedfw[1]<=0)  unmatchedf[1]=1;
        }
        if(Nj8==3){
            if (deltaRtj[0]<rmatched||deltaRantitj[0]<rmatched) {matchedt[0]=1;matchedft[0]=1;}
            if (deltaRrj[0]<rmatched) {matchedr[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched) matchedw[0]=1;
                if (deltaRwfj[i][0]<rmatched) matchedfw[0]=1;
                if (deltaRzj[i][0]<rmatched) matchedz[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched) matchedg[0]=1;}
            if(matchedt[0]<=0&&matchedr[0]<=0&&matchedw[0]<=0 && matchedz[0]<=0&&matchedg[0]<=0)  unmatched[0]=1;
            matchedtorg[0]=matchedt[0];matchedrorg[0]=matchedr[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matchedt[0]==1) matchedt[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matchedr[0]==1) matchedr[0]=0;
            if(matchedft[0]<=0&&matchedfw[0]<=0)  unmatchedf[0]=1;
            
            if (deltaRtj[1]<rmatched||deltaRantitj[1]<rmatched) {matchedt[1]=1;matchedft[1]=1;}
            if (deltaRrj[1]<rmatched) {matchedr[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][1]<rmatched) matchedw[1]=1;
                if (deltaRwfj[i][1]<rmatched) matchedfw[1]=1;
                if (deltaRzj[i][1]<rmatched) matchedz[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][1]<rmatched) matchedg[1]=1;}
            if(matchedt[1]<=0&&matchedr[1]<=0&&matchedw[1]<=0 && matchedz[1]<=0&&matchedg[1]<=0)  unmatched[1]=1;
            matchedtorg[1]=matchedt[1];matchedrorg[1]=matchedr[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[1],deltaRantitj[1])> deltaRwj2[mindrwj2[0]]&&matchedt[1]==1) matchedt[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[1]> deltaRwj2[mindrwj2[0]]&&matchedr[1]==1) matchedr[1]=0;
            if(matchedft[1]<=0&&matchedfw[1]<=0)  unmatchedf[1]=1;
            
            if (deltaRtj[2]<rmatched||deltaRantitj[2]<rmatched) {matchedt[2]=1;matchedft[2]=1;}
            if (deltaRrj[2]<rmatched) {matchedr[2]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched) matchedw[2]=1;
                if (deltaRwfj[i][2]<rmatched) matchedfw[2]=1;
                if (deltaRzj[i][2]<rmatched) matchedz[2]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[2][0]<rmatched) matchedg[2]=1;}
            if(matchedt[2]<=0&&matchedr[2]<=0&&matchedw[2]<=0 && matchedz[2]<=0&&matchedg[2]<=0)  unmatched[2]=1;
            matchedtorg[2]=matchedt[2];matchedrorg[2]=matchedr[2];
            if(IsData>0 && IsData!=10&&Mj_min<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj3[mindrwj3[0]]&&matchedt[2]==1) matchedt[2]=0;
            if(IsData==10&&Mj_min<masstorw && deltaRrj[2]> deltaRwj3[mindrwj3[0]]&&matchedr[2]==1) matchedr[2]=0;
            if(matchedft[2]<=0&&matchedfw[2]<=0)  unmatchedf[2]=1;
        }
        
        double rmatched0p4=0.4;
        if(Nj8==1){
            if (deltaRtj[0]<rmatched0p4||deltaRantitj[0]<rmatched0p4) {matched0p4t[0]=1;matched0p4ft[0]=1;}
            if (deltaRrj[0]<rmatched0p4) {matched0p4r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p4) matched0p4w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p4) matched0p4fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p4) matched0p4z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p4) matched0p4g[0]=1;}

            if(matched0p4t[0]<=0&&matched0p4r[0]<=0&&matched0p4w[0]<=0 && matched0p4z[0]<=0&&matched0p4g[0]<=0)  unmatched0p4[0]=1;
            matched0p4torg[0]=matched0p4t[0];matched0p4rorg[0]=matched0p4r[0];
            if(IsData>0 && IsData!=10&&jet_mass_puppi_corr<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p4t[0]==1) matched0p4t[0]=0;
            if(IsData==10&&jet_mass_puppi_corr<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p4r[0]==1) matched0p4r[0]=0;
            if(matched0p4ft[0]<=0&&matched0p4fw[0]<=0)  unmatched0p4f[0]=1;
        }
        if(Nj8==2){
            if (deltaRtj[0]<rmatched0p4||deltaRantitj[0]<rmatched0p4) {matched0p4t[0]=1;matched0p4ft[0]=1;}
            if (deltaRrj[0]<rmatched0p4) {matched0p4r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p4) matched0p4w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p4) matched0p4fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p4) matched0p4z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p4) matched0p4g[0]=1;}

            if(matched0p4t[0]<=0&&matched0p4r[0]<=0&&matched0p4w[0]<=0 && matched0p4z[0]<=0&&matched0p4g[0]<=0)  unmatched0p4[0]=1;
            matched0p4torg[0]=matched0p4t[0];matched0p4rorg[0]=matched0p4r[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p4t[0]==1) matched0p4t[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p4r[0]==1) matched0p4r[0]=0;
            
            if(matched0p4ft[0]<=0&&matched0p4fw[0]<=0)  unmatched0p4f[0]=1;
            
            if (deltaRtj[2]<rmatched0p4||deltaRantitj[2]<rmatched0p4) {matched0p4t[1]=1;matched0p4ft[1]=1;}
            if (deltaRrj[2]<rmatched0p4) {matched0p4r[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched0p4) matched0p4w[1]=1;
                if (deltaRwfj[i][2]<rmatched0p4) matched0p4fw[1]=1;
                if (deltaRzj[i][2]<rmatched0p4) matched0p4z[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched0p4) matched0p4g[1]=1;}

            if(matched0p4t[1]<=0&&matched0p4r[1]<=0&&matched0p4w[1]<=0 && matched0p4z[1]<=0&&matched0p4g[1]<=0)  unmatched0p4[1]=1;
            matched0p4torg[1]=matched0p4t[1];matched0p4rorg[1]=matched0p4r[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj2[mindrwj2[0]]&&matched0p4t[1]==1) matched0p4t[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[2]> deltaRwj2[mindrwj2[0]]&&matched0p4r[1]==1) matched0p4r[1]=0;
            
            if(matched0p4ft[1]<=0&&matched0p4fw[1]<=0)  unmatched0p4f[1]=1;
        }
        if(Nj8==3){
            if (deltaRtj[0]<rmatched0p4||deltaRantitj[0]<rmatched0p4) {matched0p4t[0]=1;matched0p4ft[0]=1;}
            if (deltaRrj[0]<rmatched0p4) {matched0p4r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p4) matched0p4w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p4) matched0p4fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p4) matched0p4z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p4) matched0p4g[0]=1;}

            if(matched0p4t[0]<=0&&matched0p4r[0]<=0&&matched0p4w[0]<=0 && matched0p4z[0]<=0&&matched0p4g[0]<=0)  unmatched0p4[0]=1;
            matched0p4torg[0]=matched0p4t[0];matched0p4rorg[0]=matched0p4r[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p4t[0]==1) matched0p4t[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p4r[0]==1) matched0p4r[0]=0;
            if(matched0p4ft[0]<=0&&matched0p4fw[0]<=0)  unmatched0p4f[0]=1;
            
            if (deltaRtj[1]<rmatched0p4||deltaRantitj[1]<rmatched0p4) {matched0p4t[1]=1;matched0p4ft[1]=1;}
            if (deltaRrj[1]<rmatched0p4) {matched0p4r[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][1]<rmatched0p4) matched0p4w[1]=1;
                if (deltaRwfj[i][1]<rmatched0p4) matched0p4fw[1]=1;
                if (deltaRzj[i][1]<rmatched0p4) matched0p4z[1]=1;
                if (deltaRgj[i][1]<rmatched0p4) matched0p4g[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][1]<rmatched0p4) matched0p4g[1]=1;}

            if(matched0p4t[1]<=0&&matched0p4r[1]<=0&&matched0p4w[1]<=0 && matched0p4z[1]<=0&&matched0p4g[1]<=0)  unmatched0p4[1]=1;
            matched0p4torg[1]=matched0p4t[1];matched0p4rorg[1]=matched0p4r[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[1],deltaRantitj[1])> deltaRwj2[mindrwj2[0]]&&matched0p4t[1]==1) matched0p4t[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[1]> deltaRwj2[mindrwj2[0]]&&matched0p4r[1]==1) matched0p4r[1]=0;
            if(matched0p4ft[1]<=0&&matched0p4fw[1]<=0)  unmatched0p4f[1]=1;
            
            if (deltaRtj[2]<rmatched0p4||deltaRantitj[2]<rmatched0p4) {matched0p4t[2]=1;matched0p4ft[2]=1;}
            if (deltaRrj[2]<rmatched0p4) {matched0p4r[2]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched0p4) matched0p4w[2]=1;
                if (deltaRwfj[i][2]<rmatched0p4) matched0p4fw[2]=1;
                if (deltaRzj[i][2]<rmatched0p4) matched0p4z[2]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched0p4) matched0p4g[2]=1;}

            if(matched0p4t[2]<=0&&matched0p4r[2]<=0&&matched0p4w[2]<=0 && matched0p4z[2]<=0&&matched0p4g[2]<=0)  unmatched0p4[2]=1;
            matched0p4torg[2]=matched0p4t[2];matched0p4rorg[2]=matched0p4r[2];
            if(IsData>0 && IsData!=10&&Mj_min<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj3[mindrwj3[0]]&&matched0p4t[2]==1) matched0p4t[2]=0;
            if(IsData==10&&Mj_min<masstorw && deltaRrj[2]> deltaRwj3[mindrwj3[0]]&&matched0p4r[2]==1) matched0p4r[2]=0;
            if(matched0p4ft[2]<=0&&matched0p4fw[2]<=0)  unmatched0p4f[2]=1;
        }

        double rmatched0p8=0.8;
        if(Nj8==1){
            if (deltaRtj[0]<rmatched0p8||deltaRantitj[0]<rmatched0p8) {matched0p8t[0]=1;matched0p8ft[0]=1;}
            if (deltaRrj[0]<rmatched0p8) {matched0p8r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p8) matched0p8w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p8) matched0p8fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p8) matched0p8z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p8) matched0p8g[0]=1;}

            if(matched0p8t[0]<=0&&matched0p8r[0]<=0&&matched0p8w[0]<=0 && matched0p8z[0]<=0&&matched0p8g[0]<=0)  unmatched0p8[0]=1;
            matched0p8torg[0]=matched0p8t[0];matched0p8rorg[0]=matched0p8r[0];
            if(IsData>0 && IsData!=10&&jet_mass_puppi_corr<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p8t[0]==1) matched0p8t[0]=0;
            if(IsData==10&&jet_mass_puppi_corr<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p8r[0]==1) matched0p8r[0]=0;
            if(matched0p8ft[0]<=0&&matched0p8fw[0]<=0)  unmatched0p8f[0]=1;
        }
        if(Nj8==2){
            if (deltaRtj[0]<rmatched0p8||deltaRantitj[0]<rmatched0p8) {matched0p8t[0]=1;matched0p8ft[0]=1;}
            if (deltaRrj[0]<rmatched0p8) {matched0p8r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p8) matched0p8w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p8) matched0p8fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p8) matched0p8z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p8) matched0p8g[0]=1;}

            if(matched0p8t[0]<=0&&matched0p8r[0]<=0&&matched0p8w[0]<=0 && matched0p8z[0]<=0&&matched0p8g[0]<=0)  unmatched0p8[0]=1;
            matched0p8torg[0]=matched0p8t[0];matched0p8rorg[0]=matched0p8r[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p8t[0]==1) matched0p8t[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p8r[0]==1) matched0p8r[0]=0;
            
            if(matched0p8ft[0]<=0&&matched0p8fw[0]<=0)  unmatched0p8f[0]=1;
            
            if (deltaRtj[2]<rmatched0p8||deltaRantitj[2]<rmatched0p8) {matched0p8t[1]=1;matched0p8ft[1]=1;}
            if (deltaRrj[2]<rmatched0p8) {matched0p8r[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched0p8) matched0p8w[1]=1;
                if (deltaRwfj[i][2]<rmatched0p8) matched0p8fw[1]=1;
                if (deltaRzj[i][2]<rmatched0p8) matched0p8z[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched0p8) matched0p8g[1]=1;}

            if(matched0p8t[1]<=0&&matched0p8r[1]<=0&&matched0p8w[1]<=0 && matched0p8z[1]<=0&&matched0p8g[1]<=0)  unmatched0p8[1]=1;
            matched0p8torg[1]=matched0p8t[1];matched0p8rorg[1]=matched0p8r[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj2[mindrwj2[0]]&&matched0p8t[1]==1) matched0p8t[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[2]> deltaRwj2[mindrwj2[0]]&&matched0p8r[1]==1) matched0p8r[1]=0;
            
            if(matched0p8ft[1]<=0&&matched0p8fw[1]<=0)  unmatched0p8f[1]=1;
        }
        if(Nj8==3){
            if (deltaRtj[0]<rmatched0p8||deltaRantitj[0]<rmatched0p8) {matched0p8t[0]=1;matched0p8ft[0]=1;}
            if (deltaRrj[0]<rmatched0p8) {matched0p8r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p8) matched0p8w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p8) matched0p8fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p8) matched0p8z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p8) matched0p8g[0]=1;}

            if(matched0p8t[0]<=0&&matched0p8r[0]<=0&&matched0p8w[0]<=0 && matched0p8z[0]<=0&&matched0p8g[0]<=0)  unmatched0p8[0]=1;
            matched0p8torg[0]=matched0p8t[0];matched0p8rorg[0]=matched0p8r[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p8t[0]==1) matched0p8t[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p8r[0]==1) matched0p8r[0]=0;
            if(matched0p8ft[0]<=0&&matched0p8fw[0]<=0)  unmatched0p8f[0]=1;
            
            if (deltaRtj[1]<rmatched0p8||deltaRantitj[1]<rmatched0p8) {matched0p8t[1]=1;matched0p8ft[1]=1;}
            if (deltaRrj[1]<rmatched0p8) {matched0p8r[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][1]<rmatched0p8) matched0p8w[1]=1;
                if (deltaRwfj[i][1]<rmatched0p8) matched0p8fw[1]=1;
                if (deltaRzj[i][1]<rmatched0p8) matched0p8z[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][1]<rmatched0p8) matched0p8g[1]=1;}

            if(matched0p8t[1]<=0&&matched0p8r[1]<=0&&matched0p8w[1]<=0 && matched0p8z[1]<=0&&matched0p8g[1]<=0)  unmatched0p8[1]=1;
            matched0p8torg[1]=matched0p8t[1];matched0p8rorg[1]=matched0p8r[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[1],deltaRantitj[1])> deltaRwj2[mindrwj2[0]]&&matched0p8t[1]==1) matched0p8t[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[1]> deltaRwj2[mindrwj2[0]]&&matched0p8r[1]==1) matched0p8r[1]=0;
            if(matched0p8ft[1]<=0&&matched0p8fw[1]<=0)  unmatched0p8f[1]=1;
            
            if (deltaRtj[2]<rmatched0p8||deltaRantitj[2]<rmatched0p8) {matched0p8t[2]=1;matched0p8ft[2]=1;}
            if (deltaRrj[2]<rmatched0p8) {matched0p8r[2]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched0p8) matched0p8w[2]=1;
                if (deltaRwfj[i][2]<rmatched0p8) matched0p8fw[2]=1;
                if (deltaRzj[i][2]<rmatched0p8) matched0p8z[2]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched0p8) matched0p8g[2]=1;}

            if(matched0p8t[2]<=0&&matched0p8r[2]<=0&&matched0p8w[2]<=0 && matched0p8z[2]<=0&&matched0p8g[2]<=0)  unmatched0p8[2]=1;
            matched0p8torg[2]=matched0p8t[2];matched0p8rorg[2]=matched0p8r[2];
            if(IsData>0 && IsData!=10&&Mj_min<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj3[mindrwj3[0]]&&matched0p8t[2]==1) matched0p8t[2]=0;
            if(IsData==10&&Mj_min<masstorw && deltaRrj[2]> deltaRwj3[mindrwj3[0]]&&matched0p8r[2]==1) matched0p8r[2]=0;
            if(matched0p8ft[2]<=0&&matched0p8fw[2]<=0)  unmatched0p8f[2]=1;
        }
        
        
// gen top matching reco leptonic decay t (lepw and bjet)       
        TLorentzVector wlep,bjet;
        wlep.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
        mtbtvlep=-99;        mtbmvlep=-99;deltaRbtvlep=-99;deltaRbmvlep=-99;deltaRblvlep=-99;mtblvlep=-99;mbmvlep=-99;
        deltaRantitop_wlepbm=99;deltaRtop_wlepbm=99;deltaRantitop_wlepbl=99;deltaRtop_wlepbl=99;mblvlep=-99;
        if(num_bJet_tight_deep==1){
        for(Int_t i=0; i<8; i++)  {
            deltaRAK4AK8_new[i]=0.;
            deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
            if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));
            if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.8958 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {bjet.SetPtEtaPhiE(ak4jet_pt[i],ak4jet_eta[i],ak4jet_phi[i],ak4jet_e[i]);
                mtbtvlep=Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[i]*W_pt*(1.0-cos(ak4jet_phi[i]-W_phi))));
                deltaRbtvlep=deltaR(W_eta,W_phi,ak4jet_eta[i],ak4jet_phi[i]);
                //cout<<"hhhh "<<massVlepJEC<<"   "<<bjet.Mag()<<"   "<<(wlep+bjet).Mt()<<"   "<<mtbtvlep<<"   "<<Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[i]*W_pt*(1.0-cos(ak4jet_phi[i]-W_phi))))<<"   "<<deltaRbtvlep<<endl;
            }
            
        }}
        if(num_bJet_deep>=1){
            double wlepak4[8]={99,99,99,99,99,99,99,99};
            int iclosewlepak4=0;
            for(Int_t i=0; i<8; i++)  {
                deltaRAK4AK8_new[i]=0.;
                deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
                if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));
                if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.6324 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {
                    wlepak4[i]=deltaR(W_eta,W_phi,ak4jet_eta[i],ak4jet_phi[i]);
                    //cout<<"hhhhm "<<massVlepJEC<<"   "<<bjet.Mag()<<"   "<<(wlep+bjet).Mt()<<"   "<<mtbmvlep<<"   "<<Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[i]*W_pt*(1.0-cos(ak4jet_phi[i]-W_phi))))<<"   "<<deltaRbmvlep<<endl;
                }
            }
            Int_t *indexx_wlepak4=new Int_t[8];
            TMath::Sort(8,wlepak4,indexx_wlepak4,0);
            bjet.SetPtEtaPhiE(ak4jet_pt[indexx_wlepak4[0]],ak4jet_eta[indexx_wlepak4[0]],ak4jet_phi[indexx_wlepak4[0]],ak4jet_e[indexx_wlepak4[0]]);
            mtbmvlep=Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[indexx_wlepak4[0]]*W_pt*(1.0-cos(ak4jet_phi[indexx_wlepak4[0]]-W_phi))));
            deltaRbmvlep=deltaR(W_eta,W_phi,ak4jet_eta[indexx_wlepak4[0]],ak4jet_phi[indexx_wlepak4[0]]);
            TLorentzVector wlepak4v;wlepak4v=wlep+bjet;
            mbmvlep=wlepak4v.Mag();
            deltaRtop_wlepbm=99;deltaRantitop_wlepbm=99;
            if (gentop_pt>0) deltaRtop_wlepbm = deltaR(gentop_eta,gentop_phi,wlepak4v.Eta(),wlepak4v.Phi());
            if (genantitop_pt >0) deltaRantitop_wlepbm = deltaR(genantitop_eta,genantitop_phi,wlepak4v.Eta(),wlepak4v.Phi());
        }
        if(num_bJet_loose_deep>=1){
            double wlepak4[8]={99,99,99,99,99,99,99,99};
            int iclosewlepak4=0;
            for(Int_t i=0; i<8; i++)  {
                deltaRAK4AK8_new[i]=0.;
                deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
                if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));
                if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.2219 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {
                    wlepak4[i]=deltaR(W_eta,W_phi,ak4jet_eta[i],ak4jet_phi[i]);
                    //cout<<"hhhhm "<<massVlepJEC<<"   "<<bjet.Mag()<<"   "<<(wlep+bjet).Mt()<<"   "<<mtblvlep<<"   "<<Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[i]*W_pt*(1.0-cos(ak4jet_phi[i]-W_phi))))<<"   "<<deltaRblvlep<<endl;
                }
            }
            Int_t *indexx_wlepak4=new Int_t[8];
            TMath::Sort(8,wlepak4,indexx_wlepak4,0);
            bjet.SetPtEtaPhiE(ak4jet_pt[indexx_wlepak4[0]],ak4jet_eta[indexx_wlepak4[0]],ak4jet_phi[indexx_wlepak4[0]],ak4jet_e[indexx_wlepak4[0]]);
            mtblvlep=Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[indexx_wlepak4[0]]*W_pt*(1.0-cos(ak4jet_phi[indexx_wlepak4[0]]-W_phi))));
            deltaRblvlep=deltaR(W_eta,W_phi,ak4jet_eta[indexx_wlepak4[0]],ak4jet_phi[indexx_wlepak4[0]]);
            TLorentzVector wlepak4v;wlepak4v=wlep+bjet;
            mblvlep=wlepak4v.Mag();
            deltaRtop_wlepbl=99;deltaRantitop_wlepbl=99;
            if (gentop_pt>0) deltaRtop_wlepbl = deltaR(gentop_eta,gentop_phi,wlepak4v.Eta(),wlepak4v.Phi());
            if (genantitop_pt >0) deltaRantitop_wlepbl = deltaR(genantitop_eta,genantitop_phi,wlepak4v.Eta(),wlepak4v.Phi());
        }
// gen top matching reco leptonic decay t (lepw and bjet)    

// unknown part

        for(int i=0;i<3;i++){
            matchingt[i]=-99;matchingw[i]=-99;unmatching[i]=-99;matchingz[i]=-99;matchingg[i]=-99;matchingr[i]=-99;matchingq[i]=-99;
            for (int j=0;j<9;j++){
                matchingtg[j][i]=-99;matchingtq[j][i]=-99;}
        }
        double dRtj[3]={99,99,99},dRantitj[3]={99,99,99},dRwj[5][3],dRzj[5][3],dRgj[10][3],dRrj[3]={99,99,99},dRtg[10],dRq1j[5][3],dRq2j[5][3],dRq3j[5][3],dRq4j[5][3],dRq5j[5][3],dRtq1[5],dRtq2[5],dRtq3[5],dRtq4[5],dRtq5[5];
        double dRw_q1j[5][3],dRw_q2j[5][3],dRw2j[3]={99,99,99},dRw3j[3]={99,99,99},matchingtmp[3]={-99,-99,-99};
        if (gentop_pt>0){
            dRtj[0] = deltaR(gentop_eta,gentop_phi,eta_max,phi_max);
            if(Nj8==3)dRtj[1] = deltaR(gentop_eta,gentop_phi,eta_mid,phi_mid);
            if(Nj8==2)dRtj[1] = deltaR(gentop_eta,gentop_phi,eta_min,phi_min);
            dRtj[2] = deltaR(gentop_eta,gentop_phi,eta_min,phi_min);}
        if(IsData==10&&gen_rad_pt>0){
            dRrj[0] = deltaR(gen_rad_eta,gen_rad_phi,eta_max,phi_max);
            if(Nj8==3) dRrj[1] = deltaR(gen_rad_eta,gen_rad_phi,eta_mid,phi_mid);
            if(Nj8==2) dRrj[1] = deltaR(gen_rad_eta,gen_rad_phi,eta_min,phi_min);
            dRrj[2] = deltaR(gen_rad_eta,gen_rad_phi,eta_min,phi_min);}
        if (genantitop_pt>0){
            dRantitj[0] = deltaR(genantitop_eta,genantitop_phi,eta_max,phi_max);
            if(Nj8==3) dRantitj[1] = deltaR(genantitop_eta,genantitop_phi,eta_mid,phi_mid);
            if(Nj8==2) dRantitj[1] = deltaR(genantitop_eta,genantitop_phi,eta_min,phi_min);
            dRantitj[2] = deltaR(genantitop_eta,genantitop_phi,eta_min,phi_min);}
        for(int i=0;i<5;i++){
            dRtq1[i]=99;dRtq2[i]=99;dRtq3[i]=99;dRtq4[i]=99;dRtq5[i]=99;
            for (int j=0;j<3;j++){
                dRwj[i][j]=99;
                dRzj[i][j]=99;
                dRq1j[i][j]=99;dRq2j[i][j]=99;dRq3j[i][j]=99;dRq4j[i][j]=99;dRq5j[i][j]=99;
                dRw_q1j[i][j]=99;dRw_q2j[i][j]=99;
            }
            if(ptgenwl[i]>0){
                dRwj[i][0] = deltaR(etagenwl[i],phigenwl[i],eta_max,phi_max);
                if(Nj8==3) dRwj[i][1] = deltaR(etagenwl[i],phigenwl[i],eta_mid,phi_mid);
                if(Nj8==2) dRwj[i][1] = deltaR(etagenwl[i],phigenwl[i],eta_min,phi_min);
                dRwj[i][2] = deltaR(etagenwl[i],phigenwl[i],eta_min,phi_min);}
            if(ptgenzl[i]>0){
                dRzj[i][0] = deltaR(etagenzl[i],phigenzl[i],eta_max,phi_max);
                if(Nj8==3) dRzj[i][1] = deltaR(etagenzl[i],phigenzl[i],eta_mid,phi_mid);
                if(Nj8==2) dRzj[i][1] = deltaR(etagenzl[i],phigenzl[i],eta_min,phi_min);
                dRzj[i][2] = deltaR(etagenzl[i],phigenzl[i],eta_min,phi_min);}
            if(ptgenq1l[i]>0){
                dRq1j[i][0] = deltaR(etagenq1l[i],phigenq1l[i],eta_max,phi_max);
                if(Nj8==3) dRq1j[i][1] = deltaR(etagenq1l[i],phigenq1l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq1j[i][1] = deltaR(etagenq1l[i],phigenq1l[i],eta_min,phi_min);
                dRq1j[i][2] = deltaR(etagenq1l[i],phigenq1l[i],eta_min,phi_min);}
            if(ptgenq2l[i]>0){
                dRq2j[i][0] = deltaR(etagenq2l[i],phigenq2l[i],eta_max,phi_max);
                if(Nj8==3) dRq2j[i][1] = deltaR(etagenq2l[i],phigenq2l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq2j[i][1] = deltaR(etagenq2l[i],phigenq2l[i],eta_min,phi_min);
                dRq2j[i][2] = deltaR(etagenq2l[i],phigenq2l[i],eta_min,phi_min);}
            if(ptgenq3l[i]>0){
                dRq3j[i][0] = deltaR(etagenq3l[i],phigenq3l[i],eta_max,phi_max);
                if(Nj8==3) dRq3j[i][1] = deltaR(etagenq3l[i],phigenq3l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq3j[i][1] = deltaR(etagenq3l[i],phigenq3l[i],eta_min,phi_min);
                dRq3j[i][2] = deltaR(etagenq3l[i],phigenq3l[i],eta_min,phi_min);}
            if(ptgenq4l[i]>0){
                dRq4j[i][0] = deltaR(etagenq4l[i],phigenq4l[i],eta_max,phi_max);
                if(Nj8==3) dRq4j[i][1] = deltaR(etagenq4l[i],phigenq4l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq4j[i][1] = deltaR(etagenq4l[i],phigenq4l[i],eta_min,phi_min);
                dRq4j[i][2] = deltaR(etagenq4l[i],phigenq4l[i],eta_min,phi_min);}
            if(ptgenq5l[i]>0){
                dRq5j[i][0] = deltaR(etagenq5l[i],phigenq5l[i],eta_max,phi_max);
                if(Nj8==3) dRq5j[i][1] = deltaR(etagenq5l[i],phigenq5l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq5j[i][1] = deltaR(etagenq5l[i],phigenq5l[i],eta_min,phi_min);
                dRq5j[i][2] = deltaR(etagenq5l[i],phigenq5l[i],eta_min,phi_min);}

        }
        for(int i=0;i<10;i++){
            dRtg[i]=99;
            for (int j=0;j<3;j++){dRgj[i][j]=99;}
            if(ptgengl[i]>0){
                dRgj[i][0] = deltaR(etagengl[i],phigengl[i],eta_max,phi_max);
                if(Nj8==3) dRgj[i][1] = deltaR(etagengl[i],phigengl[i],eta_mid,phi_mid);
                if(Nj8==2) dRgj[i][1] = deltaR(etagengl[i],phigengl[i],eta_min,phi_min);
                dRgj[i][2] = deltaR(etagengl[i],phigengl[i],eta_min,phi_min);}
        }
        
        double dRt_bj[3]={99,99,99},dRt_wq1j[3]={99,99,99},dRt_wq2j[3]={99,99,99},dRt_wj[3]={99,99,99},dRt_wq1q2[3]={99,99,99};
        double tmpbpt=-99,tmpbm=-99,tmpbeta=-99,tmpbphi=-99,tmpq1pt=-99,tmpq1eta=-99,tmpq1phi=-99,tmpq2pt=-99,tmpq2eta=-99,tmpq2phi=-99,tmpwpt=-99,tmpweta=-99,tmpwphi=-99,tmptpt=-99,tmpteta=-99,tmptphi=-99,tmpt_wtag=-99;
        double dRr_q1j[3]={99,99,99},dRr_q2j[3]={99,99,99},dRr_q3j[3]={99,99,99},dRr_q4j[3]={99,99,99};
        double dRt1_bj[3]={99,99,99},dRt2_bj[3]={99,99,99},dRt1_wq1j[3]={99,99,99},dRt1_wq2j[3]={99,99,99},dRt2_wq1j[3]={99,99,99},dRt2_wq2j[3]={99,99,99},dRt_wb[3]={99,99,99},dRt1_wj[3]={99,99,99},dRt2_wj[3]={99,99,99},dRt1_wb[3]={99,99,99},dRt2_wb[3]={99,99,99};
        double rmatching=0.6,rbmatching=0.8;
        TLorentzVector t_b,matchedgv,matchedq1v,matchedq2v,matchedq3v,matchedq4v,matchedq5v;
        dRWW_R=-99;dRWb=-99;dRqq_W1=-99;dRqq_W2=-99;mindRqq=-99;maxdRqq=-99;
        double mindRqq_all=-99,maxdRqq_all=-99;
        //if (jentry==100000) break;
        //cout<<jentry<<endl;
        for(int i=0;i<Nj8;i++){
            //if (Nj8==2) break;
            if(gentop_pt>0||genantitop_pt>0){
                if (dRtj[i]<rmatching) matchingtmp[i]=-1;
                if (matchingtmp[i]!=-1&&dRantitj[i]<rmatching) matchingtmp[i]=-2;
                if (matchingtmp[i]==-1&&dRantitj[i]<rmatching&&dRtj[i]<dRantitj[i]) matchingtmp[i]=-1;
                if (matchingtmp[i]==-1&&dRantitj[i]<rmatching&&dRtj[i]>dRantitj[i]) matchingtmp[i]=-2;
                //cout<<i<<" jet  "<<Nj8<<"   "<<matchingtmp[i]<<endl;
                if(matchingtmp[i]==-1||matchingtmp[i]==-2){
                    if(matchingtmp[i]==-1){tmptpt=gentop_pt;tmpteta=gentop_eta;tmptphi=gentop_phi;tmpbpt=gent_b_pt;tmpbeta=gent_b_eta;tmpbphi=gent_b_phi;tmpwpt=gent_w_pt;tmpweta=gent_w_eta;tmpwphi=gent_w_phi;tmpq1pt=gent_w_q1_pt;tmpq1eta=gent_w_q1_eta;tmpq1phi=gent_w_q1_phi;tmpq2pt=gent_w_q2_pt;tmpq2eta=gent_w_q2_eta;tmpq2phi=gent_w_q2_phi;tmpbm=gent_b_mass;tmpt_wtag=gent_w_tag;}
                    if(matchingtmp[i]==-2){tmptpt=genantitop_pt;tmpteta=genantitop_eta;tmptphi=genantitop_phi;tmpbpt=genantit_b_pt;tmpbeta=genantit_b_eta;tmpbphi=genantit_b_phi;tmpwpt=genantit_w_pt;tmpweta=genantit_w_eta;tmpwphi=genantit_w_phi;tmpq1pt=genantit_w_q1_pt;tmpq1eta=genantit_w_q1_eta;tmpq1phi=genantit_w_q1_phi;tmpq2pt=genantit_w_q2_pt;tmpq2eta=genantit_w_q2_eta;tmpq2phi=genantit_w_q2_phi;tmpbm=genantit_b_mass;tmpt_wtag=genantit_w_tag;}
                        if(tmpbpt<0) matchingt[i]=-1;
                        if(tmpbpt>0){
                            if(i==0) dRt_bj[i]=deltaR(tmpbeta,tmpbphi,eta_max,phi_max);
                            if(i==1&&Nj8!=2) dRt_bj[i]=deltaR(tmpbeta,tmpbphi,eta_mid,phi_mid);
                            if(i==1&&Nj8==2) dRt_bj[i]=deltaR(tmpbeta,tmpbphi,eta_min,phi_min);
                            if(i==2) dRt_bj[i]=deltaR(tmpbeta,tmpbphi,eta_min,phi_min);
                            if(i==0) dRt_wq1j[i]=deltaR(tmpq1eta,tmpq1phi,eta_max,phi_max);
                            if(i==1&&Nj8!=2) dRt_wq1j[i]=deltaR(tmpq1eta,tmpq1phi,eta_mid,phi_mid);
                            if(i==1&&Nj8==2) dRt_wq1j[i]=deltaR(tmpq1eta,tmpq1phi,eta_min,phi_min);
                            if(i==2) dRt_wq1j[i]=deltaR(tmpq1eta,tmpq1phi,eta_min,phi_min);
                            if(i==0) dRt_wq2j[i]=deltaR(tmpq2eta,tmpq2phi,eta_max,phi_max);
                            if(i==1&&Nj8!=2) dRt_wq2j[i]=deltaR(tmpq2eta,tmpq2phi,eta_mid,phi_mid);
                            if(i==1&&Nj8==2) dRt_wq2j[i]=deltaR(tmpq2eta,tmpq2phi,eta_min,phi_min);
                            if(i==2) dRt_wq2j[i]=deltaR(tmpq2eta,tmpq2phi,eta_min,phi_min);
                            if(i==0) dRt_wj[i]=deltaR(tmpweta,tmpwphi,eta_max,phi_max);
                            if(i==1&&Nj8!=2) dRt_wj[i]=deltaR(tmpweta,tmpwphi,eta_mid,phi_mid);
                            if(i==1&&Nj8==2) dRt_wj[i]=deltaR(tmpweta,tmpwphi,eta_min,phi_min);
                            if(i==2) dRt_wj[i]=deltaR(tmpweta,tmpwphi,eta_min,phi_min);
                            dRt_wb[i]=deltaR(tmpweta,tmpwphi,tmpbeta,tmpbphi);
                            if(i==0&&Mj_max>100){
                                dRWb=deltaR(tmpweta,tmpwphi,tmpbeta,tmpbphi);
                                dRqq_W1=deltaR(tmpq1eta,tmpq1phi,tmpq2eta,tmpq2phi);
                            }
                            //if(dRt_bj[i]<rbmatching&&dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]<rbmatching&&dRt_wj[i]<rmatching&&tmpt_wtag==4)
                            if(dRt_bj[i]<rbmatching&&dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]<rbmatching&&tmpt_wtag==4){
                                matchingt[i]=3;
                                for(int j=0;j<10;j++){
                                    dRtg[j]=deltaR(tmpteta,tmptphi,etagengl[j],phigengl[j]);
                                    if(dRtg[j]<rbmatching&&dRgj[j][i]<rbmatching) matchingt[i]=4;
                                    t_b.SetPtEtaPhiM(tmpbpt,tmpbeta,tmpbphi,tmpbm);
                                    matchedgv.SetPtEtaPhiE(ptgengl[j],etagengl[j],phigengl[j],egengl[j]);
                                    if(matchingt[i]==4&&(t_b+matchedgv).Mag()>40) matchingtg[0][i]=4;
                                    if(matchingt[i]==4&&(t_b+matchedgv).Mag()>60) matchingtg[1][i]=4;
                                    if(matchingt[i]==4&&(t_b+matchedgv).Mag()>40&&(t_b+matchedgv).Mag()<120) matchingtg[2][i]=4;
                                    if(matchingt[i]==4&&(t_b+matchedgv).Mag()>60&&(t_b+matchedgv).Mag()<100) matchingtg[3][i]=4;
                                    if(matchingt[i]==4&&ptgengl[j]>0.1*tmptpt) matchingtg[4][i]=4;
                                    if(matchingt[i]==4&&ptgengl[j]>0.15*tmptpt) matchingtg[5][i]=4;
                                    if(matchingt[i]==4&&ptgengl[j]>0.2*tmptpt) matchingtg[6][i]=4;
                                    if(matchingt[i]==4&&ptgengl[j]>0.25*tmptpt) matchingtg[7][i]=4;
                                    if(matchingt[i]==4&&(t_b+matchedgv).Pt()>0.2*tmpwpt) matchingtg[8][i]=4;
                                }
                                for(int j=0;j<5;j++){
                                    dRtq1[j]=deltaR(tmpteta,tmptphi,etagenq1l[j],phigenq1l[j]);
                                    dRtq2[j]=deltaR(tmpteta,tmptphi,etagenq2l[j],phigenq2l[j]);
                                    dRtq3[j]=deltaR(tmpteta,tmptphi,etagenq3l[j],phigenq3l[j]);
                                    dRtq4[j]=deltaR(tmpteta,tmptphi,etagenq4l[j],phigenq4l[j]);
                                    dRtq5[j]=deltaR(tmpteta,tmptphi,etagenq5l[j],phigenq5l[j]);
                                    if(matchingt[i]!=4&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching))) matchingt[i]=6;
                                    if(matchingt[i]==4&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching))) matchingt[i]=46;
                                    t_b.SetPtEtaPhiM(tmpbpt,tmpbeta,tmpbphi,tmpbm);
                                    matchedq1v.SetPtEtaPhiE(ptgenq1l[j],etagenq1l[j],phigenq1l[j],egenq1l[j]);
                                    matchedq2v.SetPtEtaPhiE(ptgenq2l[j],etagenq2l[j],phigenq2l[j],egenq2l[j]);
                                    matchedq3v.SetPtEtaPhiE(ptgenq3l[j],etagenq3l[j],phigenq3l[j],egenq3l[j]);
                                    matchedq4v.SetPtEtaPhiE(ptgenq4l[j],etagenq4l[j],phigenq4l[j],egenq4l[j]);
                                    matchedq5v.SetPtEtaPhiE(ptgenq5l[j],etagenq5l[j],phigenq5l[j],egenq5l[j]);
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>40)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>40)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>40)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>40)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>40))) matchingtq[0][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>60)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>60)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>60)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>60)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>60))) matchingtq[1][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>40&&(t_b+matchedq1v).Mag()<120)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>40&&(t_b+matchedq2v).Mag()<120)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>40&&(t_b+matchedq3v).Mag()<120)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>40&&(t_b+matchedq4v).Mag()<120)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>40&&(t_b+matchedq5v).Mag()<120))) matchingtq[2][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>60&&(t_b+matchedq1v).Mag()<100)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>60&&(t_b+matchedq2v).Mag()<100)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>60&&(t_b+matchedq3v).Mag()<100)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>60&&(t_b+matchedq4v).Mag()<100)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>60&&(t_b+matchedq5v).Mag()<100))) matchingtq[3][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.1*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.1*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.1*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.1*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.1*tmptpt))) matchingtq[4][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*tmptpt))) matchingtq[5][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.2*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.2*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.2*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.2*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.2*tmptpt))) matchingtq[6][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.25*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.25*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.25*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.25*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.25*tmptpt))) matchingtq[7][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&(((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Pt()>0.2*tmpwpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Pt()>0.2*tmpwpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Pt()>0.2*tmpwpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Pt()>0.2*tmpwpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Pt()>0.2*tmpwpt)))) matchingtq[8][i]=6;
                                }
                            }
                            if(dRt_bj[i]<rbmatching&&((dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]>rbmatching)||(dRt_wq1j[i]>rbmatching&&dRt_wq2j[i]<rbmatching))&&tmpt_wtag==4){
                                matchingt[i]=7;
                            }
                            //cout<<tmpt_wtag<<" w tag"<<endl;
                            if(dRt_bj[i]<rbmatching&&dRt_wq1j[i]>rbmatching&&dRt_wq2j[i]>rbmatching&&tmpt_wtag==4){
                                matchingt[i]=-3;
                                if(matchingq[i]!=-99) matchingq[i]+=1;
                                if(matchingq[i]==-99) matchingq[i]=1;}
                            if(dRt_bj[i]<rbmatching&&tmpt_wtag!=4&&dRt_wq1j[i]<rbmatching){
                                //cout<<"   end"<<endl;
                                //cout<<"wlnu  "<<dRt_wq1j[i]<<endl;
                                matchingt[i]=1;
                                for(int j=0;j<10;j++){
                                    dRtg[j]=deltaR(tmpteta,tmptphi,etagengl[j],phigengl[j]);
                                    if(dRtg[j]<rbmatching&&dRgj[j][i]<rbmatching) matchingt[i]=2;
                                    t_b.SetPtEtaPhiM(tmpbpt,tmpbeta,tmpbphi,tmpbm);
                                    matchedgv.SetPtEtaPhiE(ptgengl[j],etagengl[j],phigengl[j],egengl[j]);
                                    if(matchingt[i]==2&&(t_b+matchedgv).Mag()>40) matchingtg[0][i]=2;
                                    if(matchingt[i]==2&&(t_b+matchedgv).Mag()>60) matchingtg[1][i]=2;
                                    if(matchingt[i]==2&&(t_b+matchedgv).Mag()>40&&(t_b+matchedgv).Mag()<120) matchingtg[2][i]=2;
                                    if(matchingt[i]==2&&(t_b+matchedgv).Mag()>60&&(t_b+matchedgv).Mag()<100) matchingtg[3][i]=2;
                                    if(matchingt[i]==2&&ptgengl[j]>0.1*tmptpt) matchingtg[4][i]=2;
                                    if(matchingt[i]==2&&ptgengl[j]>0.15*tmptpt) matchingtg[5][i]=2;
                                    if(matchingt[i]==2&&ptgengl[j]>0.2*tmptpt) matchingtg[6][i]=2;
                                    if(matchingt[i]==2&&ptgengl[j]>0.25*tmptpt) matchingtg[7][i]=2;
                                    if(matchingt[i]==2&&(t_b+matchedgv).Pt()>0.2*tmpwpt) matchingtg[8][i]=2;
                                }
                                for(int j=0;j<5;j++){
                                    dRtq1[j]=deltaR(tmpteta,tmptphi,etagenq1l[j],phigenq1l[j]);
                                    dRtq2[j]=deltaR(tmpteta,tmptphi,etagenq2l[j],phigenq2l[j]);
                                    dRtq3[j]=deltaR(tmpteta,tmptphi,etagenq3l[j],phigenq3l[j]);
                                    dRtq4[j]=deltaR(tmpteta,tmptphi,etagenq4l[j],phigenq4l[j]);
                                    dRtq5[j]=deltaR(tmpteta,tmptphi,etagenq5l[j],phigenq5l[j]);
                                    if(matchingt[i]!=2&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching))) matchingt[i]=5;
                                    if(matchingt[i]==2&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching))) matchingt[i]=25;
                                    t_b.SetPtEtaPhiM(tmpbpt,tmpbeta,tmpbphi,tmpbm);
                                    matchedq1v.SetPtEtaPhiE(ptgenq1l[j],etagenq1l[j],phigenq1l[j],egenq1l[j]);
                                    matchedq2v.SetPtEtaPhiE(ptgenq2l[j],etagenq2l[j],phigenq2l[j],egenq2l[j]);
                                    matchedq3v.SetPtEtaPhiE(ptgenq3l[j],etagenq3l[j],phigenq3l[j],egenq3l[j]);
                                    matchedq4v.SetPtEtaPhiE(ptgenq4l[j],etagenq4l[j],phigenq4l[j],egenq4l[j]);
                                    matchedq5v.SetPtEtaPhiE(ptgenq5l[j],etagenq5l[j],phigenq5l[j],egenq5l[j]);
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>40)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>40)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>40)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>40)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>40))) matchingtq[0][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>60)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>60)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>60)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>60)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>60))) matchingtq[1][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>40&&(t_b+matchedq1v).Mag()<120)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>40&&(t_b+matchedq2v).Mag()<120)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>40&&(t_b+matchedq3v).Mag()<120)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>40&&(t_b+matchedq4v).Mag()<120)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>40&&(t_b+matchedq5v).Mag()<120))) matchingtq[2][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>60&&(t_b+matchedq1v).Mag()<100)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>60&&(t_b+matchedq2v).Mag()<100)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>60&&(t_b+matchedq3v).Mag()<100)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>60&&(t_b+matchedq4v).Mag()<100)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>60&&(t_b+matchedq5v).Mag()<100))) matchingtq[3][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.1*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.1*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.1*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.1*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.1*tmptpt))) matchingtq[4][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*tmptpt))) matchingtq[5][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.2*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.2*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.2*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.2*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.2*tmptpt))) matchingtq[6][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.25*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.25*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.25*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.25*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.25*tmptpt))) matchingtq[7][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&(((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Pt()>0.2*tmpwpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Pt()>0.2*tmpwpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Pt()>0.2*tmpwpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Pt()>0.2*tmpwpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Pt()>0.2*tmpwpt)))) matchingtq[8][i]=5;

                                }
                            }
                            dRt_wq1q2[i]=deltaR(tmpq1eta,tmpq1phi,tmpq2eta,tmpq2phi);
                            if(dRt_bj[i]>rbmatching&&dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]<rbmatching&&dRt_wj[i]<rmatching&&tmpt_wtag==4){
                                matchingt[i]=-2;
                                if(matchingw[i]!=-99) matchingw[i]+=1;
                                if(matchingw[i]==-99) matchingw[i]=1;
                                for(int j=0;j<10;j++){
                                    if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*tmpwpt) matchingw[i]=100;
                                    if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*tmpwpt) matchingw[i]=101;
                                }
                                for(int j=0;j<5;j++){
                                    if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*tmpwpt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*tmpwpt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*tmpwpt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*tmpwpt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*tmpwpt)) matchingw[i]=102;
                                    if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*tmpwpt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*tmpwpt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*tmpwpt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*tmpwpt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*tmpwpt)) matchingw[i]=103;
                                }

                            }
                            if(dRt_bj[i]>rbmatching&&((dRt_wq1j[i]>rbmatching&&dRt_wq2j[i]<rbmatching)||(dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]>rbmatching))&&tmpt_wtag==4){
                                matchingt[i]=-4;
                                if(matchingq[i]!=-99) matchingq[i]+=1;
                                if(matchingq[i]==-99) matchingq[i]=1;
                            }
                        }
                }
                if(!(matchingtmp[i]==-1||matchingtmp[i]==-2)){
                    if(i==0) dRt1_bj[i]=deltaR(gent_b_eta,gent_b_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt1_bj[i]=deltaR(gent_b_eta,gent_b_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt1_bj[i]=deltaR(gent_b_eta,gent_b_phi,eta_min,phi_min);
                    if(i==2) dRt1_bj[i]=deltaR(gent_b_eta,gent_b_phi,eta_min,phi_min);
                    if(i==0) dRt2_bj[i]=deltaR(genantit_b_eta,genantit_b_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt2_bj[i]=deltaR(genantit_b_eta,genantit_b_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt2_bj[i]=deltaR(genantit_b_eta,genantit_b_phi,eta_min,phi_min);
                    if(i==2) dRt2_bj[i]=deltaR(genantit_b_eta,genantit_b_phi,eta_min,phi_min);
                    
                    if(i==0) dRt1_wq1j[i]=deltaR(gent_w_q1_eta,gent_w_q1_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt1_wq1j[i]=deltaR(gent_w_q1_eta,gent_w_q1_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt1_wq1j[i]=deltaR(gent_w_q1_eta,gent_w_q1_phi,eta_min,phi_min);
                    if(i==2) dRt1_wq1j[i]=deltaR(gent_w_q1_eta,gent_w_q1_phi,eta_min,phi_min);
                    if(i==0) dRt2_wq1j[i]=deltaR(genantit_w_q1_eta,genantit_w_q1_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt2_wq1j[i]=deltaR(genantit_w_q1_eta,genantit_w_q1_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt2_wq1j[i]=deltaR(genantit_w_q1_eta,genantit_w_q1_phi,eta_min,phi_min);
                    if(i==2) dRt2_wq1j[i]=deltaR(genantit_w_q1_eta,genantit_w_q1_phi,eta_min,phi_min);
                    if(i==0) dRt1_wq2j[i]=deltaR(gent_w_q2_eta,gent_w_q2_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt1_wq2j[i]=deltaR(gent_w_q2_eta,gent_w_q2_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt1_wq2j[i]=deltaR(gent_w_q2_eta,gent_w_q2_phi,eta_min,phi_min);
                    if(i==2) dRt1_wq2j[i]=deltaR(gent_w_q2_eta,gent_w_q2_phi,eta_min,phi_min);
                    if(i==0) dRt2_wq2j[i]=deltaR(genantit_w_q2_eta,genantit_w_q2_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt2_wq2j[i]=deltaR(genantit_w_q2_eta,genantit_w_q2_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt2_wq2j[i]=deltaR(genantit_w_q2_eta,genantit_w_q2_phi,eta_min,phi_min);
                    if(i==2) dRt2_wq2j[i]=deltaR(genantit_w_q2_eta,genantit_w_q2_phi,eta_min,phi_min);
                    
                    if(i==0) dRt1_wj[i]=deltaR(gent_w_eta,gent_w_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt1_wj[i]=deltaR(gent_w_eta,gent_w_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt1_wj[i]=deltaR(gent_w_eta,gent_w_phi,eta_min,phi_min);
                    if(i==2) dRt1_wj[i]=deltaR(gent_w_eta,gent_w_phi,eta_min,phi_min);
                    if(i==0) dRt2_wj[i]=deltaR(genantit_w_eta,genantit_w_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt2_wj[i]=deltaR(genantit_w_eta,genantit_w_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt2_wj[i]=deltaR(genantit_w_eta,genantit_w_phi,eta_min,phi_min);
                    if(i==2) dRt2_wj[i]=deltaR(genantit_w_eta,genantit_w_phi,eta_min,phi_min);
                    dRt1_wb[i]=deltaR(gent_w_eta,gent_w_phi,gent_b_eta,gent_b_phi);
                    dRt2_wb[i]=deltaR(genantit_w_eta,genantit_w_phi,genantit_b_eta,genantit_b_phi);
                    
                    if((dRt1_wq1j[i]<rbmatching&&dRt1_wq2j[i]<rbmatching&&dRt1_bj[i]>rbmatching&&dRt1_wj[i]<rmatching&&gent_w_tag==4)) {
                        matchingt[i]=-5;
                        if(matchingw[i]!=-99) matchingw[i]+=1;
                        if(matchingw[i]==-99) matchingw[i]=1;
                        for(int j=0;j<10;j++){
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*gent_w_pt) matchingw[i]=100;
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*gent_w_pt) matchingw[i]=101;
                        }
                        for(int j=0;j<5;j++){
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*gent_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*gent_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*gent_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*gent_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*gent_w_pt)) matchingw[i]=102;
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*gent_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*gent_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*gent_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*gent_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*gent_w_pt)) matchingw[i]=103;
                        }

                    }
                    if((dRt2_wq1j[i]<rbmatching&&dRt2_wq2j[i]<rbmatching&&dRt2_bj[i]>rbmatching&&dRt2_wj[i]<rmatching&&genantit_w_tag==4)) {
                        matchingt[i]=-6;
                        if(matchingw[i]!=-99) matchingw[i]+=1;
                        if(matchingw[i]==-99) matchingw[i]=1;
                        for(int j=0;j<10;j++){
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*genantit_w_pt) matchingw[i]=100;
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*genantit_w_pt) matchingw[i]=101;
                        }
                        for(int j=0;j<5;j++){
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*genantit_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*genantit_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*genantit_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*genantit_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*genantit_w_pt)) matchingw[i]=102;
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*genantit_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*genantit_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*genantit_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*genantit_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*genantit_w_pt)) matchingw[i]=103;
                        }

                    }

                    if((dRt1_wq1j[i]<rbmatching&&dRt1_wq2j[i]>rbmatching&&dRt1_bj[i]>rbmatching&&gent_w_tag==4)||(dRt1_wq1j[i]>rbmatching&&dRt1_wq2j[i]<rbmatching&&dRt1_bj[i]>rbmatching&&gent_w_tag==4)||(dRt1_wq1j[i]>rbmatching&&dRt1_wq2j[i]>rbmatching&&dRt1_bj[i]<rbmatching&&gent_w_tag==4)) {
                        matchingt[i]=-7;
                        if(matchingq[i]!=-99) matchingq[i]+=1;
                        if(matchingq[i]==-99) matchingq[i]=1;
                    }
                    if((dRt2_wq1j[i]<rbmatching&&dRt2_wq2j[i]>rbmatching&&dRt2_bj[i]>rbmatching&&genantit_w_tag==4)||(dRt2_wq1j[i]>rbmatching&&dRt2_wq2j[i]<rbmatching&&dRt2_bj[i]>rbmatching&&genantit_w_tag==4)||(dRt2_wq1j[i]>rbmatching&&dRt2_wq2j[i]>rbmatching&&dRt2_bj[i]<rbmatching&&genantit_w_tag==4)) {
                        matchingt[i]=-8;
                        if(matchingq[i]!=-99) matchingq[i]+=1;
                        if(matchingq[i]==-99) matchingq[i]=1;
                    }

                }
                for(int j=0;j<5;j++){
                    if (dRwj[j][i]<rmatching&&ptgenwf[j]!=gent_w_pt&&ptgenwf[j]!=genantit_w_pt){
                        if(i==0) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_min,phi_min);
                        if(i==2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_min,phi_min);
                        if(i==0) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_min,phi_min);
                        if(i==2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_min,phi_min);
                        if(dRw_q1j[j][i]<rbmatching&&dRw_q2j[j][i]<rbmatching&&taggenwl[j]==4) {
                            if(matchingw[i]!=-99) matchingw[i]+=1;
                            if(matchingw[i]==-99) matchingw[i]=1;
                            for(int j=0;j<10;j++){
                                if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*gent_w_pt) matchingw[i]=100;
                                if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*gent_w_pt) matchingw[i]=101;
                            }
                            for(int j=0;j<5;j++){
                                if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*gent_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*gent_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*gent_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*gent_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*gent_w_pt)) matchingw[i]=102;
                                if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*gent_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*gent_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*gent_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*gent_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*gent_w_pt)) matchingw[i]=103;
                            }

                        }
                        if(((dRw_q1j[j][i]<rbmatching&&dRw_q2j[j][i]>rbmatching)||(dRw_q1j[j][i]>rbmatching&&dRw_q2j[j][i]<rbmatching))&&taggenwl[j]==4) {
                            if(matchingq[i]!=-99) matchingq[i]+=1;
                            if(matchingq[i]==-99) matchingq[i]=1;
                        }
                    }
                }

                    for(int j=0;j<5;j++){
                        if (dRzj[j][i]<rmatching){
                            if(matchingz[i]!=-99) matchingz[i]+=1;
                            if(matchingz[i]==-99) matchingz[i]=1;
                        }
                    }
                    for(int j=0;j<10;j++){
                        if (dRgj[j][i]<rmatching){
                            if(matchingg[i]!=-99) matchingg[i]+=1;
                            if(matchingg[i]==-99) matchingg[i]=1;
                        }
                    }
                    for(int j=0;j<5;j++){
                        if ((dRq1j[j][i]<rmatching)||(dRq2j[j][i]<rmatching)||(dRq3j[j][i]<rmatching)||(dRq4j[j][i]<rmatching)||(dRq5j[j][i]<rmatching)){
                            if(matchingq[i]!=-99) matchingq[i]+=1;
                            if(matchingq[i]==-99) matchingq[i]=1;
                        }
                    }
                //if((matchingtmp[i]==-1||matchingtmp[i]==-2)) cout<<i<<"  matched t "<<dRtj[i]<<"   "<<dRantitj[i]<<"  "<<dRt_bj[i]<<"   "<<dRt_wq1j[i]<<"   "<<dRt_wq2j[i]<<"   "<<matchingt[i]<<"   "<<matchingw[i]<<"   "<<matchingq[i]<<endl;

                //if(!(matchingtmp[i]==-1||matchingtmp[i]==-2)) cout<<i<<" matched no t  "<<dRtj[i]<<"   "<<dRantitj[i]<<"  "<<dRt1_bj[i]<<"   "<<dRt1_wq1j[i]<<"   "<<dRt1_wq2j[i]<<"  "<<dRt2_bj[i]<<"   "<<dRt2_wq1j[i]<<"   "<<dRt2_wq2j[i]<<"   "<<matchingt[i]<<"   "<<matchingw[i]<<"   "<<matchingq[i]<<"   "<<matchingg[i]<<endl;
                
                if(matchingt[i]<0&&matchingw[i]<0&&matchingg[i]<0&&matchingz[i]<0&&matchingq[i]<0&&matchingr[i]<0) unmatching[i]=1;
            }
            if(gentop_pt<0&&genantitop_pt<0){
                for(int j=0;j<5;j++){
                    if (dRwj[j][i]<rmatching){
                        if(i==0) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_min,phi_min);
                        if(i==2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_min,phi_min);
                        if(i==0) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_min,phi_min);
                        if(i==2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_min,phi_min);
                        if(dRw_q1j[j][i]<rbmatching&&dRw_q2j[j][i]<rbmatching&&taggenwl[j]==4) {
                            if(matchingw[i]!=-99) matchingw[i]+=1;
                            if(matchingw[i]==-99) matchingw[i]=1;
                            for(int j=0;j<10;j++){
                                if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*genantit_w_pt) matchingw[i]=100;
                                if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*genantit_w_pt) matchingw[i]=101;
                            }
                            for(int j=0;j<5;j++){
                                if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*genantit_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*genantit_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*genantit_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*genantit_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*genantit_w_pt)) matchingw[i]=102;
                                if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*genantit_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*genantit_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*genantit_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*genantit_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*genantit_w_pt)) matchingw[i]=103;
                            }

                        }
                        if(((dRw_q1j[j][i]<rbmatching&&dRw_q2j[j][i]>rbmatching)||(dRw_q1j[j][i]>rbmatching&&dRw_q2j[j][i]<rbmatching))&&taggenwl[j]==4) {
                            if(matchingq[i]!=-99) matchingq[i]+=1;
                            if(matchingq[i]==-99) matchingq[i]=1;
                        }
                    }
                }
                if(IsData==10){
                    double dRw2w3=99;
                    /*if(status_2==1) dRqq_W1=deltaR(gen_ele_eta_2,gen_ele_phi_2,gen_nele_eta_2,gen_nele_phi_2);
                    if(status_2==2) dRqq_W1=deltaR(gen_mu_eta_2,gen_mu_phi_2,gen_nmu_eta_2,gen_nmu_phi_2);
                    if(status_2==3) dRqq_W1=deltaR(gen_tau_eta_2,gen_tau_phi_2,gen_ntau_eta_2,gen_ntau_phi_2);
                    if(status_2==4) dRqq_W1=deltaR(etaq_2[0],phiq_2[0],etaq_2[1],phiq_2[1]);
                    if(status_3==1) dRqq_W2=deltaR(gen_ele_eta_3,gen_ele_phi_3,gen_nele_eta_3,gen_nele_phi_3);
                    if(status_3==2) dRqq_W2=deltaR(gen_mu_eta_3,gen_mu_phi_3,gen_nmu_eta_3,gen_nmu_phi_3);
                    if(status_3==3) dRqq_W2=deltaR(gen_tau_eta_3,gen_tau_phi_3,gen_ntau_eta_3,gen_ntau_phi_3);
                    if(status_3==4) dRqq_W2=deltaR(etaq_3[0],phiq_3[0],etaq_3[1],phiq_3[1]);
                    mindRqq=TMath::Min(dRqq_W1,dRqq_W2);
                    maxdRqq=TMath::Max(dRqq_W1,dRqq_W2);*/

                    if(status_2==4&&status_3==4) {
                        dRw2w3=deltaR(etaGenVhad_2,phiGenVhad_2,etaGenVhad_3,phiGenVhad_3);
                        if(i==0) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_min,phi_min);
                        if(i==2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_min,phi_min);
                        if(i==0) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_min,phi_min);
                        if(i==2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_min,phi_min);
                        dRWW_R=deltaR(etaGenVhad_2,phiGenVhad_2,etaGenVhad_3,phiGenVhad_3);
                    }
                    if(status_2==4&&status_3!=4) {
                        dRw2w3=deltaR(etaGenVhad_2,phiGenVhad_2,etaGenV_3,phiGenV_3);
                        if(i==0) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_min,phi_min);
                        if(i==2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_min,phi_min);
                        if(i==0) dRw3j[i]=deltaR(etaGenV_3,phiGenV_3,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw3j[i]=deltaR(etaGenV_3,phiGenV_3,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw3j[i]=deltaR(etaGenV_3,phiGenV_3,eta_min,phi_min);
                        if(i==2) dRw3j[i]=deltaR(etaGenV_3,phiGenV_3,eta_min,phi_min);
                        dRWW_R=deltaR(etaGenVhad_2,phiGenVhad_2,etaGenV_3,phiGenV_3);
                    }
                    if(status_2!=4&&status_3==4) {
                        dRw2w3=deltaR(etaGenV_2,phiGenV_2,etaGenVhad_3,phiGenVhad_3);
                        if(i==0) dRw2j[i]=deltaR(etaGenV_2,phiGenV_2,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw2j[i]=deltaR(etaGenV_2,phiGenV_2,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw2j[i]=deltaR(etaGenV_2,phiGenV_2,eta_min,phi_min);
                        if(i==2) dRw2j[i]=deltaR(etaGenV_2,phiGenV_2,eta_min,phi_min);
                        if(i==0) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_min,phi_min);
                        if(i==2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_min,phi_min);
                        dRWW_R=deltaR(etaGenV_2,phiGenV_2,etaGenVhad_3,phiGenVhad_3);
                    }

                    
                    if(dRrj[i]<rmatching&&dRw2j[i]<rmatching&&dRw3j[i]<rmatching){
                    //if(dRrj[i]<rmatching){
                        if(i==0) dRr_q1j[i]=deltaR(etaq_2[0],phiq_2[0],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRr_q1j[i]=deltaR(etaq_2[0],phiq_2[0],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRr_q1j[i]=deltaR(etaq_2[0],phiq_2[0],eta_min,phi_min);
                        if(i==2) dRr_q1j[i]=deltaR(etaq_2[0],phiq_2[0],eta_min,phi_min);
                        if(i==0) dRr_q2j[i]=deltaR(etaq_2[1],phiq_2[1],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRr_q2j[i]=deltaR(etaq_2[1],phiq_2[1],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRr_q2j[i]=deltaR(etaq_2[1],phiq_2[1],eta_min,phi_min);
                        if(i==2) dRr_q2j[i]=deltaR(etaq_2[1],phiq_2[1],eta_min,phi_min);
                        if(i==0) dRr_q3j[i]=deltaR(etaq_3[0],phiq_3[0],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRr_q3j[i]=deltaR(etaq_3[0],phiq_3[0],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRr_q3j[i]=deltaR(etaq_3[0],phiq_3[0],eta_min,phi_min);
                        if(i==2) dRr_q3j[i]=deltaR(etaq_3[0],phiq_3[0],eta_min,phi_min);
                        if(i==0) dRr_q4j[i]=deltaR(etaq_3[1],phiq_3[1],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRr_q4j[i]=deltaR(etaq_3[1],phiq_3[1],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRr_q4j[i]=deltaR(etaq_3[1],phiq_3[1],eta_min,phi_min);
                        if(i==2) dRr_q4j[i]=deltaR(etaq_3[1],phiq_3[1],eta_min,phi_min);
                        double dRr_qj[4]={99,99,99,99};
                        dRr_qj[0]=dRr_q1j[i];dRr_qj[1]=dRr_q2j[i];dRr_qj[2]=dRr_q3j[i];dRr_qj[3]=dRr_q4j[i];
                        Int_t *idRr_qj=new Int_t[4];
                        TMath::Sort(4,dRr_qj,idRr_qj,0);
                        //cout<<dRr_qj[idRr_qj[0]]<<"  "<<dRr_qj[idRr_qj[1]]<<"  "<<dRr_qj[idRr_qj[2]]<<"  "<<dRr_qj[idRr_qj[3]]<<endl;
                        if(dRr_qj[idRr_qj[3]]<rbmatching) {
                        matchingr[i]=4;
                        for(int j=0;j<10;j++){
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*gen_rad_pt) matchingr[i]=5;
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*gen_rad_pt) matchingr[i]=6;
                        }
                        for(int j=0;j<5;j++){
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*gen_rad_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*gen_rad_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*gen_rad_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*gen_rad_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*gen_rad_pt)) matchingr[i]=7;
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*gen_rad_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*gen_rad_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*gen_rad_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*gen_rad_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*gen_rad_pt)) matchingr[i]=8;
                        }
                        }
                        if(dRr_qj[idRr_qj[2]]<rbmatching&&dRr_qj[idRr_qj[3]]>rbmatching) matchingr[i]=3;
                        if(dRr_qj[idRr_qj[1]]<rbmatching&&dRr_qj[idRr_qj[2]]>rbmatching) matchingr[i]=2;
                        if(dRr_qj[idRr_qj[0]]<rbmatching&&dRr_qj[idRr_qj[1]]>rbmatching){
                            if(status_2==4&&status_3==4) matchingr[i]=-1;
                            if(!(status_2==4&&status_3==4)) matchingr[i]=1;
                        }
                    }
                }
                for(int j=0;j<5;j++){
                    if (dRzj[j][i]<rmatching){
                            if(matchingz[i]!=-99) matchingz[i]+=1;
                            if(matchingz[i]==-99) matchingz[i]=1;
                    }
                }
                for(int j=0;j<10;j++){
                    if (dRgj[j][i]<rmatching){
                        if(matchingg[i]!=-99) matchingg[i]+=1;
                        if(matchingg[i]==-99) matchingg[i]=1;
                    }
                }
                for(int j=0;j<5;j++){
                    if ((dRq1j[j][i]<rmatching)||(dRq2j[j][i]<rmatching)||(dRq3j[j][i]<rmatching)||(dRq4j[j][i]<rmatching)||(dRq5j[j][i]<rmatching)){
                        if(matchingq[i]!=-99) matchingq[i]+=1;
                        if(matchingq[i]==-99) matchingq[i]=1;
                    }
                }
            }
            //cout<<matchingr[i]<<" cc  "<<matchingw[i]<<"   "<<matchingq[i]<<"   "<<matchingg[i]<<endl;
            if(matchingt[i]<0&&matchingw[i]<0&&matchingg[i]<0&&matchingz[i]<0&&matchingq[i]<0&&matchingr[i]<0) unmatching[i]=1;
        }
        
        
// unknown part

 
          double etaq1=-99,etaq2=-99,etaq3=-99,etaq4=-99,phiq1=-99,phiq2=-99,phiq3=-99,phiq4=-99,ptq1=-99,ptq2=-99,ptq3=-99,ptq4=-99;
        double q4v_R[4][3]={-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99};
        for(int i=0;i<4;i++){
            nqqcr[i]=0;nqqct[i]=0;nqqcw[i]=0;
            for(int j=0;j<3;j++)    {matchnewr[i][j]=-99;matchnewt[i][j]=-99;matchneww[i][j]=-99;}
        }
        double ignoreprong=0.15;
        if(IsData==10){
            if(status_2==1){etaq1=gen_ele_eta_2;phiq1=gen_ele_phi_2;ptq1=gen_ele_pt_2;etaq2=gen_nele_eta_2;phiq2=gen_nele_phi_2;ptq2=gen_nele_phi_2;}
            if(status_2==2){etaq1=gen_mu_eta_2;phiq1=gen_mu_phi_2;etaq2=gen_nmu_eta_2;phiq2=gen_nmu_phi_2;}
            if(status_2==3){etaq1=gen_tau_eta_2;phiq1=gen_tau_phi_2;etaq2=gen_ntau_eta_2;phiq2=gen_ntau_phi_2;}
            if(status_2==4){etaq1=etaq_2[0];phiq1=phiq_2[0];etaq2=etaq_2[1];phiq2=phiq_2[1];}
            if(status_3==1){etaq3=gen_ele_eta_3;phiq3=gen_ele_phi_3;etaq4=gen_nele_eta_3;phiq4=gen_nele_phi_3;}
            if(status_3==2){etaq3=gen_mu_eta_3;phiq3=gen_mu_phi_3;etaq4=gen_nmu_eta_3;phiq4=gen_nmu_phi_3;}
            if(status_3==3){etaq3=gen_tau_eta_3;phiq3=gen_tau_phi_3;etaq4=gen_ntau_eta_3;phiq4=gen_ntau_phi_3;}
            if(status_3==4){etaq3=etaq_3[0];phiq3=phiq_3[0];etaq4=etaq_3[1];phiq4=phiq_3[1];}
            
            if(status_2==1){q4v_R[0][0]=gen_ele_pt_2;q4v_R[0][1]=gen_ele_eta_2;q4v_R[0][2]=gen_ele_phi_2;q4v_R[1][0]=gen_nele_pt_2;q4v_R[1][1]=gen_nele_eta_2;q4v_R[1][2]=gen_nele_phi_2;}
            if(status_2==2){q4v_R[0][0]=gen_mu_pt_2;q4v_R[0][1]=gen_mu_eta_2;q4v_R[0][2]=gen_mu_phi_2;q4v_R[1][0]=gen_nmu_pt_2;q4v_R[1][1]=gen_nmu_eta_2;q4v_R[1][2]=gen_nmu_phi_2;}
            if(status_2==3){q4v_R[0][0]=gen_tau_pt_2;q4v_R[0][1]=gen_tau_eta_2;q4v_R[0][2]=gen_tau_phi_2;q4v_R[1][0]=gen_ntau_pt_2;q4v_R[1][1]=gen_ntau_eta_2;q4v_R[1][2]=gen_ntau_phi_2;}
            if(status_2==4){q4v_R[0][0]=ptq_2[0];q4v_R[0][1]=etaq_2[0];q4v_R[0][2]=phiq_2[0];q4v_R[1][0]=ptq_2[1];q4v_R[1][1]=etaq_2[1];q4v_R[1][2]=phiq_2[1];}
            if(status_3==1){q4v_R[2][0]=gen_ele_pt_3;q4v_R[2][1]=gen_ele_eta_3;q4v_R[2][2]=gen_ele_phi_3;q4v_R[3][0]=gen_nele_pt_3;q4v_R[3][1]=gen_nele_eta_3;q4v_R[3][2]=gen_nele_phi_3;}
            if(status_3==2){q4v_R[2][0]=gen_mu_pt_3;q4v_R[2][1]=gen_mu_eta_3;q4v_R[2][2]=gen_mu_phi_3;q4v_R[3][0]=gen_nmu_pt_3;q4v_R[3][1]=gen_nmu_eta_3;q4v_R[3][2]=gen_nmu_phi_3;}
            if(status_3==3){q4v_R[2][0]=gen_tau_pt_3;q4v_R[2][1]=gen_tau_eta_3;q4v_R[2][2]=gen_tau_phi_3;q4v_R[3][0]=gen_ntau_pt_3;q4v_R[3][1]=gen_ntau_eta_3;q4v_R[3][2]=gen_ntau_phi_3;}
            if(status_3==4){q4v_R[2][0]=ptq_3[0];q4v_R[2][1]=etaq_3[0];q4v_R[2][2]=phiq_3[0];q4v_R[3][0]=ptq_3[1];q4v_R[3][1]=etaq_3[1];q4v_R[3][2]=phiq_3[1];}

            dRqq_W1=deltaR(q4v_R[0][1],q4v_R[0][2],q4v_R[1][1],q4v_R[1][2]);
            dRqq_W2=deltaR(q4v_R[2][1],q4v_R[2][2],q4v_R[3][1],q4v_R[3][2]);
            mindRqq=TMath::Min(dRqq_W1,dRqq_W2);
            maxdRqq=TMath::Max(dRqq_W1,dRqq_W2);
            
             double dRqq_all[6]={deltaR(q4v_R[0][1],q4v_R[0][2],q4v_R[1][1],q4v_R[1][2]),deltaR(q4v_R[0][1],q4v_R[0][2],q4v_R[2][1],q4v_R[2][2]),deltaR(q4v_R[0][1],q4v_R[0][2],q4v_R[3][1],q4v_R[3][2]),deltaR(q4v_R[1][1],q4v_R[1][2],q4v_R[2][1],q4v_R[2][2]),deltaR(q4v_R[1][1],q4v_R[1][2],q4v_R[3][1],q4v_R[3][2]),deltaR(q4v_R[2][1],q4v_R[2][2],q4v_R[3][1],q4v_R[3][2])};
             Int_t *index_all=new Int_t[6];
             TMath::Sort(6,dRqq_all,index_all,0);
             mindRqq_all=dRqq_all[index_all[0]];
             maxdRqq_all=dRqq_all[index_all[5]];

            double extraqgr[5][3]={0};
            int insideqgr=0;
            for(int j=0;j<10;j++){
                if(dRgj[j][0]<rbmatching&&ptgengl[j]>0.15*gen_rad_pt&&ptgengl[j]>extraqgr[0][0]) {extraqgr[0][0]=ptgengl[j];extraqgr[0][1]=etagengl[j];extraqgr[0][2]=phigengl[j];}
            }
            for(int j=0;j<5;j++){
                if((dRq1j[j][0]<rbmatching&&ptgenq1l[j]>0.15*gen_rad_pt&&ptgenq1l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq1l[j];extraqgr[0][1]=etagenq1l[j];extraqgr[0][2]=phigenq1l[j];}
                if((dRq2j[j][0]<rbmatching&&ptgenq2l[j]>0.15*gen_rad_pt&&ptgenq2l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq2l[j];extraqgr[0][1]=etagenq2l[j];extraqgr[0][2]=phigenq2l[j];}
                if((dRq3j[j][0]<rbmatching&&ptgenq3l[j]>0.15*gen_rad_pt&&ptgenq3l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq3l[j];extraqgr[0][1]=etagenq3l[j];extraqgr[0][2]=phigenq3l[j];}
                if((dRq4j[j][0]<rbmatching&&ptgenq4l[j]>0.15*gen_rad_pt&&ptgenq4l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq4l[j];extraqgr[0][1]=etagenq4l[j];extraqgr[0][2]=phigenq4l[j];}
                if((dRq5j[j][0]<rbmatching&&ptgenq5l[j]>0.15*gen_rad_pt&&ptgenq5l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq5l[j];extraqgr[0][1]=etagenq5l[j];extraqgr[0][2]=phigenq5l[j];}
            }
            dRr_q1j[0]=deltaR(q4v_R[0][1],q4v_R[0][2],eta_max,phi_max);
            dRr_q2j[0]=deltaR(q4v_R[1][1],q4v_R[1][2],eta_max,phi_max);
            dRr_q3j[0]=deltaR(q4v_R[2][1],q4v_R[2][2],eta_max,phi_max);
            dRr_q4j[0]=deltaR(q4v_R[3][1],q4v_R[3][2],eta_max,phi_max);
            if(status_2==4&&status_3==4) {
                dRw2j[0]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_max,phi_max);
                dRw3j[0]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_max,phi_max);
            }
            if(status_2==4&&status_3!=4) {
                dRw2j[0]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_max,phi_max);
                dRw3j[0]=deltaR(etaGenV_3,phiGenV_3,eta_max,phi_max);
            }
            if(status_2!=4&&status_3==4) {
                dRw2j[0]=deltaR(etaGenV_2,phiGenV_2,eta_max,phi_max);
                dRw3j[0]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_max,phi_max);
            }
            int tmpinr=0;
            //if(dRrj[0]<rmatching){
            if(dRrj[0]<rmatching&&dRw2j[0]<rmatching&&dRw3j[0]<rmatching){
                if(extraqgr[0][0]>0) tmpinr++;
                //if(tmpinr==1){cout<<" checkgpt "<<extraqgr[0][0]<<endl;}
                if(dRr_q1j[0]<rbmatching){extraqgr[tmpinr][0]=q4v_R[0][0];extraqgr[tmpinr][1]=q4v_R[0][1];extraqgr[tmpinr][2]=q4v_R[0][2];tmpinr++;}
                if(dRr_q2j[0]<rbmatching&&status_2==4){extraqgr[tmpinr][0]=q4v_R[1][0];extraqgr[tmpinr][1]=q4v_R[1][1];extraqgr[tmpinr][2]=q4v_R[1][2];tmpinr++;}
                if(dRr_q3j[0]<rbmatching){extraqgr[tmpinr][0]=q4v_R[2][0];extraqgr[tmpinr][1]=q4v_R[2][1];extraqgr[tmpinr][2]=q4v_R[2][2];tmpinr++;}
                if(dRr_q4j[0]<rbmatching&&status_3==4){extraqgr[tmpinr][0]=q4v_R[3][0];extraqgr[tmpinr][1]=q4v_R[3][1];extraqgr[tmpinr][2]=q4v_R[3][2];tmpinr++;}
                /*for(int i=0;i<tmpinr;i++){
                    cout<<"check inr pt"<<i<<"    "<<extraqgr[i][0];}
                cout<<endl;*/
            }
            for(int k=0;k<4;k++) {matchnewr[k][0]=tmpinr;matchnewr[k][2]=tmpinr;
                if(extraqgr[0][0]>0)matchnewr[k][2]=tmpinr-1;
            }
            
            double PTtar=gen_rad_pt;
            const int maxgroupr=tmpinr;
            int pgroupr[4][2*maxgroupr],subgroupr[4][maxgroupr][maxgroupr],nsubgroupr[4][maxgroupr],flagw1[4][2*maxgroupr],flagalone[4][maxgroupr];
            double ptgroupr[4][maxgroupr];
            for(int i=0;i<4;i++){
                for(int j=0;j<maxgroupr;j++){
                    ptgroupr[i][j]=0;nsubgroupr[i][j]=-99;flagalone[i][j]=-99;
                    for(int k=0;k<maxgroupr;k++) subgroupr[i][j][k]=-99;}
                for(int j=0;j<2*maxgroupr;j++) {pgroupr[i][j]=-99;flagw1[i][j]=-99;}
            }
            double drqqin5[2*tmpinr];
            for(int i=0;i<2*tmpinr;i++){drqqin5[i]=-99;}
            int idrqqin5=0;
            for(int m=0;m<tmpinr;m++){
                for(int n=m+1;n<tmpinr;n++){
                    drqqin5[idrqqin5]=deltaR(extraqgr[m][1],extraqgr[m][2],extraqgr[n][1],extraqgr[n][2]);
                    //cout<<m<<"   "<<n<<"   "<<extraqgr[m][1]<<"   "<<extraqgr[m][2]<<"   "<<extraqgr[n][1]<<"   "<<extraqgr[n][2]<<"   "<<drqqin5[idrqqin5]<<endl;
                    idrqqin5++;
                }
            }
            int i1sub,i2sub;
            double drwp[4]={0.1,0.15,0.2,0.25};
            int end1,end2;
            for(int i=0;i<4;i++){
                i1sub=0;i2sub=0;
            for(int j=0;j<2*tmpinr;j++){
                if(0<drqqin5[j]&&drqqin5[j]<drwp[i]){
                    end1=endpoint(j,tmpinr,1);end2=endpoint(j,tmpinr,2);
                    //cout<<" check  "<<j<<"   "<<end1<<"   "<<end2<<"   "<<flagw1[i][end1]<<"   "<<flagw1[i][end2]<<" i1sub  "<<i1sub<<"   "<<i2sub<<"   "<<endl;
                    if(flagw1[i][end1]==-99&&flagw1[i][end2]==-99) {if(subgroupr[i][i1sub][0]!=-99){i1sub++;i2sub=0;}subgroupr[i][i1sub][i2sub]=end1;subgroupr[i][i1sub][i2sub+1]=end2;
                        flagw1[i][end1]=1;flagw1[i][end2]=1;i2sub+=2;
                        for(int k=0;k<2*tmpinr;k++){
                            if(0<drqqin5[k]&&drqqin5[k]<drwp[i]){
                                end1=endpoint(k,tmpinr,1);end2=endpoint(k,tmpinr,2);
                                if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                    int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                    subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                    for(int m=0;m<2*tmpinr;m++){
                                        if(0<drqqin5[m]&&drqqin5[m]<drwp[i]){
                                            end1=endpoint(m,tmpinr,1);end2=endpoint(m,tmpinr,2);
                                            if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                for(int n=0;n<2*tmpinr;n++){
                                                    if(0<drqqin5[n]&&drqqin5[n]<drwp[i]){
                                                        end1=endpoint(n,tmpinr,1);end2=endpoint(n,tmpinr,2);
                                                        if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                            int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                            subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                        }}}
                                            }}}
                                }
                            }
                        }
                    }
                }//end
            }
                int totalgroup=0;
                for(int j=0;j<tmpinr;j++){if(subgroupr[i][j][0]!=-99)totalgroup++;}
                const int ntotal=totalgroup;
                int flagw2[ntotal][tmpinr];
                for(int j=0;j<tmpinr;j++){
                    for(int k=0;k<ntotal;k++){
                        flagw2[k][j]=-99;
                    for(int m=0;m<tmpinr;m++){
                        if(subgroupr[i][k][m]==j){flagw2[k][j]=1;}
                    }}
                    int flagw3=-99;
                    for(int k=0;k<ntotal;k++){
                        if(flagw2[k][j]==1) flagw3=1;}
                    if(flagw3==-99){subgroupr[i][totalgroup][0]=j;totalgroup++;}
                     }
            }
            for(int i=0;i<4;i++){
                for(int j=0;j<tmpinr;j++){
                    int igroup=0;
                    double tmpsubg[5][2]={0};
                    for(int k=0;k<tmpinr;k++){
                        if(subgroupr[i][j][k]!=-99) {tmpsubg[igroup][0]=extraqgr[subgroupr[i][j][k]][0];tmpsubg[igroup][1]=extraqgr[subgroupr[i][j][k]][2];igroup++;}
                    }
                    ptgroupr[i][j]=SumPt(tmpsubg,igroup);
                    if(ptgroupr[i][j]/PTtar>ignoreprong) {if(matchnewr[i][1]!=-99) matchnewr[i][1]++;if(matchnewr[i][1]==-99) matchnewr[i][1]=1;}
                    /*cout<<"groupr"<<j;
                    for(int k=0;k<tmpinr;k++)cout<<"   "<<subgroupr[i][j][k];
                    cout<<"   "<<ptgroupr[i][j]<<"    "<<ptgroupr[i][j]/PTtar<<"    "<<matchnewr[i][1]<<endl;*/
            }
            }
        }
        
        if(IsData>0&&(gentop_pt>0||genantitop_pt>0)){
            
            double extraqgt[4][3]={0};
            if (dRtj[0]<rmatching) matchingtmp[0]=-1;
            if (matchingtmp[0]!=-1&&dRantitj[0]<rmatching) matchingtmp[0]=-2;
            if (matchingtmp[0]==-1&&dRantitj[0]<rmatching&&dRtj[0]<dRantitj[0]) matchingtmp[0]=-1;
            if (matchingtmp[0]==-1&&dRantitj[0]<rmatching&&dRtj[0]>dRantitj[0]) matchingtmp[0]=-2;
            //cout<<i<<" jet  "<<Nj8<<"   "<<matchingtmp[0]<<endl;
            if(matchingtmp[0]==-1||matchingtmp[0]==-2){
                if(matchingtmp[0]==-1){tmptpt=gentop_pt;tmpteta=gentop_eta;tmptphi=gentop_phi;tmpbpt=gent_b_pt;tmpbeta=gent_b_eta;tmpbphi=gent_b_phi;tmpwpt=gent_w_pt;tmpweta=gent_w_eta;tmpwphi=gent_w_phi;tmpq1pt=gent_w_q1_pt;tmpq1eta=gent_w_q1_eta;tmpq1phi=gent_w_q1_phi;tmpq2pt=gent_w_q2_pt;tmpq2eta=gent_w_q2_eta;tmpq2phi=gent_w_q2_phi;tmpbm=gent_b_mass;tmpt_wtag=gent_w_tag;}
                if(matchingtmp[0]==-2){tmptpt=genantitop_pt;tmpteta=genantitop_eta;tmptphi=genantitop_phi;tmpbpt=genantit_b_pt;tmpbeta=genantit_b_eta;tmpbphi=genantit_b_phi;tmpwpt=genantit_w_pt;tmpweta=genantit_w_eta;tmpwphi=genantit_w_phi;tmpq1pt=genantit_w_q1_pt;tmpq1eta=genantit_w_q1_eta;tmpq1phi=genantit_w_q1_phi;tmpq2pt=genantit_w_q2_pt;tmpq2eta=genantit_w_q2_eta;tmpq2phi=genantit_w_q2_phi;tmpbm=genantit_b_mass;tmpt_wtag=genantit_w_tag;}

            int insideqgr=0;
            for(int j=0;j<10;j++){
                if(dRgj[j][0]<rbmatching&&ptgengl[j]>0.15*tmptpt&&ptgengl[j]>extraqgt[0][0]) {extraqgt[0][0]=ptgengl[j];extraqgt[0][1]=etagengl[j];extraqgt[0][2]=phigengl[j];}
            }
            for(int j=0;j<5;j++){
                if((dRq1j[j][0]<rbmatching&&ptgenq1l[j]>0.15*tmptpt&&ptgenq1l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq1l[j];extraqgt[0][1]=etagenq1l[j];extraqgt[0][2]=phigenq1l[j];}
                if((dRq2j[j][0]<rbmatching&&ptgenq2l[j]>0.15*tmptpt&&ptgenq2l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq2l[j];extraqgt[0][1]=etagenq2l[j];extraqgt[0][2]=phigenq2l[j];}
                if((dRq3j[j][0]<rbmatching&&ptgenq3l[j]>0.15*tmptpt&&ptgenq3l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq3l[j];extraqgt[0][1]=etagenq3l[j];extraqgt[0][2]=phigenq3l[j];}
                if((dRq4j[j][0]<rbmatching&&ptgenq4l[j]>0.15*tmptpt&&ptgenq4l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq4l[j];extraqgt[0][1]=etagenq4l[j];extraqgt[0][2]=phigenq4l[j];}
                if((dRq5j[j][0]<rbmatching&&ptgenq5l[j]>0.15*tmptpt&&ptgenq5l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq5l[j];extraqgt[0][1]=etagenq5l[j];extraqgt[0][2]=phigenq5l[j];}
            }
            dRt_bj[0]=deltaR(tmpbeta,tmpbphi,eta_max,phi_max);
                if(dRt_bj[0]<rbmatching){
            dRt_wq1j[0]=deltaR(tmpq1eta,tmpq1phi,eta_max,phi_max);
            dRt_wq2j[0]=deltaR(tmpq2eta,tmpq2phi,eta_max,phi_max);
            int tmpinr=0;
                if(extraqgt[0][0]>0) tmpinr++;
                if(dRt_bj[0]<rbmatching){extraqgt[tmpinr][0]=tmpbpt;extraqgt[tmpinr][1]=tmpbeta;extraqgt[tmpinr][2]=tmpbphi;tmpinr++;}
                if(dRt_wq1j[0]<rbmatching){extraqgt[tmpinr][0]=tmpq1pt;extraqgt[tmpinr][1]=tmpq1eta;extraqgt[tmpinr][2]=tmpq1phi;tmpinr++;}
                if(dRt_wq2j[0]<rbmatching&&tmpt_wtag==4){extraqgt[tmpinr][0]=tmpq2pt;extraqgt[tmpinr][1]=tmpq2eta;extraqgt[tmpinr][2]=tmpq2phi;tmpinr++;}
                int thorl=(tmpt_wtag==4)?1:-1;

            for(int k=0;k<4;k++) {matchnewt[k][0]=tmpinr;matchnewt[k][2]=tmpinr*thorl;
                if(extraqgt[0][0]>0)matchnewt[k][2]=(tmpinr-1)*thorl;
            }
                
                double PTtar=tmptpt;
                const int maxgroupr=tmpinr;
                int pgroupr[4][2*maxgroupr],subgroupr[4][maxgroupr][maxgroupr],nsubgroupr[4][maxgroupr],flagw1[4][2*maxgroupr],flagalone[4][maxgroupr];
                double ptgroupr[4][maxgroupr];
                for(int i=0;i<4;i++){
                    for(int j=0;j<maxgroupr;j++){
                        ptgroupr[i][j]=0;nsubgroupr[i][j]=-99;flagalone[i][j]=-99;
                        for(int k=0;k<maxgroupr;k++) subgroupr[i][j][k]=-99;}
                    for(int j=0;j<2*maxgroupr;j++) {pgroupr[i][j]=-99;flagw1[i][j]=-99;}
                }
                double drqqin5[2*tmpinr];
                for(int i=0;i<2*tmpinr;i++){drqqin5[i]=-99;}
                int idrqqin5=0;
                for(int m=0;m<tmpinr;m++){
                    for(int n=m+1;n<tmpinr;n++){
                        drqqin5[idrqqin5]=deltaR(extraqgt[m][1],extraqgt[m][2],extraqgt[n][1],extraqgt[n][2]);
                        //cout<<m<<"   "<<n<<"   "<<extraqgt[m][1]<<"   "<<extraqgt[m][2]<<"   "<<extraqgt[n][1]<<"   "<<extraqgt[n][2]<<"   "<<drqqin5[idrqqin5]<<endl;
                        idrqqin5++;
                    }
                }
                int i1sub,i2sub;
                double drwp[4]={0.1,0.15,0.2,0.25};
                int end1,end2;
                for(int i=0;i<4;i++){
                    i1sub=0;i2sub=0;
                    for(int j=0;j<2*tmpinr;j++){
                        if(0<drqqin5[j]&&drqqin5[j]<drwp[i]){
                            end1=endpoint(j,tmpinr,1);end2=endpoint(j,tmpinr,2);
                            //cout<<" check  "<<j<<"   "<<end1<<"   "<<end2<<"   "<<flagw1[i][end1]<<"   "<<flagw1[i][end2]<<" i1sub  "<<i1sub<<"   "<<i2sub<<"   "<<endl;
                            if(flagw1[i][end1]==-99&&flagw1[i][end2]==-99) {if(subgroupr[i][i1sub][0]!=-99){i1sub++;i2sub=0;}subgroupr[i][i1sub][i2sub]=end1;subgroupr[i][i1sub][i2sub+1]=end2;
                                flagw1[i][end1]=1;flagw1[i][end2]=1;i2sub+=2;
                                for(int k=0;k<2*tmpinr;k++){
                                    if(0<drqqin5[k]&&drqqin5[k]<drwp[i]){
                                        end1=endpoint(k,tmpinr,1);end2=endpoint(k,tmpinr,2);
                                        if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                            int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                            subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                            for(int m=0;m<2*tmpinr;m++){
                                                if(0<drqqin5[m]&&drqqin5[m]<drwp[i]){
                                                    end1=endpoint(m,tmpinr,1);end2=endpoint(m,tmpinr,2);
                                                    if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                        int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                        subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                        for(int n=0;n<2*tmpinr;n++){
                                                            if(0<drqqin5[n]&&drqqin5[n]<drwp[i]){
                                                                end1=endpoint(n,tmpinr,1);end2=endpoint(n,tmpinr,2);
                                                                if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                                    int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                                    subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                                }}}
                                                    }}}
                                        }
                                    }
                                }
                            }
                        }//end
                    }
                    int totalgroup=0;
                    for(int j=0;j<tmpinr;j++){if(subgroupr[i][j][0]!=-99)totalgroup++;}
                    const int ntotal=totalgroup;
                    int flagw2[ntotal][tmpinr];
                    for(int j=0;j<tmpinr;j++){
                        for(int k=0;k<ntotal;k++){
                            flagw2[k][j]=-99;
                            for(int m=0;m<tmpinr;m++){
                                if(subgroupr[i][k][m]==j){flagw2[k][j]=1;}
                            }}
                        int flagw3=-99;
                        for(int k=0;k<ntotal;k++){
                            if(flagw2[k][j]==1) flagw3=1;}
                        if(flagw3==-99){subgroupr[i][totalgroup][0]=j;totalgroup++;}
                    }
                }
                for(int i=0;i<4;i++){
                    for(int j=0;j<tmpinr;j++){
                        int igroup=0;
                        double tmpsubg[5][2]={0};
                        for(int k=0;k<tmpinr;k++){
                            if(subgroupr[i][j][k]!=-99) {tmpsubg[igroup][0]=extraqgt[subgroupr[i][j][k]][0];tmpsubg[igroup][1]=extraqgt[subgroupr[i][j][k]][2];igroup++;}
                        }
                        ptgroupr[i][j]=SumPt(tmpsubg,igroup);
                        if(ptgroupr[i][j]/PTtar>ignoreprong) {if(matchnewt[i][1]!=-99) matchnewt[i][1]++;if(matchnewt[i][1]==-99) matchnewt[i][1]=1;}
                        /*cout<<"groupt"<<j;
                         for(int k=0;k<tmpinr;k++)cout<<"   "<<subgroupr[i][j][k];
                         cout<<"   "<<ptgroupr[i][j]<<"    "<<ptgroupr[i][j]/PTtar<<"    "<<matchnewt[i][1]<<endl;*/
                    }
                    matchnewt[i][1]=matchnewt[i][1]*thorl;
                }
                }
        }
        }
        if(IsData>0&&(gentop_pt>0||genantitop_pt>0)){
            
            double extraqgw[3][3]={0};
            if (dRtj[0]<rmatching) matchingtmp[0]=-1;
            if (matchingtmp[0]!=-1&&dRantitj[0]<rmatching) matchingtmp[0]=-2;
            if (matchingtmp[0]==-1&&dRantitj[0]<rmatching&&dRtj[0]<dRantitj[0]) matchingtmp[0]=-1;
            if (matchingtmp[0]==-1&&dRantitj[0]<rmatching&&dRtj[0]>dRantitj[0]) matchingtmp[0]=-2;
            //cout<<i<<" jet  "<<Nj8<<"   "<<matchingtmp[0]<<endl;
            if(matchingtmp[0]==-1||matchingtmp[0]==-2){
                if(matchingtmp[0]==-1){tmptpt=gentop_pt;tmpteta=gentop_eta;tmptphi=gentop_phi;tmpbpt=gent_b_pt;tmpbeta=gent_b_eta;tmpbphi=gent_b_phi;tmpwpt=gent_w_pt;tmpweta=gent_w_eta;tmpwphi=gent_w_phi;tmpq1pt=gent_w_q1_pt;tmpq1eta=gent_w_q1_eta;tmpq1phi=gent_w_q1_phi;tmpq2pt=gent_w_q2_pt;tmpq2eta=gent_w_q2_eta;tmpq2phi=gent_w_q2_phi;tmpbm=gent_b_mass;tmpt_wtag=gent_w_tag;}
                if(matchingtmp[0]==-2){tmptpt=genantitop_pt;tmpteta=genantitop_eta;tmptphi=genantitop_phi;tmpbpt=genantit_b_pt;tmpbeta=genantit_b_eta;tmpbphi=genantit_b_phi;tmpwpt=genantit_w_pt;tmpweta=genantit_w_eta;tmpwphi=genantit_w_phi;tmpq1pt=genantit_w_q1_pt;tmpq1eta=genantit_w_q1_eta;tmpq1phi=genantit_w_q1_phi;tmpq2pt=genantit_w_q2_pt;tmpq2eta=genantit_w_q2_eta;tmpq2phi=genantit_w_q2_phi;tmpbm=genantit_b_mass;tmpt_wtag=genantit_w_tag;}
                
                int insideqgr=0;
                for(int j=0;j<10;j++){
                    if(dRgj[j][0]<rbmatching&&ptgengl[j]>0.15*tmpwpt&&ptgengl[j]>extraqgw[0][0]) {extraqgw[0][0]=ptgengl[j];extraqgw[0][1]=etagengl[j];extraqgw[0][2]=phigengl[j];}
                }
                for(int j=0;j<5;j++){
                    if((dRq1j[j][0]<rbmatching&&ptgenq1l[j]>0.15*tmpwpt&&ptgenq1l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq1l[j];extraqgw[0][1]=etagenq1l[j];extraqgw[0][2]=phigenq1l[j];}
                    if((dRq2j[j][0]<rbmatching&&ptgenq2l[j]>0.15*tmpwpt&&ptgenq2l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq2l[j];extraqgw[0][1]=etagenq2l[j];extraqgw[0][2]=phigenq2l[j];}
                    if((dRq3j[j][0]<rbmatching&&ptgenq3l[j]>0.15*tmpwpt&&ptgenq3l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq3l[j];extraqgw[0][1]=etagenq3l[j];extraqgw[0][2]=phigenq3l[j];}
                    if((dRq4j[j][0]<rbmatching&&ptgenq4l[j]>0.15*tmpwpt&&ptgenq4l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq4l[j];extraqgw[0][1]=etagenq4l[j];extraqgw[0][2]=phigenq4l[j];}
                    if((dRq5j[j][0]<rbmatching&&ptgenq5l[j]>0.15*tmpwpt&&ptgenq5l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq5l[j];extraqgw[0][1]=etagenq5l[j];extraqgw[0][2]=phigenq5l[j];}
                }
                dRt_bj[0]=deltaR(tmpbeta,tmpbphi,eta_max,phi_max);
                dRt_wq1j[0]=deltaR(tmpq1eta,tmpq1phi,eta_max,phi_max);
                dRt_wq2j[0]=deltaR(tmpq2eta,tmpq2phi,eta_max,phi_max);
                if(dRt_bj[0]>rbmatching&&tmpt_wtag==4){
                int tmpinr=0;
                if(extraqgw[0][0]>0) tmpinr++;
                if(dRt_wq1j[0]<rbmatching){extraqgw[tmpinr][0]=tmpq1pt;extraqgw[tmpinr][1]=tmpq1eta;extraqgw[tmpinr][2]=tmpq1phi;tmpinr++;}
                if(dRt_wq2j[0]<rbmatching&&tmpt_wtag==4){extraqgw[tmpinr][0]=tmpq2pt;extraqgw[tmpinr][1]=tmpq2eta;extraqgw[tmpinr][2]=tmpq2phi;tmpinr++;}
                
                for(int k=0;k<4;k++) {matchneww[k][0]=tmpinr;matchneww[k][2]=tmpinr;
                    if(extraqgw[0][0]>0)matchneww[k][2]=tmpinr-1;
                }
                int thorl=(tmpt_wtag==4)?1:-1;
                
                    double PTtar=tmpwpt;
                    const int maxgroupr=tmpinr;
                    int pgroupr[4][2*maxgroupr],subgroupr[4][maxgroupr][maxgroupr],nsubgroupr[4][maxgroupr],flagw1[4][2*maxgroupr],flagalone[4][maxgroupr];
                    double ptgroupr[4][maxgroupr];
                    for(int i=0;i<4;i++){
                        for(int j=0;j<maxgroupr;j++){
                            ptgroupr[i][j]=0;nsubgroupr[i][j]=-99;flagalone[i][j]=-99;
                            for(int k=0;k<maxgroupr;k++) subgroupr[i][j][k]=-99;}
                        for(int j=0;j<2*maxgroupr;j++) {pgroupr[i][j]=-99;flagw1[i][j]=-99;}
                    }
                    double drqqin5[2*tmpinr];
                    for(int i=0;i<2*tmpinr;i++){drqqin5[i]=-99;}
                    int idrqqin5=0;
                    for(int m=0;m<tmpinr;m++){
                        for(int n=m+1;n<tmpinr;n++){
                            drqqin5[idrqqin5]=deltaR(extraqgw[m][1],extraqgw[m][2],extraqgw[n][1],extraqgw[n][2]);
                            //cout<<m<<"   "<<n<<"   "<<extraqgw[m][1]<<"   "<<extraqgw[m][2]<<"   "<<extraqgw[n][1]<<"   "<<extraqgw[n][2]<<"   "<<drqqin5[idrqqin5]<<endl;
                            idrqqin5++;
                        }
                    }
                    int i1sub,i2sub;
                    double drwp[4]={0.1,0.15,0.2,0.25};
                    int end1,end2;
                    for(int i=0;i<4;i++){
                        i1sub=0;i2sub=0;
                        for(int j=0;j<2*tmpinr;j++){
                            if(0<drqqin5[j]&&drqqin5[j]<drwp[i]){
                                end1=endpoint(j,tmpinr,1);end2=endpoint(j,tmpinr,2);
                                //cout<<" check  "<<j<<"   "<<end1<<"   "<<end2<<"   "<<flagw1[i][end1]<<"   "<<flagw1[i][end2]<<" i1sub  "<<i1sub<<"   "<<i2sub<<"   "<<endl;
                                if(flagw1[i][end1]==-99&&flagw1[i][end2]==-99) {if(subgroupr[i][i1sub][0]!=-99){i1sub++;i2sub=0;}subgroupr[i][i1sub][i2sub]=end1;subgroupr[i][i1sub][i2sub+1]=end2;
                                    flagw1[i][end1]=1;flagw1[i][end2]=1;i2sub+=2;
                                    for(int k=0;k<2*tmpinr;k++){
                                        if(0<drqqin5[k]&&drqqin5[k]<drwp[i]){
                                            end1=endpoint(k,tmpinr,1);end2=endpoint(k,tmpinr,2);
                                            if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                for(int m=0;m<2*tmpinr;m++){
                                                    if(0<drqqin5[m]&&drqqin5[m]<drwp[i]){
                                                        end1=endpoint(m,tmpinr,1);end2=endpoint(m,tmpinr,2);
                                                        if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                            int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                            subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                            for(int n=0;n<2*tmpinr;n++){
                                                                if(0<drqqin5[n]&&drqqin5[n]<drwp[i]){
                                                                    end1=endpoint(n,tmpinr,1);end2=endpoint(n,tmpinr,2);
                                                                    if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                                        int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                                        subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                                    }}}
                                                        }}}
                                            }
                                        }
                                    }
                                }
                            }//end
                        }
                        int totalgroup=0;
                        for(int j=0;j<tmpinr;j++){if(subgroupr[i][j][0]!=-99)totalgroup++;}
                        const int ntotal=totalgroup;
                        int flagw2[ntotal][tmpinr];
                        for(int j=0;j<tmpinr;j++){
                            for(int k=0;k<ntotal;k++){
                                flagw2[k][j]=-99;
                                for(int m=0;m<tmpinr;m++){
                                    if(subgroupr[i][k][m]==j){flagw2[k][j]=1;}
                                }}
                            int flagw3=-99;
                            for(int k=0;k<ntotal;k++){
                                if(flagw2[k][j]==1) flagw3=1;}
                            if(flagw3==-99){subgroupr[i][totalgroup][0]=j;totalgroup++;}
                        }
                    }
                    for(int i=0;i<4;i++){
                        for(int j=0;j<tmpinr;j++){
                            int igroup=0;
                            double tmpsubg[5][2]={0};
                            for(int k=0;k<tmpinr;k++){
                                if(subgroupr[i][j][k]!=-99) {tmpsubg[igroup][0]=extraqgw[subgroupr[i][j][k]][0];tmpsubg[igroup][1]=extraqgw[subgroupr[i][j][k]][2];igroup++;}
                            }
                            ptgroupr[i][j]=SumPt(tmpsubg,igroup);
                            if(ptgroupr[i][j]/PTtar>ignoreprong) {if(matchneww[i][1]!=-99) matchneww[i][1]++;if(matchneww[i][1]==-99) matchneww[i][1]=1;}
                            /*cout<<"group"<<j;
                             for(int k=0;k<tmpinr;k++)cout<<"   "<<subgroupr[i][j][k];
                             cout<<"   "<<ptgroupr[i][j]<<"    "<<ptgroupr[i][j]/PTtar<<"    "<<matchneww[i][1]<<endl;*/
                        }
                    }

            }
        }
        }
        if(IsData>0&&(gentop_pt<0&&genantitop_pt<0)){
            
            double extraqgw[3][3]={0};
            double allw[5][9],dRqq_wall[5];
            for(int i=0;i<5;i++){
                for(int j=0;j<9;j++) allw[i][j]=-99;
                if(taggenwl[i]==4){
                allw[i][0]=ptgenwl[i];allw[i][1]=etagenwl[i];allw[i][2]=phigenwl[i];
                allw[i][3]=genw_q1_pt[i];allw[i][4]=genw_q1_eta[i];allw[i][5]=genw_q1_phi[i];
                    allw[i][6]=genw_q2_pt[i];allw[i][7]=genw_q2_eta[i];allw[i][8]=genw_q2_phi[i];}
                dRqq_wall[i]=deltaR(allw[i][1],allw[i][1],eta_max,phi_max);
            }
            Int_t *index_wall=new Int_t[5];
            TMath::Sort(5,dRqq_wall,index_wall,0);

                int insideqgr=0;
                for(int j=0;j<10;j++){
                    if(dRgj[j][0]<rbmatching&&ptgengl[j]>0.15*allw[index_wall[0]][0]&&ptgengl[j]>extraqgw[0][0]) {extraqgw[0][0]=ptgengl[j];extraqgw[0][1]=etagengl[j];extraqgw[0][2]=phigengl[j];}
                }
                for(int j=0;j<5;j++){
                    if((dRq1j[j][0]<rbmatching&&ptgenq1l[j]>0.15*allw[index_wall[0]][0]&&ptgenq1l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq1l[j];extraqgw[0][1]=etagenq1l[j];extraqgw[0][2]=phigenq1l[j];}
                    if((dRq2j[j][0]<rbmatching&&ptgenq2l[j]>0.15*allw[index_wall[0]][0]&&ptgenq2l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq2l[j];extraqgw[0][1]=etagenq2l[j];extraqgw[0][2]=phigenq2l[j];}
                    if((dRq3j[j][0]<rbmatching&&ptgenq3l[j]>0.15*allw[index_wall[0]][0]&&ptgenq3l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq3l[j];extraqgw[0][1]=etagenq3l[j];extraqgw[0][2]=phigenq3l[j];}
                    if((dRq4j[j][0]<rbmatching&&ptgenq4l[j]>0.15*allw[index_wall[0]][0]&&ptgenq4l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq4l[j];extraqgw[0][1]=etagenq4l[j];extraqgw[0][2]=phigenq4l[j];}
                    if((dRq5j[j][0]<rbmatching&&ptgenq5l[j]>0.15*allw[index_wall[0]][0]&&ptgenq5l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq5l[j];extraqgw[0][1]=etagenq5l[j];extraqgw[0][2]=phigenq5l[j];}
                }
            double drwq1j=deltaR(allw[index_wall[0]][4],allw[index_wall[0]][5],eta_max,phi_max);
            double drwq2j=deltaR(allw[index_wall[0]][7],allw[index_wall[0]][8],eta_max,phi_max);

                    int tmpinr=0;
                    if(extraqgw[0][0]>0) tmpinr++;
                    if(drwq1j<rbmatching){extraqgw[tmpinr][0]=allw[index_wall[0]][3];extraqgw[tmpinr][1]=allw[index_wall[0]][4];extraqgw[tmpinr][2]=allw[index_wall[0]][5];tmpinr++;}
                    if(drwq2j<rbmatching){extraqgw[tmpinr][0]=allw[index_wall[0]][6];extraqgw[tmpinr][1]=allw[index_wall[0]][7];extraqgw[tmpinr][2]=allw[index_wall[0]][8];tmpinr++;}
                    
                    for(int k=0;k<4;k++) {matchneww[k][0]=tmpinr;matchneww[k][2]=tmpinr;
                        if(extraqgw[0][0]>0)matchneww[k][2]=tmpinr-1;
                    }
                    int thorl=1;
                    
            double PTtar=allw[index_wall[0]][0];
            const int maxgroupr=tmpinr;
            int pgroupr[4][2*maxgroupr],subgroupr[4][maxgroupr][maxgroupr],nsubgroupr[4][maxgroupr],flagw1[4][2*maxgroupr],flagalone[4][maxgroupr];
            double ptgroupr[4][maxgroupr];
            for(int i=0;i<4;i++){
                for(int j=0;j<maxgroupr;j++){
                    ptgroupr[i][j]=0;nsubgroupr[i][j]=-99;flagalone[i][j]=-99;
                    for(int k=0;k<maxgroupr;k++) subgroupr[i][j][k]=-99;}
                for(int j=0;j<2*maxgroupr;j++) {pgroupr[i][j]=-99;flagw1[i][j]=-99;}
            }
            double drqqin5[2*tmpinr];
            for(int i=0;i<2*tmpinr;i++){drqqin5[i]=-99;}
            int idrqqin5=0;
            for(int m=0;m<tmpinr;m++){
                for(int n=m+1;n<tmpinr;n++){
                    drqqin5[idrqqin5]=deltaR(extraqgw[m][1],extraqgw[m][2],extraqgw[n][1],extraqgw[n][2]);
                    //cout<<m<<"   "<<n<<"   "<<extraqgw[m][1]<<"   "<<extraqgw[m][2]<<"   "<<extraqgw[n][1]<<"   "<<extraqgw[n][2]<<"   "<<drqqin5[idrqqin5]<<endl;
                    idrqqin5++;
                }
            }
            int i1sub,i2sub;
            double drwp[4]={0.1,0.15,0.2,0.25};
            int end1,end2;
            for(int i=0;i<4;i++){
                i1sub=0;i2sub=0;
                for(int j=0;j<2*tmpinr;j++){
                    if(0<drqqin5[j]&&drqqin5[j]<drwp[i]){
                        end1=endpoint(j,tmpinr,1);end2=endpoint(j,tmpinr,2);
                        //cout<<" check  "<<j<<"   "<<end1<<"   "<<end2<<"   "<<flagw1[i][end1]<<"   "<<flagw1[i][end2]<<" i1sub  "<<i1sub<<"   "<<i2sub<<"   "<<endl;
                        if(flagw1[i][end1]==-99&&flagw1[i][end2]==-99) {if(subgroupr[i][i1sub][0]!=-99){i1sub++;i2sub=0;}subgroupr[i][i1sub][i2sub]=end1;subgroupr[i][i1sub][i2sub+1]=end2;
                            flagw1[i][end1]=1;flagw1[i][end2]=1;i2sub+=2;
                            for(int k=0;k<2*tmpinr;k++){
                                if(0<drqqin5[k]&&drqqin5[k]<drwp[i]){
                                    end1=endpoint(k,tmpinr,1);end2=endpoint(k,tmpinr,2);
                                    if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                        int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                        subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                        for(int m=0;m<2*tmpinr;m++){
                                            if(0<drqqin5[m]&&drqqin5[m]<drwp[i]){
                                                end1=endpoint(m,tmpinr,1);end2=endpoint(m,tmpinr,2);
                                                if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                    int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                    subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                    for(int n=0;n<2*tmpinr;n++){
                                                        if(0<drqqin5[n]&&drqqin5[n]<drwp[i]){
                                                            end1=endpoint(n,tmpinr,1);end2=endpoint(n,tmpinr,2);
                                                            if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                                int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                                subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                            }}}
                                                }}}
                                    }
                                }
                            }
                        }
                    }//end
                }
                int totalgroup=0;
                for(int j=0;j<tmpinr;j++){if(subgroupr[i][j][0]!=-99)totalgroup++;}
                const int ntotal=totalgroup;
                int flagw2[ntotal][tmpinr];
                for(int j=0;j<tmpinr;j++){
                    for(int k=0;k<ntotal;k++){
                        flagw2[k][j]=-99;
                        for(int m=0;m<tmpinr;m++){
                            if(subgroupr[i][k][m]==j){flagw2[k][j]=1;}
                        }}
                    int flagw3=-99;
                    for(int k=0;k<ntotal;k++){
                        if(flagw2[k][j]==1) flagw3=1;}
                    if(flagw3==-99){subgroupr[i][totalgroup][0]=j;totalgroup++;}
                }
            }
            for(int i=0;i<4;i++){
                for(int j=0;j<tmpinr;j++){
                    int igroup=0;
                    double tmpsubg[5][2]={0};
                    for(int k=0;k<tmpinr;k++){
                        if(subgroupr[i][j][k]!=-99) {tmpsubg[igroup][0]=extraqgw[subgroupr[i][j][k]][0];tmpsubg[igroup][1]=extraqgw[subgroupr[i][j][k]][2];igroup++;}
                    }
                    ptgroupr[i][j]=SumPt(tmpsubg,igroup);
                    if(ptgroupr[i][j]/PTtar>ignoreprong) {if(matchneww[i][1]!=-99) matchneww[i][1]++;if(matchneww[i][1]==-99) matchneww[i][1]=1;}
                    /*cout<<"group"<<j;
                     for(int k=0;k<tmpinr;k++)cout<<"   "<<subgroupr[i][j][k];
                     cout<<"   "<<ptgroupr[i][j]<<"    "<<ptgroupr[i][j]/PTtar<<"    "<<matchneww[i][1]<<endl;*/
                }
            }
                    }

  
}

#endif