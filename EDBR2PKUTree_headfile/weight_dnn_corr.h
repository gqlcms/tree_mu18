#ifndef _weight_dnn_corr_
#define _weight_dnn_corr_

void EDBR2PKUTree::weight_dnn_corr(Int_t IsData,Double_t massdnn[5],TH1D *h_t60,TH1D *h_q60,TH1D *h_w60,TH1D *h_t100,TH1D *h_q100,TH1D *h_w100,TH1D *h_t150,TH1D *h_q150,TH1D *h_t200,TH1D *h_q200) {   
   	double dnncorr=1;//useless
	double dnncorr_all=1;        
	double dnncorr_nosft60=1;
        double dnncorr_nosfq60=1;
        double dnncorr_nosfw60=1;
        double dnncorr_nosft100=1;
        double dnncorr_nosfq100=1;
        double dnncorr_nosfw100=1;
        double dnncorr_nosft150=1;
        double dnncorr_nosfq150=1;
        double dnncorr_nosft200=1;
        double dnncorr_nosfq200=1;
	
	double dnn_corrt60[2]={1,1},dnn_corrq60[2]={1,1},dnn_corrw60[2]={1,1},dnn_corrt100[2]={1,1},dnn_corrq100[2]={1,1},dnn_corrw100[2]={1,1},dnn_corrt150[2]={1,1},dnn_corrq150[2]={1,1},dnn_corrt200[2]={1,1},dnn_corrq200[2]={1,1};
	
       if(IsData!=0){
            Double_t jet_mass[2]={Mj_max,Mj_min};
	    for(int inj8=0;inj8<Nj8&&inj8<3&&inj8<2;inj8++){
            	int ixd1=0;
	    	if(inj8==1&&Nj8==3) jet_mass[inj8]=Mj_mid;
            	double dnntagd1=dnndecw_max;
		if(inj8==1&&Nj8==2) dnntagd1=dnndecw_min;
		if(inj8==1&&Nj8==3) dnntagd1=dnndecw_mid;
            	if(jet_mass[inj8]>=100) dnntagd1=(dnndecraww_max+dnndecrawh4q_max)/(dnndecrawqcd_max+dnndecrawh4q_max+dnndecraww_max);
            	if(jet_mass[inj8]>=100&&inj8==1&&Nj8==2) dnntagd1=(dnndecraww_min+dnndecrawh4q_min)/(dnndecrawqcd_min+dnndecrawh4q_min+dnndecraww_min);
            	if(jet_mass[inj8]>=100&&inj8==1&&Nj8==3) dnntagd1=(dnndecraww_mid+dnndecrawh4q_mid)/(dnndecrawqcd_mid+dnndecrawh4q_mid+dnndecraww_mid);
            	for(ixd1=0;ixd1<4;ixd1++){
            	    if(jet_mass[inj8]<massdnn[0]) {ixd1=-1;break;}
            	    if(jet_mass[inj8]>=massdnn[ixd1]&&jet_mass[inj8]<massdnn[ixd1+1]) break;
            	}
            	bool tmatching1=(matchingt[inj8]>0)||(matchingr[inj8]==3 &&( (status_2==4&&status_3==4)))||(matchingr[inj8]>=4);
		bool wmatching1=(matchingr[inj8]<=0 &&matchingt[inj8]<=0 && matchingw[inj8]>0)||(matchingr[inj8]==2 &&(! (status_2==4&&status_3==4)));
		//bool wmatching1=(matchingr[inj8]<=0 &&matchingt[inj8]<=0 && matchingw[inj8]>0);
		bool qmatching1=matchingr[inj8]<=0 &&matchingt[inj8]<=0 && matchingw[inj8]<=0 &&(matchingg[inj8]>0 || matchingq[inj8]>0);
                if(dnntagd1>=0 && dnntagd1<=1){
		if(dnntagd1==1)dnntagd1=0.99;
            	if(ixd1==0){
			if(tmatching1) dnn_corrt60[inj8]=h_t60->GetBinContent(h_t60->FindBin(dnntagd1));
			if(wmatching1) dnn_corrw60[inj8]=h_w60->GetBinContent(h_w60->FindBin(dnntagd1));
			if(qmatching1) dnn_corrq60[inj8]=h_q60->GetBinContent(h_q60->FindBin(dnntagd1));
		}
            	if(ixd1==1){
			if(tmatching1) dnn_corrt100[inj8]=h_t100->GetBinContent(h_t100->FindBin(dnntagd1));
			if(wmatching1) dnn_corrw100[inj8]=h_w100->GetBinContent(h_w100->FindBin(dnntagd1));
			if(qmatching1) dnn_corrq100[inj8]=h_q100->GetBinContent(h_q100->FindBin(dnntagd1));
		}
            	if(ixd1==2){
			if(wmatching1&&IsData==10) dnn_corrw100[inj8]=h_w100->GetBinContent(h_w100->FindBin(dnntagd1));
			if(tmatching1) dnn_corrt150[inj8]=h_t150->GetBinContent(h_t150->FindBin(dnntagd1));
			if(qmatching1) dnn_corrq150[inj8]=h_q150->GetBinContent(h_q150->FindBin(dnntagd1));
		}
            	if(ixd1==3){
			if(wmatching1&&IsData==10) dnn_corrw100[inj8]=h_w100->GetBinContent(h_w100->FindBin(dnntagd1));
			if(tmatching1) dnn_corrt200[inj8]=h_t200->GetBinContent(h_t200->FindBin(dnntagd1));
			if(qmatching1) dnn_corrq200[inj8]=h_q200->GetBinContent(h_q200->FindBin(dnntagd1));
		}}
		dnncorr_all*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8];
		dnncorr_nosft60*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrt60[inj8];
		dnncorr_nosfq60*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrq60[inj8];
		dnncorr_nosfw60*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrw60[inj8];
		dnncorr_nosft100*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrt100[inj8];
		dnncorr_nosfq100*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrq100[inj8];
		dnncorr_nosfw100*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrw100[inj8];
		dnncorr_nosft150*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrt150[inj8];
		dnncorr_nosfq150*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrq150[inj8];
		dnncorr_nosft200*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrt200[inj8];
		dnncorr_nosfq200*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrq200[inj8];
		//if(inj8==0) cout<<" massdnn[ixd1]  "<<massdnn[ixd1]<<" Mj_max  "<<Mj_max<<" dnndecw_max "<<dnndecw_max<<" wh "<<(dnndecraww_max+dnndecrawh4q_max)/(dnndecrawqcd_max+dnndecrawh4q_max+dnndecraww_max)<<" dnntagd1 "<<dnntagd1<<" wmatching1  "<<wmatching1<<" tmatching1 "<<tmatching1<<" qmatching1 "<<qmatching1<<" sfall "<<dnncorr_all<<endl;
		//if(inj8==1) cout<<" massdnn[ixd1]  "<<massdnn[ixd1]<<" Mj_max  "<<Mj_max<<" dnndecw_max "<<dnndecw_max<<" wh "<<(dnndecraww_max+dnndecrawh4q_max)/(dnndecrawqcd_max+dnndecrawh4q_max+dnndecraww_max)<<" dnntagd1 "<<dnntagd1<<" wmatching1  "<<wmatching1<<" tmatching1 "<<tmatching1<<" qmatching1 "<<qmatching1<<" sfall "<<dnncorr_all<<endl;
	     }
	}
	//cout<<" inj8  "<<Nj8<<" dnncorr_all "<<dnncorr_all<<endl;
	//	if(jentry>=200) break;
        weight_deep_matchedcorr1=weight_deep*dnncorr;
        weight_deep_matchedcorr2=weight_deep*dnncorr*dnncorr;
        weight_deep_tig_matchedcorr1=weight_tig_deep*dnncorr;
        weight_deep_tig_matchedcorr2=weight_tig_deep*dnncorr*dnncorr;
        weight_deep_los_matchedcorr1=weight_los_deep*dnncorr;
        weight_deep_los_matchedcorr2=weight_los_deep*dnncorr*dnncorr;
	//cout<<dnncorr_all<<" dnncorr_all "<<dnncorr_nosft60<<"  "<<dnncorr_nosfq60<<"  "<<dnncorr_nosfw60<<"  "<<dnncorr_nosft100<<"  "<<dnncorr_nosfq100<<"  "<<dnncorr_nosfw100<<"  "<<dnncorr_nosft150<<"  "<<dnncorr_nosfq150<<"  "<<dnncorr_nosft200<<"  "<<dnncorr_nosfq200<<"  "<<endl;
        weight_deep_lumi_corr1=weight_deep_lumi*dnncorr_all;
	weight_deep_lumi_nosft60=weight_deep_lumi*dnncorr_nosft60;
	weight_deep_lumi_nosfq60=weight_deep_lumi*dnncorr_nosfq60;
	weight_deep_lumi_nosfw60=weight_deep_lumi*dnncorr_nosfw60;
	weight_deep_lumi_nosft100=weight_deep_lumi*dnncorr_nosft100;
	weight_deep_lumi_nosfq100=weight_deep_lumi*dnncorr_nosfq100;
	weight_deep_lumi_nosfw100=weight_deep_lumi*dnncorr_nosfw100;
	weight_deep_lumi_nosft150=weight_deep_lumi*dnncorr_nosft150;
	weight_deep_lumi_nosfq150=weight_deep_lumi*dnncorr_nosfq150;
	weight_deep_lumi_nosft200=weight_deep_lumi*dnncorr_nosft200;
	weight_deep_lumi_nosfq200=weight_deep_lumi*dnncorr_nosfq200;
        weight_deep_lumi_corr2=weight_deep_lumi*dnncorr_all*dnncorr_all;
}

#endif