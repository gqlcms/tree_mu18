#ifndef _pt_order_consider_w_
#define _pt_order_consider_w_

void EDBR2PKUTree::pt_order_consider_w() {   
           Int_t *indexx11=new Int_t[3];
        double ptwdivw[3]={-99,-99,-99};
        ptwdivw[0]=W_pt;ptwdivw[1]=jet_pt_puppi;ptwdivw[2]=jet_pt_puppi_2;
        TMath::Sort(3,ptwdivw,indexx11,1);
        PTw2divPTw1=ptwdivw[indexx11[1]]/ptwdivw[indexx11[0]];
        PTw3divPTw1=ptwdivw[indexx11[2]]/ptwdivw[indexx11[0]];
        PTw3divPTw2=ptwdivw[indexx11[2]]/ptwdivw[indexx11[1]];
        PTw[0]=ptwdivw[indexx11[0]];PTw[1]=ptwdivw[indexx11[1]];PTw[2]=ptwdivw[indexx11[2]];
        if(Nj8==1)
        {PTw3divPTw1=-99;PTw3divPTw2=-99;}
}

#endif