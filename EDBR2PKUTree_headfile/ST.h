#ifndef _ST_cal_
#define _ST_cal_

void EDBR2PKUTree::ST_cal() {   
        if(Nj8==2)
            ST=W_pt+jet_pt_puppi+jet_pt_puppi_2;
        if(Nj8==1)
            ST=W_pt+jet_pt_puppi;
        if(Nj8>=3)
            ST=W_pt+jet_pt_puppi+jet_pt_puppi_2+jet_pt_puppi_3;
}

#endif