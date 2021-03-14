#ifndef _num_of_AK8_AK4_
#define _num_of_AK8_AK4_

void EDBR2PKUTree::num_of_AK8_AK4() {   
        //    Nj8_2=1,Nj8=0;
           Nj8_2=1,Nj8=1;
        
        if(jet_pt_puppi_2>0&&abs(jetAK8puppi_eta_2)<2.4){
            Nj8=2;}
        
        // Nj4=0;
        // for(Int_t ii=0;ii<8;ii++){
        //     if(ak4jet_pt[ii]>0){
        //         Nj4++;
        //     }
        // }
        // for(Int_t ii=0;ii<4;ii++){
        //     if(jetAK8puppi_pt1[ii]>0&&abs(jetAK8puppi_eta1[ii])<2.4){
        //         Nj8++;
        //     }
        // }
}

#endif