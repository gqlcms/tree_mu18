#ifndef _Nsubjetness_
#define _Nsubjetness_

void EDBR2PKUTree::Nsubjetness() {   

        jetAK8puppi_tau31=jetAK8puppi_tau3/jetAK8puppi_tau1;
        jetAK8puppi_tau32=jetAK8puppi_tau3/jetAK8puppi_tau2;
        jetAK8puppi_tau43=jetAK8puppi_tau4/jetAK8puppi_tau3;
        jetAK8puppi_tau41=jetAK8puppi_tau4/jetAK8puppi_tau1;
        jetAK8puppi_tau42=jetAK8puppi_tau4/jetAK8puppi_tau2;

        jetAK8puppi_tau31_2=jetAK8puppi_tau3_2/jetAK8puppi_tau1_2;
        jetAK8puppi_tau32_2=jetAK8puppi_tau3_2/jetAK8puppi_tau2_2;
        jetAK8puppi_tau43_2=jetAK8puppi_tau4_2/jetAK8puppi_tau3_2;
        jetAK8puppi_tau41_2=jetAK8puppi_tau4_2/jetAK8puppi_tau1_2;
        jetAK8puppi_tau42_2=jetAK8puppi_tau4_2/jetAK8puppi_tau2_2;

        jetAK8puppi_tau31_3=jetAK8puppi_tau3_3/jetAK8puppi_tau1_3;
        jetAK8puppi_tau32_3=jetAK8puppi_tau3_3/jetAK8puppi_tau2_3;
        jetAK8puppi_tau43_3=jetAK8puppi_tau4_3/jetAK8puppi_tau3_3;
        jetAK8puppi_tau41_3=jetAK8puppi_tau4_3/jetAK8puppi_tau1_3;
        jetAK8puppi_tau42_3=jetAK8puppi_tau4_3/jetAK8puppi_tau2_3;
 
   
}
void EDBR2PKUTree::Nsubjetness2() {   
            if(Nj8==1)
        {jetAK8puppi_tau31_2=-99;jetAK8puppi_tau32_2=-99;jetAK8puppi_tau43_2=-99;jetAK8puppi_tau41_2=-99;jetAK8puppi_tau42_2=-99;jetAK8puppi_tau31_3=-99;jetAK8puppi_tau32_3=-99;jetAK8puppi_tau43_3=-99;jetAK8puppi_tau41_3=-99;jetAK8puppi_tau42_3=-99;}
        if(Nj8==2)
        {jetAK8puppi_tau32_3=-99;jetAK8puppi_tau43_3=-99;jetAK8puppi_tau41_3=-99;jetAK8puppi_tau42_3=-99;}
}

#endif