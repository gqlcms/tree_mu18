#ifndef _SumPt_
#define _SumPt_

double SumPt(double ptphi[5][2],int n){
    double pt=-99;
    double x1=0,x2=0;
    for(int i=0;i<n;i++){x1+=ptphi[i][0]*TMath::Cos(ptphi[i][1]);x2+=ptphi[i][0]*TMath::Sin(ptphi[i][1]);}
    //cout<<"  sum pt "<<ptphi[0][0]<<"   "<<x1<<"   "<<x2<<endl;
    if(n>=1)pt=sqrt(x1*x1+x2*x2);
    return pt;
}

#endif