#ifndef _endpoint_
#define _endpoint_

int endpoint(int nl,int np,int p1or2){
    int p1,p2;
    if(np==5){
        if(nl==0){p1=0;p2=1;}
        if(nl==1){p1=0;p2=2;}
        if(nl==2){p1=0;p2=3;}
        if(nl==3){p1=0;p2=4;}
        if(nl==4){p1=1;p2=2;}
        if(nl==5){p1=1;p2=3;}
        if(nl==6){p1=1;p2=4;}
        if(nl==7){p1=2;p2=3;}
        if(nl==8){p1=2;p2=4;}
        if(nl==9){p1=3;p2=4;}
        }
    if(np==4){
        if(nl==0){p1=0;p2=1;}
        if(nl==1){p1=0;p2=2;}
        if(nl==2){p1=0;p2=3;}
        if(nl==3){p1=1;p2=2;}
        if(nl==4){p1=1;p2=3;}
        if(nl==5){p1=2;p2=3;}
    }
    if(np==3){
        if(nl==0){p1=0;p2=1;}
        if(nl==1){p1=0;p2=2;}
        if(nl==2){p1=1;p2=2;}
    }
    if(np==2){
        if(nl==0){p1=0;p2=1;}
    }
    if(p1or2==2) p1=p2;
    return p1;
}

#endif
