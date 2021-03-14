#ifndef _weight_tau_or_deep_corr_
#define _weight_tau_or_deep_corr_

void EDBR2PKUTree::weight_tau_or_deep_corr(TString channelname, Int_t IsData, double weight_tmp) {   
           Double_t tau_corr1=1,tau_corr2=1;
 
        Double_t corr_tau_top[8][8]={0};
    if( channelname=="mu" ){
        Double_t corr_tau_top_tmp[8][8]={1.03402672187,1.13558084027,1.00360681557,1.00246285777,0.99502594773,0.941244969132,0.873699375981,1.31334378138,
            1.27402792202,1.20715713692,1.02959119179,0.899276775573,1.02657444441,0.884478365492,0.841185385275,0.728595008922,
            1.05702477048,1.12832937166,0.895318658766,1.23145431738,0.842625831911,0.750238436967,1.09420467154,1.91003197061,
            0.980603593355,1.02947129926,1.129511571,0.949527312819,0.868362761076,0.951283339336,0.743157927249,1.12065620372,
            1.08967172869,1.11645005347,0.993568975319,0.993556032941,0.955014915484,0.727630260733,0.855622724783,0.737436065466,
            0.402541005911,1.38553543957,1.17279430632,1.0355029671,1.04503566577,0.991063708437,0.888307528307,0.762422249715,
            1.14159534938,1.05231203514,1.02810296348,1.05506976435,0.952612431097,0.983597829862,0.811889597808,0.543743762695,
            3.38517432139,1.16865887847,1.04707550932,0.967030785779,0.970556964857,1.0061406186,0.870983443008,0.532607396537};
        for(int i=0;i<8;i++){for(int j=0;j<8;j++) corr_tau_top[i][j]=corr_tau_top_tmp[i][j];}
    }
        if( channelname=="el" ){
            Double_t corr_tau_top_tmp[8][8]={1.41025975583,1.15858307886,1.0452690859,1.07341169412,0.971667001311,0.761359522858,0.968412543795,1.11905262524,
                1.06615963681,0.991362674834,1.0832466723,1.00584033907,0.921137288222,0.911435773944,1.16787220624,1.16606543631,
                1.01413792196,1.08350052267,1.05526182615,0.884223867405,0.999574049646,0.956885484244,0.993222583722,0.624462847374,
                1.07504405837,1.10128884231,0.96863454652,1.03623323744,0.985294528506,0.740692055289,0.724261698936,0.285626769713,
                0.991639809344,1.08779087537,0.999179928681,0.989316527303,0.979710236647,0.868862648216,0.812155977336,0.0,
                0.904798527337,1.07531968905,1.06744039711,1.02579027671,0.961549277631,1.07560869706,0.936769114354,0.791346031794,
                0.670408713531,1.06802236901,1.08949541824,0.964228451804,0.905098361023,1.05957081004,0.962326062293,0.414373981633,
                0.0,1.36172365544,0.862197044607,1.09692614851,1.06289307319,0.68202643517,1.10112065357,1.73629657039};
            for(int i=0;i<8;i++){for(int j=0;j<8;j++) corr_tau_top[i][j]=corr_tau_top_tmp[i][j];}
        }

        if(IsData==1||IsData==9||IsData==11||IsData==10){
        Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
            Double_t taubin[9]={0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};
            Double_t massbin[9]={50,60,70,80,90,100,150,200,1000000};
            int ix=0;
            for(ix=0;ix<8;ix++){
                if(jet_mass1<massbin[0]) {ix=-1;break;}
                if(jet_mass1>=massbin[ix]&&jet_mass1<massbin[ix+1]) break;
                }
            int iy=0;
            for(iy=0;iy<8;iy++){
                if(jet_mass1<100&&(jetAK8puppi_tau21<0||jetAK8puppi_tau21>1)) {iy=-1;break;}
                if(jet_mass1>=100&&(jetAK8puppi_tau42<0||jetAK8puppi_tau42>1)) {iy=-1;break;}
                if(jet_mass1<100&&(jetAK8puppi_tau21>=taubin[iy]&&jetAK8puppi_tau21<taubin[iy+1])) break;
                if(jet_mass1>=100&&(jetAK8puppi_tau42>=taubin[iy]&&jetAK8puppi_tau42<taubin[iy+1])) break;
            }
            //cout<<ix<<iy<<endl;
            if(ix>=0&&iy>=0) tau_corr1=corr_tau_top[ix][iy];
            int ix2=0;
            for(ix2=0;ix2<8;ix2++){
                if(jet_mass2<massbin[0]) {ix2=-1;break;}
                if(jet_mass2>=massbin[ix2]&&jet_mass2<massbin[ix2+1]) break;
            }
            int iy2=0;
            for(iy2=0;iy2<8;iy2++){
                if(jet_mass2<100&&(jetAK8puppi_tau21_2<0||jetAK8puppi_tau21_2>1)) {iy2=-1;break;}
                if(jet_mass2>=100&&(jetAK8puppi_tau42_2<0||jetAK8puppi_tau42_2>1)) {iy2=-1;break;}
                if(jet_mass2<100&&(jetAK8puppi_tau21_2>=taubin[iy2]&&jetAK8puppi_tau21_2<taubin[iy2+1])) break;
                if(jet_mass2>=100&&(jetAK8puppi_tau42_2>=taubin[iy2]&&jetAK8puppi_tau42_2<taubin[iy2+1])) break;
            }
            if(ix2>=0&&iy2>=0) tau_corr2=corr_tau_top[ix2][iy2];
            //if(ix==1||ix2==1) cout<<jet_mass1<<"  "<<jet_mass2<<"  "<<ix<<"   "<<iy<<"   "<<ix2<<"   "<<iy2<<"  "<<tau_corr1<<"  "<<tau_corr2<<endl;
        }
        if(tau_corr1==0) tau_corr1=1;
        if(tau_corr2==0) tau_corr2=1;
        //cout<<tau_corr1<<endl;
        
        weight_tautopcorr1=weight*tau_corr1*tau_corr2;
        weight_tautopcorr2=weight*tau_corr1*tau_corr2*tau_corr1*tau_corr2;
        
        Double_t corr_tau_wjets[8][8]={0};
        if( channelname=="mu" ){
            
            Double_t corr_tau_wjets_tmp[8][8]={0.878093740102,1.35156804567,1.21240539342,1.04662781583,0.994611210851,0.810461996416,0.7266833886,0.815895651623,
                1.09906271047,1.14509245011,1.02713684287,1.14865862365,0.926744533372,0.875857159651,0.726745853445,0.706692247635,
                1.08368794285,1.07212760852,1.16308231664,0.967811047954,0.973180942193,0.89205332303,0.63844271766,0.681519372528,
                1.04284648379,1.11675705123,1.05408987696,1.04076524132,0.983030478416,0.723921597151,0.767649447996,0.326313172076,
                1.07258152449,1.11993999688,1.07836471575,1.0553622067,0.885342143572,0.832897953632,0.492360188901,0.525759746067,
                4.26894599699,0.957789410667,1.07579486091,0.989669333318,0.968463342718,1.03636625932,0.986605654291,0.924427560405,
                1.06961238475,0.961331819599,0.963633250292,0.999552852707,0.925391876451,1.08318790875,1.11442288616,0.869397442857,
                0.738869452843,1.15491886733,0.89828967094,0.846040089658,0.913241660166,1.17428972598,1.24037277721,1.10446103014};
            for(int i=0;i<8;i++){for(int j=0;j<8;j++) corr_tau_wjets[i][j]=corr_tau_wjets_tmp[i][j];}
        }
        if( channelname=="el" ){
            Double_t corr_tau_wjets_tmp[8][8]={1.17292641723,1.10891790676,1.12095389311,0.973386080896,1.05078999779,0.905922788641,0.836077714046,0.858081094185,
                1.28701282261,1.14259715244,1.14058794551,1.05418747543,0.996798379208,0.817133813464,0.666248461507,0.787187559044,
                1.13685253873,1.10369052406,1.10811962364,1.1022258448,0.938749995542,0.750294624097,0.699666712499,0.71074161622,
                0.97332603372,1.03820385269,1.16188763707,1.02336152663,0.879555321097,0.895620009769,0.67726050421,0.178077617303,
                0.559131233724,0.651974135154,0.924824246403,1.20979144196,1.61383950661,1.71306880469,2.04605542934,0.679655808092,
                0.859778527345,1.08807037348,0.997113874532,0.991518733497,1.05233073551,1.06841187914,0.911995001616,0.754012650151,
                1.96017343947,0.935163911645,0.920617472488,0.931746105761,1.06198942512,1.01536409014,1.11621738285,1.326943875,
                1.90352970457,0.919726446215,0.934034407982,0.951325045954,0.89204333586,1.26855244267,0.984714168594,1.11868970913};
            for(int i=0;i<8;i++){for(int j=0;j<8;j++) corr_tau_wjets[i][j]=corr_tau_wjets_tmp[i][j];}
        }
        
        if((IsData>=2&&IsData<=8)){
            Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
            Double_t taubin[9]={0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};
            Double_t massbin[9]={50,60,70,80,90,100,150,200,1000000};
            int ix=0;
            for(ix=0;ix<8;ix++){
                if(jet_mass1<massbin[0]) {ix=-1;break;}
                if(jet_mass1>=massbin[ix]&&jet_mass1<massbin[ix+1]) break;
            }
            int iy=0;
            for(iy=0;iy<8;iy++){
                if(jet_mass1<100&&(jetAK8puppi_tau21<0||jetAK8puppi_tau21>1)) {iy=-1;break;}
                if(jet_mass1>=100&&(jetAK8puppi_tau42<0||jetAK8puppi_tau42>1)) {iy=-1;break;}
                if(jet_mass1<100&&(jetAK8puppi_tau21>=taubin[iy]&&jetAK8puppi_tau21<taubin[iy+1])) break;
                if(jet_mass1>=100&&(jetAK8puppi_tau42>=taubin[iy]&&jetAK8puppi_tau42<taubin[iy+1])) break;
            }
            //cout<<ix<<iy<<endl;
            if(ix>=0&&iy>=0) tau_corr1=corr_tau_wjets[ix][iy];
            int ix2=0;
            for(ix2=0;ix2<8;ix2++){
                if(jet_mass2<massbin[0]) {ix2=-1;break;}
                if(jet_mass2>=massbin[ix2]&&jet_mass2<massbin[ix2+1]) break;
            }
            int iy2=0;
            for(iy2=0;iy2<8;iy2++){
                if(jet_mass2<100&&(jetAK8puppi_tau21_2<0||jetAK8puppi_tau21_2>1)) {iy2=-1;break;}
                if(jet_mass2>=100&&(jetAK8puppi_tau42_2<0||jetAK8puppi_tau42_2>1)) {iy2=-1;break;}
                if(jet_mass2<100&&(jetAK8puppi_tau21_2>=taubin[iy2]&&jetAK8puppi_tau21_2<taubin[iy2+1])) break;
                if(jet_mass2>=100&&(jetAK8puppi_tau42_2>=taubin[iy2]&&jetAK8puppi_tau42_2<taubin[iy2+1])) break;
            }
            if(ix2>=0&&iy2>=0) tau_corr2=corr_tau_wjets[ix2][iy2];
            
            
            if(tau_corr1==0) tau_corr1=1;
            if(tau_corr2==0) tau_corr2=1;
            //cout<<jet_mass1<<" a "<<jet_mass_puppi_corr<<"  "<<ix<<"  "<<jetAK8puppi_tau21<<"  "<<iy<<"  "<<tau_corr1<<endl;
            //cout<<jet_mass2<<" b "<<jet_mass_puppi_corr_2<<"  "<<ix2<<"  "<<jetAK8puppi_tau21_2<<"  "<<iy2<<"  "<<tau_corr2<<endl;
            //cout<<tau_corr1<<endl;
            weight_tautopcorr1=weight*tau_corr1*tau_corr2;
            weight_tautopcorr2=weight*tau_corr1*tau_corr2*tau_corr1*tau_corr2;
        }

        ///***** corr deepak8 ******///
        Double_t deepak8_corr1=1,deepak8_corr2=1;

        Double_t corr_deepak8_top[4][12]={0};
        if( channelname=="mu"||channelname=="el" ){
            Double_t corr_deepak8_top_tmp[4][12]={0.644693981106,0.813884479014,0.9581094444,1.10369465581,1.06317264272,1.11860474771,1.01341532195,1.11330373755,1.16134793033,1.00208593929,1.15777162098,0.971844322721,
                0.709693803338,0.796987633814,1.09747903963,1.10616039205,0.981750115444,0.880359074993,1.12636126636,1.13281039152,1.21878327287,1.1731795458,1.17219072371,1.00616642582,
                0.694304423163,0.93189847275,1.06798882602,0.917641891085,1.01675453078,0.987941134099,1.00122560345,1.02471007355,1.0754614055,1.0754614055,1.05078356126,1.05078356126,
                0.947338757927,0.947338757927,1.04751871879,1.04751871879,0.906131230598,0.906131230598,1.04013842851,1.04013842851,1.0890765813,1.0890765813,1.0890765813,1.0890765813};
            for(int i=0;i<4;i++){for(int j=0;j<12;j++) corr_deepak8_top[i][j]=corr_deepak8_top_tmp[i][j];}
        }
        if(IsData==1||IsData==9||IsData==11||IsData==10){
            Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
            Double_t deepak8bin[13]={0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
            Double_t massbin[5]={60,100,150,200,1000000};
            int ixx=0;
            for(ixx=0;ixx<4;ixx++){
                if(jet_mass1<massbin[0]) {ixx=-1;break;}
                if(jet_mass1>=massbin[ixx]&&jet_mass1<massbin[ixx+1]) break;
            }
            int iyy=0;
            for(iyy=0;iyy<12;iyy++){
                if(jet_mass1<100&&(jetAK8puppi_dnnDecorrW<0||jetAK8puppi_dnnDecorrW>1)) {iyy=-1;break;}
                if(jet_mass1>=100&&(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))<0||((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>1)) {iyy=-1;break;}
                if(jet_mass1<100&&(jetAK8puppi_dnnDecorrW>=deepak8bin[iyy]&&jetAK8puppi_dnnDecorrW<deepak8bin[iyy+1])) break;
                if(jet_mass1>=100&&(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>=deepak8bin[iyy]&&((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))<deepak8bin[iyy+1])) break;
            }
            //cout<<ixx<<iyy<<endl;
            if(ixx>=0&&iyy>=0) deepak8_corr1=corr_deepak8_top[ixx][iyy];
            int ixx2=0;
            for(ixx2=0;ixx2<4;ixx2++){
                if(jet_mass2<massbin[0]) {ixx2=-1;break;}
                if(jet_mass2>=massbin[ixx2]&&jet_mass2<massbin[ixx2+1]) break;
            }
            int iyy2=0;
            for(iyy2=0;iyy2<12;iyy2++){
                if(jet_mass2<100&&(jetAK8puppi_dnnDecorrW_2<0||jetAK8puppi_dnnDecorrW_2>1)) {iyy2=-1;break;}
                if(jet_mass2>=100&&(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))<0||((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))>1)) {iyy2=-1;break;}
                if(jet_mass2<100&&(jetAK8puppi_dnnDecorrW_2>=deepak8bin[iyy2]&&jetAK8puppi_dnnDecorrW_2<deepak8bin[iyy2+1])) break;
                if(jet_mass2>=100&&(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))>=deepak8bin[iyy2]&&((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))<deepak8bin[iyy2+1])) break;
            }
            if(ixx2>=0&&iyy2>=0) deepak8_corr2=corr_deepak8_top[ixx2][iyy2];
            //if(ixx==4&&(iyy==4)) cout<<jet_mass1<<"  "<<jetAK8puppi_dnnDecorrW<<"  "<<(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd)))<<"  "<<jet_mass2<<"  "<<ixx<<"   "<<iyy<<"   "<<ixx2<<"   "<<iyy2<<"  "<<deepak8_corr1<<"  "<<deepak8_corr2<<endl;
        
        if(deepak8_corr1==0) deepak8_corr1=1;
        if(deepak8_corr2==0) deepak8_corr2=1;

        weight_deep=btagweight_center_deep*weight_tmp;
        weight_deep_topcorr1=weight_deep*deepak8_corr1*deepak8_corr2;
        weight_deep_topcorr2=weight_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;

            weight_tig_deep=btagweight_center_tig_deep*weight_tmp;
            weight_tig_deep_topcorr1=weight_tig_deep*deepak8_corr1*deepak8_corr2;
            weight_tig_deep_topcorr2=weight_tig_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;
            weight_los_deep=btagweight_center_los_deep*weight_tmp;
            weight_los_deep_topcorr1=weight_los_deep*deepak8_corr1*deepak8_corr2;
            weight_los_deep_topcorr2=weight_los_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;

        //if(jet_mass_puppi_corr>=100&&jet_mass_puppi_corr<150&&((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>=0.8) cout<<deepak8_corr1<<"  "<<deepak8_corr2<<"  "<<(weight_deep_topcorr1/weight_deep)<<endl;
        }
        deepak8_corr1=1;deepak8_corr2=1;
        Double_t corr_deepak8_wjets[4][12]={0};
        if( channelname=="mu" ||channelname=="el" ){
            Double_t corr_deepak8_wjets_tmp[4][12]={0.654254918739,0.815210932655,0.910038032827,0.932839478126,1.05400461733,1.11882859013,1.23655146368,1.18968825811,1.19391549836,1.29924949799,1.19039910678,1.18338655006,
                0.700083721023,0.88425728413,0.874207706088,0.951375830523,1.07075137866,1.11385148772,1.13544437581,1.13784196743,1.09417735717,1.0740804446,1.24843160445,1.26607139424,
                0.912517340367,0.99157165686,1.0305264157,1.03556005263,0.984675213002,0.966808237809,0.989273515549,1.0195136825,1.01241770892,1.01241770892,1.05230576083,1.05230576083,
                1.03387166559,1.03387166559,0.865032759479,0.865032759479,0.984061741661,0.984061741661,1.00605199371,1.00605199371,1.11677164941,1.11677164941,1.11677164941,1.11677164941};
            for(int i=0;i<4;i++){for(int j=0;j<12;j++) corr_deepak8_wjets[i][j]=corr_deepak8_wjets_tmp[i][j];}
        }

        if((IsData>=2&&IsData<=8)){
            Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
            Double_t deepak8bin[13]={0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
            Double_t massbin[5]={60,100,150,200,1000000};
            int ixx=0;
            for(ixx=0;ixx<4;ixx++){
                if(jet_mass1<massbin[0]) {ixx=-1;break;}
                if(jet_mass1>=massbin[ixx]&&jet_mass1<massbin[ixx+1]) break;
            }
            int iyy=0;
            for(iyy=0;iyy<12;iyy++){
                if(jet_mass1<100&&(jetAK8puppi_dnnDecorrW<0||jetAK8puppi_dnnDecorrW>1)) {iyy=-1;break;}
                if(jet_mass1>=100&&(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))<0||((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>1)) {iyy=-1;break;}
                if(jet_mass1<100&&(jetAK8puppi_dnnDecorrW>=deepak8bin[iyy]&&jetAK8puppi_dnnDecorrW<deepak8bin[iyy+1])) break;
                if(jet_mass1>=100&&(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>=deepak8bin[iyy]&&((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))<deepak8bin[iyy+1])) break;
            }
            //cout<<ixx<<iyy<<endl;
            if(ixx>=0&&iyy>=0) deepak8_corr1=corr_deepak8_wjets[ixx][iyy];
            int ixx2=0;
            for(ixx2=0;ixx2<4;ixx2++){
                if(jet_mass2<massbin[0]) {ixx2=-1;break;}
                if(jet_mass2>=massbin[ixx2]&&jet_mass2<massbin[ixx2+1]) break;
            }
            int iyy2=0;
            for(iyy2=0;iyy2<12;iyy2++){
                if(jet_mass2<100&&(jetAK8puppi_dnnDecorrW_2<0||jetAK8puppi_dnnDecorrW_2>1)) {iyy2=-1;break;}
                if(jet_mass2>=100&&(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))<0||((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))>1)) {iyy2=-1;break;}
                if(jet_mass2<100&&(jetAK8puppi_dnnDecorrW_2>=deepak8bin[iyy2]&&jetAK8puppi_dnnDecorrW_2<deepak8bin[iyy2+1])) break;
                if(jet_mass2>=100&&(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))>=deepak8bin[iyy2]&&((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))<deepak8bin[iyy2+1])) break;
            }
            if(ixx2>=0&&iyy2>=0) deepak8_corr2=corr_deepak8_wjets[ixx2][iyy2];
            
        
        if(deepak8_corr1==0) deepak8_corr1=1;
        if(deepak8_corr2==0) deepak8_corr2=1;
        weight_deep=btagweight_center_deep*weight_tmp;
        weight_deep_topcorr1=weight_deep*deepak8_corr1*deepak8_corr2;
        weight_deep_topcorr2=weight_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;
            weight_tig_deep=btagweight_center_tig_deep*weight_tmp;
            weight_tig_deep_topcorr1=weight_tig_deep*deepak8_corr1*deepak8_corr2;
            weight_tig_deep_topcorr2=weight_tig_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;
            weight_los_deep=btagweight_center_los_deep*weight_tmp;
            weight_los_deep_topcorr1=weight_los_deep*deepak8_corr1*deepak8_corr2;
            weight_los_deep_topcorr2=weight_los_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;

        }
        if (IsData==0 ) {weight_deep=btagweight_center_deep*weight_tmp;weight_deep_topcorr1=weight_deep;weight_deep_topcorr2=weight_deep;
        weight_tig_deep=btagweight_center_tig_deep*weight_tmp;weight_tig_deep_topcorr1=weight_tig_deep;weight_tig_deep_topcorr2=weight_tig_deep;
            weight_los_deep=btagweight_center_los_deep*weight_tmp;weight_los_deep_topcorr1=weight_los_deep;weight_los_deep_topcorr2=weight_los_deep;}
        weight_deep_prefire=weight_deep*L1prefiring;
        weight_deep_prefireup=weight_deep*L1prefiringup;
        weight_deep_prefiredown=weight_deep*L1prefiringdown;
        weight_deep_lumi=weight_deep*L1prefiring*35.92;
        if (IsData==0 ){
            weight_deep_prefire=weight_deep;
            weight_deep_prefireup=weight_deep;
            weight_deep_prefiredown=weight_deep;
            weight_deep_lumi=weight_deep;
        }
}

#endif