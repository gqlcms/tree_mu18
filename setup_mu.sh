#outfile=/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/xulyu/WWW/tree_toptow_btiglos_matchedcorr
#outfile=/eos/cms/store/user/xulyu/data2016_6.27_dnn_deepb_mcw/tree_toptow_btiglos_matchedcorr
outfile=/eos/cms/store/user/xulyu/data2016_8_6_parton_match/tree
outfile=./
hadd -f ${outfile}/mu_PKUTree_SingleTop_xww.root mu_out_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root        mu_out_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root mu_out_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root  mu_out_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root mu_out_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root

cp mu_out_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root ${outfile}/mu_PKUTree_TTBARpowheg_xww.root 
hadd -f ${outfile}/mu_PKUTree_VV_xww.root mu_out_WWToLNuQQ_13TeV-powheg.root mu_out_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root  mu_out_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root
hadd -f  ${outfile}/mu_PKUTree_WJetsPt180_xww.root  mu_out_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

hadd -f ${outfile}/mu_PKUTree_Muon.root mu_out_single*.root
hadd -f ${outfile}/mu_PKUTree_MET.root mu_out_MET*.root
hadd -f ${outfile}/mu_PKUTree_2017vvv.root ${outfile}/mu_PKUTree_MET.root ${outfile}/mu_PKUTree_Muon.root
channel=mu_
hadd -f ${outfile}/${channel}PKUTree_VVV.root ${channel}out_WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root ${channel}out_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root ${channel}out_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root ${channel}out_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root
hadd -f ${outfile}/${channel}PKUTree_ttV.root ${channel}out_ttWJets_13TeV_madgraphMLM.root  ${channel}out_ttZJets_13TeV_madgraphMLM-pythia8.root
hadd -f ${outfile}/${channel}PKUTree_QCD.root ${channel}out_QCD_HT*
hadd -f ${outfile}/${channel}PKUTree_QCDPt.root ${channel}out_QCD_Pt*
rm mu_out*txt
#mv mu_out* ${outfile}
'''
cp ${outfile}/mu_out_M3000_R0-5.root ${outfile}/mu_out_case1_off_3000.root
cp ${outfile}/mu_out_M3000_R0-06.root ${outfile}/mu_out_case2_off_3000.root
cp ${outfile}/mu_out_M1500_R0-5.root ${outfile}/mu_out_case1_off_1500.root
cp ${outfile}/mu_out_M1500_R0-12.root ${outfile}/mu_out_case2_off_1500.root
