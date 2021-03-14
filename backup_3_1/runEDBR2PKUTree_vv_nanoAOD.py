import os
import random
import sys
import ctypes
import time

import ROOT

from optparse import OptionParser
parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
(options, args) = parser.parse_args()


ROOT.gSystem.Load("EDBR2PKUTree_C.so")
from ROOT import EDBR2PKUTree, TFile, TTree



if __name__ == '__main__': 
    input_dir="/eos/user/q/qiguo/NanoAOD/Ntuple/21_3_1/merge/"
    treen="tree_sf_new/"
#    input_dir="/afs/cern.ch/user/l/lewang/eos/cms/store/user/lewang/toapproval_QCD_Wprime/"
    infile=open("filelist_vv_nanoAOD.txt")

    totalsamples=infile.readlines()

    for sample in totalsamples:
        if sample[0]!="#":

            print "Loading ", sample
            sampleconfig=sample.split(';')
            samplename=sampleconfig[0]
            sampleXS=sampleconfig[1]
            sampleNum=sampleconfig[2]

            samplename=samplename.replace(' ', '').replace('\n','').replace('/n','') 
            sampleXS  =float(sampleXS.replace(' ', '').replace('\n','').replace('/n','') )
            sampleNum =int(sampleNum.replace(' ', '').replace('\n','').replace('/n','') )

            if (sampleXS-1.)==0: IsData=0
            elif "WJetsToLNu_HT-100To200" in samplename: IsData=2
            elif "WJetsToLNu_HT-200To400" in samplename: IsData=3
            elif "WJetsToLNu_HT-400To600" in samplename: IsData=4
            elif "WJetsToLNu_HT-600To800" in samplename: IsData=5
            elif "WJetsToLNu_HT-800To1200" in samplename: IsData=6
            elif "WJetsToLNu_HT-1200To2500" in samplename: IsData=7
            elif "WJetsToLNu_HT-2500ToInf" in samplename: IsData=8 
            elif "ST_" in samplename: IsData=9 
            elif "TT_" in samplename: IsData=9 
            elif "ST_tW_" in samplename: IsData=11 
            elif "_R0" in samplename: IsData=10 
            else: IsData=1

            if "MET" in samplename: IsMET=1
            else: IsMET=0

            print samplename
            print sampleXS
            print sampleNum 
            print IsData
            print IsMET

            filein =TFile(input_dir+samplename);
            treew = filein.Get("t")
            #print "Open Tree Done"

            ##channel="el"
            ##if "singleMuon" in samplename: 
            ##    print "ignore mu data in el channel"
            ##else: 
            ##    analyzer_el=EDBR2PKUTree(treein,channel+"_out_"+samplename); 
            ##    analyzer_el.Loop(channel, sampleXS, treew, IsData); 
            ##    analyzer_el.endJob();
            numofevent = filein.Get("Wgt__h_nevents").GetBinContent(1)
            numofevent = int(numofevent)
            print "numberofevent:",numofevent

            channel="mu"
            if "singleEl" in samplename: 
                print "ignore el data in mu channel"
            else:
                analyzer_mu=EDBR2PKUTree(treew,input_dir+treen+channel+"_out_"+samplename,IsData,numofevent);
                analyzer_mu.Loop(channel, sampleXS, treew, IsData,IsMET);
                analyzer_mu.endJob();
