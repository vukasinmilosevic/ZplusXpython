#import stage
import numpy as np
import ROOT
from ROOT import TCanvas, TH1F, TF1, TLegend, gPad, THStack, TColor
#import tdrstyle
import os.path
import os
from array import array
import yaml
import sys
import math

from analyzeZX import *


#fileList = ["/eos/user/y/yujil/SkimNanoAOD_2022_ZXCRv3/SkimNanoAOD_2022_ZXCRv3_noDuplicates.root",
#            "/eos/user/y/yujil/HZZRun3Share/ZX_new/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8.root",
#            "/eos/user/y/yujil/HZZRun3Share/ZX_new/TT_TuneCP5_13p6TeV_powheg-pythia8.root",
#            "/eos/user/y/yujil/HZZRun3Share/ZX_new/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8.root",
#            "/eos/user/y/yujil/HZZRun3Share/ZX_new/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8.root"
#            ]
fileList = ["/eos/user/y/yujil/HZZRun3Share/Data2022CD_noDuplicate.root",
            "/eos/user/y/yujil/HZZRun3Share/WZto3LNu-1Jets-4FS_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.root",
            "/eos/user/y/yujil/HZZRun3Share/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8.root",
            "/eos/user/y/yujil/HZZRun3Share/ZXCR/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8.root",
            "/eos/user/y/yujil/HZZRun3Share/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8.root",
           ]


#RootNickNames = ["Data","DY50","TT","WZ","ZZ"]
RootNickNames = ["Data","WZ","DY50","ZZ","TT"]
# For testing
#fileList = ["../DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_2018filter2l_new_ZX.root",]
#fRootNickNames = ["DY50"]

print("First stage of processing (FR computation and CR histogram creation) has been initiated.\n")

#for i in range(len(fileList)):
for i in range(0,5):
    inFile =  ROOT.TFile.Open(fileList[i], "READ")
    if i == 0:
        tree = inFile.Get("Events")
    else:
        tree = inFile.Get("Events")


    print ("- File: \"" + fileList[i] + "\" has been opened.\n-- It has "+str(tree.GetEntries())+" events.\n")
    analyzeZX(tree, RootNickNames[i], "mEt")
    inFile.Close()

    


