from scanClass import scanClass
from subprocess import call
import os
import argparse

def ParseOptionScan():

    parser = argparse.ArgumentParser(description='do all scan')
    parser.add_argument('--doLow', dest='doLow', action='store_true')
    parser.add_argument('--doOneBin', dest='doOneBin', action='store_true')
    parser.add_argument('--low', dest='lowCut', type=float, help='low cut')
    parser.add_argument('--high', dest='highCut', type=float, help='high cut')
    parser.add_argument('--cutIndex', dest='cutIndex', type=int, help='cut index')
    parser.add_argument('--processTag', dest='processTag', type=str, help='')
    parser.add_argument('--year', dest='year', type=str, help='')
    parser.add_argument('--date', dest='date', type=str, help='')
    parser.add_argument('--postfix', dest='postfix', type=str, help='')
    parser.add_argument('--specialCut', dest='specialCut', type=str, help='')

    args = parser.parse_args()
    return args

args = ParseOptionScan()

# inputs
namedict = \
{
## use for approval
"TTHHadronicTag":"/hadoop/cms/store/user/smay/ttH/FinalFitTrees/ttHHadronic__v3.10_8Oct2019_RunII_MVA_Presel_impute_addDNNs_addTopTag_FinalFitTree.root",
"TTHLeptonicTag":"/hadoop/cms/store/user/smay/ttH/FinalFitTrees/ttHLeptonic__v3.10_8Oct2019_RunII_MVA_Presel_addDNNs_FinalFitTree.root",
}

### to be configured begin ###

doLow = args.doLow 
doOneBin = args.doOneBin 
lowCut = args.lowCut 
highCut = args.highCut 
cutIndex = args.cutIndex 
processTag = args.processTag 
year = args.year 
date = args.date 
postfix = args.postfix 
specialCut = args.specialCut # pT splitting

#doLow = False
#doOneBin = False
#lowCut = 0
#highCut = 1
#cutIndex = 50
#processTag = "TTHLeptonicTag"
#year = "2020"
#date = "20200315"
#postfix = "pt_lep_binning_0_120_step1"
#specialCut = "(global_features[28]*mass < 120)" # pT splitting
### to be configured end ###

### hardcode begin ###
modelpathBase = os.getcwd() + "/models/"
plotpathBase = os.getcwd() + "/plots/"
#modelpathBase = "/home/users/hmei/ttH/makebdtbinstest/doscan/models/"
#plotpathBase = "/home/users/hmei/public_html/"
combineCommand = "combine -M Significance CMS-HGG_mva_13TeV_datacard.txt -t -1 --expectSignal=1"
nSplit = 100

selectSignal = "process_id == 0"
globalCuts = "(global_features[7] > 0.33 && global_features[8] > 0.25) && " + specialCut
higgsCuts = "(signal_mass_category == 127)"
bkgCuts = "(sample_id == 0 && process_id <= 10)"
massCuts = "(mass > 100 && mass < 180)"
massCuts_blind = "((mass > 100 && mass < 120) || ( mass > 130 && mass < 180))" 
mvaCuts = "(mva_score < " + str(highCut) + " && mva_score > " + str(lowCut) + ")"

procDict = \
{
        "sig": {"name":"ttH_hgg", "procCut":"process_id == 0"},
        "resBkg": {"name":"ggH_hgg", "procCut":"process_id == 14", "replaceNorm":True},
        "nonResBkg": {}
}

sigs = ["ttH_hgg"]
bkgs = ["ggH_hgg", "bkg_mass"]

# setup scan
scanConfig= {\
        "processTag":processTag,
        "baseSelection_higgs": globalCuts + " && " + higgsCuts + " && " + massCuts,
        "baseSelection_bkg": globalCuts + " && " + bkgCuts + " && " + massCuts_blind,
        "mvaSelection": mvaCuts,
        "signame":"tth_hgg",
        "filename":namedict[processTag],
        "modelpath":modelpathBase +postfix+"_"+date+"/",
        "plotpath":plotpathBase +year+"/"+date+"_"+postfix+"/",
        "var":"mass",
        "weightVar":"weight",
        "mvaName":"mva_score"
        }
### hardcode end ###

testScan = scanClass(scanConfig)
testScan.cleanDir()
testScan.getTree()
    
testScan.quantiles_to_mva_score(nSplit, selectSignal)
testScan.mvaBoundaries.insert(0,[highCut,])
testScan.mvaBoundaries.append([lowCut,])
print (len(testScan.mvaBoundaries), testScan.mvaBoundaries)

lowCuts = [testScan.mvaBoundaries[cutIndex][0], lowCut]
highCuts = [highCut, testScan.mvaBoundaries[cutIndex][0]] 

if doLow:
    # in this case, only consider cuts from range: [0, x] -> [y,x] + [x,1] 
    lowCuts = [highCut, testScan.mvaBoundaries[cutIndex][0]]
    highCuts = [1, highCut]

if doOneBin:
    if not doLow:
        lowCuts = [lowCut]
        highCuts = [highCut]
    else:
        lowCuts = [highCut]
        highCuts = [1]

call("echo mvaScore " + str(testScan.mvaBoundaries[cutIndex][0]) + " > " + scanConfig["modelpath"] + "/mvaScore_"+str(cutIndex)+".txt" , shell=True)
testScan.MakeModelsForTwoBins(procDict, "sig", lowCuts, highCuts, cutIndex)
testScan.MakeModelsForTwoBins(procDict, "resBkg", lowCuts, highCuts, cutIndex)
testScan.MakeModelsForTwoBins(procDict, "nonResBkg", lowCuts, highCuts, cutIndex)

# make datacard
if not doOneBin:
    testScan.MakeCard2Bin(sigs, bkgs, cutIndex)
else:
    testScan.MakeCard1Bin(sigs, bkgs, cutIndex)

testScan.runSignificance(combineCommand, cutIndex)

# parse result
testScan.PareseResult(cutIndex)

# set up plotdir
call("chmod -R 755 " + testScan.plotpath, shell=True)
