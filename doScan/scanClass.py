import ROOT
#import filenameDict as filenameDict
#import processIDMap as processIDMap
from cardMaker import makeCards
import numpy
import root_numpy
import subprocess
from tdrStyle import *
setTDRStyle()


ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
#gSystem.Load("$CMSSW_BASE/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so")
#ROOT.gSystem.Load("$CMSSW_BASE/lib/slc6_amd64_gcc630/libHiggsAnalysisCombinedLimit.so")
ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
ROOT.gSystem.AddIncludePath("-Iinclude/")
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.DataHandling)
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING) 

class scanClass():

    def __init__(self, config):

        self.debug = False

        self.processTag = config["processTag"]
        self.filename = config["filename"]
        self.baseSelection_higgs= config["baseSelection_higgs"]
        self.baseSelection_bkg= config["baseSelection_bkg"]
        self.mvaSelection= config["mvaSelection"]
        self.modelpath = config["modelpath"]
        self.plotpath = config["plotpath"]
        self.var = config["var"]
        self.weightVar = config["weightVar"]
        self.mvaName = config["mvaName"]
        #self.tag = config["tag"]

        self.treename = "t"
        self.config = config

        # this is an initialization, will be overwritten when needed
        self.selection = self.baseSelection_higgs

    def getTree(self):

        #self.filename = filenameDict.namedict[self.tag]
        self.file = ROOT.TFile.Open(self.filename)
        self.tree = self.file.Get(self.treename)


    def getEfficiency(self, cut_denominator, cut_numerator):

        h_denominator = ROOT.TH1F("h_denominator", "", 320, 100, 180)
        self.tree.Project(h_denominator.GetName(), self.var, self.weightVar + "*(" + cut_denominator + ")")
        h_numerator = ROOT.TH1F("h_numerator", "", 320, 100, 180)
        self.tree.Project(h_numerator.GetName(), self.var, self.weightVar + "*(" + cut_numerator + ")")
        n_denominator = h_denominator.Integral()
        n_numerator = h_numerator.Integral()

        efficiency = n_numerator/n_denominator
        print(efficiency)

    def cleanDir(self):

        pathCmd =  "mkdir -p " + self.modelpath + ";"
        pathCmd += "rm " + self.modelpath+ "*;"
        pathCmd += "mkdir -p " + self.plotpath + ";"
        pathCmd += "rm " + self.plotpath+ "*;"
        pathCmd += "cp ~/public_html/backup/index.php " + self.plotpath

        subprocess.call(pathCmd, shell=True)

    def quantiles_to_mva_score(self, n_quantiles, extraCut):
        # for given selection, return mva_scores corresponding to each quantile in n_quantiles
        self.mvaBoundaries = []
        # get a numpy array from tree
        #print ("cut: ", self.selection)
        mva_scores = (root_numpy.tree2array(self.tree, branches = [self.mvaName], selection = self.baseSelection_higgs + "&&" + self.mvaSelection + "&&" + extraCut))

        sorted_mva = numpy.flip(numpy.sort(mva_scores), 0)
        quantiles = []
        for i in range(n_quantiles):
            idx = int((float(i+1) / float(n_quantiles)) * len(sorted_mva)) - 1
            quantiles.append(float(i+1) / float(n_quantiles))
            self.mvaBoundaries.append(sorted_mva[idx])
        #print mva
        #return mva

    def set_tag_savename_sig(self, process, binIndex, cutIndex):
        '''
        binIndex = 0 or 1, each time make one cut, ends up with two bins
        cutIndex = 0 to N, used to locate the cut value, correspond to one BDT score
        process = ttH_hgg, ggH_hgg ...
        '''
        self.tag = "hggpdfsmrel_" + process + "_" + self.processTag + "_" + str(binIndex) + "_" + str(cutIndex) 
        self.savename = "CMS-HGG_sigfit_mva_" + process + "_" + self.processTag + "_" + str(binIndex) + "_" + str(cutIndex) 

    def set_tag_savename_bkg(self, binIndex, cutIndex):
        '''
        binIndex = 0 or 1, each time make one cut, ends up with two bins
        cutIndex = 0 to N, used to locate the cut value, correspond to one BDT score
        '''
        self.tag = "CMS_hgg_bkgshape_" + self.processTag + "_" + str(binIndex) + "_" + str(cutIndex) 
        self.savename = "CMS-HGG_bkg_" + self.processTag + "_" + str(binIndex) + "_" + str(cutIndex) 

    def setSelection_sig(self, processSelect, mvaCutLow, mvaCutHigh):
        '''
        processSelect is a cut that select desired process, for example process_id == 0
        '''
        self.selection = self.baseSelection_higgs + " && (" + processSelect + " && " + self.mvaName + " > " + str(mvaCutLow) + " && " + self.mvaName + " < " + str(mvaCutHigh) + ")"

    def setSelection_bkg(self, mvaCutLow, mvaCutHigh):
        self.selection = self.baseSelection_bkg + " && (" + self.mvaName + " > " + str(mvaCutLow) + " && " + self.mvaName + " < " + str(mvaCutHigh) + ")"


    def makeSignalModel(self, workspaceName, config):

        rooVar = "CMS_hgg_mass"

        replaceNorm = config["replaceNorm"]
        norm_in = config["norm_in"]
        fixParameters = config["fixParameters"]
        nofit = config["nofit"]

        w = ROOT.RooWorkspace(workspaceName)
        w.factory(rooVar + "[100,180]")
        w.factory("MH[125]")

        h_mgg = ROOT.TH1F("h_mgg", "h_mgg", 320, 100, 180)
        h_mgg.Sumw2()
        self.tree.Project(h_mgg.GetName(), self.var, self.weightVar + "*(" + self.selection + ")")
        d_mgg = ROOT.RooDataHist("roohist_data_mass_" + self.tag, "", ROOT.RooArgList(w.var(rooVar)), h_mgg, 1)
        print ("bin dataset", h_mgg.Integral(), d_mgg.sumEntries(), d_mgg.numEntries() )
        #return 1

        # normalization
        norm = d_mgg.sumEntries()
        if replaceNorm:
            norm = norm_in
        if norm <= 0:
            norm = 1e-09
        #if ("FCNC" in self.tag) and ("Leptonic" in self.tag):
        #    norm = norm*1/1.527
        if nofit:
            return norm

        rv_norm = ROOT.RooRealVar(self.tag+"_norm", "", norm)

        self.norm = norm
        # pdf
        w.factory("DoubleCB:"+self.tag+"(" + rooVar + ", mean_"+self.tag+"[125,120,130], sigma_"+self.tag+"[1,0,5], a1_"+self.tag+"[1,0,10], n1_"+self.tag+"[1,0,10], a2_"+self.tag+"[1,0,10], n2_"+self.tag+"[1,0,10])")
        exPdf = ROOT.RooExtendPdf("extend" + self.tag, "", w.pdf(self.tag), rv_norm)

        # fit
        w.pdf(self.tag).fitTo(d_mgg, ROOT.RooFit.PrintLevel(-1))

        getattr(w,'import')(rv_norm)
        getattr(w,'import')(exPdf)

        # frame
        frame = w.var("CMS_hgg_mass").frame()
        d_mgg.plotOn(frame)
        w.pdf(self.tag).plotOn(frame)

        # plot
        c1 = ROOT.TCanvas("c1", "c1", 800, 800)
        dummy = ROOT.TH1D("dummy","dummy",1,100,180)
        dummy.SetMinimum(0)
        dummy.SetMaximum(h_mgg.GetMaximum()*1.2)
        dummy.SetLineColor(0)
        dummy.SetMarkerColor(0)
        dummy.SetLineWidth(0)
        dummy.SetMarkerSize(0)
        dummy.GetYaxis().SetTitle("Events")
        dummy.GetYaxis().SetTitleOffset(1.3)
        dummy.GetXaxis().SetTitle("m_{#gamma#gamma} (GeV)")
        dummy.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.6*c1.GetTopMargin())
        latex.SetTextFont(42)
        latex.SetTextAlign(11)
        latex.SetTextColor(1)

        latex.DrawLatex(0.5, 0.85, (self.tag.split("_"))[-2] + "_" + (self.tag.split("_"))[-1])
        latex.DrawLatex(0.5, 0.78, "mean = " + str(round(w.var("mean_"+self.tag).getVal(), 3) ) + " #pm " + str(round(w.var("mean_"+self.tag).getError(), 3) ))
        latex.DrawLatex(0.5, 0.71, "sigma = " + str(round(w.var("sigma_"+self.tag).getVal(), 3) ) + " #pm " + str(round(w.var("sigma_"+self.tag).getError(), 3) ))
        latex.DrawLatex(0.5, 0.64, "a1 = " + str(round(w.var("a1_"+self.tag).getVal(), 3) ) + " #pm " + str(round(w.var("a1_"+self.tag).getError(), 3) ))
        latex.DrawLatex(0.5, 0.57, "a2 = " + str(round(w.var("a2_"+self.tag).getVal(), 3) ) + " #pm " + str(round(w.var("a2_"+self.tag).getError(), 3) ))
        latex.DrawLatex(0.5, 0.50, "n1 = " + str(round(w.var("n1_"+self.tag).getVal(), 3) ) + " #pm " + str(round(w.var("n1_"+self.tag).getError(), 3) ))
        latex.DrawLatex(0.5, 0.43, "n2 = " + str(round(w.var("n2_"+self.tag).getVal(), 3) ) + " #pm " + str(round(w.var("n2_"+self.tag).getError(), 3) ))

        latex.DrawLatex(0.5, 0.33, "norm = " + str(round(self.norm, 3)))
        frame.Draw("same")
        c1.SaveAs(self.plotpath + "/fit_sig_" + self.savename + ".png")
        c1.SaveAs(self.plotpath + "/fit_sig_" + self.savename + ".pdf")

        if fixParameters:
           w.var("mean_"+self.tag).setConstant()
           w.var("sigma_"+self.tag).setConstant()
           w.var("a1_"+self.tag).setConstant()
           w.var("a2_"+self.tag).setConstant()
           w.var("n1_"+self.tag).setConstant()
           w.var("n2_"+self.tag).setConstant()

        w.writeToFile(self.modelpath + "/" + self.savename + ".root")

    def makeBackgroundModel(self, workspaceName, datasetTag):

        rooVar = "CMS_hgg_mass"

        w = ROOT.RooWorkspace(workspaceName)
        w.factory(rooVar + "[100,180]")
        w.factory("MH[125]")

        h_mgg = ROOT.TH1F("h_mgg", "h_mgg", 320, 100, 180)
        h_mgg.Sumw2()
        self.tree.Project(h_mgg.GetName(), self.var, self.weightVar + "*(" + self.selection + ")")
        d_mgg = ROOT.RooDataHist("roohist_data_mass_" + datasetTag, "", ROOT.RooArgList(w.var(rooVar)), h_mgg, 1)
        #print "bin dataset", h_gg.Integral(), d_mgg_bin.sumEntries(), d_mgg_bin.numEntries()

        # normalization
        norm = d_mgg.sumEntries()

        # set variable range
        w.var(rooVar).setRange("SL", 100, 120)
        w.var(rooVar).setRange("SU", 130, 180)
        w.var(rooVar).setRange("full", 100, 180)
        w.var(rooVar).setRange("blind",120,130)

        # pdf
        w.factory("Exponential:"+self.tag+"(" + rooVar + ", tau[-2,-10,0])")
        w.factory("ExtendPdf:"+self.tag+"_ext("+self.tag+", nevt[100,0,10000000], 'full')")

        # fit
        w.pdf(self.tag+"_ext").fitTo(d_mgg, ROOT.RooFit.Range("SL,SU"), ROOT.RooFit.Extended(True), ROOT.RooFit.PrintLevel(-1))

        frame = w.var(rooVar).frame()
        d_mgg.plotOn(frame, ROOT.RooFit.Binning(80))
        w.pdf(self.tag+"_ext").plotOn(frame)

        c1 = ROOT.TCanvas("c1", "c1", 800, 800)
        dummy = ROOT.TH1D("dummy","dummy",1,100,180)
        dummy.SetMinimum(0)
        dummy.SetMaximum(h_mgg.GetMaximum()*1.2*4)
        dummy.SetLineColor(0)
        dummy.SetMarkerColor(0)
        dummy.SetLineWidth(0)
        dummy.SetMarkerSize(0)
        dummy.GetYaxis().SetTitle("Events")
        dummy.GetYaxis().SetTitleOffset(1.3)
        dummy.GetXaxis().SetTitle("m_{#gamma#gamma} (GeV)")
        dummy.Draw()

        frame.Draw("same")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.6*c1.GetTopMargin())
        latex.SetTextFont(42)
        latex.SetTextAlign(11)
        latex.SetTextColor(1)

        latex.DrawLatex(0.4, 0.85, (self.tag.split("_"))[-2] + "_" + (self.tag.split("_"))[-1])
        latex.DrawLatex(0.4, 0.78, "nEvents = " + str(round(w.var("nevt").getVal(), 3) ) + " #pm " + str(round(w.var("nevt").getError(), 3) ))
        latex.DrawLatex(0.4, 0.71, "#tau = " + str(round(w.var("tau").getVal(), 3) ) + " #pm " + str(round(w.var("tau").getError(), 3) ))

        c1.SaveAs(self.plotpath + "/fit_bkg_" + self.savename + ".png")
        c1.SaveAs(self.plotpath + "/fit_bkg_" + self.savename + ".pdf")

        #l=ROOT.RooArgSet(w.var("CMS_hgg_mass"))
        #frac = w.pdf(tag).createIntegral(l,l,"blind")
        #print "frac", frac.getVal()
        #norm = norm/(1-frac.getVal())
        nEvt = w.var("nevt").getVal()
        w.factory(self.tag+"_norm["+str(nEvt)+",0,"+str(3*nEvt)+"]")
        #print nEvt

        getattr(w,'import')(d_mgg, ROOT.RooCmdArg())
        w.writeToFile(self.modelpath + "/" + self.savename + ".root")

        from subprocess import call
        call("echo " + str(nEvt) + " > " + self.modelpath + "/" + self.savename + "_nbkg.txt" , shell=True)

    def MakeModelsForTwoBins(self, procDict, proc, lowCuts, highCuts, cutIndex):
        '''
        binIndex = 0 or 1, each time make one cut, ends up with two bins
        cutIndex = 0 to N, used to locate the cut value, correspond to one BDT score
        process = ttH_hgg, ggH_hgg ...
        '''
        thisDict = procDict[proc]
        sigDict = procDict["sig"]

        for binIndex in range(len(lowCuts)):
            # if only wants to one bin, then pass only one cut in lowCuts/highCuts

            if proc == "nonResBkg":
                self.set_tag_savename_bkg(binIndex, cutIndex)
                self.setSelection_bkg(lowCuts[binIndex], highCuts[binIndex])
                self.makeBackgroundModel("wbkg_13TeV", self.processTag + "_" + str(binIndex) + "_" + str(cutIndex))                

            elif "replaceNorm" not in thisDict:
                # set tag
                self.set_tag_savename_sig(thisDict["name"], binIndex, cutIndex)
                # set selection
                self.setSelection_sig(thisDict["procCut"], lowCuts[binIndex], highCuts[binIndex])
                # extra config
                sigConfig = {"replaceNorm":False, "norm_in":-1, "fixParameters":True, "nofit": False}
                self.makeSignalModel("wsig_13TeV", sigConfig)

            else:
                # set tag
                self.set_tag_savename_sig(thisDict["name"], binIndex, cutIndex)
                # set selection
                self.setSelection_sig(thisDict["procCut"], lowCuts[binIndex], highCuts[binIndex])
                # extra config
                sigConfig = {"replaceNorm":False, "norm_in":-1, "fixParameters":True, "nofit": True}
                thisNorm = self.makeSignalModel("wsig_13TeV", sigConfig)
                # set selection
                self.setSelection_sig(sigDict["procCut"], lowCuts[binIndex], highCuts[binIndex])
                # extra config
                sigConfig = {"replaceNorm":True, "norm_in":thisNorm, "fixParameters":True, "nofit": False}
                self.makeSignalModel("wsig_13TeV", sigConfig)

    def MakeCard2Bin(self, sigList, bkgList, cutIndex):
        # make datacard
        #sigList = ["ttH_hgg"]
        #bkgList = ["ggH_hgg", "bkg_mass"]
        cardName = "CMS-HGG_mva_13TeV_datacard_"+str(cutIndex)+".txt"
        card = makeCards(self.modelpath, cardName)
        card.WriteCard(sigList, bkgList, [self.processTag + "_0", self.processTag + "_1"], "_" + str(cutIndex))

    def MakeCard1Bin(self, sigList, bkgList, cutIndex):
        # make datacard
        #sigList = ["ttH_hgg"]
        #bkgList = ["ggH_hgg", "bkg_mass"]
        cardName = "CMS-HGG_mva_13TeV_datacard_"+str(cutIndex)+".txt"
        card = makeCards(self.modelpath, cardName)
        #card.WriteCard(sigList, bkgList, [self.processTag + "_0", self.processTag + "_1"], "_" + str(cutIndex))
        card.WriteCard(sigList, bkgList, [self.processTag + "_0"], "_" + str(cutIndex))

    def runSignificance(self, combineCommand, cutIndex):

        cardName = "CMS-HGG_mva_13TeV_datacard_"+str(cutIndex)+".txt"

        combineCommand = "cd " + self.modelpath + ";" + combineCommand
        combineCommand = combineCommand.replace("CMS-HGG_mva_13TeV_datacard.txt", cardName)
        combineCommand += " > combine_" + str(cutIndex) + ".txt ;"

        print (combineCommand)
        cmd = "combineCommand" + str(cutIndex) + ".sh"
        print (self.modelpath + cmd)
        with open(self.modelpath + cmd, "w") as fh:
            fh.write(combineCommand)

        subprocess.call("source " + self.modelpath + cmd, shell=True)


    def PareseResult(self, cutIndex):

        script = """
        cd modelpath
        sig=`grep "Significance:" combine_cutIndex.txt`
        nbkg0=`sed -n '1p' CMS-HGG_bkg_processTag_0_cutIndex_nbkg.txt`
        nbkg1=`sed -n '1p' CMS-HGG_bkg_processTag_1_cutIndex_nbkg.txt`
        mva=`grep mva mvaScore_cutIndex.txt`
        echo $sig" "${nbkg0}" "${nbkg1}" "${mva}
        """

        script = script.replace('processTag', self.processTag).replace('modelpath', self.modelpath).replace('cutIndex', str(cutIndex))

        scriptName ="parse_" + str(cutIndex) + ".sh" 
        with open(scriptName, "w") as fh:
            fh.write(script)
        _ = subprocess.check_output("chmod u+x " + scriptName, shell=True, executable="/bin/bash")

        cmd = "source " + scriptName + " > parseout_" + str(cutIndex) + ".txt"
        #return subprocess.check_output(cmd, shell=True, executable="/bin/bash")
        output =  subprocess.check_output(cmd, shell=True, executable="/bin/bash")
