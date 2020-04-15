import ROOT
#import filenameDict as filenameDict
#import processIDMap as processIDMap
import root_numpy

from tdrStyle import *
setTDRStyle()

ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
ROOT.gSystem.Load("$CMSSW_BASE/lib/slc6_amd64_gcc630/libHiggsAnalysisCombinedLimit.so")
ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
ROOT.gSystem.AddIncludePath("-Iinclude/")
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.DataHandling)
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING) 

class makeModel():

    def __init__(self, config):

        self.debug = False

        self.tag = config["tag"]
        self.selection = config["selection"]
        self.plotpath = config["plotpath"]
        self.modelpath = config["modelpath"]
        self.savename = config["savename"]
        self.var = config["var"] #mass
        self.weightVar = config["weightVar"] #weight
        self.norm = -1

        self.treename = "t"


    def getTree(self, t):

        self.tree = t


    #def getTreeFromFile(self):

    #    self.filename = filenameDict.namedict[self.tag]
    #    self.file = ROOT.TFile.Open(self.filename)
    #    self.tree = self.file.Get(self.treename)


    def cleanDir(self):

        pathCmd =  "mkdir -p " + self.modelpath + ";"
        pathCmd += "rm " + self.modelpath+ "*;"
        pathCmd += "mkdir -p " + self.plotpath + ";"
        pathCmd += "rm " + self.plotpath+ "*;"
        pathCmd += "cp ~/public_html/tmpFile/index.php " + self.plotpath

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
        if norm <= 0:
            norm = 1e-09
        if replaceNorm:
            norm = norm_in

        if nofit:
            return norm
        #if ("FCNC" in self.tag) and ("Leptonic" in self.tag):
        #    norm = norm*1/1.527
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
