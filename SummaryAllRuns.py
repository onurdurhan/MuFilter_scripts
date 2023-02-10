import ROOT
import os
import rootUtils as ut
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetTitleFont(102,"")
ROOT.gStyle.SetTitleFont(102,"xyz")
ROOT.gStyle.SetTextFont(10)
ROOT.gStyle.SetStatFont(102)
ROOT.gStyle.SetLabelFont(102,"xyz")
ROOT.gStyle.SetLegendFont(102)
#ROOT.gStyle.SetHistLineColor(1)
#ROOT.gStyle.SetHistLineWidth(3)


class Summary:
    def __init__(self):
        self.listOfFiles = os.listdir("/eos/user/o/onur/sndlhc_US_data/TI18/LangauFits")
        self.outf = ROOT.TFile("summary.root","RECREATE")
    def smallSiPMchannel(self,i):
        if i==2 or i==5 or i==10 or i==13: return True
        else: return False
 
    def time_evo(self):
        cv = {}
        for l in range(5):
            for bar in range(10):
                detID=int(2E4+l*1E3+bar)
                ut.bookCanvas(cv,"detID_"+str(detID))
        gr={}
        summary_graph = {}
        for l in range(5):
            for bar in range(10):
                detID=int(2E4+l*1E3+bar)
                summary_graph[detID] = ROOT.TMultiGraph()
                for c in range(16):
                    if self.smallSiPMchannel(c):continue
                    gr[str(detID)+"_"+str(c)] = ROOT.TGraph()
                    gr[str(detID)+"_"+str(c)].GetXaxis().SetTitle("Run Nr")
                    gr[str(detID)+"_"+str(c)].GetYaxis().SetTitle("MPV")
        np = 0
        for filename in self.listOfFiles:
            idx = filename.find('00')
            run = filename[idx:(idx+6)]
            f = ROOT.TFile("/eos/user/o/onur/sndlhc_US_data/TI18/LangauFits/"+filename)
            print("file is ", "/eos/user/o/onur/sndlhc_US_data/TI18/LangauFits/"+filename)
            for l in range(5):
                for bar in range(10):
                    detID=int(2E4+l*1E3+bar)
                    canvas = f.Get("langau_muon"+str(detID))
                    for c in range(16):
                        if self.smallSiPMchannel(c):continue
                        pad = canvas.GetPad(c+1)
                        h = pad.GetPrimitive("qdc_muon"+str(detID)+"_channel_"+str(c))
                        if not h :
                            print("Histogram not found", detID, c)
                            gr[str(detID)+"_"+str(c)].SetPoint(np, int(run), 0.)
                            continue
                        langau = h.GetFunction("langau")
                        if not langau  :
                            print("No Fit result for ", detID, c)
                            gr[str(detID)+"_"+str(c)].SetPoint(np, int(run), 0.)
                            continue
                        value = langau.GetParameter(1)
                        gr[str(detID)+"_"+str(c)].SetPoint(np, int(run), value)
            np+=1
        for l in range(5):
            for bar in range(10):
                detID=int(2E4+l*1E3+bar)
                for c in range(16):
                    if self.smallSiPMchannel(c):continue
                    summary_graph[detID].Add(gr[str(detID)+"_"+str(c)],"P")
                    summary_graph[detID].GetXaxis().SetTitle("Run Nr")
                    summary_graph[detID].GetYaxis().SetTitle("MPV")
                    summary_graph[detID].SetTitle(str(detID)+ "MPVs")
                    gr[str(detID)+"_"+str(c)].SetMarkerStyle(21)
        for l in range(5):
            for bar in range(10):
                detID=int(2E4+l*1E3+bar)
                cv["detID_"+str(detID)].cd()
                cv["detID_"+str(detID)].SetTitle("MPVs of " + str(detID))
                summary_graph[detID].SetTitle("MPVs of " + str(detID))
                summary_graph[detID].Draw("ALP pmc plc")

        self.outf.cd()
        for l in range(5):
            for bar in range(10):
                print("printing")
                detID= int(2E4+l*1E3+bar)
                cv["detID_"+str(detID)].Write()
        self.outf.Close()
