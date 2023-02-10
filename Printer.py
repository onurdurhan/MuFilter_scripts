import ROOT
import os
import rootUtils as ut
ROOT.gStyle.SetTitleFont(102,"")
ROOT.gStyle.SetTitleFont(102,"xyz")
ROOT.gStyle.SetStatFont(102)
ROOT.gStyle.SetLabelFont(102,"xyz")
ROOT.gStyle.SetLegendFont(102)
ROOT.gStyle.SetHistLineColor(1)
ROOT.gStyle.SetHistLineWidth(3)
ROOT.gROOT.SetBatch(ROOT.kTRUE)




class Printer:
    def __init__(self,fname):
        self.canvases = {}
        self.f = ROOT.TFile(fname)
        if fname.find('TI18')>0 : data = 'TI18'
        else : data = 'H8'
        idx = fname.find('00')
        self.runNr = fname[idx:(idx+6)]
 
   
    def smallSiPMchannel(self,i):
        if i==2 or i==5 or i==10 or i==13: return True
        else: return False
 
    def print_langau(self):
        mp_langau_gr = ROOT.TGraphErrors()
        chi2_gr =  ROOT.TGraph()
        mp_langau_gr.SetMarkerStyle(21)
        chi2_gr.SetMarkerStyle(21)
        fc = 0
        np = 0
        os.system("mkdir "+self.runNr)
        import numpy
        domain = numpy.linspace(0, 100.,251)
        cuts = {}
        Fit_Statistics_txt = open(self.runNr+"/Fit_Statistics.txt","w")
        Fit_Statistics_29 = open(self.runNr+"/Fit_Statistics_29_version2.tex", "w")
        Fit_Statistics_29.write("\\documentclass{article} \n")
        Fit_Statistics_29.write("\\usepackage{graphicx} \n")
        Fit_Statistics_29.write("\\begin{document} \n")
        Fit_Statistics_txt.write("Channel ID        MPV         Width (Scale)            Norm              Sigma             N2             chi2/Ndof        Threshold \n")
        for l in range(5):
            Fit_Statistics_29.write("\\clearpage \n")
            for bar in range(10):
                Fit_Statistics_29.write("\\begin{table}[ht] \n")
                Fit_Statistics_29.write("\\centering \n")
                Fit_Statistics_29.write("\\resizebox{\\textwidth}{!}{\\begin{tabular} { |l|| r | r | r | r | r  | r | r |} \n")
                Fit_Statistics_29.write("\\multicolumn{8}{c}{Fit Statistics : US Plane-%s Bar-%s} \\\\ \n" % (str(l),str(bar) ))
                Fit_Statistics_29.write("\\hline \n")
                Fit_Statistics_29.write("Channel ID   & $ MPV $ & $ Width (Scale) $ & $Norm$ & $\sigma$ & $N2$ & $\chi^2/ndof $ & $Threshold$\\\\ \n")
#Fit_Statistics_29.write("Channel ID   & $Width (Scale)$ & $MPV$ & $Norm$ & $\sigma$ & $N2$ & $\chi^2/ndof $ & $Threshold$\\\\ \n")                
                Fit_Statistics_29.write( "\\hline \n" )
                detID=int(2E4+l*1E3+bar)
                ut.bookCanvas(self.canvases,"threshold_"+str(detID),"Langau Distributions with Threshold",700,500,4,4)
                canvas = self.f.Get("langau_muon"+str(detID))
                for c in range(16):
                    fc += 1
                    if self.smallSiPMchannel(c):continue
                    pad = canvas.GetPad(c+1)
                    h = pad.GetPrimitive("qdc_muon"+str(detID)+"_channel_"+str(c))
                    if not h :
                        print("Histogram not fount", detID, c)
                        continue
                    langau = h.GetFunction("langau")
                    if not langau  :
                        print("No Fit result for ", detID, c)
                        continue
                    for x in range(h.GetNbinsX()):
                            threshold = h.GetBinCenter(x)
                            ratio = h.Integral(x,100)/h.Integral(0,100)
                            if ratio <= 0.9999 :
                                break
                    print("The ratio under integrals : ", ratio , detID, c , h.GetMaximum())
                    if Fit_Statistics_29!= None:
                        Fit_Statistics_29.write(format("%s" % "DetID "+str(detID)+" Channel "+str(c)))
                        Fit_Statistics_txt.write(format("%s" % str(detID)+" "+str(c)+"     " ))
                        params = [1,0,2,3,4]
#                        for par in range(langau.GetNpar()):
                        for par in params:
                            value = langau.GetParameter(par)
                            value_error = langau.GetParError(par)
                            Fit_Statistics_29.write( " & $ %.2f  \\pm  %.2f $" % (value, value_error))
                            Fit_Statistics_txt.write("%.2f +-  %.2f     " % (value, value_error))
                        chi2 = langau.GetChisquare()
                        ndof = langau.GetNDF()+1E-12
                        cndf = chi2/ndof
                        if cndf > 200 : cndf = 999
                        Fit_Statistics_29.write("& $%.2f$" % cndf)
                        Fit_Statistics_29.write("& $%.2f$" % threshold)
                        Fit_Statistics_29.write( " \\\\") ### end of row
                        Fit_Statistics_29.write( "\n")
                        Fit_Statistics_txt.write("  %.2f     " % cndf)
                        Fit_Statistics_txt.write("%.2f \n" % threshold)
                    self.canvases["threshold_"+str(detID)].cd(c+1)
                    h.GetXaxis().SetTitle("QDC")
                    ROOT.gPad.Update()
                    rc = h.Draw()
                    ROOT.gStyle.SetOptFit(1111)
                    self.canvases["threshold_"+str(detID)].Update()
                    cuts['threshold_'+str(detID)+"_"+str(c)] = ROOT.TLine(threshold, 0., threshold, h.GetMaximum()/2)    
                    cuts['threshold_'+str(detID)+"_"+str(c)].SetLineColor(ROOT.kBlue)
                    cuts['threshold_'+str(detID)+"_"+str(c)].SetLineWidth(2)
                    cuts['threshold_'+str(detID)+"_"+str(c)].Draw("same")
                    ROOT.gPad.Update()
                    ROOT.gPad.Modify()
                    self.canvases["threshold_"+str(detID)].Update()
                    self.canvases["threshold_"+str(detID)].Modify()
                    ROOT.gStyle.SetStatX(10)
                    ROOT.gStyle.SetStatY(10) 
                    mp = langau.GetParameter(1)
                    mp_langau_gr.SetPoint(np, fc+1,mp)
                    ex = 0
                    ey = langau.GetParError(1)
                    mp_langau_gr.SetPointError(np, ex, ey)
                    chi2 = langau.GetChisquare()
                    ndof = langau.GetNDF()+1E-12
                    control = chi2/ndof
                    if control > 200: control = 999
                    chi2_gr.SetPoint(int(np) , fc+1 , control)
                    np += 1

                Fit_Statistics_29.write("\\hline \n")
#                Fit_Statistics_29.write("\\caption {Fit Statistics of US Plane %s}" %l)
                Fit_Statistics_29.write("\n\\end{tabular}}\n")
                Fit_Statistics_29.write("\\end{table} \n")
 
                self.canvases["threshold_"+str(detID)].Print(self.runNr+"/threshold_"+str(detID)+".pdf")
        self.canvases["threshold_"+str(23005)].SaveAs("meeting_23003.root")
        ut.bookCanvas(self.canvases, "mp_all_channels", "Most Probable Values ", 1000,500,1,2)           
        self.canvases["mp_all_channels"].cd(1)
        mp_langau_gr.Draw('AP')
        self.canvases["mp_all_channels"].cd(2)
        chi2_gr.Draw("AP")
        Fit_Statistics_29.write( "\\end{document}")
        Fit_Statistics_29.close()
        Fit_Statistics_txt.close()
        print("Printing pdf output ...")
        os.system("pdflatex --output-directory " + self.runNr+ " "+self.runNr+"/Fit_Statistics_29_version2.tex")
        return mp_langau_gr, chi2_gr
