import ROOT,os,subprocess
import atexit
import time
from array import array
import rootUtils as ut
ROOT.gStyle.SetTitleFont(102,"")
ROOT.gStyle.SetTitleFont(102,"xyz")
ROOT.gStyle.SetStatFont(102)
ROOT.gStyle.SetLabelFont(102,"xyz")
ROOT.gStyle.SetLegendFont(102)
ROOT.gStyle.SetHistLineColor(1)
ROOT.gStyle.SetHistLineWidth(3)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

class Gaux:
    def __call__(self, x, par):
        #Fit parameters:
        #par[0]= Box_1  (constant)  box to be convoluted by gauss
        #par[1]= gaus scale par
        #par[2]= mean
        #par[3]= sigma
        #par[4]= constant (noise)
        np = 20. #50.      # number of convolution steps
        integral = 0 
#        x_low= -3
#        x_up =  3
        x_low = -3.
        x_up  = 3.
        dt = (x_up-x_low) / np
        i=0.0
        while i<=np:
            t = x_low+i*dt
            integral += par[0]*par[1]*ROOT.TMath.Gaus(x[0]-t,par[2],par[3],True)*dt
            i+=1.
        val = integral+par[4]
        if x[0]<x_up and x[0]>x_low: fit_val = val
        else :
            fit_val = par[1]*ROOT.TMath.Gaus(x[0],par[2],par[3],True)+par[4]
        return fit_val




class Langau:
    def __call__(self, x, par):
        #Fit parameters:
        #par[0]=Width (scale) parameter of Landau density
        #par[1]=Most Probable (MP, location) parameter of Landau density
        #par[2]=Total area (integral -inf to inf, normalization constant)
        #par[3]=Width (sigma) of convoluted Gaussian function
        #
        #In the Landau distribution (represented by the CERNLIB approximation),
        #the maximum is located at x=-0.22278298 with the location parameter=0.
        #This shift is corrected within this function, so that the actual
        #maximum is identical to the MP parameter.
#
        # Numeric constants
        invsq2pi = 0.3989422804014   # (2 pi)^(-1/2)
        mpshift  = -0.22278298       # Landau maximum location
#
        # Control constants
        np = 100.0      # number of convolution steps
        sc =   5.0      # convolution extends to +-sc Gaussian sigmas
#
        # Variables
        summe = 0.0
#
        # MP shift correction
        mpc = par[1] - mpshift * par[0]
#
        # Range of convolution integral
        xlow = max(0,x[0] - sc * par[3])
        xupp = x[0] + sc * par[3]
#
        step = (xupp-xlow) / np
#
        # Convolution integral of Landau and Gaussian by sum
        i=1.0
        if par[0]==0 or par[3]==0: return 9999
        while i<=np/2:
            i+=1
            xx = xlow + (i-.5) * step
            fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
            summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
#
            xx = xupp - (i-.5) * step
            fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
            summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
#
        return (par[2] * step * summe * invsq2pi / par[3])


class TwoLangaufun:
    def __call__(self,x, par):
        N1 = langaufun(x,par)
        par2 = [par[0],par[1]*2,par[4],par[3]]
        N2 = langaufun(x,par2)
        return N1+abs(N2)


gaux = Gaux()
langaufun = Langau()
twoLangaufun = TwoLangaufun()
class Analyzer:
    def __init__(self,fname):
        self.fname = fname
        self.canvases = {}
        self.f = ROOT.TFile(fname)
        if fname.find('TI18')>0 : data = 'TI18'
        else : data = 'H8'
        idx = fname.find('run')
        self.runNr = fname[idx:(idx+7)]
        self.outf = ROOT.TFile("analysis_derivative_"+data+"_"+self.runNr+".root","RECREATE")

    def fit_langau(self, hist, o, bmin,bmax,opt=''):
        if opt == 2:
            params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma',4:'N2'}
            F = ROOT.TF1('langau',langaufun,0,200,len(params))
        else: 
            params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma'}
            F = ROOT.TF1('langau',twoLangaufun,0,200,len(params))
        for p in params: F.SetParName(p,params[p])
        rc = hist.Fit('landau','S'+o,'',bmin,bmax)
        res = rc.Get()
        if not res: return res
        F.SetParameter(2,res.Parameter(0))
        F.SetParameter(1,res.Parameter(1))
        F.SetParameter(0,res.Parameter(2))
        F.SetParameter(3,res.Parameter(2))
        F.SetParameter(4,0)
        F.SetParLimits(0,0,100)
        F.SetParLimits(1,0,100)
        F.SetParLimits(3,0,10)
        rc = hist.Fit(F,'S'+o,'',bmin,bmax)
        res = rc.Get()
        return res

    def twoLangaufun(x,par):
        N1 = langaufun(x,par)
        par2 = [par[0],par[1]*2,par[4],par[3]]
        N2 = langaufun(x,par2)
        return N1+abs(N2)


    def fit_gaux(self, hist):
        params = {0:'Const(box)',1:'Scale(gaus)',2:'mean',3:'sigma',4:'Const(noise)'}
        binmax = hist.GetMaximumBin()
        binmin = hist.GetMinimumBin()
        x_max =  15.#hist.GetXaxis().GetBinCenter(binmax)
        x_min = -15#hist.GetXaxis().GetBinCenter(binmin)
        mean  =  hist.GetMean()
        sigma =  hist.GetStdDev()
        rc  = hist.Fit("gaus")
#        res = rc.Get()
        res = hist.GetFunction("gaus")
#    if not res: return res
#        gaux = Gaux()
        F = ROOT.TF1('gaux', gaux, x_min,x_max,len(params))
        print("fucntion is ", F)
        for p in params: F.SetParName(p,params[p])
        F.SetParameter(0,hist.GetMaximum())
        F.SetParameter(1,res.GetParameter(0))
        F.SetParameter(2,res.GetParameter(1))
        F.SetParameter(3,res.GetParameter(2))
        F.SetParameter(4,3.)

        F.SetParLimits(0,0,100)
        F.SetParLimits(3,0,10)
        hist.Fit("gaux")

    def smallSiPMchannel(self,i):
        if i==2 or i==5 or i==10 or i==13: return True
        else: return False
    

    def analyze_langau(self, fname):
        f = ROOT.TFile(fname)
        mp_langau_gr = ROOT.TGraphErrors()
        chi2_gr =  ROOT.TGraph()
        mp_langau_gr.SetMarkerStyle(21)
        chi2_gr.SetMarkerStyle(21)
        fc = 0
        np = 0
        os.system("mkdir "+self.runNr)
        import numpy as nmp
        domain = nmp.linspace(0,100.,101)
        cuts = {}
        Fit_Statistics_29 = open(self.runNr+"/Fit_Statistics_29_version2.tex", "w")
        Fit_Statistics_29.write("\\documentclass{article} \n")
        Fit_Statistics_29.write("\\usepackage{graphicx} \n")
        Fit_Statistics_29.write("\\begin{document} \n")
#        Fit_Statistics_29.write("\\begin{table}[ht] \n")
#        Fit_Statistics_29.write("\\centering \n")
#        Fit_Statistics_29.write("\\resizebox{\\textwidth}{!}{\\begin{tabular} { |l|| r | r | r | r | r  | r | r |} \n")
#        Fit_Statistics_29.write("\\multicolumn{8}{c}{Fit Statistics Variable 29} \\\\ \n")
#        Fit_Statistics_29.write("\\hline \n")
#        Fit_Statistics_29.write("Channel ID   & $Width (Scale)$ & $MPV$ & $Norm$ & $\sigma$ & $N2$ & $\chi^2/ndof $ & $Threshold$\\\\ \n")
#        Fit_Statistics_29.write( "\\hline \n" )

        for l in range(5):
            Fit_Statistics_29.write("\\clearpage \n")
            for bar in range(10):
                Fit_Statistics_29.write("\\begin{table}[ht] \n")
                Fit_Statistics_29.write("\\centering \n")
                Fit_Statistics_29.write("\\resizebox{\\textwidth}{!}{\\begin{tabular} { |l|| r | r | r | r | r  | r | r |} \n")
                Fit_Statistics_29.write("\\multicolumn{8}{c}{Fit Statistics : US Plane-%s Bar-%s} \\\\ \n" % (str(l),str(bar) ))
                Fit_Statistics_29.write("\\hline \n")
                Fit_Statistics_29.write("Channel ID   & $Width (Scale)$ & $MPV$ & $Norm$ & $\sigma$ & $N2$ & $\chi^2/ndof $ & $Threshold$\\\\ \n")
                Fit_Statistics_29.write( "\\hline \n" )
                detID=int(2E4+l*1E3+bar)
                ut.bookCanvas(self.canvases,"threshold_"+str(detID),"Langau Distributions with Threshold",700,500,4,4)
                canvas = f.Get("langau_muon"+str(detID))
                for c in range(16):
                    fc += 1
                    if self.smallSiPMchannel(c):continue
                    pad = canvas.GetPad(c+1)
                    h = pad.GetPrimitive("qdc_muon"+str(detID)+"_channel_"+str(c))
                    if not h : continue
                    langau = h.GetFunction("langau")
                    if not langau  :
                        print("No Fit result for ", detID, c)
                        continue
                    langau.SetNormalized(True)
                    for x in domain:
                        if langau.Derivative(x)>1E-5:
                            threshold = x
                            break
                    print("The ratio ", h.Integral(int(threshold+1),100)/h.Integral(0,100), detID, c )
                    if Fit_Statistics_29!= None:
                        Fit_Statistics_29.write(format("%s" % "DetID "+str(detID)+" Channel "+str(c)))
                        for par in range(langau.GetNpar()):
                            value = langau.GetParameter(par)
                            value_error = langau.GetParError(par)
                            Fit_Statistics_29.write( " & $ %.2f  \\pm  %.2f $" % (value, value_error))
                        chi2 = langau.GetChisquare()
                        ndof = langau.GetNDF()+1E-12
                        cndf = chi2/ndof
                        if cndf > 200 : cndf = 999
                        Fit_Statistics_29.write("& $%.2f$" % cndf)
                        Fit_Statistics_29.write("& $%.2f$" % threshold)
                        Fit_Statistics_29.write( " \\\\") ### end of row
                        Fit_Statistics_29.write( "\n")


                    langau.SetNormalized(False)
                    self.canvases["threshold_"+str(detID)].cd(c+1)
                    ROOT.gPad.Update()
                    rc = h.Draw()
                    self.canvases["threshold_"+str(detID)].Update()
                    cuts['threshold_'+str(detID)+"_"+str(c)] = ROOT.TLine(threshold, 0., threshold, 500)    
                    cuts['threshold_'+str(detID)+"_"+str(c)].SetLineColor(ROOT.kBlue)
                    cuts['threshold_'+str(detID)+"_"+str(c)].SetLineWidth(1)
                    cuts['threshold_'+str(detID)+"_"+str(c)].Draw("same")
                    ROOT.gPad.Update()
                    ROOT.gPad.Modify()
                    self.canvases["threshold_"+str(detID)].Update()
                    self.canvases["threshold_"+str(detID)].Modify()
                    mp = langau.GetParameter(1)
                    mp_langau_gr.SetPoint(np, fc+1,mp)
                    ex = 0
                    ey = langau.GetParError(1)
                    mp_langau_gr.SetPointError(np, ex, ey)
                    chi2 = langau.GetChisquare()
                    ndof = langau.GetNDF()+1E-12
                    control = chi2/ndof
                    if control > 200: control = 200
                    chi2_gr.SetPoint(int(np) , fc+1 , control)
                    np += 1


                Fit_Statistics_29.write("\\hline \n")
#                Fit_Statistics_29.write("\\caption {Fit Statistics of US Plane %s}" %l)
                Fit_Statistics_29.write("\n\\end{tabular}}\n")
                Fit_Statistics_29.write("\\end{table} \n")
 
                self.canvases["threshold_"+str(detID)].Print(self.runNr+"/threshold_"+str(detID)+".pdf")
        ut.bookCanvas(self.canvases, "mp_all_channels", "Most Probable Values ", 1000,500,1,2)           
        self.canvases["mp_all_channels"].cd(1)
        mp_langau_gr.Draw('AP')
        self.canvases["mp_all_channels"].cd(2)
        chi2_gr.Draw("AP")
#        Fit_Statistics_29.write("\\hline \n")
#        Fit_Statistics_29.write("\n\\end{tabular}}\n")
#        Fit_Statistics_29.write("\\end{table} \n")
        Fit_Statistics_29.write( "\\end{document}")
        Fit_Statistics_29.close()
        print("Printing pdf output ...")
        os.system("pdflatex --output-directory " + self.runNr+ " "+self.runNr+"/Fit_Statistics_29_version2.tex")
        return mp_langau_gr, chi2_gr


    def analyze_QDC(self):
        for l in range(5):
            for bar in range(10):
                detID=int(2E4+l*1E3+bar)
                ut.bookCanvas(self.canvases,"langau_muon"+str(detID),"Langau muon",700,500,4,4)              
        mp_langau_gr = ROOT.TGraphErrors()
        chi2_gr =  ROOT.TGraph()
        mp_langau_gr.SetMarkerStyle(21)
        chi2_gr.SetMarkerStyle(21)
        fc = 0
        np = 0
        for l in range(5):
            for bar in range(10):
                detID=int(2E4+l*1E3+bar)
                for c in range(16):
                    fc += 1
                    if self.smallSiPMchannel(c):continue
                    tmp = self.f.Get("qdc_muon"+str(detID)+"_channel_"+str(c))
                    bmin,bmax = 0,80
                    for k in range(tmp.GetMaximumBin(),1,-1):
                        if tmp.GetBinContent(k)<2:
                            bmin = k
                            break
                    print("Fitting ", detID ," channel ", c)
                    if tmp.GetEntries()<50 :
                        print("Not Enough Statistics: ",tmp.GetEntries() ,"at", detID, c, )
                        continue
                    self.canvases["langau_muon"+str(detID)].cd(c+1)
                    res = self.fit_langau(tmp,'L',0.8*tmp.GetBinCenter(bmin),1.5*tmp.GetBinCenter(bmax))
                    langau = tmp.GetFunction("langau")
                    print("Fitted")
                    tmp.Draw()
                    ROOT.gPad.Update()
                    if not langau  :
                        print("Fit Failed")
                        continue
                    mp = langau.GetParameter(1)
                    mp_langau_gr.SetPoint(np, fc+1,mp)
                    ex = 0
                    ey = langau.GetParError(1)
                    mp_langau_gr.SetPointError(np, ex, ey)
                    chi2 = langau.GetChisquare()
                    ndof = langau.GetNDF()+1E-12
                    control = chi2/ndof
                    if control > 200: control = 200
                    chi2_gr.SetPoint(np , fc+1 , control)
                    np += 1
        ut.bookCanvas(self.canvases, "mp_all_channels", "Most Probable Values ", 1000,500,1,2)           
        self.canvases["mp_all_channels"].cd(1)
        mp_langau_gr.Draw('AP')
        self.canvases["mp_all_channels"].cd(2)
        chi2_gr.Draw("AP")
        return mp_langau_gr, chi2_gr

    def analyze_us_res(self):
        ut.bookCanvas(self.canvases,"US_residuals","",1600,1000,3,2)
        ut.bookCanvas(self.canvases,"occupancies","",1600,1000,3,2)
        for l in range(5):
            self.canvases["US_residuals"].cd(l+1)
            h = self.f.Get("residual_s_"+str(l))
            h.SetLineWidth(3)
            h.GetXaxis().SetLabelSize(0.05)
            h.GetXaxis().SetTitleSize(0.05)
            h.GetYaxis().SetLabelSize(0.05)
            self.fit_gaux(h)
            h.Draw()
            ROOT.gStyle.SetOptFit(1111)
            h = self.f.Get("occupancy_"+str(l))
            self.canvases["occupancies"].cd(l+1)
            h.SetLineWidth(3)
            h.Draw()
 

    def analyze_eff(self):

        gr_bar_eff = {}
        gr_channel_eff = {}
        for i in range(5):
            tmp1 = self.f.Get("pass_bars_"+str(i))
            tmp2 = self.f.Get("total_bars_"  +str(i))
            eff = ROOT.TEfficiency(tmp1,tmp2)
            eff.Draw()
            ROOT.gPad.Update()
            gr_bar_eff[i]  = eff.GetPaintedGraph().Clone()
            gr_bar_eff[i].GetXaxis().SetTitle("bar ID")
            gr_bar_eff[i].GetXaxis().SetLabelSize(0.05)
            gr_bar_eff[i].GetXaxis().SetTitleSize(0.05)
            gr_bar_eff[i].GetXaxis().SetNdivisions(11)
            gr_bar_eff[i].GetYaxis().SetTitle("eff")
            gr_bar_eff[i].GetYaxis().SetTitleSize(0.05)
            gr_bar_eff[i].GetYaxis().SetLabelSize(0.05)
            gr_bar_eff[i].SetMarkerStyle(21)
            for j in range(gr_bar_eff[i].GetN()):
                gr_bar_eff[i].GetEXlow()[j] = 0
                gr_bar_eff[i].GetEXhigh()[j] = 0

            for bar in range(10):
                detID=int(2E4+i*1E3+bar)
                tmp1 = self.f.Get("pass_channels_"+str(detID))
                tmp2 = self.f.Get("total_channels_"+str(detID))
                eff = ROOT.TEfficiency(tmp1,tmp2)
                eff.Draw()
                ROOT.gPad.Update()
                gr_channel_eff[detID] = eff.GetPaintedGraph().Clone()
                gr_channel_eff[detID].GetXaxis().SetTitle("channel ID")
                gr_channel_eff[detID].GetXaxis().SetLabelSize(0.05)
                gr_channel_eff[detID].GetXaxis().SetTitleSize(0.05)
                gr_channel_eff[detID].GetXaxis().SetNdivisions(17)
                gr_channel_eff[detID].GetYaxis().SetTitle("eff")
                gr_channel_eff[detID].GetYaxis().SetTitleSize(0.05)
                gr_channel_eff[detID].GetYaxis().SetLabelSize(0.05)
                gr_channel_eff[detID].SetMarkerStyle(21)
                for k in range(gr_channel_eff[detID].GetN()):
                    gr_channel_eff[detID].GetEXlow()[k] = 0
                    gr_channel_eff[detID].GetEXhigh()[k] = 0

        ut.bookCanvas(self.canvases,"efficiencies_US_bars", "efficiencies",1600,1200, 3, 2)
        for i in range(5):        
            self.canvases["efficiencies_US_bars"].cd(i+1)
            gr_bar_eff[i].Draw("AP")
#            ROOT.gPad.Modified()
#            ROOT.gPad.Update()
            ut.bookCanvas(self.canvases,"efficiencies_US_channels_plane_"+str(i)," ", 1600,1000,4,3)
            for bar in range(10):
                detID = int(2E4+i*1E3+bar)
                self.canvases["efficiencies_US_channels_plane_"+str(i)].cd(bar+1)
                gr_channel_eff[detID].Draw("AP")
#                ROOT.gPad.Modified()
#                ROOT.gPad.Update()

    def write_to_file(self):    
        self.outf.cd()
        for canvas in self.canvases:
            print("printing......")
            self.canvases[canvas].Write()
        self.outf.Close()
