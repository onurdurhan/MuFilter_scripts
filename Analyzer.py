import ROOT,os,subprocess
import atexit
import time
from array import array
import rootUtils as ut
ROOT.gStyle.SetTitleFont(102,"")
ROOT.gStyle.SetTitleFont(102,"xyz")
ROOT.gStyle.SetTextFont(10)
ROOT.gStyle.SetStatFont(102)
ROOT.gStyle.SetLabelFont(102,"xyz")
ROOT.gStyle.SetLegendFont(102)
ROOT.gStyle.SetHistLineColor(1)
ROOT.gStyle.SetHistLineWidth(3)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
#ROOT.gStyle.SetPadLeftMargin(-0.5)
#ROOT.gStyle.SetPadRightMargin(-0.5)
#ROOT.gPad.SetMargin(0.02, 0.001, 0.04, 0.02)

class Gaux:
    def __call__(self, x, par):
        #Fit parameters:
        #par[0]= Box_1  (constant)  box to be convoluted by gauss
        #par[1]= gaus scale par
        #par[2]= mean
        #par[3]= sigma
        #par[4]= constant (noise)
        np = 50.      # number of convolution steps
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

class ExpoLandau:
    def __call__(self,x,par):
        mpshift  = -0.22278298
        mpc = par[1] - mpshift * par[0]
#        Landau(Double_t x, Double_t mpv = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
        land = ROOT.TMath.Landau(x[0],mpc,par[0])*par[2] #,par[0])/par[0]
        expo = ROOT.TMath.Exp(-par[3]*x[0])*par[4]
        expoland = expo+land
        return expoland

gaux = Gaux()
langaufun = Langau()
twoLangaufun = TwoLangaufun()
expoland = ExpoLandau()

class Analyzer:
    def __init__(self,fname):
        self.fname = fname
        self.canvases = {}
        self.f = ROOT.TFile(fname)
        if fname.find('TI18')>0 : data = 'TI18'
        else : data = 'H8'
        idx = fname.find('00')
        self.runNr = fname[idx:(idx+6)]
        self.outf = ROOT.TFile("analysis_"+data+"_run_"+self.runNr+".root","RECREATE")

    def fit_langau(self, hist, o, bmin,bmax):
        params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma'} #,4:'N2'}
        F = ROOT.TF1('langau',langaufun,-10.,200 , len(params))
        print(F.GetNpar(),len(params))
        for p in params: F.SetParName(p,params[p])
        rc = hist.Fit('landau','S'+o,'',bmin,bmax)
        res = rc.Get()
        if not res: return res
        F.SetParameter(2,res.Parameter(0))
        F.SetParameter(1,res.Parameter(1))
        F.SetParameter(0,res.Parameter(2))
        F.SetParameter(3,res.Parameter(2))
        F.SetParLimits(0,0,100)
        F.SetParLimits(1,0,100)
        F.SetParLimits(3,0,10)
        rc = hist.Fit(F,'S'+o,'',bmin,bmax)

        ChiSqr = F.GetChisquare()
        NDF = F.GetNDF()
        
        print(ChiSqr,NDF, ChiSqr/NDF)

        res = rc.Get()
        return res


    def fit_expoland(self,hist,o,bmin,bmax):
        params = {0:'Width(scale)',1:'mostProbable',2:'NormLandau',3:'Tau',4:'NormExpo'} 
        F = ROOT.TF1('langau',expoland,-20.,200.,len(params))
        for p in params: F.SetParName(p,params[p])
        F.SetParameter(0,3)
        F.SetParameter(1,hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()))
        F.SetParameter(2,1)
        F.SetParameter(3,1)
        F.SetParameter(4,1)
        F.SetParLimits(0,0,10)
        F.SetParLimits(1,-10, 200)
        F.SetParLimits(3,0,100)
#        F.SetParLimits(4,0,1E6)
#        F.SetParLimits(3,0,hist.GetEntries()) 
        rc = hist.Fit(F,'S'+o,'',bmin,bmax)

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
    

    def fit_convolution(self, hist, o, bmin,bmax):
#        f_conv.SetNofPointsFFT(1000)
        F = ROOT.TF1("langau", f_conv, -10., 200., f_conv.GetNpar())
        params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma',4:'N2'}
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

    def analyze_QDC(self):
        hist= {}
        ut.readHists(hist,self.fname)             
        mp_langau_gr = ROOT.TGraphErrors()
        chi2_gr =  ROOT.TGraph()
        mp_langau_gr.SetMarkerStyle(21)
        chi2_gr.SetMarkerStyle(21)
        fc = 0
        np = 0
        for histogram in hist:
                    if histogram.find("qdcDS")<0 : continue
                    ut.bookCanvas(self.canvases, histogram,histogram)
        for histogram in hist:
                    if histogram.find("qdcDS")<0 : continue
                    fc+=1
                    tmp = self.f.Get(histogram)
                    tmp.SetTitle(histogram)
                    print("nbr of entries in this ", tmp.GetEntries())
                    bmin,bmax = 0,200
                    print(tmp.GetMaximumBin())
                    for k in range(tmp.GetMaximumBin(),1,-1):
                        if tmp.GetBinContent(k)<2:
                            bmin = k
                            break
                    print("Fitting ", histogram)
                    if tmp.GetEntries() < 100 :
                        print("Not Enough Statistics: ",tmp.GetEntries() ,"at", histogram )
                        continue
                    self.canvases[histogram].cd()
                    if  histogram.find('qdcDS') > -1:
                        print("found histogram ", histogram)
                        res =self.fit_expoland(tmp,'LM',0.8*tmp.GetBinCenter(bmin),1.5*tmp.GetBinCenter(bmax))
                        langau = tmp.GetFunction("langau")
                    else :                    
                        rc = tmp.Fit('landau',"SLM",'',0.8*tmp.GetBinCenter(bmin),1.5*tmp.GetBinCenter(bmax))
                        langau = tmp.GetFunction("landau")
                    print("Fitted")
                    tmp.Draw()
#                    ROOT.gPad.Modify()
                    ROOT.gPad.Update()
 #                   self.canvases[histogram].Update()
                    if not langau  :
                        print("Fit Failed for histogram ", histogram )
                        continue
                    mp = langau.GetParameter(1)
                    mp_langau_gr.SetPoint(np, fc+1,mp)
                    ex = 0
                    ey = langau.GetParError(1)
                    mp_langau_gr.SetPointError(np, ex, ey)
                    chi2 = langau.GetChisquare()
                    ndof = langau.GetNDF() #+1E-12
                    control = chi2/ndof
#                    if control > 200: control = 200
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
            gr_bar_eff[i].GetXaxis().SetLabelSize(0.035)
            gr_bar_eff[i].GetXaxis().SetTitleSize(0.035)
            gr_bar_eff[i].GetXaxis().SetNdivisions(12)
            gr_bar_eff[i].GetXaxis().SetRangeUser(0,10.9)
            gr_bar_eff[i].GetYaxis().SetTitle("eff")
            gr_bar_eff[i].GetYaxis().SetTitleSize(0.035)
            gr_bar_eff[i].GetYaxis().SetLabelSize(0.035)
            gr_bar_eff[i].GetYaxis().SetTitleOffset(1)
#            gr_bar_eff[i].GetYaxis().SetRangeUser(0.9,1.01)
            gr_bar_eff[i].SetMarkerStyle(21)
            gr_bar_eff[i].SetTitle("US Plane " + str(i))
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
                gr_channel_eff[detID].GetXaxis().SetLabelSize(0.035)
                gr_channel_eff[detID].GetXaxis().SetTitleSize(0.035)
                gr_channel_eff[detID].GetXaxis().SetNdivisions(17)
                gr_channel_eff[detID].GetXaxis().SetRangeUser(0,16.9)
                gr_channel_eff[detID].GetYaxis().SetTitle("eff")
                gr_channel_eff[detID].GetYaxis().SetTitleSize(0.035)
                gr_channel_eff[detID].GetYaxis().SetLabelSize(0.035)
                gr_channel_eff[detID].GetYaxis().SetTitleOffset(1)
                gr_channel_eff[detID].SetMarkerStyle(21)
                gr_channel_eff[detID].SetTitle("US Plane " + str(i)+" bar "+str(detID))
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

            self.canvases["efficiencies_US_channels_plane_"+str(i)].cd(11)        
            t = ROOT.TLatex()
            t.SetNDC()
            t.SetTextFont( 102 )
            t.SetTextColor( 1 )
            t.SetTextSize( 0.08 )
            t.SetTextAlign( 12 )
            t.DrawLatex( 0.4, 0.6, 'RUN '+str(self.runNr))
     

        self.canvases["efficiencies_US_bars"].cd(6)
        t = ROOT.TLatex()
        t.SetNDC()
        t.SetTextFont( 102 )
        t.SetTextColor( 1 )
        t.SetTextSize( 0.08 )
        t.SetTextAlign( 12 )
        t.DrawLatex( 0.4, 0.6, 'RUN '+str(self.runNr))
 
    def write_to_file(self):    
        self.outf.cd()
        for canvas in self.canvases:
            print("printing......")
            self.canvases[canvas].Write()
        self.outf.Close()
