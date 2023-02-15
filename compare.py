import rootUtils as ut
import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)

danger = []

def analyze_QDC():
        hist_old= {}
        hist_recalib={}
        canvases={}
        fname_old = "histograms_oldcalib_TI18_run004819_DS.root"
        fname_new = "histograms_recalib_TI18_run004819_DS.root"
        f_old= ROOT.TFile("histograms_oldcalib_TI18_run004819_DS.root")
        f_new=ROOT.TFile("histograms_recalib_TI18_run004819_DS.root")
        fout = ROOT.TFile("old-recalib.root","RECREATE")
        ut.readHists(hist_old, fname_old)
        ut.readHists(hist_recalib, fname_new)
#        for histogram in hist_old:
#                    if histogram.find("qdc")<0 : continue
#                    ut.bookCanvas(canvases,histogram,"compare"+histogram)
        res = list(hist_old.keys())
        for histogram in res:
            if histogram.find('qdc')<0:
                hist_old.pop(histogram)
                hist_recalib.pop(histogram)
#            if hist_old[histogram].GetEntries()< 10 or hist_recalib[histogram].GetEntries()<10:
#                hist_old.pop(histogram)
#                hist_recalib.pop(histogram)

        c1 = ROOT.TCanvas("c1","sub data",200,10,700,500)
        for i, histogram in enumerate(hist_old):
#                    if histogram.find('qdc')<0 : continue
#                    old = f_old.Get(histogram)
#                    new = f_new.Get(histogram)
                    old = hist_old[histogram]
                    new = hist_recalib[histogram]
#                    canvases[histogram].cd()
#                    print("Fitted")
                    old.Draw("HIST")
                    new.Draw("HISTSAME")
                    old.SetLineColor(2)
                    new.SetLineStyle(2)
                    leg = ROOT.TLegend(.73,.32,.97,.53)
                    leg.SetBorderSize(0)
                    leg.SetFillColor(0)
                    leg.SetFillStyle(0)
                    leg.SetTextFont(42)
                    leg.SetTextSize(0.035)
                    leg.AddEntry(old,"old calib","L")
                    leg.AddEntry(new,"new calib","L")
#                    ROOT.gPad.Update()
#                    leg.Draw()
                    ROOT.gPad.Update()

                    c1.Update()
                    thresholds = {}
                    for x in range(old.GetNbinsX()):
                        if old.GetEntries()< 10:
                            threshold_old = 0 
                            break
                        ratio = old.Integral(x,200)/old.Integral(0,200)
                        if ratio <= 0.9999 :
                            threshold_old = old.GetBinCenter(x)
                            thresholds[old]=threshold_old
                            break
                    for x in range(new.GetNbinsX()):
                        if new.GetEntries()<10:
                            threshold_new = 0.
                            break
                        ratio = new.Integral(x,200)/new.Integral(0,200)
                        if ratio <= 0.9999 :
                            threshold_new = new.GetBinCenter(x)
                            thresholds[new]=threshold_new
                            break
                    if old.GetEntries()>1000 and abs(threshold_new-threshold_old)<1. : danger.append(histogram)
                    threshold_line = {}
                    j = 0
                    for key in thresholds:
                        threshold_line[key] = ROOT.TLine(thresholds[key], 0., thresholds[key], key.GetMaximum()/2)    
                        if j<1: 
                            threshold_line[key].SetLineColor(ROOT.kBlue)
                            threshold_line[key].SetLineStyle(1)
                        if j>0:
                            threshold_line[key].SetLineColor(ROOT.kMagenta)
                            threshold_line[key].SetLineStyle(2)
                        threshold_line[key].SetLineWidth(2)
                        threshold_line[key].Draw("same")

                        if key == old: leg.AddEntry(threshold_line[key],"old threshold","L")
                        if key == new: leg.AddEntry(threshold_line[new],"new threshold","L")
                        leg.Draw()
 
                        ROOT.gPad.Update()
                        ROOT.gPad.Modify()
                        c1.Update()
                        c1.Modify()
                        j+=1


#                    ROOT.gPad.Modify()
#                    canvases[histogram].Update()
#                    canvases[histogram].Modify()
                    if i == 0:
                        c1.Print("h1.pdf(","pdf")
                        print(i, "first")
                    elif i !=0 and i != len(hist_old)-1 :
                        c1.Print("h1.pdf","pdf")
                        print(i,"mid")
                    else : 
                        c1.Print("h1.pdf)","pdf")
                        print(i,"last")

#        fout.cd()
        for canvas in canvases:
            canvases[canvas].cd()
            leg = ROOT.TLegend(.73,.32,.97,.53)
            leg.SetBorderSize(0)
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.035)
            leg.AddEntry(old,"old calib","L")
            leg.AddEntry(new,"new calib","L")
#                    ROOT.gPad.Update()
            leg.Draw()
            fout.cd()
            canvases[canvas].Write()
        fout.Close()
        with open('danger.txt', 'w') as f:
            for thr in danger:
                f.write(thr)
                f.write('\n')
