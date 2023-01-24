import ROOT,os,subprocess
import atexit
import time
import ctypes
from array import array
import os, psutil
import gc
ROOT.gStyle.SetTitleFont(102,"")
ROOT.gStyle.SetTitleFont(102,"xyz")
ROOT.gStyle.SetStatFont(102)
ROOT.gStyle.SetLabelFont(102,"xyz")
ROOT.gStyle.SetLegendFont(102)
ROOT.gStyle.SetHistLineColor(1)
ROOT.gStyle.SetHistLineWidth(3)

#ROOT.gROOT.SetBatch(ROOT.kTRUE)

# for fixing a root bug,  will be solved in the forthcoming 6.26 release.
ROOT.gInterpreter.Declare("""
#include "MuFilterHit.h"
#include "AbsMeasurement.h"
#include "TrackPoint.h"

void fixRoot(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllSignals(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRootT(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllTimes(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRoot(MuFilterHit& aHit, std::vector<TString>& key,std::vector<float>& value, bool mask) {
   std::map<TString, float> m = aHit.SumOfSignals();
   std::map<TString, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}

void fixRoot(std::vector<genfit::TrackPoint*>& points, std::vector<int>& d,std::vector<int>& k, bool mask) {
      for(std::size_t i = 0; i < points.size(); ++i) {
        genfit::AbsMeasurement*  m = points[i]->getRawMeasurement();
        d.push_back( m->getDetId() );
        k.push_back( int(m->getHitId()/1000) );
    }
}
""")
def pyExit():
       print("Make suicide until solution found for freezing")
       os.system('kill '+str(os.getpid()))
atexit.register(pyExit)

Tkey  = ROOT.std.vector('TString')()
Ikey   = ROOT.std.vector('int')()
Value = ROOT.std.vector('float')()

def map2Dict(aHit,T='GetAllSignals',mask=True):
     if T=="SumOfSignals":
        key = Tkey
     elif T=="GetAllSignals" or T=="GetAllTimes":
        key = Ikey
     else: 
           print('use case not known',T)
           1/0
     key.clear()
     Value.clear()
     if T=="GetAllTimes": ROOT.fixRootT(aHit,key,Value,mask)
     else:                         ROOT.fixRoot(aHit,key,Value,mask)
     theDict = {}
     for k in range(key.size()):
         if T=="SumOfSignals": theDict[key[k].Data()] = Value[k]
         else: theDict[key[k]] = Value[k]
     return theDict

import rootUtils as ut
import shipunit as u
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=True)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=False, default=False)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")
parser.add_argument("-npj", "--nPerJob",dest = "nPerJob",help="nbr of events per one job", default = 1E6)
parser.add_argument("-j","--nJob" ,dest = "nJob",help = "nbr of jobs",default=0)
parser.add_argument("-platform","--platform", dest = "platform", help = "lxplus or HTCondor",default = "lxplus" )
options = parser.parse_args()

runNr   = str(options.runNumber).zfill(6)
path     = options.path
partitions = 0
if options.runNumber > 0: 
 if path.find('eos')>0:
    path     = os.environ['EOSSHIP']+options.path
    dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+options.path,shell=True) )
# single file, before Feb'22
    data = "sndsw_raw_"+runNr+".root"
    if  dirlist.find(data)<0:
# check for partitions
       dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+options.path+"run_"+runNr,shell=True) )
       while 1>0:
        data = "raw-"+ str(partitions).zfill(4)
        if dirlist.find(data)>0:
            partitions+=1
        else: break
 else:
# check for partitions
       data = "sndsw_raw_"+runNr+".root"
       dirlist = os.listdir(options.path)
       if  not data in dirlist:
          dirlist  = os.listdir(options.path+"run_"+runNr)
          for x in dirlist:
            data = "raw-"+ str(partitions).zfill(4)
            if x.find(data)>0:
               partitions+=1

import SndlhcGeo
if not options.geoFile:
    if options.runNumber < 4575:
        options.geoFile =  "geofile_sndlhc_TI18_V3_08August2022.root"
    elif options.runNumber < 4855:
        options.geoFile =  "geofile_sndlhc_TI18_V5_14August2022.root"
    else:
        options.geoFile =  "geofile_sndlhc_TI18_V6_08October2022.root"

print("geofile :", options.geoFile)
if (options.geoFile).find('../')<0: geo = SndlhcGeo.GeoInterface(path+options.geoFile)
else:                                         geo = SndlhcGeo.GeoInterface(options.geoFile[3:])
MuFilter = geo.modules['MuFilter']
Scifi       = geo.modules['Scifi']
nav = ROOT.gGeoManager.GetCurrentNavigator()

A,B,locA,locB,globA,globB    = ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3()
latex = ROOT.TLatex()


# initialize 


if options.runNumber>0:
              eventChain = ROOT.TChain('rawConv')
              if partitions==0:
                   eventChain.Add(path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
              else:
                   for p in range(partitions):
                       eventChain.Add(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')

else:
# for MC data and other files
              f=ROOT.TFile.Open(options.fname)
              if f.Get('rawConv'):   eventChain = f.rawConv
              else:                        eventChain = f.cbmsim
#if options.remakeScifiClusters: eventChain.SetBranchStatus("Cluster_Scifi*",0)
rc = eventChain.GetEvent(0)
run      = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()
ioman.SetTreeName(eventChain.GetName())
outFile = ROOT.TMemFile('dummy','CREATE')
source = ROOT.FairFileSource(eventChain.GetCurrentFile())

if partitions>0:
    for p in range(1,partitions):
                       source.AddFile(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')
run.SetSource(source)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

houghTransform = False # under construction, not yet tested
if houghTransform:
   muon_reco_task = SndlhcMuonReco.MuonReco()
   muon_reco_task.SetName("houghTransform")
   run.AddTask(muon_reco_task)
else:
   import SndlhcTracking
   trackTask = SndlhcTracking.Tracking()
   trackTask.SetName('simpleTracking')
   run.AddTask(trackTask)

#avoiding some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

run.Init()
if partitions>0:  eventTree = ioman.GetInChain()
else:                 eventTree = ioman.GetInTree()

OT = ioman.GetSink().GetOutTree()
Reco_MuonTracks = trackTask.fittedTracks
Cluster_Scifi         = trackTask.clusScifi
# wait for user action 
   
class MuFilterEngine:
    "generate histograms for Muon Filter"
    ##### generate histograms for Muon system #####
    def __init__(self, options =  parser.parse_args()):
        print("Initializing Muon Filter Engine ...")
        print("Total nbr of events in this run : ", eventTree.GetEntries()/1E6, "M" )
        if options.path.find('TI18') or options.path.find('physics')>0: self.data = 'TI18'
        else : self.data = 'H8'
        self.nPerJob = int(options.nPerJob)
        self.nJob = int(options.nJob)
        self.nTot = eventTree.GetEntries() ### Total nbr of events
        self.nStart = self.nPerJob*self.nJob
        self.nEnd   = min(self.nTot,self.nStart + self.nPerJob)
        if options.platform == 'HTCondor' : self.loopOver = range(self.nStart,self.nEnd)
        else: self.loopOver = range(options.nEvents)
        self.hist = {}
        self.systemAndPlanes   = {1:2,2:5,3:7}
        self.systemAndBars     = {1:7,2:10,3:60}
        self.systemAndChannels = {1:16,2:16,3:2}
        self.sdict             = {1:'Veto',2:'US',3:'DS'}
        for system in self.systemAndPlanes:
            for plane in range(self.systemAndPlanes[system]):
                if system == 2:
                    ut.bookHist(self.hist,"pass-bars-"+str(plane),   " " , 22,0.,11)
                    ut.bookHist(self.hist,"total-bars-"  +str(plane)," " , 22,0.,11)
                    ut.bookHist(self.hist,"residuals-"+str(plane),"Residuals plane"+str(plane), 60,-30.,30.)
                    self.hist["residuals-"+str(plane)].GetXaxis().SetTitle("res [cm]")
                    ut.bookHist(self.hist,"nbrHits-"+str(plane)," nbr of hits ",10,0.,10.)
                    self.hist["nbrHits-"+str(plane)].GetXaxis().SetTitle("nbr of hits")
                orientation = self.systemAndOrientation(system,plane)
                for bar in range(self.systemAndBars[system]):
                    bar_key = self.sdict[system]+"plane"+str(plane)+orientation+"bar"+str(bar)
                    if system == 2:
                        ut.bookHist(self.hist,"pass-channels-"+bar_key , " ", 34, 0.,17)
                        ut.bookHist(self.hist,"total-channels-"+bar_key, " ", 34, 0.,17)
                        self.hist["pass-channels-"+bar_key].GetXaxis().SetTitle("fired_sipm")
                    for channel in range(self.systemAndChannels[system]):
                        if orientation == "vertical": side = ""
                        elif self.systemAndChannels[system]/2 > channel : side = "L"
                        else : side = "R" 
                        channel_key = bar_key+side+"largesipm"+str(channel)
                        if system == 2 and self.smallSiPMchannel(channel): channel_key = bar_key+side+"smallsipm"+str(channel)
                        ut.bookHist(self.hist, "qdc"+channel_key," ",200,-50.,200)
                        if orientation == "vertical" : break
        ut.bookHist(self.hist,"chi2Ndof", "",400,0,800)
        ut.bookHist(self.hist,"track_angle_x","", 100,0.,1.)
        ut.bookHist(self.hist,"track_angle_y","", 100,0.,1.)
        ut.bookHist(self.hist,"xEx_diff","",100,-25.,25.)
        ut.bookHist(self.hist,"yEx_diff","",100,-25.,25.)
        ut.bookHist(self.hist,"xvsy","",100,-25.,25.,100,-25.,25.)

        for key in self.hist:
            print(key)
        print("--------------------- Histograms booked -----------------------")


    def systemAndOrientation(self,s,plane):
        if s==1 or s==2: return "horizontal"
        if plane%2==1 or plane == 6: return "vertical"
        return "horizontal"

    def residual(self, theTrack, detID):
        MuFilter.GetPosition(detID,A,B)
        Z = 0.5*(A[2]+B[2])
        Y = 0.5*(A[1]+B[1])
        xEx, yEx = self.extrapolate(theTrack,Z)
        res = yEx-Y
        return res

    def expected_bars(self, theTrack):
        expected_detIDs = []
        for l in range(5):
            Y_in_planes={}
            for bar in range(10):
                detID=int(2E4+l*1E3+bar)
                MuFilter.GetPosition(detID,A,B)
                Z = (A[2]+B[2])/2.
                Y = (A[1]+B[1])/2.
                xEx,yEx = self.extrapolate(theTrack, Z) 
                Y_in_planes[abs(Y-yEx)]=detID
            sorted_keys = sorted(Y_in_planes.items())
            expected_detIDs.append(sorted_keys[0][1])
        return expected_detIDs

    def extrapolate(self, theTrack,z_mid):
        state = theTrack.getFittedState()
        pos   = state.getPos()
        mom = state.getMom()
        slope_x = mom.x()/mom.z()
        slope_y = mom.y()/mom.z()
        x=pos.x()+slope_x*(z_mid-pos.z())
        y=pos.y()+slope_y*(z_mid-pos.z())
        return x,y


    def check_scoring(self,theTrack,system,plane):
        reaching  = True
        orientation = self.systemAndOrientation(system,plane)
        bars = range(self.systemAndBars[system])
        if orientation == 'horizontal':
            MuFilter.GetPosition(int(system*10**4+plane*1E3+bars[0]),A,B)
            x_start, x_end = A[0], B[0]
            y_start = (A[1]+B[1])/2
            MuFilter.GetPosition(int(system*10**4+plane*1E3+bars[-1]),A,B)
            y_end = (A[1]+B[1])/2
        if orientation == 'vertical':
            MuFilter.GetPosition(int(system*10**4+plane*1E3+bars[0]),A,B)
            y_start, y_end = A[1], B[1]
            x_start = (A[0]+B[0])/2
            MuFilter.GetPosition(int(system*10**4+plane*1E3+bars[-1]),A,B)
            x_end = (A[0]+B[0])/2
        x_det = [x_start,x_end]
        y_det = [y_start,y_end]
        z_mid = (A[2]+B[2])/2
        x, y = self.extrapolate(theTrack,z_mid)
        reaching_x = (min(x_det) <= x and x <= max(x_det))
        reaching_y = (min(y_det) <= y and y <= max(y_det))
        reaching = reaching_x*reaching_y
        return reaching

    def chris(self, test_plane, observed_detIDs, expected_detIDs):
        detIDs_per_layer = {0:[],1:[],2:[],3:[],4:[]}
        chris = True
        if len(observed_detIDs)> 5 : chris = False
        for the_detID in observed_detIDs:
            l = (the_detID%10000)//1000
            detIDs_per_layer[l].append(the_detID)
        for key in detIDs_per_layer:
            if key == test_plane : continue
            if len(detIDs_per_layer[key])!=1 : chris = False
            if chris == False : break  ## break it already !
        return chris
      
    def check_veto(self):
        veto=False
        for hit in eventTree.Digi_MuFilterHits:
            if not hit.isValid() : continue
            detID = hit.GetDetectorID()
            system = hit.GetSystem()
            if system == 1 : 
                veto=True
                break
        return veto

    def smallSiPMchannel(self,i):
        if i==2 or i==5 or i==10 or i==13: return True
        else: return False

    def check_lr(self,aHit):       ### valid for US stations
        nSiPMs = aHit.GetnSiPMs()
        nSides = aHit.GetnSides()
        left_right={}
        allChannels= map2Dict(aHit,'GetAllSignals')
        for c in allChannels:
            Nleft = 0
            Nright = 0
            if nSiPMs>c : Nleft+=1
            else:Nright+=1
        left_right['left']=Nleft
        left_right['right']=Nright
        return left_right

    def isChannel_efficient(self,channel, aHit):
        allChannels = map2Dict(aHit,'GetAllSignals')
        eff = False
        if channel in allChannels : eff = True
        return eff
    
    def is_bar_efficient(self,expected_detID,observed_detIDs):
        tmp = []
        eff = -999999
        for i in range(-2,3):
            tmp.append(i+expected_detID)
        for bar in tmp:
            if bar in observed_detIDs : eff = bar
        if eff<0: print("inefficency at event ", eventTree.GetReadEvent(),"missing hit at: ", expected_detID, "hits at bars:", observed_detIDs)
        return eff

    def ringing_timestamp(self,aTree):
        isRinging = False
        h = aTree.EventHeader
        time_stamp = h.GetEventTime()
        t = time_stamp%(4*3564)/4
        if t > 1525 and t < 1535 : isRinging = True
        return isRinging

    def veto_eff(self,tree):
        for event in tree:
            Reco_MuonTracks.Clear()
            for aTrack in Reco_MuonTracks: aTrack.Delete()
            rc = trackTask.ExecuteTask("ScifiDS")
            if Reco_MuonTracks.GetEntries()>2 : continue
            TRACKS={}
            for track in Reco_MuonTracks:
                if track.GetUniqueID()==1 : TRACKS['Scifi']=track
                if track.GetUniqueID()==3 : TRACKS['DS']=track
            if 'Scifi' not in TRACKS :continue
            if  'DS'   not in TRACKS :continue
            both_on_plane = True
            for key in TRACKS:
                reach = self.check_scoring(TRACKS[key],1,1)
                if not reach :
                    both_on_plane = False
                    break
            if not both_on_plane : continue
            fitStatus = TRACKS['DS'].getFitStatus()
            chi2Ndof = fitStatus.getChi2()/(fitStatus.getNdf()+1E-10)
            MuFilter.GetPosition(int(1*10**4+1*1E3+4),A,B)
            z_mid = (A[2]+B[2])/2
            xEx_Scifi, yEx_Scifi = self.extrapolate(TRACKS['Scifi'],z_mid)
            xEx_DS, yEx_DS = self.extrapolate(TRACKS['DS'],z_mid)
            x_diff=xEx_Scifi-xEx_DS
            y_diff=yEx_Scifi-yEx_DS
            self.hist['xEx_diff'].Fill(x_diff)
            self.hist['yEx_diff'].Fill(y_diff)
            self.hist['xvsy'].Fill(x_diff,y_diff)
                                                
    def us_eff(self):
        background_region = {0:15.,1:15.,2:12.5,3:10.,4:10.}
        ROOT.EnableImplicitMT()
        newtree = eventTree.CloneTree(0)
        for event in self.loopOver: 
            if event > eventTree.GetEntries(): break # don't loopover the same event in condor !
            Reco_MuonTracks.Clear()
            for aTrack in Reco_MuonTracks: aTrack.Delete()
            eventTree.GetEvent(event)
            process = psutil.Process(os.getpid())
            if eventTree.GetReadEvent()%100000==0 : print("now event at",eventTree.GetReadEvent(), process.memory_info().rss/1024**2)
            optionTrack = options.trackType
            if optionTrack == 'ScifiDS' :rc = trackTask.ExecuteTask("ScifiDS")
            if optionTrack=='DS': rc = trackTask.ExecuteTask("DS")
            else: rc = trackTask.ExecuteTask("Scifi")
            if not Reco_MuonTracks.GetEntries() == 1 : continue ### dont forget this
            theTrack = Reco_MuonTracks[0]
            if not hasattr(theTrack,"getFittedState"): continue
            if not theTrack.getFitStatus().isFitConverged() and optionTrack!='DS':  continue # for H8 where only two planes / proj were avaiable
            if self.data == 'TI18' and not theTrack.getFitStatus().isFitConverged(): continue
            state = theTrack.getFittedState()
            pos   = state.getPos()
            mom = state.getMom()
            fitStatus = theTrack.getFitStatus()
            chi2Ndof = fitStatus.getChi2()/(fitStatus.getNdf()+1E-10)
            self.hist["chi2Ndof"].Fill(chi2Ndof)
            if fitStatus.getNdf()<2 : continue
            slope_x=mom.x()/mom.z()
            slope_y=mom.y()/mom.z()
            slope=ROOT.TMath.Sqrt(slope_x**2+slope_y**2)
            self.hist["track_angle_x"].Fill(slope_x)
            self.hist["track_angle_y"].Fill(slope_y)
            if abs(mom.x()/mom.z())>0.1: continue
            if abs(mom.y()/mom.z())>0.2: continue
            nHit_per_stat  = {0:0,1:0,2:0,3:0,4:0}
            reach=self.check_scoring(theTrack,2,0)
            if reach == False : continue
            veto_hits = {}
            us_hits =   {}
            ds_hits =   {}
            is_noisy = False
            for hit in eventTree.Digi_MuFilterHits:
                if not hit.isValid() : continue
                detID = hit.GetDetectorID()
                system = hit.GetSystem()
                l  = (detID%10000)//1000
                bar = detID%1000
                if system == 1: continue
                if system == 3:
                    pointsWithMeasurement=theTrack.getPointsWithMeasurement()
                    for point in pointsWithMeasurement:
                        rawM=point.getRawMeasurement()
                        if rawM.getDetId()==detID : ds_hits[hit]=detID
                if system == 2:
                    us_hits[hit]=detID
                    nHit_per_stat[l]+=1
                res = self.residual(theTrack,detID)
                self.hist["residuals-"+str(l)].Fill(res)
                if abs(res) > background_region[l] : is_noisy = True
            expected_detIDs = self.expected_bars(theTrack)
            if len(expected_detIDs) < 5 : print("a bug is here ", expected_detIDs)
            if is_noisy :continue
            observed_detIDs = us_hits.values()
            clean_mu_event = True
            for expected_detID in expected_detIDs :
                l = (expected_detID%10000)//1000
                chris_condition = self.chris(l, observed_detIDs, expected_detIDs)  ### where l is the test plane
                if chris_condition == False :
                    clean_mu_event = False
                    break
#            if sum(nHit_per_stat.values())<4 : continue
                self.hist["nbrHits-"+str(l)].Fill(nHit_per_stat[l])
                which_bar = self.is_bar_efficient(expected_detID,observed_detIDs)
                if which_bar > 0 :
                    bar = which_bar%1000
                    self.hist["pass-bars-"+str(l)].Fill(bar+1)
                    self.hist["total-bars-"+str(l)].Fill(bar+1)
                else :
                    clean_mu_event = False
                    bar = expected_detID%1000   #### blame the expected
                    self.hist["total-bars-"+str(l)].Fill(bar+1)
            if not clean_mu_event : continue
            us_hits.update(ds_hits)
            newtree.Fill()
            for hit in us_hits:
                detID = hit.GetDetectorID()
                system = hit.GetSystem()
                l  = (detID%10000)//1000
                bar = detID%1000
                if system>2:
                    l=2*l
                    if bar>59:
                        bar=bar-60
                        if l<6: l+=1
                signals = map2Dict(hit,"GetAllSignals",mask = True)
                times   = map2Dict(hit, "GetAllTimes", mask= True)
                orientation = self.systemAndOrientation(system,l)
                bar_key = self.sdict[system]+"plane"+str(l)+orientation+"bar"+str(bar)
                nSiPMs = hit.GetnSiPMs()
                nSides = hit.GetnSides()
                for key in signals:
                    if nSides > 1 and nSiPMs>key : side="L"
                    if nSides > 1 and nSiPMs<=key: side="R"
                    if nSides < 2 : side = ""
                    channel_key = bar_key+side+"largesipm"+str(key)
                    if system == 2 and self.smallSiPMchannel(key): channel_key = bar_key+side+"smallsipm"+str(key)
                    self.hist["qdc"+channel_key].Fill(signals[key])
                    if system != 2 : continue
                    if self.smallSiPMchannel(key):continue
                    if self.isChannel_efficient(key, hit) :
                        self.hist["pass-channels-"+bar_key].Fill(key+1)
                        self.hist["total-channels-"+bar_key].Fill(key+1)
                    else:
                        self.hist["total-channels-"+bar_key].Fill(key+1)
            theTrack.Delete()
            gc.collect()
        return newtree


    def write_to_file(self):
        if options.platform == "HTCondor":
            ut.writeHists(self.hist, "histograms_"+str(self.nJob)+"_"+str(self.data)+"_run"+str(options.runNumber)+"_"+str(options.trackType)+".root")
        else:
            ut.writeHists(self.hist, "histograms_"+str(self.data)+"_run"+str(options.runNumber)+"_"+str(options.trackType)+".root")
