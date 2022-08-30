import ROOT,os,subprocess
import atexit
import time
import ctypes
from array import array
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
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
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
# for MC data
              f=ROOT.TFile.Open(options.fname)
              eventChain = f.cbmsim
eventChain.GetEvent(0)

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

houghTransform = False#True # under construction, not yet tested
if houghTransform:
   import SndlhcMuonReco
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
# backward compatbility for early converted events
eventTree.GetEvent(0)
if eventTree.GetBranch('Digi_MuFilterHit'): eventTree.Digi_MuFilterHits = eventTree.Digi_MuFilterHit

OT = ioman.GetSink().GetOutTree()
Reco_MuonTracks = trackTask.fittedTracks

   
class EngineUS:
    "generate histograms for Upstream Muon Detector"
    ##### generate histograms related to the Upstream Muon Detector #####
    def __init__(self, options =  parser.parse_args()):
        print( "Total nbr of events in this run : ", eventTree.GetEntries())
        if options.path.find('TI18')>0: self.data = 'TI18'
        else : self.data = 'H8'
        self.nPerJob = int(options.nPerJob)
        self.nJob = int(options.nJob)
        self.nEvent = eventTree.GetEntries()
        self.nStart = self.nPerJob*self.nJob
        self.nEnd   = min(self.nEvent,self.nStart + self.nPerJob)
        if options.platform == 'HTCondor' : self.loopOver = range(self.nStart,self.nEnd)
        else: self.loopOver = range(options.nEvents)
        self.hist = {}
        for l in range(5):
            ut.bookHist(self.hist,"pass_bars_"+str(l),   " " , 22,0.,11)
            ut.bookHist(self.hist,"total_bars_"  +str(l)," " , 22,0.,11)
            ut.bookHist(self.hist,"residual_s_"+str(l),"plane "+str(l), 20,-30.,30.)
            self.hist["residual_s_"+str(l)].GetXaxis().SetTitle("res [cm]")
            ut.bookHist(self.hist,"occupancy_"+str(l)," nbr of hits ",10,0.,10.)
            self.hist["occupancy_"+str(l)].GetXaxis().SetTitle("nbr of hits")
            for bar in range(10):
                detID=int(2E4+l*1E3+bar)
                ut.bookHist(self.hist,"residual_"+str(detID),"residual station "+str(l)+"_bar_"+str(bar),40,-20.,20.)
                self.hist["residual_"+str(detID)].GetXaxis().SetTitle("residual [cm]")
                ut.bookHist(self.hist,"pass_channels_"+str(detID) , " ", 34, 0.,17)
                ut.bookHist(self.hist,"total_channels_"+str(detID), " ", 34, 0.,17)
                self.hist["pass_channels_"+str(detID)].GetXaxis().SetTitle("fired_sipm")
                for channel in range(16):
                    ut.bookHist(self.hist,"qdc_muon"+str(detID)+"_channel_"+ str(channel),"QDC Distribution at "+str(detID)+" " +str(channel),250,0.,100)
        print("----------histograms booked----------")

    
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

    def check_reach(self, theTrack,trackType,data):
        reaching  = True
        if trackType == 'DS' and data == 'H8':
            l, bar = 0, 0
            MuFilter.GetPosition(int(2E4+l*1E3+bar),A,B)
            z_0 = (A[2]+B[2])/2
            x, y = extrapolate(theTrack,z_0)
            x_det = [-74.36500144004822, 8.160000085830688]
            y_det = [0.9049997925758362, 54.99499979056418]        
        ####for TI18#####
        if trackType == 'DS' and data == 'TI18':
            z_0 = 379.4010
            x,y = self.extrapolate(theTrack,z_0)
            x_det = [-78.9938,3.5294]
            y_det = [15., 65.]
        if trackType == 'Scifi' and data == 'TI18': 
            z_0=470.3488
            x,y = self.extrapolate(theTrack,z_0)
            x_det = [-79.3248,3.1984]
            y_det = [7.9,68.0699]
        if x < x_det[0] or x > x_det[1]: reaching = False
        if y < y_det[0] or y > y_det[1]: reaching = False
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
        eff = -99999
        for i in range(-2,3):
            tmp.append(i+expected_detID)
        for bar in tmp:
            if bar in observed_detIDs : eff = bar
        if eff<0: print( "inefficency at event  ", eventTree.GetReadEvent(),"missing hit at bar: ", expected_detID, "hits at bars:", observed_detIDs )
        return eff

    def ringing_timestamp(self,aTree):
        isRinging = False
        h = aTree.EventHeader
        time_stamp = h.GetEventTime()
        t = time_stamp%(4*3564)/4
        if t > 1525 and t < 1535 : isRinging = True
        return isRinging


    def us_eff(self):
        background_region = {0:15.,1:15.,2:12.5,3:10.,4:10.}
        for event in self.loopOver:
            if event > eventTree.GetEntries(): break
            Reco_MuonTracks.Clear()
            for aTrack in Reco_MuonTracks: aTrack.Delete()
            if event > eventTree.GetEntries(): break
            eventTree.GetEvent(event)
#            if platform == 'lxplus' and eventTree.GetReadEvent() > options.nEvents : break
            if eventTree.GetReadEvent()%100000==0 : print("now event at",eventTree.GetReadEvent())
            optionTrack = options.trackType
            if optionTrack=='DS': rc = trackTask.ExecuteTask("DS")
            else: rc = trackTask.ExecuteTask("Scifi")
            if not Reco_MuonTracks.GetEntries() == 1 : continue
            theTrack = Reco_MuonTracks[0]
            if not hasattr(theTrack,"getFittedState"): continue
            if not theTrack.getFitStatus().isFitConverged() and optionTrack!='DS':  continue # for H8 where only two planes / proj were avaiable
            if self.data == 'TI18' and not theTrack.getFitStatus().isFitConverged(): continue
            state = theTrack.getFittedState()
            pos   = state.getPos()
            mom = state.getMom()
            fitStatus = theTrack.getFitStatus()
            chi2Ndof = fitStatus.getChi2()/(fitStatus.getNdf()+1E-10)
            if chi2Ndof>9: continue
            if abs(mom.x()/mom.z())>0.25:continue
            if abs(mom.y()/mom.z())>0.1: continue
            nHit_per_stat  = {0:0,1:0,2:0,3:0,4:0}
            reach=self.check_reach(theTrack,options.trackType, self.data)
            if reach == False : continue
            all_hits = {}  ### means all US hits
            is_noisy = False
            for hit in eventTree.Digi_MuFilterHits:
                if not hit.isValid() : continue
                detID = hit.GetDetectorID()
                system = hit.GetSystem()
                if system != 2 : continue
                l  = (detID%10000)//1000
                bar = detID%1000
                all_hits[hit]=detID
                nHit_per_stat[l]+=1
                theDict = map2Dict(hit,"GetAllSignals",mask = True)
                all_signals = hit.GetAllSignals()
                res = self.residual(theTrack,detID)
                self.hist["residual_s_"+str(l)].Fill(res)
                if abs(res) > background_region[l] : is_noisy = True
            expected_detIDs = self.expected_bars(theTrack)
            if len(expected_detIDs) < 5 : print("a bug is here ", expected_detIDs)
            if is_noisy :continue
            clean_mu_event = True
            observed_detIDs = all_hits.values()
            for expected_detID in expected_detIDs :
                l = (expected_detID%10000)//1000
                chris_condition = self.chris(l, observed_detIDs, expected_detIDs)  ### where l is the test plane
                if chris_condition == False :
                    clean_mu_event = False
                    break
#            if sum(nHit_per_stat.values())<4 : continue
                self.hist["occupancy_"+str(l)].Fill(nHit_per_stat[l])
                which_bar = self.is_bar_efficient(expected_detID,observed_detIDs)
                if which_bar > 0 :
                    bar = which_bar%1000
                    self.hist["pass_bars_"+str(l)].Fill(bar+1)
                    self.hist["total_bars_"+str(l)].Fill(bar+1)
                else :
                    clean_mu = False
                    bar = expected_detID%1000   #### blame the expected
                    self.hist["total_bars_"+str(l)].Fill(bar+1)
            if not clean_mu_event : continue        
            for hit in all_hits:
                detID = hit.GetDetectorID()
                l  = (detID%10000)//1000
                bar = detID%1000
                theDict = map2Dict(hit,"GetAllSignals",mask = True)
                for channel in range(16):
                    if self.smallSiPMchannel(channel):continue
                    if self.isChannel_efficient(channel, hit) : 
                        self.hist["pass_channels_"+str(detID)].Fill(channel+1)
                        self.hist["total_channels_"+str(detID)].Fill(channel+1)
                    else:
                        self.hist["total_channels_"+str(detID)].Fill(channel+1)
                for key in theDict:
                    self.hist["qdc_muon"+str(detID)+"_channel_"+str(key)].Fill(theDict[key])
            theTrack.Delete()


    def write_to_file(self):
        if options.platform == "HTCondor":
            ut.writeHists(self.hist, "histograms_"+str(self.nJob)+"_"+str(self.data)+"_run"+str(options.runNumber)+"_"+str(options.trackType)+".root")
        else:
            ut.writeHists(self.hist, "histograms_"+str(self.data)+"_run"+str(options.runNumber)+"_"+str(options.trackType)+".root")
