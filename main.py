from MuFilterEngine import MuFilterEngine
#from Analyzer import Analyzer
from argparse import ArgumentParser

#parser = ArgumentParser()
#parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
#options = parser.parse_args()
#fname = options.fname


def main():

    produce_histos  = MuFilterEngine()
    tree=produce_histos.us_eff()
    rc = produce_histos.veto_eff(tree)
    produce_histos.write_to_file()
#    ana = Analyzer(fname)
#    fit_qdc = ana.analyze_QDC()
#    eff = ana.analyze_eff()
#   rc = ana.write_to_file() 


if __name__ == "__main__":
    print("creating histograms for Muon Filter System")
    rc  = main()
    print("Finished taks")
