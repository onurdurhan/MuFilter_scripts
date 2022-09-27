from US_Engine import EngineUS


def main():

    produce_histos  = EngineUS()
    produce_histos.us_eff()
    produce_histos.write_to_file()

if __name__ == "__main__":
    print("creating histograms for US Muon detector")
    rc  = main()
