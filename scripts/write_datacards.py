#!/bin/env python


from __future__ import absolute_import
from __future__ import print_function

import argparse
import sys
import os
from errors import E
import datacard_writer as dw
from six.moves import range
from six.moves import zip


parser = argparse.ArgumentParser(description="Plotter for the WVZ analysis")
parser.add_argument('-b' , '--baseline_tag'    , dest='baseline_tag'    , help='baseline tag (e.g. test, test1. test2, etc.)' , required=True)
parser.add_argument('-t' , '--ntuple_type'     , dest='ntuple_type'     , help='WVZ, Trilep, etc.'                            , required=True)
parser.add_argument('-v' , '--ntuple_version'  , dest='ntuple_version'  , help='v0.1.6, v0.1.7, etc.'                         , required=True)
parser.add_argument('-1' , '--emu_one_bin'     , dest='emu_one_bin'     , help='write emu in one bin'                         , action='store_true', default=False)
parser.add_argument('-4' , '--emu_four_bin'    , dest='emu_four_bin'    , help='write emu in four bin'                        , action='store_true', default=True)
parser.add_argument('-3' , '--eemm_three_bin'  , dest='eemm_three_bin'  , help='write eemm in three bin'                      , action='store_true', default=False)
parser.add_argument('-y' , '--print_yields'    , dest='print_yields'    , help='to print wysiwyg yields'                      , action='store_true', default=False)
parser.add_argument('-d' , '--print_detail'    , dest='print_detail'    , help='to print wysiwyg uncertainty detail'          , action='store_true', default=False)
parser.add_argument('-w' , '--wwz_only'        , dest='wwz_only'        , help='to write wwz only signal datacards'           , action='store_true', default=False)
parser.add_argument('-s' , '--nonh_vh_split'   , dest='nonh_vh_split'   , help='to write with vh and vvv split'               , action='store_true', default=False)

args = parser.parse_args()

import pyrootutil as pr
import plottery_wrapper as p
import sys
import os
from errors import E
import datacard_writer as dw

import ROOT as r

def main():

    sample2016 = "{}2016_{}".format(args.ntuple_type, args.ntuple_version)
    sample2017 = "{}2017_{}".format(args.ntuple_type, args.ntuple_version)
    sample2018 = "{}2018_{}".format(args.ntuple_type, args.ntuple_version)
    tag2016 = "y2016_" + args.baseline_tag
    tag2017 = "y2017_" + args.baseline_tag
    tag2018 = "y2018_" + args.baseline_tag
    write_datacards(sample2016, tag2016)
    write_datacards(sample2017, tag2017)
    write_datacards(sample2018, tag2018)
    write_datacards("_".join([ sample2016, sample2017, sample2018 ]), "_".join([ tag2016, tag2017, tag2018 ]))

def write_datacards(ntuple_version, tag):

    # ntuple_version = args.sample_set_name
    # tag = args.tag

    wwzonlysuffix = ""

    if args.wwz_only:
        wwzonlysuffix = "wwzonly"
        if args.nonh_vh_split:
            fname_sig     = "outputs/{}/{}/nonh_wwz.root".format(ntuple_version, tag)
        else:
            fname_sig     = "outputs/{}/{}/wwz.root".format(ntuple_version, tag)
    else:
        fname_sig     = "outputs/{}/{}/sig.root".format(ntuple_version, tag)
    fname_wwz     = "outputs/{}/{}/wwz.root".format(ntuple_version, tag)
    fname_wzz     = "outputs/{}/{}/wzz.root".format(ntuple_version, tag)
    fname_zzz     = "outputs/{}/{}/zzz.root".format(ntuple_version, tag)
    fname_ttz     = "outputs/{}/{}/ttz.root".format(ntuple_version, tag)
    fname_zz      = "outputs/{}/{}/zz.root".format(ntuple_version, tag)
    fname_wz      = "outputs/{}/{}/wz.root".format(ntuple_version, tag)
    fname_twz     = "outputs/{}/{}/twz.root".format(ntuple_version, tag)
    fname_rare    = "outputs/{}/{}/rare.root".format(ntuple_version, tag)
    fname_dyttbar = "outputs/{}/{}/dyttbar.root".format(ntuple_version, tag)
    fname_higgs   = "outputs/{}/{}/higgs.root".format(ntuple_version, tag)
    fname_othernoh= "outputs/{}/{}/othernoh.root".format(ntuple_version, tag)
    fname_data    = "outputs/{}/{}/data.root".format(ntuple_version, tag)
    fname_nonh_wwz= "outputs/{}/{}/nonh_wwz.root".format(ntuple_version, tag)
    fname_nonh_wzz= "outputs/{}/{}/nonh_wzz.root".format(ntuple_version, tag)
    fname_nonh_zzz= "outputs/{}/{}/nonh_zzz.root".format(ntuple_version, tag)
    fname_zh_wwz  = "outputs/{}/{}/zh_wwz.root".format(ntuple_version, tag)
    fname_wh_wzz  = "outputs/{}/{}/wh_wzz.root".format(ntuple_version, tag)
    fname_zh_zzz  = "outputs/{}/{}/zh_zzz.root".format(ntuple_version, tag)

    year = "2" + ntuple_version.split("_")[0].split("2")[1]
    if "2016" in ntuple_version and "2017" in ntuple_version:
        year = "All"
    prefix = "{}/{}".format(ntuple_version, tag)

    procs = ["data_obs", "sig", "ttz", "zz", "wz", "twz", "rare", "dyttbar", "higgs"]
    mcprocs = procs[1:]
    bkgprocs = procs[2:]
    fnames = [fname_data, fname_sig, fname_ttz, fname_zz, fname_wz, fname_twz, fname_rare, fname_dyttbar, fname_higgs]
    nonzzbkg = [fname_sig, fname_ttz, fname_wz, fname_twz, fname_rare, fname_dyttbar, fname_higgs]
    nonttzbkg = [fname_sig, fname_zz, fname_wz, fname_twz, fname_rare, fname_dyttbar, fname_higgs]

    procs = ["data_obs", "sig", "ttz", "zz", "wz", "twz", "higgs", "other"]
    mcprocs = procs[1:]
    bkgprocs = procs[2:]
    fnames = [fname_data, fname_sig, fname_ttz, fname_zz, fname_wz, fname_twz, fname_higgs, fname_othernoh]
    nonzzbkg = [fname_sig, fname_ttz, fname_wz, fname_twz, fname_higgs, fname_othernoh]
    nonttzbkg = [fname_sig, fname_zz, fname_wz, fname_twz, fname_higgs, fname_othernoh]

    if args.wwz_only:
        if args.nonh_vh_split:
            procs = ["data_obs", "sig", "zhwwz", "nonhwzz", "whwzz", "nonhzzz", "zhzzz", "zz", "ttz", "twz", "wz", "higgs", "other"]
            mcprocs = procs[1:]
            bkgprocs = procs[2:]
            fnames =    [ fname_data , fname_nonh_wwz , fname_zh_wwz , fname_nonh_wzz , fname_wh_wzz , fname_nonh_zzz , fname_zh_zzz , fname_zz  , fname_ttz , fname_twz , fname_wz  , fname_higgs , fname_othernoh]
            nonzzbkg =  [              fname_nonh_wwz , fname_zh_wwz , fname_nonh_wzz , fname_wh_wzz , fname_nonh_zzz , fname_zh_zzz ,             fname_ttz , fname_twz , fname_wz  , fname_higgs , fname_othernoh]
            nonttzbkg = [              fname_nonh_wwz , fname_zh_wwz , fname_nonh_wzz , fname_wh_wzz , fname_nonh_zzz , fname_zh_zzz , fname_zz  ,             fname_twz , fname_wz  , fname_higgs , fname_othernoh]
        else:
            procs = ["data_obs", "sig", "wzz", "zzz", "zz", "ttz", "twz", "wz", "higgs", "other"]
            mcprocs = procs[1:]
            bkgprocs = procs[2:]
            fnames =    [ fname_data , fname_wwz , fname_wzz , fname_zzz , fname_zz  , fname_ttz , fname_twz , fname_wz  , fname_higgs , fname_othernoh]
            nonzzbkg =  [              fname_wwz , fname_wzz , fname_zzz ,             fname_ttz , fname_twz , fname_wz  , fname_higgs , fname_othernoh]
            nonttzbkg = [              fname_wwz , fname_wzz , fname_zzz , fname_zz  ,             fname_twz , fname_wz  , fname_higgs , fname_othernoh]


    systcategs = ["BTagHF", "BTagLF", "JES", "Pileup", "Qsq", "PDF", "AlphaS", "MET", "JER", "METPileup","ElLepSF","MuLepSF"] # Null string is the nominal variation
    systnames = ["Nominal"] # Nominal always exist
    for systcateg in systcategs:
        systnames.append(systcateg+"Up")
        systnames.append(systcateg+"Down")

    #############
    # Open TFiles
    #############
    tfiles = {}
    for proc, fname in zip(procs, fnames):
        tfiles[proc] = r.TFile(fname)

    r.gROOT.cd()

    ###########################
    # OnZ Control region yields
    ###########################

    onz_zz_h = pr.get_summed_histogram([fname_zz], "ChannelOnZ__Yield")
    onz_data_h = pr.get_summed_histogram([fname_data], "ChannelOnZ__Yield")
    onz_nonzz_h = pr.get_summed_histogram(nonzzbkg, "ChannelOnZ__Yield")
    zz_sf = pr.get_sf(onz_zz_h, onz_data_h, onz_nonzz_h).GetBinContent(1)
    zz_sferr = pr.get_sf(onz_zz_h, onz_data_h, onz_nonzz_h).GetBinError(1)
    expected_nevt_zz = onz_data_h.GetBinContent(1)

    ############################
    # BTag Control region yields
    ############################

    bcr_ttz_h = pr.get_summed_histogram([fname_ttz], "ChannelBTagEMu__Yield")
    bcr_data_h = pr.get_summed_histogram([fname_data], "ChannelBTagEMu__Yield")
    bcr_nonttz_h = pr.get_summed_histogram(nonttzbkg, "ChannelBTagEMu__Yield")
    ttz_sf = pr.get_sf(bcr_ttz_h, bcr_data_h, bcr_nonttz_h).GetBinContent(1)
    ttz_sferr = pr.get_sf(bcr_ttz_h, bcr_data_h, bcr_nonttz_h).GetBinError(1)
    expected_nevt_ttz = bcr_data_h.GetBinContent(1)

    if not args.print_yields:
        print(year, "ttz_sf", "{:.2f} +/- {:.2f}".format(ttz_sf, ttz_sferr), expected_nevt_ttz)
        print(year, "zz_sf", "{:.2f} +/- {:.2f}".format(zz_sf, zz_sferr), expected_nevt_zz)

    ###############################
    # EMu channel data card writing
    ###############################

    # number of bins
    fitreg = "EMuHighMT"
    if args.emu_one_bin:
        nbins = 1
        fitvar = "Yield"
    else:
        if args.emu_four_bin:
            nbins = 4
            fitvar = "emuSR"
        else:
            nbins = 5
            fitvar = "MllNom"
        # nbins = 5
        # fitvar = "pt_zeta"

    # Main data base to hold all the histograms
    hists_db = {}

    # Loop over the processes
    for proc in procs:

        # Retrieve the tfile
        tfile = tfiles[proc]

        # For each processes create another map to hold various histograms
        hists_db[proc] = {}

        # Loop over the systematic variations
        for syst in systnames:

            if syst == "Nominal":
                if nbins == 5:
                    print(fitreg, fitvar)
                    print(tfile)
                    h = rebin36(tfile.Get("Channel{}__{}".format(fitreg, fitvar)).Clone())
                else:
                    h = tfile.Get("Channel{}__{}".format(fitreg, fitvar)).Clone()
            else:
                systhacked = syst
                if proc == "NONE":
                    systhacked = ""
                if nbins == 5:
                    h = rebin36(tfile.Get("Channel{}{}__{}".format(fitreg, systhacked, fitvar)).Clone())
                else:
                    h = tfile.Get("Channel{}{}__{}".format(fitreg, systhacked, fitvar)).Clone()
                # if nbins == 5:
                #     h = rebin36(tfile.Get("Channel{}{}__{}".format(fitreg, syst, fitvar)).Clone())
                # else:
                #     h = tfile.Get("Channel{}{}__{}".format(fitreg, syst, fitvar)).Clone()

            h.SetTitle("emu{}_{}".format(year, proc))

            if proc == "ttz":
                before_scale = h.Integral()
                h.Scale(ttz_sf)
                after_scale = h.Integral()
                if syst == "Nominal":
                    print(year, "ttz", before_scale, after_scale)
            if proc == "zz":
                before_scale = h.Integral()
                h.Scale(zz_sf)
                after_scale = h.Integral()
                if syst == "Nominal":
                    print(year, "zz", before_scale, after_scale)

            hists_db[proc][syst] = h

    systs = []

    # ZZ CR systematic line
    onz_cr_hist = r.TH1F("onz_cr", "", nbins, 0, nbins)
    for i in range(1, nbins+1):
        onz_cr_hist.SetBinContent(i, expected_nevt_zz)
    alpha = hists_db["zz"]["Nominal"].Clone("alpha")
    alpha.Divide(onz_cr_hist)
    thissyst = {}
    for proc in mcprocs:
       if proc == "zz":
           thissyst["emu{}_".format(year) + proc] = [ "{:4f}".format(alpha.GetBinContent(i)) for i in range(1,nbins+1) ]
       else:
           thissyst["emu{}_".format(year) + proc] = 0
    systs.append(("CRZZ{}".format(year), "gmN", [onz_cr_hist], thissyst))

    # ttZ CR systematic line
    btag_cr_hist = r.TH1F("btag_cr", "", nbins, 0, nbins)
    for i in range(1, nbins+1):
        btag_cr_hist.SetBinContent(i, expected_nevt_ttz)
    alpha = hists_db["ttz"]["Nominal"].Clone("alpha")
    alpha.Divide(btag_cr_hist)
    thissyst = {}
    for proc in mcprocs:
       if proc == "ttz":
           thissyst["emu{}_".format(year) + proc] = [ "{:4f}".format(alpha.GetBinContent(i)) for i in range(1,nbins+1) ]
       else:
           thissyst["emu{}_".format(year) + proc] = 0
    systs.append(("CRTTZ{}".format(year), "gmN", [btag_cr_hist], thissyst))

    # Experimental systematics
    for systcateg in systcategs:
        thissyst = {}
        for proc in mcprocs:
            # if proc not in ["zz", "ttz", "wz", "twz"]:
            if proc not in ["zz", "ttz", "wz"]:
                thissyst["emu{}_".format(year) + proc] = [hists_db[proc][systcateg+"Up"], hists_db[proc][systcateg+"Down"]]
            else:
                thissyst["emu{}_".format(year) + proc] = 0
        systs.append( (systcateg+year, "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     if proc == "ttz": thissyst["emu{}_".format(year) + proc] = "1.11"
    #     else: thissyst["emu{}_".format(year) + proc] = 0
    # systs.append( ("FlatSystTFNbTTZ{}".format(year), "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     if proc == "ttz": thissyst["emu{}_".format(year) + proc] = "1.02"
    #     else: thissyst["emu{}_".format(year) + proc] = 0
    # systs.append( ("FlatSystMTexpTTZ{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        if proc == "ttz": thissyst["emu{}_".format(year) + proc] = get_tf_uncertainty_from_txt_file("exp/ttz_emu_tf.txt")
        else: thissyst["emu{}_".format(year) + proc] = 0
    systs.append( ("FlatSystTFEMuTTZ{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        if proc == "zz": thissyst["emu{}_".format(year) + proc] = get_tf_uncertainty_from_txt_file("exp/zz_emu_tf.txt")
        else: thissyst["emu{}_".format(year) + proc] = 0
    systs.append( ("FlatSystTFEMuZZ{}".format(year), "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #t     if proc == "zz": thissyst["emu{}_".format(year) + proc] = "1.05"
    #     else: thissyst["emu{}_".format(year) + proc] = 0
    # systs.append( ("FlatSystMTexpZZ{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        if proc == "wz": thissyst["emu{}_".format(year) + proc] = "1.3" # Fake Syst
        else: thissyst["emu{}_".format(year) + proc] = 0
    systs.append( ("FlatSystWZ{}".format(year), "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     if proc == "twz": thissyst["emu{}_".format(year) + proc] = "1.47" # TWZ Syst
    #     else: thissyst["emu{}_".format(year) + proc] = 0
    # systs.append( ("FlatSystTWZ{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        thissyst["emu{}_".format(year) + proc] = "1.025"
    systs.append( ("FlatSystLumi{}".format(year), "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     thissyst["emu{}_".format(year) + proc] = "1.03"
    # systs.append( ("FlatSystsIP3D{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        thissyst["emu{}_".format(year) + proc] = "1.032"
    systs.append( ("FlatSystsTrigSF{}".format(year), "lnN", [], thissyst) )

    # Now create data card writer
    sig = hists_db["sig"]["Nominal"]
    bgs = [ hists_db[proc]["Nominal"] for proc in bkgprocs ]
    data = hists_db["data_obs"]["Nominal"]
    d = dw.DataCardWriter(sig=sig, bgs=bgs, data=data, systs=systs, no_stat_procs=["emu{}_zz".format(year), "emu{}_ttz".format(year)])

    finalyields = []
    if nbins > 1:
        for i in xrange(1, nbins+1):
            d.set_bin(i)
            d.set_region_name("bin{}".format(i))
            d.write("stats/{}/emu_datacard{}_bin{}{}.txt".format(prefix, "_"+wwzonlysuffix if wwzonlysuffix != "" else "", i, "_split" if args.nonh_vh_split else ""))
            if args.print_yields:
                vals = d.print_yields(detail=args.print_detail)
                if vals:
                    if args.nonh_vh_split:
                        print_yield_table(vals[0], vals[1], "textable/emuSplit{}{}{}".format(wwzonlysuffix, year, i), vals[2])
                        finalyields.append(vals)
                    else:
                        print_yield_table(vals[0], vals[1], "textable/emu{}{}{}".format(wwzonlysuffix, year, i), vals[2])
                        finalyields.append(vals)
    elif nbins == 1:
        d.set_bin(1)
        d.set_region_name("bin{}".format(1))
        d.write("stats/{}/emu_datacard_singlebin{}.txt".format(prefix, 1))
        if args.print_yields:
            vals = d.print_yields(detail=args.print_detail)
            if vals:
                if args.nonh_vh_split:
                    print_yield_table(vals[0], vals[1], "textable/emuSplit{}{}".format(wwzonlysuffix, year), vals[2])
                else:
                    print_yield_table(vals[0], vals[1], "textable/emu{}{}".format(wwzonlysuffix, year), vals[2])

    # colors = [2005, 2001, 2003, 2007, 920, 2012, 2011, 2002]
    # p.plot_hist(data=None, bgs=bgs, sigs=[sig], options={"bkg_sort_method":"ascending", "yaxis_range":[0.,2.5]}, colors=colors, sig_labels=["sig"], legend_labels=bkgprocs)

    ################################
    # OffZ channel data card writing
    ################################

    # number of bins
    if not args.eemm_three_bin:
        nbins = 1
        fitvar = "Yield"
        fitreg = "OffZSR"
    else:
        nbins = 3
        fitvar = "eemmSR"
        fitreg = "OffZSR"


    # Main data base to hold all the histograms
    hists_db = {}

    # Loop over the processes
    for proc in procs:

        # Retrieve the tfile
        tfile = tfiles[proc]

        # For each processes create another map to hold various histograms
        hists_db[proc] = {}

        # Loop over the systematic variations
        for syst in systnames:

            if syst == "Nominal":
                h = tfile.Get("Channel{}__{}".format(fitreg, fitvar)).Clone()
            else:
                systhacked = syst
                if proc == "NONE":
                    systhacked = ""
                h = tfile.Get("Channel{}{}__{}".format(fitreg, systhacked, fitvar)).Clone()
                # h = tfile.Get("ChannelOffZHighMET{}__Yield".format(syst)).Clone()

            h.SetTitle("offz{}_{}".format(year, proc))

            if proc == "ttz":
                before_scale = h.Integral()
                h.Scale(ttz_sf)
                after_scale = h.Integral()
                if syst == "Nominal":
                    print(year, "ttz", before_scale, after_scale)
            if proc == "zz":
                before_scale = h.Integral()
                h.Scale(zz_sf)
                after_scale = h.Integral()
                if syst == "Nominal":
                    print(year, "zz", before_scale, after_scale)

            hists_db[proc][syst] = h

    systs = []

    # ZZ CR systematic line
    onz_cr_hist = r.TH1F("onz_cr", "", nbins, 0, nbins)
    for i in range(1, nbins+1):
        onz_cr_hist.SetBinContent(i, expected_nevt_zz)
    alpha = hists_db["zz"]["Nominal"].Clone("alpha")
    alpha.Divide(onz_cr_hist)
    thissyst = {}
    for proc in mcprocs:
       if proc == "zz":
           thissyst["offz{}_".format(year) + proc] = [ "{:4f}".format(alpha.GetBinContent(i)) for i in range(1,nbins+1) ]
       else:
           thissyst["offz{}_".format(year) + proc] = 0
    systs.append(("CRZZ{}".format(year), "gmN", [onz_cr_hist], thissyst))

    # ttZ CR systematic line
    btag_cr_hist = r.TH1F("btag_cr", "", nbins, 0, nbins)
    for i in range(1, nbins+1):
        btag_cr_hist.SetBinContent(i, expected_nevt_ttz)
    alpha = hists_db["ttz"]["Nominal"].Clone("alpha")
    alpha.Divide(btag_cr_hist)
    thissyst = {}
    for proc in mcprocs:
       if proc == "ttz":
           thissyst["offz{}_".format(year) + proc] = [ "{:4f}".format(alpha.GetBinContent(i)) for i in range(1,nbins+1) ]
       else:
           thissyst["offz{}_".format(year) + proc] = 0
    systs.append(("CRTTZ{}".format(year), "gmN", [btag_cr_hist], thissyst))

    # Experimental systematics
    for systcateg in systcategs:
        thissyst = {}
        for proc in mcprocs:
            # if proc not in ["zz", "ttz", "wz", "twz"]:
            if proc not in ["zz", "ttz", "wz"]:
                thissyst["offz{}_".format(year) + proc] = [hists_db[proc][systcateg+"Up"], hists_db[proc][systcateg+"Down"]]
            else:
                thissyst["offz{}_".format(year) + proc] = 0
        systs.append( (systcateg+year, "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     if proc == "ttz": thissyst["offz{}_".format(year) + proc] = "1.10"
    #     else: thissyst["offz{}_".format(year) + proc] = 0
    # systs.append( ("FlatSystTFeemmTTZ{}".format(year), "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     if proc == "ttz": thissyst["offz{}_".format(year) + proc] = "1.03"
    #     else: thissyst["offz{}_".format(year) + proc] = 0
    # systs.append( ("FlatSystMETexpTTZ{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        if proc == "ttz": thissyst["offz{}_".format(year) + proc] = get_tf_uncertainty_from_txt_file("exp/ttz_eemm_tf.txt")
        else: thissyst["offz{}_".format(year) + proc] = 0
    systs.append( ("FlatSystTFEEMMTTZ{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        if proc == "zz": thissyst["offz{}_".format(year) + proc] = get_tf_uncertainty_from_txt_file("exp/zz_eemm_tf.txt")
        else: thissyst["offz{}_".format(year) + proc] = 0
    systs.append( ("FlatSystTFEEMMZZ{}".format(year), "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     if proc == "zz": thissyst["offz{}_".format(year) + proc] = "1.23"
    #     else: thissyst["offz{}_".format(year) + proc] = 0
    # systs.append( ("FlatSystMETexpZZ{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        if proc == "wz": thissyst["offz{}_".format(year) + proc] = "1.3" # Fake Syst
        else: thissyst["offz{}_".format(year) + proc] = 0
    systs.append( ("FlatSystWZ{}".format(year), "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     if proc == "twz": thissyst["offz{}_".format(year) + proc] = "1.47" # TWZ Syst
    #     else: thissyst["offz{}_".format(year) + proc] = 0
    # systs.append( ("FlatSystTWZ{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        thissyst["offz{}_".format(year) + proc] = "1.025"
    systs.append( ("FlatSystLumi{}".format(year), "lnN", [], thissyst) )

    # # Flat additional systematics
    # thissyst = {}
    # for proc in mcprocs:
    #     thissyst["offz{}_".format(year) + proc] = "1.03"
    # systs.append( ("FlatSystsIP3D{}".format(year), "lnN", [], thissyst) )

    # Flat additional systematics
    thissyst = {}
    for proc in mcprocs:
        thissyst["offz{}_".format(year) + proc] = "1.032"
    systs.append( ("FlatSystsTrigSF{}".format(year), "lnN", [], thissyst) )

    # Now create data card writer
    sig = hists_db["sig"]["Nominal"]
    bgs = [ hists_db[proc]["Nominal"] for proc in bkgprocs ]
    data = hists_db["data_obs"]["Nominal"]
    d = dw.DataCardWriter(sig=sig, bgs=bgs, data=data, systs=systs, no_stat_procs=["offz{}_zz".format(year), "offz{}_ttz".format(year)])

    if nbins > 1:
        for i in xrange(1, nbins+1):
            d.set_bin(i)
            d.set_region_name("bin{}".format(i))
            d.write("stats/{}/offz_datacard{}_bin{}{}.txt".format(prefix, "_"+wwzonlysuffix if wwzonlysuffix != "" else "", i, "_split" if args.nonh_vh_split else ""))
            if args.print_yields:
                vals = d.print_yields(detail=args.print_detail)
                if vals:
                    if args.nonh_vh_split:
                        print_yield_table(vals[0], vals[1], "textable/offzSplit{}{}{}".format(wwzonlysuffix, year, i), vals[2])
                        finalyields.append(vals)
                    else:
                        print_yield_table(vals[0], vals[1], "textable/offz{}{}{}".format(wwzonlysuffix, year, i), vals[2])
                        finalyields.append(vals)
    elif nbins == 1:
        d.set_bin(1)
        d.set_region_name("bin{}".format(1))
        d.write("stats/{}/offz_datacard_singlebin{}.txt".format(prefix, 1))
        if args.print_yields:
            vals = d.print_yields(detail=args.print_detail)
            if vals:
                if args.nonh_vh_split:
                    print_yield_table(vals[0], vals[1], "textable/offzSplit{}{}".format(wwzonlysuffix, year), vals[2])
                else:
                    print_yield_table(vals[0], vals[1], "textable/offz{}{}".format(wwzonlysuffix, year), vals[2])


    if len(finalyields) > 5:

        procs = ["sig", "wzz", "zzz", "zz", "ttz", "twz", "wz", "higgs", "other"]
        if args.nonh_vh_split:
            procs = ["sig", "zhwwz", "nonhwzz", "whwzz", "nonhzzz", "zhzzz", "zz", "ttz", "twz", "wz", "higgs", "other"]

        histsdict = {}
        for proc in procs:
            # if proc == "sig":
            #     h = r.TH1F("WWZ", "", 7, 0, 7)
            # elif proc == "wzz":
            #     h = r.TH1F("WZZ", "", 7, 0, 7)
            # elif proc == "zzz":
            #     h = r.TH1F("ZZZ", "", 7, 0, 7)
            # else:
            #     h = r.TH1F("Fit{}".format(proc), "", 7, 0, 7)
            h = r.TH1F("Fit{}".format(proc), "", 7, 0, 7)
            h.GetXaxis().SetBinLabel(1, "e#mu 1")
            h.GetXaxis().SetBinLabel(2, "e#mu 2")
            h.GetXaxis().SetBinLabel(3, "e#mu 3")
            h.GetXaxis().SetBinLabel(4, "e#mu 4")
            h.GetXaxis().SetBinLabel(5, "ee/#mu#mu A")
            h.GetXaxis().SetBinLabel(6, "ee/#mu#mu B")
            h.GetXaxis().SetBinLabel(7, "ee/#mu#mu C")
            # h.GetXaxis().SetBinLabel(8, "ee/#mu#mu Bin C")
            histsdict[proc] = h
        
        for index, item in enumerate(finalyields):
            for procfullname, rate in zip(item[0], item[1]):
                procname = procfullname.split("_")[1]
                print(index, procname, rate)
                histsdict[procname].SetBinContent(index+1, rate.val)
                histsdict[procname].SetBinError(index+1, rate.err)
        if args.nonh_vh_split:
            bkghists = [ histsdict[proc].Clone() for proc in procs[6:] ]
            sighists = [ histsdict[proc].Clone() for proc in procs[:6] ]
        else:
            bkghists = [ histsdict[proc].Clone() for proc in procs[3:] ]
            sighists = [ histsdict[proc].Clone() for proc in procs[:3] ]

        lumi = 137
        if "2016" in year: lumi = 35.9
        if "2017" in year: lumi = 41.3
        if "2018" in year: lumi = 59.74

        p.plot_hist(bgs=bkghists,
                sigs=sighists,
                options={
                "output_name": "fitplot/fit{}{}.pdf".format(wwzonlysuffix, year),
                "print_yield":True,
                "signal_scale": 1,
                "legend_scalex":1.8,
                "legend_scaley":1.0,
                "legend_ncolumns": 3,
                "legend_smart": True,
                "yaxis_log":False,
                "ymax_scale": 1.2,
                "lumi_value":lumi,
                # "no_overflow": True,
                "remove_underflow": True,
                "xaxis_ndivisions":505,
                "ratio_range":[0.,2.],
                "xaxis_label":"Fit regions",
                "ratio_xaxis_title":"Fit regions",
                "no_ratio": True,
                },
                colors = [2001, 2005, 2007, 2003, 2011, 920, 2012, 2011, 2002],
                legend_labels = ["ZZ", "t#bar{t}Z", "tWZ", "WZ", "Higgs", "Other"],
                # sig_labels = ["WWZ","WZZ","ZZZ"]
                )

def rebin36(h):

    # To rebin emu channels 180 bins to 5 bin

    h.Rebin(36)
    bin5 = E(h.GetBinContent(5), h.GetBinError(5))
    bin6 = E(h.GetBinContent(6), h.GetBinError(6))
    bin5plus6 = bin5 + bin6
    h.SetBinContent(6, 0)
    h.SetBinError(6, 0)
    h.SetBinContent(5, bin5plus6.val)
    h.SetBinError(5, bin5plus6.err)

    return h

def print_yield_table(procs, rates, output_name, datayield):

    hists = []
    bkgh = r.TH1F("Total", "Total", 1, 0, 1)
    total_rate = E(0, 0)
    for proc, rate in zip(procs, rates):
        procname = proc.split("_")[1]
        h = r.TH1F(procname, procname, 1, 0, 1)
        if args.nonh_vh_split:
            if procname != "sig" and procname != "whwzz" and procname != "nonhwzz" and procname != "nonhzzz" and procname != "zhzzz" and procname != "zhwwz" and procname != "nonhwwz":
                total_rate += rate
        else:
            if procname != "sig" and procname != "wzz" and procname != "wzz" and procname != "zzz":
                total_rate += rate
        h.SetBinContent(1, rate.val)
        h.SetBinError(1, rate.err)
        hists.append(h)
    bkgh.SetBinContent(1, total_rate.val)
    bkgh.SetBinError(1, total_rate.err)
    hists.insert(0, bkgh)
    obsh = bkgh.Clone("obs")
    obsh.SetBinContent(1, datayield)
    hists.insert(0, obsh)

    p.print_yield_table_from_list(hists, output_name + ".txt", prec=2, binrange=[1])
    p.print_yield_tex_table_from_list(hists, output_name + ".tex", prec=2)

def get_tf_uncertainty_from_txt_file(txtpath):

    # e.g. get_tf_uncertainty_from_txt_file("exp/ttz_emu_tf.txt")

    # +-------+---------+----------+----------+---------+----------+----------+---------+---------+---------+---------+---------+---------+---------+---------+
    # | Bin#  |  Ratio  |  Yield   |  Total   |  Stat   | ElLepSF  | MuLepSF  |   JES   | Pileup  | BTagHF  | BTagLF  |   MET   |   PDF   |   Qsq   | AlphaS  |
    # |-------+---------+----------+----------+---------+----------+----------+---------+---------+---------+---------+---------+---------+---------+---------+
    # | Bin1  | 0.0680  |  0.0000  | 11.1660  | 4.4533  |  1.1700  |  0.0076  | 3.8027  | 3.7372  | 8.4400  | 1.1282  | 0.0000  | 1.0664  | 1.1861  | 0.0284  |
    # | Bin2  | 0.0000  |  3.2383  | 12.3090  | 4.3122  |  6.5930  |  1.0479  | 3.5228  | 3.2690  | 7.6211  | 1.0312  | 0.0000  | 0.6620  | 1.3699  | 1.9480  |
    # | Bin3  | 0.0000  | 47.5963  |  9.5567  | 1.1119  |  5.4230  |  1.0403  | 0.2795  | 7.2161  | 0.8178  | 0.0970  | 0.0000  | 1.7295  | 0.0821  | 1.9431  |
    # +-------+---------+----------+----------+---------+----------+----------+---------+---------+---------+---------+---------+---------+---------+---------+

    f = open(txtpath)
    lines = f.readlines()

    result = ""

    for line in lines:
        if "Bin1" not in line:
            continue
        result = line.split()[7]

    print("from TF file ", txtpath, "got TF uncertainty of", float(result))

    return str(1 + float(result) / 100.)

if __name__ == "__main__":

    main()

