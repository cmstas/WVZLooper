#!/bin/env python

import ROOT as r
import plottery_wrapper as p

dirpath = "outputs/WVZMVA2016_v0.1.21_WVZMVA2017_v0.1.21_WVZMVA2018_v0.1.21/y2016_syst_nosmear_20191106_y2017_syst_nosmear_20191106_y2018_syst_nosmear_20191106/"

fnames = [
"data.root",
"dy.root",
"dyttbar.root",
"higgs.root",
"nonh_www.root",
"nonh_wwz.root",
"nonh_wzz.root",
"nonh_zzz.root",
"other.root",
"othernoh.root",
"othervvv.root",
"rare.root",
"rarevvv.root",
"sig.root",
"triother.root",
"ttbar.root",
"ttz.root",
"twz.root",
"wh_www.root",
"wh_wzz.root",
"www.root",
"wwz.root",
"wz.root",
"wzz.root",
"zh_wwz.root",
"zh_zzz.root",
"zz.root",
"zzz.root",
]

histname = "ChannelBTagEMuOneNb__Njet"

tfiles = {}
thists = {}

for fname in fnames:

    tfiles[fname] = r.TFile(dirpath + fname)
    thists[fname] = tfiles[fname].Get(histname)

hother = thists["othernoh.root"].Clone("other")
hother.Add(thists["higgs.root"])
hother.Add(thists["wz.root"])
hother.Add(thists["sig.root"])

htwz = thists["twz.root"].Clone("twz")
httz = thists["ttz.root"].Clone("ttz")
hzz = thists["zz.root"].Clone("zz")

hdata = thists["data.root"].Clone("data")
bgs = [httz, hzz, hother, htwz] 

p.plot_hist(bgs=[httz.Clone(), hzz.Clone(), hother.Clone(), htwz.Clone()], data=hdata, options={"output_name":"fit_onenb.pdf", "ratio_range":[0.,2.]})

njet = r.RooRealVar("njet", "njet", 0., 6.)

datahists = {}
pdfhists = {}
normhists = {}

for bg in bgs + [hdata]:
    name = bg.GetName()
    datahists[name] = r.RooDataHist(name, name, r.RooArgList(njet), bg)
    pdfhists[name] = r.RooHistPdf(name+"_pdf", name+"_pdf", r.RooArgSet(njet), datahists[name])
    if name not in ["twz", "ttz"]:
        normhists[name] = r.RooRealVar(name+"_n", name+"_n", bg.Integral(), bg.Integral() * 0.999, bg.Integral() * 1.001)
    else:
        normhists[name] = r.RooRealVar(name+"_n", name+"_n", bg.Integral(), bg.Integral() * 0.1, bg.Integral() * 10.0)

model = r.RooAddPdf("model", "model", r.RooArgList(pdfhists["ttz"], pdfhists["zz"], pdfhists["other"], pdfhists["twz"]), r.RooArgList(normhists["ttz"], normhists["zz"], normhists["other"], normhists["twz"]))
fitres = model.fitTo(datahists["data"], r.RooFit.SumW2Error(r.kFALSE), r.RooFit.Extended(), r.RooFit.Save(r.kTRUE), r.RooFit.Range(0., 6.))

c1 = r.TCanvas()

mesframe = njet.frame()
datahists["data"].plotOn(mesframe)
model.plotOn(mesframe)
model.plotOn(mesframe, r.RooFit.Components("twz_pdf"), r.RooFit.LineStyle(r.kDashed))
# model.plotOn(mesframe, r.RooFit.Components("ttz_pdf"), r.RooFit.LineColor(r.kRed), r.RooFit.LineStyle(r.kDashed))
mesframe.Draw()

c1.SaveAs("test.pdf")

for bg in bgs + [hdata]:

    print bg.GetName(), normhists[bg.GetName()].getValV() / bg.Integral(), normhists[bg.GetName()].getError() / bg.Integral()

httz.Scale(normhists["ttz"].getValV() / httz.Integral())
htwz.Scale(normhists["twz"].getValV() / htwz.Integral())

p.plot_hist(bgs=[httz, hzz, hother, htwz], data=hdata, options={"output_name":"fit_onenb_postfit.pdf", "ratio_range":[0.,2.]})
