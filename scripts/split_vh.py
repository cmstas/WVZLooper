#!/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import ROOT as r
import sys
import os

user_data_path = os.getenv("WVZ_DATA_PATH")

try:
    version = sys.argv[1]
except:
    print("Usage:")
    print("")
    print("   {} VERSION".format(sys.argv[0]))
    print("")
    sys.exit(-1)

f = r.TFile(user_data_path + "/babies/{}/vh_nonbb_amcatnlo_1.root".format(version))

t = f.Get("t")
h = f.Get("h_neventsinfile")

f_www = r.TFile(user_data_path + "/babies/{}/wh_ww_amcatnlo_1.root".format(version), "recreate")
t_www = t.CopyTree("VHchannel==24&&Higgschannel==24")
t_www.Write()
h.Write()

f_zww = r.TFile(user_data_path + "/babies/{}/zh_ww_amcatnlo_1.root".format(version), "recreate")
t_zww = t.CopyTree("VHchannel==23&&Higgschannel==24")
t_zww.Write()
h.Write()

f_wzz = r.TFile(user_data_path + "/babies/{}/wh_zz_amcatnlo_1.root".format(version), "recreate")
t_wzz = t.CopyTree("VHchannel==24&&Higgschannel==23")
t_wzz.Write()
h.Write()

f_zzz = r.TFile(user_data_path + "/babies/{}/zh_zz_amcatnlo_1.root".format(version), "recreate")
t_zzz = t.CopyTree("VHchannel==23&&Higgschannel==23")
t_zzz.Write()
h.Write()

f_rest = r.TFile(user_data_path + "/babies/{}/vh_nonbbwwzz_amcatnlo_1.root".format(version), "recreate")
t_rest = t.CopyTree("Higgschannel!=23&&Higgschannel!=24")
t_rest.Write()
h.Write()

