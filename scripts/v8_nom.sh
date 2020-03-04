
# # sh runall.sh -v v0.1.21 -t WVZMVA -b v8_nom -1 -2
# sh runall.sh -v v0.1.21 -t WVZMVA -b v8_nom -a

# Z peak plots
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__lepZPt0 -u -n 45 -x \"#it{p}_{T}^{Z-tag1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__lepZPt1 -u -n 45 -x \"#it{p}_{T}^{Z-tag2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__lepNPt0 -u -n 45 -x \"#it{p}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__lepNPt1 -u -n 45 -x \"#it{p}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__MllZCandZoom -u -n 45 -x \"#it{m}_{ll}^{Z-tag} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__MllNomZoom -u -n 45 -x \"#it{m}_{ll} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__MET -u -n 45 -x \"MET [GeV]\" -y -r 0.01,1e6
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__Njet -u -x \"n_{j}\" -y
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d zpeak -p Cut4LepLeptonPt__Nbjet -u -x \"n_{b}\" -y

# Signal discriminants
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d signal_discriminant -p ChannelEMu__MT2 -n 30 -x \"#it{m}_{T2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d signal_discriminant -p ChannelEMu__MllNom -n 30 -x \"#it{m}_{ll} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d signal_discriminant -p ChannelOffZ__MET -n 30 -y -x \"MET [GeV]\" -r 0.01,1e6
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d signal_discriminant -p ChannelOffZ__Pt4l -n 30 -y -x \"#it{p}_{T,4l} [GeV]\" -r 0.01,1e6

# MC only table yields
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuHighMT__emuSR -n 4 -x \"e#mu SR bins\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelOffZSR__eemmSR -n 3 -x \"ee/#mu#mu SR bins\" -a

# BDT discriminant output
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuBDT__emuBDTTTZScore -n 30 -x \"e#mu channel t#bar{t}Z BDT score\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuBDT__emuBDTZZScore -n 30 -x \"e#mu channel ZZ BDT score\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelOffZBDTPre__offzBDTScore -n 30 -x \"ee/#mu#mu channel ZZ BDT score\" -a -S 10
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuBDT__emuBDT_Nominal -n 30 -x \"e#mu channel BDT bins\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelOffZBDT__offzBDT_Nominal -n 30 -x \"ee/#mu#mu channel BDT bins\" -a -y
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuBDT0__Yield -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuBDT1__Yield -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuBDT2__Yield -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuBDT3__Yield -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelEMuBDT4__Yield -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelOffZBDTA__Yield -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelOffZBDTB__Yield -a

# ZZ control region table
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelOnZ__Yield -n 1 -x \"Yield\" -a -u -i

# TTZ control region table
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d fits -p ChannelBTagEMu__Yield -n 1 -x \"Yield\" -a -u -i

# TWZ validation region
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d twzcr -p ChannelBTagEMuOneNb__Njet  -u -x \"#it{n}_{jet}\"

# HZZ validation region
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d hzzcr -p ChannelHZZ4l__MZZ4lZoom -u -x \"#it{m}_{4l} [GeV]\" -n 15
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d hzzcr -p ChannelHZZ4l__MllZCand -u -x \"#it{m}_{ll,1st} [GeV]\" -n 45
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d hzzcr -p ChannelHZZ4l__MllZ2Cand -u -x \"#it{m}_{ll,2nd} [GeV]\" -n 45

# On Z MET distribution modeling check
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d met_correction -u -p ChannelOnZ__MET -n 45 -x \"MET [GeV]\" -y -r 0.01,1e7
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d met_correction -u -p ChannelOnZ__OrigMET -n 45 -x \"MET before smearing [GeV]\" -y -r 0.01,1e7

# Btag CR plots
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d btagcr -p ChannelBTagEMu__MllNom -n 15 -u -x \"#it{m}_{ll} [GeV]\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d btagcr -p ChannelBTagEMu__MT2 -n 15 -u -x \"#it{m}_{T2} [GeV]\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d btagcr -p ChannelBTagEMu__MET -n 15 -u -x \"MET [GeV]\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d btagcr -p ChannelBTagEMu__Pt4l -n 15 -u -x \"#it{p}_{T,4l} [GeV]\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d btagcr -p ChannelBTagEMu__emuBDTTTZScore -n 15 -u -x \"e#mu channel t#bar{t}Z BDT score\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d btagcr -p ChannelBTagEMu__emuBDTZZScore -n 15 -u -x \"e#mu channel ZZ BDT score\" -a
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d btagcr -p ChannelBTagEMu__offzBDTScore -n 15 -u -x \"ee/#mu#mu channel ZZ BDT score\" -a

# On Z CR plots
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZ__MT2 -n 30 -u -x \"#it{m}_{T2} [GeV]\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZ__MET -n 45 -u -x \"MET [GeV]\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZ__Pt4l -n 45 -u -x \"#it{p}_{T,4l} [GeV]\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZ__emuBDT_Nominal  -u -x \"e#mu channel BDT bins\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZ__offzBDTScore -n 30 -x \"ee/#mu#mu channel ZZ BDT score\" -a -u
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZBDT0__Yield  -u -x \"e#mu channel BDT bin 0\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZBDT1__Yield  -u -x \"e#mu channel BDT bin 1\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZBDT2__Yield  -u -x \"e#mu channel BDT bin 2\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZBDT3__Yield  -u -x \"e#mu channel BDT bin 3\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOnZBDT4__Yield  -u -x \"e#mu channel BDT bin 4\" -a 
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOffZBDTCR__offzBDTScore -n 30 -x \"ee/#mu#mu channel ZZ BDT score\" -a -u
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d onzcr -p ChannelOffZBDTCR__Yield -n 30 -x \"ee/#mu#mu channel ZZ BDT score\" -a -u

# WZ AR plots
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d wzar -p ChannelAREMuHighMT__lepFrelIso03EA  -u -x \"I_{rel,0.3,EA,e}\" -n 5
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d wzar -p ChannelAREMuHighMT__lepFrelIso04DB  -u -x \"I_{rel,0.4,DB,#mu}\" -n 5

# BDT input plots
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_minDRJetToLep3 -u -x \"min#DeltaR(jet,l1)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_minDRJetToLep4 -u -x \"min#DeltaR(jet,l2)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_m_4l -u -x \"#it{m}_{4l}\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_vecsum_pt_4l -u -x \"#it{p}_{T}^{4l} vector sum [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_scalarsum_pt_4l -u -x \"#Sigma_{i}p_{T}^{i} scalar sum of four-leptons [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_jet1Pt -u -x \"#it{p}_{T}^{j1} [GeV]\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_lep3Pt -u -x \"#it{p}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_lep4Pt -u -x \"#it{p}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_met_pt -u -x \"MET [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_mt2 -u -x \"#it{m}_{T2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_lep3MT -u -x \"#it{m}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_lep4MT -u -x \"#it{m}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_MllN -u -x \"#it{m}_{ll} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_pt_zeta -u -x \"P_{#zeta} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_pt_zeta_vis -u -x \"P_{#zeta}^{vis} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOnZ__BDTInput_ZPt -u -x \"#it{p}_{T,ll}^{Z-tag} [GeV]\"

sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_minDRJetToLep3 -u -x \"min#DeltaR(jet,l1)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_minDRJetToLep4 -u -x \"min#DeltaR(jet,l2)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_m_4l -u -x \"#it{m}_{4l}\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_vecsum_pt_4l -u -x \"#it{p}_{T}^{4l} vector sum [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_scalarsum_pt_4l -u -x \"#Sigma_{i}p_{T}^{i} scalar sum of four-leptons [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_jet1Pt -u -x \"#it{p}_{T}^{j1} [GeV]\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_lep3Pt -u -x \"#it{p}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_lep4Pt -u -x \"#it{p}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_met_pt -u -x \"MET [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_mt2 -u -x \"#it{m}_{T2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_lep3MT -u -x \"#it{m}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_lep4MT -u -x \"#it{m}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_MllN -u -x \"#it{m}_{ll} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_pt_zeta -u -x \"P_{#zeta} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_pt_zeta_vis -u -x \"P_{#zeta}^{vis} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 45 -p ChannelOffZBDTCR__BDTInput_ZPt -u -x \"#it{p}_{T,ll}^{Z-tag} [GeV]\"

sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_minDRJetToLep3 -u -x \"min#DeltaR(jet,l1)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_minDRJetToLep4 -u -x \"min#DeltaR(jet,l2)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_m_4lLarge -u -x \"#it{m}_{4l}\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_vecsum_pt_4l -u -x \"#it{p}_{T}^{4l} vector sum [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_scalarsum_pt_4l -u -x \"#Sigma_{i}p_{T}^{i} scalar sum of four-leptons [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_jet1Pt -u -x \"#it{p}_{T}^{j1} [GeV]\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_lep3Pt -u -x \"#it{p}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_lep4Pt -u -x \"#it{p}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_met_pt -u -x \"MET [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_mt2 -u -x \"#it{m}_{T2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_lep3MT -u -x \"#it{m}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_lep4MT -u -x \"#it{m}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_MllN -u -x \"#it{m}_{ll} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_pt_zeta -u -x \"P_{#zeta} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_pt_zeta_vis -u -x \"P_{#zeta}^{vis} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 20 -p ChannelBTagEMu__BDTInput_ZPt -u -x \"#it{p}_{T,ll}^{Z-tag} [GeV]\"

# BDT input plots in SR
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_minDRJetToLep3 -x \"min#DeltaR(jet,l1)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_minDRJetToLep4 -x \"min#DeltaR(jet,l2)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_m_4l -x \"#it{m}_{4l}\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_vecsum_pt_4l -x \"#it{p}_{T}^{4l} vector sum [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_scalarsum_pt_4l -x \"#Sigma_{i}p_{T}^{i} scalar sum of four-leptons [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_jet1Pt -x \"#it{p}_{T}^{j1} [GeV]\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_lep3Pt -x \"#it{p}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_lep4Pt -x \"#it{p}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_met_pt -x \"MET [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_mt2 -x \"#it{m}_{T2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_lep3MT -x \"#it{m}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_lep4MT -x \"#it{m}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_MllN -x \"#it{m}_{ll} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_pt_zeta -x \"P_{#zeta} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_pt_zeta_vis -x \"P_{#zeta}^{vis} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelEMuBDT__BDTInput_ZPt -x \"#it{p}_{T,ll}^{Z-tag} [GeV]\"

sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_minDRJetToLep3 -x \"min#DeltaR(jet,l1)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_minDRJetToLep4 -x \"min#DeltaR(jet,l2)\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_m_4l -x \"#it{m}_{4l}\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_vecsum_pt_4l -x \"#it{p}_{T}^{4l} vector sum [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_scalarsum_pt_4l -x \"#Sigma_{i}p_{T}^{i} scalar sum of four-leptons [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_jet1Pt -x \"#it{p}_{T}^{j1} [GeV]\" -U
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_lep3Pt -x \"#it{p}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_lep4Pt -x \"#it{p}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_met_pt -x \"MET [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_mt2 -x \"#it{m}_{T2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_lep3MT -x \"#it{m}_{T}^{l1} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_lep4MT -x \"#it{m}_{T}^{l2} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_MllN -x \"#it{m}_{ll} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_pt_zeta -x \"P_{#zeta} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_pt_zeta_vis -x \"P_{#zeta}^{vis} [GeV]\"
sh runall.sh -b v8_nom -t WVZMVA -v v0.1.21 -d bdtinput -n 30 -p ChannelOffZBDT__BDTInput_ZPt -x \"#it{p}_{T,ll}^{Z-tag} [GeV]\"

