
# sh runall.sh -v v0.1.21 -t WVZMVA -b v2 -s -1 -2
# sh runall.sh -v v0.1.21 -t WVZMVA -b v2 -a

# extrapolation table
# python scripts/extrapolation_uncertainty.py

# systematic table
mkdir -p systs
python scripts/write_datacards.py -b v8 -t WVZMVA -v v0.1.21  -y -1 -d -w > systs/syst.log
python scripts/write_datacards.py -b v8 -t WVZMVA -v v0.1.21  -y -1 -d > systs/syst_onesig.log
python scripts/write_datacards_bdt.py -b v8 -t WVZMVA -v v0.1.21  -y -1 -d -w > systs/syst_bdt.log
python scripts/write_datacards_bdt.py -b v8 -t WVZMVA -v v0.1.21  -y -1 -d > systs/syst_bdt_onesig.log

# Yield table with estimates
python scripts/write_datacards.py -b v8 -t WVZMVA -v v0.1.21  -y -w -3 -4
python scripts/write_datacards.py -b v8 -t WVZMVA -v v0.1.21  -y -w -1
python scripts/write_datacards_bdt.py -b v8 -t WVZMVA -v v0.1.21  -y -w
python scripts/write_datacards_bdt.py -b v8 -t WVZMVA -v v0.1.21  -y -w -1
python scripts/write_datacards.py -b v8 -t WVZMVA -v v0.1.21  -y -3 -4
python scripts/write_datacards_bdt.py -b v8 -t WVZMVA -v v0.1.21  -y
python scripts/write_datacards.py -b v8 -t WVZMVA -v v0.1.21  -y -w -3 -4 -s
python scripts/write_datacards.py -b v8 -t WVZMVA -v v0.1.21  -y -w -1 -s
python scripts/write_datacards_bdt.py -b v8 -t WVZMVA -v v0.1.21  -y -w -s
python scripts/write_datacards_bdt.py -b v8 -t WVZMVA -v v0.1.21  -y -w -1 -s

# 5/6L
python scripts/write_5l_datacards.py -b v8 -t WVZMVA -v v0.1.21 -w
python scripts/write_6l_datacards.py -b v8 -t WVZMVA -v v0.1.21 -w
