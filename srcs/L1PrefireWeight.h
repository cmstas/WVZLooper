#ifndef L1PrefireWeight_h
#define L1PrefireWeight_h

#include "rooutil.h"
#include "wvztree.h"

class L1PrefireWeight
{
    public:

        RooUtil::HistMap* histmap_l1wgt_2016;
        RooUtil::HistMap* histmap_l1wgt_2017;

        L1PrefireWeight()
        {
            histmap_l1wgt_2016 = new RooUtil::HistMap("scale_factors/l1prefire//L1prefiring_jetpt_2016BtoH.root:L1prefiring_jetpt_2016BtoH");
            histmap_l1wgt_2017 = new RooUtil::HistMap("scale_factors/l1prefire//L1prefiring_jetpt_2017BtoF.root:L1prefiring_jetpt_2017BtoF");
        }

        ~L1PrefireWeight()
        {
            delete histmap_l1wgt_2016;
            delete histmap_l1wgt_2017;
        }

        float l1wgt(int year, bool isData)
        {
            if (isData)
                return 1;
            if (year == 2018)
                return 1;
            float weight = 1;
            for (unsigned int i = 0; i < wvz.jets_p4().size(); ++i)
            {
                if (fabs(wvz.jets_p4()[i].Eta()) <= 2.0) continue;
                if (fabs(wvz.jets_p4()[i].Eta()) >= 3.0) continue;
                if (year == 2016)
                    weight *= (1. - (histmap_l1wgt_2016->eval(wvz.jets_p4()[i].Eta(), std::max(std::min(float(499.), wvz.jets_p4()[i].Pt()), float(16.)))));
                else
                    weight *= (1. - (histmap_l1wgt_2017->eval(wvz.jets_p4()[i].Eta(), std::max(std::min(float(499.), wvz.jets_p4()[i].Pt()), float(16.)))));
            }
            return weight;
        }
};

#endif
