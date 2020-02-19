#include "ElectronScaleFactors.h"

#include <stdexcept>

namespace {
  struct IDInfo {
    LeptonID id;       // the enum value in this class
    std::string flag;  // the flag in the TnP tool which appears in the filename of the histogram file
  };
}  // namespace

ElectronScaleFactors::ElectronScaleFactors(std::string const& electronScaleFactorsPath) {
  const std::vector<IDInfo> idInfos{
      {LeptonID::commonVeto, "passingVVVCommonVetoID"},
      {LeptonID::zCandidate, "passingWVZCandidateZID"},
      {LeptonID::wCandidate, "passingWVZCandidateWID"},
      {LeptonID::threeLepton, "passingWWW3LID"},
      {LeptonID::sameSign, "passingWWWSameSignID"},
  };

  std::string idPath = electronScaleFactorsPath + "/identification/";

  for (int year : {2016, 2017, 2018}) {
    for (auto const& idInfo : idInfos) {
      histMaps_.emplace_back(idPath + std::to_string(year) + "_" + idInfo.flag + ".root:EGamma_SF2D");
    }
  }
}

float ElectronScaleFactors::operator()(int year, LeptonID id, float eta, float pt, int vare) {
  // basic eta and pt checks to make sure they are not accidentally swapped
  if (pt < 0) {
    throw std::runtime_error("The pt passed to ElectronScaleFactors is negative. Did you swap pt and eta?");
  }
  if (pt < 5.) {
    throw std::runtime_error(
        "The pt passed to ElectronScaleFactors is smaller than 5 GeV. This is ridiculous. Did you swap pt and eta?");
  }
  if (vare == 0) {
    return histMaps_[getHistMapIndex(year, id)].eval(eta, pt);
  }
  if (vare == 1) {
    return histMaps_[getHistMapIndex(year, id)].eval_up(eta, pt);
  }
  if (vare == -1) {
    return histMaps_[getHistMapIndex(year, id)].eval_down(eta, pt);
  }
}

// Let's check at compile time that the indices in the histMap array turn out fine, just in case
static_assert(ElectronScaleFactors::getHistMapIndex(2016, LeptonID::commonVeto) == 0);
static_assert(ElectronScaleFactors::getHistMapIndex(2016, LeptonID::zCandidate) == 1);
static_assert(ElectronScaleFactors::getHistMapIndex(2016, LeptonID::wCandidate) == 2);
static_assert(ElectronScaleFactors::getHistMapIndex(2016, LeptonID::threeLepton) == 3);
static_assert(ElectronScaleFactors::getHistMapIndex(2016, LeptonID::sameSign) == 4);
static_assert(ElectronScaleFactors::getHistMapIndex(2017, LeptonID::commonVeto) == 5);
static_assert(ElectronScaleFactors::getHistMapIndex(2017, LeptonID::zCandidate) == 6);
static_assert(ElectronScaleFactors::getHistMapIndex(2017, LeptonID::wCandidate) == 7);
static_assert(ElectronScaleFactors::getHistMapIndex(2017, LeptonID::threeLepton) == 8);
static_assert(ElectronScaleFactors::getHistMapIndex(2017, LeptonID::sameSign) == 9);
static_assert(ElectronScaleFactors::getHistMapIndex(2018, LeptonID::commonVeto) == 10);
static_assert(ElectronScaleFactors::getHistMapIndex(2018, LeptonID::zCandidate) == 11);
static_assert(ElectronScaleFactors::getHistMapIndex(2018, LeptonID::wCandidate) == 12);
static_assert(ElectronScaleFactors::getHistMapIndex(2018, LeptonID::threeLepton) == 13);
static_assert(ElectronScaleFactors::getHistMapIndex(2018, LeptonID::sameSign) == 14);
