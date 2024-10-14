#include <string>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

class MtdTracksHarvester : public DQMEDHarvester {
public:
  explicit MtdTracksHarvester(const edm::ParameterSet& iConfig);
  ~MtdTracksHarvester() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

private:
  void computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result);
  void normalize(MonitorElement* h, double scale);

  const std::string folder_;

  // --- Histograms
  MonitorElement* meBtlEtaEff_;
  MonitorElement* meBtlPhiEff_;
  MonitorElement* meBtlPtEff_;
  MonitorElement* meEtlEtaEff_[2];
  MonitorElement* meEtlEtaEffLowPt_[2];
  MonitorElement* meEtlPhiEff_[2];
  MonitorElement* meEtlPtEff_[2];
  MonitorElement* meEtlEtaEff2_[2];
  MonitorElement* meEtlEtaEff2LowPt_[2];
  MonitorElement* meEtlPhiEff2_[2];
  MonitorElement* meEtlPtEff2_[2];
  MonitorElement* meTPPtSelEff_;
  MonitorElement* meTPEtaSelEff_;
  MonitorElement* meTPPtMatchEff_;
  MonitorElement* meTPEtaMatchEff_;
  MonitorElement* meTPPtMatchEtl2Eff_;
  MonitorElement* meTPEtaMatchEtl2Eff_;
  MonitorElement* meTPmtdPtSelEff_;
  MonitorElement* meTPmtdEtaSelEff_;
  MonitorElement* meTPmtdPtMatchEff_;
  MonitorElement* meTPmtdEtaMatchEff_;
  MonitorElement* meNoTimeFraction_;
  MonitorElement* meExtraPtEff_;
  MonitorElement* meExtraPtEtl2Eff_;
  MonitorElement* meExtraEtaEff_;
  MonitorElement* meExtraEtaEtl2Eff_;
  MonitorElement* meExtraPhiAtBTLEff_;
  MonitorElement* meExtraMTDfailExtenderEtaEff_;
  MonitorElement* meExtraMTDfailExtenderPtEff_;

  // 1/total
  MonitorElement* meBTLTrackTPSimDirectMatchedEtaEff_;
  MonitorElement* meBTLTrackTPSimDirectMatchedPtEff_;
  MonitorElement* meBTLTrackTPSimOtherMatchedEtaEff_;
  MonitorElement* meBTLTrackTPSimOtherMatchedPtEff_;

  MonitorElement* meETLTrackTPSimMatchedPtposEff_;
  MonitorElement* meETL2TrackTPSimMatchedPtposEff_;
  MonitorElement* meETLTrackTPSimMatchedEtaposEff_;
  MonitorElement* meETL2TrackTPSimMatchedEtaposEff_;

  // 2/1

  MonitorElement* meBTLRecoDirectEtaEff_;
  MonitorElement* meBTLRecoDirectPtEff_;
  MonitorElement* meBTLRecoOtherEtaEff_;
  MonitorElement* meBTLRecoOtherPtEff_;

  MonitorElement* meETLRecoEtaposEff_;
  MonitorElement* meETL2RecoEtaposEff_;
  MonitorElement* meETLRecoPtposEff_;
  MonitorElement* meETL2RecoPtposEff_;

  // 2/Total

  MonitorElement* meBTLRecoDirectEtaTotEff_;
  MonitorElement* meBTLRecoDirectPtTotEff_;
  MonitorElement* meBTLRecoOtherEtaTotEff_;
  MonitorElement* meBTLRecoOtherPtTotEff_;

  MonitorElement* meETLRecoEtaposTotEff_;
  MonitorElement* meETL2RecoEtaposTotEff_;
  MonitorElement* meETLRecoPtposTotEff_;
  MonitorElement* meETL2RecoPtposTotEff_;

  MonitorElement* meBTLRecoNoRecoDirectEtaTotEff_;
  MonitorElement* meBTLRecoNoRecoDirectPtTotEff_;
  MonitorElement* meBTLRecoNoRecoOtherEtaTotEff_;
  MonitorElement* meBTLRecoNoRecoOtherPtTotEff_;

  MonitorElement* meETLRecoNoRecoEtaposTotEff_;
  MonitorElement* meETL2RecoNoRecoEtaposTotEff_;
  MonitorElement* meETLRecoNoRecoPtposTotEff_;
  MonitorElement* meETL2RecoNoRecoPtposTotEff_;

  // 3/2

  MonitorElement* meBTLCorrectSimDirectEtaEff_;
  MonitorElement* meBTLCorrectSimDirectPtEff_;
  MonitorElement* meBTLCorrectSimOtherEtaEff_;
  MonitorElement* meBTLCorrectSimOtherPtEff_;

  MonitorElement* meETLCorrectSimEtaposEff_;
  MonitorElement* meETL2CorrectSimEtaposEff_;
  MonitorElement* meETLCorrectSimPtposEff_;
  MonitorElement* meETL2CorrectSimPtposEff_;

  // 3/total

  MonitorElement* meBTLCorrectSimDirectEtaTotEff_;
  MonitorElement* meBTLCorrectSimDirectPtTotEff_;
  MonitorElement* meBTLCorrectSimOtherEtaTotEff_;
  MonitorElement* meBTLCorrectSimOtherPtTotEff_;

  MonitorElement* meETLCorrectSimEtaposTotEff_;
  MonitorElement* meETL2CorrectSimEtaposTotEff_;
  MonitorElement* meETLCorrectSimPtposTotEff_;
  MonitorElement* meETL2CorrectSimPtposTotEff_;

  // 3/reco

  MonitorElement* meBTLCorrectSimDirectEtaMtdEff_;
  MonitorElement* meBTLCorrectSimDirectPtMtdEff_;
  MonitorElement* meBTLCorrectSimOtherEtaMtdEff_;
  MonitorElement* meBTLCorrectSimOtherPtMtdEff_;

  MonitorElement* meETLCorrectSimEtaposMtdEff_;
  MonitorElement* meETL2CorrectSimEtaposMtdEff_;
  MonitorElement* meETLCorrectSimPtposMtdEff_;
  MonitorElement* meETL2CorrectSimPtposMtdEff_;

  // 5/reco

  MonitorElement* meBTLIncorrectSimEtaMtdEff_;
  MonitorElement* meBTLIncorrectSimPtMtdEff_;

  MonitorElement* meETLIncorrectSimEtaposMtdEff_;
  MonitorElement* meETLIncorrectSimPtposMtdEff_;

  MonitorElement* meBTLIncorrectSimEtaTotEff_;
  MonitorElement* meBTLIncorrectSimDirectEtaTotEff_;
  MonitorElement* meBTLIncorrectSimSecondEtaTotEff_;
  MonitorElement* meBTLIncorrectSimBackEtaTotEff_;
  MonitorElement* meBTLIncorrectSimLoopEtaTotEff_;
  MonitorElement* meBTLIncorrectSimPtTotEff_;
  MonitorElement* meBTLIncorrectSimDirectPtTotEff_;
  MonitorElement* meBTLIncorrectSimSecondPtTotEff_;
  MonitorElement* meBTLIncorrectSimBackPtTotEff_;
  MonitorElement* meBTLIncorrectSimLoopPtTotEff_;

  MonitorElement* meETLIncorrectSimEtaposTotEff_;
  MonitorElement* meETLIncorrectSimPtposTotEff_;

  MonitorElement* meBTLTrackTPNoSimEtaEff_;
  MonitorElement* meBTLTrackTPNoSimPtEff_;

  MonitorElement* meETLTrackTPNoSimEtaposEff_;
  MonitorElement* meETLTrackTPNoSimPtposEff_;
};

// ------------ constructor and destructor --------------
MtdTracksHarvester::MtdTracksHarvester(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")) {}

MtdTracksHarvester::~MtdTracksHarvester() {}

// auxiliary method to compute efficiency from the ratio of two 1D MonitorElement
void MtdTracksHarvester::computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result) {
  for (int ibin = 1; ibin <= den->getNbinsX(); ibin++) {
    double eff = num->getBinContent(ibin) / den->getBinContent(ibin);
    double bin_err = sqrt((num->getBinContent(ibin) * (den->getBinContent(ibin) - num->getBinContent(ibin))) /
                          pow(den->getBinContent(ibin), 3));
    if (den->getBinContent(ibin) == 0) {
      eff = 0;
      bin_err = 0;
    }
    result->setBinContent(ibin, eff);
    result->setBinError(ibin, bin_err);
  }
}

void MtdTracksHarvester::normalize(MonitorElement* h, double scale) {
  double integral = h->getTH1F()->Integral();
  double norma = (integral > 0.) ? scale / integral : 0.;
  for (int ibin = 1; ibin <= h->getNbinsX(); ibin++) {
    double eff = h->getBinContent(ibin) * norma;
    double bin_err = h->getBinError(ibin) * norma;
    h->setBinContent(ibin, eff);
    h->setBinError(ibin, bin_err);
  }
}

// ------------ endjob tasks ----------------------------
void MtdTracksHarvester::dqmEndJob(DQMStore::IBooker& ibook, DQMStore::IGetter& igetter) {
  // --- Get the monitoring histograms
  MonitorElement* meBTLTrackEffEtaTot = igetter.get(folder_ + "TrackBTLEffEtaTot");
  MonitorElement* meBTLTrackEffPhiTot = igetter.get(folder_ + "TrackBTLEffPhiTot");
  MonitorElement* meBTLTrackEffPtTot = igetter.get(folder_ + "TrackBTLEffPtTot");
  MonitorElement* meBTLTrackEffEtaMtd = igetter.get(folder_ + "TrackBTLEffEtaMtd");
  MonitorElement* meBTLTrackEffPhiMtd = igetter.get(folder_ + "TrackBTLEffPhiMtd");
  MonitorElement* meBTLTrackEffPtMtd = igetter.get(folder_ + "TrackBTLEffPtMtd");
  MonitorElement* meETLTrackEffEtaTotZneg = igetter.get(folder_ + "TrackETLEffEtaTotZneg");
  MonitorElement* meETLTrackEffEtaTotLowPt0 = igetter.get(folder_ + "TrackETLEffEtaTotLowPt0");
  MonitorElement* meETLTrackEffEtaTotLowPt1 = igetter.get(folder_ + "TrackETLEffEtaTotLowPt1");
  MonitorElement* meETLTrackEffPhiTotZneg = igetter.get(folder_ + "TrackETLEffPhiTotZneg");
  MonitorElement* meETLTrackEffPtTotZneg = igetter.get(folder_ + "TrackETLEffPtTotZneg");
  MonitorElement* meETLTrackEffEtaMtdZneg = igetter.get(folder_ + "TrackETLEffEtaMtdZneg");
  MonitorElement* meETLTrackEffEtaMtdLowPt0 = igetter.get(folder_ + "TrackETLEffEtaMtdLowPt0");
  MonitorElement* meETLTrackEffEtaMtdLowPt1 = igetter.get(folder_ + "TrackETLEffEtaMtdLowPt1");
  MonitorElement* meETLTrackEffPhiMtdZneg = igetter.get(folder_ + "TrackETLEffPhiMtdZneg");
  MonitorElement* meETLTrackEffPtMtdZneg = igetter.get(folder_ + "TrackETLEffPtMtdZneg");
  MonitorElement* meETLTrackEffEta2MtdZneg = igetter.get(folder_ + "TrackETLEffEta2MtdZneg");
  MonitorElement* meETLTrackEffEta2MtdLowPt0 = igetter.get(folder_ + "TrackETLEffEta2MtdLowPt0");
  MonitorElement* meETLTrackEffEta2MtdLowPt1 = igetter.get(folder_ + "TrackETLEffEta2MtdLowPt1");
  MonitorElement* meETLTrackEffPhi2MtdZneg = igetter.get(folder_ + "TrackETLEffPhi2MtdZneg");
  MonitorElement* meETLTrackEffPt2MtdZneg = igetter.get(folder_ + "TrackETLEffPt2MtdZneg");
  MonitorElement* meETLTrackEffEtaTotZpos = igetter.get(folder_ + "TrackETLEffEtaTotZpos");
  MonitorElement* meETLTrackEffPhiTotZpos = igetter.get(folder_ + "TrackETLEffPhiTotZpos");
  MonitorElement* meETLTrackEffPtTotZpos = igetter.get(folder_ + "TrackETLEffPtTotZpos");
  MonitorElement* meETLTrackEffEtaMtdZpos = igetter.get(folder_ + "TrackETLEffEtaMtdZpos");
  MonitorElement* meETLTrackEffPhiMtdZpos = igetter.get(folder_ + "TrackETLEffPhiMtdZpos");
  MonitorElement* meETLTrackEffPtMtdZpos = igetter.get(folder_ + "TrackETLEffPtMtdZpos");
  MonitorElement* meETLTrackEffEta2MtdZpos = igetter.get(folder_ + "TrackETLEffEta2MtdZpos");
  MonitorElement* meETLTrackEffPhi2MtdZpos = igetter.get(folder_ + "TrackETLEffPhi2MtdZpos");
  MonitorElement* meETLTrackEffPt2MtdZpos = igetter.get(folder_ + "TrackETLEffPt2MtdZpos");
  MonitorElement* meTrackPtTot = igetter.get(folder_ + "TrackPtTot");
  MonitorElement* meExtraPtMtd = igetter.get(folder_ + "ExtraPtMtd");
  MonitorElement* meExtraPtEtl2Mtd = igetter.get(folder_ + "ExtraPtEtl2Mtd");
  MonitorElement* meTrackMatchedTPEffPtTot = igetter.get(folder_ + "MatchedTPEffPtTot");
  MonitorElement* meTrackMatchedTPEffPtTotLV = igetter.get(folder_ + "MatchedTPEffPtTotLV");
  MonitorElement* meTrackMatchedTPEffPtMtd = igetter.get(folder_ + "MatchedTPEffPtMtd");
  MonitorElement* meTrackMatchedTPEffPtEtl2Mtd = igetter.get(folder_ + "MatchedTPEffPtEtl2Mtd");
  MonitorElement* meTrackMatchedTPmtdEffPtTot = igetter.get(folder_ + "MatchedTPmtdEffPtTot");
  MonitorElement* meTrackMatchedTPmtdEffPtMtd = igetter.get(folder_ + "MatchedTPmtdEffPtMtd");
  MonitorElement* meTrackEtaTot = igetter.get(folder_ + "TrackEtaTot");
  MonitorElement* meExtraEtaMtd = igetter.get(folder_ + "ExtraEtaMtd");
  MonitorElement* meExtraEtaEtl2Mtd = igetter.get(folder_ + "ExtraEtaEtl2Mtd");
  MonitorElement* meTrackMatchedTPEffEtaTot = igetter.get(folder_ + "MatchedTPEffEtaTot");
  MonitorElement* meTrackMatchedTPEffEtaTotLV = igetter.get(folder_ + "MatchedTPEffEtaTotLV");
  MonitorElement* meTrackMatchedTPEffEtaMtd = igetter.get(folder_ + "MatchedTPEffEtaMtd");
  MonitorElement* meTrackMatchedTPEffEtaEtl2Mtd = igetter.get(folder_ + "MatchedTPEffEtaEtl2Mtd");
  MonitorElement* meTrackMatchedTPmtdEffEtaTot = igetter.get(folder_ + "MatchedTPmtdEffEtaTot");
  MonitorElement* meTrackMatchedTPmtdEffEtaMtd = igetter.get(folder_ + "MatchedTPmtdEffEtaMtd");
  MonitorElement* meTrackNumHits = igetter.get(folder_ + "TrackNumHits");
  MonitorElement* meTrackNumHitsNT = igetter.get(folder_ + "TrackNumHitsNT");
  MonitorElement* meExtraPhiAtBTL = igetter.get(folder_ + "ExtraPhiAtBTL");
  MonitorElement* meExtraPhiAtBTLmatched = igetter.get(folder_ + "ExtraPhiAtBTLmatched");
  MonitorElement* meExtraBTLeneInCone = igetter.get(folder_ + "ExtraBTLeneInCone");
  MonitorElement* meExtraMTDfailExtenderEta = igetter.get(folder_ + "ExtraMTDfailExtenderEta");
  MonitorElement* meExtraMTDfailExtenderPt = igetter.get(folder_ + "ExtraMTDfailExtenderPt");

  // 1
  MonitorElement* meBTLTrackTPSimDirectMatchedEta = igetter.get(folder_ + "BTLTrackTPSimDirectMatchedEta");
  MonitorElement* meBTLTrackTPSimOtherMatchedEta = igetter.get(folder_ + "BTLTrackTPSimOtherMatchedEta");
  MonitorElement* meBTLTrackTPSimDirectMatchedPt = igetter.get(folder_ + "BTLTrackTPSimDirectMatchedPt");
  MonitorElement* meBTLTrackTPSimOtherMatchedPt = igetter.get(folder_ + "BTLTrackTPSimOtherMatchedPt");
  MonitorElement* meETLTrackTPSimMatchedEtapos = igetter.get(folder_ + "ETLTrackTPSimMatchedEtapos");
  MonitorElement* meETLTrackTPSimMatchedPtpos = igetter.get(folder_ + "ETLTrackTPSimMatchedPtpos");
  MonitorElement* meETL2TrackTPSimMatchedEtapos = igetter.get(folder_ + "ETL2TrackTPSimMatchedEtapos");
  MonitorElement* meETL2TrackTPSimMatchedPtpos = igetter.get(folder_ + "ETL2TrackTPSimMatchedPtpos");

  // 2
  MonitorElement* meBTLTrackTPSimRecoDirectMatchedEta = igetter.get(folder_ + "BTLTrackTPSimRecoDirectMatchedEta");
  MonitorElement* meBTLTrackTPSimRecoOtherMatchedEta = igetter.get(folder_ + "BTLTrackTPSimRecoOtherMatchedEta");
  MonitorElement* meBTLTrackTPSimRecoDirectMatchedPt = igetter.get(folder_ + "BTLTrackTPSimRecoDirectMatchedPt");
  MonitorElement* meBTLTrackTPSimRecoOtherMatchedPt = igetter.get(folder_ + "BTLTrackTPSimRecoOtherMatchedPt");
  MonitorElement* meETLTrackTPSimRecoMatchedEtapos = igetter.get(folder_ + "ETLTrackTPSimRecoMatchedEtapos");
  MonitorElement* meETLTrackTPSimRecoMatchedPtpos = igetter.get(folder_ + "ETLTrackTPSimRecoMatchedPtpos");
  MonitorElement* meETL2TrackTPSimRecoMatchedEtapos = igetter.get(folder_ + "ETL2TrackTPSimRecoMatchedEtapos");
  MonitorElement* meETL2TrackTPSimRecoMatchedPtpos = igetter.get(folder_ + "ETL2TrackTPSimRecoMatchedPtpos");

  MonitorElement* meBTLTrackTPSimNoRecoDirectMatchedEta = igetter.get(folder_ + "BTLTrackTPSimNoRecoDirectMatchedEta");
  MonitorElement* meBTLTrackTPSimNoRecoDirectMatchedPt = igetter.get(folder_ + "BTLTrackTPSimNoRecoDirectMatchedPt");
  MonitorElement* meBTLTrackTPSimNoRecoOtherMatchedEta = igetter.get(folder_ + "BTLTrackTPSimNoRecoOtherMatchedEta");
  MonitorElement* meBTLTrackTPSimNoRecoOtherMatchedPt = igetter.get(folder_ + "BTLTrackTPSimNoRecoOtherMatchedPt");
  MonitorElement* meETLTrackTPSimNoRecoMatchedEtapos = igetter.get(folder_ + "ETLTrackTPSimNoRecoMatchedEtapos");
  MonitorElement* meETLTrackTPSimNoRecoMatchedPtpos = igetter.get(folder_ + "ETLTrackTPSimNoRecoMatchedPtpos");
  MonitorElement* meETL2TrackTPSimNoRecoMatchedEtapos = igetter.get(folder_ + "ETL2TrackTPSimNoRecoMatchedEtapos");
  MonitorElement* meETL2TrackTPSimNoRecoMatchedPtpos = igetter.get(folder_ + "ETL2TrackTPSimNoRecoMatchedPtpos");

  // 3
  MonitorElement* meBTLTrackSimToSimDirectMatchedEta = igetter.get(folder_ + "BTLTrackSimToSimDirectMatchedEta");
  MonitorElement* meBTLTrackSimToSimOtherMatchedEta = igetter.get(folder_ + "BTLTrackSimToSimOtherMatchedEta");
  MonitorElement* meBTLTrackSimToSimDirectMatchedPt = igetter.get(folder_ + "BTLTrackSimToSimDirectMatchedPt");
  MonitorElement* meBTLTrackSimToSimOtherMatchedPt = igetter.get(folder_ + "BTLTrackSimToSimOtherMatchedPt");
  MonitorElement* meETLTrackSimToSimMatchedEtapos = igetter.get(folder_ + "ETLTrackSimToSimMatchedEtapos");
  MonitorElement* meETLTrackSimToSimMatchedPtpos = igetter.get(folder_ + "ETLTrackSimToSimMatchedPtpos");
  MonitorElement* meETL2TrackSimToSimMatchedEtapos = igetter.get(folder_ + "ETL2TrackSimToSimMatchedEtapos");
  MonitorElement* meETL2TrackSimToSimMatchedPtpos = igetter.get(folder_ + "ETL2TrackSimToSimMatchedPtpos");

  // 4

  MonitorElement* meBTLTrackTPNoSimEta = igetter.get(folder_ + "BTLTrackTPNoSimEta");
  MonitorElement* meBTLTrackTPNoSimPt = igetter.get(folder_ + "BTLTrackTPNoSimPt");
  MonitorElement* meETLTrackTPNoSimEtapos = igetter.get(folder_ + "ETLTrackTPNoSimEtapos");
  MonitorElement* meETLTrackTPNoSimPtpos = igetter.get(folder_ + "ETLTrackTPNoSimPtpos");

  // 5

  MonitorElement* meBTLTrackNoSimToSimMatchedEta = igetter.get(folder_ + "BTLTrackNoSimToSimMatchedEta");
  MonitorElement* meBTLTrackNoSimToSimDirectMatchedEta = igetter.get(folder_ + "BTLTrackNoSimToSimDirectMatchedEta");
  MonitorElement* meBTLTrackNoSimToSimSecondMatchedEta = igetter.get(folder_ + "BTLTrackNoSimToSimSecondMatchedEta");
  MonitorElement* meBTLTrackNoSimToSimBackMatchedEta = igetter.get(folder_ + "BTLTrackNoSimToSimBackMatchedEta");
  MonitorElement* meBTLTrackNoSimToSimLoopMatchedEta = igetter.get(folder_ + "BTLTrackNoSimToSimLoopMatchedEta");

  MonitorElement* meBTLTrackNoSimToSimMatchedPt = igetter.get(folder_ + "BTLTrackNoSimToSimMatchedPt");
  MonitorElement* meBTLTrackNoSimToSimDirectMatchedPt = igetter.get(folder_ + "BTLTrackNoSimToSimDirectMatchedPt");
  MonitorElement* meBTLTrackNoSimToSimSecondMatchedPt = igetter.get(folder_ + "BTLTrackNoSimToSimSecondMatchedPt");
  MonitorElement* meBTLTrackNoSimToSimBackMatchedPt = igetter.get(folder_ + "BTLTrackNoSimToSimBackMatchedPt");
  MonitorElement* meBTLTrackNoSimToSimLoopMatchedPt = igetter.get(folder_ + "BTLTrackNoSimToSimLoopMatchedPt");

  MonitorElement* meETLTrackNoSimToSimMatchedEtapos = igetter.get(folder_ + "ETLTrackNoSimToSimMatchedEtapos");
  MonitorElement* meETLTrackNoSimToSimMatchedPtpos = igetter.get(folder_ + "ETLTrackNoSimToSimMatchedPtpos");

  if (!meBTLTrackEffEtaTot || !meBTLTrackEffPhiTot || !meBTLTrackEffPtTot || !meBTLTrackEffEtaMtd ||
      !meBTLTrackEffPhiMtd || !meBTLTrackEffPtMtd || !meETLTrackEffEtaTotZneg || !meETLTrackEffPhiTotZneg ||
      !meETLTrackEffPtTotZneg || !meETLTrackEffEtaTotLowPt0 || !meETLTrackEffEtaTotLowPt1 || !meETLTrackEffEtaMtdZneg ||
      !meETLTrackEffEtaMtdLowPt0 || !meETLTrackEffEtaMtdLowPt1 || !meETLTrackEffPhiMtdZneg || !meETLTrackEffPtMtdZneg ||
      !meETLTrackEffEta2MtdZneg || !meETLTrackEffEta2MtdLowPt0 || !meETLTrackEffEta2MtdLowPt1 ||
      !meETLTrackEffPhi2MtdZneg || !meETLTrackEffPt2MtdZneg || !meETLTrackEffEtaTotZpos || !meETLTrackEffPhiTotZpos ||
      !meETLTrackEffPtTotZpos || !meETLTrackEffEtaMtdZpos || !meETLTrackEffPhiMtdZpos || !meETLTrackEffPtMtdZpos ||
      !meETLTrackEffEta2MtdZpos || !meETLTrackEffPhi2MtdZpos || !meETLTrackEffPt2MtdZpos || !meTrackMatchedTPEffPtTot ||
      !meTrackMatchedTPEffPtTotLV || !meTrackMatchedTPEffPtMtd || !meTrackMatchedTPEffPtEtl2Mtd ||
      !meTrackMatchedTPmtdEffPtTot || !meTrackMatchedTPmtdEffPtMtd || !meTrackMatchedTPEffEtaTot ||
      !meTrackMatchedTPEffEtaTotLV || !meTrackMatchedTPEffEtaMtd || !meTrackMatchedTPEffEtaEtl2Mtd ||
      !meTrackMatchedTPmtdEffEtaTot || !meTrackMatchedTPmtdEffEtaMtd || !meTrackNumHits || !meTrackNumHitsNT ||
      !meTrackPtTot || !meTrackEtaTot || !meExtraPtMtd || !meExtraPtEtl2Mtd || !meExtraEtaMtd || !meExtraEtaEtl2Mtd ||
      !meExtraPhiAtBTL || !meExtraPhiAtBTLmatched || !meExtraBTLeneInCone || !meExtraMTDfailExtenderEta ||
      !meExtraMTDfailExtenderPt || !meBTLTrackTPSimDirectMatchedPt || !meBTLTrackTPSimDirectMatchedEta ||
      !meBTLTrackTPSimOtherMatchedPt || !meBTLTrackTPSimOtherMatchedEta || !meETLTrackTPSimMatchedEtapos ||
      !meETLTrackTPSimMatchedPtpos || !meETL2TrackTPSimMatchedEtapos || !meETL2TrackTPSimMatchedPtpos ||
      !meBTLTrackTPSimRecoDirectMatchedEta || !meBTLTrackTPSimRecoDirectMatchedPt ||
      !meBTLTrackTPSimRecoOtherMatchedEta || !meBTLTrackTPSimRecoOtherMatchedPt || !meETLTrackTPSimRecoMatchedEtapos ||
      !meETLTrackTPSimRecoMatchedPtpos || !meETL2TrackTPSimRecoMatchedEtapos || !meETL2TrackTPSimRecoMatchedPtpos ||
      !meBTLTrackTPSimNoRecoDirectMatchedEta || !meBTLTrackTPSimNoRecoDirectMatchedPt ||
      !meBTLTrackTPSimNoRecoOtherMatchedEta || !meBTLTrackTPSimNoRecoOtherMatchedPt ||
      !meETLTrackTPSimNoRecoMatchedEtapos || !meETLTrackTPSimNoRecoMatchedPtpos ||
      !meETL2TrackTPSimNoRecoMatchedEtapos || !meETL2TrackTPSimNoRecoMatchedPtpos || !meETLTrackTPNoSimEtapos ||
      !meETLTrackTPNoSimPtpos || !meBTLTrackTPNoSimEta || !meBTLTrackTPNoSimPt || !meBTLTrackNoSimToSimMatchedEta ||
      !meBTLTrackNoSimToSimDirectMatchedEta || !meBTLTrackNoSimToSimSecondMatchedEta ||
      !meBTLTrackNoSimToSimBackMatchedEta || !meBTLTrackNoSimToSimLoopMatchedEta || !meBTLTrackNoSimToSimMatchedPt ||
      !meBTLTrackNoSimToSimDirectMatchedPt || !meBTLTrackNoSimToSimSecondMatchedPt ||
      !meBTLTrackNoSimToSimBackMatchedPt || !meBTLTrackNoSimToSimLoopMatchedPt || !meETLTrackNoSimToSimMatchedEtapos ||
      !meETLTrackNoSimToSimMatchedPtpos) {
    edm::LogError("MtdTracksHarvester") << "Monitoring histograms not found!" << std::endl;
    return;
  }

  // --- Book  histograms
  ibook.cd(folder_);

  meBtlEtaEff_ = ibook.book1D("BtlEtaEff",
                              " Track Efficiency VS Eta;#eta;Efficiency",
                              meBTLTrackEffEtaTot->getNbinsX(),
                              meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                              meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBtlEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffEtaMtd, meBTLTrackEffEtaTot, meBtlEtaEff_);

  meBtlPhiEff_ = ibook.book1D("BtlPhiEff",
                              "Track Efficiency VS Phi;#phi [rad];Efficiency",
                              meBTLTrackEffPhiTot->getNbinsX(),
                              meBTLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmin(),
                              meBTLTrackEffPhiTot->getTH1()->GetXaxis()->GetXmax());
  meBtlPhiEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffPhiMtd, meBTLTrackEffPhiTot, meBtlPhiEff_);

  meBtlPtEff_ = ibook.book1D("BtlPtEff",
                             "Track Efficiency VS Pt;Pt [GeV];Efficiency",
                             meBTLTrackEffPtTot->getNbinsX(),
                             meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                             meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBtlPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackEffPtMtd, meBTLTrackEffPtTot, meBtlPtEff_);

  meEtlEtaEff_[0] = ibook.book1D("EtlEtaEffZneg",
                                 " Track Efficiency VS Eta (-Z);#eta;Efficiency",
                                 meETLTrackEffEtaTotZneg->getNbinsX(),
                                 meETLTrackEffEtaTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffEtaTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtdZneg, meETLTrackEffEtaTotZneg, meEtlEtaEff_[0]);

  meEtlEtaEffLowPt_[0] = ibook.book1D("EtlEtaEffLowPt0",
                                      " Track Efficiency VS Eta, 0.2 < pt < 0.45;#eta;Efficiency",
                                      meETLTrackEffEtaTotLowPt0->getNbinsX(),
                                      meETLTrackEffEtaTotLowPt0->getTH1()->GetXaxis()->GetXmin(),
                                      meETLTrackEffEtaTotLowPt0->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEffLowPt_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtdLowPt0, meETLTrackEffEtaTotLowPt0, meEtlEtaEffLowPt_[0]);

  meEtlPhiEff_[0] = ibook.book1D("EtlPhiEffZneg",
                                 "Track Efficiency VS Phi (-Z);#phi [rad];Efficiency",
                                 meETLTrackEffPhiTotZneg->getNbinsX(),
                                 meETLTrackEffPhiTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffPhiTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlPhiEff_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPhiMtdZneg, meETLTrackEffPhiTotZneg, meEtlPhiEff_[0]);

  meEtlPtEff_[0] = ibook.book1D("EtlPtEffZneg",
                                "Track Efficiency VS Pt (-Z);Pt [GeV];Efficiency",
                                meETLTrackEffPtTotZneg->getNbinsX(),
                                meETLTrackEffPtTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                meETLTrackEffPtTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlPtEff_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPtMtdZneg, meETLTrackEffPtTotZneg, meEtlPtEff_[0]);

  meEtlEtaEff_[1] = ibook.book1D("EtlEtaEffZpos",
                                 " Track Efficiency VS Eta (+Z);#eta;Efficiency",
                                 meETLTrackEffEtaTotZpos->getNbinsX(),
                                 meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtdZpos, meETLTrackEffEtaTotZpos, meEtlEtaEff_[1]);

  meEtlEtaEffLowPt_[1] = ibook.book1D("EtlEtaEffLowPt1",
                                      " Track Efficiency VS Eta, 0.45 < pt < 0.7;#eta;Efficiency",
                                      meETLTrackEffEtaTotLowPt1->getNbinsX(),
                                      meETLTrackEffEtaTotLowPt1->getTH1()->GetXaxis()->GetXmin(),
                                      meETLTrackEffEtaTotLowPt1->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEffLowPt_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEtaMtdLowPt1, meETLTrackEffEtaTotLowPt1, meEtlEtaEffLowPt_[1]);

  meEtlPhiEff_[1] = ibook.book1D("EtlPhiEffZpos",
                                 "Track Efficiency VS Phi (+Z);#phi [rad];Efficiency",
                                 meETLTrackEffPhiTotZpos->getNbinsX(),
                                 meETLTrackEffPhiTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffPhiTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlPhiEff_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPhiMtdZpos, meETLTrackEffPhiTotZpos, meEtlPhiEff_[1]);

  meEtlPtEff_[1] = ibook.book1D("EtlPtEffZpos",
                                "Track Efficiency VS Pt (+Z);Pt [GeV];Efficiency",
                                meETLTrackEffPtTotZpos->getNbinsX(),
                                meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlPtEff_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPtMtdZpos, meETLTrackEffPtTotZpos, meEtlPtEff_[1]);

  meEtlEtaEff2_[0] = ibook.book1D("EtlEtaEff2Zneg",
                                  " Track Efficiency VS Eta (-Z, 2 hit);#eta;Efficiency",
                                  meETLTrackEffEtaTotZneg->getNbinsX(),
                                  meETLTrackEffEtaTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                  meETLTrackEffEtaTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff2_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEta2MtdZneg, meETLTrackEffEtaTotZneg, meEtlEtaEff2_[0]);

  meEtlEtaEff2LowPt_[0] = ibook.book1D("EtlEtaEff2LowPt0",
                                       " Track Efficiency VS Eta (2 hits), 0.2 < pt < 0.45;#eta;Efficiency",
                                       meETLTrackEffEtaTotLowPt0->getNbinsX(),
                                       meETLTrackEffEtaTotLowPt0->getTH1()->GetXaxis()->GetXmin(),
                                       meETLTrackEffEtaTotLowPt0->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff2LowPt_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEta2MtdLowPt0, meETLTrackEffEtaTotLowPt0, meEtlEtaEff2LowPt_[0]);

  meEtlPhiEff2_[0] = ibook.book1D("EtlPhiEff2Zneg",
                                  "Track Efficiency VS Phi (-Z, 2 hits);#phi [rad];Efficiency",
                                  meETLTrackEffPhiTotZneg->getNbinsX(),
                                  meETLTrackEffPhiTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                  meETLTrackEffPhiTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlPhiEff2_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPhi2MtdZneg, meETLTrackEffPhiTotZneg, meEtlPhiEff2_[0]);

  meEtlPtEff2_[0] = ibook.book1D("EtlPtEff2Zneg",
                                 "Track Efficiency VS Pt (-Z, 2 hits);Pt [GeV];Efficiency",
                                 meETLTrackEffPtTotZneg->getNbinsX(),
                                 meETLTrackEffPtTotZneg->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffPtTotZneg->getTH1()->GetXaxis()->GetXmax());
  meEtlPtEff2_[0]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPt2MtdZneg, meETLTrackEffPtTotZneg, meEtlPtEff2_[0]);

  meEtlEtaEff2_[1] = ibook.book1D("EtlEtaEff2Zpos",
                                  "Track Efficiency VS Eta (+Z, 2 hits);#eta;Efficiency",
                                  meETLTrackEffEtaTotZpos->getNbinsX(),
                                  meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                  meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff2_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEta2MtdZpos, meETLTrackEffEtaTotZpos, meEtlEtaEff2_[1]);

  meEtlEtaEff2LowPt_[1] = ibook.book1D("EtlEtaEff2LowPt1",
                                       " Track Efficiency VS Eta (2 hits), 0.45 < pt < 0.7;#eta;Efficiency",
                                       meETLTrackEffEtaTotLowPt1->getNbinsX(),
                                       meETLTrackEffEtaTotLowPt1->getTH1()->GetXaxis()->GetXmin(),
                                       meETLTrackEffEtaTotLowPt1->getTH1()->GetXaxis()->GetXmax());
  meEtlEtaEff2LowPt_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffEta2MtdLowPt1, meETLTrackEffEtaTotLowPt1, meEtlEtaEff2LowPt_[1]);

  meEtlPhiEff2_[1] = ibook.book1D("EtlPhiEff2Zpos",
                                  "Track Efficiency VS Phi (+Z, 2 hits);#phi [rad];Efficiency",
                                  meETLTrackEffPhiTotZpos->getNbinsX(),
                                  meETLTrackEffPhiTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                  meETLTrackEffPhiTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlPhiEff2_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPhi2MtdZpos, meETLTrackEffPhiTotZpos, meEtlPhiEff2_[1]);

  meEtlPtEff2_[1] = ibook.book1D("EtlPtEff2Zpos",
                                 "Track Efficiency VS Pt (+Z, 2 hits);Pt [GeV];Efficiency",
                                 meETLTrackEffPtTotZpos->getNbinsX(),
                                 meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                 meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meEtlPtEff2_[1]->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackEffPt2MtdZpos, meETLTrackEffPtTotZpos, meEtlPtEff2_[1]);

  meExtraPtEff_ =
      ibook.book1D("ExtraPtEff",
                   "MTD matching efficiency wrt extrapolated track associated to LV VS Pt;Pt [GeV];Efficiency",
                   meTrackMatchedTPEffEtaTotLV->getNbinsX(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmin(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraPtMtd, meTrackMatchedTPEffPtTotLV, meExtraPtEff_);

  meExtraPtEtl2Eff_ =
      ibook.book1D("ExtraPtEtl2Eff",
                   "MTD matching efficiency (2 ETL) wrt extrapolated track associated to LV VS Pt;Pt [GeV];Efficiency",
                   meTrackMatchedTPEffEtaTotLV->getNbinsX(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmin(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraPtEtl2Eff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraPtEtl2Mtd, meTrackMatchedTPEffPtTotLV, meExtraPtEtl2Eff_);

  meExtraEtaEff_ = ibook.book1D("ExtraEtaEff",
                                "MTD matching efficiency wrt extrapolated track associated to LV VS Eta;Eta;Efficiency",
                                meTrackMatchedTPEffEtaTotLV->getNbinsX(),
                                meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmin(),
                                meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraEtaMtd, meTrackMatchedTPEffEtaTotLV, meExtraEtaEff_);

  meExtraEtaEtl2Eff_ =
      ibook.book1D("ExtraEtaEtl2Eff",
                   "MTD matching efficiency (2 ETL) wrt extrapolated track associated to LV VS Eta;Eta;Efficiency",
                   meTrackMatchedTPEffEtaTotLV->getNbinsX(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmin(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraEtaEtl2Eff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraEtaEtl2Mtd, meTrackMatchedTPEffEtaTotLV, meExtraEtaEtl2Eff_);

  meTPPtSelEff_ = ibook.book1D("TPPtSelEff",
                               "Track selected efficiency TP VS Pt;Pt [GeV];Efficiency",
                               meTrackPtTot->getNbinsX(),
                               meTrackPtTot->getTH1()->GetXaxis()->GetXmin(),
                               meTrackPtTot->getTH1()->GetXaxis()->GetXmax());
  meTPPtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffPtTot, meTrackPtTot, meTPPtSelEff_);

  meTPEtaSelEff_ = ibook.book1D("TPEtaSelEff",
                                "Track selected efficiency TP VS Eta;Eta;Efficiency",
                                meTrackEtaTot->getNbinsX(),
                                meTrackEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                meTrackEtaTot->getTH1()->GetXaxis()->GetXmax());
  meTPEtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffEtaTot, meTrackEtaTot, meTPEtaSelEff_);

  meTPPtMatchEff_ = ibook.book1D("TPPtMatchEff",
                                 "Track matched to TP efficiency VS Pt;Pt [GeV];Efficiency",
                                 meTrackMatchedTPEffPtTot->getNbinsX(),
                                 meTrackMatchedTPEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                 meTrackMatchedTPEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meTPPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffPtMtd, meTrackMatchedTPEffPtTot, meTPPtMatchEff_);

  meTPEtaMatchEff_ = ibook.book1D("TPEtaMatchEff",
                                  "Track matched to TP efficiency VS Eta;Eta;Efficiency",
                                  meTrackMatchedTPEffEtaTot->getNbinsX(),
                                  meTrackMatchedTPEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                  meTrackMatchedTPEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meTPEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffEtaMtd, meTrackMatchedTPEffEtaTot, meTPEtaMatchEff_);

  meTPPtMatchEtl2Eff_ = ibook.book1D("TPPtMatchEtl2Eff",
                                     "Track matched to TP efficiency VS Pt, 2 ETL hits;Pt [GeV];Efficiency",
                                     meTrackMatchedTPEffPtTot->getNbinsX(),
                                     meTrackMatchedTPEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                     meTrackMatchedTPEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meTPPtMatchEtl2Eff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffPtEtl2Mtd, meTrackMatchedTPEffPtTot, meTPPtMatchEtl2Eff_);

  meTPEtaMatchEtl2Eff_ = ibook.book1D("TPEtaMatchEtl2Eff",
                                      "Track matched to TP efficiency VS Eta, 2 ETL hits;Eta;Efficiency",
                                      meTrackMatchedTPEffEtaTot->getNbinsX(),
                                      meTrackMatchedTPEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                      meTrackMatchedTPEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meTPEtaMatchEtl2Eff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPEffEtaEtl2Mtd, meTrackMatchedTPEffEtaTot, meTPEtaMatchEtl2Eff_);

  meTPmtdPtSelEff_ = ibook.book1D("TPmtdPtSelEff",
                                  "Track selected efficiency TP-mtd hit VS Pt;Pt [GeV];Efficiency",
                                  meTrackPtTot->getNbinsX(),
                                  meTrackPtTot->getTH1()->GetXaxis()->GetXmin(),
                                  meTrackPtTot->getTH1()->GetXaxis()->GetXmax());
  meTPmtdPtSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPmtdEffPtTot, meTrackPtTot, meTPmtdPtSelEff_);

  meTPmtdEtaSelEff_ = ibook.book1D("TPmtdEtaSelEff",
                                   "Track selected efficiency TPmtd hit VS Eta;Eta;Efficiency",
                                   meTrackEtaTot->getNbinsX(),
                                   meTrackEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                   meTrackEtaTot->getTH1()->GetXaxis()->GetXmax());
  meTPmtdEtaSelEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPmtdEffEtaTot, meTrackEtaTot, meTPmtdEtaSelEff_);

  meTPmtdPtMatchEff_ = ibook.book1D("TPmtdPtMatchEff",
                                    "Track matched to TP-mtd hit efficiency VS Pt;Pt [GeV];Efficiency",
                                    meTrackMatchedTPmtdEffPtTot->getNbinsX(),
                                    meTrackMatchedTPmtdEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                    meTrackMatchedTPmtdEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meTPmtdPtMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPmtdEffPtMtd, meTrackMatchedTPmtdEffPtTot, meTPmtdPtMatchEff_);

  meTPmtdEtaMatchEff_ = ibook.book1D("TPmtdEtaMatchEff",
                                     "Track matched to TP-mtd hit efficiency VS Eta;Eta;Efficiency",
                                     meTrackMatchedTPmtdEffEtaTot->getNbinsX(),
                                     meTrackMatchedTPmtdEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                     meTrackMatchedTPmtdEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meTPmtdEtaMatchEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackMatchedTPmtdEffEtaMtd, meTrackMatchedTPmtdEffEtaTot, meTPmtdEtaMatchEff_);

  meNoTimeFraction_ = ibook.book1D("NoTimeFraction",
                                   "Fraction of tracks with MTD hits and no time associated; Num. of hits",
                                   meTrackNumHits->getNbinsX(),
                                   meTrackNumHits->getTH1()->GetXaxis()->GetXmin(),
                                   meTrackNumHits->getTH1()->GetXaxis()->GetXmax());
  meNoTimeFraction_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meTrackNumHitsNT, meTrackNumHits, meNoTimeFraction_);

  meBtlEtaEff_->getTH1()->SetMinimum(0.);
  meBtlPhiEff_->getTH1()->SetMinimum(0.);
  meBtlPtEff_->getTH1()->SetMinimum(0.);
  for (int i = 0; i < 2; i++) {
    meEtlEtaEff_[i]->getTH1()->SetMinimum(0.);
    meEtlEtaEffLowPt_[i]->getTH1()->SetMinimum(0.);
    meEtlPhiEff_[i]->getTH1()->SetMinimum(0.);
    meEtlPtEff_[i]->getTH1()->SetMinimum(0.);
    meEtlEtaEff2_[i]->getTH1()->SetMinimum(0.);
    meEtlEtaEff2LowPt_[i]->getTH1()->SetMinimum(0.);
    meEtlPhiEff2_[i]->getTH1()->SetMinimum(0.);
    meEtlPtEff2_[i]->getTH1()->SetMinimum(0.);
  }

  meExtraPhiAtBTLEff_ = ibook.book1D("ExtraPhiAtBTLEff",
                                     "Efficiency to match hits at BTL surface of extrapolated tracks associated to LV",
                                     meExtraPhiAtBTL->getNbinsX(),
                                     meExtraPhiAtBTL->getTH1()->GetXaxis()->GetXmin(),
                                     meExtraPhiAtBTL->getTH1()->GetXaxis()->GetXmax());
  meExtraPhiAtBTLEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraPhiAtBTLmatched, meExtraPhiAtBTL, meExtraPhiAtBTLEff_);

  normalize(meExtraBTLeneInCone, 1.);

  meExtraMTDfailExtenderEtaEff_ =
      ibook.book1D("ExtraMTDfailExtenderEtaEff",
                   "Track associated to LV extrapolated at MTD surface no extender efficiency VS Eta;Eta;Efficiency",
                   meTrackMatchedTPEffEtaTotLV->getNbinsX(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmin(),
                   meTrackMatchedTPEffEtaTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraMTDfailExtenderEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraMTDfailExtenderEta, meTrackMatchedTPEffEtaTotLV, meExtraMTDfailExtenderEtaEff_);

  meExtraMTDfailExtenderPtEff_ = ibook.book1D(
      "ExtraMTDfailExtenderPtEff",
      "Track associated to LV extrapolated at MTD surface no extender efficiency VS Pt;Pt [GeV];Efficiency",
      meTrackMatchedTPEffPtTotLV->getNbinsX(),
      meTrackMatchedTPEffPtTotLV->getTH1()->GetXaxis()->GetXmin(),
      meTrackMatchedTPEffPtTotLV->getTH1()->GetXaxis()->GetXmax());
  meExtraMTDfailExtenderPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meExtraMTDfailExtenderPt, meTrackMatchedTPEffPtTotLV, meExtraMTDfailExtenderPtEff_);

  meBTLTrackTPSimDirectMatchedEtaEff_ = ibook.book1D("BTLTrackTPSimDirectMatchedEtaEff",
                                                     "Track Efficiency VS Eta;#eta;Efficiency",
                                                     meBTLTrackEffEtaTot->getNbinsX(),
                                                     meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                     meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTrackTPSimDirectMatchedEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimDirectMatchedEta, meBTLTrackEffEtaTot, meBTLTrackTPSimDirectMatchedEtaEff_);

  meBTLTrackTPSimDirectMatchedPtEff_ = ibook.book1D("BTLTrackTPSimDirectMatchedPtEff",
                                                    "Track Efficiency VS Pt;Pt;Efficiency",
                                                    meBTLTrackEffPtTot->getNbinsX(),
                                                    meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                                    meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTrackTPSimDirectMatchedPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimDirectMatchedPt, meBTLTrackEffPtTot, meBTLTrackTPSimDirectMatchedPtEff_);

  meBTLTrackTPSimOtherMatchedEtaEff_ = ibook.book1D("BTLTrackTPSimOtherMatchedEtaEff",
                                                    "Track Efficiency VS Eta;#eta;Efficiency",
                                                    meBTLTrackEffEtaTot->getNbinsX(),
                                                    meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                    meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTrackTPSimOtherMatchedEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimOtherMatchedEta, meBTLTrackEffEtaTot, meBTLTrackTPSimOtherMatchedEtaEff_);

  meBTLTrackTPSimOtherMatchedPtEff_ = ibook.book1D("BTLTrackTPSimOtherMatchedPtEff",
                                                   "Track Efficiency VS Pt;Pt;Efficiency",
                                                   meBTLTrackEffPtTot->getNbinsX(),
                                                   meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                                   meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTrackTPSimOtherMatchedPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimOtherMatchedPt, meBTLTrackEffPtTot, meBTLTrackTPSimOtherMatchedPtEff_);

  meETLTrackTPSimMatchedEtaposEff_ = ibook.book1D("ETLTrackTPSimMatchedEtaposEff",
                                                  "Track Efficiency VS Eta (positive);#eta;Efficiency",
                                                  meETLTrackEffEtaTotZpos->getNbinsX(),
                                                  meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                                  meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLTrackTPSimMatchedEtaposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPSimMatchedEtapos, meETLTrackEffEtaTotZpos, meETLTrackTPSimMatchedEtaposEff_);

  meETL2TrackTPSimMatchedEtaposEff_ = ibook.book1D("ETL2TrackTPSimMatchedEtaposEff",
                                                   "Two hits Track Efficiency VS Eta (positive);#eta;Efficiency",
                                                   meETLTrackEffEtaTotZpos->getNbinsX(),
                                                   meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                                   meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2TrackTPSimMatchedEtaposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackTPSimMatchedEtapos, meETLTrackEffEtaTotZpos, meETL2TrackTPSimMatchedEtaposEff_);

  meETLTrackTPSimMatchedPtposEff_ = ibook.book1D("ETLTrackTPSimMatchedPtposEff",
                                                 "Track Efficiency VS Pt (positive);Pt;Efficiency",
                                                 meETLTrackEffPtTotZpos->getNbinsX(),
                                                 meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                                 meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLTrackTPSimMatchedPtposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPSimMatchedPtpos, meETLTrackEffPtTotZpos, meETLTrackTPSimMatchedPtposEff_);

  meETL2TrackTPSimMatchedPtposEff_ = ibook.book1D("ETL2TrackTPSimMatchedPtposEff",
                                                  "Two hits Track Efficiency VS Pt (positive);Pt;Efficiency",
                                                  meETLTrackEffPtTotZpos->getNbinsX(),
                                                  meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                                  meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2TrackTPSimMatchedPtposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackTPSimMatchedPtpos, meETLTrackEffPtTotZpos, meETL2TrackTPSimMatchedPtposEff_);

  // 2/1
  meBTLRecoDirectEtaEff_ = ibook.book1D("BTLRecoDirectEtaEff",
                                        "Reconstruction Efficiency VS Eta;#eta;Efficiency",
                                        meBTLTrackTPSimDirectMatchedEta->getNbinsX(),
                                        meBTLTrackTPSimDirectMatchedEta->getTH1()->GetXaxis()->GetXmin(),
                                        meBTLTrackTPSimDirectMatchedEta->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoDirectEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimRecoDirectMatchedEta, meBTLTrackTPSimDirectMatchedEta, meBTLRecoDirectEtaEff_);

  meBTLRecoDirectPtEff_ = ibook.book1D("BTLRecoDirectPtEff",
                                       "Reconstruction Efficiency VS Pt;Pt;Efficiency",
                                       meBTLTrackTPSimDirectMatchedPt->getNbinsX(),
                                       meBTLTrackTPSimDirectMatchedPt->getTH1()->GetXaxis()->GetXmin(),
                                       meBTLTrackTPSimDirectMatchedPt->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoDirectPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimRecoDirectMatchedPt, meBTLTrackTPSimDirectMatchedPt, meBTLRecoDirectPtEff_);

  meBTLRecoOtherEtaEff_ = ibook.book1D("BTLRecoOtherEtaEff",
                                       "Reconstruction Efficiency VS Eta;#eta;Efficiency",
                                       meBTLTrackTPSimOtherMatchedEta->getNbinsX(),
                                       meBTLTrackTPSimOtherMatchedEta->getTH1()->GetXaxis()->GetXmin(),
                                       meBTLTrackTPSimOtherMatchedEta->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoOtherEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimRecoOtherMatchedEta, meBTLTrackTPSimOtherMatchedEta, meBTLRecoOtherEtaEff_);

  meBTLRecoOtherPtEff_ = ibook.book1D("BTLRecoOtherPtEff",
                                      "Reconstruction Efficiency VS Pt;Pt;Efficiency",
                                      meBTLTrackTPSimOtherMatchedPt->getNbinsX(),
                                      meBTLTrackTPSimOtherMatchedPt->getTH1()->GetXaxis()->GetXmin(),
                                      meBTLTrackTPSimOtherMatchedPt->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoOtherPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimRecoOtherMatchedPt, meBTLTrackTPSimOtherMatchedPt, meBTLRecoOtherPtEff_);

  meETLRecoEtaposEff_ = ibook.book1D("ETLRecoEtaposEff",
                                     "Reconstruction Efficiency VS Eta (positive);#eta;Efficiency",
                                     meETLTrackTPSimMatchedEtapos->getNbinsX(),
                                     meETLTrackTPSimMatchedEtapos->getTH1()->GetXaxis()->GetXmin(),
                                     meETLTrackTPSimMatchedEtapos->getTH1()->GetXaxis()->GetXmax());
  meETLRecoEtaposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPSimRecoMatchedEtapos, meETLTrackTPSimMatchedEtapos, meETLRecoEtaposEff_);

  meETL2RecoEtaposEff_ = ibook.book1D("ETL2RecoEtaposEff",
                                      "Two hits Reconstruction Efficiency VS Eta (positive);#eta;Efficiency",
                                      meETL2TrackTPSimMatchedEtapos->getNbinsX(),
                                      meETL2TrackTPSimMatchedEtapos->getTH1()->GetXaxis()->GetXmin(),
                                      meETL2TrackTPSimMatchedEtapos->getTH1()->GetXaxis()->GetXmax());
  meETL2RecoEtaposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackTPSimRecoMatchedEtapos, meETL2TrackTPSimMatchedEtapos, meETL2RecoEtaposEff_);

  meETLRecoPtposEff_ = ibook.book1D("ETLRecoPtposEff",
                                    "Reconstruction Efficiency VS Pt (positive);Pt;Efficiency",
                                    meETLTrackTPSimMatchedPtpos->getNbinsX(),
                                    meETLTrackTPSimMatchedPtpos->getTH1()->GetXaxis()->GetXmin(),
                                    meETLTrackTPSimMatchedPtpos->getTH1()->GetXaxis()->GetXmax());
  meETLRecoPtposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPSimRecoMatchedPtpos, meETLTrackTPSimMatchedPtpos, meETLRecoPtposEff_);

  meETL2RecoPtposEff_ = ibook.book1D("ETL2RecoPtposEff",
                                     "Two hits Reconstruction Efficiency VS Pt (positive);Pt;Efficiency",
                                     meETL2TrackTPSimMatchedPtpos->getNbinsX(),
                                     meETL2TrackTPSimMatchedPtpos->getTH1()->GetXaxis()->GetXmin(),
                                     meETL2TrackTPSimMatchedPtpos->getTH1()->GetXaxis()->GetXmax());
  meETL2RecoPtposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackTPSimRecoMatchedPtpos, meETL2TrackTPSimMatchedPtpos, meETL2RecoPtposEff_);

  // 2/total
  meBTLRecoDirectEtaTotEff_ = ibook.book1D("BTLRecoDirectEtaTotEff",
                                           "Reconstruction Efficiency VS Eta;#eta;Efficiency",
                                           meBTLTrackEffEtaTot->getNbinsX(),
                                           meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                           meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoDirectEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimRecoDirectMatchedEta, meBTLTrackEffEtaTot, meBTLRecoDirectEtaTotEff_);

  meBTLRecoDirectPtTotEff_ = ibook.book1D("BTLRecoDirectPtTotEff",
                                          "Reconstruction Efficiency VS Pt;Pt;Efficiency",
                                          meBTLTrackEffPtTot->getNbinsX(),
                                          meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                          meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoDirectPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimRecoDirectMatchedPt, meBTLTrackEffPtTot, meBTLRecoDirectPtTotEff_);

  meBTLRecoOtherEtaTotEff_ = ibook.book1D("BTLRecoOtherEtaTotEff",
                                          "Reconstruction Efficiency VS Eta;#eta;Efficiency",
                                          meBTLTrackEffEtaTot->getNbinsX(),
                                          meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                          meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoOtherEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimRecoOtherMatchedEta, meBTLTrackEffEtaTot, meBTLRecoOtherEtaTotEff_);

  meBTLRecoOtherPtTotEff_ = ibook.book1D("BTLRecoOtherPtTotEff",
                                         "Reconstruction Efficiency VS Pt;Pt;Efficiency",
                                         meBTLTrackEffPtTot->getNbinsX(),
                                         meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                         meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoOtherPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimRecoOtherMatchedPt, meBTLTrackEffPtTot, meBTLRecoOtherPtTotEff_);

  meETLRecoEtaposTotEff_ = ibook.book1D("ETLRecoEtaposTotEff",
                                        "Reconstruction Efficiency VS Eta (positive);#eta;Efficiency",
                                        meETLTrackEffEtaTotZpos->getNbinsX(),
                                        meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                        meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLRecoEtaposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPSimRecoMatchedEtapos, meETLTrackEffEtaTotZpos, meETLRecoEtaposTotEff_);

  meETL2RecoEtaposTotEff_ = ibook.book1D("ETL2RecoEtaposTotEff",
                                         "Two hits Reconstruction Efficiency VS Eta (positive);#eta;Efficiency",
                                         meETLTrackEffEtaTotZpos->getNbinsX(),
                                         meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                         meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2RecoEtaposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackTPSimRecoMatchedEtapos, meETLTrackEffEtaTotZpos, meETL2RecoEtaposTotEff_);

  meETLRecoPtposTotEff_ = ibook.book1D("ETLRecoPtposTotEff",
                                       "Reconstruction Efficiency VS Pt (positive);Pt;Efficiency",
                                       meETLTrackEffPtTotZpos->getNbinsX(),
                                       meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                       meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLRecoPtposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPSimRecoMatchedPtpos, meETLTrackEffPtTotZpos, meETLRecoPtposTotEff_);

  meETL2RecoPtposTotEff_ = ibook.book1D("ETL2RecoPtposTotEff",
                                        "Two hits Reconstruction Efficiency VS Pt (positive);Pt;Efficiency",
                                        meETLTrackEffPtTotZpos->getNbinsX(),
                                        meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                        meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2RecoPtposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackTPSimRecoMatchedPtpos, meETLTrackEffPtTotZpos, meETL2RecoPtposTotEff_);

  meBTLRecoNoRecoDirectEtaTotEff_ = ibook.book1D("BTLRecoNoRecoDirectEtaTotEff",
                                                 "Efficiency (No Reco) VS Eta;#eta;Efficiency",
                                                 meBTLTrackEffEtaTot->getNbinsX(),
                                                 meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                 meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoNoRecoDirectEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimNoRecoDirectMatchedEta, meBTLTrackEffEtaTot, meBTLRecoNoRecoDirectEtaTotEff_);

  meBTLRecoNoRecoDirectPtTotEff_ = ibook.book1D("BTLRecoNoRecoDirectPtTotEff",
                                                "Efficiency (No Reco) VS Pt;Pt;Efficiency",
                                                meBTLTrackEffPtTot->getNbinsX(),
                                                meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                                meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoNoRecoDirectPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimNoRecoDirectMatchedPt, meBTLTrackEffPtTot, meBTLRecoNoRecoDirectPtTotEff_);

  meBTLRecoNoRecoOtherEtaTotEff_ = ibook.book1D("BTLRecoNoRecoOtherEtaTotEff",
                                                "Efficiency (No Reco) VS Eta;#eta;Efficiency",
                                                meBTLTrackEffEtaTot->getNbinsX(),
                                                meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoNoRecoOtherEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimNoRecoOtherMatchedEta, meBTLTrackEffEtaTot, meBTLRecoNoRecoOtherEtaTotEff_);

  meBTLRecoNoRecoOtherPtTotEff_ = ibook.book1D("BTLRecoNoRecoOtherPtTotEff",
                                               "Efficiency (No Reco) VS Pt;Pt;Efficiency",
                                               meBTLTrackEffPtTot->getNbinsX(),
                                               meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                               meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLRecoNoRecoOtherPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPSimNoRecoOtherMatchedPt, meBTLTrackEffPtTot, meBTLRecoNoRecoOtherPtTotEff_);

  meETLRecoNoRecoEtaposTotEff_ = ibook.book1D("ETLRecoNoRecoEtaposTotEff",
                                              "Efficiency (No Reco) VS Eta (positive);#eta;Efficiency",
                                              meETLTrackEffEtaTotZpos->getNbinsX(),
                                              meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                              meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLRecoNoRecoEtaposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPSimNoRecoMatchedEtapos, meETLTrackEffEtaTotZpos, meETLRecoNoRecoEtaposTotEff_);

  meETL2RecoNoRecoEtaposTotEff_ = ibook.book1D("ETL2RecoNoRecoEtaposTotEff",
                                               "Two hits Efficiency (No Reco) VS Eta (positive);#eta;Efficiency",
                                               meETLTrackEffEtaTotZpos->getNbinsX(),
                                               meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                               meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2RecoNoRecoEtaposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackTPSimNoRecoMatchedEtapos, meETLTrackEffEtaTotZpos, meETL2RecoNoRecoEtaposTotEff_);

  meETLRecoNoRecoPtposTotEff_ = ibook.book1D("ETLRecoNoRecoPtposTotEff",
                                             "Efficiency (No Reco) VS Pt (positive);Pt;Efficiency",
                                             meETLTrackEffPtTotZpos->getNbinsX(),
                                             meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                             meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLRecoNoRecoPtposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPSimNoRecoMatchedPtpos, meETLTrackEffPtTotZpos, meETLRecoNoRecoPtposTotEff_);

  meETL2RecoNoRecoPtposTotEff_ = ibook.book1D("ETL2RecoNoRecoPtposTotEff",
                                              "Two hits Efficiency (No Reco) VS Pt (positive);Pt;Efficiency",
                                              meETLTrackEffPtTotZpos->getNbinsX(),
                                              meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                              meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2RecoNoRecoPtposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackTPSimNoRecoMatchedPtpos, meETLTrackEffPtTotZpos, meETL2RecoNoRecoPtposTotEff_);

  // 3/2
  meBTLCorrectSimDirectEtaEff_ = ibook.book1D("BTLCorrectSimDirectEtaEff",
                                              "SimToSim Efficiency VS Eta;#eta;Efficiency",
                                              meBTLTrackTPSimRecoDirectMatchedEta->getNbinsX(),
                                              meBTLTrackTPSimRecoDirectMatchedEta->getTH1()->GetXaxis()->GetXmin(),
                                              meBTLTrackTPSimRecoDirectMatchedEta->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimDirectEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(
      meBTLTrackSimToSimDirectMatchedEta, meBTLTrackTPSimRecoDirectMatchedEta, meBTLCorrectSimDirectEtaEff_);

  meBTLCorrectSimDirectPtEff_ = ibook.book1D("BTLCorrectSimDirectPtEff",
                                             "SimToSim Efficiency VS Pt;Pt;Efficiency",
                                             meBTLTrackTPSimRecoDirectMatchedPt->getNbinsX(),
                                             meBTLTrackTPSimRecoDirectMatchedPt->getTH1()->GetXaxis()->GetXmin(),
                                             meBTLTrackTPSimRecoDirectMatchedPt->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimDirectPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(
      meBTLTrackSimToSimDirectMatchedPt, meBTLTrackTPSimRecoDirectMatchedPt, meBTLCorrectSimDirectPtEff_);

  meBTLCorrectSimOtherEtaEff_ = ibook.book1D("BTLCorrectSimOtherEtaEff",
                                             "SimToSim Efficiency VS Eta;#eta;Efficiency",
                                             meBTLTrackTPSimRecoOtherMatchedEta->getNbinsX(),
                                             meBTLTrackTPSimRecoOtherMatchedEta->getTH1()->GetXaxis()->GetXmin(),
                                             meBTLTrackTPSimRecoOtherMatchedEta->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimOtherEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(
      meBTLTrackSimToSimOtherMatchedEta, meBTLTrackTPSimRecoOtherMatchedEta, meBTLCorrectSimOtherEtaEff_);

  meBTLCorrectSimOtherPtEff_ = ibook.book1D("BTLCorrectSimOtherPtEff",
                                            "SimToSim Efficiency VS Pt;Pt;Efficiency",
                                            meBTLTrackTPSimRecoOtherMatchedPt->getNbinsX(),
                                            meBTLTrackTPSimRecoOtherMatchedPt->getTH1()->GetXaxis()->GetXmin(),
                                            meBTLTrackTPSimRecoOtherMatchedPt->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimOtherPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimOtherMatchedPt, meBTLTrackTPSimRecoOtherMatchedPt, meBTLCorrectSimOtherPtEff_);

  meETLCorrectSimEtaposEff_ = ibook.book1D("ETLCorrectSimEtaposEff",
                                           "SimToSim Efficiency VS Eta (positive);#eta;Efficiency",
                                           meETLTrackTPSimRecoMatchedEtapos->getNbinsX(),
                                           meETLTrackTPSimRecoMatchedEtapos->getTH1()->GetXaxis()->GetXmin(),
                                           meETLTrackTPSimRecoMatchedEtapos->getTH1()->GetXaxis()->GetXmax());
  meETLCorrectSimEtaposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackSimToSimMatchedEtapos, meETLTrackTPSimRecoMatchedEtapos, meETLCorrectSimEtaposEff_);

  meETL2CorrectSimEtaposEff_ = ibook.book1D("ETL2CorrectSimEtaposEff",
                                            "Two hits SimToSim Efficiency VS Eta (positive);#eta;Efficiency",
                                            meETL2TrackTPSimRecoMatchedEtapos->getNbinsX(),
                                            meETL2TrackTPSimRecoMatchedEtapos->getTH1()->GetXaxis()->GetXmin(),
                                            meETL2TrackTPSimRecoMatchedEtapos->getTH1()->GetXaxis()->GetXmax());
  meETL2CorrectSimEtaposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackSimToSimMatchedEtapos, meETL2TrackTPSimRecoMatchedEtapos, meETL2CorrectSimEtaposEff_);

  meETLCorrectSimPtposEff_ = ibook.book1D("ETLCorrectSimPtposEff",
                                          "SimToSim Efficiency VS Pt (positive);Pt;Efficiency",
                                          meETLTrackTPSimRecoMatchedPtpos->getNbinsX(),
                                          meETLTrackTPSimRecoMatchedPtpos->getTH1()->GetXaxis()->GetXmin(),
                                          meETLTrackTPSimRecoMatchedPtpos->getTH1()->GetXaxis()->GetXmax());
  meETLCorrectSimPtposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackSimToSimMatchedPtpos, meETLTrackTPSimRecoMatchedPtpos, meETLCorrectSimPtposEff_);

  meETL2CorrectSimPtposEff_ = ibook.book1D("ETL2CorrectSimPtposEff",
                                           "Two hits SimToSim Efficiency VS Pt (positive);Pt;Efficiency",
                                           meETL2TrackTPSimRecoMatchedPtpos->getNbinsX(),
                                           meETL2TrackTPSimRecoMatchedPtpos->getTH1()->GetXaxis()->GetXmin(),
                                           meETL2TrackTPSimRecoMatchedPtpos->getTH1()->GetXaxis()->GetXmax());
  meETL2CorrectSimPtposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackSimToSimMatchedPtpos, meETL2TrackTPSimRecoMatchedPtpos, meETL2CorrectSimPtposEff_);

  // 3/Total

  meBTLCorrectSimDirectEtaTotEff_ = ibook.book1D("BTLCorrectSimDirectEtaTotEff",
                                                 "SimToSim Efficiency VS Eta;#eta;Efficiency",
                                                 meBTLTrackEffEtaTot->getNbinsX(),
                                                 meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                 meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimDirectEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimDirectMatchedEta, meBTLTrackEffEtaTot, meBTLCorrectSimDirectEtaTotEff_);

  meBTLCorrectSimDirectPtTotEff_ = ibook.book1D("BTLCorrectSimDirectPtTotEff",
                                                "SimToSim Efficiency VS Pt;Pt;Efficiency",
                                                meBTLTrackEffPtTot->getNbinsX(),
                                                meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                                meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimDirectPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimDirectMatchedPt, meBTLTrackEffPtTot, meBTLCorrectSimDirectPtTotEff_);

  meBTLCorrectSimOtherEtaTotEff_ = ibook.book1D("BTLCorrectSimOtherEtaTotEff",
                                                "SimToSim Efficiency VS Eta;#eta;Efficiency",
                                                meBTLTrackEffEtaTot->getNbinsX(),
                                                meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimOtherEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimOtherMatchedEta, meBTLTrackEffEtaTot, meBTLCorrectSimOtherEtaTotEff_);

  meBTLCorrectSimOtherPtTotEff_ = ibook.book1D("BTLCorrectSimOtherPtTotEff",
                                               "SimToSim Efficiency VS Pt;Pt;Efficiency",
                                               meBTLTrackEffPtTot->getNbinsX(),
                                               meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                               meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimOtherPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimOtherMatchedPt, meBTLTrackEffPtTot, meBTLCorrectSimOtherPtTotEff_);

  meETLCorrectSimEtaposTotEff_ = ibook.book1D("ETLCorrectSimEtaposTotEff",
                                              "SimToSim Efficiency VS Eta (positive);#eta;Efficiency",
                                              meETLTrackEffEtaTotZpos->getNbinsX(),
                                              meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                              meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLCorrectSimEtaposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackSimToSimMatchedEtapos, meETLTrackEffEtaTotZpos, meETLCorrectSimEtaposTotEff_);

  meETL2CorrectSimEtaposTotEff_ = ibook.book1D("ETL2CorrectSimEtaposTotEff",
                                               "Two hits SimToSim Efficiency VS Eta (positive);#eta;Efficiency",
                                               meETLTrackEffEtaTotZpos->getNbinsX(),
                                               meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                               meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2CorrectSimEtaposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackSimToSimMatchedEtapos, meETLTrackEffEtaTotZpos, meETL2CorrectSimEtaposTotEff_);

  meETLCorrectSimPtposTotEff_ = ibook.book1D("ETLCorrectSimPtposTotEff",
                                             "SimToSim Efficiency VS Pt (positive);Pt;Efficiency",
                                             meETLTrackEffPtTotZpos->getNbinsX(),
                                             meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                             meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLCorrectSimPtposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackSimToSimMatchedPtpos, meETLTrackEffPtTotZpos, meETLCorrectSimPtposTotEff_);

  meETL2CorrectSimPtposTotEff_ = ibook.book1D("ETL2CorrectSimPtposTotEff",
                                              "Two hits SimToSim Efficiency VS Pt (positive);Pt;Efficiency",
                                              meETLTrackEffPtTotZpos->getNbinsX(),
                                              meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                              meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2CorrectSimPtposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackSimToSimMatchedPtpos, meETLTrackEffPtTotZpos, meETL2CorrectSimPtposTotEff_);

  // 3/ Reco

  meBTLCorrectSimDirectEtaMtdEff_ = ibook.book1D("BTLCorrectSimDirectEtaMtdEff",
                                                 "SimToSim Efficiency VS Eta;#eta;Efficiency",
                                                 meBTLTrackEffEtaMtd->getNbinsX(),
                                                 meBTLTrackEffEtaMtd->getTH1()->GetXaxis()->GetXmin(),
                                                 meBTLTrackEffEtaMtd->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimDirectEtaMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimDirectMatchedEta, meBTLTrackEffEtaMtd, meBTLCorrectSimDirectEtaMtdEff_);

  meBTLCorrectSimDirectPtMtdEff_ = ibook.book1D("BTLCorrectSimDirectPtMtdEff",
                                                "SimToSim Efficiency VS Pt;Pt;Efficiency",
                                                meBTLTrackEffPtMtd->getNbinsX(),
                                                meBTLTrackEffPtMtd->getTH1()->GetXaxis()->GetXmin(),
                                                meBTLTrackEffPtMtd->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimDirectPtMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimDirectMatchedPt, meBTLTrackEffPtMtd, meBTLCorrectSimDirectPtMtdEff_);

  meBTLCorrectSimOtherEtaMtdEff_ = ibook.book1D("BTLCorrectSimOtherEtaMtdEff",
                                                "SimToSim Efficiency VS Eta;#eta;Efficiency",
                                                meBTLTrackEffEtaMtd->getNbinsX(),
                                                meBTLTrackEffEtaMtd->getTH1()->GetXaxis()->GetXmin(),
                                                meBTLTrackEffEtaMtd->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimOtherEtaMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimOtherMatchedEta, meBTLTrackEffEtaMtd, meBTLCorrectSimOtherEtaMtdEff_);

  meBTLCorrectSimOtherPtMtdEff_ = ibook.book1D("BTLCorrectSimOtherPtMtdEff",
                                               "SimToSim Efficiency VS Pt;Pt;Efficiency",
                                               meBTLTrackEffPtMtd->getNbinsX(),
                                               meBTLTrackEffPtMtd->getTH1()->GetXaxis()->GetXmin(),
                                               meBTLTrackEffPtMtd->getTH1()->GetXaxis()->GetXmax());
  meBTLCorrectSimOtherPtMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackSimToSimOtherMatchedPt, meBTLTrackEffPtMtd, meBTLCorrectSimOtherPtMtdEff_);

  meETLCorrectSimEtaposMtdEff_ = ibook.book1D("ETLCorrectSimEtaposMtdEff",
                                              "SimToSim Efficiency VS Eta (positive);#eta;Efficiency",
                                              meETLTrackEffEtaMtdZpos->getNbinsX(),
                                              meETLTrackEffEtaMtdZpos->getTH1()->GetXaxis()->GetXmin(),
                                              meETLTrackEffEtaMtdZpos->getTH1()->GetXaxis()->GetXmax());
  meETLCorrectSimEtaposMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackSimToSimMatchedEtapos, meETLTrackEffEtaMtdZpos, meETLCorrectSimEtaposMtdEff_);

  meETL2CorrectSimEtaposMtdEff_ = ibook.book1D("ETL2CorrectSimEtaposMtdEff",
                                               "Two hits SimToSim Efficiency VS Eta (positive);#eta;Efficiency",
                                               meETLTrackEffEta2MtdZpos->getNbinsX(),
                                               meETLTrackEffEta2MtdZpos->getTH1()->GetXaxis()->GetXmin(),
                                               meETLTrackEffEta2MtdZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2CorrectSimEtaposMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackSimToSimMatchedEtapos, meETLTrackEffEta2MtdZpos, meETL2CorrectSimEtaposMtdEff_);

  meETLCorrectSimPtposMtdEff_ = ibook.book1D("ETLCorrectSimPtposMtdEff",
                                             "SimToSim Efficiency VS Pt (positive);Pt;Efficiency",
                                             meETLTrackEffPtMtdZpos->getNbinsX(),
                                             meETLTrackEffPtMtdZpos->getTH1()->GetXaxis()->GetXmin(),
                                             meETLTrackEffPtMtdZpos->getTH1()->GetXaxis()->GetXmax());
  meETLCorrectSimPtposMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackSimToSimMatchedPtpos, meETLTrackEffPtMtdZpos, meETLCorrectSimPtposMtdEff_);

  meETL2CorrectSimPtposMtdEff_ = ibook.book1D("ETL2CorrectSimPtposMtdEff",
                                              "Two hits SimToSim Efficiency VS Pt (positive);Pt;Efficiency",
                                              meETLTrackEffPt2MtdZpos->getNbinsX(),
                                              meETLTrackEffPt2MtdZpos->getTH1()->GetXaxis()->GetXmin(),
                                              meETLTrackEffPt2MtdZpos->getTH1()->GetXaxis()->GetXmax());
  meETL2CorrectSimPtposMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETL2TrackSimToSimMatchedPtpos, meETLTrackEffPt2MtdZpos, meETL2CorrectSimPtposMtdEff_);

  // 5/ Reco

  meBTLIncorrectSimEtaMtdEff_ = ibook.book1D("BTLIncorrectSimEtaMtdEff",
                                             "NoSimToSim Efficiency VS Eta;#eta;Efficiency",
                                             meBTLTrackEffEtaMtd->getNbinsX(),
                                             meBTLTrackEffEtaMtd->getTH1()->GetXaxis()->GetXmin(),
                                             meBTLTrackEffEtaMtd->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimEtaMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimMatchedEta, meBTLTrackEffEtaMtd, meBTLIncorrectSimEtaMtdEff_);

  meBTLIncorrectSimPtMtdEff_ = ibook.book1D("BTLIncorrectSimPtMtdEff",
                                            "NoSimToSim Efficiency VS Pt;Pt;Efficiency",
                                            meBTLTrackEffPtMtd->getNbinsX(),
                                            meBTLTrackEffPtMtd->getTH1()->GetXaxis()->GetXmin(),
                                            meBTLTrackEffPtMtd->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimPtMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimMatchedPt, meBTLTrackEffPtMtd, meBTLIncorrectSimPtMtdEff_);

  meETLIncorrectSimEtaposMtdEff_ = ibook.book1D("ETLIncorrectSimEtaposMtdEff",
                                                "NoSimToSim Efficiency VS Eta (positive);#eta;Efficiency",
                                                meETLTrackEffEtaMtdZpos->getNbinsX(),
                                                meETLTrackEffEtaMtdZpos->getTH1()->GetXaxis()->GetXmin(),
                                                meETLTrackEffEtaMtdZpos->getTH1()->GetXaxis()->GetXmax());
  meETLIncorrectSimEtaposMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackNoSimToSimMatchedEtapos, meETLTrackEffEtaMtdZpos, meETLIncorrectSimEtaposMtdEff_);

  meETLIncorrectSimPtposMtdEff_ = ibook.book1D("ETLIncorrectSimPtposMtdEff",
                                               "NoSimToSim Efficiency VS Pt (positive);Pt;Efficiency",
                                               meETLTrackEffPtMtdZpos->getNbinsX(),
                                               meETLTrackEffPtMtdZpos->getTH1()->GetXaxis()->GetXmin(),
                                               meETLTrackEffPtMtdZpos->getTH1()->GetXaxis()->GetXmax());
  meETLIncorrectSimPtposMtdEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackNoSimToSimMatchedPtpos, meETLTrackEffPtMtdZpos, meETLIncorrectSimPtposMtdEff_);

  meBTLIncorrectSimEtaTotEff_ = ibook.book1D("BTLIncorrectSimEtaTotEff",
                                             "NoSimToSim Efficiency VS Eta;#eta;Efficiency",
                                             meBTLTrackEffEtaTot->getNbinsX(),
                                             meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                             meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimMatchedEta, meBTLTrackEffEtaTot, meBTLIncorrectSimEtaTotEff_);

  meBTLIncorrectSimDirectEtaTotEff_ = ibook.book1D("BTLIncorrectSimDirectEtaTotEff",
                                                   "NoSimToSim Efficiency VS Eta;#eta;Efficiency",
                                                   meBTLTrackEffEtaTot->getNbinsX(),
                                                   meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                   meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimDirectEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimDirectMatchedEta, meBTLTrackEffEtaTot, meBTLIncorrectSimDirectEtaTotEff_);

  meBTLIncorrectSimSecondEtaTotEff_ = ibook.book1D("BTLIncorrectSimSecondEtaTotEff",
                                                   "NoSimToSim Efficiency VS Eta;#eta;Efficiency",
                                                   meBTLTrackEffEtaTot->getNbinsX(),
                                                   meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                   meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimSecondEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimSecondMatchedEta, meBTLTrackEffEtaTot, meBTLIncorrectSimSecondEtaTotEff_);

  meBTLIncorrectSimBackEtaTotEff_ = ibook.book1D("BTLIncorrectSimBackEtaTotEff",
                                                 "NoSimToSim Efficiency VS Eta;#eta;Efficiency",
                                                 meBTLTrackEffEtaTot->getNbinsX(),
                                                 meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                 meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimBackEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimBackMatchedEta, meBTLTrackEffEtaTot, meBTLIncorrectSimBackEtaTotEff_);

  meBTLIncorrectSimLoopEtaTotEff_ = ibook.book1D("BTLIncorrectSimLoopEtaTotEff",
                                                 "NoSimToSim Efficiency VS Eta;#eta;Efficiency",
                                                 meBTLTrackEffEtaTot->getNbinsX(),
                                                 meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                                 meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimLoopEtaTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimLoopMatchedEta, meBTLTrackEffEtaTot, meBTLIncorrectSimLoopEtaTotEff_);

  meBTLIncorrectSimPtTotEff_ = ibook.book1D("BTLIncorrectSimPtTotEff",
                                            "NoSimToSim Efficiency VS Pt;Pt;Efficiency",
                                            meBTLTrackEffPtTot->getNbinsX(),
                                            meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                            meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimMatchedPt, meBTLTrackEffPtTot, meBTLIncorrectSimPtTotEff_);

  meBTLIncorrectSimDirectPtTotEff_ = ibook.book1D("BTLIncorrectSimDirectPtTotEff",
                                                  "NoSimToSim Efficiency VS Pt;Pt;Efficiency",
                                                  meBTLTrackEffPtTot->getNbinsX(),
                                                  meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                                  meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimDirectPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimDirectMatchedPt, meBTLTrackEffPtTot, meBTLIncorrectSimDirectPtTotEff_);

  meBTLIncorrectSimSecondPtTotEff_ = ibook.book1D("BTLIncorrectSimSecondPtTotEff",
                                                  "NoSimToSim Efficiency VS Pt;Pt;Efficiency",
                                                  meBTLTrackEffPtTot->getNbinsX(),
                                                  meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                                  meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimSecondPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimSecondMatchedPt, meBTLTrackEffPtTot, meBTLIncorrectSimSecondPtTotEff_);

  meBTLIncorrectSimBackPtTotEff_ = ibook.book1D("BTLIncorrectSimBackPtTotEff",
                                                "NoSimToSim Efficiency VS Pt;Pt;Efficiency",
                                                meBTLTrackEffPtTot->getNbinsX(),
                                                meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                                meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimBackPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimBackMatchedPt, meBTLTrackEffPtTot, meBTLIncorrectSimBackPtTotEff_);

  meBTLIncorrectSimLoopPtTotEff_ = ibook.book1D("BTLIncorrectSimLoopPtTotEff",
                                                "NoSimToSim Efficiency VS Pt;Pt;Efficiency",
                                                meBTLTrackEffPtTot->getNbinsX(),
                                                meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                                meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLIncorrectSimLoopPtTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackNoSimToSimLoopMatchedPt, meBTLTrackEffPtTot, meBTLIncorrectSimLoopPtTotEff_);

  meETLIncorrectSimEtaposTotEff_ = ibook.book1D("ETLIncorrectSimEtaposTotEff",
                                                "NoSimToSim Efficiency VS Eta (positive);#eta;Efficiency",
                                                meETLTrackEffEtaTotZpos->getNbinsX(),
                                                meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                                meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLIncorrectSimEtaposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackNoSimToSimMatchedEtapos, meETLTrackEffEtaTotZpos, meETLIncorrectSimEtaposTotEff_);

  meETLIncorrectSimPtposTotEff_ = ibook.book1D("ETLIncorrectSimPtposTotEff",
                                               "NoSimToSim Efficiency VS Pt (positive);Pt;Efficiency",
                                               meETLTrackEffPtTotZpos->getNbinsX(),
                                               meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                               meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLIncorrectSimPtposTotEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackNoSimToSimMatchedPtpos, meETLTrackEffPtTotZpos, meETLIncorrectSimPtposTotEff_);

  meBTLTrackTPNoSimEtaEff_ = ibook.book1D("BTLTrackTPNoSimEtaEff",
                                          "Track Efficiency VS Eta;#eta;Efficiency",
                                          meBTLTrackEffEtaTot->getNbinsX(),
                                          meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmin(),
                                          meBTLTrackEffEtaTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTrackTPNoSimEtaEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPNoSimEta, meBTLTrackEffEtaTot, meBTLTrackTPNoSimEtaEff_);

  meBTLTrackTPNoSimPtEff_ = ibook.book1D("BTLTrackTPNoSimPtEff",
                                         "Track Efficiency VS Pt;Pt;Efficiency",
                                         meBTLTrackEffPtTot->getNbinsX(),
                                         meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmin(),
                                         meBTLTrackEffPtTot->getTH1()->GetXaxis()->GetXmax());
  meBTLTrackTPNoSimPtEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meBTLTrackTPNoSimPt, meBTLTrackEffPtTot, meBTLTrackTPNoSimPtEff_);

  meETLTrackTPNoSimEtaposEff_ = ibook.book1D("ETLTrackTPNoSimEtaposEff",
                                             "Track Efficiency VS Eta;#eta;Efficiency",
                                             meETLTrackEffEtaTotZpos->getNbinsX(),
                                             meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                             meETLTrackEffEtaTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLTrackTPNoSimEtaposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPNoSimEtapos, meETLTrackEffEtaTotZpos, meETLTrackTPNoSimEtaposEff_);

  meETLTrackTPNoSimPtposEff_ = ibook.book1D("ETLTrackTPNoSimPtposEff",
                                            "Track Efficiency VS Pt;#pt;Efficiency",
                                            meETLTrackEffPtTotZpos->getNbinsX(),
                                            meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmin(),
                                            meETLTrackEffPtTotZpos->getTH1()->GetXaxis()->GetXmax());
  meETLTrackTPNoSimPtposEff_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meETLTrackTPNoSimPtpos, meETLTrackEffPtTotZpos, meETLTrackTPNoSimPtposEff_);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ----------
void MtdTracksHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks/");

  descriptions.add("MtdTracksPostProcessor", desc);
}

DEFINE_FWK_MODULE(MtdTracksHarvester);
