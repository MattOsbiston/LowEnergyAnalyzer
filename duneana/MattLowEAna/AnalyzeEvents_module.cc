////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEvents_module.cc
//
// Generated at Sat Apr 20 06:39:03 2024 by Matthew Osbiston using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include <TTree.h>
#include "art/Utilities/ToolMacros.h"

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <tuple>

namespace lowe {
  class AnalyzeEvents;
}


class lowe::AnalyzeEvents : public art::EDAnalyzer {
public:
  explicit AnalyzeEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEvents(AnalyzeEvents const&) = delete;
  AnalyzeEvents(AnalyzeEvents&&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents const&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reset(bool deepClean=false);

private:

  // Declare member data here.
  double CalculateEnergy(const detinfo::DetectorClocksData& clockData,
                           const detinfo::DetectorPropertiesData& detProp,
                           const std::vector<art::Ptr<recob::Hit>>& hits,
                           const geo::PlaneID::PlaneID_t plane) const;

  TTree *fTree;
  unsigned int fEventID;
  int fBestView;
  unsigned int fBestViewHitCount;
  unsigned int fBestUViewHitCount;
  unsigned int fBestVViewHitCount;
  unsigned int fBestWViewHitCount;
  double fBestUViewEnergy;
  double fBestVViewEnergy;
  double fBestWViewEnergy;
  double fPFPVertexX;
  double fPFPVertexY;
  double fPFPVertexZ;
  unsigned int fNMCParticles;
  unsigned int fNPFPs;

  static const int kNMaxMCParticles = 200000;

  int fMCParticleIsPrimary[kNMaxMCParticles];
  double fMCParticleTrueEnergy[kNMaxMCParticles];
  double fMCNeutrinoTrueEnergy;
  double fMCNeutronEnergy[5];
  double fMCElectronTrueEnergy;
  double fMCNeutrinoTrueVertexX;
  double fMCNeutrinoTrueVertexY;
  double fMCNeutrinoTrueVertexZ;
  int fMCParticlePdgCode[kNMaxMCParticles];
  int fMCNeutrinoPdgCode;
  
  art::ServiceHandle<geo::Geometry> fGeometry;

  std::string fLabels;
  double fRecombinationFactor;
  std::string fHitLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fTruthLabel;
};


lowe::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fLabels(p.get<std::string>("MarleyLabel"))
  , fRecombinationFactor(p.get<double>("RecombinationFactor"))
  , fHitLabel(p.get<std::string>("HitLabel"))
  , fPFParticleLabel(p.get<std::string>("ParticleLabel"))
  , fTrackLabel(p.get<std::string>("TrackLabel"))
  , fTruthLabel(p.get<std::string>("TruthLabel"))
{
}

void lowe::AnalyzeEvents::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  reset(true);
  fEventID = e.id().event();
 
  std::cout << "*--------------------------- Running LowE Analysis --------------------------------------*" << std::endl; 
  //get all hits
  art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
  std::vector<art::Ptr<recob::PFParticle>> pfpVector;
  std::vector<art::Ptr<recob::Hit>> allHits;
  
  if (pfpHandle.isValid())
    art::fill_ptr_vector(pfpVector, pfpHandle);
  art::FindManyP<recob::Track> pfpTrackAssociation(pfpHandle, e, fTrackLabel);
  float Max(-std::numeric_limits<float>::max());
  fNPFPs = pfpVector.size();
  std::cout << "Number of reconstructed particles:  " << fNPFPs << std::endl;

  for (auto const &pfp : pfpVector)
  { 
    std::vector<art::Ptr<recob::Track>> pfpTracks(pfpTrackAssociation.at(pfp.key()));

    art::FindManyP<recob::Hit> trackHitAssociation(trackHandle, e, fHitLabel);

    for (const art::Ptr<recob::Track> &track : pfpTracks)
    {
	float size = trackHitAssociation.at(track.key()).size();
        std::cout << "This pfp has size:  " << size << std::endl;
	if (size > Max)
	{
            std::cout << "Added hit list length: " << size << " | Max was set to: " << Max << std::endl;
            Max = size;
	    fPFPVertexX = track->Vertex().X();
	    fPFPVertexY = track->Vertex().Y();
	    fPFPVertexZ = track->Vertex().Z();
            allHits.clear();
	    allHits.insert(allHits.end(), trackHitAssociation.at(track.key()).begin(), trackHitAssociation.at(track.key()).end());
	}

    }
  }

  //get mc information
  if (!e.isRealData())
  {
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    if (mcParticles.isValid())
    {
      fNMCParticles = mcParticles->size();
      fMCElectronTrueEnergy=0;
      int index=0;
      for (unsigned int i = 0 ; i < fNMCParticles ; i++)
      {
        const simb::MCParticle trueParticle = mcParticles->at(i);
        fMCParticleTrueEnergy[i]=trueParticle.E();
        fMCParticlePdgCode[i]=trueParticle.PdgCode();
        fMCParticleIsPrimary[i] = trueParticle.Process() == "primary" ? 1 : 0;
        if ((trueParticle.PdgCode() == 2112) && (index < 5))
	{
            fMCNeutronEnergy[index]=trueParticle.E();
	    std::cout << "True Neutron Energy: " << trueParticle.E() * 1000 << std::endl;
            ++index;
	}
	if ((trueParticle.PdgCode() == 11) && (trueParticle.Process() == "primary"))
	    fMCElectronTrueEnergy=trueParticle.E();	
      }
    std::cout << "True Electron Energy: " << fMCElectronTrueEnergy * 1000 << std::endl;
    std::cout << "Completed MCParticle Loop" << std::endl;
    }
    
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    e.getByLabel(fLabels,ThisHandle);
    
    if (ThisHandle)
    {
      auto Marley = e.getValidHandle<std::vector<simb::MCTruth>>(fLabels); //Get MCTruth from MARLEY
      for (const auto &MarleyTruth : *Marley)
      { 
        const simb::MCNeutrino &nue = MarleyTruth.GetNeutrino();
        fMCNeutrinoTrueEnergy=nue.Nu().E();
	fMCNeutrinoTrueVertexX=nue.Nu().Vx();
	fMCNeutrinoTrueVertexY=nue.Nu().Vy();
	fMCNeutrinoTrueVertexZ=nue.Nu().Vz();
        fMCNeutrinoPdgCode=nue.Nu().PdgCode();
	std::cout << "True Neutrino Energy: " << fMCNeutrinoTrueEnergy * 1000 << std::endl;
      }
    }
  }

  
  std::cout << " Size of all hits: " << allHits.size() << std::endl;
  std::map<int, std::vector<art::Ptr<recob::Hit>>> viewHits;

  //loop to get hits in plane
  for (auto const& hit : allHits)
  {
    const int view(hit->View());
    viewHits[view].emplace_back(hit);
  }
  
  //calculate energy in each plane
  fBestView = std::numeric_limits<int>::max();
  fBestViewHitCount = 0;
  fBestUViewHitCount = 0;
  fBestUViewEnergy = 0;
  fBestVViewHitCount = 0;
  fBestVViewEnergy = 0;
  fBestWViewHitCount = 0;
  fBestWViewEnergy = 0; 

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);
 
  std::vector<double> energyVector(fGeometry->Nplanes(), -999.);

  //check best plane
  for (auto const& [view, hits] : viewHits)
  {
    unsigned int viewHitCount = hits.size();
    double energy = AnalyzeEvents::CalculateEnergy(clockData, detProp, hits, view);
   
    std::cout << " Looking at view: " << view << " | Number of hits in view: " << viewHitCount << std::endl;
    if (energy > 0)
      energyVector.at(view) = energy;
  
    if (viewHitCount > fBestViewHitCount)
    {
      fBestView = view;
      fBestViewHitCount = viewHitCount;
    }
    if (view == 0)
    {
      fBestUViewHitCount = viewHitCount;
      fBestUViewEnergy = energy;
    }
    if (view == 1)
    {
      fBestVViewHitCount = viewHitCount;
      fBestVViewEnergy = energy;
    }
    if (view == 2)
    {
      fBestWViewHitCount = viewHitCount;
      fBestWViewEnergy = energy;
    }

  }

  //Fill tree
  fTree->Fill();
}

void lowe::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  // Get TFileService to output a TTree
  // Deep clean variables
  reset(true);
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");

  //Add branches to TTree
  fTree->Branch("EventID", &fEventID, "EventID/i");
  fTree->Branch("BestView", &fBestView, "BestView/i");
  fTree->Branch("BestUViewHitCount", &fBestUViewHitCount, "BestUViewHitCount/i");
  fTree->Branch("BestUViewEnergy", &fBestUViewEnergy, "BestUViewEnergy/d");
  fTree->Branch("BestVViewHitCount", &fBestVViewHitCount, "BestVViewHitCount/i");
  fTree->Branch("BestVViewEnergy", &fBestVViewEnergy, "BestVViewEnergy/d");
  fTree->Branch("BestWViewHitCount", &fBestWViewHitCount, "BestWViewHitCount/i");
  fTree->Branch("BestWViewEnergy", &fBestWViewEnergy, "BestWViewEnergy/d");
  fTree->Branch("nMCParticles", &fNMCParticles, "nMCParticles/i");
  fTree->Branch("nPfPs", &fNPFPs, "nPFPs/i");
  fTree->Branch("PFPVertexX", &fPFPVertexX, "PFPVertexX/d");
  fTree->Branch("PFPVertexY", &fPFPVertexY, "PFPVertexY/d");
  fTree->Branch("PFPVertexZ", &fPFPVertexZ, "PFPVertexZ/d");
  fTree->Branch("MCNeutrinoTrueVertexX", &fMCNeutrinoTrueVertexX, "MCNeutrinoTrueVertexX/d");
  fTree->Branch("MCNeutrinoTrueVertexY", &fMCNeutrinoTrueVertexY, "MCNeutrinoTrueVertexY/d");
  fTree->Branch("MCNeutrinoTrueVertexZ", &fMCNeutrinoTrueVertexZ, "MCNeutrinoTrueVertexZ/d");
  fTree->Branch("mcParticleTrueEnergy", &fMCParticleTrueEnergy, "fMCParticleTrueEnergy[nMCParticles]/D");
  fTree->Branch("mcNeutrinoTrueEnergy", &fMCNeutrinoTrueEnergy, "fMCNeutrinoTrueEnergy/d");
  fTree->Branch("mcNeutronEnergy", &fMCNeutronEnergy, "fMCNeutronEnergy[5]/D");
  fTree->Branch("mcElectronTrueEnergy", &fMCElectronTrueEnergy, "fMCElectronTrueEnergy/d");
  fTree->Branch("mcNeutrinoPdgCode", &fMCNeutrinoPdgCode, "fMCNeutrinoPdgCode/i");
  fTree->Branch("mcParticlePdgCode", &fMCParticlePdgCode, "fMCParticlePdgCode[nMCParticles]/I");
  fTree->Branch("mcParticleIsPrimary", &fMCParticleIsPrimary, "fMCParticleIsPrimary[nMCParticles]/I");

}

void lowe::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
}




double lowe::AnalyzeEvents::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
                                                   const detinfo::DetectorPropertiesData& detProp,
                                                   const std::vector<art::Ptr<recob::Hit>>& hits,
                                                   const geo::PlaneID::PlaneID_t plane) const
  {
    double totalCharge = 0;
    double totalEnergy = 0;
    double correctedTotalCharge = 0;
    double nElectrons = 0;
    std::vector<double> CalAreaConstants = {5.346e-3, 5.339e-3, 5.292e-3};       

    for (auto const& hit : hits) 
    {
      // obtain charge and correct for lifetime
      double const T0 = 0;
      double const t = hit->PeakTime() - trigger_offset(clockData);
      double const timetick = sampling_rate(clockData) * 1.e-3;
      double const adjustedtime = (t * timetick) - (T0 * 1.e-3);
      double const tau = detProp.ElectronLifetime();
      double const correction = exp(adjustedtime / tau);
      std::cout << "hit initial : " << hit->Integral() << " | hit correction : " << correction << std::endl;
      totalCharge +=  hit->Integral() * correction;
    }
     
    std::cout << " *** Total Charge : " << totalCharge << std::endl;
    
    // correct charge due to recombination
    correctedTotalCharge = totalCharge / fRecombinationFactor;
    
    // calculate # of electrons and the corresponding energy
    nElectrons = correctedTotalCharge / CalAreaConstants[plane];
    
    // energy in MeV
    totalEnergy = (nElectrons / util::kGeVToElectrons) * 1000;
    
    return totalEnergy;
  }

  void lowe::AnalyzeEvents::reset(bool deepClean)
  {
    for ( unsigned i=0 ; i < (deepClean ? kNMaxMCParticles : fNMCParticles); i++)
    {
      fMCParticleIsPrimary[i]=0;
      fMCParticlePdgCode[i]=0;
      fMCParticleTrueEnergy[i]=-9999999;
    }
    for ( unsigned i=0 ; i < 5; i++)
    {
      fMCNeutronEnergy[i]=0;
    }
    fMCNeutrinoTrueEnergy=-9999999;
    fMCElectronTrueEnergy=-9999999;
    fMCNeutrinoPdgCode=0;
    fNMCParticles=0;
    fNPFPs=0;
    fPFPVertexX=0;
    fPFPVertexY=0;
    fPFPVertexZ=0;
    fMCNeutrinoTrueVertexX=0;
    fMCNeutrinoTrueVertexY=0;
    fMCNeutrinoTrueVertexZ=0;
  }
DEFINE_ART_MODULE(lowe::AnalyzeEvents)
