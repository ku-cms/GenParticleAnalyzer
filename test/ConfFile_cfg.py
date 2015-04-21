import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00D2084E-C869-E411-A5E2-0025B3E063FA.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/24B6B427-DB69-E411-9490-001E67398CA0.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A877FC1-E069-E411-B1D5-002481E150C2.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7218631F-D069-E411-ABD3-0025B3E05D0A.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8094C03C-E969-E411-B347-002590A370DC.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/860C971D-CD69-E411-8E13-002590A8882C.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C70F6BE-D469-E411-BE95-0025901248FA.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A05AE83E-D969-E411-B5D0-0025902008E4.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AE442967-D769-E411-8BB2-001E67397CBA.root',
      '/store/mc/Phys14DR/TprimeJetToTH_M800GeV_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B00DB2F6-E469-E411-9A3B-002590A88830.root',

    )
)

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("GenPlots.root")
       )

process.gen = cms.EDAnalyzer('GenParticleAnalyzer'
)


process.p = cms.Path(process.gen)
