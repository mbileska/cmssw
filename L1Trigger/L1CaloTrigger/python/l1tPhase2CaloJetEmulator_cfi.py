import FWCore.ParameterSet.Config as cms

l1tPhase2CaloJetEmulator = cms.EDProducer("Phase2L1CaloJetEmulator",
						gctFullTowers = cms.InputTag("l1tPhase2L1CaloEGammaEmulator","GCTFullTowers"),
						hgcalTowers = cms.InputTag("l1tHGCalTowerProducer","HGCalTowerProcessor"),
                                                hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis"),
                                                nHits_to_nvtx_params = cms.VPSet(
                                                     cms.PSet(
                                                         fit = cms.string( "hgcalEM" ),
                                                         params = cms.vdouble( 157.522, 0.090 )
                                                     ),
                                                     cms.PSet(
                                                         fit = cms.string( "hgcalHad" ),
                                                         params = cms.vdouble( 159.295, 0.178 )
                                                     ),
                                                     cms.PSet(
                                                         fit = cms.string( "hf" ),
                                                         params = cms.vdouble( 165.706, 0.153 )
                                                     ),
                                                ),
                                                nvtx_to_PU_sub_params = cms.VPSet(
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalEM" ),
                                                         iEta = cms.string( "er1p4to1p8" ),
                                                         params = cms.vdouble( -0.011772, 0.004142 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalEM" ),
                                                         iEta = cms.string( "er1p8to2p1" ),
                                                         params = cms.vdouble( -0.015488, 0.005410 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalEM" ),
                                                         iEta = cms.string( "er2p1to2p4" ),
                                                         params = cms.vdouble( -0.021150, 0.006078 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalEM" ),
                                                         iEta = cms.string( "er2p4to2p7" ),
                                                         params = cms.vdouble( -0.015705, 0.005339 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalEM" ),
                                                         iEta = cms.string( "er2p7to3p1" ),
                                                         params = cms.vdouble( -0.018492, 0.005620 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalHad" ),
                                                         iEta = cms.string( "er1p4to1p8" ),
                                                         params = cms.vdouble( 0.005675, 0.000615 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalHad" ),
                                                         iEta = cms.string( "er1p8to2p1" ),
                                                         params = cms.vdouble( 0.004560, 0.001099 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalHad" ),
                                                         iEta = cms.string( "er2p1to2p4" ),
                                                         params = cms.vdouble( 0.000036, 0.001608 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalHad" ),
                                                         iEta = cms.string( "er2p4to2p7" ),
                                                         params = cms.vdouble( 0.000869, 0.001754 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hgcalHad" ),
                                                         iEta = cms.string( "er2p7to3p1" ),
                                                         params = cms.vdouble( -0.006574, 0.003134 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hf" ),
                                                         iEta = cms.string( "er29to33" ),
                                                         params = cms.vdouble( -0.203291, 0.044096 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hf" ),
                                                         iEta = cms.string( "er34to37" ),
                                                         params = cms.vdouble( -0.210922, 0.045628 )
                                                     ),
                                                     cms.PSet(
                                                         calo = cms.string( "hf" ),
                                                         iEta = cms.string( "er38to41" ),
                                                         params = cms.vdouble( -0.229562, 0.050560 )
                                                     ),
                                                ),
)
