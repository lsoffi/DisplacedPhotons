import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_390.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_389.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_388.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_387.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_386.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_385.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_384.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_383.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_382.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_381.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_380.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_210.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_269.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_268.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_267.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_266.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_265.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_264.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_263.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_262.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_261.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_260.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_230.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_189.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_188.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_187.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_186.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_185.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_184.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_183.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_182.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_181.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_180.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_310.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_289.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_288.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_287.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_286.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_285.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_284.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_283.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_282.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_281.root',
        'file:/eos/cms/store/group/phys_exotica/displacedPhotons/UpgradeStudies/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_cff_py_STEP0/171002_010948/0000/GMSB_L100TeV_Ctau200cm_Pythia8_13TeV_GENSIM_280.root'
    )
)


process.TFileService = cms.Service("TFileService", 
                                                      fileName = cms.string("output.root")
)

process.demo = cms.EDAnalyzer('GMSBGenAnalyzer',
genparts = cms.InputTag("prunedGenParticles")
)


process.p = cms.Path(process.demo)
