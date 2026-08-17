"""Pixel hit efficiency from the stored track hit pattern, on AOD.

Usage
-----
  cmsRun hitPatternFromAod_cfg.py \
      inputFiles=root://cms-xrd-global.cern.ch//store/data/Run2026D/ZeroBias/AOD/... \
      lumiMask=Cert_Collisions2026_Golden.json \
      outputFileName=hitPattern.root maxEvents=-1
"""

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

opt = VarParsing.VarParsing("analysis")
opt.register("outputFileName", "hitPattern.root",
             VarParsing.VarParsing.multiplicity.singleton,
             VarParsing.VarParsing.varType.string, "output tree file")
opt.register("lumiMask", "",
             VarParsing.VarParsing.multiplicity.singleton,
             VarParsing.VarParsing.varType.string,
             "path to the JSON file")
opt.register("tracks", "generalTracks",
             VarParsing.VarParsing.multiplicity.singleton,
             VarParsing.VarParsing.varType.string, "input track collection")
opt.register("d0Max", 0.1,
             VarParsing.VarParsing.multiplicity.singleton,
             VarParsing.VarParsing.varType.float,
             "preselection |d0| to the leading PV, cm; negative disables")
opt.register("dzMax", 0.5,
             VarParsing.VarParsing.multiplicity.singleton,
             VarParsing.VarParsing.varType.float,
             "preselection |dz| to the leading PV, cm; negative disables")
opt.register("nThreads", 1,
             VarParsing.VarParsing.multiplicity.singleton,
             VarParsing.VarParsing.varType.int, "number of threads")
opt.setDefault("maxEvents", -1)
opt.parseArguments()

process = cms.Process("HITPATTERN")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.MessageLogger.cerr.threshold = "WARNING"

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opt.maxEvents))
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(opt.inputFiles),
                            secondaryFileNames=cms.untracked.vstring())

if opt.lumiMask:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(
        filename=opt.lumiMask).getVLuminosityBlockRange()

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(opt.nThreads),
    numberOfStreams=cms.untracked.uint32(opt.nThreads),
    wantSummary=cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(opt.outputFileName))

process.hitPattern = cms.EDAnalyzer(
    "SiPixelTrackHitPatternTree",
    tracks=cms.InputTag(opt.tracks),
    primaryVertices=cms.InputTag("offlinePrimaryVertices"),
    trackPtCut=cms.double(1.0),
    trackNStripCut=cms.int32(10),
    requireHighPurity=cms.bool(True),
    pvMinNdof=cms.double(4.0),
    pvMaxAbsZ=cms.double(24.0),
    pvMaxRho=cms.double(2.0),
    pvMinNTrack=cms.int32(10),
    d0Max=cms.double(opt.d0Max),
    dzMax=cms.double(opt.dzMax),
)

process.p = cms.Path(process.hitPattern)
