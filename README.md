# SiPixelTools/TrackHitPattern

Pixel hit efficiency from the track hit pattern stored in AOD.
Important to be aware of the bug fixed by this PR: https://github.com/cms-sw/cmssw/pull/51478
It affects only muon tracks built with out-in tracking algorithm.
It would not store L1 (most likely missing) hits in secondary hit patterns
biasing the results towards higher efficiencies.
The recommendation is to run the tool on general tracks because of this.

## Install

```sh
cmsrel CMSSW_16_0_6 && cd CMSSW_16_0_6/src && cmsenv
git clone git@github.com:CMSTrackerDPG/SiPixelTools-TrackHitPattern.git SiPixelTools/TrackHitPattern
scram b -j 8
```

## Run

```sh
cd SiPixelTools/TrackHitPattern/test
cmsRun hitPatternFromAod_cfg.py \
    inputFiles_load=filelist.txt \
    lumiMask=Cert_Collisions2026_Golden.json \
    outputFileName=hitPattern.root maxEvents=-1
```

## Batch

Slurm, one job per run.

```sh
cd SiPixelTools/TrackHitPattern/batch
./submit.py tasks/Run2026D.json --create
./submit.py tasks/Run2026D.json --submit
./submit.py tasks/Run2026D.json --status
```

`--create` queries DAS, shuffles the file list (fixed seed per run, so it is
reproducible but not in lumisection order), and writes `job.sh`. A run is
~0.2 CPU-hours and ~0.3 GB at 100k events.

`tasks/` holds one config per era, 2025C–G and 2026B–D. Golden JSONs are in
`certs/`; 2025 must use the combined `Cert_Collisions2025_Golden.json`, since
the per-era files are missing runs.

## Output

Two trees.

### `hitPattern` — one row per selected track

| branch | meaning |
|---|---|
| `run`, `ls`, `event`, `bx` | event id |
| `npv` | good primary vertices in the event (pileup proxy) |
| `pv_ntrk`, `pv_z` | tracks on, and z of, the leading PV |
| `trk_pt`, `trk_eta`, `trk_phi` | track kinematics |
| `trk_d0`, `trk_dz` | distance to the **leading PV** in cm |
| `trk_algo` | `reco::TrackBase::TrackAlgorithm` (for debugging) |
| `nvalid[7]`, `nmissing[7]`, `ninactive[7]` | per region: BPIX L1–L4, then FPIX D1–D3 |

Counts hit statuses per layer so that (in an unlikely case) 
a track crossing a layer twice through a ladder overlap is not collapsed.

Efficiency to be measured is `valid / (valid + missing)`.

### `jobInfo` — one row per job

Bookkeeping with the cuts actually applied and counters.

## Selection

Applied here:

| cut | default |
|---|---|
| `highPurity`, `pt > 1`, `nstrip > 10` | on | 
| good leading PV, `ndof >= 4`, \|z\| < 24, rho < 2, `ntrk > 10` | on |
| \|d0\| < `d0Max` = 0.1 cm | on | 
| \|dz\| < `dzMax` = 0.5 cm | on |

The d0/dz defaults are looser than any working point the analysis is expected to
use, 0.5 cm in dz is the loosest of the DQM reference selections (0.1 in BPIX).

Everything layer-dependent, including the requirement of valid pixel hits on the
*other* layers, is left to analysis.
