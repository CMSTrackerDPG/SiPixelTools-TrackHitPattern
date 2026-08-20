#!/usr/bin/env python3
"""Slurm / HTCondor submission, one job per run.

  ./submit.py tasks/Run3_1fb.json --create
  ./submit.py tasks/Run3_1fb.json --submit
  ./submit.py tasks/Run3_1fb.json --status

Runs come from an explicit "runs" list, from every run in the golden JSON
if "allcertified" is set, or -- if "step" is set -- from the brilcalc table in
"lumi", taking the first run past every "step" /fb of delivered luminosity.

Runs whose dataset has no complete disk replica are dropped

Task config fields
  taskname, maxevents, outdir
  cert          golden JSON, or {year: path}
  backend       "slurm" (default) or "condor"
  runs          explicit run list                        (explicit mode)
  sample        DAS dataset for those runs               (explicit mode)
  allcertified  take every run in cert                   (all-golden mode)
  lumi, step    brilcalc table and spacing in /fb        (stepped mode)
  runmin/runmax, dataset, datasetFilter, events_per_file, maxconcurrent
  walltime      slurm time limit;  flavour  condor +JobFlavour
"""

import argparse
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CFG = os.path.join(os.path.dirname(HERE), "test", "hitPatternFromAod_cfg.py")

# Measured on one 2026 run, might be out of date
SEC_PER_EVENT = 0.0075
KB_PER_EVENT = 3.2

JOB = """#!/bin/bash -e
#SBATCH --account=t3
#SBATCH --partition={partition}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time={walltime}
#SBATCH --nodes=1
#SBATCH --job-name={taskname}
#SBATCH --array=1-{njobs}{concurrent}
#SBATCH --output={taskdir}/logs/%x_%A_%a.out
#SBATCH --error={taskdir}/logs/%x_%A_%a.err

set -e
read RUN CERT <<< $(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {taskdir}/runlist.txt)
OUTPUT=hitPattern_run${{RUN}}.root

export SCRAM_ARCH={scram_arch}
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd {cmsswbase}/src
eval `scramv1 runtime -sh`

# Work on node-local /scratch
mkdir -p /scratch/$USER
WORKDIR=$(mktemp -d /scratch/$USER/{taskname}_XXXXXXXX)
trap 'cd /; rm -rf "$WORKDIR"' EXIT
cd $WORKDIR
export TMPDIR=$WORKDIR

cmsRun {cmsrun_cfg} \\
    lumiMask=$CERT \\
    maxEvents={maxevents} \\
    outputFileName=$OUTPUT \\
    inputFiles=$(paste -sd, {taskdir}/filelists/input_run${{RUN}}.txt)

if [[ "{outdir}" == /pnfs/* ]]; then
    ( eval `scram unsetenv -sh`; gfal-mkdir -p root://t3dcachedb03.psi.ch:1094/{outdir} || true )
    xrdcp -f -N $OUTPUT root://t3dcachedb03.psi.ch:1094//{outdir}/$OUTPUT
else
    mkdir -p {outdir}
    cp $OUTPUT {outdir}/$OUTPUT
fi
"""

CONDOR_JOB = """#!/bin/bash -e
set -e
read RUN CERT <<< $(sed -n "$(($1 + 1))p" {taskdir}/runlist.txt)
OUTPUT=hitPattern_run${{RUN}}.root

export SCRAM_ARCH={scram_arch}
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd {cmsswbase}/src
eval `scramv1 runtime -sh`

# Condor cleans up its own scratch dir
cd $_CONDOR_SCRATCH_DIR
export TMPDIR=$_CONDOR_SCRATCH_DIR

cmsRun {cmsrun_cfg} \\
    lumiMask=$CERT \\
    maxEvents={maxevents} \\
    outputFileName=$OUTPUT \\
    inputFiles=$(paste -sd, {taskdir}/filelists/input_run${{RUN}}.txt)

if [[ "{outdir}" == /eos/* ]]; then
    xrdcp -f -N $OUTPUT root://eoscms.cern.ch/{outdir}/$OUTPUT
else
    mkdir -p {outdir}
    cp $OUTPUT {outdir}/$OUTPUT
fi
"""

CONDOR_SUB = """executable              = {taskdir}/job.sh
arguments               = $(ProcId)
output                  = {taskdir}/logs/{taskname}_$(ClusterId).$(ProcId).out
error                   = {taskdir}/logs/{taskname}_$(ClusterId).$(ProcId).err
log                     = {taskdir}/logs/{taskname}_$(ClusterId).log
+JobFlavour             = "{flavour}"
request_cpus            = 1
request_memory          = 4000M
request_disk            = 4000M
MY.WantOS               = "el9"
transfer_output_files   = ""
use_x509userproxy       = true
MY.SendCredential       = true
max_retries             = 3
{concurrent}queue {njobs}
"""


def das(query):
    try:
        out = subprocess.run(["dasgoclient", "-query=" + query],
                             capture_output=True, text=True, timeout=1800)
    except FileNotFoundError:
        sys.exit("dasgoclient not found -- run cmsenv first")
    if out.returncode != 0:
        print("  DAS failed [%s]: %s" % (query, out.stderr.strip()))
        return []
    return [l.strip() for l in out.stdout.splitlines() if l.strip()]


def das_json(query):
    out = subprocess.run(["dasgoclient", "-query=" + query, "-json"],
                         capture_output=True, text=True, timeout=1800)
    if out.returncode != 0 or not out.stdout.strip():
        return []
    try:
        return json.loads(out.stdout)
    except ValueError:
        return []


def certified_ls(cert_path, run):
    """Certified lumisections of one run, from the golden JSON."""
    return {ls for lo, hi in json.load(open(cert_path)).get(str(run), [])
            for ls in range(lo, hi + 1)}


def pick_files(dataset, run, certls, want):
    """Files holding certified lumis, closest to the middle of the run first.

    Reading uncertified files wastes an open for no events, and the middle of
    the run is preffered because charge loss stabilizes after a few LS in a fill.
    """
    cand = []
    for rec in das_json("file,lumi dataset=%s run=%d" % (dataset, run)):
        try:
            name = rec["file"][0]["name"]
            lumis = sorted(set(rec["lumi"][0]["number"]) & certls)
        except (KeyError, IndexError, TypeError):
            continue
        if name.endswith(".root") and lumis:
            cand.append((abs(lumis[len(lumis) // 2] - middle(certls)), name))
    cand.sort()
    return [name for _, name in (cand[:want] if want else cand)]


def middle(certls):
    ls = sorted(certls)
    return ls[len(ls) // 2] if ls else 0


def parse_lumi(path):
    """(run, year, delivered) from a brilcalc ASCII table."""
    rows = []
    for line in open(path):
        line = line.strip()
        if not line.startswith("|"):
            continue
        p = [c.strip() for c in line.strip("|").split("|")]
        if len(p) != 6 or ":" not in p[0] or not p[0].split(":")[0].isdigit():
            continue
        rows.append((int(p[0].split(":")[0]),
                     "20" + p[1].split("/")[2].split()[0],
                     float(p[4])))
    rows.sort()
    return rows


def pick_runs(rows, step, runmin, runmax):
    """First run past every `step` /fb, with the cumulative reached there."""
    picked, cum, target = [], 0.0, step
    for run, year, dlv in rows:
        cum += dlv
        while cum >= target:
            if runmin <= run <= runmax:
                picked.append({"run": run, "year": year, "cum": cum})
            target += step
    seen = set()
    return [p for p in picked if not (p["run"] in seen or seen.add(p["run"]))]


def drop_uncertified(picked, cfg):
    """Runs absent from the golden JSON yield an empty tree, so skip them."""
    certs = cfg["cert"] if isinstance(cfg["cert"], dict) else None
    if not certs:
        return picked
    good = []
    for year, path in certs.items():
        runs = set(json.load(open(resolve(path))))
        good += [p for p in picked if p["year"] == year and str(p["run"]) in runs]
    dropped = len(picked) - len(good)
    if dropped:
        print("  %d runs not in the golden JSON, dropped" % dropped)
    return sorted(good, key=lambda p: p["run"])


def all_certified(cfg):
    """Every run in the golden JSON(s), with cumulative delivered lumi if known."""
    if not isinstance(cfg["cert"], dict):
        sys.exit("allcertified needs 'cert' as a {year: path} map")
    cum, tot = {}, 0.0
    if "lumi" in cfg:
        for run, _, dlv in parse_lumi(resolve(cfg["lumi"])):
            tot += dlv
            cum[run] = tot
    lo, hi = int(cfg.get("runmin", 0)), int(cfg.get("runmax", 10 ** 9))
    picked = [{"run": int(r), "year": year, "cum": cum.get(int(r), 0.0)}
              for year, path in cfg["cert"].items()
              for r in json.load(open(resolve(path))) if lo <= int(r) <= hi]
    return sorted(picked, key=lambda p: p["run"])


def on_disk(dataset):
    """True if a DISK site holds a complete replica; tape-only reads would stall."""
    out = subprocess.run(["dasgoclient", "-query=site dataset=" + dataset, "-json"],
                         capture_output=True, text=True, timeout=1800)
    if out.returncode != 0 or not out.stdout.strip():
        return False
    for entry in json.loads(out.stdout):
        for site in entry.get("site", []):
            if site.get("kind") != "DISK":
                continue
            try:
                if float(site.get("replica_fraction", "0").strip().rstrip("%")) > 99.9:
                    return True
            except ValueError:
                pass
    return False


def drop_tape_only(picked):
    ondisk = {ds for ds in sorted({p["dataset"] for p in picked}) if on_disk(ds)}
    kept = [p for p in picked if p["dataset"] in ondisk]
    if len(kept) != len(picked):
        print("  %d runs with no disk replica, dropped" % (len(picked) - len(kept)))
    return kept


def drop_done(picked, outdir):
    """Runs already written to outdir, so --create can resume a partial task."""
    have = set(os.listdir(outdir)) if os.path.isdir(outdir) else set()
    kept = [p for p in picked
            if "hitPattern_run%d.root" % p["run"] not in have]
    if len(kept) != len(picked):
        print("  %d runs already done, skipped" % (len(picked) - len(kept)))
    return kept


def resolve_datasets(picked, pattern, filt):
    """run -> dataset, one DAS call per candidate dataset rather than per run."""
    out = {}
    for year in sorted({p["year"] for p in picked}):
        datasets = [d for d in das("dataset dataset=" + pattern.format(year=year))
                    if filt in d]
        if not datasets:
            print("  %s: no dataset matching %r" % (year, filt))
            continue
        for ds in sorted(datasets):
            for r in das("run dataset=" + ds):
                if r.isdigit():
                    out[int(r)] = ds
        print("  %s: %d datasets" % (year, len(datasets)))
    return out


def create(cfg, taskdir):
    os.makedirs(os.path.join(taskdir, "filelists"), exist_ok=True)
    os.makedirs(os.path.join(taskdir, "logs"), exist_ok=True)

    if cfg.get("allcertified"):
        picked = all_certified(cfg)
        print("  %d certified runs" % len(picked))
    elif "step" in cfg:
        rows = parse_lumi(resolve(cfg["lumi"]))
        picked = pick_runs(rows, float(cfg["step"]),
                           int(cfg.get("runmin", 0)), int(cfg.get("runmax", 10 ** 9)))
        print("  %d runs at every %g /fb" % (len(picked), float(cfg["step"])))
        picked = drop_uncertified(picked, cfg)
    else:
        picked = [{"run": r, "year": None, "cum": 0.0, "dataset": cfg["sample"]}
                  for r in cfg["runs"]]

    if not picked:
        sys.exit("no runs selected")

    if "dataset" not in picked[0]:  # explicit mode already knows it
        ds_of = resolve_datasets(picked, cfg.get("dataset", "/ZeroBias/Run{year}*/AOD"),
                                 cfg.get("datasetFilter", "PromptReco"))
        picked = [p for p in picked if p["run"] in ds_of]
        for p in picked:
            p["dataset"] = ds_of[p["run"]]

    picked = drop_tape_only(picked)

    minls = int(cfg.get("minls", 0))
    if minls:
        n = len(picked)
        picked = [p for p in picked
                  if len(certified_ls(cert_for(cfg, p["year"]), p["run"])) >= minls]
        if n - len(picked):
            print("  %d runs with < %d certified LS, dropped" % (n - len(picked), minls))

    if cfg.get("resume"):
        picked = drop_done(picked, cfg["outdir"])
    if not picked:
        sys.exit("no runs selected")

    per_file = float(cfg.get("events_per_file", 3000))
    want = max(3, int(cfg["maxevents"] / per_file * 3) + 1) if cfg["maxevents"] > 0 else 0

    kept = []
    for i, p in enumerate(picked, 1):
        cached = os.path.join(taskdir, "filelists", "input_run%d.txt" % p["run"])
        if os.path.exists(cached) and os.path.getsize(cached):
            # cached lines carry the redirector prefix; keep only the LFN
            files = ["/store/" + l.strip().split("/store/", 1)[1]
                     for l in open(cached) if "/store/" in l]
        else:
            files = pick_files(p["dataset"], p["run"],
                               certified_ls(cert_for(cfg, p["year"]), p["run"]), want)
        if not files:
            print("  [%d/%d] run %d: no files" % (i, len(picked), p["run"]))
            continue
        files = files[:want] if want else files
        with open(os.path.join(taskdir, "filelists",
                               "input_run%d.txt" % p["run"]), "w") as fh:
            fh.writelines("root://cms-xrd-global.cern.ch/" + f + "\n" for f in files)
        p["nfiles"] = len(files)
        kept.append(p)
        print("  [%d/%d] run %d: %d files" % (i, len(picked), p["run"], len(files)))

    # spread concurrent jobs over storage sites; runs of one era typically share one site
    # running many jobs that read from the same site may signicantly slow down jobs
    random.Random(0).shuffle(kept)

    with open(os.path.join(taskdir, "runlist.txt"), "w") as fh:
        for p in kept:
            fh.write("%d %s\n" % (p["run"], cert_for(cfg, p["year"])))

    with open(os.path.join(taskdir, "runs.csv"), "w") as fh:
        fh.write("run,year,cumulative_delivered_fb,nfiles,dataset\n")
        for p in kept:
            fh.write("%d,%s,%.4f,%d,%s\n"
                     % (p["run"], p["year"] or "", p["cum"], p["nfiles"], p["dataset"]))

    if not cfg["outdir"].startswith("/pnfs/"):
        os.makedirs(cfg["outdir"], exist_ok=True)

    nmax = cfg.get("maxconcurrent")
    common = dict(taskname=cfg["taskname"], taskdir=taskdir, njobs=len(kept),
                  scram_arch=os.environ.get("SCRAM_ARCH", "el9_amd64_gcc13"),
                  cmsswbase=os.environ["CMSSW_BASE"], cmsrun_cfg=CFG,
                  maxevents=cfg["maxevents"], outdir=cfg["outdir"])
    path = os.path.join(taskdir, "job.sh")
    with open(path, "w") as fh:
        if backend(cfg) == "condor":
            fh.write(CONDOR_JOB.format(**common))
        else:
            fh.write(JOB.format(concurrent="%%%d" % int(nmax) if nmax else "",
                                partition=cfg.get("partition", "standard"),
                                walltime=cfg.get("walltime", "3:00:00"), **common))
    os.chmod(path, 0o755)

    if backend(cfg) == "condor":
        with open(os.path.join(taskdir, "job.sub"), "w") as fh:
            fh.write(CONDOR_SUB.format(
                flavour=cfg.get("flavour", "workday"),
                concurrent="max_materialize        = %d\n" % int(nmax) if nmax else "",
                **common))

    n = len(kept) * cfg["maxevents"]
    print("\n%d jobs, %d events each" % (len(kept), cfg["maxevents"]))
    print("~%.1f CPU-hours, ~%.1f GB" % (n * SEC_PER_EVENT / 3600.0, n * KB_PER_EVENT * 1e-6))


def status(cfg, taskdir):
    runlist = os.path.join(taskdir, "runlist.txt")
    if not os.path.exists(runlist):
        sys.exit("no runlist.txt -- run --create first")
    runs = [int(l.split()[0]) for l in open(runlist) if l.strip()]

    if cfg["outdir"].startswith("/pnfs/"):
        have = set(das_ls(cfg["outdir"]))
    else:
        have = set(os.listdir(cfg["outdir"])) if os.path.isdir(cfg["outdir"]) else set()
    done = sum("hitPattern_run%d.root" % r in have for r in runs)

    if backend(cfg) == "condor":
        q = subprocess.run(["condor_q", os.environ.get("USER", ""), "-af", "JobStatus"],
                           capture_output=True, text=True)
        states = q.stdout.split() if q.returncode == 0 else []
        pending, running = states.count("1"), states.count("2")
    else:
        q = subprocess.run(["squeue", "-u", os.environ.get("USER", ""), "-h", "-o", "%t"],
                           capture_output=True, text=True)
        states = q.stdout.split() if q.returncode == 0 else []
        pending, running = states.count("PD"), states.count("R")
    print("%s: %d runs, %d done, %d pending, %d running"
          % (cfg["taskname"], len(runs), done, pending, running))


def das_ls(outdir):
    out = subprocess.run(["gfal-ls", "root://t3dcachedb03.psi.ch:1094/" + outdir],
                         capture_output=True, text=True)
    return out.stdout.split() if out.returncode == 0 else []


def backend(cfg):
    b = cfg.get("backend", "slurm")
    if b not in ("slurm", "condor"):
        sys.exit("backend must be slurm or condor, not %r" % b)
    return b


def resolve(p):
    return p if os.path.isabs(p) else os.path.join(HERE, p)


def cert_for(cfg, year):
    c = cfg["cert"]
    if isinstance(c, dict):
        if year not in c:
            sys.exit("no cert configured for year %s" % year)
        return resolve(c[year])
    return resolve(c)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("config")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--create", action="store_true")
    g.add_argument("--submit", action="store_true")
    g.add_argument("--status", action="store_true")
    args = p.parse_args()

    if "CMSSW_BASE" not in os.environ:
        sys.exit("run cmsenv first")

    cfg = json.load(open(args.config))
    for k in ("taskname", "maxevents", "cert", "outdir"):
        if k not in cfg:
            sys.exit("config is missing: %s" % k)
    if not any(k in cfg for k in ("step", "runs", "allcertified")):
        sys.exit("config needs one of 'runs', 'step', 'allcertified'")
    for c in (cfg["cert"].values() if isinstance(cfg["cert"], dict) else [cfg["cert"]]):
        if not os.path.exists(resolve(c)):
            sys.exit("golden JSON not found: %s" % resolve(c))

    taskdir = os.path.join(HERE, cfg["taskname"])
    if args.create:
        create(cfg, taskdir)
    elif args.submit:
        job = os.path.join(taskdir, "job.sh")
        if not os.path.exists(job):
            sys.exit("no job.sh -- run --create first")
        if backend(cfg) == "condor":
            subprocess.call(["condor_submit", os.path.join(taskdir, "job.sub")])
        else:
            subprocess.call(["sbatch", job])
    else:
        status(cfg, taskdir)


if __name__ == "__main__":
    main()
