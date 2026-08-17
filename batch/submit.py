#!/usr/bin/env python3
"""Slurm submission, one job per run.

  ./submit.py tasks/Run2026D.json --create
  ./submit.py tasks/Run2026D.json --submit
  ./submit.py tasks/Run2026D.json --status

Task config fields: taskname, sample, runs, maxevents, cert, outdir.
Optional: walltime, events_per_file.
"""

import argparse
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CFG = os.path.join(os.path.dirname(HERE), "test", "hitPatternFromAod_cfg.py")
REQUIRED = ["taskname", "sample", "runs", "maxevents", "cert", "outdir"]

# Measured on 2026D AOD.
SEC_PER_EVENT = 0.0075
KB_PER_EVENT = 3.2

JOB = """#!/bin/bash -e
#SBATCH --account=t3
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time={walltime}
#SBATCH --nodes=1
#SBATCH --job-name={taskname}
#SBATCH --array=1-{njobs}
#SBATCH --output={taskdir}/logs/%x_%A_%a.out
#SBATCH --error={taskdir}/logs/%x_%A_%a.err

set -e
RUN=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {taskdir}/runlist.txt)
OUTPUT=hitPattern_run${{RUN}}.root

export SCRAM_ARCH={scram_arch}
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd {cmsswbase}/src
eval `scramv1 runtime -sh`
cd $TMPDIR

cmsRun {cmsrun_cfg} \\
    lumiMask={cert} \\
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


def das_files(dataset, run):
    q = "file dataset=%s run=%d" % (dataset, run)
    try:
        out = subprocess.run(["dasgoclient", "-query=" + q],
                             capture_output=True, text=True, timeout=900)
    except FileNotFoundError:
        sys.exit("dasgoclient not found -- run cmsenv first")
    if out.returncode != 0:
        print("  run %d: DAS query failed: %s" % (run, out.stderr.strip()))
        return []
    return [l.strip() for l in out.stdout.splitlines() if l.strip().endswith(".root")]


def create(cfg, taskdir, cert):
    os.makedirs(os.path.join(taskdir, "filelists"), exist_ok=True)
    os.makedirs(os.path.join(taskdir, "logs"), exist_ok=True)

    # Files come back from DAS in lumisection order, so taking the first N
    # events would sample only the start of the run. Fixed seed per run keeps
    # the choice reproducible.
    per_file = float(cfg.get("events_per_file", 3000))
    want = max(3, int(cfg["maxevents"] / per_file * 3) + 1) if cfg["maxevents"] > 0 else 0

    kept = []
    for run in cfg["runs"]:
        files = das_files(cfg["sample"], run)
        if not files:
            continue
        random.Random(run).shuffle(files)
        if want:
            files = files[:want]
        with open(os.path.join(taskdir, "filelists", "input_run%d.txt" % run), "w") as fh:
            fh.writelines("root://cms-xrd-global.cern.ch/" + f + "\n" for f in files)
        kept.append(run)
        print("  run %d: %d files" % (run, len(files)))

    if not kept:
        sys.exit("no runs with files")

    with open(os.path.join(taskdir, "runlist.txt"), "w") as fh:
        fh.writelines("%d\n" % r for r in kept)

    path = os.path.join(taskdir, "job.sh")
    with open(path, "w") as fh:
        fh.write(JOB.format(
            taskname=cfg["taskname"],
            taskdir=taskdir,
            njobs=len(kept),
            walltime=cfg.get("walltime", "3:00:00"),
            scram_arch=os.environ.get("SCRAM_ARCH", "el9_amd64_gcc13"),
            cmsswbase=os.environ["CMSSW_BASE"],
            cmsrun_cfg=CFG,
            cert=cert,
            maxevents=cfg["maxevents"],
            outdir=cfg["outdir"]))
    os.chmod(path, 0o755)

    n = len(kept) * cfg["maxevents"]
    print("\n%d runs, %d events each" % (len(kept), cfg["maxevents"]))
    print("~%.2f CPU-hours, ~%.2f GB" % (n * SEC_PER_EVENT / 3600.0, n * KB_PER_EVENT * 1e-6))


def status(cfg, taskdir):
    runlist = os.path.join(taskdir, "runlist.txt")
    if not os.path.exists(runlist):
        sys.exit("no runlist.txt -- run --create first")
    runs = [int(l) for l in open(runlist) if l.strip()]

    done = 0
    for run in runs:
        name = "hitPattern_run%d.root" % run
        if cfg["outdir"].startswith("/pnfs/"):
            rc = subprocess.run(
                ["gfal-stat", "root://t3dcachedb03.psi.ch:1094/%s/%s" % (cfg["outdir"], name)],
                capture_output=True)
            done += rc.returncode == 0
        else:
            done += os.path.exists(os.path.join(cfg["outdir"], name))

    q = subprocess.run(["squeue", "-u", os.environ.get("USER", ""), "-h", "-o", "%t"],
                       capture_output=True, text=True)
    states = q.stdout.split() if q.returncode == 0 else []
    print("%s: %d runs, %d done, %d pending, %d running"
          % (cfg["taskname"], len(runs), done, states.count("PD"), states.count("R")))


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
    missing = [k for k in REQUIRED if k not in cfg]
    if missing:
        sys.exit("config is missing: %s" % ", ".join(missing))

    cert = cfg["cert"] if os.path.isabs(cfg["cert"]) else os.path.join(HERE, cfg["cert"])
    if not os.path.exists(cert):
        sys.exit("golden JSON not found: %s" % cert)

    taskdir = os.path.join(HERE, cfg["taskname"])
    if args.create:
        create(cfg, taskdir, cert)
    elif args.submit:
        job = os.path.join(taskdir, "job.sh")
        if not os.path.exists(job):
            sys.exit("no job.sh -- run --create first")
        subprocess.call(["sbatch", job])
    else:
        status(cfg, taskdir)


if __name__ == "__main__":
    main()
