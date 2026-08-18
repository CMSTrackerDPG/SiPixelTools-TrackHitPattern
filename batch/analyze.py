#!/usr/bin/env python3
"""Hit efficiency trend vs integrated luminosity, and L1 efficiency vs eta.

  ./analyze.py # defaults below
  ./analyze.py --task Run3_1fb --etaruns 8

Reads hitPattern_run<RUN>.root from --indir and the run -> lumi mapping from
<task>/runs.csv. Efficiency is valid/(valid+missing), requiring >=3 valid pixel
hits on regions other than the one under test.

Per-run counts are cached in <task>/measurements.npz, so a second call only
reads the ROOT files of runs that are new; restyling the plots is then instant.
Pass --refresh to re-read everything.
"""

import argparse
import csv
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

HERE = os.path.dirname(os.path.abspath(__file__))

REGIONS = ["BPIX L1", "BPIX L2", "BPIX L3", "BPIX L4", "FPIX D1", "FPIX D2", "FPIX D3"]
COLORS = ["red", "blue", "green", "purple", "red", "blue", "green"]
MARKERS = ["o", "s", "^", "D", "o", "s", "^"]
BPIX = [0, 1, 2, 3]
FPIX = [4, 5, 6]
ETA_COLORS = ["red", "blue", "green", "purple", "orange", "brown", "magenta", "black"]
NPIXHIT = 3

# Year boundaries are annotated in place; the reconstruction changes go in the legend.
YEAR_MARKS = [(198.0, "2025"), (323.0, "2026")]
CHANGE_MARKS = [(171.98, "HV increase", "blue"),
                (197.48, "Clu.thr. decrease", "violet"),
                (266.4, "CPE update", "black"),
                (296.0, "Digi morphing", "darkorange"),
                (341.0, "Generic-only", "teal")]

# Gets evaluated to how many pixel layers / disks have hits 
NPIX = "+".join("(nvalid[%d]>0)" % i for i in range(7))


def measure(path, edges):
    """Per-region (num, den) and, for L1, the same binned in eta."""
    df = ROOT.RDataFrame("hitPattern/tracks", path).Define("npix", NPIX)
    books, hists = {}, {}
    for r in range(7):
        d = (df.Define("v%d" % r, "1.0*nvalid[%d]" % r) # Valid
               .Define("t%d" % r, "1.0*(nvalid[%d]+nmissing[%d])" % (r, r)) # Total
               .Filter("npix-(nvalid[%d]>0) >= %d" % (r, NPIXHIT))) # hits in >=NPIXHIT regions other than that particular layer / disk
        books[r] = (d.Sum("v%d" % r), d.Sum("t%d" % r))
        if r == 0: # Measure eta dependency in L1
            nb = len(edges) - 1
            hists["n"] = d.Histo1D(("n", "", nb, edges[0], edges[-1]), "trk_eta", "v0")
            hists["d"] = d.Histo1D(("d", "", nb, edges[0], edges[-1]), "trk_eta", "t0")

    num = np.array([books[r][0].GetValue() for r in range(7)])
    den = np.array([books[r][1].GetValue() for r in range(7)])
    hn = np.array([hists["n"].GetBinContent(i + 1) for i in range(len(edges) - 1)])
    hd = np.array([hists["d"].GetBinContent(i + 1) for i in range(len(edges) - 1)])
    return num, den, hn, hd


def eff(num, den):
    with np.errstate(divide="ignore", invalid="ignore"):
        e = np.where(den > 0, num / den, np.nan)
        err = np.where(den > 0, np.sqrt(np.abs(e * (1 - e)) / np.maximum(den, 1)), 0.0)
    return e, err


def load_cache(path, edges):
    """Per-run measurements from a previous pass, keyed by run number.

    The eta binning is baked into the stored histograms, so a cache made with a
    different --bin is discarded rather than silently mixed in."""
    if not os.path.exists(path):
        return {}
    z = np.load(path)
    if z["edges"].shape != edges.shape or not np.allclose(z["edges"], edges):
        print("cache %s uses a different eta binning, ignoring it" % path)
        return {}
    return {int(r): (z["num"][i], z["den"][i], z["etan"][i], z["etad"][i])
            for i, r in enumerate(z["runs"])}


def save_cache(path, edges, data):
    runs = sorted(data)
    np.savez_compressed(
        path, edges=edges, runs=np.array(runs, dtype=np.int64),
        num=np.array([data[r][0] for r in runs]),
        den=np.array([data[r][1] for r in runs]),
        etan=np.array([data[r][2] for r in runs]),
        etad=np.array([data[r][3] for r in runs]))


def plot_trend(x, num, den, regions, title, path):
    """Efficiency vs delivered luminosity for the given region indices."""
    fig, ax = plt.subplots(figsize=(9, 6))
    for r in regions:
        e, err = eff(num[:, r], den[:, r])
        ax.errorbar(x, e, yerr=err, color=COLORS[r], marker=MARKERS[r], ms=5,
                    ls="none", label=REGIONS[r])
    for xv, lab in YEAR_MARKS:
        ax.axvline(xv, color="gray", ls="--", lw=1.2)
        ax.text(xv + 0.005 * (max(x) - min(x)), 0.98, lab, transform=ax.get_xaxis_transform(),
                ha="left", va="top", color="gray", fontsize=9)
    for xv, lab, col in CHANGE_MARKS:
        ax.axvline(xv, color=col, ls="--", lw=1.2, label=lab)
    ax.set_xlabel("Delivered integrated luminosity [fb$^{-1}$]")
    ax.set_ylabel("Hit efficiency")
    ax.set_title(title)
    ax.grid(True)
    # errorbar containers would otherwise be listed after the plain vlines
    h, l = ax.get_legend_handles_labels()
    order = sorted(range(len(l)), key=lambda i: l[i] not in REGIONS)
    ax.legend([h[i] for i in order], [l[i] for i in order], loc="lower left", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_eta(ctr, sel, etan, etad, lumi, xlim, path):
    """L1 efficiency vs eta, one curve per selected run."""
    fig, ax = plt.subplots(figsize=(9, 6))
    for i, run in enumerate(sel):
        e, err = eff(etan[run], etad[run])
        e = np.where(etad[run] >= 200, e, np.nan)
        ax.errorbar(ctr, e, yerr=err, color=ETA_COLORS[i % len(ETA_COLORS)],
                    marker="o", ms=3, lw=1.5,
                    label="run %d  (%.0f /fb)" % (run, lumi[run]))
    ax.set_xlabel(r"track $\eta$")
    ax.set_ylabel("BPIX L1 hit efficiency")
    ax.set_title(r"BPIX L1 hit efficiency vs $\eta$")
    ax.set_xlim(*xlim)
    ax.grid(True)
    ax.legend(fontsize=9, loc="center left", bbox_to_anchor=(1.015, 0.5))
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--task", default="Run3_10fb")
    p.add_argument("--indir", default="/t3home/mrogulji/store/PixelHitPattern/Run3_10fb")
    p.add_argument("--outdir", default=os.path.join(HERE, "plots"))
    p.add_argument("--etaruns", type=int, default=6)
    p.add_argument("--bin", type=float, default=0.1)
    p.add_argument("--cache", default=None,
                   help="measurement cache (default <task>/measurements.npz)")
    p.add_argument("--refresh", action="store_true",
                   help="re-read every ROOT file instead of using the cache")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lumi = {}
    with open(os.path.join(HERE, args.task, "runs.csv")) as fh:
        for row in csv.DictReader(fh):
            lumi[int(row["run"])] = float(row["cumulative_delivered_fb"])

    files = {}
    for f in sorted(glob.glob(os.path.join(args.indir, "hitPattern_run*.root"))):
        m = re.search(r"run(\d+)\.root$", f)
        if m and int(m.group(1)) in lumi:
            files[int(m.group(1))] = f
    if not files:
        raise SystemExit("no files in %s matching %s/runs.csv" % (args.indir, args.task))

    edges = np.arange(-2.6, 2.6 + 0.5 * args.bin, args.bin)
    ctr = 0.5 * (edges[:-1] + edges[1:])

    # Reading the trees is the slow part, so each run's numbers are kept in an
    # .npz and only runs that are not in it yet are measured. Delete the file
    # (or pass --refresh) after changing anything in measure().
    cache_path = args.cache or os.path.join(HERE, args.task, "measurements.npz")
    data = {} if args.refresh else load_cache(cache_path, edges)
    todo = [r for r in sorted(files) if r not in data]
    print("%d run(s) cached, %d to measure" % (len(files) - len(todo), len(todo)))
    for i, run in enumerate(todo):
        data[run] = measure(files[run], edges)
        print("  [%d/%d] measured run %d" % (i + 1, len(todo), run))
    if todo:
        save_cache(cache_path, edges, data)
        print("cache -> %s" % cache_path)

    # Skip empty TTrees
    runs, nums, dens, etan, etad = [], [], [], {}, {}
    for run in sorted(files):
        n, d, hn, hd = data[run]
        if d[0] == 0:
            print("  run %d  %6.1f /fb  empty, skipped" % (run, lumi[run]))
            continue
        runs.append(run)
        nums.append(n), dens.append(d)
        etan[run], etad[run] = hn, hd
        e, _ = eff(n, d)
        print("  run %d  %6.1f /fb  L1 %.4f  L2 %.4f  D1 %.4f  (%d L1 meas)"
              % (run, lumi[run], e[0], e[1], e[4], d[0]))
    if not runs:
        raise SystemExit("every file was empty")
    num, den = np.array(nums), np.array(dens)

    x = np.array([lumi[r] for r in runs])

    plot_trend(x, num, den, BPIX, "BPIX hit efficiency vs luminosity",
               os.path.join(args.outdir, "eff_trend_vs_lumi_bpix.png"))
    plot_trend(x, num, den, FPIX, "FPIX hit efficiency vs luminosity",
               os.path.join(args.outdir, "eff_trend_vs_lumi_fpix.png"))

    sel = [runs[i] for i in np.linspace(0, len(runs) - 1, min(args.etaruns, len(runs))).astype(int)]
    plot_eta(ctr, sel, etan, etad, lumi, (edges[0], edges[-1]),
             os.path.join(args.outdir, "eff_l1_vs_eta.png"))

    print("\nplots -> %s" % args.outdir)


if __name__ == "__main__":
    main()
