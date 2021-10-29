#!/usr/bin/env python
import numpy as np
import pylab as plt
import argparse
import os
import glob
import time
from collections import OrderedDict

import common
import onescan

depends = [__file__, "common.py"]
DIR = "/Users/fowlerj/data/tomcat2021/slab/"


def make_scan_plots(filename, force=False):
    path, base = os.path.split(filename)
    partialname = base.replace("Results_", "").replace(".hdf5", "")

    plots_needed = (
        ("EDSPlots", "eds"),
        ("EDSPlots", "scan"),
        ("TESrates", "rate"),
        ("ScanPlots", "scan"),
        ("Radiographs", "rg"),
    )

    needany = force
    if not force:
        for (dir, prefix) in plots_needed:
            output_plot = ("../plots/{}/{}_{}.png".format(dir, prefix, partialname))
            dependencies = depends+[filename, "onescan.py", "radiograph.py"]
            if common.needs_updating(output_plot, dependencies):
                needany = True
                break
    if not needany:
        print("Plots for '{}' are up to date...skipping.".format(filename))
        return

    try:
        print("Updating scan plots for '{}'.".format(filename))
        scan = onescan.OneScan(filename)
        basename = os.path.splitext(os.path.basename(filename))[0]
        key = basename.replace("Results_", "")

        output_plot = "../plots/ScanPlots/scan_{}.png".format(key)
        dependencies = depends+[filename, "onescan.py"]
        if common.needs_updating(output_plot, dependencies):
            common.myfigure((9.4, 10), 21)
            onescan.summarize_scan(scan)
            plt.tight_layout()
            plt.savefig(output_plot)

        output_plot = "../plots/EDSPlots/eds_{}.png".format(key)
        dependencies = depends+[filename, "onescan.py"]
        if common.needs_updating(output_plot, dependencies):
            common.myfigure((6, 7), 22)
            scan.plot_eds_summary()
            plt.tight_layout()
            plt.savefig(output_plot)

        output_plot = "../plots/Targets/target_{}.png".format(key)
        dependencies = depends+[filename, "onescan.py"]
        if common.needs_updating(output_plot, dependencies):
            common.myfigure((7, 2.6), 24)
            scan.plot_targets()
            plt.tight_layout()
            plt.savefig(output_plot)

        output_plot = "../plots/TimeWeightedDrift/twd_{}.png".format(key)
        dependencies = depends+[filename, "onescan.py"]
        if common.needs_updating(output_plot, dependencies):
            common.myfigure((5.5, 3.5), 25)
            scan.plot_timewt_drift()
            plt.tight_layout()
            plt.savefig(output_plot)

        output_plot = "../plots/TESrates/rate_{}.png".format(key)
        dependencies = depends+[filename, "onescan.py"]
        if common.needs_updating(output_plot, dependencies):
            common.myfigure((5.5, 4.5), 26)
            scan.plot_tes_rate()
            plt.tight_layout()
            plt.savefig(output_plot)

    except Exception as e:
        print("Failed on {} with ".format(filename), e)


def make_scan_page(filename, force=False, prev=None, next=None):
    path, base = os.path.split(filename)
    basebackslashed = base.replace("_", "\_")
    partialname = base.replace("Results_", "").replace(".hdf5", "")
    output_page = "../slab2021/scans/scan_{}.md".format(partialname)
    dependencies = depends+[filename]
    if common.needs_updating(output_page, dependencies) or force:
        _, date, index, angle, outin, _ = base.split("_")
        try:
            angle = float(angle)
        except Exception as e:
            print("Problem in file {} with angle {}".format(base, angle))
            raise e
        intangle = int(angle)
        strangle = "{:+5.1f}".format(angle)
        outin = outin.capitalize()
        eds_plot = "eds_{}.png".format(partialname)
        scan_plot = "scan_{}.png".format(partialname)
        targets_plot = "target_{}.png".format(partialname)
        radiograph = "rg_{}.png".format(partialname)
        twd_plot = "twd_{}.png".format(partialname)
        prevpage = nextpage = ""
        if prev is not None:
            pbase = os.path.split(prev)[1]
            prevpage = "prev: scan_{}.html\n".format(
                pbase.replace("Results_", "").replace(".hdf5", ""))
        if next is not None:
            nbase = os.path.split(next)[1]
            nextpage = "next: scan_{}.html\n".format(
                nbase.replace("Results_", "").replace(".hdf5", ""))

        plotfiles = OrderedDict()
        plotfiles["Radiograph"] = "/plots/Radiographs/rg_{}.png".format(partialname)
        plotfiles["Offset squares"] = "/plots/Targets/{}".format(targets_plot)
        plotfiles["Scan analysis"] = "/plots/ScanPlots/{}".format(scan_plot)
        plotfiles["Time-weighted drift"] = "/plots/TimeWeightedDrift/{}".format(twd_plot)
        plotfiles["TES rates"] = "/plots/TESrates/rate_{}.png".format(partialname)
        plotfiles["EDS analysis"] = "/plots/EDSPlots/{}".format(eds_plot)
        allplots_markdown = "\n\n".join(["![]({})".format(f) for f in plotfiles.values()])

        pagedate = time.ctime()
        print("Updating scan page '{}'.".format(output_page))
        L = locals()
        contents = """---
layout: page
{prevpage}{nextpage}title: {outin} area, {strangle} degrees, {date}
---
[Up: Index of all scans](/slab2021/scans)

Scan file name: {basebackslashed}

{allplots_markdown}

# Page produced

| What | Date |
| ---- | ---- |
| Page text | {pagedate} |
""".format(**L)
        plottimes = []
        for k in plotfiles:
            file = "..{}".format(plotfiles[k])
            try:
                mtime = os.stat(file).st_mtime
            except FileNotFoundError:
                continue
            date = time.ctime(mtime)
            t = "| {} | {} |".format(k, date)
            plottimes.append(t)
        contents += "\n".join(plottimes)

        with open(output_page, "w") as fp:
            fp.write(contents)

    else:
        print("Page '{}' is up to date...skipping.".format(output_page))


def make_radiographs_page(files, force=False):
    output_page = "../slab2021/radiographs.md"
    if not common.needs_updating(output_page, depends):
        print("All-radiographs page '{}' is up to date...skipping".format(output_page))
        return
    print("Updating all-radiographs page '{}'".format(output_page))
    with open(output_page, "w") as fp:
        header = """---
layout: page
date: {}
title: All radiographs
---

[Up: slab target](/slab2021/)

""".format(time.ctime())
        fp.write(header)
        for filename in files:
            path, base = os.path.split(filename)
            partialname = base.replace("Results_", "").replace(".hdf5", "")
            line = "![rg](/plots/Radiographs/rg_{}.png)\n".format(partialname)
            fp.write(line)


class ScanSummary():
    def __init__(self, scan=None):
        if scan is None:
            return
        self.angle = scan.angle
        self.file = scan.filename
        self.nscans = 1
        if "inner" in scan.filename:
            self.nouter, self.ninner = 0, 1
        else:
            self.nouter, self.ninner = 1, 0
        self.duration = scan.duration.sum()
        self.nphot_total = 0
        self.nphot_pt = 0
        for cnum in scan.xdet.keys():
            chan = scan.h5["tes/chan{:03d}".format(cnum)]
            self.nphot_total += np.sum(chan["counts"])
            self.nphot_pt += np.sum(chan["sig"], axis=0)[0]
        self.bytes = os.stat(scan.h5.filename).st_size

    def __add__(self, other):
        sum = ScanSummary()
        if self.angle == other.angle:
            sum.angle = self.angle
        else:
            sum.angle = "All angles"
        sum.file = ", ".join((self.file, other.file))
        sum.nscans = self.nscans + other.nscans
        sum.nouter = self.nouter + other.nouter
        sum.ninner = self.ninner + other.ninner
        sum.duration = self.duration + other.duration
        sum.nphot_pt = self.nphot_pt + other.nphot_pt
        sum.nphot_total = self.nphot_total + other.nphot_total
        sum.bytes = self.bytes + other.bytes
        return sum

    def tableline(self):
        if isinstance(self.angle, str):
            angle = self.angle
        elif self.angle == 0.0:
            angle = "0.0&deg;"
        else:
            angle = "{:+.1f}&deg;".format(self.angle)
            angle = angle.replace("-", "&minus;")
        return "| {:17s} | {:d} | {:d}+{:d} | {:.2f} | {:.1f} | {:.0f} | {:.1f} MB |".format(
            angle, self.nscans, self.nouter, self.ninner, self.duration/3600,
            self.nphot_pt/1e6, self.nphot_total/1e6, self.bytes/1e6)


def make_summary_table(files, force=False):
    N = len(files)
    results = {}
    for file in files:
        print("Analyzing {} for summary".format(file))
        try:
            s = onescan.OneScan(file)
            summary = ScanSummary(s)
            a = s.angle
            if a in results:
                results[a] = results[a] + summary
            else:
                results[a] = summary
        except Exception as e:
            print("...failed with:")
            print(e)

    keys = list(results.keys())
    keys.sort()
    allscans = results[keys[0]]
    for k in keys[1:]:
        allscans = allscans + results[k]

    with open("../slab2021/data_summary.md", "w") as fp:
        lines = []
        lines.append(
            "| Angle | Scans<br/>total | Scans<br/>out+in | Time (hr) | Photons/1M<br/>Pt L&alpha; | Photons/1M<br/>TES Total | File size |")
        lines.append("| ----: | ----: | -------: | --------: | ----------: | ----: | -----: |")
        lines.append(allscans.tableline())
        for k in keys:
            lines.append(results[k].tableline())

        fp.write("\n".join(lines))


def make_good_dets_table():
    files = common.get_all_hdf5s_sorted()
    ratios = []
    with open("good_dets_table.csv", "w") as fp:
        for f in files:
            scan = onescan.OneScan(f)
            r = np.array([v for v in scan.medianrates.values()])
            ratios.append(r/scan.medianrate)

            # Report
            ratio = np.array([scan.medianrates.get(c, 0)/scan.medianrate for c in range(241)])
            good = np.abs(ratio-1.0) < 0.1
            goodlist = ",".join([str(i) for i in np.arange(241)[good]])
            fp.write("{},{},{}\n".format(scan.basename, good.sum(), goodlist))
            print("{:40s}  {:3d} good sensors with ±10% of median rate of {:.2f} cps".format(
                scan.basename, good.sum(), scan.medianrate))

    return np.hstack(ratios)


def report_eds_values(stop_after=None):
    results = OrderedDict()
    for angle in np.linspace(-45, 45, 13):
        results[angle] = []

    files = common.get_all_hdf5s_sorted()
    print("Filename                                    Angle   Conv factor")
    print("--------                                   ------   -----------")
    for f in files[:stop_after]:
        scan = onescan.OneScan(f, verbose=False)
        esf = scan.compute_eds_scale_factor()
        print("{:43s} {:5.1f}   {:.6f}".format(scan.basename, scan.angle, esf))
        results[scan.angle].append(esf)

    def median_abs_dev(x):
        return np.median(np.abs(x-np.median(x)))
    print()
    print("Angle      Med(conv)  MedAbsDev(conv)")
    print("---------  ---------  ---------------")
    for a, r in results.items():
        if len(r) == 0:
            continue
        print("{:5.1f} deg  {:.6f}      {:.6f}".format(a, np.median(r), median_abs_dev(r)))


def main():
    parser = argparse.ArgumentParser(description="Make scan plots and web pages.")
    parser.add_argument("-f", "--force-update", action="store_true",
                        help="update all plots, even if they are newer than the code")
    parser.add_argument("-n", "--num-files", default=-1, type=int,
                        help="max # of files to process (-1=all)")
    parser.add_argument("--no_pages", action="store_true", help="skip making page markdown files")
    parser.add_argument("--no_plots", action="store_true", help="skip making plots")
    args = parser.parse_args()
    files = common.get_all_hdf5s_sorted()
    if args.num_files >= 0:
        files = files[:args.num_files]

    print("Running with {} files, no_pages={}, no_plots={}".format(
        len(files), args.no_pages, args.no_plots))

    if not args.no_plots:
        plt.ioff()
        for i, f in enumerate(files):
            prev = None
            if i > 0:
                prev = files[i-1]
            next = None
            if i < len(files)-1:
                next = files[i+1]
            make_scan_plots(f, force=args.force_update)

    if not args.no_pages:
        for i, f in enumerate(files):
            prev = None
            if i > 0:
                prev = files[i-1]
            next = None
            if i < len(files)-1:
                next = files[i+1]
            make_scan_page(f, force=args.force_update, prev=prev, next=next)
        make_radiographs_page(files, force=args.force_update)
        make_summary_table(files, force=args.force_update)


if __name__ == "__main__":
    main()
