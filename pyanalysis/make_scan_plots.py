#!/usr/bin/env python
import numpy as np
import pylab as plt
import argparse
import os
import glob
import time

import common
import onescan

depends = [__file__, "common.py"]
DIR = "/Volumes/Data2017/Data_Tomcat/tomcat2021/slab/"


def get_all_hdf5s_sorted():
    files = glob.glob(DIR+"LaLB/Results*.hdf5")
    files += glob.glob(DIR+"La/Results*.hdf5")

    def sortkey(path):
        fname = os.path.split(path)[1]
        return fname.lower()
    files.sort(key=sortkey)
    return files


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
            scan.eds_summary()
            plt.tight_layout()
            plt.savefig(output_plot)

        output_plot = "../plots/Targets/target_{}.png".format(key)
        dependencies = depends+[filename, "onescan.py"]
        if common.needs_updating(output_plot, dependencies):
            common.myfigure((7, 2.6), 24)
            scan.plot_targets()
            plt.tight_layout()
            plt.savefig(output_plot)

        output_plot = "../plots/TESrates/rate_{}.png".format(key)
        dependencies = depends+[filename, "onescan.py"]
        if common.needs_updating(output_plot, dependencies):
            common.myfigure((5.5, 4.5), 25)
            scan.tes_rate()
            plt.tight_layout()
            plt.savefig(output_plot)

        output_plot = "../plots/Radiographs/rg_{}.png".format(key)
        dependencies = depends+[filename, "onescan.py", "radiograph.py"]
        if common.needs_updating(output_plot, dependencies):
            rg = onescan.compute_radiograph(scan, voxsize_nm=60)
            common.myfigure((9.5, 4), 26)
            rg.plot()
            plt.tight_layout()
            plt.savefig(output_plot)

    except Exception as e:
        print("Failed on {} with ".format(filename), e)


def make_scan_page(filename, force=False, prev=None, next=None):
    path, base = os.path.split(filename)
    partialname = base.replace("Results_", "").replace(".hdf5", "")
    output_page = "../slab2021/scans/scan_{}.md".format(partialname)
    dependencies = depends+[filename]
    if common.needs_updating(output_page, dependencies) or force:
        _, date, index, angle, outin, _ = base.split("_")
        angle = float(angle)
        intangle = int(angle)
        strangle = "{:+5.1f}".format(angle)
        outin = outin.capitalize()
        eds_plot = "eds_{}.png".format(partialname)
        scan_plot = "scan_{}.png".format(partialname)
        targets_plot = "target_{}.png".format(partialname)
        radiograph = "rg_{}.png".format(partialname)
        prevpage = nextpage = ""
        if prev is not None:
            pbase = os.path.split(prev)[1]
            prevpage = "prev: scan_{}.html\n".format(
                pbase.replace("Results_", "").replace(".hdf5", ""))
        if next is not None:
            nbase = os.path.split(next)[1]
            nextpage = "next: scan_{}.html\n".format(
                nbase.replace("Results_", "").replace(".hdf5", ""))

        print("Updating scan page '{}'.".format(output_page))
        L = locals()
        contents = """---
layout: page
{prevpage}{nextpage}title: {outin} area, {strangle} degrees, {date}
---
[Up: Index of all scans](/slab2021/scans)

Scan file name: {base}
![](/plots/cartoons/angle{intangle}.png)

![](/plots/Radiographs/rg_{partialname}.png)

![](/plots/TESrates/rate_{partialname}.png)

![](/plots/Targets/{targets_plot})

![](/plots/ScanPlots/{scan_plot})

![](/plots/EDSPlots/{eds_plot})
""".format(**L)
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


def main():
    parser = argparse.ArgumentParser(description="Make scan plots and web pages.")
    parser.add_argument("-f", "--force-update", action="store_true",
                        help="update all plots, even if they are newer than the code")
    parser.add_argument("-n", "--num-files", default=-1, type=int,
                        help="max # of files to process (-1=all)")
    parser.add_argument("--no_pages", action="store_true", help="skip making page markdown files")
    parser.add_argument("--no_plots", action="store_true", help="skip making plots")
    args = parser.parse_args()
    files = get_all_hdf5s_sorted()
    if args.num_files >= 0:
        files = files[:args.num_files]

    print("Running with {} files, no_pages={}, no_plots={}".format(
        len(files), args.no_pages, args.no_plots))

    if not args.no_plots:
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


if __name__ == "__main__":
    main()
