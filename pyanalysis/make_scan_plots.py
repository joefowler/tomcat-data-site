#!/usr/bin/env python
import numpy as np
import pylab as plt
import argparse
import os
import glob
import time
import common

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


def make_scan_page(filename, force=False):
    path, base = os.path.split(filename)
    partialname = base.replace("Results_", "").replace(".hdf5", "")
    output_page = "../slab2021/scans/scan_{}.md".format(partialname)
    dependencies = depends+[filename]
    if common.needs_updating(output_page, dependencies):
        _, date, index, angle, outin, _ = base.split("_")
        angle = float(angle)
        intangle = int(angle)
        outin = outin.capitalize()
        eds_plot = "eds_{}.png".format(partialname)
        scan_plot = "scan_{}.png".format(partialname)
        targets_plot = "target_{}.png".format(partialname)
        radiograph = "rg_{}.png".format(partialname)

        print("Updating scan page '{}'.".format(output_page))
        L = locals()
        contents = """---
layout: page
title: {outin} area, {angle} degrees, {date}
---
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
        print("Sll-radiographs page '{}' is up to date...skipping".format(output_page))
        return
    print("Updating all-radiographs page '{}'".format(output_page))
    with open(output_page, "w") as fp:
        header = """---
layout: page
date: {}
title: All radiographs
---

""".format(time.ctime())
        fp.write(header)
        for filename in files:
            path, base = os.path.split(filename)
            partialname = base.replace("Results_", "").replace(".hdf5", "")
            line = "![rg](/plots/Radiographs/rg_{}.png)\n".format(partialname)
            fp.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--force-update", action="store_true",
                        help="update all plots, even if they are newer than the code")
    args = parser.parse_args()
    files = get_all_hdf5s_sorted()
    for f in files:
        make_scan_page(f, force=args.force_update)
    make_radiographs_page(files, force=args.force_update)


if __name__ == "__main__":
    main()
