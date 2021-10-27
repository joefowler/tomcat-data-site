#!/usr/bin/env python
import numpy as np
import pylab as plt
import argparse
import os
# import glob
import time
# from collections import OrderedDict
import h5py
import tempfile
#
import common
import onescan


def generate_alldata_summary(files, force=False):
    target = os.path.join(common.DIR, "alldata_summary.hdf5")
    dependencies = [__file__, "common.py", "onescan.py"]
    if not (force or common.needs_updating(target, dependencies)):
        print("Summary HDF5 '{}' is up to date. Skipping...".format(target))
        return

    # daynum = 0
    eds = []
    durations = []
    starts = []
    x, y = [], []
    idx = []
    names = []
    nsteps = []
    for i, file in enumerate(files):
        try:
            scan = onescan.OneScan(file)
            print("Analyzing {}".format(file))
        except Exception:
            print("Analyzing {} ...failed. Skipping to next file.".format(file))
            continue
        eds.append(scan.eds)
        durations.append(scan.duration)
        starts.append(scan.h5["times/start/sec"][:])
        x.append(scan.drift[:, 0])
        y.append(scan.drift[:, 1])
        N = len(scan.eds)
        nsteps.append(N)
        idx.append(np.repeat(i, N))
        names.append(scan.basename)

    # Store results to HDF5
    _fd, fname = tempfile.mkstemp(suffix=".hdf5")
    with h5py.File(fname, "a") as h:
        h.create_dataset("eds", data=np.hstack(eds))
        h.create_dataset("durations", data=np.hstack(durations))
        h.create_dataset("starts", data=np.hstack(starts))
        h.create_dataset("x", data=np.hstack(x))
        h.create_dataset("y", data=np.hstack(y))
        h.create_dataset("nsteps", data=np.array(nsteps))
        h.create_dataset("idx", data=np.hstack(idx))
        h.create_dataset("names", data=names)

    os.rename(fname, target)


def plot_eds(h=None):
    if h is None:
        target = os.path.join(common.DIR, "alldata_summary.hdf5")
        with h5py.File(target, "r") as h:
            plot_eds(h)
    t0 = time.mktime((2021, 8, 31, 12, 0, 0, 0, 0, 1))  # Noon MDT on Aug 31
    t = h["starts"][:]-t0
    day = t//86400
    plt.clf()
    rate = h["eds"][:]/h["durations"][:]
    plt.scatter(t/86400, rate, c=day % 12, cmap="Paired", s=3)
    plt.xlabel("Date (1.00=noon MDT Sept 1; 2.00=noon Sept 2)")
    plt.ylabel("EDS Rate")
    plt.title("EDS rate over time")
    plt.ylim([0, rate.max()*1.01])
    plt.xlim([2, t.max()/86400+0.1])
    plt.grid(True)


def main():
    parser = argparse.ArgumentParser(description="Make overview plots and web pages.")
    parser.add_argument("-f", "--force-update", action="store_true",
                        help="update all plots, even if they are newer than the code")
    parser.add_argument("-n", "--num-files", default=-1, type=int,
                        help="max # of files to process (-1=all)")
    parser.add_argument("--no_pages", action="store_true", help="skip making page markdown files")
    parser.add_argument("--no_plots", action="store_true", help="skip making plots")
    args = parser.parse_args()
    allrawfiles = common.get_all_hdf5s_sorted()
    if args.num_files >= 0:
        allrawfiles = allrawfiles[:args.num_files]

    print("Running with {} files, no_pages={}, no_plots={}".format(
        len(allrawfiles), args.no_pages, args.no_plots))

    generate_alldata_summary(allrawfiles, force=args.force_update)

    # if not args.no_plots:
    #     plt.ioff()
    #     for i, f in enumerate(files):
    #         prev = None
    #         if i > 0:
    #             prev = files[i-1]
    #         next = None
    #         if i < len(files)-1:
    #             next = files[i+1]
    #         make_scan_plots(f, force=args.force_update)
    #
    # if not args.no_pages:
    #     for i, f in enumerate(files):
    #         prev = None
    #         if i > 0:
    #             prev = files[i-1]
    #         next = None
    #         if i < len(files)-1:
    #             next = files[i+1]
    #     #     make_scan_page(f, force=args.force_update, prev=prev, next=next)
    #     # make_radiographs_page(files, force=args.force_update)
    #     # make_summary_table(files, force=args.force_update)


if __name__ == "__main__":
    main()
