import numpy as np
import scipy as sp
import pylab as plt
import h5py
import mass
import os
import glob

import radiograph


class OneScan():
    """
    Represent one raster scan on the target.
    """

    def __init__(self, filename):
        f = h5py.File(filename, "r")
        self.h5 = f
        self.filename = "/".join(filename.split(os.path.sep)[-3:])
        self.basename = os.path.basename(filename)
        self.angle = float(self.basename.split("_")[3])

        self.good = f["detInfo/goodChans"][:]
        try:
            f["eds"]
            self.eds = f["eds/sig"][:, 0]
            self.eds_unc = f["eds/unc"][:, 0]
        except KeyError:
            print("No EDS data in {}. Yikes".format(filename))
        self.target = f["steps/corrected_sampleFixed"][:, :2]
        self.targetend = f["steps/endPos_sampleFixed"][:, :2]
        self.command = f["steps/command_sampleFixed"][:, :2]
        self.offset = self.command-self.target
        self.duration = f["times/duration"][:]
        self.medianrates = {}
        self.xdet = {}
        self.ydet = {}
        for k in f["tes"].keys():
            if not k.startswith("chan"):
                continue
            cnum = int(k.replace("chan", ""))
            self.medianrates[cnum] = np.median(f["tes/{}/sig".format(k)][:, 0]/self.duration)
            self.xdet[cnum] = f["tes/{}/xdet".format(k)][0, 0]
            self.ydet[cnum] = f["tes/{}/ydet".format(k)][0, 0]
        self.medianrate = np.median([x for x in self.medianrates.values()])
        self.Nx, self.Xcoords = findgrid(self.target[:, 0])
        self.Ny, self.Ycoords = findgrid(self.target[:, 1])

    def tes_rate(self):
        plt.clf()
        plt.subplot(111, aspect="equal")
        plt.xlabel("TES X position (µm)")
        plt.ylabel("TES Y position (µm)")
        plt.title("Pt L$\\alpha$ rate (cps) {}".format(self.basename[:-5]))
        cm = plt.cm.inferno
        for k, v in self.medianrates.items():
            plt.plot(self.xdet[k], self.ydet[k], "s", ms=14, color=cm(.5*v/self.medianrate))
            if v > self.medianrate*1.5:
                color = "k"
            else:
                color = "w"
            plt.text(self.xdet[k], self.ydet[k], "{:.1f}".format(
                v), color=color, ha="center", va="center", fontsize="x-small")

    def eds_summary(self):
        "Plot a summary of the EDS rate vs time for this scan."
        plt.clf()
        ax = plt.subplot(211)
        idx = np.arange(1, 1+len(self.eds))
        dur = self.duration
        e = self.eds/dur
        de = self.eds_unc/dur
        scale = np.median(e)
        e /= scale
        de /= scale

        j = 0
        y = self.target[:, 1]
        plt.plot(idx[[0, -1]], [1, 1], "k")
        while True:
            use = np.abs(y-y[j]) < 1e-4
            # plt.errorbar(idx[use], e[use], de[use], fmt=".")
            plt.plot(idx[use], e[use], ".")
            plt.plot(idx[use], 1+de[use], idx[use], 1-de[use], color="k", lw=1, alpha=0.5)

            # Fix `use` for next time
            use[:j] = True
            if np.all(use):
                break
            j = (~use).argmax()

        plt.xlabel("Scan step #")
        plt.ylabel("EDS rate relative to median")
        plt.title("SEM-side energy-dispersive sensor (EDS) rate")
        msg = "Horizontal lines show median EDS rate ± 1σ: {:.1f}±{:.1f} cps".format(
            scale, np.median(de)*scale)
        plt.text(.1, .1, msg, transform=ax.transAxes)

        plt.subplot(212, aspect="equal")
        x, y = self.target.T
        cm = plt.cm.inferno
        rate = self.eds/dur
        vmin, vmax = np.percentile(rate, [5, 95])
        plt.scatter(x, y, s=1400/self.Nx, marker="s", c=rate, vmin=vmin, vmax=vmax, cmap=cm)
        plt.xlabel("X location (mm)")
        plt.ylabel("Y location (mm)")
        plt.title("EDS rate vs scan position")
        plt.colorbar(location="bottom", label="EDS counts per second")

    def plot_targets(self, color=None, axis=None, showEnd=True):
        x, y = self.target.T
        if axis is None:
            plt.clf()
            axis = plt.subplot(111, aspect="equal")
            plt.xlim([x.min()-.0007, x.max()+.0007])
            plt.ylim([y.min()-.0007, y.max()+.0007])
            plt.title("Scan target locations: {}".format(self.basename))
            plt.xlabel("X location (mm)")
            plt.ylabel("Y location (mm)")

        if color is None:
            if "outer_1" in self.basename:
                color = "r"
            elif "outer_2" in self.basename:
                color = "orange"
            elif "inner_1" in self.basename:
                color = "b"
            elif "inner_2" in self.basename:
                color = "c"
            elif "inner_3" in self.basename:
                color = "g"
            else:
                color = "k"

        plt.plot(x, y, "s", color=color, ms=10, fillstyle="none")
        if showEnd:
            x, y = self.targetend.T
            plt.plot(x, y, "o", color=color, ms=4)


def findgrid(x, mindist=2e-4):
    """
    Given a vector `x`
    """
    vmin = np.median(x[np.abs(x-x.min()) < mindist])
    v = [vmin]
    while True:
        use = x > v[-1]+mindist
        if use.sum() == 0:
            break
        v2 = (x[use]).min()
        use = np.abs(x-v2) < mindist
        if use.sum() == 0:
            break
        v2 = np.median(x[use])
        v.append(v2)
    N = len(v)
    v = np.array(v)
    return N, v


def summarize_scan(scan):
    plt.clf()
    plt.gcf().suptitle("File: {}".format(scan.filename))
    ax = plt.subplot(411)
    ax.set_aspect("equal")
    plt.title("Scan targets and commands needed to hit them")
    plt.plot(scan.target[:, 0], scan.target[:, 1], "xk", label="Targets")
    plt.plot(scan.command[:, 0], scan.command[:, 1], "o-", color="r", ms=4, label="Commands")
    # plt.plot(scan.command[:,0], scan.command[:,1], "-", color="r")
    plt.xlabel("X location (mm)")
    plt.ylabel("Y location (mm)")
    plt.legend()

    dx = (scan.offset[1:, 0] - scan.offset[:-1, 0])*1e6
    dy = (scan.offset[1:, 1] - scan.offset[:-1, 1])*1e6
    d = np.sqrt(dx**2+dy**2)
    ax = plt.subplot(412)

    fx = np.poly1d(np.polyfit(scan.Xcoords, np.arange(scan.Nx)+0.5, 1))
    fy = np.poly1d(np.polyfit(scan.Ycoords, np.arange(scan.Ny)+0.5, 1))
    ix = np.asarray(fx(scan.target[:, 0]), dtype=int)
    iy = np.asarray(fy(scan.target[:, 1]), dtype=int)
    map = np.zeros((scan.Ny, scan.Nx), dtype=float)
    for i, (x, y) in enumerate(zip(dx, dy)):
        map[iy[i], ix[i]] = (x**2+y**2)**.5
    xstep = (scan.Xcoords[-1]-scan.Xcoords[0])/scan.Nx
    ystep = (scan.Ycoords[-1]-scan.Ycoords[0])/scan.Ny
    lft = scan.Xcoords[0]-.5*xstep
    rgt = scan.Xcoords[-1]+.5*xstep
    bot = scan.Ycoords[0]-.5*ystep
    top = scan.Ycoords[-1]+.5*ystep
    plt.subplot(412)
    plt.title("Drift during each step but the last")
    plt.imshow(map, extent=(lft, rgt, bot, top), origin="lower", aspect="equal", cmap=plt.cm.magma)
    for i, (x, y) in enumerate(zip(dx, dy)):
        tx, ty = scan.target[i, :2]
        plt.plot([tx, tx+x/1e6], [ty, ty+y/1e6], color="white")
        plt.plot(tx, ty, "o", color="w", ms=2)

    plt.subplot(413)
    plt.title("Change in offset with each step")
    plt.plot(dx, "-o", ms=3, label="x")
    plt.plot(dy, "-o", ms=3, label="y")
    plt.plot(np.zeros_like(dx), "k")
    plt.ylabel("Δoffset (nm)")
    plt.legend()

    dist = (dx**2+dy**2)**.5
    titles = ("X offset (nm)", "Y offset (nm)", "Offset (nm)")
    for i, d in enumerate((dx, dy, dist)):
        ax = plt.subplot(4, 4, 13+i)
        if d.min() >= 0:
            hrange = [0, 500]
        else:
            hrange = [-250, 250]
        plt.hist(d, 100, hrange, histtype="step", color="C{}".format(i))
        rms = (d**2).mean()**.5
        mad = np.median(np.abs(d))
        text = "RMS: {:.2f} nm".format(rms)
        plt.text(.9, .9, text, ha="right", transform=ax.transAxes)
        text = "MAD: {:.2f} nm".format(mad)
        plt.text(.9, .8, text, ha="right", transform=ax.transAxes)
        plt.title(titles[i])
        plt.xlim(hrange)

    ax = plt.subplot(4, 4, 16)
    plt.title("Drift offsets")
    ax.set_aspect("equal")
    good = dist < 100
    xrange = dx.max()-dx.min()
    yrange = dy.max()-dy.min()
    if yrange > xrange and False:
        plt.plot(dy[good], dx[good], "ok", ms=3)
        plt.plot(dy[~good], dx[~good], "or", ms=3)
        plt.plot(.1, .4, marker=r'$\circlearrowleft$', ms=20, transform=ax.transAxes, color="g")
        plt.text(.03, .2, "Rotate for xy", color="g", transform=ax.transAxes)
        plt.xlabel("ΔY (nm)")
        plt.ylabel("ΔC (nm)")
    else:
        plt.plot(dx[good], dy[good], "ok", ms=3)
        plt.plot(dx[~good], dy[~good], "or", ms=3)
        plt.xlabel("ΔX (nm)")
        plt.ylabel("ΔY (nm)")
    plt.tight_layout()


def compute_radiograph(scan, parity=False, rotation=0.0, mag=6340, voxsize_nm=50.0, std_extent=True):
    """
    Convert a scan to its 2d radiograph
    """
    cosang = np.cos(scan.angle*np.pi/180.)
    magX = mag * cosang**2
    magY = mag * cosang
    print("Magnifications: {:.2f}, {:.2f}".format(magX, magY))
    voxsize_mm = voxsize_nm/1e6
    xnom = np.array(scan.target[:, 0])
    ynom = np.array(scan.target[:, 1])
    # Correct to middle
    dx = (scan.offset[1:, 0] - scan.offset[:-1, 0])
    dy = (scan.offset[1:, 1] - scan.offset[:-1, 1])
    xnom[:-1] += 0.5*dx
    ynom[:-1] += 0.5*dy
    Nstep = len(xnom)
    assert Nstep == len(ynom)

    if std_extent:
        xshift = np.tan(scan.angle*np.pi/180.)*.01
        x0 = -.067+xshift
        x9 = -.033+xshift
        y0, y9 = -.009, -.002
    else:
        # Pixelize a region that's (6.5 mm/mag) around the nominal positions
        border = 6.5
        x0 = xnom.min()-border/magX
        x9 = xnom.max()+border/magX
        y0 = ynom.min()-border/magY
        y9 = ynom.max()+border/magY
    nx = int((x9-x0)/voxsize_mm+0.5)
    ny = int((y9-y0)/voxsize_mm+0.5)
    dx = nx*voxsize_mm
    dy = ny*voxsize_mm
    x9 = x0+dx
    y9 = y0+dy
    print("Making a ({:}x{:}) map with dx={:.2f} µm dy={:.2f} µm".format(nx, ny, dx*1e3, dy*1e3))
    if nx > 2000 or ny > 1000:
        raise ValueError("too many bins: {}x{}".format(nx, ny))

    Wt = np.zeros((nx, ny), dtype=float)
    Vt = np.zeros((nx, ny), dtype=float)
    Ct = np.zeros((nx, ny), dtype=float)
    xidx, yidx = [a for a in np.indices((nx+1, ny+1))]
    xedges = voxsize_mm*xidx+x0  # In mm
    yedges = voxsize_mm*yidx+y0  # In mm
    xpix = voxsize_mm*0.5+xedges[:-1]
    ypix = voxsize_mm*0.5+yedges[:-1]
    xpixctr = voxsize_mm*np.linspace(.5, nx-.5, nx)+x0
    ypixctr = voxsize_mm*np.linspace(.5, ny-.5, ny)+y0

    f = scan.h5
    sintheta = np.sin(rotation/57.2958)
    costheta = np.cos(rotation/57.2958)
    sig = 50
    sig2 = sig**2
    plt.figure(2)
    plt.clf()
    ax = plt.subplot(111)
    ax.set_aspect("equal")
    plt.title("Pixel positions on sample (µm) $\\theta={:.1f}^\\circ$".format(scan.angle))
    total_photons = 0.0
    for channum in range(1, 241):
        k = "chan{:03d}".format(channum)
        if k not in scan.h5["tes"]:
            continue
        mr = scan.medianrates[channum]
        if abs(mr/scan.medianrate-1.0) > 0.3:
            print("Skipping channel {} with mr={:.2f} standard={:.2f}".format(
                channum, mr, scan.medianrate))
            continue
        print(k)
        g = f["tes/{}".format(k)]
        try:
            X = g["xdet"][0, 0]/1e3  # In mm
            Y = g["ydet"][0, 0]/1e3
            X, Y = costheta*X-sintheta*Y, sintheta*X+costheta*Y
            if parity:
                X = -X
        except IndexError:
            continue
        plt.plot(X/magX*1000, Y/magY*1000, "or")
        # continue
        xtes = xnom + X/magX  # In mm
        ytes = ynom + Y/magY
        for idx in range(Nstep):
            this_signal = g["sig"][idx, 0]/scan.eds[idx]
            if this_signal <= 0:
                # print(k, idx, this_signal)
                continue
            this_unc = g["sig_unc"][idx, 0]/scan.eds[idx]

            ix = int((xtes[idx]-x0)/voxsize_mm)
            i0 = ix-2
            i9 = ix+3
            if i9 < 0:
                continue
            if i0 < 0:
                i0 = 0
            if i9 >= nx+1:
                i9 = nx
            iy = int((ytes[idx]-y0)/voxsize_mm)
            j0 = iy-2
            j9 = iy+3
            if j9 < 0:
                continue
            if j0 < 0:
                j0 = 0
            if j9 >= ny+1:
                j9 = ny
            # print(xtes[idx], xpix[ix, iy], ytes[idx], ypix[ix, iy])
            # print(xtes[idx]-xpix[ix, iy], ytes[idx]-ypix[ix, iy])
            # print((xtes[idx]-x0)/voxsize_mm, (ytes[idx]-y0)/voxsize_mm)
            # print("i,j: ", ix, i0, i9, iy, j0, j9)
            pixd2 = (xtes[idx]-xpix[i0:i9, j0:j9])**2+(ytes[idx]-ypix[i0:i9, j0:j9])**2
            wt_spot = np.exp(-0.5*(pixd2*1e12/sig2))
            wt_sn = this_unc**-2
            weight = wt_spot*wt_sn
            Vt[i0:i9, j0:j9] += weight*this_unc
            Wt[i0:i9, j0:j9] += weight
            Ct[i0:i9, j0:j9] += weight*this_signal
            total_photons += this_signal
    Wt[Wt <= 0] = Wt.min()/10
    rg = radiograph.Radiograph(Ct/Wt, Vt/Wt, Wt, xpixctr, ypixctr,
                               (xedges, yedges), scan.basename, scan.angle)
    rg.total_photons = total_photons
    return rg
