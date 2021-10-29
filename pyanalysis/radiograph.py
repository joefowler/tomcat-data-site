import numpy as np
import scipy as sp
import pylab as plt


class Radiograph:
    def __init__(self, rg, uncertainty, weight, xctr, yctr, xyedges, basename, angle):
        self.rg = rg
        self.uncert = uncertainty
        self.wt = weight
        self.xctr = xctr
        self.yctr = yctr
        self.xyedges = xyedges
        self.basename = basename
        self.angle = angle
        self.xpixsize = np.median(np.diff(xctr))
        self.ypixsize = np.median(np.diff(yctr))
        self.Nx, self.Ny = rg.shape

    def copy(self):
        return Radiograph(self.rg, self.uncert, self.wt, self.xctr, self.yctr, self.xyedges,
                          self.basename, self.angle)

    def _combine(self, other, shiftx=0, shifty=0, mode="add"):
        if self.angle != other.angle:
            raise ValueError("Cannot combine radiographs from unequal angles: {}, {}".format(
                self.angle, other.angle))

        def shiftb(B):
            B = B.copy()
            if shiftx > 0:
                B[:-shiftx, :] = B[shiftx:, :]
                B[-shiftx:, :] = 0
            elif shiftx < 0:
                B[-shiftx:, :] = B[:shiftx, :]
                B[:-shiftx, :] = 0
            if shifty > 0:
                B[:, :-shifty] = B[:, shifty:]
                B[:, -shifty] = 0
            elif shifty < 0:
                B[:, -shifty:] = B[:, :shifty]
                B[:, :-shifty] = 0
            return B

        def wtmean(A, wa, B, wb):
            A = A.copy()
            B = shiftb(B)
            if not np.isscalar(wb):
                wb = shiftb(wb)
            A[np.isnan(A) & ~np.isnan(B)] = 0
            B[np.isnan(B) & ~np.isnan(A)] = 0
            return (A*wa+B*wb)/(wa+wb)

        def add(A, B):
            return A + shiftb(B)

        if mode in ("average", ):
            w1 = self.uncert**-2
            w2 = other.uncert**-2
            w1[np.isnan(w1)] = 0
            w2[np.isnan(w2)] = 0

            rg = wtmean(self.rg, w1, other.rg, w2)
            wt = add(self.wt, other.wt)
            uncert = add(w1, w2)**-0.5
            uncert[np.isinf(uncert)] = np.nan
            basename = "Avg ({}, {})".format(self.basename, other.basename)
        elif mode in ("subtract", "difference"):
            rg = add(self.rg, -other.rg)
            wt = add(self.wt, other.wt)
            uncert = add(self.uncert**2, other.uncert**2)**0.5
            basename = "Diff ({}, {})".format(self.basename, other.basename)
        else:
            raise ValueError("Radiograph._combine() requires mode='average' or mode='subtract'")
        return Radiograph(rg, uncert, wt, self.xctr, self.yctr, self.xyedges, basename,
                          self.angle)

    def average(self, other, shiftx=0, shifty=0):
        return self._combine(other, shiftx, shifty, "average")

    def subtract(self, other, shiftx=0, shifty=0):
        return self._combine(other, shiftx, shifty, "subtract")

    def plot(self, name="rg", vmin=None, vmax=None, midpct=95):
        plt.clf()
        plt.subplot(111, aspect="equal")
        plt.title("Radiograph of {}".format(self.basename))
        xpix, ypix = self.xyedges
        z = self.rg

        if name == "rg":
            z = self.rg
            label = "Radiograph value"
        elif name == "wt":
            z = self.wt
            label = "Weights"
        elif name.startswith("un"):
            z = self.uncert
            label = "Uncertainty"
        elif name.startswith("sig"):
            z = self.rg/self.uncert
            label = "Signal to noise ratio"

        # If the vmin and/or vmax aren't given, we need to choose them.
        # Instead of the default (set equal to the lowest and highest z values),
        # We'll set the scale so that the middle `midpct` values use the
        # middle `midpct` of the color scale.
        if midpct > 99:
            midpct = 99
        elif midpct < 1:
            midpct = 1
        plo, phi = np.nanpercentile(z.ravel(), [50-midpct/2, 50+midpct/2])
        dp = phi-plo
        if vmin is None:
            vmin = plo-dp*5/midpct
        if vmax is None:
            vmax = phi+dp*5/midpct
        plt.pcolor(xpix, -ypix, z, vmin=vmin, vmax=vmax, cmap=plt.cm.inferno)
        plt.colorbar(fraction=.07, location="bottom", label=label)
        plt.xlabel("X location (mm)")
        plt.ylabel("$-$Y location (mm)")

    def smoother(self):
        raise NotImplementedError("This smoother isn't working yet, I think")
        x = self.xyedges[0][:-1, :-1].ravel()+self.xpixsize*0.5
        y = self.xyedges[1][:-1, :-1].ravel()+self.ypixsize*0.5
        weight = 1.0/self.uncert.ravel()
        weight[np.isnan(weight)] = 0
        weight[weight < 0] = 0
        poswgt = weight > 0
        return sp.interpolate.SmoothBivariateSpline(x[poswgt], y[poswgt], self.rg.ravel()[poswgt], weight[poswgt])

    def spline(self):
        "Return a bivariate spline approximation to the map"
        z = self.rg.copy()
        z[np.isnan(z)] = 0.002
        return sp.interpolate.RectBivariateSpline(self.xctr, self.yctr, z)


def rgs_var(g, h, xoff, yoff):
    """Radiographs and variance with xoff and yoff offsets"""
    gx0, hx0 = xoff, 0
    if xoff < 0:
        gx0, hx0 = 0, -xoff
    nx = min(g.rg.shape[0]-gx0, h.rg.shape[0]-hx0)
    gy0, hy0 = yoff, 0
    if yoff < 0:
        gy0, hy0 = 0, -yoff
    ny = min(g.rg.shape[1]-gy0, h.rg.shape[1]-hy0)
    mg = g.rg[gx0:gx0+nx, gy0:gy0+ny]
    mh = h.rg[hx0:hx0+nx, hy0:hy0+ny]
    vg = g.uncert[gx0:gx0+nx, gy0:gy0+ny]**2
    vh = h.uncert[hx0:hx0+nx, hy0:hy0+ny]**2
    return mg, mh, vg, vh


def compare_graphs(g, h, xoffctr=0, yoffctr=0):
    chi2x, chi2y = [], []
    nsteps = 7

    mg, vg = g.rg, g.uncert**2
    mh = h.spline()
    z = h.uncert**2
    z[np.isnan(z)] = np.nanmedian(z)
    vh = sp.interpolate.RectBivariateSpline(h.xctr, h.yctr, z)
    x = g.xctr
    y = g.yctr
    print("Try one now")
    print(x.shape, y.shape)
    print(mh(x, y).shape)

    allxoff = np.linspace(xoffctr-nsteps, xoffctr+nsteps, 5*nsteps+1)
    for xoff in allxoff:
        yoff = yoffctr
        # mg, mh, vg, vh = rgs_var(g, h, xoff, yoff)
        # chisq = (mg-mh)**2/(vg+vh)
        chisq = (mg-mh(x, y+xoff*g.xpixsize))**2/(vg+vh(x, y+yoff*g.ypixsize))
        chisq[np.isnan(chisq)] = 0
        chi2x.append(chisq[chisq > 0].mean())
    chi2x = np.array(chi2x)

    allyoff = np.linspace(yoffctr-nsteps, yoffctr+nsteps, 4*nsteps+1)
    for yoff in allyoff:
        xoff = xoffctr
        # mg, mh, vg, vh = rgs_var(g, h, xoff, yoff)
        # chisq = (mg-mh)**2/(vg+vh)
        chisq = (mg-mh(x, y+xoff*g.xpixsize))**2/(vg+vh(x, y+yoff*g.ypixsize))
        chisq[np.isnan(chisq)] = 0
        chi2y.append(chisq[chisq > 0].mean())
    chi2y = np.array(chi2y)

    plt.clf()
    plt.subplot(211)
    plt.plot(allxoff, chi2x, "o-")
    fituse = np.abs(allxoff-xoffctr) <= 2
    f = np.poly1d(np.polyfit(allxoff[fituse], chi2x[fituse], 2))
    plt.plot(allxoff, f(allxoff), "--")
    plt.plot(allxoff[fituse], f(allxoff[fituse]), color="C1")
    plt.xlabel("Pixel shift in x")
    plt.ylabel("<$\\chi^2$> non-empty pix")
    xoff = f.deriv(1).roots[0]

    plt.subplot(212)
    plt.plot(allyoff, chi2y, "o-")
    fituse = np.abs(allyoff-yoffctr) <= 2
    f = np.poly1d(np.polyfit(allyoff[fituse], chi2y[fituse], 2))
    plt.plot(allyoff, f(allyoff), "--")
    plt.plot(allyoff[fituse], f(allyoff[fituse]), color="C1")
    plt.xlabel("Pixel shift in y")
    plt.ylabel("<$\\chi^2$> non-empty pix")
    yoff = f.deriv(1).roots[0]

    print("Best x offset: {:7.3f} pix = {:8.1f} nm".format(xoff, xoff*g.xpixsize*1e6))
    print("Best y offset: {:7.3f} pix = {:8.1f} nm".format(yoff, yoff*g.ypixsize*1e6))
