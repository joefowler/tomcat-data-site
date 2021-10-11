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
        # self.edges =

    def copy(self):
        return Radiograph(self.rg, self.uncert, self.wt, self.xctr, self.yctr, self.xyedges,
                          self.basename, self.angle)

    def plot(self):
        plt.clf()
        plt.title("Radiograph of {}".format(self.basename))
        ax = plt.subplot(111)
        xpix, ypix = self.xyedges
        plt.pcolor(xpix, ypix, self.rg, vmin=.0014, vmax=.003, cmap=plt.cm.inferno)
        ax.set_aspect("equal")
        plt.colorbar(fraction=.07, location="bottom")
        plt.xlabel("X location (mm)")
        plt.ylabel("Y location (mm)")

    def smoother(self):
        x = self.xyedges[0][:-1, :-1].ravel()+self.xpixsize*0.5
        y = self.xyedges[1][:-1, :-1].ravel()+self.ypixsize*0.5
        weight = 1.0/self.uncert.ravel()
        weight[np.isnan(weight)] = 0
        weight[weight < 0] = 0
        poswgt = weight > 0
        raise NotImplementedError("This smoother isn't working yet, I think")
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
    chi2x = []
    nsteps = 7
    allxoff = np.arange(xoffctr-nsteps, xoffctr+nsteps+1)
    for xoff in allxoff:
        yoff = yoffctr
        mg, mh, vg, vh = rgs_var(g, h, xoff, yoff)
        chisq = (mg-mh)**2/(vg+vh)
        chisq[np.isnan(chisq)] = 0
        chi2x.append(chisq[chisq > 0].mean())

    chi2y = []
    allyoff = np.arange(yoffctr-nsteps, yoffctr+nsteps+1)
    for yoff in allyoff:
        xoff = xoffctr
        mg, mh, vg, vh = rgs_var(g, h, xoff, yoff)
        chisq = (mg-mh)**2/(vg+vh)
        chisq[np.isnan(chisq)] = 0
        chi2y.append(chisq[chisq > 0].mean())

    plt.clf()
    plt.subplot(211)
    plt.plot(allxoff, chi2x, "o-")
    f = np.poly1d(np.polyfit(allxoff[nsteps-3:nsteps+4], chi2x[nsteps-3:nsteps+4], 2))
    plt.plot(allxoff, f(allxoff), "--")
    plt.plot(allxoff[nsteps-3:nsteps+4], f(allxoff[nsteps-3:nsteps+4]), color="C1")
    plt.xlabel("Pixel shift in x")
    plt.ylabel("<$\\chi^2$> non-empty pix")
    xoff = f.deriv(1).roots[0]

    plt.subplot(212)
    plt.plot(allyoff, chi2y, "o-")
    f = np.poly1d(np.polyfit(allyoff[nsteps-3:nsteps+4], chi2y[nsteps-3:nsteps+4], 2))
    plt.plot(allyoff, f(allyoff), "--")
    plt.plot(allyoff[nsteps-3:nsteps+4], f(allyoff[nsteps-3:nsteps+4]), color="C1")
    plt.xlabel("Pixel shift in y")
    plt.ylabel("<$\\chi^2$> non-empty pix")
    yoff = f.deriv(1).roots[0]

    print("Best x offset: {:7.3f} pix = {:8.1f} nm".format(xoff, xoff*g.xpixsize*1e6))
    print("Best y offset: {:7.3f} pix = {:8.1f} nm".format(yoff, yoff*g.ypixsize*1e6))
