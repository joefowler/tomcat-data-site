#!/usr/bin/env python
import numpy as np
import pylab as plt


def makeplots():
    plt.figure(2, (4.5, 1.5))
    for angle in np.linspace(-45, 45, 13):
        print("Making plot at angle {:+5.1f} degrees".format(angle))
        makeplot(angle)
        plt.tight_layout()
        plt.savefig("../plots/cartoons/angle{:.0f}.png".format(angle))


def makeplot(angle):
    c = np.cos(angle*np.pi/180)
    s = np.sin(angle*np.pi/180)
    plt.clf()
    plt.subplot(111, aspect="equal")
    plt.title("Top view, $\\theta={:.1f}^\\circ$".format(angle))

    # electrons
    for y in np.linspace(-16, 16, 7):
        plt.plot([-60, 0], [y, 0], color="gray", lw=1)
    plt.text(-50, 0, "SEM", va="center")
    # x rays
    for y in np.linspace(-60, 60, 49):
        if np.abs(y) < 22.5:
            plt.plot([0, 357], [0, y], color="gold", lw=.5)
        else:
            plt.plot([0, 350], [0, y], ":", color="gold", lw=.5)

    # TES circle
    phi = np.linspace(0, 2*np.pi, 101)
    x = np.cos(phi)*12+350
    y = np.sin(phi)*22.5
    plt.plot(x, y, "k")
    y = np.array([22.5, 35, 60])
    for a in (y, -y):
        plt.plot([350, 350, 400], a,  "k")

    targety = np.array([-45, 45])
    targetx = targety*s/c
    plt.fill_betweenx(targety, targetx, targetx+8*c, color="yellow", alpha=.1)
    plt.fill_betweenx(targety, targetx+8*c, targetx+11*c, color="r", alpha=.4)
    plt.plot(targetx, targety, "k")

    if angle > 0:
        plt.text(-5+25*s/c, 25, "Pt", color="k", ha="right")
        plt.text(20+25*s/c, 25, "Circuit", color="r", ha="left")
        plt.text(20+25*s/c, -35, "X rays", color="gold", ha="left")
    else:
        plt.text(-5-30*s/c, -30, "Pt", color="k", ha="right", va="top")
        plt.text(20-30*s/c, -30, "Circuit", color="r", ha="left", va="top")
        plt.text(20-33*s/c, 25, "X rays", color="gold", ha="left")

    plt.text(370, 0, "TES", rotation=90, va="center")
    plt.plot(0, 0, "o", color="gold", ms=4)
    plt.xlim([-50, 385])
    plt.ylim([-45, 45])
    plt.xticks(())
    plt.yticks(())


if __name__ == "__main__":
    makeplots()
