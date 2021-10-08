import os
import pylab as plt


def needs_updating(target, dependencies):
    """
    Does the file `target` need to be updated, in light of the collection `dependencies`?

    Return True if target doesn't exist, has zero size, or is older than any dependency.
    Otherwise, return False.
    """
    if not os.path.exists(target):
        return True
    if os.path.getsize(target) <= 0:
        return True
    target_mtime = os.path.getmtime(target)
    for d in dependencies:
        if target_mtime < os.path.getmtime(d):
            return True
    return False


def myfigure(size, fignum=9):
    """
    Generate a pylab figure number `fignum` to replace the current figure.

    Do so only if the requested size tuple (width,height) is different by more
    than 0.01 inches in either dimension from the current figure's size. If the
    figure exists and is the same size as `size`, then no need to close it.
    """
    curr_size = plt.gcf().get_size_inches()
    if abs(curr_size[0] - size[0]) + abs(curr_size[1] - size[1]) > .01:
        plt.close(fignum)
    return plt.figure(fignum, figsize=size)
