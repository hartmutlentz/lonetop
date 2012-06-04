import os, sys

# add parent directory to PythonPath (here: src)
sys.path.insert(1, os.path.join(sys.path[0], os.pardir))


def twoUp():
    """returns path to directory two layers above"""
    p = os.getcwd()
    return p.rsplit("/", 2)[0]

def dataPath(fname):
    """
    returns string to <fname> in data directory.
    File <fname> does not necessarily have to exist.
    """
    p = twoUp()
    p = "/".join([p,"data",fname])
    return p

def plotPath(fname):
    """returns string to <fname> in plot directory"""
    p = twoUp()
    p = "/".join([p,"plot",fname])
    return p

def resultsPath(fname):
    """returns string to <fname> in results directory"""
    p = twoUp()
    p = "/".join([p,"results",fname])
    return p