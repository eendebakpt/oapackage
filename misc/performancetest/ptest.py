"""Script to test performance of components of the OA package

Log with:

    > python ptest.py > test0-ddmmyyy.txt
"""

# %% Load necessary packages
import getopt
import os
import platform
import subprocess
import sys
import time

import numpy as np

import oapackage
import oapackage.oahelper as oahelper
from oapackage import oalib

tmp = oalib.log_print(-oalib.SYSTEM, "")

# %% Define test functions


def pytest7(verbose=1):
    """Reduced array with random transformtion to LMC form"""
    al = oalib.exampleArray(3)
    adata = oalib.arraylink2arraydata(al)
    reduction = oalib.LMCreduction_t(adata)

    for ii in range(10):
        if verbose >= 2:
            print("pytest7: randomized reduction of array %d" % ii)
        reduction.transformation.randomize()
        al2 = reduction.transformation.apply(al)

        alr = al2.reduceLMC()
        c = al == alr
        if not c:
            print("pytest7: error: reduction of randomized array failed!")


# %%


def pytest6(verbose=1):
    """Check generation of OA(32, t, 2^a)"""
    N = 32
    k = 8
    t = 3
    l = [2] * 8
    rr = []
    gtrr = [3, 5, 10, 17, 33]
    s = oalib.intVector(l)
    adata = oalib.arraydata_t(s, N, t, k)

    if verbose:
        print("pytest6: run different algorithms on the same case")
    algs = [oalib.MODE_ORIGINAL, oalib.MODE_J4]
    for _, alg in enumerate(algs):
        algname = oalib.algnames(alg)
        if verbose >= 2:
            print(f"pytest6: running {adata.fullidstr()}, alg {algname}")
        rr = []
        oahelper.runExtend(N, k, t, l, verbose=verbose, nums=rr, algorithm=alg)
        if not rr == gtrr:
            print(f"pytest6: case {adata.fullidstr()}")
            print(f"   algorithm {algname}: error: incorrect number of arrays! {rr} -> {gtrr}")


# %%


def pytest3(verbose=2):
    """Test rank and Defficiency properties"""
    k = 6
    strength = 3
    N = 40
    lll = oahelper.runExtend(N, k, strength, 2, verbose=0)

    rnk0 = np.array([19, 19, 21, 21, 19, 22, 22, 22, 22])  # ground-truth
    D0 = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.1703475035645109,
            0.0,
            0.8107439985672059,
            0.8969911485707478,
            0.8953919248315042,
            0.9007069559712273,
        ]
    )  # gt

    rnk = [oalib.array2xf(al).rank() for al in lll]
    D = [al.Defficiency() for al in lll]

    if (rnk - rnk0).any():
        if verbose:
            print("  rank test FAILED!")
            return False
    if not (np.abs(D0 - D) < 1e-6).all():
        if verbose:
            print("  Defficiency test FAILED!")
        return False
    return True


def pytest4(verbose=1):
    """Check generation of OA(16,2, 4^2 2^a)"""
    N = 16
    k = 6
    t = 2
    l = [4, 4, 2]
    rr = []
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr)
    if not rr == [2, 7, 17, 27]:
        print(f"ERROR: incorrect number of arrays! {rr}")


def pytest2level(verbose=1):
    """Check generation of OA(16,2, 4^2 2^a)"""
    N = 32
    k = 10
    t = 3
    l = [2]
    rr = []
    t0 = time.time()
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr, algorithm=oalib.MODE_ORIGINAL)
    dt = time.time() - t0
    t0 = time.time()
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr, algorithm=oalib.MODE_LMC_2LEVEL)
    dt2level = time.time() - t0

    ll = l * k
    s = oalib.intVector(ll)

    adata = oalib.arraydata_t(s, N, t, k)

    if verbose:
        print(f"case {adata.idstr()}: 2-level method {dt:.2f} [s] -> {dt2level:.2f} [s]")
    # cases.append('runExtend(N=%d, k=%d, t=%d, l=%s)' %(N,k,t, str(l)) )


def pytest_analysis(oadir, verbose=1):
    """Performance of analysis algorithms"""
    cmd = "export OMP_NUM_THREADS=1; oaanalyse --gwp -A {}".format(os.path.join(oadir, "testdata", "test64.oa"))
    if verbose:
        print("pytest_analysis: running oaanalyse")
    t0 = time.time()
    return_value = os.system(cmd)
    assert return_value == 0
    ta = time.time() - t0
    return ta


def pytest2(oadir, verbose=1):
    afile = os.path.join(oadir, "testdata", "unstable-rank.oa")
    sols = oalib.readarrayfile(afile)
    afile = os.path.join(oadir, "testdata", "unstable-rank-extended.oa")
    sols = oalib.readarrayfile(afile)

    for ii, al in enumerate(sols):
        vv = oalib.doubleVector(3)
        rnk = oalib.array2rank_Deff_Beff(al, vv)

        al2 = oalib.array2xf(al)
        rnk2 = al2.rank()
        if (not rnk2 == rnk) or (not rnk < 29):
            print("warning: problem?!: %d %d %d" % (ii, rnk, rnk2))
        if verbose >= 2:
            print("array %d: rank: %d %d" % (ii, rnk, rnk2))


def pytest(verbose=1):
    t0 = time.time()
    tt = []
    cases = []

    N = 48
    k = 6
    t = 3
    l = 2
    rr = []
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr)
    tt.append(time.time() - t0)
    cases.append("runExtend(N=%d, k=%d, t=%d, l=%d)" % (N, k, t, l))
    if not rr == [4, 10, 45]:
        print(f"ERROR: incorrect number of arrays! {rr}")

    N = 32
    k = 6
    t = 3
    l = 2
    rr = []
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr)
    tt.append(time.time() - t0)
    cases.append("runExtend(N=%d, k=%d, t=%d, l=%d)" % (N, k, t, l))
    if not rr == [3, 5, 10]:
        print(f"ERROR: incorrect number of arrays! {rr}")

    N = 16
    k = 16
    t = 2
    l = 2
    rr = []
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr)
    tt.append(time.time() - t0)
    cases.append("runExtend(N=%d, k=%d, t=%d, l=%d)" % (N, k, t, l))
    if not rr == [3, 5, 11, 27, 55, 80, 87, 78, 58, 36, 18, 10, 5, 0]:
        print(f"ERROR: incorrect number of arrays! {rr}")

    N = 96
    k = 8
    t = 4
    l = 2
    rr = []
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr)
    tt.append(time.time() - t0)
    cases.append("runExtend(N=%d, k=%d, t=%d, l=%d)" % (N, k, t, l))
    if not rr == [4, 9, 4, 0]:
        print(f"ERROR: incorrect number of arrays! {rr}")

    dt = time.time() - t0
    if verbose:
        print(f"total time {dt:.3f} [s]")

    if verbose:
        print("extra testing cases")
    N = 24
    k = 5
    t = 2
    l = [3, 2]
    rr = []
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr)
    tt.append(time.time() - t0)
    cases.append("runExtend(N=%d, k=%d, t=%d, l=%s)" % (N, k, t, str(l)))
    if not rr == [4, 29, 573]:
        print(f"ERROR: incorrect number of arrays! {rr}")

    N = 18
    k = 9
    t = 2
    # old test: l=[2,3]
    l = [3] * 8 + [2]
    rr = []
    oahelper.runExtend(N, k, t, l, verbose=1, nums=rr)
    tt.append(time.time() - t0)
    cases.append("runExtend(N=%d, k=%d, t=%d, l=%s)" % (N, k, t, str(l)))
    if not rr == [4, 12, 10, 8, 3, 0, 0]:
        print(f"ERROR: incorrect number of arrays! {rr}")

    dt2 = time.time() - t0
    return (dt, tt, cases, dt2)


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


# %%


def testExtendBinary(verbose=1):
    repos = dict()
    repos["oatest1.txt"] = "result-16.2-2-2-2-2-2-2-2.oa"
    repos["oatest2.txt"] = "result-18.3-3-3-3-3.oa"
    repos["oatest3.txt"] = "result-8.4-2-1-4-2-2.oa"
    repos["oatest4.txt"] = "result-28.7-2-2-2.oa"
    repos["oatest5.txt"] = "result-16.2-2-2-2-2-2-2-2-2.oa"
    repos["oatest6.txt"] = "result-25.5-5-5-5.oa"
    repos["oatest7.txt"] = "result-64.2-2-2-2-2-2.oa"
    repos["oatest8.txt"] = "result-56.2-2-2-2-2.oa"


# testExtendBinary()


# %%
def main(argv=None):
    """Main testing function"""
    if platform.system() == "Windows":
        oadir = os.path.join(os.path.split(oapackage.__file__)[0], "..")
    else:
        oadir = os.path.join(os.path.expanduser("~"), "misc/oa/oacode/")

    print("OA performance testing: version 2.0")
    ss = oapackage.version()
    print(f"OAlib: version {ss}")
    ss = oalib.compile_information()
    print(ss)
    print("System:")
    print("  Python: " + " ".join(sys.version.split("\n")))
    print("  Machine: " + " ".join([f"{x}" for x in platform.uname()]))
    print(f"  Processor: {platform.processor()}")

    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        # more code, unchanged
    except Usage as err:
        print >> sys.stderr, err.msg
        print >> sys.stderr, "for help use --help"
        return 2

    print("\nRunning tests:\n")
    (dt, tt, cases, dt2) = pytest()

    pytest2(oadir)
    pytest3(verbose=0)
    pytest4(verbose=1)
    ta1 = pytest_analysis(oadir, verbose=1)
    pytest6(verbose=0)
    pytest7(verbose=0)

    pytest2level(verbose=1)

    print("-----")

    dtx = ta1

    print("\nResults:\n")
    for ii, t in enumerate(tt):
        print(f"{cases[ii]}: time {t:.2f} [s]")
    version = oapackage.__version__
    print(f"Total time: {dt:.3f} [s], {dt2:.3f} [s], {dtx:.3f} [s] ({platform.node()}, {version})")
    print("   should be of order 3.408 [s], 3.494 [s], 7.313 [s] (woelmuis, 2.4.8)")
    print("   should be of order 4.4 [s], 4.6 [s], 5.9 [s] (woelmuis)")
    print("   should be of order 3.2 [s], 3.3 [s], 4.7 [s] (marmot) [v 2.0.0]")
    print("   should be of order 2.96 [s], 3.0 [s], 4.6 [s] (marmot) [v 2.0.24]")


if __name__ == "__main__":
    main()


def timeconfig(
    configfile="oaconfig.txt", timebin="/usr/bin/time", pdir="/home/eendebakpt/misc/oa/oacode/performancetest"
):
    os.chdir(pdir)
    res = subprocess.check_output(
        [(f"cd {pdir};ls;") + timebin, '--format="%%E %%S %%U"', "./oaextendsingle", f"-c {configfile} -l 1"]
    )
    return res
