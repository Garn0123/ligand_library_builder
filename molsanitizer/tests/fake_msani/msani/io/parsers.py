import argparse


def parseArguments(args=None, batch_mode=False):
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", dest="input_files", nargs="+")
    ap.add_argument("--prefix")
    ap.add_argument("--no-neutralize", dest="neutralize", action="store_false", default=None)
    ap.add_argument("--no-stereoisomers", dest="stereoisomers", action="store_false", default=None)
    ap.add_argument("--metal")
    ap.add_argument("-3d", dest="gen3d", action="store_true")
    ap.add_argument("-f", dest="format", nargs="+")
    for f in ("--numconfs", "--energywindow", "--rmsd", "--nringconfs", "--timeout", "--randomSeed", "-j"):
        ap.add_argument(f)
    ap.add_argument("--timing", action="store_true")
    return ap.parse_args(args)


def Sanitycheck(args):
    # the real one turns stereo/neutralize ON with -3d unless told not to
    assert args.gen3d and args.stereoisomers is False and args.neutralize is False and args.metal == "off"
    return args
