import collections
import datetime


def _flatten(x):
    """Replace deprecated ``compiler.ast.flatten``"""
    for e in x:
        if not isinstance(e, collections.abc.Iterable) or isinstance(e, str):
            yield e
        else:
            yield from _flatten(e)


def logger(log, logstring):
    with open(log, "a") as ofile:
        ofile.write(str(datetime.datetime.now()) + " " + logstring)