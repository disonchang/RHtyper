#__all__=['BGClf','RHtyper']

from os.path import dirname, basename, isfile, join
import glob, os, sys, pkgutil

if sys.version_info >= (3,0):
    __all__ = list(module for _, module, _ in pkgutil.iter_modules([os.path.dirname(__file__)]))
else:
    modules = glob.glob(join(dirname(__file__), "*.py"))
    __all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]
#print(__all__)
