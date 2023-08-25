### import all modules at the same time
#__all__=["coordinates","database"]
#from os.path import dirname, basename, isfile
#import glob, sys

### python 2
#modules = glob.glob(dirname(__file__)+"/*.py")
#__all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]


#if (sys.version_info > (3, 0)):
#    ### python 3
#    import os, pkgutil; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
#    __all__ = list(module for _, module, _ in pkgutil.iter_modules([os.path.dirname(__file__)]))



from os.path import dirname, basename, isfile, join
import sys, glob, os, pkgutil

if sys.version_info >= (3,0):
    sys.path.append(os.path.dirname(os.path.realpath(__file__)))
    __all__ = list(module for _, module, _ in pkgutil.iter_modules([os.path.dirname(__file__)]))
else:
    modules = glob.glob(join(dirname(__file__), "*.py"))
    __all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]
#print(__all__)
