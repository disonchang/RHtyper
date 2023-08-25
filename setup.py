import sys, setuptools
from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()
#logging not added
install_reqPy3=['numpy','pandas','scipy','pysam','matplotlib','fpdf','plotly',"xlsxwriter"]
install_reqPy2=['numpy<=1.16.6','pandas<=0.24.0','scipy','pysam<=0.15.3','matplotlib<=2.1.0','fpdf==1.7.2','plotly<=2.1.0',"xlsxwriter"]

if sys.version_info >= (3,0):
    install_requires = install_reqPy3
else:
    install_requires = install_reqPy2



setup(name='RHtyper',
      version='0.1',
      description='RH typing using whole genome sequencing data',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Biodata :: Polymorphism',
      ],
      #url='http://github.com/storborg/funniest',
      keywords='RH typing',
      author='Ti-Cheng Chang',
      author_email='ti-cheng.chang@stjude.org',
      license='MIT',
      #packages=['RHtyper'],
      packages=setuptools.find_packages(),
      install_requires=install_requires,
      scripts=['bin/RHtyper'],
      #entry_points = {
      #  'console_scripts': ['RHtyper=RHtyper.RHtyper:main'],
      #},
      package_data={'RHtyper': ['BGClf/*.jpg','BGClf/*.png', 'BGClf/coordinates/*.*', 'BGClf/database/*.*']},
      include_package_data=True,
      zip_safe=False)



