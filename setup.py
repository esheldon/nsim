import os
import glob
from distutils.core import setup

scripts=['nsim-average-outputs',
         'nsim-average-runs',
         'nsim-combine-trials',
         'nsim-make-condor',
         'nsim-plot',
         'nsim-run']


scripts=[os.path.join('bin',s) for s in scripts]

conf_files=glob.glob('config/*.yaml')

data_files=[]
for f in conf_files:
    d=os.path.dirname(f)
    n=os.path.basename(f)
    d=os.path.join('share',d)

    data_files.append( (d,[f]) )


setup(name="nsim", 
      version="0.1.0",
      description="Run shear simulations using ngmix",
      license = "GPL",
      author="Erin Scott Sheldon",
      author_email="erin.sheldon@gmail.com",
      scripts=scripts,
      data_files=data_files,
      packages=['nsim'])
