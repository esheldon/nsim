import os
import glob
from distutils.core import setup

scripts=['nsim-average-outputs',
         'nsim-combine-averaged',
         'nsim-average-runs',
         'nsim-combine-trials',
         'nsim-make-condor',
         'nsim-make-wq',
         'nsim-plot',
         'nsim-make-paper-plots',
         'nsim-run']


scripts=[os.path.join('bin',s) for s in scripts]

conf_files=glob.glob('config/*.yaml')

data_files=[]
for f in conf_files:
    data_files.append( ('share/nsim_config',[f]) )


setup(name="nsim", 
      version="0.1.0",
      description="Run shear simulations using ngmix",
      license = "GPL",
      author="Erin Scott Sheldon",
      author_email="erin.sheldon@gmail.com",
      scripts=scripts,
      data_files=data_files,
      packages=['nsim'])
