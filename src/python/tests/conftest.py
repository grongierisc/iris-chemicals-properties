import sys
from os.path import dirname as d
from os.path import abspath, join
import os
os.environ['IRISINSTALLDIR'] = '/opt/intersystems/iris'
os.environ['IRISUSERNAME'] = 'SuperUser'
os.environ['IRISPASSWORD'] = 'SYS'
os.environ['IRISNAMESPACE'] = 'IRISAPP'
root_dir = d(d(abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(join(root_dir, 'opm'))