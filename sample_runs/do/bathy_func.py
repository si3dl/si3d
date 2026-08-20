import numpy as np
import pandas as pd
import sys
import os


# Load Si3D functions
if os.name == 'nt':
    root = "C:/Users/"
    name = os.getlogin()
    root += name
else:
    import pwd
    name = pwd.getpwuid(os.getuid()).pw_name
    root = '/Users/' + name

func = '/Documents/Github/si3dInputs/pythonlibrary/'
FuncPath = root + func
sys.path.append(FuncPath)
from bathy4si3d import bathy4si3d # type: ignore

cwd = os.path.dirname(os.path.abspath(__file__))

BasinType=2
dx = 10
L = 500
B = 500
H = 15

PathSave = cwd

dxsave = ' (dx= '+str(dx)+'),'
header = 'RectangularLak'
if len(header + dxsave) != 27:
    while len(header + dxsave) < 27:
        header += ' '
    while len(header + dxsave) > 27:
        header = header[0:-1]



[X,Y,Z] = bathy4si3d(BasinType,header,dx,PathSave,L,B,H)