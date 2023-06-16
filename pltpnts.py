import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

#====================================================
# INPUT PARAMETERS
#====================================================
scale = 5
area = np.pi*2  #controls the size of the particle
var = "phi"
#====================================================

path = "pics/"+var+"_jpg"
try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed" % path)
else:
    print ("Successfully created the directory %s " % path)

#====================================================
# SETTING FIGURE SIZE AND LIMITS
#====================================================    
pltctrl = 'out/gmtry/pltctrl.txt'
pltctrltype=[int, int, float, float, float, float]
pltctrlnames = ["A", "B", "C", "D", "E", "F"]
ary1 = np.genfromtxt(pltctrl, names=pltctrlnames, dtype=pltctrltype)
print(ary1)

nt = ary1['A']
inc = ary1['B']
num = inc
xmin = ary1['C']
xmax = ary1['D']
ymin = ary1['E']
ymax = ary1['F']
L = 1.25*scale*(xmax-xmin)
W = scale*(ymax-ymin)
figure(figsize=(L, W), dpi=80)

#====================================================
# PLOTTING DATA AND EXPORTING TO JPG
#====================================================
n = 1
while (num <= nt):

    if num < 10:
        fname = 'out/'+var+'_txt/'+var+'00'+str(num)+'.txt'
    if num < 100 and num >= 10:
        fname = 'out/'+var+'_txt/'+var+'0'+str(num)+'.txt'
    if num >= 100:
        fname = 'out/'+var+'_txt/'+var+str(num)+'.txt'

    ndtype=[float, float, float, float]
    names = ["A", "B", "V", "D"]
    ary = np.genfromtxt(fname, names=names, dtype=ndtype)
    x = ary['A']
    y = ary['B']
    val = ary['D']
    Nx = np.size(ary['A'])

    plt.scatter(x, y, s=area, c=val, alpha=1, cmap='jet');
    plt.title(var)
    plt.colorbar()

    plt.savefig('pics/'+var+'_jpg/'+var+str(n)+".jpg")
    #if n < 10:
    #    plt.savefig('pics/'+var+'_jpg/'+var+'0000'+str(n)+".jpg")
    #if n < 100 and n >= 10:
    #    plt.savefig('pics/'+var+'_jpg/'+var+'000'+str(n)+".jpg")
    #if n < 1000 and n >= 100:
    #    plt.savefig('pics/'+var+'_jpg/'+var+'00'+str(n)+".jpg")
    #if n < 10000 and n >= 1000:
    #    plt.savefig('pics/'+var+'_jpg/'+var+'0'+str(n)+".jpg")
    #if n >= 10000:
    #    plt.savefig('pics/'+var+'_jpg/'+var+str(n)+".jpg")
    plt.clf()
    n = n + 1
    num = num + inc





