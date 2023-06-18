import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

#-------------------------------------------------------
# INPUT PARAMETERS
#-------------------------------------------------------
scale = float(sys.argv[1]) #parameter that scales the plot
show_mesh = int(sys.argv[2])

v = 3
while (v < len(sys.argv)):
    var = sys.argv[v]

    #-------------------------------------------------------
    # CREATING A PATH DIRECTORY TO EXPORT PICS
    #-------------------------------------------------------
    path = "out/"+var+"_jpg"
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

    #-------------------------------------------------------
    # SETTING FIGURE SIZE AND LIMITS
    #-------------------------------------------------------
    pltctrl = 'out/gmtry/pltctrl.txt'
    pltctrltype = [int, int, float, float, float, float]
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
    L = 1.25 * scale * abs(xmax - xmin)
    W = scale * abs(ymax - ymin)
    figure(figsize=(L, W), dpi=80)

    #-------------------------------------------------------
    # READING ELEMENTS
    #-------------------------------------------------------
    fname = 'out/gmtry/elems.txt'
    ndtype = [int, int, int]
    names = ["A", "B", "C"]
    ary = np.genfromtxt(fname, names=names, dtype=ndtype)
    p1 = ary['A']
    p2 = ary['B']
    p3 = ary['C']
    elems = np.vstack((p1, p2, p3)).T

    #-------------------------------------------------------
    # PLOTTING DATA AND EXPORTING TO JPG
    #-------------------------------------------------------
    n = 1
    while (num <= nt):

        if num < 10:
            fname = 'out/'+var+'_txt/'+var+'00'+str(num)+'.txt'
        if num < 100 and num >= 10:
            fname = 'out/'+var+'_txt/'+var+'0'+str(num)+'.txt'
        if num >= 100:
            fname = 'out/'+var+'_txt/'+var+str(num)+'.txt'

        ndtype = [float, float, float, float]
        names = ["A", "B", "V", "D"]
        ary = np.genfromtxt(fname, names=names, dtype=ndtype)
        x = ary['A']
        y = ary['B']
        val = ary['D']
        Nx = np.size(ary['A'])

        if show_mesh == 1:
            for element in elems:
                px = [x[element[i]] for i in range(len(element))]
                py = [y[element[i]] for i in range(len(element))]
                plt.fill(px, py, edgecolor='black', fill=False)

        tcf = plt.tricontourf(x, y, elems, val, cmap='jet')
        #tcf = plt.tricontour(x, y, elems, val, cmap='RdBu_r', levels=[0.5])
        plt.colorbar()
        plt.title(var)

        # OPTION 2.
        #fig = plt.figure(figsize=plt.figaspect(1))
        #ax = fig.add_subplot(projection='3d')
        #ax.plot_trisurf(x, y, val, triangles=elems, cmap='RdBu_r')
        #ax.view_init(90, -90)

        plt.savefig('out/'+var+'_jpg/'+var+str(n)+".jpg")

        #if n < 10:
        #    plt.savefig('pictures/'+var+'_jpg/'+var+'0000'+str(n)+".jpg")
        #if n < 100 and n >= 10:
        #    plt.savefig('pictures/'+var+'_jpg/'+var+'000'+str(n)+".jpg")
        #if n < 1000 and n >= 100:
        #    plt.savefig('pictures/'+var+'_jpg/'+var+'00'+str(n)+".jpg")
        #if n < 10000 and n >= 1000:
        #    plt.savefig('pictures/'+var+'_jpg/'+var+'0'+str(n)+".jpg")
        #if n >= 10000:
        #    plt.savefig('pictures/'+var+'_jpg/'+var+str(n)+".jpg")
        plt.clf()
        
        n = n + 1
        num = num + inc
    v = v + 1





