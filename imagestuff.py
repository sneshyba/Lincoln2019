# Functions for reading in new data set
from scipy import misc
import os
import copy
import numpy as np
from PIL import Image, ImageFont, ImageDraw
from scipy.interpolate import griddata


def getval2(Filename): #gets pixel values for a single bmp image
    value = misc.imread(Filename)
    Nx,Nz = value.shape
    return value, Nx, Nz, Filename
    
def getc2(folder,namebase,imageroot,mydet='D'): #applies getval to all four images in a set
    detectors = 'A', 'B', 'C', 'D'
    for det in detectors:
        Filename = folder+namebase+imageroot + '-' + det + '.bmp'
        print(Filename)
        value, Nx, Ny, Filename = getval2(Filename)
        if det == 'A':
            cA = value #raw bmp data
            
        if det == 'B':
            cB = value
            
        if det == 'C':
            cC = value
                
        if det == 'D':
            cD = value   

    #Getting pixel size
    textfile = folder+namebase+imageroot+'-A'+'.txt'
    fileref = open(textfile,"r")
    dummy = fileref.read().split()
    for i in range(len(dummy)):
        test = dummy[i].find('PixelSize')
        #print (i, test)
        if test!=-1:
            pixelsize = float(dummy[i][10:])/1000.
            break
    #pixelsize = float(dummy[16][10:])/1000 #microns
    fileref.close()
    dx = dy = pixelsize
    Filename = folder+namebase+imageroot + '-' + mydet + '.bmp'
    
    return dx,dy,cA,cB,cC,cD,Filename

def getc2tif(folder,namebase,imageroot): #Same as getc2, but picking up .tif files instead of .bmp files
    detectors = 'A', 'B', 'C', 'D'
    for det in detectors:
        Filename = folder+namebase+imageroot + '-' + det + '.tif'
        #print(Filename)
        value, Nx, Ny, Filename = getval2(Filename)
        if det == 'A':
            cA = value #raw tif data
            
        if det == 'B':
            cB = value
            
        if det == 'C':
            cC = value
                
        if det == 'D':
            cD = value   

    #Getting pixel size
    textfile = folder+namebase+imageroot+'-A'+'.txt'
    fileref = open(textfile,"r")
    dummy = fileref.read().split()
    for i in range(len(dummy)):
        test = dummy[i].find('PixelSize')
        #print (i, test)
        if test!=-1:
            pixelsize = float(dummy[i][10:])/1000.
            break
    #pixelsize = float(dummy[16][10:])/1000 #microns
    fileref.close()
    dx = dy = pixelsize
    
    return dx,dy,cA,cB,cC,cD,Filename

def getval(upperpath, image, folder): #gets pixel values for a single bmp image
    path = upperpath + folder
    Filename = os.path.join(path,image + '.bmp')
    value = misc.imread(Filename)
    Nx,Nz = value.shape
    return value, Nx, Nz, Filename

def getc(upperpath,namebase,folder): #applies getval to all four images in a set
    detectors = 'A', 'B', 'C', 'D'
    for det in detectors:
        image = namebase + '-' + det
        value, Nx, Ny, Filename = getval(upperpath,image,folder)
        if det == 'A':
            cA = value #raw bmp data
            
        if det == 'B':
            cB = value
            
        if det == 'C':
            cC = value
                
        if det == 'D':
            cD = value   
    #Getting pixel size
    fileref = open(os.path.join(upperpath+folder,image +'.txt'),"r")
    dummy = fileref.read().split()
    pixelsize = float(dummy[16][10:])/1000 #microns
    fileref.close()
    dx = dy = pixelsize
    
    return dx,dy,cA,cB,cC,cD,Filename

# Function for calculating mean normals
def getmeannormal(surf_dzgrid_dy, surf_dzgrid_dx, ixstart, ixstop, iystart, iystop):

    # Get normal vectors and do some averaging (meant for checking smooth surfaces)
    fy = surf_dzgrid_dy[:,:-1]
    fx = surf_dzgrid_dx[:-1,:]
    nynormal, nxnormal = np.shape(fx); #print nxnormal, nynormal
    normalgrid = np.zeros((nxnormal, nynormal,3))
    for ix in range(nxnormal):
        for iy in range(nynormal):
            normalgrid[ix,iy] = [-fx[iy,ix],-fy[iy,ix],1]

    mean = [0,0,0]
    count = 0
    for ix in range(ixstart, ixstop):
        for iy in range(iystart, iystop):
            mean += normalgrid[ix,iy]
            count += 1
    mean /= count
    mean /= sum(mean**2)**.5
    return mean

# Graphics functions
def myrectangle(draw,a,b,width=2):
    #fnt = ImageFont.truetype('Keyboard.ttf', 18)
    width = 2
    draw.line(((a[0],a[1]),(a[0],b[1]),(b[0],b[1]),(b[0],a[1]),(a[0],a[1])),width=width)

    
# Graphics functions
def myrectanglelabel(draw,a,b,label=''):
    fnt = ImageFont.truetype('Keyboard.ttf', 24)
    width = 4
    draw.line(((a[0],a[1]),(a[0],b[1]),(b[0],b[1]),(b[0],a[1]),(a[0],a[1])),width=width)
    if label!='':
       draw.text(a,' '+label,font=fnt)
       #draw.text(a,' '+label)
        
def linearFit(y, z):
    # Fitting with linearly generated sequence
    A = np.array([y, np.ones(y.size)])
    w = np.linalg.lstsq(A.T, z)[0]  # obtaining the parameters
    zline = w[0]*y+w[1]
    zfixed = z-zline  # substracting baseline from every point
    return zfixed

# Function for getting dot products
def mygets(nvec,theta):
    # This scales the normal vector so that nz=1 
    nxigrid = nvec[0]/nvec[2]
    nyigrid = nvec[1]/nvec[2]
    sA = (-nxigrid*np.sin(theta)+np.cos(theta)-1)/np.sqrt((1+nxigrid**2+nyigrid**2))
    sB = -(-nyigrid*np.sin(theta)+np.cos(theta)-1)/np.sqrt((1+nxigrid**2+nyigrid**2))
    sC = (+nxigrid*np.sin(theta)+np.cos(theta)-1)/np.sqrt((1+nxigrid**2+nyigrid**2))
    sD = -(+nyigrid*np.sin(theta)+np.cos(theta)-1)/np.sqrt((1+nxigrid**2+nyigrid**2))
    return sA, sB, sC, sD

def myrotation_matrix(axis, theta_deg):
    """
    Return the rotation matrix associated with clockwise rotation of an object
    (clockwise as seen by looking along the rotation axis toward the origin)
    about the given axis by theta degrees
    """
    import math
    theta = theta_deg*np.pi/180
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    #axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.asmatrix(np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]]))

import itertools
def polyfit2d(x, y, z, order=3, linear=False):
    """Two-dimensional polynomial fit. Based uppon code provided by 
    Joe Kington.

    References:
        http://stackoverflow.com/questions/7997152/
            python-3d-polynomial-surface-fit-order-dependent/7997925#7997925

    """
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
        if linear & (i != 0.) & (j != 0.):
            G[:, k] = 0
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    """Values to two-dimensional polynomial fit. Based upon code 
        provided by Joe Kington.
    """
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


def flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid, Rotx):
    # Rotates every point in the dataset
    thisshape = np.shape(surf_xseggrid)
    surf_xseggridp = np.zeros(thisshape)
    surf_yseggridp = np.zeros(thisshape)
    surf_zseggridp = np.zeros(thisshape)
    for ix in range (thisshape[1]):
        for iy in range (thisshape[0]):
            vec = np.matrix([surf_xseggrid[iy,ix],surf_yseggrid[iy,ix],surf_zseggrid[iy,ix]]).T
            vecp = Rotx*vec
            surf_xseggridp[iy,ix] = vecp[0]
            surf_yseggridp[iy,ix] = vecp[1]
            surf_zseggridp[iy,ix] = vecp[2]
    return surf_xseggridp, surf_yseggridp, surf_zseggridp

def extractflat(npzfile,dx,dy):
    
    # Extracting each segment, flattened
    nx1list = npzfile['nx1list']
    nx2list = npzfile['nx2list']
    ny1list = npzfile['ny1list']
    ny2list = npzfile['ny2list']
    solution = npzfile['solution']
    nsegments = len(nx1list)
    xseggridtot = []
    yseggridtot = []
    zseggridtot = []
    surf_xseggridtot = []
    surf_yseggridtot = []
    surf_zseggridtot = []
    thetazpptot = []
    
    for isegment in range(0,nsegments):

        # Extract this segment
        nx1=nx1list[isegment]; nx2=nx2list[isegment]; nxsegment = nx2-nx1+1
        ny1=ny1list[isegment]; ny2=ny2list[isegment]; nysegment = ny2-ny1+1
        surf_xseg = np.linspace(0,(nxsegment-1)*dx,nxsegment); 
        surf_yseg = np.linspace(0,(nysegment-1)*dy,nysegment); 
        surf_xseggrid, surf_yseggrid = np.meshgrid(surf_xseg,surf_yseg) # 1st index is y, 2nd is x
        surf_zseggrid = copy.copy(np.flipud(solution[ny1:ny2+1,nx1:nx2+1])) # This flips the y-coordinate

        # Fit a plane to the data and adjust data to start at the origin
        m = polyfit2d(\
                      surf_xseggrid.reshape(nysegment*nxsegment), \
                      surf_yseggrid.reshape(nysegment*nxsegment), \
                      surf_zseggrid.reshape(nysegment*nxsegment), \
                      linear=True,order=1)

        # Get the angles of the plane
        dzdy = m[1]; thetay = np.arctan(dzdy)*180/np.pi; #print 'y:', thetay

        # Get rotation matrix & flatten in one direction
        Roty = myrotation_matrix([1,0,0], -thetay)
        surf_xseggridp, surf_yseggridp, surf_zseggridp = \
            flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid, Roty)

        # Fit a plane to the data and adjust data to start at the origin
        mp = polyfit2d(\
                      surf_xseggridp.reshape(nysegment*nxsegment), \
                      surf_yseggridp.reshape(nysegment*nxsegment), \
                      surf_zseggridp.reshape(nysegment*nxsegment), \
                      linear=True,order=1)

        # Get the angle of the plane in another direction
        dzdx = mp[2]; thetaxp = np.arctan(dzdx)*180/np.pi; #print 'x:', thetaxp

        # Get rotation matrix & flatten in another direction
        Rotxp = myrotation_matrix([0,1,0], thetaxp)
        surf_xseggridpp, surf_yseggridpp, surf_zseggridpp = \
            flatten(surf_xseggridp, surf_yseggridp, surf_zseggridp, Rotxp)

        # Trying out the polyval2d
        surf_zseggrid_theory_long = polyval2d(\
                      surf_xseggrid.reshape(nysegment*nxsegment), \
                      surf_yseggrid.reshape(nysegment*nxsegment), \
                      m)
        surf_zseggrid_theory = surf_zseggrid_theory_long.reshape(nysegment,nxsegment)
        surf_xseggridp_theory, surf_yseggridp_theory, surf_zseggridp_theory = \
            flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid_theory, Roty)
        surf_xseggridpp_theory, surf_yseggridpp_theory, surf_zseggridpp_theory = \
            flatten(surf_xseggridp_theory, surf_yseggridp_theory, surf_zseggridp_theory, Rotxp)

        # Now rotate
        deltay = surf_yseggridpp_theory[0,-1]-surf_yseggridpp_theory[0,0]
        deltax = surf_xseggridpp_theory[0,-1]-surf_xseggridpp_theory[0,0]
        thetazpp = -np.arctan(deltay/deltax)*180/np.pi;
        Rotzpp = myrotation_matrix([0,0,1], thetazpp)
        surf_xseggridppp, surf_yseggridppp, surf_zseggridppp = \
            flatten(surf_xseggridpp, surf_yseggridpp, surf_zseggridpp, Rotzpp)
        surf_xseggridppp_theory, surf_yseggridppp_theory, surf_zseggridppp_theory = \
            flatten(surf_xseggridpp_theory, surf_yseggridpp_theory, surf_zseggridpp_theory, Rotzpp)

        # Now we have to extract an orthogonal subset
        dxsub = dysub = dx
        xsubstart = np.max(surf_xseggridppp_theory[[0,-1],0])+dxsub*2
        xsubstop = np.min(surf_xseggridppp_theory[[0,-1],-1])-dxsub*2
        ysubstart = np.max(surf_yseggridppp_theory[0,[0,-1]])+dysub*2
        ysubstop = np.min(surf_yseggridppp_theory[-1,[0,-1]])-dysub*2
        xsub = np.arange(xsubstart,xsubstop,dxsub)
        ysub = np.arange(ysubstart,ysubstop,dysub)
        sub_xseggrid, sub_yseggrid = np.meshgrid(xsub,ysub) # 1st index is y, 2nd is x
        nsuby, nsubx = np.shape(sub_xseggrid)
        surf_xseggridppp_theory_long = np.reshape(surf_xseggridppp_theory,nysegment*nxsegment)
        surf_yseggridppp_theory_long = np.reshape(surf_yseggridppp_theory,nysegment*nxsegment)
        points = np.vstack((surf_xseggridppp_theory_long,surf_yseggridppp_theory_long)).T # rows are x,y pairs
        values = np.reshape(surf_zseggridppp,nysegment*nxsegment)
        sub_zseggrid_long = griddata(points, values, (sub_xseggrid, sub_yseggrid), method='cubic')
        sub_zseggrid = np.reshape(sub_zseggrid_long,(nsuby, nsubx))
        
        xseggridtot.append(sub_xseggrid)
        yseggridtot.append(sub_yseggrid)
        zseggridtot.append(sub_zseggrid)
        surf_xseggridtot.append(surf_xseggrid)
        surf_yseggridtot.append(surf_yseggrid)
        surf_zseggridtot.append(surf_zseggrid)
        
    return xseggridtot, yseggridtot, zseggridtot, surf_xseggridtot, surf_yseggridtot, surf_zseggridtot


