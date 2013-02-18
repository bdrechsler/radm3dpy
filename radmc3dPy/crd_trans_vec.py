import numpy as np
import scipy.linalg as la

# ----------------------------------------------------------------------------------------------------------------------------
def meshgrid3d_vec(a=None, b=None, c=None):
    """
    Function to create a 3D mesh from three vectors (1d ndarrays)
    This function works in a brute force way (i.e. using for loops) 
    with no optimization and it might be slow for long vectors

    INPUT:
    ------
        a,b,c : 1d ndarray vectors

    OUTPUT:
    -------
        The function returns an ndarray with [Na, Nb, Nc, 3] dimensions

    """

    if a=None:
        print 'ERROR' 
        print 'Number of input vectors should be three'
        return -1
    if b=None:
        print 'ERROR' 
        print 'Number of input vectors should be three'
        return -1
    if c=None:
        print 'ERROR' 
        print 'Number of input vectors should be three'
        return -1

    na = a.shape[0]
    nb = b.shape[0]
    nc = c.shape[0]
    res = np.zeros([na,nb,nc,3], dtype=float64)
    for ia in range(na):
        for ib in range(nb):
            for ic in range(nc):
                res[ia,ib,ic,:] = np.array([a[ia], b[ib], c[ic]])

    return res

# ----------------------------------------------------------------------------------------------------------------------------
def ctrans_sph2cart(crd=[[0.], [0,], [0.]], reverse=False):
    """
    Function to transform coordinates between spherical to cartesian systems

    INPUT : 
    -------
            crd      : Three element array containing the input
                       coordinates [x,y,z] or [r,phi,theta] by default
                       the coordinates assumed to be in the cartesian system
    KEYWORDS :
    ----------
            reverse=False : Calculates the inverse trasnformation
                       (cartesian -> spherical). In this case crd should be [r,phi,theta]

    OUTPUT : 
    --------
            result   : A three element array containg the output
                       coordinates [r,phi,theta] or [x,y,z]
    """
    if (reverse==False):
        # First create the coordinate array
        vcrd = meshgrid3d_vec(crd[0], crd[1], crd[2])
        # Flatten the arrays and separate the different coordinate and vector components
        r = vcrd[:,:,:,0].flatten('F')
        phi = vcrd[:,:,:,1].flatten('F')
        # I do not remember why I have to add 1e-50 to theta put I probably had a reason.. :-)
        theta = vcrd[:,:,:,2].flatten('F') + 1e-50

        x = np.sin(theta) * np.cos(phi) * r
        y = np.sin(theta) * np.sin(phi) * r
        z = np.cos(theta) * r
    
        nx = x.shape[0]
        ny = y.shape[0]
        nz = z.shape[0]


        vcrd_new = array(vcrd) * 0.
        vcrd_new[:,:,:,0] = x.reshape(nx,ny,nz, order='F')
        vcrd_new[:,:,:,1] = y.reshape(nx,ny,nz, order='F')
        vcrd_new[:,:,:,2] = z.reshape(nx,ny,nz, order='F')

    else:

        # First create the coordinate array
        vcrd = meshgrid3d_vec(crd[0], crd[1], crd[2])
        # Flatten the arrays and separate the different coordinate and vector components
        x = vcrd[:,:,:,0].flatten('F')
        y = vcrd[:,:,:,1].flatten('F')
        z = vcrd[:,:,:,2].flatten('F') 
        
        r     = np.sqrt(x**2 + y**2 + z**2)
        phi   = np.arccos(x / np.sqrt(x**2 + y**2) + 1e-90)
        theta = np.arccos(z / r)
        
        nr = r.shape[0]
        nphi = phi.shape[0]
        ntheta = theta.shape[0]

        phi[y<0.] = 2.0 *np.pi - phi[y<0.]
       
        vcrd_new = array(vcrd) * 0.
        vcrd_new[:,:,:,0] = r.reshape(nx,ny,nz, order='F')
        vcrd_new[:,:,:,1] = phi.reshape(nx,ny,nz, order='F')
        vcrd_new[:,:,:,2] = theta.reshape(nx,ny,nz, order='F')
        

    return vcrd_new
# ----------------------------------------------------------------------------------------------------------------------------
def vtrans_sph2cart(crd=None, v=None, reverse=False):
    """
    Function to transform velocities between spherical to cartesian systems

    INPUT : 
    -------
            crd      : Three element array containing the input
                       coordinates [x,y,z] or [r,phi,theta] by default
                       the coordinates assumed to be in the cartesian system

            v        : Three element array containing the input
                       velocities in the same coordinate system as crd


    KEYWORDS :
    ----------
            reverse=False : Calculates the inverse trasnformation (cartesian -> spherical)

    OUTPUT : 
    --------
            result   : A three element array containg the output
                       velocities [vr,vphi,vtheta] or [vx,vy,vz]


    NOTE!!!!! The velocities in the spherical system are not angular velocities!!!!
    v[1] = dphi/dt * r
    v[2] = dtheta/dt * r
    """
    if (reverse==False):
        # First create the coordinate array
        vcrd = meshgrid3d_vec(crd[0], crd[1], crd[2])
        # Flatten the arrays and separate the different coordinate and vector components
        r = vcrd[:,:,:,0].flatten('F')
        phi = vcrd[:,:,:,1].flatten('F')
        theta = vcrd[:,:,:,2].flatten('F')
        nr = r.shape[0]
        nphi = phi.shape[0]
        ntheta = theta.shape[0]
        
        vr     = v[:,:,:,0]
        vphi   = v[:,:,:,1]
        vtheta = v[:,:,:,2]
        
        vx     = vr*np.sin(theta)*np.cos(phi) - vphi*np.sin(phi) + vtheta*np.cos(theta)*np.cos(phi)
        vy     = vr*np.sin(theta)*np.sin(phi) + vphi*np.cos(phi) + vtheta*np.cos(theta)*np.sin(phi)
        vz     = vr*np.cos(theta) - vtheta*np.sin(theta)

        v_new = array(v) * 0.
        v_new[:,:,:,0] = vx.reshape(nr,nphi,ntheta, order='F')
        v_new[:,:,:,1] = vy.reshape(nr,nphi,ntheta, order='F')
        v_new[:,:,:,2] = vz.reshape(nr,nphi,ntheta, order='F')

    else:
        
        crd_sph = ctrans_sph2cart(crd, reverse=True)
        r       = crd_sph[:,:,:,0].flatten('F')
        phi     = crd_sph[:,:,:,1].flatten('F')
        theta   = crd_sph[:,:,:,2].flatten('F')
   
        nr = r.shape[0]
        nphi = phi.shape[0]
        ntehta = theta.shape[0]

        a       = [[np.sin(theta)*np.cos(phi), -np.sin(phi), np.cos(theta)*np.cos(phi)],\
                   [np.sin(theta)*np.sin(phi), np.cos(phi), np.cos(theta)*np.sin(phi)],\
                   [np.cos(theta), 0., -np.sin(theta)]]
    
        a       = np.matrix(np.array(a, dtype=np.float64))
        ai      = a.I
        ai      = np.array(ai)

        vr     = vx * ai[0,0] + vy * ai[0,1] + vz * ai[0,2]
        vphi   = vx * ai[1,0] + vy * ai[2,1] + vz * ai[1,2]
        vtheta = vx * ai[2,0] + vy * ai[1,1] + vz * ai[2,2]

        v_new = array(v) * 0.
        v_new[:,:,:,0] = vr.reshape(nr,nphi,ntheta, order='F')
        v_new[:,:,:,1] = vphi.reshape(nr,nphi,ntheta, order='F')
        v_new[:,:,:,2] = vtheta.reshape(nr,nphi,ntheta, order='F')
        
    return v_new
# ----------------------------------------------------------------------------------------------------------------------------
def csrot(crd=None, ang=None, deg=False):
    """
    Procedure to make coordinate system rotation
 
    INPUT : 
    -------
           crd  : three element vector containing the coordinates of a
                  given point in a cartesian system

           ang  : three element array, angles of rotation around the x,y,z axes

    KEYWORDS : 
    ----------
           deg=True : if this keyword is set angles should be given in
                  angles instead of radians (as by default)
 
    Rotation matrices :
    -------------------
    X-axis

     |      1               0            0        |
     |      0          cos(alpha)     -sin(alpha)  | 
     |      0         sin(alpha)     cos(alpha)  |

    Y-axis

     |   cos(beta)          0        sin(beta)   |
     |      0               1            0        |
     |   -sin(beta)          0         cos(beta)   |

    Z-axis
 
     |   cos(gamma)     -sin(gamma)       0        |
     |  sin(gamma)     cos(gamma)       0        |
     |      0               0            1        |

    """

    # First create the coordinate array
    vcrd = meshgrid3d_vec(crd[0], crd[1], crd[2])
    # Flatten the arrays and separate the different coordinate and vector components
    x = vcrd[:,:,:,0].flatten('F')
    y = vcrd[:,:,:,1].flatten('F')
    z = vcrd[:,:,:,2].flatten('F')

   
    # Rotation angles
    xang = ang[0]
    yang = ang[1]
    zang = ang[2]

#
# Convert degree into radian if the angles are given in degree
#

    if (deg==True):
        xang = xang / 180.0 * np.pi
        yang = yang / 180.0 * np.pi
        zang = zang / 180.0 * np.pi


#
# Rotation around the x axis
#

    if (xang!=0.0):
        x = x
        y = np.cos(xang)*y  - np.sin(xang)*z
        z = np.sin(xang)*y + np.cos(xang)*z 
    
#
# Rotation around the y axis
#

    if (yang!=0.0):
        x = np.cos(yang)*x + np.sin(yang)*z
        y = y
        z = -np.sin(yang)*x + np.cos(yang)*z
        

#
# Rotation around the z axis
#
  
    if (zang!=0.0):
        x = np.cos(zang)*x - np.sin(zang)*y + 0.0
        y = np.sin(zang)*x + np.cos(zang)*y + 0.0
        z = z

    dumx = x.reshape(nx,ny,nz, order='F')
    dumy = y.reshape(nx,ny,nz, order='F')
    dumz = z.reshape(nx,ny,nz, order='F')

    vcrd_new = array(vcrd) * 0.
    vcrd_new[:,:,:,0] = x.reshape(nx,ny,nz, order='F')
    vcrd_new[:,:,:,1] = y.reshape(nx,ny,nz, order='F')
    vcrd_new[:,:,:,2] = z.reshape(nx,ny,nz, order='F')

    return vcrd_new
# ----------------------------------------------------------------------------------------------------------------------------
def vrot(crd=None, v=None, ang=None):
    """
    Procedure to rotate a vector in spherical coordinate system

    First transform the vector to cartesian coordinate system do the rotation then make the
     inverse transformation

    INPUT : 
    -------
           crd  : three element vector containing the coordinates of a
                  given point in a cartesian system

           v    : three element array, angles of rotation around the x,y,z axes

           ang  : angle around the x, y, z, axes with which the vector should be rotated

    """
# Convert the position vector to cartesian coordinate system
    crd_xyz   = ctrans_sph2cart(crd=crd)
# Convert the velocity vector to cartesian coordinate system
    v_xyz     = vtrans_sph2cart(crd=crd, v=v)
# Rotate the vecto
    v_xyz_rot = csrot(crd=v_xyz, ang=ang)
# Transform the rotated vector back to the spherical coordinate system
    v_rot     = vtrans_sph2cart(crd=crd_xyz, v=v_xyz_rot, reverse=True) 
    
    return v_rot
