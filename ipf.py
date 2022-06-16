from xfab import tools, symmetry
import numpy as np
import os
from ImageD11.grain import read_grain_file
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


class RodriguezRotator(object):
    """Object for rotating vectors in the plane described by yhe unit normal rotation_axis.

    Args:
        rotation_axis (:obj:`numpy array`): A unit vector in 3d euclidean space (``shape=(3,)``)

    Attributes:
        rotation_axis (:obj:`numpy array`): A unit vector in 3d euclidean space (``shape=(3,)``)
        K (:obj:`numpy array`): (``shape=(3,3)``)
        K2 (:obj:`numpy array`): (``shape=(3,3)``)
        I (:obj:`numpy array`): (``shape=(3,3)``)

    """

    def __init__(self, rotation_axis):
        assert np.allclose(np.linalg.norm(rotation_axis),
                           1), "The rotation axis must be length unity."
        self.rotation_axis = rotation_axis
        rx, ry, rz = self.rotation_axis
        self.K = np.array([[0, -rz, ry],
                           [rz, 0, -rx],
                           [-ry, rx, 0]])
        self.K2 = self.K.dot(self.K)

    def get_rotation_matrix(self, rotation_angle):
        return np.eye(3, 3) + np.sin(rotation_angle) * self.K + \
            (1 - np.cos(rotation_angle)) * self.K2

    def __call__(self, vectors, rotation_angle):
        """Rotate a vector in the plane described by v1 and v2 towards v2 a fraction s=[0,1].

        Args:
            vectors (:obj:`numpy array`): A set of vectors in 3d euclidean space to be rotated (``shape=(3,N)``)
            rotation_angle (:obj:`float`): Radians to rotate vectors around the rotation_axis (positive rotation).

        Returns:
            Rotated vectors (:obj:`numpy array`) of ``shape=(3,N)``.

        """
        R = self.get_rotation_matrix(rotation_angle)
        return R.dot(vectors)

def _backtrace(coord, south_pole):
    n = coord - south_pole
    n = n / np.linalg.norm(n)
    #print(n)
    # k = coord + s * n
    # k*k == 1
    # coord*coord + 2 * s * n * coord + s*s - 1 == 0
    q =  coord.dot(coord) - 1
    p =  -2 * n.dot(coord)
    s1 = (p/2.) + np.sqrt( (p*p/4.) - q )
    s2 = (p/2.) - np.sqrt( (p*p/4.) - q )
    vector = coord + np.max([s1,s2])*n
    #print('norm', np.linalg.norm(vector))
    return vector

def plot_cmap( triangle_corners ):
    south_pole = np.array([0,0,-1])

    ti = np.linalg.inv(triangle_corners.T)

    xb, yb = ipf_bounds(triangle_corners, number_of_points=100)
    xlow  = np.min(xb) - 0.02
    xhigh = np.max(xb) + 0.02
    ylow  = np.min(yb) - 0.02
    yhigh = np.max(yb) + 0.02

    normals, points = [],[]
    c1,c2,c3 = triangle_corners
    for cc in [(c1,c2),(c2,c3),(c3,c1)]:
        normals.append( np.cross(cc[0],cc[1]) )
        points.append( cc[0] )

    resolution = 250
    cmap = np.ones((resolution, resolution ,3))
    for i,y in enumerate(np.linspace(yhigh, ylow, resolution)):
        for j,x in enumerate(np.linspace(xlow, xhigh, resolution)):
            coord = np.array([x,y,0])
            vector = _backtrace(coord, south_pole)

            in_symmetry_triangle = True
            for nn,pp in zip(normals, points):
                if (vector - pp).dot(nn) < -1e-8:
                    in_symmetry_triangle = False
                    break
            
            if in_symmetry_triangle:
                c = np.abs( ti.dot(vector) )
                cmap[i,j,:] = c / np.max(c)

    fig,ax = plt.subplots(1,1)
    ax.imshow(cmap)
    xt, yt = ipf_bounds(triangle_corners, number_of_points=25)
    xt = (xt -xlow) / (xhigh - xlow)
    xt = (resolution * xt)
    yt = (-(yt - ylow) / (yhigh - ylow) ) + 1.0 - (0.5/resolution)
    yt = (resolution * yt)
    ax.plot(xt, yt, 'k-', linewidth=2.5)
    xc, yc = ipf_corners(triangle_corners)
    xc = (xc -xlow) / (xhigh - xlow)
    xc = (resolution * xc)
    yc = (-(yc -ylow) / (yhigh - ylow) ) + 1.0 - (0.5/resolution)
    yc = (resolution * yc)
    ax.text( xc[0]-0.05*resolution, yc[0]+0.1*resolution, str(triangle_corners[0,:]).replace(' ', '').replace(']', '').replace('[', ''), fontsize=24 )
    ax.text( xc[1]-0.03*resolution, yc[1]+0.1*resolution, str(triangle_corners[1,:]).replace(' ', '').replace(']', '').replace('[', ''), fontsize=24 )
    ax.text( xc[2]-0.033*resolution, yc[2]-0.035*resolution, str(triangle_corners[2,:]).replace(' ', '').replace(']', '').replace('[', ''), fontsize=24 )
    ax.plot(xc, yc, 'ko', markersize=10)
    ax.axis('off')

    return cmap, fig, ax

def ipf_corners(triangle_corners):
    """Return the symmetry triangle corners of an inverse pole figure.
    """
    X, Y = [], []
    for c in triangle_corners:
        x,y = stereographic_projection(c)
        X.append(x)
        Y.append(y)
    return X, Y

def ipf_bounds(triangle_corners, number_of_points):
    """Return gridded points on the inverse pole figure symmetry triangle.
    """
    norm_triangle_corners = triangle_corners / np.linalg.norm( triangle_corners, axis=1 ).reshape(3,1)
    X, Y = [], []
    for i in range(norm_triangle_corners.shape[0]):
        c1 = norm_triangle_corners[i,:]
        if i==2:
            c2 = norm_triangle_corners[0,:]
        else:
            c2 = norm_triangle_corners[i+1,:]
        
        axis = np.cross(c1, c2)
        axis = axis / np.linalg.norm(axis)
        angle = np.arccos( c1.dot(c2) )
        rotator = RodriguezRotator(axis)
        for ang in np.linspace(0, angle, number_of_points):
            rax = rotator(c1, ang)
            x,y = stereographic_projection(rax)
            X.append(x)
            Y.append(y)
    return X, Y

def stereographic_projection(vector):
    """Steriographically project a vector unto the unit sphere equator plane (z=0).
    ref: https://en.wikipedia.org/wiki/Stereographic_projection

        Args:
            vector (:obj:`numpy array`): shape=(3,)
        
        Returns:
            x,y (:obj:`tuple` of `float`), projected coordinates in equator plane.
    """
    normalised_vector = vector / np.linalg.norm( vector )

    if np.abs(normalised_vector[2])<=1e-8:
        intersection = normalised_vector
    else:
        south_pole = np.array([0,0,-1])
        n = normalised_vector - south_pole
        s = -south_pole[2] / n[2]
        intersection = south_pole + s*n

    assert np.abs(intersection[2])<1e-8, str(intersection)

    return intersection[0], intersection[1]

def cartesian_to_polar(x, y):
    th = np.arctan2(y, x)
    r = np.sqrt(x*x + y*y)
    return th, r

def in_symmetry_triangle(vector, triangle_corners):
    c1,c2,c3 = triangle_corners
    for cc in [(c1,c2),(c2,c3),(c3,c1)]:
        normal = np.cross(cc[0],cc[1])
        point = cc[0]
        if (vector - point).dot(normal) < -1e-8:
            return False
    return True


def inverse_pole_figure(ubi_matrices, crystal_system, triangle_corners):
    """
    """
    sample_axis = np.array([0,0,1])
    X, Y, colors = [], [], []
    ti = np.linalg.inv(triangle_corners.T)
    symmetry_rotations = symmetry.rotations(crystal_system)

    for ubi in ubi_matrices:
        
        u, b = tools.ub_to_u_b( np.linalg.inv(ubi) )
        #u, b = tools.ubi_to_u_b(ubi)

        for symrot in symmetry_rotations:
            
            #u_new = symrot.dot(u)

            u_new = u.dot(symrot)

            #ubi_new = symrot.dot( ubi )

            ubi_new = np.linalg.inv(u_new.dot(b))
            sample_axis_in_crystal = ubi_new.dot( sample_axis )
            sample_axis_in_crystal = sample_axis_in_crystal/ np.linalg.norm(sample_axis_in_crystal)

            if sample_axis_in_crystal[2]<0:
                continue

            if not in_symmetry_triangle(sample_axis_in_crystal, triangle_corners):
                continue

            x, y = stereographic_projection(sample_axis_in_crystal)
            # th, r = cartesian_to_polar(x, y)
        
            X.append( x )
            Y.append( y )

            c = np.abs( ti.dot(sample_axis_in_crystal) )
            colors.append( c / np.max(c) )

    return np.array(X), np.array(Y), colors

if __name__=='__main__':

    # ubi_matrices = []
    # base_path = os.path.dirname(__file__)
    # for layer in range(15,16):
    #     gr = os.path.join(base_path, 'ubi_files', 'matrices_'+str(layer).zfill(3)+'.ubi')
    #     ubi_matrices.extend( [ g.ubi for g in read_grain_file(gr)] )

    #print(len(ubi_matrices))

    B = tools.form_b_mat([1., 1., 1., 90., 90., 90.])

    # u = np.array([[1,0,0],[0,1,0],[0,0,1]])
    # ubi_matrices = [ np.linalg.inv(u.dot(B)) ]

    ubi_matrices=[]
    zhat = np.array([0,0,1])
    U = Rotation.random(100).as_matrix()
    UBIs = np.linalg.inv( U.dot(B) )
    for ubi in UBIs: ubi_matrices.append(ubi)
    # for ubi in UBIs:
    #     hkl = ubi.dot(zhat)
    #     n = hkl / np.linalg.norm(hkl)
    #     hkl_dir = np.array([1,0,0])
    #     hkl_dir =  hkl_dir / np.linalg.norm(hkl_dir)
    #     if np.abs(n.dot(hkl_dir)-1)<0.01: # z-axis close to 111
    #         ubi_matrices.append(ubi)

    #u = np.array([[-1,0,0],[0,-1,0],[0,0,1]])
    #ubi_matrices = [ np.linalg.inv(B) ]
 
    sample_axis = np.array([0,0,1])
    triangle_corners = np.array([[0,0,1], [1,1,1], [0,1,0] ])
    x, y, colors = inverse_pole_figure(ubi_matrices, crystal_system=7, triangle_corners=triangle_corners)
    print(len(colors), len(colors)//24)
    #cmap, fig_cmap, ax_cmap = plot_cmap( triangle_corners )

    fig = plt.figure()
    ax = fig.add_subplot()
    c = ax.scatter(x, y, c=colors, s=35, alpha=0.5)

    xt, yt = ipf_bounds(triangle_corners, number_of_points=25)
    c = ax.plot(xt, yt, 'k-')

    xc, yc = ipf_corners(triangle_corners)

    ax.text( xc[0]-0.02, yc[0]-0.04, str(triangle_corners[0,:]).replace(' ', '').replace(']', '').replace('[', ''), fontsize=24 )
    ax.text( xc[1]-0.01, yc[1]-0.04, str(triangle_corners[1,:]).replace(' ', '').replace(']', '').replace('[', ''), fontsize=24 )
    ax.text( xc[2]-0.015, yc[2]+0.015, str(triangle_corners[2,:]).replace(' ', '').replace(']', '').replace('[', ''), fontsize=24 )

    #c = ax.plot(xc, yc, 'ko', markersize=12)
    ax.axis('off')

    plt.show()