import numpy as np
import matplotlib.pyplot as plt
import xfab.symmetry
from shapely import Polygon, Point
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation

class inverse_pole_figure(object):

    """An inverse pole figure (IPF).

    The IPF computes stereographic projections of reciprocal points and plots these in a fundamental zone.

    Stereographic projection is here defined with z=[0,0,1] as pole and z=-1 as projection plane.

    Args:
        crystal_system (:obj:`str`): One of triclinic, monoclinic, orthorhombic, tetragonal, trigonal
            hexagonal or cubic. The symmetry of the crystal_system is used in the creation of the ipf.
        view_axis (:obj:`numpy array`): The real space axis to which the ipf is related, shape=(3,)
        fundamental_triangle (:obj:`numpy array`): Corners of the zone in the z=-1 plane inside which the
            any reciprocal point will be stereographically projected exactly once - when taking the symmetry of
            the crystal into account.

    Attributes:
        crystal_system (:obj:`str`): One of triclinic, monoclinic, orthorhombic, tetragonal, trigonal
            hexagonal or cubic. The symmetry of the crystal_system is used in the creation of the ipf.
        view_axis (:obj:`numpy array`): The real space axis to which the ipf is related, shape=(3,)
        fundamental_triangle (:obj:`numpy array`): Corners of the zone in the z=-1 plane inside which the
            any reciprocal point will be stereographically projected exactly once - when taking the symmetry of
            the crystal into account.

    """

    def __init__(self, crystal_system, view_axis):
        crystal_dict = { 'triclinic': 1,
                        'monoclinic': 2,
                        'orthorhombic': 3,
                        'tetragonal': 4,
                        'trigonal': 5,
                        'hexagonal': 6,
                        'cubic': 7    }

        fundamental_dict = {
            'triclinic': [(-1,0,0), (1, 1, 0), (1, -1, 0)],
            'monoclinic': [(-1,0,0), (0, 1, 0), (1, 0, 0)],
            'orthorhombic': [(0, 0, 1), (0, 1, 0), (1, 0, 0)],
            'tetragonal': [(0, 0, 1), (0, 1, 0), (1, 1, 0)],
            'trigonal': [(0,0,1), (0,1,0), (1,0,0)],
            'hexagonal': [(0,0,1), (1,0,0), (2,-1,0)],
            'cubic':[(0,0,1), (0,1,1), (1,1,1)]
            }

        self.crystal_system = crystal_system
        self._P = xfab.symmetry.permutations(crystal_dict[crystal_system])
        self.view_axis = view_axis / np.linalg.norm(view_axis)
        self.fundamental_triangle = np.array(fundamental_dict[crystal_system]).T

        self._p = Polygon( self.get_fundamental_triangle(resolution=150).T )

    def colorbar(self):
        """Make and plot colorbar in the fundamental zone.
        """
        ori = [ Rotation.random().as_matrix().T for _ in range(2000) ]
        points = self.get_reciprocal_points( ori, filter_by_fundamental_zone=False )
        xy = self.stereographic_projection( points )
        edges = self.get_fundamental_triangle()
        rgb = self.get_colors(points)
        xmin,xmax = np.min(edges[0,:]), np.max(edges[0,:])
        ymin,ymax = np.min(edges[1,:]), np.max(edges[1,:])
        grid_x, grid_y = np.meshgrid(np.linspace(xmin, xmax, 128), np.linspace(ymin, ymax, 128), indexing='ij')

        mask = np.zeros( grid_x.shape,  dtype=bool)
        for i in range(grid_x.shape[0]):
            for j in range(grid_y.shape[1]):
                mask[i,j] = Point( [grid_x[i,j], grid_y[i,j]] ).within( self._p )

        grid_z1 = griddata(xy.T, rgb.T, (grid_x, grid_y), method='linear')
        grid_z1[~mask] = np.nan

        fig,ax = plt.subplots(1,1, figsize=(9,9))
        ax.axis('equal')
        ax.plot(edges[0], edges[1], 'gray', linewidth=int(10*64./grid_y.shape[0]))
        ax.axis('off')
        plt.pcolormesh(grid_x, grid_y, grid_z1)
        ft = self.fundamental_triangle
        xyf = self.stereographic_projection(self.fundamental_triangle)
        ax.annotate( str(ft[0,0])+' '+str(ft[1,0])+' '+str(ft[2,0]), (xyf[0,0]+0.01,xyf[1,0]-0.01), fontsize=20, annotation_clip=False )
        ax.annotate( str(ft[0,1])+' '+str(ft[1,1])+' '+str(ft[2,1]), (xyf[0,1]-0.01,xyf[1,1]+0.01), fontsize=20, annotation_clip=False)
        ax.annotate( str(ft[0,2])+' '+str(ft[1,2])+' '+str(ft[2,2]), (xyf[0,2]+0.01,xyf[1,2]+0.01), fontsize=20 , annotation_clip=False)
        plt.show()


    def show( self, ubi_matrices, point_size=50, alpha=1 ):
        """Plot the IPF given a set of ubi matrices.

        Args:
            ubi_matrices (:obj:`list` of `numpy array`): 3x3 matrices mapping from reciprocal to real space.
                The must map a reciprocal vector v_hkl to the same coordinate system as the view_axis.
            point_size (int, optional): scatter point point-size in plot. Defaults to 50.
            alpha (int, optional): scatter point transparency in plot. Defaults to 1.
        """
        points = self.get_reciprocal_points( ubi_matrices )
        xy = self.stereographic_projection( points )
        edges = self.get_fundamental_triangle()
        rgb = self.get_colors(points)

        fig,ax = plt.subplots(1,1, figsize=(8,8))
        ax.axis('equal')
        ax.scatter(xy[0], xy[1], c=rgb.T, s=point_size, alpha=alpha, linewidth=0)
        ax.plot(edges[0], edges[1], 'gray')
        ax.axis('off')
        ft = self.fundamental_triangle
        xyf = self.stereographic_projection(self.fundamental_triangle)
        ax.annotate( str(ft[0,0])+' '+str(ft[1,0])+' '+str(ft[2,0]), (xyf[0,0]+0.01,xyf[1,0]-0.01), fontsize=20 , annotation_clip=False)
        ax.annotate( str(ft[0,1])+' '+str(ft[1,1])+' '+str(ft[2,1]), (xyf[0,1]-0.01,xyf[1,1]+0.01), fontsize=20, annotation_clip=False)
        ax.annotate( str(ft[0,2])+' '+str(ft[1,2])+' '+str(ft[2,2]), (xyf[0,2]+0.01,xyf[1,2]+0.01), fontsize=20 , annotation_clip=False)
        va = self.view_axis
        plt.show()

    def get_fundamental_triangle(self, resolution = 100 ):
        """Generate points on the boundary of the fundamental zone

        Args:
            resolution (int, optional): Number of points to generate on fundamental zone boundary. Defaults to 100.

        Returns:
            :obj:`numpy array`: 2xN array of boundary points in the z=-1 plane.
        """
        t = np.linspace(0, 1, resolution)
        U = self.fundamental_triangle[:,0:1]
        V = self.fundamental_triangle[:,1:2]
        W = self.fundamental_triangle[:,2:3]
        U = U / np.linalg.norm(U)
        V = V / np.linalg.norm(V)
        W = W / np.linalg.norm(W)
        l1 = U*(1 - t) + V*t
        l2 = V*(1 - t) + W*t
        l3 = W*(1 - t) + U*t
        points = np.concatenate( (l1, l2, l3), axis=1 )
        q = self.stereographic_projection(l3)
        return self.stereographic_projection(points)

    def stereographic_projection(self, points):
        """Project 3D points unto the z=-1 plane using z=[0,0,1] as pole.

        Args:
            points (:obj:`numpy array`): 3,N array of 3D points.

        Returns:
            :obj:`numpy array`: 2xN array of boundary points in the z=-1 plane.
        """
        plane_z_offset = 0
        pole = np.array([0, 0, -1]).reshape(3, 1)
        norm = np.linalg.norm(points, axis=0)
        norm[norm==0]=1
        pn = points /norm
        s = (plane_z_offset - pole[2]) / (pole[2] - pn[2,:])
        projection = (pole + s*(pole - pn))
        return projection

    def get_colors(self, points):
        """Compute the IPF color of a set of 2D points.

        The color is selected as the normalised solution to

            p = r*t1 + g*t2 + c*t3

        where p is the input point, r,g,b the un-normalised color value,
        and t1,t2,t3 are the corners of the fundamental triangle/zone.

        Args:
            points (:obj:`numpy array`): 3xN array of 3D points.

        Returns:
            :obj:`numpy array`: 3xN array of RGB-color values.
        """
        uvw = np.linalg.lstsq( self.fundamental_triangle, points, rcond=-1 )[0]
        return ( uvw - np.min(uvw, axis=0)  ) / ( np.max(uvw- np.min(uvw, axis=0), axis=0))

    def _filter_fundamental_zone(self, points, number_of_orientations):
        """Remove projected points not in the fundamental zone.
        """

        xy = self.stereographic_projection(points)

        edges = self.get_fundamental_triangle( resolution=60 )
        xmin,xmax = np.min(edges[0,:]), np.max(edges[0,:])
        ymin,ymax = np.min(edges[1,:]), np.max(edges[1,:])
        mask = ~( (xy[0,:] < xmin) + (xy[1,:] < ymin) + (xy[0,:] > xmax) +  (xy[1,:] > ymax) )

        for i in range(len(mask)):
            if mask[i]:
                mask[i] = Point( xy[:,i] ).within( self._p )

        if np.sum(mask)!=number_of_orientations:
            self._p = Polygon( self.get_fundamental_triangle(resolution=900).T )
            for i in range(len(mask)):
                if mask[i]:
                    mask[i] = Point( xy[:,i] ).within( self._p )
        if np.sum(mask)!=number_of_orientations:
            print('WARNING: Some points on the edge of the fundamental zone could not be resolved... ')

        return points[:, mask]

    def get_reciprocal_points(self, ubi_matrices, filter_by_fundamental_zone=True):
        """_summary_

        Args:
            ubi_matrices (:obj:`list` of `numpy array`): 3x3 matrices mapping from reciprocal to real space.
                The must map a reciprocal vector v_hkl to the same coordinate system as the view_axis.
            filter_by_fundamental_zone (bool, optional): If False points outside the fundamental zone are preserved. 
                Defaults to True.

        Returns:
            :obj:`numpy array`: 3xN array of reciprocal points. These correspond to the points gerenated by applying
            the ubi matrices to the view_axis. I.e they are the reciprocal system view_axes.
        """
        points = []
        for ubi in ubi_matrices:
            hkl = ubi @ self.view_axis
            points.extend( list(self._P @ hkl) )
        points = np.array(points).T
        points = np.concatenate( (points, -points), axis=1 )
        if filter_by_fundamental_zone:
            return self._filter_fundamental_zone(points, len(ubi_matrices))
        else:
            return points


if __name__ == "__main__":
    # E X A M P L E:
    #plt.style.use('dark_background')
    #ipf = inverse_pole_figure(crystal_system='trigonal', view_axis= np.array([0,0,1]))
    #ipf.show( [Rotation.random().as_matrix().T for _ in range(1000)] )
    #ipf.colorbar()