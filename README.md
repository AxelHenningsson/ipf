Inverse Pole Figures
-------------------------------

An inverse pole figure is a great tool for visualising crystal orientation data. The essential idea is to map each crystal orientation to a series rgb color values and planar points to illustrate the texture of the entire orientation-set. A simple example of 100o orientations illustrated in an ipf for trigonal symmetry below.

````python
    from scipy.spatial.transform import Rotation
    ipf = inverse_pole_figure(crystal_system='trigonal', view_axis= np.array([0,0,1]))
    ipf.show( [Rotation.random().as_matrix().T for _ in range(1000)] )
````

![image](https://github.com/AxelHenningsson/ipf/assets/31615210/c77359e8-3869-44c0-ae88-937f8e15c331)
