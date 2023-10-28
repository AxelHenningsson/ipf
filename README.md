A bare-bones implementation of - Inverse Pole Figures
-------------------------------------

An inverse pole figure (IPF) is a way to plot rotation data (elements in SO3, i.e., 3x3 rotation matrices, tuples of euler angles, etc.)
in 2D while taking care to consider the symmetry of an attached point-group. It is a great and widely used tool for visualisation of
crystallographic textures. It is also a wildly confusing concept wing to the fact that higher dimensional rotation elements are
compressed into 2D coordinates and rgb-color values. This document is an attempt of me to clarify the interpretation and
mathematical definition of the IPF. In this quest it seem use-full to implement a bare-bones version of the IPF which is what
the code in this repository aims at. 

*Discalimer: this is my interpretation of these concepts - if you find any misstakes please let me know! -*

My interpretation of the stereographic-fundamental-zones
-------------------------------------
![image](https://github.com/AxelHenningsson/ipf/assets/31615210/75b22698-96cb-4256-863b-066a67ac1dc8)

````python
    from scipy.spatial.transform import Rotation
    ipf = inverse_pole_figure(crystal_system='trigonal', view_axis= np.array([0,0,1]))
    ipf.show( [Rotation.random().as_matrix().T for _ in range(1000)] )
````

![image](https://github.com/AxelHenningsson/ipf/assets/31615210/c77359e8-3869-44c0-ae88-937f8e15c331)
