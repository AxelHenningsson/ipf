A bare-bones implementation of - Inverse Pole Figures
-------------------------------------

An inverse pole figure (IPF) is a way to plot rotation data (elements in SO3, i.e., 3x3 rotation matrices, tuples of euler angles, etc.)
in 2D while taking care to consider the symmetry of an attached point-group. It is a great and widely used tool for visualisation of
crystallographic textures. It is also a wildly confusing concept wing to the fact that higher dimensional rotation elements are
compressed into 2D coordinates and rgb-color values. This document is an attempt of me to clarify the interpretation and
mathematical definition of the IPF. In this quest it seem use-full to implement a bare-bones version of the IPF which is what
the code in this repository aims at. 

*Discalimer: this is my interpretation of these concepts - if you find any misstakes please let me know! -*

Example
-------------------------------------
Lets create 100 random orientation matrices, just to have somethign to plot
````python
from scipy.spatial.transform import Rotation
random_orientations = [Rotation.random().as_matrix().T for _ in range(1000)]
````
Next, we create an inverse pole figure object in the trigonal crystal system. We specify the sample normal direction with the `view_axis` keyword. In this example we are looking at the texture with respect tot he sample **z**-axis.
````python
    from ipf import inverse_pole_figure
    ipf = inverse_pole_figure(crystal_system='trigonal', view_axis= np.array([0,0,1]))
````
We can now visualise the data we created in an IPF using the `show()` command. Since the texure is uniform, the IPF has data points scattered over the full area of the *"stereographic-fundamental zone"*.
````python
    ipf.show( random_orientations )
````
![image](https://github.com/AxelHenningsson/ipf/assets/31615210/c77359e8-3869-44c0-ae88-937f8e15c331)
We may also produce a colorbar that covers the *"stereographic-fundamental zone"* using the  `colorbar()` command
````python
    ipf.colorbar()
````
IN the next section the result of this `colorbar()` command on each of the 7 crystal systems is shown.

My interpretation of the stereographic-fundamental-zones
-------------------------------------
![image](https://github.com/AxelHenningsson/ipf/assets/31615210/75b22698-96cb-4256-863b-066a67ac1dc8)



Installation
-------------------------------------
You may install from source using `pip` as
````python
    git clone https://github.com/AxelHenningsson/ipf.git
    cd ipf
    pip install -e .
````
Also, this is a single module repo with a handfull of dependecies so you will probably figure it out 😙🤟
