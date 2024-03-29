**Inverse Pole Figures** - *A Bare-bones Educational Implementation*
-------------------------------------

An inverse pole figure (IPF) is a way to plot rotation data (elements in SO3, i.e., 3x3 rotation matrices, tuples of euler angles, etc.)
in 2D while taking care to consider the symmetry of an attached point-group. It is a great and widely used tool for visualisation of
crystallographic textures. It is also a wildly confusing concept wing to the fact that higher dimensional rotation elements are
compressed into 2D coordinates and rgb-color values. This document is an attempt of me to clarify the interpretation and
mathematical definition of the IPF. In this quest it seem use-full to implement a bare-bones version of the IPF which is what
the code in this repository aims at. 

*Discalimer: this is my interpretation of these concepts - if you find any misstakes please let me know! -*

**Note** *the color coding is currently non-optimal, future persepectives should go in this direction:*

[Area-preserving colour coding of inverse pole figure domain. KARTHIKEYAN, T. (2017) Journal of Microscopy, 267: 107-113.](https://doi.org/10.1111/jmi.12578)

Example
-------------------------------------
Lets create 100 random crystal cell matrices, just to have something to plot
````python
from scipy.spatial.transform import Rotation
random_ubis = [Rotation.random().as_matrix().T for _ in range(1000)]
````
Next, we create an inverse pole figure object in the trigonal crystal system. We specify the sample normal direction with the `view_axis` keyword. In this example we are looking at the texture with respect tot he sample **z**-axis.
````python
    from ipf.pole_figures import inverse_pole_figure
    crystal_system ='cubic'
    view_axis = np.array([0,0,1])
    ipf = inverse_pole_figure(crystal_system, view_axis)
````
We can now visualise the data we created in an IPF using the `show()` command. Since the texure is uniform, the IPF has data points scattered over the full area of the *"stereographic-fundamental zone"*.
````python
    fig, ax = ipf.show( random_ubis )
    plt.show()
````
![image](https://github.com/AxelHenningsson/ipf/assets/31615210/7d5e9365-f494-46b6-907a-1e5e6230cf46)
We may also produce a colorbar that covers the *"stereographic-fundamental zone"* using the  `colorbar()` command
````python
    fig, ax = ipf.colorbar()
    plt.show(
````
In the next section the result of this `colorbar()` command on each of the 7 crystal systems is shown. Note that the color coding is non-optimal in this implementation, future persepectives should go in this direction:

[Area-preserving colour coding of inverse pole figure domain. KARTHIKEYAN, T. (2017) Journal of Microscopy, 267: 107-113.](https://doi.org/10.1111/jmi.12578)


More demos can be found as notebooks in the demos/ folder.

My interpretation of the stereographic-fundamental-zones
-------------------------------------
![image](https://github.com/AxelHenningsson/ipf/assets/31615210/aeeb3580-023a-489f-ae42-837521aae91d)

Algorithm & Interpretation
-------------------------------------
To create the inverse pole figures there are five essental steps taken:

1. An  input ubi matrix mapps a fixed real space point to a reciprocal point. The ubi matrix is just that the 3x3 matrix that maps from the reciprocal space to the real. We call the fixed point the "view axis". The choice of view axis determines how the IPF should be interpreted.
2. The Reciprocal point is now normalised and mapped to a 2D plane at a fixed z=-1 position [using a stereographc projection](https://en.wikipedia.org/wiki/Stereographic_projection). (The pole is z=1 in this implementation.)
3. If the point is in the fundamental zone (this is the corners of a spherical triangle) defined by the crystal symmetry we keep the point around. (Due to the crystal symmetry the ubi matrix is not unique)
4. The 2D projected point is mapped to an RGB color value. In this project we use a simple, but non optimal map as: p = r x t1 + g x t2 + b x t3, where p is the input point, r,g,b the un-normalised color value and t1,t2,t3 are the spherical-corners of the fundamental triangle. we normalise the rgb value by subtracting the minimum and then dividing with the maximum. Better colormaps exists: [Area-preserving colour coding of inverse pole figure domain. KARTHIKEYAN, T. (2017) Journal of Microscopy, 267: 107-113.](https://doi.org/10.1111/jmi.12578). Not that when the r=1 and g=b=0 we have p=t1 and so red will mean that the ubi matrix is mapping the view axis to the first corner of the fundamental zone.
5. We may now plot the 2D points and their color along with the boundaries of the fundamental zone. That is an IPF. The position and color in the IPF is displaying which reciprocal point the the ubi matrix is mappining the view axis to. In crystalography this can be though of as a measure of set of Miller planes that have a normal that is aligned with the view axis.

Installation
-------------------------------------
You may install from source using `pip` as
````python
    git clone https://github.com/AxelHenningsson/ipf.git
    cd ipf
    pip install -e .
````
Also, this is a single module repo with a handfull of dependecies so you will probably figure it out 😙🤟
