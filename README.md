# Evanescent Plane Wave Approximation of Helmholtz Solutions in Spherical Domains

## Wave Goodbye to Instability!

<div><img src="/images/README_image.png" width="250px" align="right"></div>

The recent results presented within the paper [*Stable approximation of Helmholtz solutions by evanescent plane waves*](https://arxiv.org/abs/2202.05658) have led to significant developments in achieving stable Helmholtz solution approximations by plane wave superposition. The study shows that the numerical instability and ill-conditioning inherent in plane wave-based Trefftz methods can be effectively overcome with regularization techniques, provided there exist accurate approximations in the form of expansions with bounded coefficients. While propagative plane waves fail to yield stable approximations due to the exponential growth of the expansion coefficients, evanescent plane waves, which contain high Fourier modal content, provide both accurate and stable results. The developed numerical approach &mdash; available at [https://github.com/EmileParolin/evanescent-plane-wave-approx](https://github.com/EmileParolin/evanescent-plane-wave-approx) &mdash; results in substantial improvements when compared to conventional propagative plane wave schemes.

The numerical recipe proposed here has been tailored for the 3D setting and extended with new sampling strategies involving extremal system of points (see [https://web.maths.unsw.edu.au/~rsw/Sphere/](https://web.maths.unsw.edu.au/~rsw/Sphere/)). The results show the desired accuracy and bounded-coefficient stability, in line with the two-dimensional case.

## Play with the code!

The MATLAB code is available in [src](src). For this to work properly, it is necessary to download the extremal point sets and the related weights from [https://web.maths.unsw.edu.au/~rsw/Sphere/MaxDet/MaxDet1.html#ExtSys](https://web.maths.unsw.edu.au/~rsw/Sphere/MaxDet/MaxDet1.html#ExtSys). Moreover, the [example](example) section contains a MATLAB live script that enables running some code and playing the parameters, as well as a static version of the same script. This provides a means of navigating through the construction of both propagative and evanescent plane wave approximation sets.
