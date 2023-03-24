# Stable Helmholtz Solution Approximation by Evanescent Plane Waves in Spherical Domains

The recent results presented within the paper [*Stable approximation of Helmholtz solutions by evanescent plane waves*](https://arxiv.org/abs/2202.05658) have led to significant developments in achieving stable Helmholtz solution approximations by plane wave superposition. The study shows that the numerical instability and ill-conditioning inherent in plane wave-based Trefftz methods can be effectively overcome with regularization techniques, provided there exist accurate approximations in the form of expansions with bounded coefficients. While propagative plane waves fail to yield stable approximations due to the exponential growth of the expansion coefficients, evanescent plane waves, which contain high-frequency modal content, provide both accurate and stable results. The developed numerical approach -- available at [https://github.com/EmileParolin/evanescent-plane-wave-approx](https://github.com/EmileParolin/evanescent-plane-wave-approx) -- results in substantial improvements when compared to conventional propagative plane wave schemes.

The numerical recipe proposed here has been tailored for the 3D setting and extended with new sampling strategies involving extremal system of points (see [https://web.maths.unsw.edu.au/~rsw/Sphere/](https://web.maths.unsw.edu.au/~rsw/Sphere/)). The results show the desired accuracy and bounded-coefficient stability, in line with the two-dimensional case.

# Documentation

The MATLAB
