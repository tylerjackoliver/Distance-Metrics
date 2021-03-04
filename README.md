# Distance-Metrics

In supporting ongoing research into dynamical systems theory, particularly in the use of Lagrangian Coherent Structures for flow analysis and automatic orbit classification in astrodynamics systems, efficient metrics for determining distances between two trajectories are required.

This repository stores two such commonly-used metrics: the Hausdorff distance ![(_here_)]([https://en.wikipedia.org/wiki/Hausdorff_distance]) and the Frechet distance ![(_here_)](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance). 

## Hausdorff Distance

The Hausdorff algorithm is a direct port of the SciPy directed Hausdorff distance ![found here](https://github.com/scipy/scipy/blob/v1.6.1/scipy/spatial/distance.py#L365-L462), and it is the greatest of all the distances from a point in one set to the closest point in the other set. This algorithm has been slightly optimised for Intel-based CPUs by writing the program in a manner for targeting CPU vectorisation using the Intel C++ compiler, but no other optimisations have been made.

## Frechet Distance

The Frechet distance is an implementation of the following two papers:

   * Devogele, T., Esnault, M., Etienne, L., & Lardy, F. (2017). Optimized Discrete Fréchet Distance between trajectories. BigSpatial 2017 - Proceedings of the 6th ACM SIGSPATIAL International Workshop on Analytics for Big Geospatial Data, (November), 11–19. https://doi.org/10.1145/3150919.3150924
   * João Paulo Figueira. (2021). Fast Discrete Fréchet Distance.

Both papers present highly optimised versions of the Frechet distance computation; further improvements have been added to support automatic CPU-directed vectorisation. The C++17 standard and the C++ STL is used exclusively.
