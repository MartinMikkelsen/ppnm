My student number ends with 71 which means I have to do exam problem number 5 (71 mod 22 = 5).

Adaptive 1D integrator with random nodes.

Implement an adaptive one-dimensional integrator with random abscissas. Reuse points. Note that you don't need to pass the points to the next level of recursion, only statistics. 

As mentioned in the lecture notes "Yet Another Introduction to Numerical Methods". Classical quadratures are mostly of historical interest. Alternative methods could be quadratues with optimized abscissas, adaptive or variable transformation variables. In this exam project I will consider another case with random abscissas closely related to Monte Carlo integration -- a quadrature where the abscinnas are chosen randomly. The adaptive part implies that the distribution of points is chosen non-uniformly to reduce integration errors. 

An example of a recursive adaptive integration algorithm is the Stratified Sampling Algorithm (SSA.c) which is implemented here. GSL also has an stratified sampling algorithm called MISER. This technique aims to reduce the overall integration error by concentrating integration points in the regions of highest variance.

In this exam project I compare my implementation of the Stratified Sampling Algorithm to other routines related to one-dimensional integrators. I consider the Recursive Adaptive Integrator (Recursive_Adaptive_Integrator.c), Clenshaw–Curtis variable transformation (Clenshaw–Curtis.c),  plain Monte Carlo (pMonteCarlo.) and quasi-random Monte-Carlo integrator (qMonteCarlo.o). 

I the main.c I have included some different integrals, but the output.txt prints data using 

\begin{equation}
\int_{0}^{1} dx 4\sqrt{1-x²} = \pi,
\end{equation}

