## Target of the approximation problem

We consider the Helmholtz solution in the unit ball with *wavenumber* `k` such as

````
k = 5;
````

We define the *maximum mode number* `L` in the approximation target:
````
L=5*k;
````
The *solution surrogate* `u` is defined as a superposition of *spherical waves* $b_{\ell}^m$, namely
```math
u=\sum_{\ell=0}^{\infty}\sum_{m=-\ell}^{\ell}u_{\ell}^m\,b_{\ell}^m.
```
We need to define the modal content of the solution surrogate `u` by specifying the coefficients in the spherical wave basis $b_{\ell}^m$ for $0\leq|m|\leq \ell$. This is done by defining a *modal matrix* `M` whose lines are of the form $(\ell,m,u_{\ell}^m)$, where $\ell$ is the *degree*, $m$ is the *order*, and $u_{\ell}^m$ is the associated coefficient.
````
M=[];
for l=0:L
    ulm=randn(2*l+1,1)/max(1,l-k);
    M=[M;l*ones(2*l+1,1),(-l:l)',ulm];
end
````
The coefficients $u_{\ell}^m$, for $0\leq|m|\leq \ell$, within the expansion are products of normally-distributed random numbers (with mean 0 and standard deviation 1) and the scaling factors $(\text{max}(1,\ell-\kappa))^{-1}$. Since the coefficients of any Helmholtz solution decay in modulus as $o(\ell^{-1})$ for $\ell \rightarrow \infty$, this scenario is quite challenging.


The target of the approximation problem `u` can then be constructed as:
````
u=@(x)sol_surrogate(M,k,x);	
````

