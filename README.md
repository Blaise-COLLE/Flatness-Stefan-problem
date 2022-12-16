
> Dirichlet boundary control for the two phases Stefan problem using flatness

+ [Context](#Context)
+ [How to tune your own simulation](#How-to-tune-your-own-simulation)
+ [Documentation](#Documentation)
+ [Licence](#Licence)
+ [Authors](#Authors)

---

# Context

This provides codes related to the paper [Controllability of the Stefan problem by the flatness approach](https://hal.archives-ouvertes.fr/hal-03721544) by Blaise Colle, Jérôme Lohéac and Takéo Takahashi.

The scripts have been tested on MATLAB Version 9.9 (R2020b) and uses the "Parallel Computing" and "Symbolic Math" Toolboxes of MATLAB.

---

# How to tune your own simulation

Run `example` on MATLAB.

Parameters on interest are:
- `cs`, `cl`: diffusivity coefficients in the solid and liquid phases 
- `ntrap`:    number of points for some integral evaluations (using the trapeze method) 
- `N`:        truncation order in the sums 
- `nt`:       number of points for the time evaluation of the solution 
- `nx`:       number of points for the spacial evaluation of the temperature 
- `T`:        controllability time (**T>0**) 
- `sigma`:    Gevrey order (**1<sigma<2**) 
- `b0`, `v0`: define the initial state (theta0(x)=v0(x-b0)) 
- `b1`, `v1`: define the final (target) state (theta1(x)=v1(x-b1)) 

Once the simulation has ended, the results are stored in:
- `t_`:       the time evaluations
- `us`, `ul`: the two controls at times `t_`
- `bs`:       the position at of the solid/liquid interface at times `t_`
- `xs`, `xl`: the abscissa for the spatial evaluation at times `t_`, `xs(:,i)=linspace(0,bs(i),nx)'` and `xl(:,i)=linspace(bs(i),1,nx)'`
- `thetass`, `thetals`: the temperatures in the two phases, for instance `thetass(i,j)=theta(t_(j),xs(i,j))`

---

# Documentation

A brief html documentation can be obtained by running `Contents` (this requires [makehtmldoc](http://www.fast.u-psud.fr/~moisy/ml/)).
The documentation will then be available at `./Documentation/Contents.html`.

---

# Licence

CeCILL-C Free Software License Agreement

---

# Authors

Blaise Colle, Jérôme Lohéac and Takéo Takahashi
