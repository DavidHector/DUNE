# Notes
## Basic understanding
### Finite Element Method
- Used for solving parial differential equations
- Works as follows:
    1. Subdivide large system into smaller parts (finite Elements)
    2. Approximates unknown function over space by minimizing associated error function (common example: Poisson ($ - \Delta u = f \text{ in } \Omega$) equation with Dirichlet boundary condition ($ u = g \text{ in } \partial \Omega$). Inhomogenous Dirichlet: $u = g$, homogeneous: $g=0$ )
        1. convert function to weak formulation: multiply both sides with smooth function $v$ that satisfies boundary condition and integrate over the domain. Then do partial integration and obtain a formula like 
        $$ \int_{\Omega} \nabla u \cdot \nabla v \ dx = \int_{\Omega} fv \ dx \\ 
        a(u, v) = I(v)
        $$
        These are the bilinear and linear forms, also called the weak formulation. To work u and v must be part of the [Sobolev space](https://en.wikipedia.org/wiki/Sobolev_space) $H^1$. An example for v is a piecewise linear function. Also the space of the functions u is typically represented by Lagrange basis functions: 
        $$
        \Phi_h = \{\phi_1, ..., \phi_N\}, \quad \forall i, j \in \mathcal{I}_h: \phi_i(x_j) = \delta_{i, j}
        $$
        where $ \mathcal{I}$ is the index space of the mesh. If this basis representation is inserted in the bilinear form, the integral becomes a sum with coefficients $z_j$ (that form the vector z), by expressing u as:
        $$
        u_h = \sum^{N}_{j=1} z_j\phi_j
        $$
        The same is done for v, and putting this in the bilinear form, we eventually get a matrix problem:
        $$
        Az = b\\
         A_{i, j} = a(\phi_j, \phi_i), \ b_i = I(\phi_i) \text{ inside the space}\\
        A_{i,j} = \delta_{i, j}, \ b_i = z_i = g(x_i) \text{ on the boundary}
        $$
        This can then be solved using several methods, like Gaussian elimination or Iterative methods like *Richardson's* iteration:
        $$
        z^{k+1} = z^k + w(b-Az^k), \quad w>0
        $$
        Other solvers: Krylov, etc. 

### Mesh
Always remember the __integrationFactor__, i.e. the shape function transformation when integrating: Don't just use $f(x) \cdot w $, use $f(\hat{x}) \cdot w \cdot \mu(\hat{x})$, where $\hat{x}$ is the quadrature point.  