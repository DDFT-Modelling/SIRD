# Exploring nonpharmaceutical controls



---


[`InitialCondition.m`](InitialCondition.m) contains three different initial conditions, can be selected using the keyword `choice`.

[`Kernels.m`](Kernels.m) defines two Gaussian kernels.

[`State.m`](State.m) returns the vector field $\rho(x;t)$ from initial data and parameters.

[`Adjoint.m`](Adjoint.m) returns the vector field $q(x;t)$ from terminal data, $\rho(x;t)$, and parameters.

[`FISTA_SIRD_DDFT_box.mlx`](FISTA_SIRD_DDFT_box.mlx) 