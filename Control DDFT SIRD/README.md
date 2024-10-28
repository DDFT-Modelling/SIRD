# Exploring nonpharmaceutical controls

$$
\begin{align}
	\rho &= \big( \begin{smallmatrix} S \\ I \\ R \end{smallmatrix} \big)
	\\
	\frac{\partial S}{\partial t} &=
            \nabla \cdot \Big( D_S \nabla S - \Gamma_S S \nabla \big( K_{\text{sd}} \star (S+R) +  K_{\text{si}}
                                \star I\big)  \Big) - \beta SI
	\\
	\frac{\partial I}{\partial t} &= D \Delta I + \beta SI - \gamma I
	\\
	\frac{\partial R}{\partial t} &= D \Delta R + \gamma I
\end{align}
$$

---
## Main files

[`FISTA_SIRD_DDFT_box.mlx`](FISTA_SIRD_DDFT_box.mlx)  runs FISTA and determines a pair of controls.


[`InitialCondition.m`](InitialCondition.m) contains three different initial conditions, can be selected using the keyword `choice`.

[`Kernels.m`](Kernels.m) defines two Gaussian kernels.

[`State.m`](State.m) returns the vector field $\rho(x;t)$ from initial data and parameters.

[`Adjoint.m`](Adjoint.m) returns the vector field $q(x;t)$ from terminal data, $\rho(x;t)$, and parameters.

