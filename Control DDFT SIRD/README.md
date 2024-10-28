# Exploring nonpharmaceutical controls

$$
\begin{align}
	\rho &= \big( \begin{smallmatrix} S \\ I \\ R \end{smallmatrix} \big)
	\\
	\frac{\partial S}{\partial t} &=
            \nabla \cdot \Big( D_S \nabla S - \Gamma_S S \nabla \big( u K_{\text{sd}} \star (S+R) + v K_{\text{si}}
                                \star I\big)  \Big) - \beta SI
	\\
	\frac{\partial I}{\partial t} &=
           \nabla \cdot \Big( D_I \nabla I - \Gamma_I I \nabla \big( v K_{\text{si}} \star (S+I+R) \big)  \Big)
                + \beta SI - [\gamma+m] I
	\\
	\frac{\partial R}{\partial t} &=
           \nabla \cdot \Big( D_R \nabla R - \Gamma_R R \nabla \big( u K_{\text{sd}} \star (S+R) +  v K_{\text{si}}
                                 \star I\big)  \Big) + \gamma I
\end{align}
$$

---
## Main files

[`FISTA_SIRD_DDFT_box.mlx`](FISTA_SIRD_DDFT_box.mlx)  runs FISTA and determines a pair of controls.


[`InitialCondition.m`](InitialCondition.m) contains three different initial conditions, can be selected using the keyword `choice`.

[`Kernels.m`](Kernels.m) defines two Gaussian kernels.

[`State.m`](State.m) returns the vector field $\rho(x;t)$ from initial data and parameters.

[`Adjoint.m`](Adjoint.m) returns the vector field $q(x;t)$ from terminal data, $\rho(x;t)$, and parameters.


---

## Outputs

[`Panels_Plot.m`](Panels_Plot.m) returns a density plot at user selected times.

[`plots_in_box.m`](plots_in_box.m) generates an animation for a vector field of dimension 3.

[`Plot_SIR_Mean_Curves.m`](Plot_SIR_Mean_Curves.m) aggregates each compartment in space and plots it.

