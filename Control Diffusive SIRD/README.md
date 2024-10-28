# Parameter identification for the diffusive SIR model

$$
\begin{align}
	\rho &= \big( \begin{smallmatrix} S \\ I \\ R \end{smallmatrix} \big)
	\\
	\frac{\partial S}{\partial t} &= D \Delta S - \beta SI
	\\
	\frac{\partial I}{\partial t} &= D \Delta I + \beta SI - \gamma I
	\\
	\frac{\partial R}{\partial t} &= D \Delta R + \gamma I
\end{align}
$$

---
## Main files

[`FISTA_DiffusiveSIR_nb.mlx`](FISTA_DiffusiveSIR_nb.mlx) runs FISTA and determines a pair of controls.

[`Forward_Diffusive_SIR.m`](Forward_Diffusive_SIR.m) computes the forward system and returns the vector field $\rho(x;t)$. It has an option for plotting. A simplified version is included in [`State.m`](State.m).

[`Adjoint_Diffusive_SIR.m`](Adjoint_Diffusive_SIR.m) computes the adjoint system and returns the vector field $q(x;t)$. It has an option for plotting. A simplified version is included in [`Adjoint.m`](Adjoint.m).

[`InitialCondition.m`](InitialCondition.m) contains three different initial conditions, can be selected using the keyword `choice`.

[`Objective.m`](Objective.m) evaluates the objective $j(\alpha) = \frac{1}{2} \| \rho - \widehat{\rho}\|^2$.

[`Gradient.m`](Gradient.m) evaluates $\nabla j(\alpha)$.

---

## Outputs

[`Panels_Plot.m`](Panels_Plot.m) returns a density plot at user selected times.

[`plots_in_box.m`](plots_in_box.m) generates an animation for a vector field of dimension 3.

[`Plot_SIR_Mean_Curves.m`](Plot_SIR_Mean_Curves.m) aggregates each compartment in space and plots it.

[`Contour_Plots.m`](Contour_Plots.m) generates a coarse and a fine contour plot.
