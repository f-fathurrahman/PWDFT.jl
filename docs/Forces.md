# Calculation of ionic forces

These notes are adapted from Kohanoff2006.

The forces on the nuceli of atomic species $s$ are given by
$$
\mathbf{F}_{I}^{s} = -\frac{\partial E_{\mathrm{KS}}}{\partial\mathbf{R}_{I}^{s}}
$$
Specifically:
$$
\mathbf{F}_{I}^{s} = \int \mathbf{\rho}(\mathbf{r})\left(V^{\mathrm{Ps,loc}}_{s}\right)'
\left(\left| \mathbf{r} - \mathbf{R}_{I} \right|\right)
\frac{\mathbf{r} - \mathbf{R}_{I}}{\left| \mathbf{r} - \mathbf{R}_{I} \right|}
\,\mathrm{d}\mathbf{r} +
\frac{Z_{I}}{2} \sum_{J \neq I} Z_{J}
\frac{\mathbf{R}_{I} - \mathbf{R}_{J}}{\left| \mathbf{R}_{I} - \mathbf{R}_{J} \right|^3} +
\mathbf{F}^{\mathrm{Ps,nloc}}_{s,I}
$$

$\left( V^{\mathrm{Ps,loc}}_{s} \right)'$ is the derivative of the local part of the pseudopotential

$\mathbf{F}_{\mathrm{Ps,nloc}}^{s}$: nonlocal pseudopotential contribution to the force

Using plane wave basis set.

Nuclear-nuclear contribution:
$$
F^{NN}_{s,I} = \frac{Z_{I}}{2}\sum_{J \neq I}^{P} Z_{J}
\sum_{L=-L_{\mathrm{max}}}^{L_{\mathrm{max}}}
\left( \mathbf{R}_{I} + l - \mathbf{R}_{J} \right) \times
\left[
\frac{\mathrm{erfc}\left(\left|\mathbf{R}_{I} + l - \mathbf{R}_{J} \right|\eta\right)}{\left|\mathbf{R}_{I} + l - \mathbf{R}_{J} \right|^3} +
\frac{\eta \mathrm{e}^{-\eta^2}\left|\mathbf{R}_{I} + l - \mathbf{R}_{J} \right|^2}
{\left|\mathbf{R}_{I} + l - \mathbf{R}_{J} \right|}
\right]
$$
Local pseudopotential contribution:
$$
F^{\mathrm{Ps,loc}}_{s,I} = \Omega \sum_{\mathbf{G} \neq \mathbf{0}}
\mathrm{i}\mathbf{G}\,\mathrm{e}^{\mathbf{G}\cdot\mathbf{R}_{I}}
V^{\mathrm{Ps,loc}}_{\mathrm{s}}(G) \rho(\mathbf{G})
$$
