---
theme: presentation
marp: true
math: katex
---

<!-- _class: titlepage-->

<div style="flex-basis: 100%">

# Scientific Computing for Geophysical Problems

### Introduction to observational Seismology

by Léonard Seydoux
`seydoux@ipgp.fr`

<!-- _footer: Institut de physique du globe de Paris, Université Paris Cité -->

</div>
<div>

![bg brightness:0.6 sepia](images/cover/synthetic_waves.jpg)

</div>
<span class="logo"></span>

---

## Simple oscillator

<div>

__Hooke's law__ models the force $F$ needed to move a spring by a displacement $u$ by

$$F = ku$$

where $k$ is the spring constant (i.e. stiffness).

__This model is valid for small displacements (linear elasticity).__

</div>
<div style="flex-basis: 20%">

![](images/waves/oscillator.gif)

</div>

---

## Simple oscillator

<div>

__Newton's law of motion__ models the position $u$ of a mass $m$ attached to a spring of constant $k$ with

$$m \ddot{u} = -ku.$$

This gives the oscillator equation

$$\ddot{u} + \omega^2 u = 0$$

where $\omega = \sqrt{k/m}$ is the resonnant frequency, with simple solution $u(t) = a \cos(\omega t)$ valid only for __small deformations__ and __no attenuation__.

</div>
<div style="flex-basis: 20%">

![](images/waves/oscillator.gif)

</div>

---

## 1D wave equation

<div>

Newton's law applied to the mass $m_i$ of a __spring-coupled mass array__ gives
$$\ddot{u}_i = \frac{k}{m}(u_{i+1} - 2u_i + u_{i-1}).$$
Considering a mass distance $h$, mass density $\rho = m/h$ and elastic constant $E = kh$ leads to
$$\ddot{u} = \frac{E}{\rho}\frac{u_{i+1} - 2u_i + u_{i-1}}{h^2}.$$

</div>

<div>

![width:700px](images/waves/oscillator_coupled.png)

</div>

---

## 1D wave equation

<!-- ![bg right 90%](Images/spring_couple_masses.png) -->

<div>

Now considering $h \rightarrow 0$, we identify the second-order space derivative of the displacement:

$$\lim_{h \rightarrow 0} \frac{u_{i+1} - 2u_i + u_{i-1}}{h^2} = u''(x).$$

Back into the equation, we obtain the 1D wave equation:

$$\ddot{u}(x, t) = c^2 u''(x, t)$$

where $c = \sqrt{E / \rho}$ is the wavespeed.

</div>
<div>

![width:700px](images/waves/oscillator_coupled.png)

</div>

---

## 1D wave equation

<div>

The displacement $u(x, t)$ at a given position $x$ and time $t$ is given by

$$\frac{\partial^2 u(x, t)}{\partial x^2} = \frac{1}{c^2}\frac{\partial^2 u(x, t)}{\partial t^2}$$

This equation admits solutions of the form $u(x, t) = f(x \pm vt)$, which describe a propagating pulse __in both directions__ along the space dimension.

</div>
<div>

_One-sided solution_

![](stein_and_wysession/chapter2-figures/2_2_02.JPG)

</div>

<!-- _footer: Stein and Wysession, Figure 2.2-2.-->

---

## 1D wave equation

<div>

Depending on the boundary conditions in space and time, one solution to the 1D wave equation is

$$u(x, t) = A e^{i(\omega t \pm k x)}$$

where $\omega$ is the angular frequency ($2\pi f$), and $k$ is the wavenumber. We also have $v = \omega / k$.

__This solution is known as the harmonic solution.__

</div>
<div>

![](stein_and_wysession/chapter2-figures/2_2_03.JPG)

</div>

---

## 1D wave equation

<div style="font-size:smaller">

__Relationship between wave variables__

| Quantity | Units | Relationship |
|:-|:-|:-|
| Velocity | m/s | $v = \omega / k$|
| Period | s | $T = 2 \pi / \omega$|
| Angular frequency | rad/s | $\omega = 2 \pi f$|
| Frequency | Hz | $f= 1 / T$|
| Wavelength | m | $\lambda= v / f$|
| Wavenumber | rad/m | $k = 2 \pi / \lambda$|

</div>
<div>

![](stein_and_wysession/chapter2-figures/2_2_04.JPG)

</div>

---

## 3D wave equation

<div>

We need to turn everything into three dimensions:

- Stress tensor
- Strain tensor
- Equations of motion
- Constitutive equations

</div>
<div>

![width:800px](images/waves/3D_wave.gif)

</div>

---

<div style="padding: 40px; background-color: #FFFFFF99;">

In continuum mechanics, the basic concept of a continuous homogeneous media is an __approximation__. Real Earth __is__ heterogeneous.

</div>

![bg sepia](images/cover/earth_is_heterogeneous.png)

---

## 3D wave equation

<div>

We transform the quantities __per unit volume__, such as the mass and the applied forces:

$$ \mathbf F = m \mathbf a \longrightarrow  \mathbf f ( \mathbf x, t) = \rho \ddot{\mathbf{u}}( \mathbf x, t)$$

where $\rho$ is the volumic mass, and the __bold typeface__ denotes a vector quantity.

</div>
<div>

![](stein_and_wysession/chapter2-figures/2_3_02.JPG)

</div>

---

## Stress tensor

<div>

The equation of motion writes

$$\partial_j \sigma_{ij}(\mathbf x, t) + f_i(\mathbf x, t) = \rho \ddot{u}_i(\mathbf x, t)$$

where $i, j = 1 \ldots 3$, and $\partial_j$ denotes the partial derivative along direction $j$, as in $\partial/\partial x_j$, and $\partial^2_t = \partial^2 / \partial t^2$.

Finally, $\sigma_{ij}$ denotes the __stress__ tensor.

</div>
<div>

![](stein_and_wysession/chapter2-figures/2_3_08.JPG)

</div>

<!-- _footer: Enstein convention used with implicit summation over repeated indices.-->

---

## Strain tensor

<div>

The __strain__ tensor describes the deformation resulting from the differential motion within a body $u_i(\mathbf{x} + \delta \mathbf{x}) \approx u_i(\mathbf{x}) + \delta u_j$, with the deformation

$$\delta u_i = \partial_j u_i(\mathbf{x})\delta x_j = (\epsilon_{ij} + \omega_{ij})\delta x_j$$

where $\epsilon_{ij}$ is the __strain__ tensor, and $\omega_{ij}$ is the rigid-body rotation without deformation. Valid for __small displacements__

</div>
<div>

![](stein_and_wysession/chapter2-figures/2_3_11.JPG)

</div>

---

## Strain tensor

<div>

$$
\epsilon_{ij} = \begin{pmatrix} \partial_1 u_1 & \dfrac{1}{2}(\partial_2 u_1 + \partial_1 u_2) & \dfrac{1}{2}(\partial_3 u_1 + \partial_1 u_3) \\ \dfrac{1}{2}(\partial_2 u_1 + \partial_1 u_2) & \partial_2 u_2 & \dfrac{1}{2}(\partial_3 u_2 + \partial_2 u_3) \\ \dfrac{1}{2}(\partial_3 u_1 + \partial_1 u_3) & \dfrac{1}{2}(\partial_3 u_2 + \partial_2 u_3) &\partial_3 u_3  \end{pmatrix}
$$

</div>
<div>

![](stein_and_wysession/chapter2-figures/2_3_11.JPG)

</div>

---

## Constitutive equations

<div>

Constitutive equations relate the stress $\sigma_{ij}$ to the strain $\epsilon_{ij}$. We distinguish:

- The __linear elasticity__, which implies a linear relationship between stress and strain.
- The __non-linear elasticity__, which can be viscous, visco-elastic, elastic-plastic, etc.

__We will consider linear elasticity.__

</div>
<div>

![](images/waves/elasticity.webp)

</div>

---

## Linear elasticity

<div>

The constitutive equation is the __Hooke's law__, again valid for small displacements

$$\sigma_{ij} = c_{ijkl} \epsilon_{kl}.$$

Here, the constants in the tensor $c_{ijlk}$ are the __elastic moduli__ describing the material's properties. In 3D, $c_{ijkl}$ contains 3⁴ = 81 independant coefficients.  

> 81 independant constants?

The simmerties imply $c_{ijkl} = c_{jikl} = c_{jilk}$. The number of independant constants is 36.

> 36 independant constants, still?

Wait, there is a last symmetry relation based on the idea of strain energy: $c_{ijkl} = c_{klij}.$

> We need 21 independant components to describe an anisotropic elastic medium.

</div>

---

## 3D wave equation

<div align="center">

Finally, we can write the elastodynamic equation of motion as

$$(c_{ijkl} \partial_l u_k)_j (\mathbf x, t) + f_i (\mathbf x, t) = \rho \ddot{u}_i (\mathbf x, t).$$

</div>

---

## 3D wave equation for istropic media

<div>

If the material behaves the same regardless of orientation, is it __isotropic__. This reduces the number of inependant elastic constants to __2__.

### Lamé constants

$c_{ijlk}=\lambda \delta_{ij}\delta_{kl} + \mu (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})$
$\sigma_{ij}=\lambda e_{kk}\delta_{ij} + 2 \mu e_{ij}$

</div>

<div>

We re-formulate the 3D wave equation in isotropic media

$$(\lambda + \mu)\nabla(\nabla \cdot \mathbf{u}(\mathbf{x}, t)) + \mu \nabla^2 \mathbf{u}(\mathbf{x}, t) = \rho \ddot{\mathbf{u}}(\mathbf{x}, t)$$

</div>

---

## Wave equations for _P_ and _S_ waves

Rewriting $\mathbf{u}(\mathbf{x}, t) = \nabla\phi(\mathbf{x}, t) + \nabla \times \psi(\mathbf{x}, t)$, we get two different independant equations

<div>

P wave equation:
$\nabla^2\phi(\mathbf{x}, t) = \frac{1}{\alpha^2} \ddot{\phi}(\mathbf{x}, t)$

S wave equation:
$\nabla^2\psi(\mathbf{x}, t) = \frac{1}{\beta^2} \ddot{\psi}(\mathbf{x}, t)$

with _P_ wave velocity $\alpha = \sqrt{(\lambda + 2\mu)/\rho}$ and _S_ wave velocity $\beta = \sqrt\frac{\mu}{\rho}$

</div>
<div align="center">

![](images/waves/seismic_waves.gif)

</div>

---

## Plane-wave solution

<div>

Correspond to the solution of the 1D scalar wave equation:

$$u(x, t) = Ae^{i(\omega t \pm kx)}$$

</div>

<div>

![](stein_and_wysession/chapter2-figures/2_4_01.JPG)

</div>

---

## Spherical solution

<div>

We can write the scalar wave equation in spherical coordinates

$$\frac{1}{v^2}\ddot{\phi}(r, t) = \nabla^2 \phi(r, t) = \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial \phi(r, t)}{\partial r}\right)$$

A solution is any function of the form

$$\phi(r, t) = \frac{f(t \pm r/v)}{r}$$

Geometric spreading: amplitudes decay as $1/r$ energy decays as $1/r^2$.

</div>

<div>

![](stein_and_wysession/chapter2-figures/2_4_02.JPG)

</div>

---

## Single force source

<!-- _footer: From Aki and Richards (chap. 4) -->

<div>

$$u_i({\bf x}, t) = \frac{1}{4\pi \rho}(3\gamma_i \gamma_j - \delta_{ij})\frac{1}{r^3}\int_{r/\alpha}^{r/\beta}\tau x_0(t - \tau) d\tau + \frac{1}{4\pi \rho \alpha^2}\textcolor{orange}{\gamma_i\gamma_j}\frac{1}{r}x_0\left(t - \frac{r}{\alpha}\right) - \frac{1}{4\pi \rho \beta^2}\textcolor{orange}{(\gamma_i\gamma_j - \delta_{ij})}\frac{1}{r}x_0\left(t - \frac{r}{\beta}\right) $$

![](images/waves/radiation_pattern.png)

</div>

---

## Homogeneous travel times of the far-field P and S waves

<div>

$$T(x) = \frac{r}{c}$$

</div>

<div>

![](stein_and_wysession/chapter2-figures/2_4_02.JPG)

</div>

---

## Travel times of the far-field P and S waves

<div>

1. Wavefront still exist
2. Local plane wave:\
$\phi(x, t) = A(x)e^{i\omega(t - T(x))}$
3. Put it into the scalar wave equation leads to
$$|\nabla T|^2 = \frac{1}{c^2} = \mathbf{s} $$

We call it the __Eikonal equation__, where $\mathbf{s}$ is the slowness vector normal to the wavefront (locally) $|\mathbf{s}| = 1/c$. Together with __ray tracing__, we can predict the travel times from a given source!

</div>

<div>

![](images/waves/ray_tracing.png)

</div>

---

## Finally, observed seismic phases on seismograms

![width:800px](images/seismograms/seismicwaves.webp)

---

## Finally, observed seismic phases

![](stein_and_wysession/chapter1-figures/1_1_03.JPG)

---

<!-- _class: titlepage-->
<!-- _backgroundColor: black-->

# Time to locate
