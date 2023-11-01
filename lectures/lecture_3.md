---
theme: presentation
marp: true
math: katex
---

<!-- _class: titlepage-->

<div style="flex-basis: 100%">

# Scientific Computing for Geophysical Problems

## Noise-based virtual seismograms

Léonard Seydoux
Institut de physique du globe de Paris
Université Paris Cité
seydoux@ipgp.fr

</div>

<div>

![bg brightness:0.6 sepia](images/cover/hurricane_isabel_from_iss.jpg)

</div>
<span class="logo">

</span>

---

<!-- paginate: true -->

## Foreword

| __When__ | __Who__ | __What__ | __Topic__|
| :- | :- | :- | :- |
| 9/14 | AF | Lecture | Planetary magnetic fields |
| 9/21 | AF | Lab | Basic signal processing with python and application to geomagnetic data |
| 16/21 | LS | Lecture & lab | Introduction to observational seismology |
| 10/5 | LS + AF | Lecture & lab | Seismograms, earthquakes, travel times, hypocenter |
| 10/12 | AF | Lab | Models of the secular variation of the Earth’s magnetic field part I |
| 10/19 | AF | Lab | Models of the secular variation of the Earth’s magnetic field part II |
| 10/26 | LS | Lecture & lab | Earth’s structure and seismic sources |
| Today | LS | Lecture & lab | Basic properties of earthquake sources, seismic tremors, and noise |
| TBD | LS + AF | Final exam (2h) | Observational seismology |

<!-- _footer: __AF__: Alexandre Fournier <br> __LS__: Leonard Seydoux -->

---

## Digital seismological observations

<div>

### Nowadays seismic data

- 10k+ permanent
- 1000k+ temporary
- High-sensitivity, big data
- Open-access$^*$

![width:400px](images/networks/network_sjf.png)

</div>
<div>

![width:400px](images/networks/network_gsn.png)
![width:400px](images/networks/network_usarray.png)

</div>

<!-- _footer: $^*$ Depending on funding organizations. Global Seismogarphic Network, US Transportable Array, and San Jacinto seismic array.-->

---

## Digital seismological observations

![](images/seismograms/subduction_eq.png)

---

## Digital seismological observations

![](images/seismograms/subduction_no.png)

---

## Digital seismological observations

![](images/seismograms/subduction_tr.png)

---

## Digital seismological observations

![](images/seismograms/subduction_vo.png)

---

## Summary of seismological observations

<div>

![width:400px](images/seismograms/subduction_eq_signal.png)
![width:400px](images/seismograms/subduction_all.png)

</div>
<div>

### External sources

- Oceanic noise (continuous)
- Atmospheric noise (high frequency)
- Antropogenic noise (modulated)

### Internal sources

- Earthquakes (<1% of data)
- Tectonic tremors (> 1 Hz)
- Volcanic tremors (> 1 Hz)

</div>

---

## Forward and inverse problems in seismology

<div>

### Forward problem

A displacement seismogram $u(t)$ convolves a source $s(t)$ with the medium $g(t)$ and a instrument $i(t)$, such as

$$u(t) = (s * g * i )(t).$$

</div>
<div>

### Inverse problems: seismic imaging

Find $g$ from a source $u$ assuming earthquake properties $s$ and instrument response $i$, such as

$$g = u * (s * i)^{-1}$$

</div>

---

## Issues with earthquake-based imaging

<div style="flex-basis: 45%;">

- Earthquakes do not occur everywhere
- Earthquakes do not occur continuously
- Earthquakes rarely occur at the same place

### → Limited resolution and monitoring

### → What about using ambient noise

</div>
<div style="flex-basis:30%;">

![](images/seismograms/earthquake_location.png)

</div>

---

## Noise correlation theorem

<div align="center" style="font-size: smaller;">

Assuming a __random__ wavefield issued by sources distributed __homogeneously__ in space in the medium,

$$\frac{d}{d\tau}C_{ij}(\tau, \mathbf{r}_i, \mathbf{r}_j) \approx -\frac{\sigma^2}{4a} \left( G({\tau, \mathbf{r}_i, \mathbf{r}_j}) - G({-\tau, \mathbf{r}_i, \mathbf{r}_j}) \right)$$

where $C_{ij}(\tau)$ is the noise cross-correlation between the two locations $\mathbf{r}_i$ and $\mathbf{r}_j$, and $G({\tau, \mathbf{r}_i, \mathbf{r}_j})$ is the medium's __Green's function__ between those locations.

![width:600px](images/seismograms/crosscorrelation.png)

</div>

---

## New way for seismic imaging

<div>

### Forward problem

A displacement seismogram $u(t)$ convolves a source $s(t)$ with the medium $g(t)$ and a instrument $i(t)$, such as

$$u(t) = (s * g * i )(t).$$

</div>
<div>

### Inverse problems: seismic imaging

~~Find $g$ from a source $u$ assuming earthquake properties $s$ and instrument response $i$, such as~~

The data correlation converges towards the imaginary part of the Green's function __under noise properties assumtions__

$$\frac{1}{2}(g - g^*) \approx u u^*$$

</div>

---

## New way for seismic imaging

<div>

- We can turn every seismometer into a virtual source
- Noise sources are continuous in time
- We now have tons of useful data!

![width:400px](images/networks/network_sjf.png)

</div>
<div>

![width:400px](images/networks/network_gsn.png)
![width:400px](images/networks/network_usarray.png)

</div>

---

## First noise-based high-resolution surface-wave tomography

![](images/tomography/shapiro_2005.jpeg)

<!-- _footer: Shapiro et al. (2005)-->

---

## Example: USArray noise-based surface-wave tomography

![](images/tomography/fan-chi-lin_2009.png)

<!-- _footer: Fan-Chi Lin et al. (2009)-->

---

## Example: global noise-based surface-wave tomography

![](images/tomography/nishida_2009.jpg)

<!-- _footer: Nishida et al. (2009)-->

---

## Example: high-resolution tomography at Valhall

![height:460px](images/tomography/valhall_map.jpeg)

![height:500px](images/tomography/valhall_tomo.jpeg)

<!-- _footer: Mordret et al. (2013)-->

---

## Example: anisotropic tomography of the Toba volcano, Indonesia

![](images/tomography/jaxibulatov_2014.jpg)

<!-- _footer: Jaxybulatov et al. (2014)-->

---

## Example: monitoring of the Parkfield postseismic relaxation

![height:500px](images/tomography/brenguier_1.jpeg)

![width:600px](images/tomography/brenguier_2.jpeg)

<!-- _footer: Brenguier et al. (2008)-->

---

## Example: monitoring during the Tohoku-Oki earthquake

<div>

![](images/tomography/brenuiger_3.gif)

</div>
<div>

How would you proceed without a continuous observable?

</div>

<!-- _footer: Brenguier et al. (2014)-->

---

## Noise correlation theorem assumptions

<div align="center">

Assuming a __random__ wavefield issued by sources distributed __homogeneously__ in space in the medium,

$$\frac{d}{d\tau}C_{ij}(\tau, \mathbf{r}_i, \mathbf{r}_j) \approx -\frac{\sigma^2}{4a} \left( G({\tau, \mathbf{r}_i, \mathbf{r}_j}) - G({-\tau, \mathbf{r}_i, \mathbf{r}_j}) \right)$$

where $C_{ij}(\tau)$ is the noise cross-correlation between the two locations $\mathbf{r}_i$ and $\mathbf{r}_j$, and $G({\tau, \mathbf{r}_i, \mathbf{r}_j})$ is the medium's __Green's function__ between those locations.

### __What if we do not meet the assumtions?__

</div>

---

## Noise correlation theorem assumptions

<div style="flex-basis:30%;">

Assuming a __random__ wavefield issued by sources distributed __homogeneously__ in space in the medium.

### __A random noise should have a flat spectrum.__

😬

</div>
<div>

![height:800px](images/seismograms/ppsd.png)

</div>

---

## Noise correlation theorem assumptions

<div style="flex-basis:30%;">

Assuming a __random__ wavefield issued by sources distributed __homogeneously__ in space in the medium.

### __A random noise should be stationnary.__

😬

</div>
<div>

![height:300px](images/seismograms/signal_pet.jpg)

</div>

---

## Noise correlation theorem assumptions

<div style="flex-basis:30%;">

Assuming a __random__ wavefield issued by sources distributed __homogeneously__ in space in the medium.

### __Sources are distributed in coastal areas.__

😬

</div>
<div>

![](images/seismograms/generation.png)

</div>

<!-- _footer: Seydoux et al. (2016)-->

---

## Noise correlation theorem assumptions

<div style="font-size:smaller;" align="center">

Assuming a __random__ wavefield issued by sources distributed __homogeneously__ in space in the medium.

### __Effect of source distribution.__

![](images/seismograms/distribution.png)

</div>

<!-- _footer: from Stehly et al. (2007)-->

---

## Pre-processing data to meet the requirements

![](images/seismograms/Illustration-of-the-pre-processing-of-seismic-records-a-Raw-data-b-Raw-data.png)

<!-- _footer: Bensen et al. (2007) / Seydoux et al. (2016)-->

---

## Pre-processing data to meet the requirements

![](images/seismograms/bensen_2007.png)

<!-- _footer: Bensen et al. (2007) -->

---

<!-- _class: titlepage-->
<!-- _backgroundColor: black-->

# Notebook time
