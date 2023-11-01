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

![bg brightness:0.6](images/cover/computer_4.gif)

</div>
<span class="logo"></span>

---

<!-- paginate: true -->

## Foreword



<img src="images/cover/me.jpeg"/>

<div style="flex-basis: 20%">

__Leonard Seydoux__
Assistant Professor
Institut de Physique du Globe de Paris
Université Paris Cité
seydoux@ipgp.fr

</div>
<div style="flex-basis: 65%">

__Research keywords__
Seismology, Volcanology, Artificial Intelligence, Geodesy, Applied Mathematics.

__Public websites__
[`Google Scholar`](https://scholar.google.fr/citations?user=TqLdn9YAAAAJ&hl) [`ResearchGate`](https://www.researchgate.net/profile/Leonard-Seydoux) [`GitHub`](https://github.com/leonard-seydoux)
</div>

---

## Foreword

| Class | __Who__ | __What__ | __Topic__|
| :- | :- | :- | :- |
| 1 | AF | Lecture | Planetary magnetic fields |
| 2| AF | Lab | Basic signal processing with python: application to geomagnetic data |
| 3 | AF | Lab | Models of the secular variation of the Earth’s magnetic field part I |
| 4 | AF | Lab | Models of the secular variation of the Earth’s magnetic field part II |
| 5 | LS | Lecture & lab | Introduction to observational seismology |
| 6 | LS | Lecture & lab | Seismograms, earthquakes, travel times, hypocenter |
| 7 | LS | Lecture & lab | Earth’s structure and seismic sources |
| 8 | LS | Lecture & lab | Basic properties of earthquake sources, seismic tremors, and noise |
| 9 | LS + AF | __Final exam__  | Observational seismology |

<!-- _footer: __AF__: Alexandre Fournier <br> __LS__: Leonard Seydoux -->

---

## Foreword: evaluation

<div>

- __50%__ of the grade is on observational seismology

- __25%__ on four in-class notebooks,\
5+1 points each, 20 points max

- __25%__ on final exam, 2 hours in-class with material

</div>
<div style="width: 80%;">

![width:600px](images/flowcharts/grading.svg)

</div>

<!-- _footer: __AF__: Alexandre Fournier <br> __LS__: Leonard Seydoux -->

---

<!-- _class: titlepage-->

# 1. Introduction

![bg brightness:0.6 grayscale](images/cover/synthetic_waves.jpg)

<!-- _footer: Synthetic displacement after an earthquake in Northern Italy. © Matthias Meschede.-->

---

## Objectives

<div>

- Get a __scientific method__ to pose and solve scientific problems

- Get __reference__ textbooks, papers, web resources, and softwares
- Get access to and analyze __data__
- Have tutorial-like __programs__ to help you later in your carrier
- Benefit from our experience by __asking questions__ (email, class)

</div>

<div>

![](images/flowcharts/scientific_method_geophysics.png)

</div>

---

## About __self-learning__

<div>

### Advantages

- Community-based
- Rapid understanding

</div>
<div>

### Disadvantages

- Can be executed "as-is"
- No learning involved

</div>
<div style="flex-basis:90%">

### Recommendations

- Read the __documentation__ and __additional information__ (websites, textbooks)
- Read until you can __explain it yourself__ to seomeone else
- Ask __questions__ (teachers, classmates) if something remains not clear

</div>

---

## References

<div>

### For a complete introduction

Seth Stein and Michael Wysession, An introduction to seismology, earthquakes, and earth structure. John Wiley & Sons, 2009.

</div>
<div>

![width:350px drop-shadow](images/references/cover_stein_and_wysession.png)

</div>

<!-- _footer: Free content at epsc.wustl.edu/seismology/book -->

---

## References

<div>

### In case you still have some appetite

Keiiti Aki and Paul G. Richards. Quantitative seismology. University Science Books, 2002.

> The overcomplete reference textbook. Here, we call it _"LE"_ Aki.

</div>
<div style="flex-basis: 30%">

![width:320px drop-shadow](images/references/cover_aki_and_richards.jpg)

</div>

<!-- _footer: Homepage at www.ldeo.columbia.edu/~richards/Aki_Richards.html -->

---

## References

<div  >

### To deal with seismological data

We will use ObsPy, an open-source Python library for processing seismological data:

- Parsers for common formats
- Clients to access data
- Tools for signal processing
- Tutorials and examples

</div>
<div  >

![](images/references/homepage_obspy.png)

</div>

<!-- _footer: Documentation at obspy.org -->

---

## References

<div  >

### To plot results

We will use Matplotlib, an open-source Python library for plotting data

- Matlab-like syntax
- Publication-quality plots
- Tutorials and examples

> Also consider other libraries! `Seaborn` `Plotly` `Altair` `Bokeh`

</div>
<div >

![](images/references/homepage_matplotlib.png)

</div>

<!-- _footer: Documentation at matplotlib.org -->

---

## References

<div  >

### For calculating and data handling

We will use NumPy, an open-source Python library for array manipulation and linear algebra

- Powerful (written in C)
- Well-documented and easy-to-use

</div>
<div >

![](images/references/homepage_numpy.png)

</div>

<!-- _footer: Documentation at numpy.org -->

---

## References

<div  >

### The suggested interface

I recommend using Jupyter, an open-source software for running code (not only Python) inside notebooks.

- Cell-style coding
- Report-oriented

> You have also Atom, VSCode that can execute Jupyter notebooks.

</div>
<div >

![](images/references/homepage_jupyter.png)

</div>

<!-- _footer: Documentation at jupyter.org -->

---

## References

<div  >

### For managing packages

I recommand using Anaconda, a freemium distribution of Python and R

- Package management
- Virtual environments
- Stable-ish releases

> Alternatives include pip, docker or virtualenv (non-exhaustive).

</div>
<div >

![](images/references/homepage_conda.png)

</div>

<!-- _footer: Documentation at www.anaconda.com -->

---

<!-- _class: titlepage-->

# 2. Forward and inverse problems

![bg brightness:0.6 grayscale](images/cover/synthetic_waves.jpg)

<!-- _footer: Synthetic displacement after an earthquake in Northern Italy. © Matthias Meschede.-->

---

## Scientific method (hypothesis testing)

<div class="columns-center">

<br/>

![](images/flowcharts/scientific_method.svg)

</div>
<div>

- Get __observations__.
- Make __hypotheses__, and __predictions__.
- __Test__ against observations
- Report __conclusions__, improve hypotheses

The set of valid hypotheses allow us to __model__ the observations.

</div>

<!-- _footer: Diagram from [en.wikipedia.org/wiki/Scientific_method](http://en.wikipedia.org/wiki/Scientific_method). -->

---

<!-- class: smaller-->

## Observations in seismology

<div>

__Seismology__ (/saɪzˈmɒlədʒi, saɪs-/; from Ancient Greek σεισμός (seismós) meaning "earthquake" and -λογία (-logía) meaning "study of") is the scientific <mark>study of earthquakes</mark> and the <mark>propagation of elastic waves through the Earth</mark> or through other planet-like bodies. It also includes studies of earthquake environmental effects such as tsunamis as well as diverse seismic sources such as volcanic, tectonic, glacial, fluvial, oceanic, atmospheric, and artificial processes such as explosions. A related field that uses geology to infer information regarding past earthquakes is paleoseismology. <mark>A recording of Earth motion as a function of time is called a __seismogram__.</mark> A seismologist is a scientist who does research in seismology.

![](images/seismograms/seismogram.png)

</div>

<!-- _footer: Definition from [en.wikipedia.org/wiki/Seismology](http://en.wikipedia.org/wiki/Seismology). <br> Daylong seismogram from the Black Forest Observatory from [examples.obspy.org](http://examples.obspy.org).-->

---

## Observations in seismology

<div style="text-align: center;">

Daylong seismogram recorded by the seismic station BFO. A lot of __knowledge__ was inferred from this type of recordings (Earth interior, earthquake mechanism).

![](images/seismograms/seismogram_long.png)

</div>

<!-- _footer: Daylong seismogram from the Black Forest Observatory from [examples.obspy.org](examples.obspy.org).-->

---

<!-- _backgroundColor: white -->

## Knowledge inferred from seismograms

<div style="text-align:center">

### Seismic wave types

![](images/waves/seismic_waves.gif)

</div>

<div style="text-align:center">

### Seismic phases and Earth structure

![width:450px](images/waves/seismic_phases.png)

</div>

<!-- _footer: Shear and compressional waves from `user194703` on [Stack Exchange](http://tex.stackexchange.com/questions/528579/draw-p-and-s-waves-illustration) <br> Sample seismic phases from the [Encyclopedia of Solid Earth Geophysics](http://link.springer.com/referenceworkentry/10.1007/978-3-030-58631-7_11)-->

---

## Forward problem in seismology

<div style="flex-basis:8%;">

A displacement seismogram $u(t)$ is the convolution of a source $s(t)$ with the medium $g(t)$ and a instrument $i(t)$, such as

$$u(t) = (s * g * i )(t),$$

where $*$ stands for convolution

$$(f * g)(t) = \int_{-\infty}^{+\infty}f(\tau)g(t - \tau)d\tau.$$

</div>
<div>

![](images/stein_and_wysession_transparent/1_1_01.png)

</div>

<!-- _footer: Figure 1.1-1 from Stein and Wysession.-->

---

## Forward problem in seismology

<div style="flex-basis:8%;">

A displacement seismogram $u(t)$ is the convolution of a source $s(t)$ with the medium $g(t)$ and a instrument $i(t)$, such as

$$u(t) = (s * g * i )(t),$$

where $*$ stands for convolution

$$(f * g)(t) = \int_{-\infty}^{+\infty}f(\tau)g(t - \tau)d\tau.$$

</div>
<div>

![width:600px](images/flowcharts/convolution.png)

</div>

<!-- _footer: Figure 1.1-1 from Stein and Wysession.-->

---

## Forward and inverse problems in seismology

<div>

### Forward problem

A displacement seismogram $u(t)$ is the convolution of a source $s(t)$ with the medium $g(t)$ and a instrument $i(t)$, such as

$$u(t) = (s * g * i )(t).$$

![height:300px](images/flowcharts/inverse.svg)

</div>
<div>

### Inverse problems

The inverse problem aims at finding the parameters that best fit the data.
<br>

- _Seismic imaging_ (tomography): find $g$ from earthquake signals $u$ assuming earthquake properties $s$ and instrument $i$.

$$g = u * (s * i)^{-1}$$

- _Earthquake characterization_ (location, mechanism, source time function): find $s$ from $u$ assuming $g$ and $i$.

$$s = u * (g * i)^{-1}$$

</div>

---

## Forward and inverse problems in seismology

<div>

### Forward problem

A displacement seismogram $u(t)$ is the convolution of a source $s(t)$ with the medium $g(t)$ and a instrument $i(t)$, such as

$$u(t) = (s * g * i )(t).$$

### Inverse problems

The inverse problem aims at finding the parameters that best fit the data.

</div>
<div>

### Inversion strategy

We try to _minimize_ the error between the __data__ and the __synthetics__:

![height:300px](images/flowcharts/inversion.svg)

We can also solve this problem by maximizing the likelihood, like in Bayesian inference.

</div>

---

## Inverse problems in seismology

<div style="flex-basis:30%">

__Seismic tomography__: find $g$ from earthquake signals $u$ assuming earthquake properties $s$ and instrument $i$.

</div>
<div>

![](images/tomography/seismic_tomography.png)

</div>

<!-- _footer: Tomographic paths from [Koulakov & Shapiro (2021)](https://link.springer.com/referenceworkentry/10.1007/978-3-642-36197-5_51-1) -->

---

## Inverse problems in seismology

<div style="flex-basis:1%">

__Seismic tomography__: find $g$ from earthquake signals $u$ assuming earthquake properties $s$ and instrument $i$.

<mark>__Depending on the method, the data, the instruments, the results may differ!__</mark>

</div>
<div>

![](images/tomography/models_tomography.png)

</div>

<!-- _footer: Comparison of global tomographic model comparison from [Durand et al. (2017)](https://academic.oup.com/gji/article/211/3/1628/4222797) -->

---

## Inverse problems in seismology

<div style="flex-basis:40%">

__Earthquake characterization__: (location, mechanism, source time function): find $s$ from $u$ assuming $g$ and $i$.

<mark>Again, __depending on the method, the data, the instruments, the results may differ!__</mark>

</div>
<div>

![](images/catalogs/catalog_compare.png)

</div>

<!-- _footer: Comparison of earthquake catalogs in Greece from [Duverger et al. (2018)](https://academic.oup.com/gji/article/215/1/196/5046732) -->

---

## Inverse problems in seismology

<div style="flex-basis: 40%">  

__Inverse problems are always inaccurate.__

- Forward problem always __approximate__ reality
- Inverse problems rarely find __exact and unique__ solutions
- Observations can only be predicted with __limited accuracy__.
- Improving accurracy involves improving models, method or instruments.

- __Altogether, this leads to new discoveries__.

</div>
<div>

![](images/seismograms/synthetic_seismogram_blom_et_al.png)

</div>

<!-- _footer: Predicted versus observed seismograms from [Blom et al. (2019)](https://se.copernicus.org/preprints/se-2019-152/se-2019-152.pdf)-->

---

<!-- _class: titlepage-->

# 3. Seismological data

![bg brightness:0.6](images/cover/ring_of_fire.jpg)

---

## Seismic instruments

<div style="flex-basis:55%;">

### (Maybe) the first seismic instrument

Seismoscope made in China in 132 by Zhang Heng. This instrument possibly led to the first unfelt ground motion detection in 143.

<br>

![width:700px drop-shadow](images/seismometers/seismograph_action.png)

</div>
<div>

![](images/seismometers/seismograph.png)

</div>

---

## Seismic instruments

<div>

### First instruments

Early seismometer constructed by Wiechert in 1904 with a pendulum weight of _17 tons_.

<br>

![width:300px drop-shadow](images/seismometers/sismometre.png)

</div>
<div>

![](images/seismometers/wiechert_seismograph_1904.png)

</div>

<!-- _footer: _Left_: Seismometer scheme from www.seis-insight.eu <br> _Right_: Seismometer original sketch from [Ritter et al. (2000)](https://www.sciencedirect.com/science/article/abs/pii/S0040195100002122?via%3Dihub) -->

---

## Seismic instruments

<div>

### First record of a distance earthquake

Historical first recording of an earthquake in Japan from a seismometer (a tiltmeter) in Potsdam (in 1895).

<br>

![width:400px](images/seismograms/first_record.png)

</div>
<div>

![width:500px; drop-shadow](images/seismograms/Nature%201889.png)

</div>

<!-- _footer: _Left and right_: image from [Nature, July 25, 1889](https://www.nature.com/articles/040294e0#citeas)-->

---

## Seismic instruments

<div style="flex-basis: 30%">

![height:500px](images/stein_and_wysession_transparent/6_6_05.png)

</div>
<div>

### Modern digital seismometers

Nowadays, we mesure the feedback force voltage $V(t)$ required to keep the mass in place. The data is then turn into digital data and stored in datacenters.

![width:600px](images/flowcharts/seismograph_circuit.svg)

</div>

<!-- _footer: Figure 6.6-5 from Stein and Wysession-->

---

## Seismic instruments

![bg](images/seismometers/station_NEA_site.jpg)

<!-- _footer: <span style="color:black; background-color: rgba(255, 255, 255, 0.5); padding:4px;"> Example seismic station site in Northern California, with the solar panels on the right and the vault entrance on the left. <br> © _David Croker_, United States Geological Survey.</span>-->

---

## Digital seismic data

<div style="flex-basis:40%">

The continuous displacement $u(t)$ is measured every sampling period $\Delta t$, such as

$$u[n] = u(n \Delta t)$$

where $n$ is the sample index $n = 1 \ldots N$. __We call $u[n]$ the digital signal__, with the following quantities:

- Number of samples: $N$
- Sampling rate: $f_s = N / T = 1 / \Delta t$ in Hz (or s.p.s. for samples per seconds)
- Duration = $N\Delta t = N/f_s$ in seconds.

</div>
<div>

![width:600px](images/seismograms/digital.svg)

</div>

---

## Digital seismic networks

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

## Digital seismic networks

<div>

### Reference website for all networks

- All networks designated with a unique Digital Object Identifier (DOI)
- Map, coordinates, and operating range of all networks
- Other metadata (sensitivity, operating range, sensor type)

</div>
<div>

![](images/references/homepage_fdsn.png)

</div>

<!-- _footer: Online access at www.fdsn.org-->

---

<div>

## Digital seismic data

### Reference website for downloading data

- The IRIS (Incorporated Research Institutions for Seismology) dataceter is one of the biggest datacenter in seismology.
- Tools for seeking waveform availability against time and space
- Other datacenters (among others) include: IPGP, RESIF, SCEDC, NCEDC, GFZ

</div>
<div>

![](images/references/homepage_iris.png)

</div>

<!-- _footer: Online access at ds.iris.edu-->

---

<!-- _class: titlepage-->

# 4. Cutting-edge seismology

![bg brightness:0.5](images/cover/das.jpg)

---

## Cutting-edge seismology

- Build new type of seismic sensors
- Build modern seismic networks
- Build data acquisition, archiving and distribution systems
- __Explore the big data__:
  - Digital signal analysis
  - Big data challenges
  - Open access challenges
  - High-performance computing

![bg right cover](images/cover/das.jpg)

---

## Cutting-edge instruments

<div style="text-align: center">

### Towards future seismic networks with distributed acoustic sensing

![width:960px](images/seismometers/das.png)

</div>

<!-- _footer: Fiber optics scheme modified from [motionsignaltechnologies.com](https://motionsignaltechnologies.com/what-is-das-and-what-is-it-measuring/)-->

---

## Cutting-edge methods

<div style="text-align: center">

### Using seismic noise to infer seismic images

![width:900px](images/seismograms/crosscorrelation.png)

</div>

<!-- _footer: Ambient noise virtual source from [Weaver (2005)](https://www.science.org/doi/10.1126/science.1109834)-->

---

## Cutting-edge methods

<div style="text-align: center">

### Deep-learning-based earthquake signal denoising

![width:1000px](images/seismograms/denoising.png)

</div>

<!-- _footer: From [Zhu et al. (2018)](https://www.semanticscholar.org/reader/248caa1cd0f2bbb6b5b4fc2550e85b622ce51ef7)-->

---

## Cutting-edge methods

<div style="text-align: center">

### Deep-learning-based earthquake signal detection

![width:900px](images/seismograms/detector.png)

</div>

<!-- _footer: From [Ross et al. (2018)](https://www.semanticscholar.org/paper/Generalized-Seismic-Phase-Detection-with-Deep-Ross-Meier/e178d94a0601f0f395cf6d81b884a238331fa869)-->

---

## Cutting-edge methods

<div style="text-align: center">

### Machine-learning-based pattern recognition in continuous data

![width:1000px](images/seismograms/clustering.png)

</div>

<!-- _footer: From [Beyreuther et al. (2012)](https://www.semanticscholar.org/paper/Generalized-Seismic-Phase-Detection-with-Deep-Ross-Meier/e178d94a0601f0f395cf6d81b884a238331fa869)-->

---

<!-- _class: titlepage-->

<div style="text-align:center;">

# 5. Notebook time

Now it is time to switch on your machine, and launch Jupyter. We will first work on `lecture_1_observational_seismology.ipynb`

</div>

![bg brightness:0.3](images/cover/tonga_eruption.png)

---

## Just before we go

<div>

### About questions during coding sessions

- You can ask me all the questions you want
- __But__, for the sake of the evaluation, I will not all your questions
- Those questions I did not answer, I will do afterwards

</div>

---

## The Zen of Python

<div style="font-size:70%; margin-top: 15px;">

Beautiful is better than ugly.
Explicit is better than implicit.
Simple is better than complex.
Complex is better than complicated.
Flat is better than nested.
Sparse is better than dense.
Readability counts.
Special cases aren't special enough to break the rules.
Although practicality beats purity.
Errors should never pass silently.
Unless explicitly silenced.
In the face of ambiguity, refuse the temptation to guess.
There should be one-- and preferably only one --obvious way to do it.
Although that way may not be obvious at first unless you're Dutch.
Now is better than never.
Although never is often better than _right_ now.
If the implementation is hard to explain, it's a bad idea.
If the implementation is easy to explain, it may be a good idea.
Namespaces are one honking great idea -- let's do more of those!

</div>

---

## Coding with style is not an option

- Please, read [The PEP8 guidelines](https://peps.python.org/pep-0008/) before naming as  `f` you variable for the frequency, and `t` for the time.

- Also, consider [transforming code into beautiful, idiomatic Python](https://gist.github.com/0x4D31/f0b633548d8e0cfb66ee3bea6a0deff9).

- Although, I won't base the evaluation on the code's readability.

---

<!-- _class: titlepage-->

![bg brightness:0.3](images/cover/tonga_eruption.png)

<div style="position:absolute; left: 200px; right: 200px; font-size: 15pt; color:white; text-align: left">

__© 2022 Leonard Seydoux.__ This class was created for the class _Scientific Computing for Geophysical Problem_  driven by Alexandre Fournier, professor at the Institut de Physique du Globe de Paris, France. It is intended to be freely used for personal purposes only. This slide set was created with the Marpit framework, with a custom CSS theme. The present slides are avaialable on GitHub in my respositories (user leonard-seydoux) in both HTML and PDF formats. The exercices were originally created for the present class. Other tutorials on seismology are available online, and I do recommend to test them.

</div>
