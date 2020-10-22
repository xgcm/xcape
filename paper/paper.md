---
title: 'xcape: : Fast convective parameters for numpy, dask, and xarray'
tags:
  - Python
  - Fortran
authors:
  - name: Chiara Lepore
    orcid: 0000-0002-1606-6982
    affiliation: 1 
  - name: Ryan Abernathey
    affiliation: "1, 2"
  - name: John T. Allen
    affiliation: 3
affiliations:
 - name: Lamonth-Doherty Earth Observatory
   index: 1
 - name: Columbia University
   index: 2
 - name: Central Michigan University
   index: 3
date: 22 September 2020
bibliography: paper.bib

---

# Summary and Statement of need 


Atmospheric scientists rely on profiles of the atmosphere from observations and modeled snapshots to assess the potential for thunderstorms and other hazards. Model (those data produced to forecast into the future) and reanalysis (fixed snapshots of the historical atmosphere) output parameters are rarely consistent between different datasets, outside of the standard parameters such as temperature (T), specific humidity (Q), pressure (P), and vector winds (U,V). These data are typically available on at least a four-dimensional grid of time, vertical level, and spatial grid, and in some applications can also involve multiple realizations, necessitating calculations in an n-dimensional framework. This leads to a significant challenge when working with these data to evaluate parameters that describe vertical profiles, compare parameters between datasets, or even calculate these parameters for a single dataset. In the field of severe convective storms, this issue arises often as a result of many of the quantities being calculated by integration, and therefore requiring significant computational resources, and requiring application of approaches that can take as long as a second per profile to calculate (Allen 2018). This has served to hamper development and intercomparison of model forecasts and historical reanalysis profiles. The high cost of calculations has also served to exclude users, a problem that has only increased in response to exponential growth in model and reanalysis dataset size. This necessitates an approach to calculate desired parameters that is a) applicable to any n-dimensional model or reanalysis dataset and b) scalable and efficient in its calculations to allow users to perform these calculations from personal machines through high-performance computing centers.  To address this challenge we have developed `xcape`, an open-source Python package for the calculation of two numerically intensive parameters commonly used in convective storms meteorology as well as other fields: Convective Available Potential Energy (CAPE) and Storm Relative Helicity (SRH). 


CAPE is a measure of the potential for instability of an atmospheric parcel [@emanuel1994atmospheric] . It corresponds to the integrated virtual temperature difference between a theoretical parcel and the environmental profile while the parcel is warmer than the environment, while CIN corresponds to the integration of the negative area. 

\begin{equation}
\text{CAPE} = g \int_{LFC}^{EL} \frac{(\Theta_{v,parcel} - \Theta_{v,env})}{  \
              \Theta_{v,env}} d\text{dz}
\end{equation}

\begin{equation}
\text{CIN} = g \int_{SFC}^{LFC} \frac{(\Theta_{v,parcel} - \Theta_{v,env})}{  \
              \Theta_{v,env}} d\text{dz}
\end{equation}

Where SFC is surface, LFC is the Level of free convection, and EL is the equilibrium level, $Theta_{v,parcel}$ is the virtual potential temperature of the parcel, $Theta_{v,env}$ is the virtual potential temperature of the environment, both in K (true? necessary to say?), $g$ the gravitational acceleration, and $z$ the height above ground. 

SRH is a measure of the potential for updraft rotation in supercells, and is often applied to anticipating tornadic storms [@markowski2011mesoscale]. It is a quantity proportional to streamwise vorticity and storm-relative winds, defined as: 

\begin{equation}
  SRH=\int_{0}^{h}(V-C)\cdot\omega dz   
\end{equation}

Where $V$ is the ground relative wind vector $ms^{-1}$, $C$ is the storm motion $ms^{-1}$, $\omega$ is the horizontal vorticity vector $s^{-1}$, and the resulting parameter SRH has units of $m^{2}s^{-2}$. SRH is typically calculated by an integration of some depth depending on the application,


It allows for post-processing of 1D observed profile data to 4D gridded 
model data. By wrapping low-level language routines (Fortran) for speed, and `dask` for rapid and efficient scaling of computation, `xcape` allows for fast computation of such intensive calculations..


# State of the field 
Current alternative to `xcape`, such as metpy and sharpy, . (if we say that our is much faster we have to include an example to test it)

Given the significant resources that are usually necessary to calculate CAPE and SRH, `xcape` ability to calculate them in a fast and scalable way, this package will be of interest to the global convective storm and tropical meteorology communities, by forecasters, researcher, and students as well.




# Citations


# Acknowledgements

Chiara Lepore and Ryan Abernathey acknowledge funding support by grant NSF OCE 1740648: “Collaborative Research:  EarthCube Integration:  Pangeo:  An Open Source Big Data Climate Science Platform”.

# References
