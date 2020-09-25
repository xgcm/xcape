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

# Summary


# Statement of need 

`xcape` is an open-source Python package for the calculation of parameters 
commonly used in meteorology: Convective Available Potential Energy (CAPE) and 
Storm Relative Helicity (SRH). 

CAPE is a measure of the potential for instability of an atmospheric parcel [ref]. It corresponds to the integrated virtual temperature difference between a theoretical parcel and the environmental profile while the parcel is warmer than the environment, while CIN corresponds to the integration of the negative area. 

\begin{equation}

 \text{CAPE} = g \int_{LFC}^{EL} \frac{(\Theta_v_{parcel} - \Theta_v_{env})}{  \
              \Theta_v_{env}} d\text{dz}

\end{equation}

\begin{equation}

\text{CIN} = g \int_{SFC}^{LFC} \frac{(\Theta_v_{parcel} - \Theta_v_{env})}{  \
              \Theta_v_{env}} d\text{dz}

\end{equation}

Where SFC is surface, LFC is the Level of free convection, and EL is the equilibrium level, $Theta_v_{parcel}$ is the virtual potential temperature of the parcel, $Theta_v_{env}$ is the virtual potential temperature of the environment, both in K (?), $g$ the gravitational acceleration, and $z$ the height above ground. 

SRH is a measure of the potential for updraft rotation in supercells, and is often applied to anticipating tornadic storms [ref] (Davies-Jones 1984, Davies-Jones et al. 1990, Kerr and Darkow 1996). It is a quantity proportional to streamwise vorticity and storm-relative winds, defined as: 

\begin{equation}

  SRH=\int_{0}^{h}(V-C).\omega dz   

\end{equation}

Where $V$ is the ground relative wind vector $ms^{-1}$, $C$ is the storm motion $ms^{-1}$, $\omega$ is the horizontal vorticity vector $s^{-1}$, and the resulting parameter SRH has units of $m^{2}s^{-2}$. SRH is typically calculated by an integration of some depth depending on the application,


It allows for post-processing of 1D observed profile data to 4D gridded 
model data. By wrapping low-level language routines (Fortran) for speed, and 'dask' for rapid and efficient scaling of computation, 'xcape' allows for fast computation of such intensive calculations…..

The API for 'xcape' was designed… ryan

# State of the field 
Current alternative to 'xcape', such as metpy and sharpy, …. (if we say that our is much faster we have to include an example to test it)

Given the significant resources that are usually necessary to calculate CAPE and SRH, 'xcape' ability to calculate them in a fast and scalable way, this package will be of interest to the global convective storm and tropical meteorology communities, by forecasters, researcher, and students as well.




# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"


# Acknowledgements

Chiara Lepore and Ryan Abernathey acknowledge funding support by grant # NSF - OCE 17-40648: “Collaborative Research:  EarthCube Integration:  Pangeo:  An Open Source Big Data Climate Science Platform”.

# References
