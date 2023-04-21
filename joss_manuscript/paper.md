---
title: 'MolPainter: A Tool for Painting and Solvating Layered Molecular Systems'
tags:
  - python
  - molecularDynamics
  - lipidBilayers
authors:
  - name: George A. Pantelopulos
    orcid: 0000-0002-4373-1677
    affiliation: "1,2"
    corresponding: true
  - name: Aaron Liberatore
    affiliation: 3
affiliations:
 - name: Laboratory of Chemical Physics, National Institute of Diabetes and Digestive and Kidney Diseases, National Institutes of Health, Bethesda, Maryland, USA
   index: 1
 - name: Department of Chemistry, Boston University, Boston, Massachusetts, USA
   index: 2
 - name: Independent Researcher, USA
   index: 3
date: 25 December 2022
bibliography: paper.bib

---

# Summary

MolPainter presents the metaphor of painting as an abstraction for the construction of complex layered molecular systems. In particular, systems whose equilibrium states exhibit distinctive visual artifacts can be intentionally crafted through interactive, visually-oriented molecule placement. This is in contrast to non-interactive methods, which allow neither manual placement of molecules nor visual confirmation of desired states. The principal use-case is presented for modelling lipid bilayers with molecular dynamics simulation, though MolPainter may be used to construct any layered molecular system provided appropriate PDB-format files to define molecules (\autoref{fig:Figure1}).

![MolPainter uses xy-plane grids at defined z-positions and allows the user to paint user-supplied PDB files into these grids to define a molecular system. MolSolvator solvates these systems after construction.\label{fig:Figure1}](figures/MolPainterGraphic.png)

# Statement of Need

Complexity in systems investigated using molecular dynamics (MD) and similar simulation methods employing MD force fields is increasing in terms of composition, initial condition, and system size. Due to limitations in time scales accessible even in highly coarse-grained force fields, careful selection of molecular composition and initial positions in complex systems is essential to designing useful simulations. With current tools for preparing molecular systems, it is straightforward to construct systems of complex composition but with random placement of molecules in the system. From such random initial distributions, systems would typically need very long, potentially intractable time scales to reach equilibrium. Construction of complex systems of realistic compositions, such as the *Mycoplasma genitalium* cell [@Yu2016], are impressive, though it is difficult to construct any non-random initial positions for such systems.

# Software with similar functionalities

Several useful tools have been made available for the construction of complex molecular systems. Tools for off-lattice random insertion of molecules such as PACKMOL [@Martinez2009] and the insert-molecules program in GROMACS [@VanDerSpoel2005], enable construction of off-lattice molecular systems by placement of randomly translated and rotated molecules into user-defined volumes with optional restraints on specific atom positions. Tools for phase-specific preparation of simulations such as CHARMM-GUI [@Jo2008], PACKMOL-MemGen [@Schott-Verdugo2019], and insane.py [@Wassenaar2015] for all-atom and coarse-grained lipid bilayer simulations are popular and widely used, and handle complex lipid mixtures by allowing the user to control the number and identity of lipids of lipids randomly, laterally placed within two flat leaflets. These programs also provide support for the inclusion of solutes, such as transmembrane proteins. Similar tools for constructing arbitrary, flat bilayers such as MemGen [@Knight2015], vesicles with insane.py and CHARMM-GUI, and more complex curved geometries with BUMPy [@Boyd2018] and LipidWrapper [@Durrant2014], which also insert lipids and other molecules in a random, lateral distribution, provide for more flexibility and creativity in system construction though with less complete support in preparing the system for simulation. However, there is a lack of tools to assist in the fine-tuning of the spatial placement of molecules in lipid bilayers, and more generally into any layered geometry.


# Applications

Lipid membranes are believed to exhibit specific, functional clusters of lipids, cholesterol, and proteins, often generally referred to as lipid rafts [@Sych2021]. Lipid rafts are characterized by the liquid ordered ($\mathrm{L_o}$) lipid phase, unique from the liquid disordered phase ($\mathrm{L_d}$) of the bulk membrane. The formation of the $\mathrm{L_o}$ lipid phase and phase separation from the $\mathrm{L_d}$ phase is an example of a long time scale phenomenon which is computationally intractable in all-atom MD simulation. However, as the phase diagrams of ternary lipid mixtures and cholesterol have been characterized, it is possible to estimate *a priori* whether the lipid bilayer should be phase separated, and it is also possible to infer the relative concentration of cholesterol in $\mathrm{L_o}$ and $\mathrm{L_d}$ phases. With an appropriate tool for the task, it is possible to model phase-separated states with a near-equilibrium initial spatial placement of lipids which can then be equilibrated using conventional all-atom MD simulation.

The shape of such molecular phase separations and clusters can also be unique. The $\mathrm{L_o}$ phase at high cholesterol concentrations has been proposed to form "threads" of cholesterol, in which each cholesterol has no more and no less than two cholesterol in the first solvation shell, manifesting a maze-like pattern [@Huang1999; @Miao2002; @Pantel2018]. This has been proposed to either be due to unfavorable penetration of water introduced by loss of protective lipid headgroups or specific cholesterol-cholesterol interactions involving the $\alpha$ and $\beta$ faces of cholesterol [@Bandara2017]. The shape and size of lipid phase separations can also vary significantly, manifesting as stripes or dots in vesicles [@Pantel2017] and molecular simulations with varying system dimensions and line tensions between the $\mathrm{L_o}$ and $\mathrm{L_d}$ phase [@Kwon2022].

Along the z-axis, multi-layered lipid systems, such as the long-period-periodicity (LPP) phase of lipids [@Bouwstra2001], have been an interesting subject of lipid biophysics for their importance in mediating permeation of small molecules through skin. Large scale molecular simulations of model multilamellar skin lipids have been performed to better understand these complex environments [@MacDermaid2020]. MolPainter is uniquely well-suited to construction of systems with complex lateral distributions as well as such multi-layered lipid systems.

# Functions

![MolPainter components and functions. (A) *n* layers defined at assigned z-axial coordinates in angstrom units. (B) the grid used to define the xy-plane position of each lipid in terms of a cell spacing in angstrom units and number of cells in the x- and y-dimension. (C) The molecules are defined by PDB files supplied by the user. (D) Anchor points, defined in PDB files using B-factors of 1.00, define where the centroid of all anchor points z-axially align into the grid layer pained upon. (E) The canvas holds labels of each molecule painted onto xy-positions in each grid on each layer via various painting tools, such as the spray can. (G) Solutes can be inserted to the canvas, such as the C99 protein monomer shown here [@Pantel2022], obfuscating cells across multiple grid layers to prevent clashes with painted molecules.\label{fig:Figure2}](figures/MolPainterObjectsAndFunctions.png)

Through the use of user-provided PDB files, MolPainter can be used to paint molecules onto a canvas, a lattice in the xy-plane, aligning each molecule to a specific xy-coordinate and specific z-positions defined for any number of canvases (\autoref{fig:Figure2}). MolPainter also supports the inclusion of a "solute" to the painting, allowing for ease in constructing unique environments defined by complex mixtures of solvents. With MolPainter we include a light companion script, MolSolvator, for the addition of bulk solvent to systems painted using MolPainter.

# Acknowledgements

We acknowledge insightful conversations with John E. Straub (Boston University) regarding potential MolPainter functions and applications. GAP thanks the NIH Postdoctoral Intramural Research Training Awards program for support.

# References
