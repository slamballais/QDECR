---
title: About
description: 'Who makes QDECR, where it came from, and the conference material from 2019 and 2020.'
---

QDECR is free software for vertex-wise analysis of FreeSurfer data in R, released under the [GPL-3 licence](https://www.gnu.org/licenses/gpl-3.0.html). Its code is on [GitHub](https://github.com/slamballais/QDECR).

## The team

### Sander Lamballais

Wrote QDECR's code, made it fast, and maintains it. Sander is a bioinformatician at the Erasmus MC in Rotterdam, the Netherlands, building and validating the data pipelines behind long-read rare-disease diagnostics. Before that came a PhD from Erasmus University Rotterdam in 2022, with the thesis [Shaping the brain: Causes and consequences of the changing brain across the lifespan](https://pure.eur.nl/en/publications/shaping-the-brain-causes-and-consequences-of-the-changing-brain-a/), and postdoctoral research in neurogenetics. [ORCID 0000-0003-3118-6330](https://orcid.org/0000-0003-3118-6330).

### Ryan Muetzel

Conceived QDECR. Ryan is an assistant professor and principal investigator at the [Department of Child and Adolescent Psychiatry/Psychology](https://www.erasmusmc.nl/en/research/researchers/muetzel-ryan) of the Erasmus MC in Rotterdam, leading the Integrative and Precision Neuroimaging group and co-directing Generation R Neuroimaging. Ryan's research uses longitudinal population neuroimaging to study typical and atypical brain development, after a PhD from Erasmus University Rotterdam in 2016 with the thesis [The Connections Within: Pediatric population-based neuroimaging of brain development](https://repub.eur.nl/pub/94642). [ORCID 0000-0003-3215-1287](https://orcid.org/0000-0003-3215-1287).

The conference work below was written with Henning Tiemeier, Meike W. Vernooij and M. Arfan Ikram.

## History

QDECR began at the Erasmus MC as a side project: a way to make FreeSurfer's group analysis, the QDEC tool and `mri_glmfit`, easier to use for colleagues who worked in R rather than in Bash. The name is QDEC with an R. Researchers in epidemiology wanted what their other analyses already had: adjustment for many covariates, imputed data for missing values, and samples of thousands. So the analysis was rewritten in R, around R's own formulas.

- **2019.** Version 0.7.0, the first public release, was presented at OHBM 2019 in Rome with a poster and a software demonstration.
- **2020.** At OHBM 2020, held online, a prototype for vertex-wise mixed models. It was never released; [verywise](https://github.com/SereDef/verywise) now fits those models.
- **2021.** The paper describing QDECR appeared in *Frontiers in Neuroinformatics*. [How to cite it](/cite).
- **2022.** Version 0.9.0, "Lausanne", added weighted regression. It is the current release.
- **2026.** This site replaced the 2019 one, and work on 0.10.0 began.

## Archive

The material from the two OHBM meetings, kept as it was presented. It describes QDECR as it was then: some of what it says has changed since, and the tutorials describe the current release.

### OHBM 2019

- [Poster: The QDECR package, a flexible, extensible vertex-wise analysis framework in R](/archive/qdecr_ohbm2019_poster.pdf) (PDF, 1.4 MB).
- [Brochure: A quick guide to QDECR](/archive/qdecr_ohbm_brochure.pdf), handed out at the software demonstration (PDF, 8 pages, 2.0 MB).

### OHBM 2020

- [Poster: Vertex-wise mixed modeling using QDECR](/archive/qdecr_lmm.pdf) (PDF, 1.3 MB). The `qdecr_lmm` function it shows was a prototype and is not part of QDECR.
- [Video: the poster presented](https://www.youtube.com/watch?v=utqh73BSTqA), on YouTube.

The questions answered at the 2020 meeting are now part of [Troubleshooting and FAQ](/tutorials/troubleshooting#frequently-asked-questions).

### Code Ocean

- [A Code Ocean capsule of QDECR](https://doi.org/10.24433/CO.2177760.v1) from 2020, which the 2020 poster cites.
