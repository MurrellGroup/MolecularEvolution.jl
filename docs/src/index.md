```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "MolecularEvolution.jl"
  text: "Phylogenetic Modeling Framework"
  tagline: "A Julia package for flexible development of phylogenetic models"
  image:
    src: '/logo.png'
  actions:
    - theme: brand
      text: Intro
      link: /intro
    - theme: alt
      text: API reference
      link: /api
    - theme: alt
      text: View on Github
      link: https://github.com/MurrellGroup/MolecularEvolution.jl
features:
  - title: Flexible Model Development
    details: Create custom evolutionary models by defining how your new (or already existing) probability distributions evolve along a branch, and get likelihood calculations (and much more) for free! • Mix and match different models on the same phylogeny.
    link: /framework
  - title: Rich Model Library
    details: Codon, AA, nucleotide and generic discrete character models as well as continuous Brownian motion in 1D. Site- and branch-wise mixture models and more.
    link: /models
  - title: Powerful Tree Operations
    details: Optimize model parameters, branch lengths and root state/location. Perform tree rearrangements (NNI). The above can also be made in the Bayesian framework by sampling over posterior tree and parameter spaces with MCMC. • Infer ancestral states.
    link: /optimization
  - title: Simulation Tools
    details: Generate phylogenies using e.g. coalescent and logistic growth processes, and simulate both sequence evolution and continuous trait evolution on those trees.
    link: /simulation
---
```