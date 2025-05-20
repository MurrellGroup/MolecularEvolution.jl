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
    details: Create custom evolutionary models by defining how a process evolves along a branch and get likelihood calculations (and much more) for free! • Mix and match different models on the same phylogeny.
    link: /framework
  - title: Model Library
    details: Nucleotide, AA, Codon, and generic discrete character models • Continuous models (eg. Brownian motion) • Site- and branch-wise mixture models and more.
    link: /models
  - title: Tree Operations
    details: Optimize model parameters, branch lengths, tree topology (NNI), and root position for maximum likelihood inference • Sample from the model and tree posterior with MCMC for Bayesian inference • Infer ancestral states.
    link: /optimization
  - title: Simulation Tools
    details: Sample from a range of realistic phylogenies using a flexible coalescent process • Simulate discrete and continuous data under any model over a simulated or imported phylogeny.
    link: /simulation
---
```
