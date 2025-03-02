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
    details: Create custom evolutionary models with your own data types and probability distributions while leveraging the existing infrastructure for tree operations and likelihood calculations.
    link: /framework
  - title: Rich Model Library
    details: Includes standard models like JC69, K80, GTR for nucleotides, JTT and WAG for amino acids, MG94 and GY94 for codons, plus Brownian motion for continuous traits.
    link: /models
  - title: Powerful Tree Operations
    details: Optimize branch lengths, infer ancestral states, perform tree rearrangements (NNI), and sample topologies with MCMC.
    link: /optimization
  - title: Simulation Tools
    details: Generate phylogenies using both coalescent and birth-death processes, and simulate sequence evolution on those trees.
    link: /simulation
---
```

```@meta
CurrentModule = MolecularEvolution
```
