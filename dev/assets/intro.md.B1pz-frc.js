import{_ as s,c as a,o as e,aA as n}from"./chunks/framework.Dyz83CXA.js";const t="/MurrellGroup.github.io/MolecularEvolution.jl/dev/assets/quick_example.BPXtjq6w.svg",u=JSON.parse('{"title":"Intro","description":"","frontmatter":{},"headers":[],"relativePath":"intro.md","filePath":"intro.md","lastUpdated":null}'),l={name:"intro.md"};function o(h,i,p,r,d,k){return e(),a("div",null,i[0]||(i[0]=[n(`<h1 id="intro" tabindex="-1">Intro <a class="header-anchor" href="#intro" aria-label="Permalink to &quot;Intro&quot;">​</a></h1><p>MolecularEvolution.jl exploits Julia&#39;s multiple dispatch, implementing a fully generic suite of likelihood calculations, branchlength optimization, topology optimization, and ancestral inference. Users can construct trees using already-defined data types and models. But users can define probability distributions over their own data types, and specify the behavior of these under their own model types, and can mix and match different models on the same phylogeny.</p><p>If the behavior you need is not already available in <code>MolecularEvolution.jl</code>:</p><ul><li><p>If you have a new data type:</p><ul><li><p>A <code>Partition</code> type that represents the uncertainty over your state.</p></li><li><p><code>combine!()</code> that merges evidence from two <code>Partition</code>s.</p></li></ul></li><li><p>If you have a new model:</p><ul><li><p>A <code>BranchModel</code> type that stores your model parameters.</p></li><li><p><code>forward!()</code> that evolves state distributions over branches, in the root-to-tip direction.</p></li><li><p><code>backward!()</code> that reverse-evolves state distributions over branches, in the tip-to-root direction.</p></li></ul></li></ul><p>And then sampling, likelihood calculations, branch-length optimization, ancestral reconstruction, etc should be available for your new data or model.</p><h3 id="Design-principles" tabindex="-1">Design principles <a class="header-anchor" href="#Design-principles" aria-label="Permalink to &quot;Design principles {#Design-principles}&quot;">​</a></h3><p>In order of importance, we aim for the following:</p><ul><li><p>Flexibility and generality</p><ul><li><p>Where possible, we avoid design decisions that limit the development of new models, or make it harder to develop new models.</p></li><li><p>We do not sacrifice flexibility for performance.</p></li></ul></li><li><p>Scalability</p><ul><li>Analyses implemented using <code>MolecularEvolution.jl</code> should scale to large, real-world datasets.</li></ul></li><li><p>Performance</p><ul><li>While the above take precedence over speed, it should be possible to optimize your <code>Partition</code>, <code>combine!()</code>, <code>BranchModel</code>, <code>forward!()</code> and <code>backward!()</code> functions to obtain competative runtimes.</li></ul></li></ul><h3 id="authors" tabindex="-1">Authors: <a class="header-anchor" href="#authors" aria-label="Permalink to &quot;Authors:&quot;">​</a></h3><p>Venkatesh Kumar and Ben Murrell, with additional contributions by Sanjay Mohan, Alec Pankow, Hassan Sadiq, and Kenta Sato.</p><h3 id="Quick-example:-Likelihood-calculations-under-phylogenetic-Brownian-motion:" tabindex="-1">Quick example: Likelihood calculations under phylogenetic Brownian motion: <a class="header-anchor" href="#Quick-example:-Likelihood-calculations-under-phylogenetic-Brownian-motion:" aria-label="Permalink to &quot;Quick example: Likelihood calculations under phylogenetic Brownian motion: {#Quick-example:-Likelihood-calculations-under-phylogenetic-Brownian-motion:}&quot;">​</a></h3><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> MolecularEvolution, Plots</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#First simulate a tree, using a coalescent process</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">tree </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sim_tree</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(n</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">200</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">internal_message_init!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(tree, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GaussianPartition</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#Simulate brownian motion over the tree</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">bm_model </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> BrownianMotion</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sample_down!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(tree, bm_model)</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#And plot the log likelihood as a function of the parameter value</span></span>
<span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">ll</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> log_likelihood!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(tree,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">BrownianMotion</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,x))</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.7</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.001</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.6</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,ll, xlabel </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;variance per unit time&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ylabel </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;log likelihood&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="`+t+'" alt=""></p>',13)]))}const g=s(l,[["render",o]]);export{u as __pageData,g as default};
