import{_ as o,C as r,c as p,o as d,j as e,a,aA as n,G as i,w as l}from"./chunks/framework.ByC1u5oh.js";const C=JSON.parse('{"title":"Input/Output","description":"","frontmatter":{},"headers":[],"relativePath":"IO.md","filePath":"IO.md","lastUpdated":null}'),u={name:"IO.md"},h={class:"jldocstring custom-block",open:""},k={class:"jldocstring custom-block",open:""},c={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""},f={class:"jldocstring custom-block",open:""},b={class:"jldocstring custom-block",open:""};function m(E,s,_,y,v,F){const t=r("Badge");return d(),p("div",null,[s[25]||(s[25]=e("h1",{id:"input-output",tabindex:"-1"},[a("Input/Output "),e("a",{class:"header-anchor",href:"#input-output","aria-label":'Permalink to "Input/Output"'},"​")],-1)),e("details",h,[e("summary",null,[s[0]||(s[0]=e("a",{id:"MolecularEvolution.write_nexus",href:"#MolecularEvolution.write_nexus"},[e("span",{class:"jlbinding"},"MolecularEvolution.write_nexus")],-1)),s[1]||(s[1]=a()),i(t,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[3]||(s[3]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">write_nexus</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fname</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">String</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,tree</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">FelNode</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Writes the tree as a nexus file, suitable for opening in eg. FigTree. Data in the <code>node_data</code> dictionary will be converted into annotations. Only tested for simple <code>node_data</code> formats and types.</p>',2)),i(t,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[2]||(s[2]=[e("a",{href:"https://github.com/MurrellGroup/MolecularEvolution.jl/blob/2b710268f705c279490460497de47001e6d12bc2/src/utils/misc.jl#L281-L287",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),e("details",k,[e("summary",null,[s[4]||(s[4]=e("a",{id:"MolecularEvolution.newick",href:"#MolecularEvolution.newick"},[e("span",{class:"jlbinding"},"MolecularEvolution.newick")],-1)),s[5]||(s[5]=a()),i(t,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[7]||(s[7]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">newick</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(root)</span></span></code></pre></div><p>Returns a newick string representation of the tree.</p>',2)),i(t,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[6]||(s[6]=[e("a",{href:"https://github.com/MurrellGroup/MolecularEvolution.jl/blob/2b710268f705c279490460497de47001e6d12bc2/src/core/nodes/AbstractTreeNode.jl#L626-L630",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),e("details",c,[e("summary",null,[s[8]||(s[8]=e("a",{id:"MolecularEvolution.read_newick_tree",href:"#MolecularEvolution.read_newick_tree"},[e("span",{class:"jlbinding"},"MolecularEvolution.read_newick_tree")],-1)),s[9]||(s[9]=a()),i(t,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[11]||(s[11]=e("p",null,"read_newick_tree(treefile)",-1)),s[12]||(s[12]=e("p",null,"Reads in a tree from a file, of type FelNode",-1)),i(t,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[10]||(s[10]=[e("a",{href:"https://github.com/MurrellGroup/MolecularEvolution.jl/blob/2b710268f705c279490460497de47001e6d12bc2/src/utils/misc.jl#L256-L260",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),e("details",g,[e("summary",null,[s[13]||(s[13]=e("a",{id:"MolecularEvolution.populate_tree!",href:"#MolecularEvolution.populate_tree!"},[e("span",{class:"jlbinding"},"MolecularEvolution.populate_tree!")],-1)),s[14]||(s[14]=a()),i(t,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[16]||(s[16]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">populate_tree!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(tree</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">FelNode</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, starting_message, names, data; init_all_messages </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, tolerate_missing </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, leaf_name_transform </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x)</span></span></code></pre></div><p>Takes a tree, and a <code>starting_message</code> (which will serve as the memory template for populating messages all over the tree). <code>starting_message</code> can be a message (ie. a vector of Partitions), but will also work with a single Partition (although the tree) will still be populated with a length-1 vector of Partitions. Further, as long as <code>obs2partition</code> is implemented for your Partition type, the leaf nodes will be populated with the data from <code>data</code>, matching the names on each leaf. When a leaf on the tree has a name that doesn&#39;t match anything in <code>names</code>, then if</p><ul><li><p><code>tolerate_missing = 0</code>, an error will be thrown</p></li><li><p><code>tolerate_missing = 1</code>, a warning will be thrown, and the message will be set to the uninformative message (requires identity!(::Partition) to be defined)</p></li><li><p><code>tolerate_missing = 2</code>, the message will be set to the uninformative message, without warnings (requires identity!(::Partition) to be defined)</p></li></ul><p>A renaming function that can eg. strip tags from the tree when matching leaf names with <code>names</code> can be passed to <code>leaf_name_transform</code></p>',4)),i(t,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[15]||(s[15]=[e("a",{href:"https://github.com/MurrellGroup/MolecularEvolution.jl/blob/2b710268f705c279490460497de47001e6d12bc2/src/utils/misc.jl#L110-L122",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),e("details",f,[e("summary",null,[s[17]||(s[17]=e("a",{id:"MolecularEvolution.read_fasta",href:"#MolecularEvolution.read_fasta"},[e("span",{class:"jlbinding"},"MolecularEvolution.read_fasta")],-1)),s[18]||(s[18]=a()),i(t,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[20]||(s[20]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">read_fasta</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(filepath</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">String</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Reads in a fasta file and returns a tuple of (seqnames, seqs).</p>',2)),i(t,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[19]||(s[19]=[e("a",{href:"https://github.com/MurrellGroup/MolecularEvolution.jl/blob/2b710268f705c279490460497de47001e6d12bc2/src/utils/fasta_io.jl#L4-L8",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),e("details",b,[e("summary",null,[s[21]||(s[21]=e("a",{id:"MolecularEvolution.write_fasta",href:"#MolecularEvolution.write_fasta"},[e("span",{class:"jlbinding"},"MolecularEvolution.write_fasta")],-1)),s[22]||(s[22]=a()),i(t,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[24]||(s[24]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">write_fasta</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(filepath</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">String</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, sequences</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Vector{String}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; seq_names </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> nothing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Writes a fasta file from a vector of sequences, with optional seq_names.</p>',2)),i(t,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[23]||(s[23]=[e("a",{href:"https://github.com/MurrellGroup/MolecularEvolution.jl/blob/2b710268f705c279490460497de47001e6d12bc2/src/utils/fasta_io.jl#L18-L22",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})])])}const w=o(u,[["render",m]]);export{C as __pageData,w as default};
