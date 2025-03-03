import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/MurrellGroup.github.io/MolecularEvolution.jl/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Core Concepts', collapsed: false, items: [
{ text: 'Intro', link: '/intro' },
{ text: 'The MolecularEvolution.jl Framework', link: '/framework' },
{ text: 'Models', link: '/models' }]
 },
{ text: 'Methods & Algorithms', collapsed: false, items: [
{ text: 'Simulation', link: '/simulation' },
{ text: 'Optimization', link: '/optimization' },
{ text: 'Ancestral Reconstruction', link: '/ancestors' }]
 },
{ text: 'Extensions & Utilities', collapsed: false, items: [
{ text: 'Input/Output', link: '/IO' },
{ text: 'Visualization', link: '/generated/viz' },
{ text: 'Updating a phylogenetic tree', link: '/generated/update' }]
 },
{ text: 'Examples', link: '/examples' },
{ text: 'Full API', link: '/api' }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/MurrellGroup.github.io/MolecularEvolution.jl/',// TODO: replace this in makedocs!
  title: 'MolecularEvolution.jl',
  description: 'Documentation for MolecularEvolution.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../final_site', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: '/favicon.ico' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  ignoreDeadLinks: true,
  vite: {
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.png', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Core Concepts', collapsed: false, items: [
{ text: 'Intro', link: '/intro' },
{ text: 'The MolecularEvolution.jl Framework', link: '/framework' },
{ text: 'Models', link: '/models' }]
 },
{ text: 'Methods & Algorithms', collapsed: false, items: [
{ text: 'Simulation', link: '/simulation' },
{ text: 'Optimization', link: '/optimization' },
{ text: 'Ancestral Reconstruction', link: '/ancestors' }]
 },
{ text: 'Extensions & Utilities', collapsed: false, items: [
{ text: 'Input/Output', link: '/IO' },
{ text: 'Visualization', link: '/generated/viz' },
{ text: 'Updating a phylogenetic tree', link: '/generated/update' }]
 },
{ text: 'Examples', link: '/examples' },
{ text: 'Full API', link: '/api' }
]
,
    editLink: { pattern: "https://https://github.com/MurrellGroup/MolecularEvolution.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/MurrellGroup/MolecularEvolution.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
