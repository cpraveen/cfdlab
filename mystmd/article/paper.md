---
title: This is the main title
subject: Subject
subtitle: This is the subtitle
short_title: This is the short title
authors:
  - name: Praveen Chandrashekar
    affiliations:
      - Centre for Applicable Mathematics
      - Tata Institute of Fundamental Research
      - Bangalore - 560065, India
    orcid: 0000-0003-1903-4107
    email: praveen@mail.com
license: CC-BY-4.0
keywords: myst, markdown, open-science
numbering:
  code: false
  equation: true
  title: false
  headings: true
exports:
  - format: pdf
    template: arxiv_nips
    output: paper.pdf
  - format: tex
    template: arxiv_nips
math:
  '\re': '\mathbb{R}'
  '\complex': '\mathbb{C}'
  '\float': '\mathbb{F}'
  '\poly': '\mathbb{P}'
  '\half': '\frac{1}{2}'
  '\shalf': '\tfrac{1}{2}'
  '\ud': '\textrm{d}'
  '\od': '\frac{\ud #1}{\ud #2}'
  '\pd': '\frac{\partial #1}{\partial #2}'
  '\cts': '\mathcal{C}'
  '\ii': '\mathfrak{i}'
  '\imh': '{i-\tfrac{1}{2}}'
  '\norm': '\|#1\|'
  '\ip': '\left(#1\right)'
  '\ee': '\textrm{e}'
  '\limplies': '\qquad \Longrightarrow \qquad'
  '\clr': '{\color{#1}#2}'
abstract: |
  This is the abstract.
---

# Introduction

Single line equation

$$
\label{eq:ode}
\od{y}{t} = f(y,t)
$$

A PDE

$$
\pd{u}{t} + \pd{f}{x} = 0
$$

Navier-Stokes equations

$$
v_t + \clr{red}{v \cdot \nabla v} + \nabla p = \clr{blue}{\nu \Delta v}
$$

Multiline equation using align

\begin{align}
y &= a_1 x_1 + a_2 x_2 \\
  &= b_1 z_1 + b_2 z_2
\end{align}

Using gather

\begin{gather*}
a_1=b_1+c_1\\
a_2=b_2+c_2-d_2+e_2
\end{gather*}

We see from [](#eq:ode) that we have an ODE.

# This is another section

See this for some tips: https://curvenote.com/docs/publish/authoring-in-myst

Definition
: This is a definition.

## This is a sub-section

This is itemize.

* first item
* second item

This is enumerate.

1. first item
1. second item

Create a note, see more here https://mystmd.org/guide/admonitions

:::{note} Watch out
This is a note.
:::

This is an exercise, it does not work in latex/pdf. https://mystmd.org/guide/exercises

:::{exercise} Exercise
This is an exercise.
:::

Sphinx proof provides may environments, see https://mystmd.org/guide/proofs-and-theorems

:::{prf:remark} Remark
This is a remark.
:::

:::{prf:theorem} Mean-value theorem
:label: thm:mvt
This is a theorem.
:::

See [](#thm:mvt).

:::{prf:proof}
This is the proof.
:::

:::{prf:lemma} Mean-value theorem
This is a lemma.
:::

:::{prf:definition} Definition
This is a definition.
:::

Put a figure

:::{figure} https://www.math.tifrbng.res.in/logo.png
:label: fig:sunset
:align: center
:width: 50%
This is is the figure caption.
:::

For tabular figures, see https://sphinx-subfigure.readthedocs.io

Adding table is like this https://mystmd.org/guide/tables

:::{table} Table caption
:label: table
:align: center

| left | center | right |
| :--- | :---:  | ---:  |
| baz  | bim    | xyz   |
| baz  | bim    | abc   |

:::
