#import "@preview/hydra:0.6.2": hydra

#set page(
  paper: "a4",
  margin: (x: 2cm, y: auto),
  numbering: "-1-",
  header: context{
    if calc.odd(here().page()) {
    align(right, emph(hydra(1)))
  } else {
    align(left, emph(hydra(2)))
  }
  line(length: 100%)
  }
)
#set text(font: "New Computer Modern", size: 12pt)
#set par(justify: true)
#set heading(numbering: "1.")

#show heading.where(level: 1): it => pagebreak(weak: true) + it

#page(numbering: none)[
  #set align(center)
  
  // Logo ou Nom de l'Institution
  #table(
    stroke:none,
    columns: (1fr, 1fr),
    align: (left+horizon, right+horizon),
    image("logocomp.png", width: 60%),
    image("logo.png", width: 74%)
  )
  #v(1cm)
  #text(size: 18pt, weight: "bold")[MASTER THESIS]
  #line(length: 70%, stroke: 0.5pt)
  #v(1cm)
  
  // Titre du Document
  #text(size: 26pt, weight: "bold")[Evolution of Quantum Open System under Continuous Measurement]
  
  #v(0.5cm)
  #text(size: 14pt, style: "italic", fill: gray.darken(20%))[Keywords: \ \ Stochastic Master Equation \ Quantum Trajectories \ Photodetection \ Homodyne detection]
  
  #v(5fr)

  // Informations sur l'auteur et l'encadrant
  #grid(
    columns: (1fr, 1fr),
    align(left)[
      #text(weight: "bold", size: 15pt)[STUDENT:] \
      Jules STEVENOT
    ],
    align(right)[
      #text(weight: "bold", size: 15pt)[SUPERVISOR:] \
      Bruno BELLOMO
    ]
  )
  
  // Date et lieu
  #v(1fr)
  #align(center, text(size: 13pt, [January -- June \ 2026]))
]

// On réinitialise le compteur de pages pour que le contenu commence à 1
#counter(page).update(0)

#pagebreak()

#outline()
#counter(page).update(0)

#pagebreak()

= Introduction

