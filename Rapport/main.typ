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
#set heading(numbering: "I.1")
#show heading.where(level: 1): it => pagebreak(weak: true) + it

// Page de Garde
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
  #v(10fr)
  
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
  
  #table(
    stroke:none,
    columns: (1fr, 1fr,1fr),
    align: (left+horizon, center+bottom, right+horizon),
    image("logoPTH.png", width: 50%),
    align(center, text(size: 13pt, [January -- June \ 2026])) ,
    image("logoUTI.jpeg", width: 70%)
  )
]

// Remerciements
#page(numbering: none)[
  #set align(center)
  #v(0.2cm)
  #text(size: 18pt, weight: "bold")[ACKNOWLEDGEMENTS]
  #v(0.2cm)

  #text(size: 12pt)[
    #h(20pt)I would like to express my deepest gratitude to my supervisor, Bruno Bellomo, for his guidance and support throughout this research. He is a man who always took the time to explore ideas, ensure our understanding, and seek out details. His expertise in Quantum Mechanics has been fundamental to the successful completion of this Master’s thesis. I am grateful to conclude my short career as a physicist with such a demanding and rewarding journey, thanks to Bruno.

    #h(20pt)I would like to thank the members of the jury for their time and interest in this work, as well as for the opportunity to present my work in front of such esteemed peers.

    #h(20pt)I would also like to thank my colleagues for their insightful discussions and moral support during this adventure: Silvio Da Silva, Jana El Badawi, Anjana Mohan, and Léo Jacquerot-Legros, for welcoming me into your office. Thank you Léo, in particular, for your time and for the discussions we had, which greatly helped me.  

    #h(20pt)Finally, I would like to thank my friends and family for their never-ending support and love: Erwan LeDoeuff, Titouan Alphonse, and Nabyan Arkan for these two years of Master’s studies and friendship together; my parents Catherine and Lionel for my entire life and education; my brother Pierre for his righteousness and compassion; and finally my other half, without whom I would not be who I am today, Alexiane.
  ]
  #line(length: 60%, stroke: 0.5pt)
  #text(size: 12pt)[
    #h(20pt)Je souhaite exprimer ma plus profonde gratitude envers mon tuteur, Bruno Bellomo, pour son accompagnement et son soutien tout au long de ce travail. Il a toujours pris le temps d’explorer les idées, de s’assurer de notre compréhension et d’en examiner les moindres détails. Son expertise en mécanique quantique a été déterminante pour la réussite de ce mémoire de Master. Je suis reconnaissant de conclure ma courte carrière de physicien par un parcours aussi exigeant qu’enrichissant, grâce à Bruno.

    #h(20pt)Je tiens également à remercier les membres du jury pour le temps et l’intérêt qu’ils ont accordés à ce travail, ainsi que pour l’opportunité de présenter mes recherches devant des pairs aussi estimés.

    #h(20pt)Je souhaite aussi remercier mes collègues pour leurs discussions enrichissantes et leur soutien moral tout au long de cette aventure : Silvio Da Silva, Jana El Badawi, Anjana Mohan et Léo Jacquerot-Legros, pour m’avoir accueilli dans votre bureau. Merci en particulier à Léo pour ton temps et pour les discussions que nous avons eues, qui m’ont beaucoup apporté.

    #h(20pt)Enfin, je voudrais remercier mes amis et ma famille pour leur soutien et leur affection sans faille : Erwan LeDoeuff, Titouan Alphonse et Nabyan Arkan pour ces deux années de Master et d’amitié partagées ; mes parents Catherine et Lionel pour toute ma vie et mon éducation ; mon frère Pierre pour sa droiture et sa compassion ; et enfin mon autre moitié, sans qui je ne serais pas celui que je suis aujourd’hui, Alexiane.
  ]

]

// Table des matières
#set page(numbering: none)
#outline()

// Introduction
#counter(page).update(0)
#set page(numbering: "-1-")
= Motivations
#v(1cm)
#text(size: 14pt)[
  Salut
]
