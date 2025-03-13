#import "@preview/physica:0.9.4": super-T-as-transpose
#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *
#import "utils.typ": mathformatter, mref

// Font fix because FruityFeedback is trash
#set text(font: "STIX Two Text")
#show math.equation: set text(font: "STIX Two Math")

#show: codly-init.with()
#codly(languages: codly-languages)


#set heading(numbering: "A.1.1")
#set math.equation(numbering: "(1)", supplement: [Eq.])
#show: super-T-as-transpose
#set table.header()

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: true, upright: true)

///////////////////////////////////////////////////////////////////////////////

= 02417 Times Series Analysis - Assignment 1

Authors:
- Jeppe Klitgaard `<s250250@dtu.dk>`
- Yunis Wirkus `<s250700@dtu.dk>`

Date: 2025-03-02

#pagebreak()

== Plot Data

#bibliography("report.bib")
