# glycanPredict

## Overview
R package for calculation of possible glycan compositions and annotation in mass spectrometry data.

This package is intended for “calculation” of all possible glycan compositions within a set of constraining parameters. It dynamically builds a database based on a set of combinatorial equations. Specifically, the user indicates which monomer types (hexose only or hexose and pentose), degree of polymerisation (length) range and modification types (e.g. sulphate, carboxylic acid, deoxy - full list of options is described below) should be included, the desired maximum for the average number of modifications per monomer and whether mono-/oligosaccharides are labelled or not (see label options below). There is also the option for whether or not double sulphation of a single monomer is allowed. The function “builds” names (any of IUPAC, GlycoCT or Oxford), formulas and masses for all glycan compositions possible within the constraining parameters. The user also provides  parameters related to mass spectrometry: ionisation type, polarity and scan range. Formulae and *m/z* values of ions are calculated depending on the ionisation mode and modifications. Compositions which give no ions with *m/z* values within the given scan range are removed. The final output is returned as a (wide or long format) dataframe. 
