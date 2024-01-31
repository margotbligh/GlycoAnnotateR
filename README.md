# GlycoAnnotateR

## Overview
GlycoAnnotateR is an R package for data-base free annotation of glycan compositions in mass spectrometry data. The package is designed to be flexible and work with many different types of mass spectrometry data (e.g. LC-MS, MALDI, direct injection), as well as the output of many different data processing pipelines (e.g. [XMCS](https://github.com/sneumann/xcms), [Cardinal](https://github.com/kuwisdelu/Cardinal)). It currently consists of two main functionalities:

`glycoPredict`: this function uses combinatorial mathematics to 'predict' and build every combination of monomers and modifications possible based on the input parameters.

`glycoAnnotate`: this function annotates mass spectrometry peaks or features by comparison to the predicted, theoretical compositions with the specified error allowed.

## Installation

This package can be installed directly from Github using devtools:

```
library(devtools)
devtools::install_github('margotbligh/GlycoAnnotateR')
```

**Please note that python is required for the package to function**. If you do not have a local version of python available, please follow [instructions](https://wiki.python.org/moin/BeginnersGuide/Download) to download and install.

## Prediction parameters

The 'prediction' or 'calculation' of glycan compositions is the core utility of this package. Therefore a detailed description of the arguments is provided here.

### Glycan composition parameters

* Degree of polymerisation, `dp`

  This is **always** a range from the lowest to highest DP desired (e.g. `c(1,10)` for DPs from 1 to 10)- if you need only a single DP, provide that DP twice (e.g. `c(2,2)` for only DP 2).

* Should pentose be included in addition to hexose, `pent_option`

  This is a logical argument for whether pentose monomers should be included in compositions in addition to hexose monomers.

* Maximum number of modifications per monomer on average, `nmod_max`

  Calculated by the number of modifications over the number of monomers. Does not take into account unsaturated, alditol or dehydrated. For example, for a tetramer (DP4) of deoxyhexoses with four sulphate groups, the average number of modifications per monomer is 2. By default `nmod_max` is 1, and the maximum allowed value is 3. Consider carefully whether you need to increase this value above the default.

* Label, `label`

  Are sugars labelled by reductive amination? Current supported labels are: "none", "procainamide","2-aminobenzoic acid", "2-aminobenzamide", "1-phenyl-3-methyl-5-pyrazolone". Common abbreviations or notations for these labels are generally accepted (e.g. 'pmp' or 'PMP' for the latter).

* Modifications, `modifications`

  By default, each modification can occur once per monomer, and it is possible to have all modifications selected present on one monomer. After calculation of modified monomers they are filtered by the `nmod_max` term before output is returned. So, for example, for `modifications = c('deoxy', 'sulphate', 'carboxylicacid')`, the program will generate as one possible composition all three modifications on one monomer (i.e. 'DeoxyHex1 CarboxylicAcid1 Sulphate1'). If `nmod_max` is at the default 1, this composition will be filtered out before output is returned (as the `nmod` = 3). Sulphate is the only modification which is allowed to occur twice per mononer. For this, you need to set `double_sulphate=TRUE` and `nmod_max` to be at least 2.

  The different modifications and their namings are summarised below:


| **Modification**  | **Definition / description**                                                                                                                                                                                  | **IUPAC naming** | **GlycoCT naming** | **Oxford naming** |
|-------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|--------------------|-------------------|
| carboxylicacid    | Effective loss of two hydrogens and gain of one oxygen to form a carboxylic acid group on C6. The modified monomer is commonly called a 'uronic acid'                                                         | CarboxylicAcid   | COOH               | A                 |
| sialicacid        | Effect addition of C11H19N1O9 to hexose. Here, sialic acid only refers to N-Acetylneuraminic acid (Neu5Ac), the most common sialic acid. Predominantly found in complex mammalian glycans.                    | NeuAc            | SIA                | SA                |
| phosphate         |                                                                                                                                                                                                               | Phosphate        | PO4                | P                 |
| sulphate          | Addition of SO3. Only modification allowed to occur twice per monomer (see options for `double_sulphate`)                                                                                                     | Sulfate          | SO4                | S                 |
| amino             | Gain of NH and loss of of O - result ofreplacing a hydroxyl group with an amino group.                                                                                                                        | Amino            | NH2                | Am                |
| deoxy             | One hydroxyl group is replaced by an H atom. Fucose and rhamnose are two common deoxyhexoses. NB: GlycoAnnotateR currently only considers deoxyhexoses and not deoxypentoses.                                 | DeoxyHex         | DHEX               | D                 |
| nacetyl           | Addition of an N-acetyl group (net change = +C2H3N) . Common example of N-acetylated hexose is N-acetylglucosamine. Note that here, N-acetylglucosamine would be termed in e.g. IUPAC naming Hex1 N-Acetyl1.  | N-Acetyl         | NAc                | N                 |
| oacetyl           | Acetylation of a hydroxyl group (net change = +C2H2O).                                                                                                                                                        | O-Acetyl         | Ac                 | Ac                |
| omethyl           | Addition of CH2 to an hydroxyl group. Natural modification, but can also be generated by permethylation.                                                                                                      | O-Methyl         | OMe                | M                 |
| anhydrobridge     | Water loss formed by bridge between two hydroxyl groups. Occurs from C6 to C3, C2 or C1. Seen in e.g. carrageenans.                                                                                           | AnhydroBridge    | ANH                | B                 |
| unsaturated       | Water loss to form a C-C double bond inside a ring. Seen for example in ulvans and are the target of polysaccharide lyases.                                                                                   | Unsaturated      | UNS                | U                 |
| dehydrated        | Water loss that occurs during ionisation or other reactions.                                                                                                                                                  | Dehydrated       | Y                  | Y                 |
| alditol           | Reducing end monomer is opened and the aldehyde reduced to an alcohol. Commonly done before PGC-LC to reduce anomer splitting of peaks. Refers to an alditol 'modification' not a monomer here.               | Alditol          | ALD                | o                 |
| aminopentyllinker | Functional group used in synthetic chemistry. Can occur once per composition.                                                                                                                                 | NH2Pent1         | NH2Pent1           | NH2Pent1          |





