### exforUncertainty - R package

This package contains a collection of functions to
extract systematic and statistical uncertainties 
from entries in the EXFOR database.

It also contains functions to derive (in lack of 
a better word) uncertainties from those present
in EXFOr subentries. Namely, inconsistent 
uncertainty specifications are penalized by the
introduction of extra uncertainty and conservative
assumptions are made about missing uncertainties.

The assessment of uncertainties of nuclear experiments
is a non-trivial undertaking and hence the functions
for the automatic construction and correction of 
uncertainties should be seen as a first blueprint 
and are likely to change in the future.

## General notes

This package operates on nested lists representing 
EXFOR subentries.
To obtain the nested list of correct structure,
the R package [exforParser](https://github.com/gschnabel/exforParser)
has been used to convert EXFOR given in text format to
nested lists. Specifically, first the function `parseEntry`
has been applied and subsequently the function `transformSubent`.
See the repo for [createExforDb](https://github.com/gschnabel/createExforDb)
for an example of this conversion.
