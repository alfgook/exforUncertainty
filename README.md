### exforUncertainty - R package

This package contains a collection of functions to
extract systematic and statistical uncertainties 
from entries in the EXFOR database.
It also contains functions to derive (in lack of 
a better word) uncertainties from those present
in EXFOR subentries. Namely, inconsistent 
uncertainty specifications are penalized by the
introduction of an extra uncertainty and conservative
assumptions are made about missing uncertainties.

The assessment of uncertainties of nuclear experiments
is a non-trivial undertaking and hence the functions
in this package for the automatic construction and 
correction of uncertainties should be seen as a first
attempt and are likely subject to change in the future.

This package operates on nested lists representing 
EXFOR subentries.
To obtain nested lists of correct structure,
the R package [exforParser](https://github.com/gschnabel/exforParser)
can be used to convert EXFOR entries given in text format to
nested lists. Specifically, first the function `parseEntry`
can be applied and subsequently the function `transformSubent`.
See the repo [createExforDb](https://github.com/gschnabel/createExforDb)
for an example of this conversion.
