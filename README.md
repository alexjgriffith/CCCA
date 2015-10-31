Cross Condition ChIP Analysis (CCCA)
==
Alexander Griffith, University of Ottawa
--

This system is composed of two portions


1. First the combination of data sets
   - Generate a unified peak set
   - determine the read pile up count
     under each peak
2. Extraction of relevent subsets

```r
library(CCCA)
bedData<-loadBedFile("test.bed")
bedData[,1]
```

This file is part of CCCA, http://github.com/alexjgriffith/CCCA/, and is Copyright (C) University of Ottawa, 2015. It is Licensed under the three-clause BSD License; see LICENSE.txt.

