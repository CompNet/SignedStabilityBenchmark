#!/bin/bash

echo '.libPaths(c("~/libs/R-server", .libPaths()))' > ~/.Rprofile

/usr/bin/Rscript src/main2.R
