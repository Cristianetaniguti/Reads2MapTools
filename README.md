[![Build Status](https://travis-ci.com/Cristianetaniguti/onemap.svg?branch=master)](https://travis-ci.org/Cristianetaniguti/onemap) 
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)

# OneMap Utils

Contains utilities functions for [OneMap](https://github.com/augusto-garcia/onemap) package usage:

- Uses as inputs VCF files to run updog, polyRAD and SuperMASSA and converts the output to onemap object;
- Simulates mapping populations using PedigreeSim software and converts the output to onemap object;
- Draws samples haplotypes from phased VCF files.

# How to install

```R
install.packages("devtools")
devtools::install_github("Cristianetaniguti/onemapUTILS")
```

This package is one of the tools used in [Reads2Map]() workflows. You can find it in the docker image [cristaniguti/reads2map](). Usage through docker:

```bash
docker pull cristaniguti/reads2map:latest
```

# Tutorials


