[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)

# Reads2Map Tools

Contains utilities functions for [Reads2Map](https://github.com/Cristianetaniguti/Reads2Map) workflows usage:

- Uses as inputs VCF files to run [updog](https://github.com/dcgerard/updog), [polyRAD](https://github.com/lvclark/polyRAD) and [SuperMASSA](https://bitbucket.org/orserang/supermassa.git/src) and converts the output to onemap object or VCF files;
- Simulates mapping populations using PedigreeSim software and converts the output to onemap object or VCF files;
- Draws samples haplotypes from phased VCF files.

# How to install

```R
install.packages("devtools")
devtools::install_github("Cristianetaniguti/Reads2MapTools")
```

This package is one of the tools used in [Reads2Map](https://github.com/Cristianetaniguti/Reads2Map) workflows. You can find it in the docker image [cristaniguti/reads2map](https://hub.docker.com/repository/docker/cristaniguti/reads2map). Usage through docker:

```bash
docker pull cristaniguti/reads2map:0.0.2
```

# Documentation

Check features description and example of usage [here](https://cristianetaniguti.github.io/Tutorials/Reads2MapTools/Simulations.html).

