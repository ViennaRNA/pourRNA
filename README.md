[![GitHub release](https://img.shields.io/github/release/ViennaRNA/pourRNA.svg)](https://github.com/ViennaRNA/pourRNA/releases)
[![Build Status](https://travis-ci.org/ViennaRNA/pourRNA.svg?branch=master)](https://travis-ci.org/ViennaRNA/pourRNA)
[![Github All Releases](https://img.shields.io/github/downloads/ViennaRNA/pourRNA/total.svg)](https://github.com/ViennaRNA/pourRNA/releases)
[![Conda](https://img.shields.io/conda/v/bioconda/pourRNA.svg)](https://anaconda.org/bioconda/pourRNA)

# pourRNA

A tool that computes local minima and respective transition rates of an RNA energy landscape.

## Installation

To configure, compile and install execute the following commands on your command line:
```
  - autoreconf -i
  - ./configure [--help for additional configuration options]
  - make
  - make install
```
Dependencies:
  - [ViennaRNA library (>= 2.4.11)](https://www.tbi.univie.ac.at/RNA/#download)
  
## Execute
The minimal input is an RNA sequence.
```
pourRNA --seq="ACGUACGUACGUACGU"
```
The output consists of the representative structures of the local minima and the transition rates between adjacent local minima.
You can adjust the output and speed up the computation by using additional command line parameters. All parameters are shown by
```
pourRNA --help
```

In case you run into problems, please contact us!



&copy; Gregor Entzian, 2019

