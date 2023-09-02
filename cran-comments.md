## Test environments
* local macOS R installation, R 4.3.1
* macOS latest release (via [R-CMD-check](https://github.com/astamm/nevada/actions/workflows/R-CMD-check.yaml) github action)
* windows latest release (via [R-CMD-check](https://github.com/astamm/nevada/actions/workflows/R-CMD-check.yaml) github action)
* ubuntu 20.04 latest both release and devel (via [R-CMD-check](https://github.com/astamm/nevada/actions/workflows/R-CMD-check.yaml) github action)
* [win-builder](https://win-builder.r-project.org/) (release and devel)
* [R-hub](https://builder.r-hub.io)
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran

## R CMD check results
There were no ERRORs or WARNINGs.

There was 2 NOTEs:

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Aymeric Stamm <aymeric.stamm@math.cnrs.fr>'

    New submission

    Possibly misspelled words in DESCRIPTION:
      Lovato (35:74, 37:60)
      Multiscale (38:48)
      Pini (36:9, 37:72)
      Stamm (36:19, 38:5)
      Taquet (38:16)
      Vantini (36:30, 38:28)
      flipr (25:67)

    * checking installed package size ... NOTE
      installed size is  5.2Mb
      sub-directories of 1Mb or more:
        libs   4.7Mb
