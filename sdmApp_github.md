Species Distribution Models Application (**sdmApp**)
================
Aboubacar HEMA, Babacar NDAO, Louise LEROUX, Abdoul Aziz DIOUF
16/09/2020

# Introduction and Main Features

Package *sdmApp* contains a shiny app that should help users that are
non-experts in R (command-line) to apply species distribution models
techniques. For this reason, users may upload data files (raster images and species
occurence) from different software products into the app and then
start by working through an interactive graphical user interface (GUI).
This document give you an overview of the functionalities of the
graphical user interface which can be started with sdmApp(). The main
functionality of the GUI is:

  - Upload data (raster and species occurence files)
  - View correlation between raster images
  - Use Climate Ecological Niche Factors Analysis (CENFA) to select
    species predictors
  - Assignment of the spatial blocks to cross-validation folds that can
    be done in two different ways: randomly and systematically using the
    blockCV package
  - Use of the two methods: spatial blocking  and no blocking
  - Export results
  - Keep reproducibility (R code) to be able to download the
    underlying code from sdmApp

The following lines describe the main features of the interactive graphical user interface. The GUI is separated into 6 main categories, which can be
selected from the navigation bar at the top of the screen. Initially,
some of these pages will be empty until data
(both raster and species occurence files) have been uploaded.

<!--![caption](sticker_vstone_07_1.png)-->

# Help/About

This is the first page that is shown once the graphical user interface
has been started using sdmApp() after loading package sdmApp. On this
page, users will found information about the package and
contact or feedback.

<img src="image1.png" align="center" width="500" />

# Data Upload

On this page, the user can either upload data sets stored as files on
the hard drive into the GUI or to select data frames that exist in the
users’ workspace before working the graphical user interface was
started.

We note that the content of this page changes depending on whether data
have already been uploaded or not. This is described in chapter Upload
data below.

## Select data source

<img src="image_data_source.png" align="center" width="200" />

### Testdata/internal data

This screen allows the user to select data that are included in sdmApp.
Pressing the action button below the drop-down selection input will make
the GUI use the selected data.In Niakhar data are included rasters and
species occurence.

## SPSS-file (.sav)

Here users can opt to upload a file exported from SPSS. Users can change
the options if character vectors should be automatically converted to
factors and if variables that contain only missing-values (‘NA’) only
should be dropped. By clicking on the Browse button the user needs to
select a sav-file on disk which he wants to upload.

## SAS-file (.sasb7dat)

Here users can opt to upload a file exported from SAS. Users can change
the options if character vectors should be automatically converted to
factors and if variables that contain only missing-values (‘NA’) only
should be dropped. By clicking on the Browse button the user needs to
select a sas7bdat-file on disk which he wants to upload.

## CSV-file (.csv, .txt)

Here users can opt to upload a text file where variables are separated
by some characters. Typically these data would be exported from software
such as Excel. It is crucial that users indicate if the data file has
variable names in the first row and how variables are separated. At this
point, users have the option to have character vectors automatically
converted to factor or have variables that contain only missing-values
(‘NA’) dropped when the data are read into the GUI. For columns read
as character (text), the character " is ignored as quoting character and
not imported. By clicking on the Browse button the user needs to select
a txt or csv-file on disk which he wants to upload.

### Select the field separator (Comma, Semicolon, Tab)

This option is only available when a text/csv file is imported. The
radio button input has three possible choices, Comma (the default value)
Semicolon and Tab defining the value that is used to separate variables
in the input file.

  - Comma: the “,” character is used as separator
  - Semicolon: the ; character is used as separator
  - Tab: tabulators ( are used as separators

## STATA-file (.dta)

Here users can opt to upload a file exported from Stata. Users can
change the options if character vectors should be automatically
converted to factors and if variables that contain only missing-values
(‘NA’) only should be dropped. By clicking on the Browse button the
user needs to select a dta-file on disk which he wants to upload.

# Data Preparation

## View Data

### Species Data

### Species distribution

## Spatial Analysis

### Correlation between rasters (predictors)

### ENFA

# Spatial blocking

## Barchart

## Mapplot

## Range

## Variogram

# Modeling

## Profils method

### Bioclim

### Domain

### Mahalanobis distance

## Classical regression models

### GLM

### GAM

## Machine learning methods

### MaxEnt

### BRT

### Random Forest

### SVM

## Combining model predictions

# R-code
