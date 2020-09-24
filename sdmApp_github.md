Species Distribution Models Application (**sdmApp**)
================
Aboubacar HEMA, Babacar NDAO, Louise LEROUX, Abdoul Aziz DIOUF
16/09/2020

# Introduction and Main Features

Package *sdmApp* contains a shiny app that should help users that are
non-experts in R (command-line) to apply species distribution models
techniques. For this reason, users may upload data (rasters and spcies
occurence) files from different software products into the app and then
start by working within the interactive, graphical user interface (GUI).
This document will give an overview of the functionalities of the
graphical user interface which can be started with sdmApp(). The main
functionality of the GUI is:

  - Uploading data (raster and species occurence files)
  - View correlation between raster
  - Using Climate Ecological Niche Factors Analysis (CENFA) to select
    species predictors
  - The assignment of the spatial blocks to cross-validation folds can
    be done in two different ways: random and systematic by using
    blockCV package
  - Using two methods both spatial blocking strategy and no blocking
  - Export results
  - Keep reproducibility (R code) by being able do download the
    underlying code from sdmApp

We now describe the features of the interactive graphical user interface
in detail. The GUI is separated into 6 main categories, which can be
selected from the navigation bar at the top of the screen. Initially,
some of these pages will be empty and their content changes once data
(both raster and species occurence files) have been uploaded.


<!--![caption](sticker_vstone_07_1.png)-->

# Help/About <img src="image1.png" align="right" width="120" />

This is the first page that is shown once the graphical user interface
has been started using sdmApp() after loading package sdmApp. On this
page, the user is presented the information on how to open this package
vignette which contains extensive information on how to use the GUI.

# Data Upload

On this page, the user can either upload data sets stored as files on
the hard drive into the GUI or to select data frames that exist in the
users’ workspace before working the graphical user interface was
started. This allows to perform common data manipulation steps directly
in R before continuing to anonymize the dataset using the GUI.

We note that the content of this page changes depending on whether
microdata have already been uploaded or not. In the former case, the
user can view, modify or reset variables from the uploaded dataset as
described in chapter Modify microdata. In the latter case, the user is
asked to upload data in the GUI. This is described in chapter Upload
microdata below.

## Testdata/internal data

This screen allows the user to select data.frames that are available in
the users-workspace when starting the user interface. Two test-data sets
(testdata and testdata2, information on which is available from
?testdata) that are included in sdmApp are always available. Pressing
the action button below the drop-down selection input will make the GUI
use the selected data frame.

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
