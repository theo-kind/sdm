# ================================= Notes =================================
#
# Setup script called at the start of every pipeline script.
# Loads R packages, creates the project and folder structure
# folder structure via envimaR, and sources all helper functions.
#
# ================================= Set up =================================

# list of packages to load
packagesToLoad <- c(
  'terra',
  'png',
  'sf',
  'rsample',
  'purrr',
  'stars',
  'magick',
  'fs',
  'tidyverse',
  'dplyr',
  'virtualspecies',
  'predicts',
  'blockCV',
  'USE',
  'maxnet',
  'geodata'
)

# mandantory folder structure
projectDirList <- c(
  'data/',
  'data/raster/',
  'data/virtual_species/',
  'docs/',
  'docs/figures/',
  'run/',
  'tmp/',
  'src/',
  'src/functions/'
)


# append additional folders if defined by calling script
if (exists("appendProjectDirList") && appendProjectDirList[[1]] != "") {
  projectDirList <- append(projectDirList, appendProjectDirList)
}

# append additional packages if defined by calling script
if (exists("appendPackagesToLoad") && appendPackagesToLoad[[1]] != "") {
  packagesToLoad <- append(packagesToLoad, appendPackagesToLoad)
}

library(envimaR)

# Automatically set root direcory, folder structure and load libraries
envrmt <- envimaR::createEnvi(
  root_folder = rootDir,             # if it exists, this is the root folder
  folders = projectDirList,          # mandantory folder structure
  path_prefix = "path_",             # prefix to all path variables that are created
  libs = packagesToLoad,             # list of R packages to be loaded
  alt_env_id = "COMPUTERNAME",       # check the environment variable "COMPUTERNAME"
  alt_env_value = "PCRZP",           # check if it contains the string "PCRZP" (e.g. local PC pools)
  alt_env_root_folder = "F:/BEN/edu" # use the alternative root folder
)

# set temp path to speed up terra package operations
terraOptions(tempdir = envrmt$path_tmp)

# load functions from the functions folder
lapply(
  list.files(
    file.path(envrmt$path_functions),
    pattern = "\\.[Rr]$",
    full.names = TRUE,
    recursive = TRUE
  ),
  source
)

