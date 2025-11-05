## ------------------------------------------------------- ##
# Above-Below Synchrony - Setup Tasks
## ------------------------------------------------------- ##
# Purpose:
## Do generally useful setup stuff that multiple downstream scripts will benefit from

# Clear environment
rm(list = ls()); gc()

## -------------------------------------- ##
# Create Folders ----
## -------------------------------------- ##

# Make needed folder(s)
dir.create(path = file.path("data", "raw_neon"), showWarnings = F, recursive = T)

# End ----
