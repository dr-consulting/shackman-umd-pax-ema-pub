#######################################################################################################################
# Study 1 Psychometrics Script: Dispositional Negativity Measurement Properties

# Description: 
#   Analysis that explores the psychometric properties of dispositional negativity", with a focus on omega. Analysis
#   relies on Revelle's implementation of omega available in the psych package. 

#######################################################################################################################

load("~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study1_data.RData")

# "R" denotes a reverse-scored item 
DN_items <- c(
  "Battery_BFINeurot_1",
  "Battery_BFINeurot_2R",
  "Battery_BFINeurot_3",
  "Battery_BFINeurot_4",
  "Battery_BFINeurot_5R",
  "Battery_BFINeurot_6",
  "Battery_BFINeurot_7R",
  "Battery_BFINeurot_8",
  "Battery_IPIPAnx_9R",
  "Battery_IPIPAnx_10",
  "Battery_IPIPAnx_11",
  "Battery_IPIPAnx_12",
  "Battery_IPIPAnx_13R",
  "Battery_IPIPAnx_14R",
  "Battery_IPIPAnx_15",
  "Battery_IPIPAnx_16",
  "Battery_IPIPAnx_17R",
  "Battery_IPIPAnx_18")

# Get item-level statistics
psych::alpha(dat.study1_DN_items[DN_items])

# Calculate omega 
psych::omega((dat.study1_DN_items[DN_items]))
