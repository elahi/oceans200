# Just for  working with Rmarkdown; not needed for analysis: "knitr"
# Load tidyverse to make things simple: "tidyverse", "broom", "broom.mixed"
# Common lm add-ons: "car", "MASS", "Rmisc", "ez", "afex", "emmeans", "lme4", "lmerTest", "nlme", "VCA", "MuMIn", "lattice", "effectsize", "lmtest", 
# Graphics enhancements: "ggExtra", "patchwork", "viridis", "ggsci", "ggforce"

# Package names
packages <- c("knitr", "tidyverse", "broom", "broom.mixed", "sjstats", "pwr", 
              "car", "MASS", "Rmisc", "ez", "afex", "emmeans", "lme4", "lmerTest", "nlme", "VCA", "MuMIn", "lattice", "effectsize", "lmtest", 
              "ggExtra", "patchwork", "viridis", "ggsci", "ggforce") 

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {  install.packages(packages[!installed_packages]) }

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
