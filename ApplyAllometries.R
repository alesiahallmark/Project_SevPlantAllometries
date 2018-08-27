### This script matches quad observations with known and inferred allometries to estimate biomass 

# Load required libraries
library(ggplot2)
library(RColorBrewer)
library(lubridate)
library(plyr)
library(reshape2)
library(gridExtra)
library(cowplot)
library(knitr)
library(markdown)
library(xtable)
library(zoo)
library(tidyr)
library(stats)
library(car)
library(nlme)
library(lme4)
library(visreg)
library(data.table)

# Read in all quad data
allquads <- read.csv("/Users/alesia/Desktop/unsorted_AJdata/AJ_AllQuads.csv", strip.white = T)

allquads$kartez <- as.factor(toupper(allquads$kartez))

allquads$cover[allquads$cover < 0] <- NA
allquads$height[allquads$height < 0] <- NA

allquads$volume <- allquads$cover * allquads$height

allquads$site <- revalue(allquads$site, c(
  "BURN03" = "B03",
  "BURN09" = "B09"))

allquads$SiteCluster <- revalue(allquads$site, c(
  "B" = "BlueGrama",
  "BG" = "DeepWell",
  "B03" = "DeepWell",
  "B09" = "DeepWell",
  "C" = "FivePoints",
  "CG" = "FivePoints",
  "EB" = "BlueGrama",
  "EG" = "FivePoints",
  "F" = "DeepWell",
  "G" = "FivePoints",
  "GG" = "FivePoints",
  "M" = "FivePoints",
  "MG" = "DeepWell",
  "MS" = "FivePoints",
  "N" = "DeepWell",
  "P" = "CerroMontosa",
  "W" = "DeepWell"))

# Fix some species names
allquads$kartez <- revalue(allquads$kartez, c(
  "SATR2" = "SATR12",
  "SPP06" = "SPPO6", 
  "OPUN" = "OPUNT",
  "ASMIM" = "ASFE2",
  "CHGR2" = "DYGR",
  "LEDED" = "LEDE",
  "GACO5" = "OESU3",
  "GASUN" = "OESUN",
  "GUWR" = "GUSP",
  "PHIN" = "PHCR",
  "SPWR" = "SPPO6"))

# Match site names with nearest met stations
allquads$MetStation <- 
  revalue(allquads$SiteCluster,
          c(DeepWell = 40, # Deep well
            CerroMontosa = 42, # Cerro Montosa PJ
            FivePoints = 49, # Five Points
            BlueGrama = 50)) # Blue grama

allquads$season <- revalue(as.factor(allquads$season), c(
  "1" = "winter", "2" = "spring", "3" = "fall"))

# Get rid of winter observations for now
allquads <- allquads[allquads$season %in% c("spring", "fall"),]
corequads <- allquads[allquads$site %in% c("B", "G", "C") &
                        allquads$dataset == "sev129",]

unique(corequads$kartez[is.na(corequads$genus)])


# Read in allometry data
allos <- read.csv("~/Desktop/SevPlantBiomass_06Jul18_highlightsAPPLY.csv", header = T)
invar.allos <- allos[allos$Best.Invar.Model == "x", c("kartez", "genus", "sp.epithet", "family", "LifeHistory", "PhotoPath", "FunctionalGroup", "season", "model.name", "main.eff.full", "mainclim.eff.full", "inter.eff.full")]

# replicate allometries for missing seasons
ex.sp.se <- expand.grid(unique(invar.allos$kartez), as.factor(c("spring", "fall")))
colnames(ex.sp.se) <- c("kartez", "season")
ex.sp.se <- unique(merge(ex.sp.se, invar.allos[,c("kartez", "genus", "sp.epithet", "family", "LifeHistory", "PhotoPath", "FunctionalGroup")], all = T))

invar.allos <- merge(invar.allos, ex.sp.se, all = T)

for (i in 1:length(unique(invar.allos$kartez))) {
  invar.allos$model.name[
    invar.allos$kartez == invar.allos$kartez[i] & 
      is.na(invar.allos$model.name)] <- 
    invar.allos$model.name[
      invar.allos$kartez == invar.allos$kartez[i] & 
        !is.na(invar.allos$model.name)]
  invar.allos$main.eff.full[
    invar.allos$kartez == invar.allos$kartez[i] & 
      is.na(invar.allos$main.eff.full)] <- 
    invar.allos$main.eff.full[
      invar.allos$kartez == invar.allos$kartez[i] & 
        !is.na(invar.allos$main.eff.full)]
}

# Read in taxonomy data and merge
taxo <- read.csv("/Users/alesia/Documents/Project_SevPublicPlantRCode/SevilletaSpeciesList_AJH.csv", strip.white = T, na.strings = c("NA",""))
# Revalue a_p and g_f columns
taxo$a_p <- revalue(taxo$a_p, c(
  "a" = "annual",
  "p" = "perennial",
  "a/p" = "ann/peren"))
taxo$g_f <- revalue(taxo$g_f, c(
  "f" = "forb",
  "t" = "tree", 
  "g" = "grass",
  "s" = "shrub"))

# Create temporary substitution allometries (sister species, other season, average of larger taxa)
# Find average slope for each PFG
PFG.means <- ddply(allos[allos$model.name %in% c("Volume"),], c("LifeHistory", "PhotoPath", "FunctionalGroup", "season", "model.name"), function(x) data.frame(
  main.eff.full = mean(x$main.eff.full, na.rm = T)
))
# Apply other-season allometry from missing PFG-season combos
PFG.means[dim(PFG.means)[1] + 1,] <- PFG.means[3,]
PFG.means$season[dim(PFG.means)[1]] <- "spring"
PFG.means[dim(PFG.means)[1] + 1,] <- PFG.means[4,]
PFG.means$season[dim(PFG.means)[1]] <- "spring"
PFG.means[dim(PFG.means)[1] + 1,] <- PFG.means[4,]
PFG.means$PhotoPath[dim(PFG.means)[1]] <- "C3"
PFG.means[dim(PFG.means)[1] + 1,] <- PFG.means[20,]
PFG.means$PhotoPath[dim(PFG.means)[1]] <- "C3"

# Find observations with no matching allometry
find.miss <- merge(corequads, invar.allos, all.x = T)
find.miss <- find.miss[is.na(find.miss$main.eff.full),]
unique(find.miss$kartez)

substitution.allos <- expand.grid(unique(find.miss$kartez), c("spring", "fall"))
colnames(substitution.allos) <- c("kartez", "season")
substitution.allos <- unique(merge(substitution.allos, taxo[,c("kartez", "family", "genus", "species", "path", "a_p", "g_f")], all.x = T))
colnames(substitution.allos) <- c("kartez", "season", "family", "genus", "sp.epithet", "PhotoPath", "LifeHistory", "FunctionalGroup")

# Apply average slope for that functional group
substitution.allos <- merge(substitution.allos, PFG.means, by = c("LifeHistory", "PhotoPath", "FunctionalGroup", "season"), all.x = T)


# Apply a few special allometries - when we can do better than PFG mean


# Find groups that still don't have allometry and apply a substitute
substitution.allos[is.na(substitution.allos$main.eff.full),]

# Combine real and imagined allometries
all.invar.allos <- merge(invar.allos, substitution.allos, all = T)

# Merge quad data with allometries
corequads <- merge(corequads, taxo[,c("kartez", "family", "genus", "species", "path", "a_p", "g_f")], all.x = T)
quad.allo <- merge(corequads, invar.allos, all.x = T)


# Calculate 3 versions of biomass data
quad.allo$biomass <- NA
quad.allo$biomass[quad.allo$model.name %in% "Volume"] <- 
  quad.allo$volume[quad.allo$model.name %in% "Volume"] * 
  quad.allo$main.eff.full[quad.allo$model.name %in% "Volume"]
quad.allo$biomass[quad.allo$model.name %in% "Cover"] <- 
  quad.allo$cover[quad.allo$model.name %in% "Cover"] * 
  quad.allo$main.eff.full[quad.allo$model.name %in% "Cover"]

# Find illegal or outlier values

# Count number of quads at each site
num.quads <- ddply(quad.allo, c("site", "year"), function(x) data.frame(
  total.quads = length(unique(paste(x$web, x$plot, x$subplot)))
))
num.quads$total.quads[num.quads$total.quads %in% c(38, 39, 41)] <- 41
num.quads$total.quads[num.quads$total.quads %in% c(79)] <- 80

quad.allo <- merge(quad.allo, num.quads, all.x = T)

# Summarize biomass by quad and species
quad.sp.summ <- as.data.frame(data.table(quad.allo)[, list(
  volume = sum(volume, na.rm = T),
  cover = sum(cover, na.rm = T),
  biomass = sum(biomass, na.rm = T)),
  by = list(year, site, season, kartez, web, plot, treatment, dataset, subplot, SiteCluster, MetStation, genus, sp.epithet, family, LifeHistory, PhotoPath, FunctionalGroup, total.quads)])

avg_dates <- as.data.frame(data.table(quad.allo)[, list(
  date = mean(as.Date(date), na.rm = T)),
  by = list(year, site, season)])

quad.sp.summ <- merge(quad.sp.summ, avg_dates, all.x = T)


# Save for phenomass
write.csv(quad.sp.summ, "~/Desktop/unsorted_AJdata/corebiomass_TEMP.csv", row.names = F)

# Summarize biomass by site and species
site.sp.summ <- as.data.frame(data.table(quad.sp.summ)[, list(
  volume = sum(volume, na.rm = T) / total.quads,
  cover = sum(cover, na.rm = T) / total.quads,
  biomass = sum(biomass, na.rm = T) / total.quads),
  by = list(year, site, season, date, kartez, treatment, dataset, SiteCluster, MetStation, genus, sp.epithet, family, LifeHistory, PhotoPath, FunctionalGroup, total.quads)])


# Summarize biomass by site and PFG
site.PFG.summ <- as.data.frame(data.table(site.sp.summ)[, list(
  volume = sum(volume, na.rm = T),
  cover = sum(cover, na.rm = T),
  biomass = sum(biomass, na.rm = T)),
  by = list(year, date, site, season, treatment, dataset, SiteCluster, MetStation, LifeHistory, PhotoPath, FunctionalGroup)])
site.PFG.summ$PFG <- paste(site.PFG.summ$LifeHistory, site.PFG.summ$PhotoPath, site.PFG.summ$FunctionalGroup)
site.PFG.summ$PFG.PF <- paste(site.PFG.summ$PhotoPath, site.PFG.summ$FunctionalGroup)


ggplot(site.sp.summ, aes(x = date, y = biomass, group = kartez)) + 
  geom_col(width = 100, aes(fill = kartez), show.legend = F) + 
  facet_grid(site~.)

ggplot(site.PFG.summ[site.PFG.summ$site %in% "G",], aes(x = date, y = biomass, group = PFG.PF)) + 
  geom_col(width = 110, aes(fill = PFG.PF)) + 
  facet_grid(site~., scales = "free_y")

