#Code from: "A global meta-analysis of ecological functions and regulating ecosystem services of freshwater bivalves"
#Alexandra Zieritz, et al., Limnology and Oceanography, 2025

#Code prepared by Joshua I. Brian, Ana Sofia Vaz

# The code performs the following:
#   1. Calculates pooled effect sizes across all case studies using random-effects model (REMA) 
#   2. Tests for publication bias on REMA
#   3. Calculates a multilevel-meta-analyis (MLMA) to account for remaining non-independence and compares REMA and MLMA
#   4. Performs structures meta-analyses to explain heterogeneity by adding moderators to REMA and MLMA

#NOTE: This code runs through the full procedure using the complete data (Zieritz_et_al_data.csv).
#Results for each EFRES category individually can be obtained by filtering the full database per 
#category and repeating this code.

library("metafor")
library("plyr")

cst<-read.csv("Zieritz_et_al_data.csv", header = T)

# Create a unique id for each effect size
cst$id<-as.factor(seq(from = 1, to = nrow(cst), by = 1))
str(cst)

####### 1. GRAND EFFECT: RANDOM-EFFECTS META-ANALYSIS (REMA) ##############

res<-rma(yi, vi, data=cst) 
summary(res) 
permutest(res, exact = FALSE, iter=1000, retpermdist = TRUE) #apply permutation test for significance 

#Forest plot
forest(res,
       back="lightgray", shade="white", hlines="white",
       pch=19, pch.fill=21, main = "Chemicals")

#Funnel plot
funnel(res, yaxis="sei", xlab="effect size", ylab="standard error",
       back="lightgray", shade="white", hlines="white",
       pch=19, pch.fill=21, main = "Biological control") # Save figure in a folder

#Check forest and funnel results, If some case studies emerge as outliers, try removing them, and re-run the REMA

#qqnorm(residuals(res))
#qqline(residuals(res))

####### 2. Testing for Publication bias tests on REMA ###################

#Egger's regression test
regtest(res, model = "rma", predictor = "sei", ret.fit=TRUE) #Copy results to a table

#Trim and fill
tf.res<-trimfill(res)
structure(tf.res)
funnel(tf.res, main = "alldata") 

########3. MULTILEVEL META-ANALYSIS (MLMA) TO ACCOUNT FOR SOURCES OF NON-INDEPENDENCY #########

# Get within study variance (needed below for I2 partition)
wsv<-((sum(1/cst$vi)) * (length(cst$vi)-1)) / (((sum(1/cst$vi))^2) - (sum((1/cst$vi)^2)))

# Tests for a source of non-independence (random factor), use the following moderators:
#paper
#country
#country2
#continent
#experiment.type
#system.type
#experiment.design
#other.fb
#taxa
#genera
#order
#origin
#significant

#Test each random factor
mlma<-rma.mv(yi, vi, random=~1|paper/id, data=cst) 
summary(mlma) 

mlma<-rma.mv(yi, vi, random=~1|country/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|country2/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|continent/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|experiment.type/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|system.type/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|experiment.design/id, data=cst)
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|other.fb/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|taxa/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|genera/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|order/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|origin/id, data=cst) 
summary(mlma)

mlma<-rma.mv(yi, vi, random=~1|significant/id, data=cst) 
summary(mlma)

########4. STRUCTURED META-ANALYSIS TO EXPLAIN HETEROGENEITY WITH MODERATORS #########

## REMA with categorical moderators (repeat for each moderator)
#country
#country2
#continent
#experiment.type
#system.type
#experiment.design
#other.fb
#taxa
#genera
#order
#origin
#significant

#res<-rma(yi, vi, mods = ~
         #+factor(continent)
         #+factor(order)
         #+factor(origin)
         #+factor(experiment.type)
         #+factor(system.type)
         #- 1, data = cst)

res<-rma(yi, vi, mods = ~ factor(country) - 1, data = cst) 
summary(res)
count(cst$country) 

res<-rma(yi, vi, mods = ~ factor(country2) - 1, data = cst)
summary(res)
count(cst$country2)

res<-rma(yi, vi, mods = ~ factor(continent) - 1, data = cst) 
summary(res)
count(cst$continent)

res<-rma(yi, vi, mods = ~ factor(experiment.type) - 1, data = cst) 
summary(res)
count(cst$experiment.type)

res<-rma(yi, vi, mods = ~ factor(system.type) - 1, data = cst)
summary(res)
count(cst$system.type)

res<-rma(yi, vi, mods = ~ factor(experiment.design) - 1, data = cst) 
summary(res)
count(cst$experiment.design)

res<-rma(yi, vi, mods = ~ factor(other.fb) - 1, data = cst) 
summary(res)
count(cst$other.fb)

res<-rma(yi, vi, mods = ~ factor(taxa) - 1, data = cst) 
summary(res)
count(cst$taxa)

res<-rma(yi, vi, mods = ~ factor(genera) - 1, data = cst) 
summary(res)
count(cst$genera)

res<-rma(yi, vi, mods = ~ factor(order) - 1, data = cst) 
summary(res)
count(cst$order)

res<-rma(yi, vi, mods = ~ factor(origin) - 1, data = cst) 
summary(res)
count(cst$origin)

res<-rma(yi, vi, mods = ~ factor(significant) - 1, data = cst) 
summary(res)
count(cst$significant)

## 5. REMA with quantitative moderators (repeat for each moderator)

### Quantitative moderators
#size.area
#tank.size
#water.volume
#water.depth
#sediment.depth
#chlorophyll
#nitrogen
#phosphorous
#carbon
#o2
#temperature
#pH
#soft.tissue
#dry.weight
#weight.g.m2
#total.weight
#number.fb
#density.fb.volume
#density.fb.area
#length

#first, ensure moderators are numeric

cst$size.area<-as.numeric(cst$size.area) 
res <- rma(yi, vi, mods = ~ size.area, data=cst) 
summary(res)

cst$tank.size<-as.numeric(cst$tank.size) 
res <- rma(yi, vi, mods = ~ tank.size, data=cst) 
summary(res)

cst$water.volume<-as.numeric(cst$water.volume) 
res <- rma(yi, vi, mods = ~ water.volume, data=cst) 
summary(res)

cst$water.depth<-as.numeric(cst$water.depth) 
res <- rma(yi, vi, mods = ~ water.depth, data=cst) 
summary(res)

cst$sediment.depth<-as.numeric(cst$sediment.depth) 
res <- rma(yi, vi, mods = ~ sediment.depth, data=cst) 
summary(res)

cst$chlorophyll<-as.numeric(cst$chlorophyll) 
res <- rma(yi, vi, mods = ~ chlorophyll, data=cst) 
summary(res)

cst$nitrogen<-as.numeric(cst$nitrogen) 
res <- rma(yi, vi, mods = ~ nitrogen, data=cst) 
summary(res)

cst$phosphorous<-as.numeric(cst$phosphorous) 
res <- rma(yi, vi, mods = ~ phosphorous, data=cst) 
summary(res)

cst$carbon<-as.numeric(cst$carbon) 
res <- rma(yi, vi, mods = ~ carbon, data=cst) 
summary(res)

cst$o2<-as.numeric(cst$o2) 
res <- rma(yi, vi, mods = ~ o2, data=cst) 
summary(res)

cst$temperature<-as.numeric(cst$temperature) 
res <- rma(yi, vi, mods = ~ temperature, data=cst) 
summary(res)

cst$pH<-as.numeric(cst$pH) 
res <- rma(yi, vi, mods = ~ pH, data=cst) 
summary(res)

cst$soft.tissue<-as.numeric(cst$soft.tissue) 
res <- rma(yi, vi, mods = ~ soft.tissue, data=cst) 
summary(res)

cst$dry.weight<-as.numeric(cst$dry.weight) 
res <- rma(yi, vi, mods = ~ dry.weight, data=cst) 
summary(res)

cst$weight.g.m2<-as.numeric(cst$weight.g.m2) 
res <- rma(yi, vi, mods = ~ weight.g.m2, data=cst) 
summary(res)

cst$total.weight<-as.numeric(cst$total.weight) 
res <- rma(yi, vi, mods = ~ total.weight, data=cst) 
summary(res)

cst$number.fb<-as.numeric(cst$number.fb) 
res <- rma(yi, vi, mods = ~ number.fb, data=cst) 
summary(res)

cst$density.fb.volume<-as.numeric(cst$density.fb.volume) 
res <- rma(yi, vi, mods = ~ density.fb.volume, data=cst) 
summary(res)

cst$density.fb.area<-as.numeric(cst$density.fb.area) 
res <- rma(yi, vi, mods = ~ density.fb.area, data=cst) 
summary(res)

cst$length<-as.numeric(cst$length) 
res <- rma(yi, vi, mods = ~ length, data=cst) 
summary(res)



