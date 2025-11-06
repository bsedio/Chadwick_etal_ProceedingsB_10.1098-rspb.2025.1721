# Chadwick_etal_ProceedingsB_10.1098-rspb.2025.1721
# R code and supplementary datasets for Chadwick, Sierra E., David Henderson, Dale L. Forrister, Leslie Cayola, Alfredo F. Fuentes, Belén Alvestegui, Nathan Muchhala, J. Sebastián Tello, Martin Volf, Jonathan A. Myers, and Brian E. Sedio. "Chemical properties of foliar metabolomes represent a key axis of functional trait variation in forests of the tropical Andes" Proceedings B. https://doi.org/10.1098/rspb.2025.1721


#### R code used for Chadwick, Sierra E., David Henderson, Dale L. Forrister, 
#### Leslie Cayola, Alfredo F. Fuentes, Belén Alvestegui, Nathan Muchhala,
#### J. Sebastián Tello, Martin Volf, Jonathan A. Myers, and Brian E. Sedio 
#### Proceedings B 
#### https://doi.org/10.1098/rspb.2025.1721

## Set working directory to folder that contains supplementary datasets S1-S4
setwd("~/Documents/Madidi_Project")

## Load required R packages
library(vegan)
library(ggplot2)
library(ggrepel)
library(corrplot)
library(phytools)
library(picante)
library(phylosignal)
library(V.PhyloMaker2)
library(phylomatic)
library(phangorn)
library(ape)
library(dplyr)
library(ggtree)

## Read supplmentary datasets
metaprops = read.csv("~/Documents/Madidi_Project/Chadwick_ProcB_10.1098:rspb.2025.1721_DatasetS1.csv")

metaplotsprops = read.csv("~/Documents/Madidi_Project/Chadwick_ProcB_10.1098:rspb.2025.1721_DatasetS2.csv")
row.names(metaplotsprops) = metaplotsprops$Plot

chemtree = read.tree("~/Documents/Madidi_Project/Chadwick_ProcB_10.1098:rspb.2025.1721_DatasetS3.tre")

mastertable = read.csv("~/Documents/Madidi_Project/Chadwick_ProcB_10.1098:rspb.2025.1721_DatasetS4.csv")

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
#### For Fig. 1: Species-level correlation plot of chemical properties and morphological traits, accounting for phylogeny
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

## Simplify long variable names and extract matrix that only contains trait data
names(metaprops)[which(names(metaprops) == "TwigBarkThickness_Relative")] = "TwigBarkThick"
names(metaprops)[which(names(metaprops) == "LeafThickness")] = "LeafThick"
metaprops.foc = metaprops[,c("nAtomP", "ALogP", "TopoPSA", "Fsp3", "MW", "SLA", "LeafArea", "LeafThick", "TwigBarkThick", "TwigSpecDens")]

## Create empty matrix for recording species-level mean trait values
sppmeanprops = as.data.frame(matrix(NA, nrow = length(unique(metaprops$Species_binomial)[!is.na(unique(metaprops$Species_binomial))]), ncol = 10))
names(sppmeanprops) = c("nAtomP", "TopoPSA", "MW", "ALogP", "Fsp3", "SLA", "LeafArea", "LeafThick", "TwigBarkThick", "TwigSpecDens")
names(sppmeanprops) = c("nAtomP", "ALogP", "TopoPSA", "Fsp3", "MW", "SLA", "LeafArea", "LeafThick", "TwigBarkThick", "TwigSpecDens")
row.names(sppmeanprops) = unique(metaprops$Species_binomial)[!is.na(unique(metaprops$Species_binomial))]
dim(sppmeanprops)
# [1] 472  10

## Calculate species-level mean trait values
for(i in 1:nrow(sppmeanprops)){
	spp = row.names(sppmeanprops)[i]
	if(length(which(metaprops$Species_binomial == spp)) == 1){
		sppmeanprops[i,] = metaprops.foc[which(metaprops$Species_binomial == spp),]
	}
	if(length(which(metaprops$Species_binomial == spp)) > 1){
		sppmeanprops[i,] = colMeans(metaprops.foc[which(metaprops$Species_binomial == spp),], na.rm = T)
	}
}
sppmeanprops = sppmeanprops[!is.na(sppmeanprops$SLA),]
sppmeanprops = sppmeanprops[!is.na(sppmeanprops$TwigBarkThick),]
sppmeanprops = sppmeanprops[!is.na(sppmeanprops$TwigSpecDens),]

dim(sppmeanprops)
# [1] 437  10

## Calculate and plot simple correlation matrix
sppcorrm = cor(sppmeanprops)
corrplot(corr = sppcorrm)


## Prune Madidi species from a megaphylogeny using V.phylomaker to create 'phylo', a phylogeny of Madidi species
mad.tax = as.data.frame(cbind(metaprops$Species_binomial, metaprops$Genus, metaprops$Family))
names(mad.tax) = c("Species", "Genus", "Family")
phy = phylo.maker(mad.tax)
phylo <- phy$scenario.3

## Create a bifurcating phylogeny that only contains Madidi tree species for which we have trait data
traitvec = (t(sppmeanprops))[1,]
names(traitvec) = gsub(" ", "_", names(traitvec))
phylo.props = prune.missing(x = traitvec, phylo)$tree
phylo.props = multi2di(phylo.props)
plot(phylo.props, cex = 0.1)

## Reorganize species mean chemical properties so it is in the same order as the Madidi megaphylogeny
sppmeanprops.phylo = sppmeanprops
row.names(sppmeanprops.phylo) = gsub(" ", "_", row.names(sppmeanprops.phylo))
sppmeanprops.phylo = sppmeanprops.phylo[phylo.props$tip.label,]
dim(sppmeanprops.phylo)
phylo.props

## Create phylogenetic covariance matrix
obj<-phyl.vcv(as.matrix(sppmeanprops.phylo),as.matrix(vcv(phylo.props)),1)
obj$R

                     # nAtomP         ALogP       TopoPSA          Fsp3            MW           SLA      LeafArea
# nAtomP         6.817337e-01 -5.385959e-02    3.62378873 -2.097556e-02  2.719911e+00  -0.565412277 -2.827712e+02
# ALogP         -5.385959e-02  4.209581e-02   -0.62852439  3.241475e-03  8.192659e-01   0.033713204 -8.007950e+01
# TopoPSA        3.623789e+00 -6.285244e-01   37.27276427 -1.068895e-01  2.930752e+01  -3.708408557 -2.598950e+02
# Fsp3          -2.097556e-02  3.241475e-03   -0.10688949  9.456165e-04  3.043788e-02   0.016058604  4.968763e+00
# MW             2.719911e+00  8.192659e-01   29.30752258  3.043788e-02  1.037952e+02  -3.446238062 -2.740424e+03
# SLA           -5.654123e-01  3.371320e-02   -3.70840856  1.605860e-02 -3.446238e+00   1.844104999  1.120813e+02
# LeafArea      -2.827712e+02 -8.007950e+01 -259.89503140  4.968763e+00 -2.740424e+03 112.081279321  2.715873e+07
# LeafThick      3.664556e-03 -1.921571e-04    0.02557271 -9.457283e-05  2.826367e-02  -0.006765656 -7.692860e-01
# TwigBarkThick -1.131145e-02  4.340104e-04   -0.09480847  2.386330e-04 -1.273920e-01   0.019062563 -1.094200e+01
# TwigSpecDens  -5.814591e-02  2.852046e-03   -0.43137888  1.363493e-03 -5.311015e-01   0.101666130 -4.823629e+01
                  # LeafThick TwigBarkThick  TwigSpecDens
# nAtomP         3.664556e-03 -1.131145e-02  -0.058145910
# ALogP         -1.921571e-04  4.340104e-04   0.002852046
# TopoPSA        2.557271e-02 -9.480847e-02  -0.431378876
# Fsp3          -9.457283e-05  2.386330e-04   0.001363493
# MW             2.826367e-02 -1.273920e-01  -0.531101473
# SLA           -6.765656e-03  1.906256e-02   0.101666130
# LeafArea      -7.692860e-01 -1.094200e+01 -48.236292176
# LeafThick      4.838582e-05 -1.664693e-04  -0.000831091
# TwigBarkThick -1.664693e-04  1.055697e-03   0.003420497
# TwigSpecDens  -8.310910e-04  3.420497e-03   0.017936567


## Create empty matrices for recording phylogenetic correlations, t-values, and p-values
phycor.props = as.data.frame(matrix(1,nrow = 10, ncol = 10))
names(phycor.props) = c("nAtomP", "ALogP", "TopoPSA",  "Fsp3", "MW", "SLA", "LeafArea", "LeafThick", "TwigBarkThick", "TwigSpecDens")
row.names(phycor.props) = c("nAtomP", "ALogP", "TopoPSA",  "Fsp3", "MW", "SLA", "LeafArea", "LeafThick", "TwigBarkThick", "TwigSpecDens")
phycor.props.t = phycor.props.p = phycor.props

## calculate phylogenetic correlations, t-values, and p-values
for (i in 1:nrow(phycor.props)){
	for(j in 2:ncol(phycor.props)){
		r.xy<-cov2cor(obj$R)[i,j]
		t.xy<-r.xy*sqrt((Ntip(phylo.props)-2)/(1-r.xy^2))
		P.xy<-2*pt(abs(t.xy),df=Ntip(phylo.props)-2,lower.tail=F)
		phycor.props[i,j] = phycor.props[j,i] = r.xy
		phycor.props.t[i,j] = phycor.props.t[j,i] = t.xy
		phycor.props.p[i,j] = phycor.props.p[j,i] = P.xy
	}
}

## Plot Figure 1 correlation plot
corrplot(corr = as.matrix(phycor.props))


## Correlations between metabolome chemical properties
tree = phylo.props

## Correlation between nAtomP and TopoPSA
r.xy<-cov2cor(obj$R)["nAtomP","TopoPSA"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.7188854
P.xy
# [1] 3.513331e-05

## Visualize the correlation between nAtomP and TopoPSA
phylomorphospace(phylo.props,sppmeanprops.phylo[,c("nAtomP","TopoPSA")],node.size=c(0,0))
points(sppmeanprops.phylo[,c("nAtomP","TopoPSA")],pch=21,cex=0.7,bg="grey")
text(3,165,paste("evolutionary correlation\nr = ",round(r.xy,4),
    ", P = ",round(P.xy,4),sep="" ),pos=4)


## Correlation between nAtomP and TopoPSA
r.xy<-cov2cor(obj$R)["nAtomP","TopoPSA"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.7188854
P.xy
# [1] 3.513331e-05

## Visualize the correlation between nAtomP and TopoPSA
phylomorphospace(phylo.props,sppmeanprops.phylo[,c("nAtomP","TopoPSA")],node.size=c(0,0))
points(sppmeanprops.phylo[,c("nAtomP","TopoPSA")],pch=21,cex=0.7,bg="grey")
text(3,165,paste("evolutionary correlation\nr = ",round(r.xy,4),
    ", P < 0.0001",sep="" ),pos=4)

## Correlation between nAtomP and MW
r.xy<-cov2cor(obj$R)["nAtomP","MW"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.3233394
P.xy
# [1] 0.107125

## Correlation between nAtomP and ALogP
r.xy<-cov2cor(obj$R)["nAtomP","ALogP"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] -0.3179335
P.xy
# [1] 0.1134659


## Correlation between nAtomP and Fsp3
r.xy<-cov2cor(obj$R)["nAtomP","Fsp3"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] -0.8261305
P.xy
# [1] 2.013028e-07


## Correlation between TopoPSA and MW
r.xy<-cov2cor(obj$R)["TopoPSA","MW"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.4711884
P.xy
# [1] 0.0151092

## Correlation between TopoPSA and ALogP
r.xy<-cov2cor(obj$R)["TopoPSA","ALogP"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] -0.5017726
P.xy
# [1] 0.009007835

## Correlation between TopoPSA and Fsp3
r.xy<-cov2cor(obj$R)["TopoPSA","Fsp3"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] -0.5693531
P.xy
# [1] 0.002399307


## Correlation between MW and ALogP
r.xy<-cov2cor(obj$R)["MW","ALogP"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.3919374
P.xy
# [1] 0.04767533

## Correlation between MW and Fsp3
r.xy<-cov2cor(obj$R)["MW","Fsp3"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.09715566
P.xy
# [1] 0.6368169

## Correlation between ALogP and Fsp3
r.xy<-cov2cor(obj$R)["ALogP","Fsp3"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.5137665
P.xy
# [1] 0.007260083





## Correlations between morphological traits


## Correlation between SLA and LeafArea
r.xy<-cov2cor(obj$R)["SLA","LeafArea"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.01583746
P.xy
# [1] 0.9387918


## Correlation between SLA and LeafThick
r.xy<-cov2cor(obj$R)["SLA","LeafThick"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] -0.716239
P.xy
# [1] 3.871827e-05

## Correlation between SLA and TwigBarkThick
r.xy<-cov2cor(obj$R)["SLA","TwigBarkThick"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.4320348
P.xy
# [1] 0.02751911

## Correlation between SLA and TwigSpecDens
r.xy<-cov2cor(obj$R)["SLA","TwigSpecDens"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] 0.5590026
P.xy
# [1] 0.002991626



## Correlation between LeafArea and LeafThick
r.xy<-cov2cor(obj$R)["LeafArea","LeafThick"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] -0.0212214
P.xy
# [1] 0.9180441

## Correlation between LeafArea and LeafThick
r.xy<-cov2cor(obj$R)["LeafArea","TwigBarkThick"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] -0.06462082
P.xy
# [1] 0.7538055

## Correlation between LeafArea and LeafThick
r.xy<-cov2cor(obj$R)["LeafArea","TwigSpecDens"]
## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
r.xy
# [1] -0.06911138
P.xy
# [1] 0.7372704





##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
#### For Fig. 2a: Forest plot-level correlation plot of chemical properties and morphological traits
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

## Read DatasetS2 and use the forest plot names as row names
metaplotsprops = read.csv("~/Documents/Madidi_Project/Chadwick_ProcB_10.1098:rspb.2025.1721_DatasetS2.csv")
row.names(metaplotsprops) = metaplotsprops$Plot

## Simplify long variable names
names(metaplotsprops)[which(names(metaplotsprops) == "TwigBarkThickness_Relative")] = "TwigBarkThick"
names(metaplotsprops)[which(names(metaplotsprops) == "LeafThickness")] = "LeafThick"
row.names(metaplotsprops) = gsub("PP_", "", row.names(metaplotsprops))

## Calculate correlations among forest plot-level variables
plotcorr = metaplotsprops[,c("Elevation", "Climate_PC1", "nAtomP", "ALogP", "TopoPSA", "Fsp3", "MW", "SLA", "LeafArea", "LeafThick", "TwigBarkThick", "TwigSpecDens")]
plotcorrm = cor(plotcorr)

## Plot correlation plot for Fig. 2a
corrplot(corr = plotcorrm)



##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
#### For Fig. 2b: Redundancy analysis (RDA) of trait and environmental variation among forest plots
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

## Here we control for the effect of phylogenetic similarity in species composition of the plots. 
## We calculate Phylogenetic Community Similarity (PCD) using the SQ function from 'adiv' and represent it as a distance matrix.
## We use this as a covariate in RDA of the relationship between community plant traits, climate, and elevation

## Requires sppmeanprops.phylo from Fig 1 section
dim(sppmeanprops.phylo)
[1] 437  10 # 437 species

## Create an empty matrix in which to record the abundances of each of 437 tree species (columns) in each of 16 forest plots (rows)
madcomm = as.data.frame(matrix(0,nrow = 16, ncol = 437))
row.names(madcomm) = metaplotsprops$Plot
names(madcomm) = gsub("_", " ", row.names(sppmeanprops.phylo))

## Record species abundances in the forest plots to create a community composition matrix
for(i in 1:nrow(madcomm)){
	madplot = droplevels(metaprops[which(metaprops$Plot == row.names(madcomm)[i]),])
	for(j in 1:ncol(madcomm)){
		spp = names(madcomm)[j]
		if(spp %in% madplot$Species_binomial){
			madcomm[i,j] = madplot$Abundance[which(madplot$Species_binomial == spp)]
		}		
	}
}

## Remove any columns representing species not found in the 16 forest plots
colSums(madcomm)[which(colSums(madcomm) < 1)]
madcomm = madcomm[which(colSums(madcomm) > 0)]
names(madcomm) = gsub(" ", "_", names(madcomm))

## Create a phylogenetic distance matrix (cophenetic) for all the tree species
coph.props = as.data.frame(1-(cophenetic(phylo.props)/max(cophenetic(phylo.props))))
coph.props = coph.props[which(row.names(coph.props) %in% names(madcomm)),which(names(coph.props) %in% names(madcomm))]

## Calculate the phylogenetic similarity of every pair of 16 forest plots using dsimcom from the 'adiv' package
madsq = SQ(madcomm, Sigma = coph.props, type = "similarity")
madsimcom = dsimcom(as.matrix(madcomm), Sigma = as.matrix(coph.props), method = 3, option = "relative")
madpcoa <- cmdscale (1-madsimcom, eig = TRUE)

# ordiplot(madpcoa, display = 'sites', type = 'text')

## Simplify long variable names and extract matrix that only contains trait data
names(metaplotsprops)[which(names(metaplotsprops) == "TwigBarkThickness_Relative")] = "TwigBarkThick"
names(metaplotsprops)[which(names(metaplotsprops) == "LeafThickness")] = "LeafThick"

plotcorr = metaplotsprops[,c("Elevation", "Climate_PC1", "nAtomP", "TopoPSA", "MW", "ALogP", "Fsp3", "SLA", "LeafArea", "LeafThick", "TwigBarkThick", "TwigSpecDens")]

madsimcom = as.data.frame(madsimcom)
names(madsimcom) = gsub("PP_", "", names(madsimcom))
row.names(madsimcom) = gsub("PP_", "", row.names(madsimcom))

## Add positions on two phylogenetic similarity axes for each of the 16 Madidi forest plots
metaplotsprops$PhyloPCoA1 = madpcoa$points[,1]
metaplotsprops$PhyloPCoA2 = madpcoa$points[,2]
row.names(metaplotsprops) = metaplotsprops$Plot

## Select only the columns representing the community weighted-mean chemical properties from the data matrix
traitsforrdachemsim = metaplotsprops[,c("nAtomP", "TopoPSA", "MW", "ALogP", "Fsp3", "SLA", "LeafArea", "LeafThick", "TwigSpecDens", "TwigBarkThick")]

## Try an RDA with Elevation, Climate, and two axes of phylogenetic similarity
cdk.elevclimphy = rda(scale(traitsforrdachemsim) ~ Elevation + Climate_PC1 + PhyloPCoA1 + PhyloPCoA2, madplotsprops)
plot(cdk.elevclimphy)

cdk.elevclimphy
# Call: rda(formula = scale(traitsforrdachemsim) ~ Elevation + Climate_PC1 + PhyloPCoA1 + PhyloPCoA2, data =
# madplotsprops)

# -- Model Summary --

              # Inertia Proportion Rank
# Total         10.0000     1.0000     
# Constrained    3.6912     0.3691    4
# Unconstrained  6.3088     0.6309   10

# Inertia is variance

# -- Eigenvalues --

# Eigenvalues for constrained axes:
  # RDA1   RDA2   RDA3   RDA4 
# 2.6119 0.8225 0.1894 0.0674 

# Eigenvalues for unconstrained axes:
   # PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10 
# 2.8432 1.2425 0.7901 0.6787 0.3657 0.2902 0.0832 0.0084 0.0051 0.0017 


## Evaluate which variables contribute significantly to variation among the forest plots
cdk.mod0 = rda(scale(traitsforrdachemsim) ~ 1, madplotsprops)
cdk.elevclim12phy12div = rda(scale(traitsforrdachemsim) ~ Elevation + invSimpson + Climate_PC1 + Climate_PC2 + PhyloPCoA1 + PhyloPCoA2, madplotsprops)
ordiR2step(cdk.mod0, cdk.elevclim12phy12div)

# Step: R2.adj= 0 
# Call: scale(traitsforrdachemsim) ~ 1 
 
                  # R2.adjusted
# + Elevation      0.1989150233
# + PhyloPCoA1     0.1825948033
# + Climate_PC1    0.1465954716
# <All variables>  0.1407054095
# + invSimpson     0.0505244460
# + Climate_PC2    0.0305078189
# <none>           0.0000000000
# + PhyloPCoA2    -0.0005325269

# Call: rda(formula = scale(traitsforrdachemsim) ~ 1, data = madplotsprops)

# -- Model Summary --

              # Inertia Rank
# Total              10     
# Unconstrained      10   10

# Inertia is variance

# -- Eigenvalues --

# Eigenvalues for unconstrained axes:
  # PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8   PC9  PC10 
# 4.730 2.357 0.890 0.813 0.575 0.331 0.221 0.064 0.012 0.008 

## Conclusion: the best model includes Elevation, PhyloPCoA1, and ClimatePC1

## Perform the final RDA with climatic variables as explanatory variables and the significant phylogenetic similarity axes as covariates 
## This accounts for the effect of phylogenetic similarity of the forest plots

cdk.elevclim1phy1 = rda(scale(traitsforrdachemsim) ~ Elevation + Climate_PC1, Z = PhyloPCoA1, madplotsprops)

cdk.elevclim1phy1
# Call: rda(formula = scale(traitsforrdachemsim) ~ Elevation + Climate_PC1, data = madplotsprops, Z =
# PhyloPCoA1)

# -- Model Summary --

              # Inertia Proportion Rank
# Total         10.0000     1.0000     
# Constrained    3.3545     0.3354    2
# Unconstrained  6.6455     0.6646   10

# Inertia is variance

# -- Eigenvalues --

# Eigenvalues for constrained axes:
  # RDA1   RDA2 
# 2.5341 0.8204 

# Eigenvalues for unconstrained axes:
   # PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10 
# 2.8931 1.3498 0.8398 0.6797 0.3960 0.3264 0.1225 0.0217 0.0095 0.0069 

## Plot Figure 2b: RDA of Elevation, Climate_PC1 and traits using phylogenetic similarity as a conditioning matrix
plot(cdk.elevclim1phy1)

#Note: 
## high-elevation forests: 						Kanupa_44, Chaqui_32, Titiri_42
## low-elevation dry forests: 					Pintat_5, Resina_12, Yarimi_9
## low-elevation species-rich wet forests:		Lomasa_39, Lomaka_40, Sumpul_34, Tintay_25, Tintay_24


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
#### For Fig. 3: Linear regression of community weighted mean chemical properties: nAtomP, ALogP, TopoPSA, Fsp3, and MW
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

## Linear regressions of community-weighted mean metabolome chemical properties, elevation, and climate PC1


mod.nAtomP.elev = lm(metaplotsprops$nAtomP ~ metaplotsprops$Elevation)
summary(mod.nAtomP.elev)
# Call:
# lm(formula = madplotsprops$nAtomP ~ madplotsprops$elevation)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -2.68842 -0.60311  0.04005  0.63768  2.71097 

# Coefficients:
                         # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             1.022e+01  8.030e-01  12.726 4.38e-09 ***
# madplotsprops$elevation 8.104e-04  3.902e-04   2.077   0.0567 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.28 on 14 degrees of freedom
# Multiple R-squared:  0.2356,	Adjusted R-squared:  0.1809 
# F-statistic: 4.314 on 1 and 14 DF,  p-value: 0.05669

plot(metaplotsprops$nAtomP ~ madplotsprops$Elevation)
abline(mod.nAtomP.elev)

mod.nAtomP.clim1 = lm(metaplotsprops$nAtomP ~ metaplotsprops$Climate_PC1)
summary(mod.nAtomP.clim1)
# Call:
# lm(formula = metaplotsprops$nAtomP ~ metaplotsprops$Climate_PC1)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -2.7216 -0.5829  0.1030  0.4242  3.0081 

# Coefficients:
                           # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 11.7494     0.3323   35.36 4.29e-15 ***
# metaplotsprops$Climate_PC1  -0.3518     0.2033   -1.73    0.106    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.329 on 14 degrees of freedom
# Multiple R-squared:  0.1762,	Adjusted R-squared:  0.1173 
# F-statistic: 2.994 on 1 and 14 DF,  p-value: 0.1056

mod.ALogP.elev = lm(metaplotsprops$ALogP ~ metaplotsprops$Elevation)
summary(mod.ALogP.elev)
# Call:
# lm(formula = madplotsprops$ALogP ~ madplotsprops$elevation)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -0.38823 -0.22238 -0.02858  0.11253  0.70511 

# Coefficients:
                          # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              3.532e+00  1.800e-01  19.630 1.39e-11 ***
# madplotsprops$elevation -1.912e-04  8.744e-05  -2.187   0.0462 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.2869 on 14 degrees of freedom
# Multiple R-squared:  0.2547,	Adjusted R-squared:  0.2015 
# F-statistic: 4.784 on 1 and 14 DF,  p-value: 0.04619

plot(madplotsprops$ALogP ~ madplotsprops$Elevation)
abline(mod.ALogP.elev)

mod.ALogP.clim1 = lm(metaplotsprops$ALogP ~ metaplotsprops$Climate_PC1)
summary(mod.ALogP.clim1)
# Call:
# lm(formula = metaplotsprops$ALogP ~ metaplotsprops$Climate_PC1)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -0.34065 -0.18433 -0.06396  0.07815  0.61538 

# Coefficients:
                           # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 3.17137    0.06635  47.795   <2e-16 ***
# metaplotsprops$Climate_PC1  0.11447    0.04060   2.819   0.0137 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.2654 on 14 degrees of freedom
# Multiple R-squared:  0.3621,	Adjusted R-squared:  0.3165 
# F-statistic: 7.947 on 1 and 14 DF,  p-value: 0.01366

plot(metaplotsprops$ALogP ~ metaplotsprops$Climate_PC1)
abline(mod.ALogP.clim1)


mod.TopoPSA.elev = lm(metaplotsprops$TopoPSA ~ metaplotsprops$Elevation)
summary(mod.TopoPSA.elev)
# Call:
# lm(formula = madplotsprops$TopoPSA ~ madplotsprops$elevation)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -16.247  -2.673   0.738   4.369  15.935 

# Coefficients:
                         # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             90.169254   5.266095  17.123 8.73e-11 ***
# madplotsprops$elevation  0.007376   0.002559   2.883    0.012 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 8.395 on 14 degrees of freedom
# Multiple R-squared:  0.3725,	Adjusted R-squared:  0.3277 
# F-statistic:  8.31 on 1 and 14 DF,  p-value: 0.01204

plot(metaplotsprops$TopoPSA ~ metaplotsprops$Elevation)
abline(mod.TopoPSA.elev)

mod.TopoPSA.clim1 = lm(metaplotsprops$TopoPSA ~ metaplotsprops$Climate_PC1)
summary(mod.TopoPSA.clim1)
# Call:
# lm(formula = metaplotsprops$TopoPSA ~ metaplotsprops$Climate_PC1)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -18.1029  -1.7349   0.5621   2.6429  17.9390 

# Coefficients:
                           # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 104.092      2.040  51.032   <2e-16 ***
# metaplotsprops$Climate_PC1   -3.872      1.248  -3.102   0.0078 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 8.159 on 14 degrees of freedom
# Multiple R-squared:  0.4073,	Adjusted R-squared:  0.365 
# F-statistic: 9.622 on 1 and 14 DF,  p-value: 0.007802

plot(metaplotsprops$TopoPSA ~ metaplotsprops$Climate_PC1)
abline(mod.TopoPSA.clim1)

mod.Fsp3.elev = lm(metaplotsprops$Fsp3  ~ metaplotsprops$Elevation)
summary(mod.Fsp3.elev)
# Call:
# lm(formula = madplotsprops$Fsp3 ~ madplotsprops$elevation)

# Residuals:
      # Min        1Q    Median        3Q       Max 
# -0.113024 -0.026957  0.008663  0.025174  0.112317 

# Coefficients:
                          # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              5.035e-01  3.589e-02   14.03 1.23e-09 ***
# madplotsprops$elevation -3.923e-05  1.744e-05   -2.25   0.0411 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.05721 on 14 degrees of freedom
# Multiple R-squared:  0.2655,	Adjusted R-squared:  0.2131 
# F-statistic: 5.062 on 1 and 14 DF,  p-value: 0.04106

mod.Fsp3.clim1 = lm(metaplotsprops$Fsp3  ~ metaplotsprops$Climate_PC1)
summary(mod.Fsp3.clim1)
# Call:
# lm(formula = madplotsprops$Fsp3 ~ madplotsprops$PCA1)

# Residuals:
      # Min        1Q    Median        3Q       Max 
# -0.128123 -0.028334  0.008627  0.026018  0.112333 

# Coefficients:
                   # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.429470   0.015086  28.467  8.6e-14 ***
# madplotsprops$PCA1 0.016346   0.009232   1.771   0.0984 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.06035 on 14 degrees of freedom
# Multiple R-squared:  0.183,	Adjusted R-squared:  0.1246 
# F-statistic: 3.135 on 1 and 14 DF,  p-value: 0.09839

plot(madplotsprops$Fsp3  ~ metaplotsprops$Climate_PC1)
abline(mod.Fsp3.clim1)

mod.MW.elev = lm(metaplotsprops$MW ~ metaplotsprops$Elevation)
summary(mod.MW.elev)
# Call:
# lm(formula = madplotsprops$MW ~ madplotsprops$elevation)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -16.6955  -3.0921  -0.1878   1.0582  15.9055 

# Coefficients:
                         # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             3.716e+02  5.544e+00  67.027   <2e-16 ***
# madplotsprops$elevation 2.546e-03  2.694e-03   0.945    0.361    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 8.839 on 14 degrees of freedom
# Multiple R-squared:  0.05996,	Adjusted R-squared:  -0.007185 
# F-statistic: 0.893 on 1 and 14 DF,  p-value: 0.3607

mod.MW.clim1 = lm(metaplotsprops$MW ~ metaplotsprops$Climate_PC1)
summary(mod.MW.clim1)
# Call:
# lm(formula = metaplotsprops$MW ~ metaplotsprops$Climate_PC1)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -16.348  -3.099  -1.268   1.473  17.631 

# Coefficients:
                           # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                376.4290     2.2654 166.164   <2e-16 ***
# metaplotsprops$Climate_PC1  -0.5727     1.3863  -0.413    0.686    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 9.062 on 14 degrees of freedom
# Multiple R-squared:  0.01204,	Adjusted R-squared:  -0.05852 
# F-statistic: 0.1707 on 1 and 14 DF,  p-value: 0.6858




## create Figure 3

plotdata = metaplotsprops

colvec = rep("black", 16)
# colvec[which(plotdata$ForestType == "dry forest")] = "red"

shpvec = rep(16, 16)
# shpvec[which(plotdata$ForestType == "dry forest")] = 1

par(mfrow = c(5, 5), oma = c(2,2,2,2))

par(mar = c(1,0,0,0))

plot.new()

plot(plotdata$nAtomP ~ plotdata$Elevation, pch = shpvec, cex = 2,  xaxt = 'n', ylab = "Aromaticity/Light Absorption (nAtomP)")
segments(x0 = min(plotdata$Elevation), y0 = min(plotdata$Elevation)*(8.104e-04) + 1.022e+01, x1 = max(plotdata$Elevation), y1 = max(plotdata$Elevation)*(8.104e-04) + 1.022e+01, lwd = 2, lty = 2)
text(x = 724, y = 15, labels = "a")
# text(x = 1400, y = 14.5, labels = c(expression(italic("R")^2), "= 0.18"))
# text(x = 1400, y = 14, labels = c(expression("p"), "= 0.056")
text(x = 1400, y = 14.5, labels = "R2 = 0.18")
text(x = 1400, y = 14, labels = "p = 0.056")

plot(plotdata$nAtomP ~ plotdata$Climate_PC1, pch = shpvec,, cex = 2, yaxt = 'n', xaxt = 'n')
segments(x0 = min(plotdata$Climate_PC1), y0 = min(plotdata$Climate_PC1)*(0.012575) + 0.215372, x1 = max(plotdata$Climate_PC1), y1 = max(plotdata$Climate_PC1)*(0.012575) + 0.215372, lwd = 2, lty = 1)
# segments(x0 = min(plotdata$PCA1), y0 = min(plotdata$PCA1)*(0.012900) + 0.217415, x1 = max(plotdata$PCA1), y1 = max(plotdata$PCA1)*(0.012900) + 0.217415, lwd = 2, lty = 2)
text(x = -3.5, y = 15, labels = "b")
# text(x = 1.4, y = 0.18, labels = "R2 = 0.45")
# text(x = 1.4, y = 0.17, labels = "p < 0.01")
# text(x = 2500, y = 0.55, labels = "R2 = 0.48")
# text(x = 2500, y = 0.50, labels = "p < 0.01")

# plot(plotdata$nAtomP ~ plotdata$PCA2, pch = shpvec, cex = 2, yaxt = 'n', xaxt = 'n')
# text(x = -1.85, y = 15, labels = "d")

# plot(plotdata$nAtomP ~ plotdata$invSimpson, pch = shpvec, cex = 2, yaxt = 'n', xaxt = 'n')
# text(x = 3.5, y = 15, labels = "c")

plot.new()
plot.new()
plot.new()

plot(plotdata$ALogP ~ plotdata$Elevation, pch = shpvec, cex = 2, xaxt = 'n', ylab = "Topological Polar Surface Area")
segments(x0 = min(plotdata$Elevation), y0 = min(plotdata$Elevation)*(-1.912e-04) + 3.532e+00, x1 = max(plotdata$Elevation), y1 = max(plotdata$Elevation)*(-1.912e-04) + 3.532e+00, lwd = 2, lty = 1)
text(x = 724, y = 3.86, labels = "c")
text(x = 2500, y = 3.86, labels = "R2 = 0.20")
text(x = 2500, y = 3.6, labels = "p < 0.05")

plot(plotdata$ALogP ~ plotdata$Climate_PC1, pch = shpvec, cex = 2,  yaxt = 'n', xaxt = 'n')
segments(x0 = min(plotdata$Climate_PC1), y0 = min(plotdata$Climate_PC1)*(0.11447) + 3.17137, x1 = max(plotdata$Climate_PC1), y1 = max(plotdata$PCA1)*(0.11447) + 3.17137, lwd = 2, lty = 1)
text(x = -3.5, y = 3.86, labels = "d")
text(x = -1, y = 3.86, labels = "R2 = 0.32")
text(x = -1, y = 3.6, labels = "p = 0.01")

# plot(plotdata$ALogP ~ plotdata$PCA2, pch = shpvec, cex = 2, yaxt = 'n', xaxt = 'n')
# text(x = -1.85, y = 3.86, labels = "p")

# plot(plotdata$ALogP ~ plotdata$invSimpson, pch = shpvec,, cex = 2, xaxt = 'n', yaxt = 'n')
# # segments(x0 = min(plotdata$invSimpson), y0 = min(plotdata$invSimpson)*(0.0022101) + 0.1833859, x1 = max(plotdata$invSimpson), y1 = max(plotdata$invSimpson)*(0.0022101) + 0.1833859, lwd = 2, lty = 1)
# text(x = 3.5, y = 3.86, labels = "l")

plot.new()
plot.new()
plot.new()


plot(plotdata$TopoPSA ~ plotdata$Elevation, pch = shpvec, cex = 2, xaxt = 'n', ylab = "Topological Polar Surface Area")
segments(x0 = min(plotdata$Elevation), y0 = min(plotdata$Elevation)*(0.007376) + 90.169254, x1 = max(plotdata$Elevation), y1 = max(plotdata$Elevation)*(0.007376) + 90.169254, lwd = 2, lty = 1)
text(x = 724, y = 124, labels = "e")
# text(x = 1600, y = 124, labels = expression(italic("R")^2 = 0.33))
# text(x = 1600, y = 120, labels = expression(italic("p") = 0.01))
text(x = 1600, y = 124, labels = "R2 = 0.33")
text(x = 1600, y = 120, labels = "p = 0.01")

plot(plotdata$TopoPSA ~ plotdata$Climate_PC1, pch = shpvec, cex = 2,  yaxt = 'n', xaxt = 'n')
segments(x0 = min(plotdata$Climate_PC1), y0 = min(plotdata$Climate_PC1)*(-3.872) + 104.092, x1 = max(plotdata$Climate_PC1), y1 = max(plotdata$Climate_PC1)*(-3.872) + 104.092, lwd = 2, lty = 1)
text(x = -3.5, y = 124, labels = "f")
text(x = 1.4, y = 124, labels = "R2 = 0.37")
text(x = 1.4, y = 120, labels = "p < 0.01")

# plot(plotdata$TopoPSA ~ plotdata$PCA2, pch = shpvec, cex = 2, yaxt = 'n', xaxt = 'n')
# text(x = -1.85, y = 124, labels = "h")

# plot(plotdata$TopoPSA ~ plotdata$invSimpson, pch = shpvec, cex = 2, yaxt = 'n',, xaxt = 'n')
# # segments(x0 = min(plotdata$invSimpson), y0 = min(plotdata$invSimpson)*(0.0022101) + 0.1833859, x1 = max(plotdata$invSimpson), y1 = max(plotdata$invSimpson)*(0.0022101) + 0.1833859, lwd = 2, lty = 1)
# text(x = 3.5, y = 124, labels = "f")

plot.new()
plot.new()
plot.new()

plot(plotdata$Fsp3 ~ plotdata$Elevation, pch = shpvec, cex = 2, ylab = "Carbon Fsp3 Bond Saturation", xaxt = 'n')
segments(x0 = min(plotdata$Elevation), y0 = min(plotdata$Elevation)*(-3.923e-05) + 5.035e-01, x1 = max(plotdata$Elevation), y1 = max(plotdata$Elevation)*(-3.923e-05) + 5.035e-01, lwd = 2, lty = 1)
text(x = 724, y = 0.51, labels = "g")
text(x = 1600, y = 0.33, labels = "R2 = 0.21")
text(x = 1600, y = 0.30, labels = "p = 0.04")

plot(plotdata$Fsp3 ~ plotdata$Climate_PC1, pch = shpvec, cex = 2,  yaxt = 'n', xaxt = 'n')
# segments(x0 = min(plotdata$PCA1), y0 = min(plotdata$PCA1)*(0.016346) + 0.429470, x1 = max(plotdata$PCA1), y1 = max(plotdata$PCA1)*(0.016346) + 0.429470, lwd = 2, lty = 1)
text(x = -3.5, y = 0.51, labels = "h")
# text(x = -1, y = 0.51, labels = "R2 = 0.12")
# text(x = -1, y = 0.49, labels = "p = 0.10")

# plot(plotdata$Fsp3 ~ plotdata$PCA2, pch = shpvec, cex = 2, yaxt = 'n', xaxt = 'n')
# text(x = -1.85, y = 0.51, labels = "t")

# plot(plotdata$Fsp3 ~ plotdata$invSimpson, pch = shpvec,, cex = 2, xlab = "Diversity (Inverse Simpson)", yaxt = 'n')
# # segments(x0 = min(plotdata$invSimpson), y0 = min(plotdata$invSimpson)*(0.0022101) + 0.1833859, x1 = max(plotdata$invSimpson), y1 = max(plotdata$invSimpson)*(0.0022101) + 0.1833859, lwd = 2, lty = 1)
# text(x = 3.5, y = 0.51, labels = "o")


plot.new()
plot.new()
plot.new()

plot(plotdata$MW ~ plotdata$Elevation, pch =shpvec, cex = 2, ylab = "Molecular Weight", xlab = "Elevation (m)")
text(x = 724, y = 395, labels = "i")

plot(plotdata$MW ~ plotdata$Climate_PC1, pch = shpvec, cex = 2, yaxt = 'n', xlab = "Climate PC1")
text(x = -3.5, y = 395, labels = "j")

# plot(plotdata$MW ~ plotdata$PCA2, pch = shpvec, cex = 2, yaxt = 'n', xlab = "Climate PC2")
# text(x = -1.8, y = 395, labels = "l")
# # text(x = -1.2, y = 0.26, labels = "R2 = 0.24")
# # text(x = -1.2, y = 0.24, labels = "p = 0.03")

# plot(plotdata$MW ~ plotdata$invSimpson, pch = shpvec, cex = 2, xaxt = 'n', yaxt = 'n')
# text(x = 3.5, y = 395, labels = "i")


## end Fig. 3



##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
#### For Fig. 4: Linear regressions of chemical similarity for upper/lower quartile nAtomP, ALogP, TopoPSA, Fsp3 vs elevation and climate
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

## Calculate upper and lower quartiles of metabolties for CDK varaibles representing metabolite chemical properties
nAtomPquant = quantile(mastertable$nAtomP, probs = c(0.25,0.75), na.rm = T)
TopoPSAquant = quantile(mastertable$TopoPSA, probs = c(0.25,0.75), na.rm = T)
ALogPquant = quantile(mastertable$ALogP, probs = c(0.25,0.75), na.rm = T)
Fsp3quant = quantile(mastertable$Fsp3, probs = c(0.25,0.75), na.rm = T)

## Define input variables for the calcCSCS function
spp = FALSE
date = 20251104
class = "total"
proj = "madidi_cscs_"

## Define 'calcCSCS' function to calculate the chemical structural-compositional similarity (CSCS) metric
## following Sedio et al. 2017 Ecology https://doi.org/10.1002/ecy.1689
calcCSCS = function(spp = "TRUE", class = "total", date = 20240814, proj = "MyProjectCSCS"){
  outfile = paste(proj, class, date, ".RData",sep="")
	if(class == "nAtomP_75"){
  		heat.class = mastertable[which(mastertable$nAtomP > nAtomPquant[2]),42:ncol(mastertable)]
	}
	if(class == "ALogP_25"){
		heat.class = mastertable[which(mastertable$ALogP < ALogPquant[1]),42:ncol(mastertable)]
	}
	if(class == "TopoPSA_75"){
		heat.class = mastertable[which(mastertable$TopoPSA > TopoPSAquant[2]),42:ncol(mastertable)]
	}
	if(class == "TopoPSA_25"){
  		heat.class = mastertable[which(mastertable$TopoPSA < TopoPSAquant[1]),42:ncol(mastertable)]
	}
	if(class == "Fsp3_75"){
		heat.class = mastertable[which(mastertable$Fsp3 > Fsp3quant[2]),42:ncol(mastertable)]
	}	
	sampsByCompounds = as.data.frame(t(as.data.frame(heat.class)))
	names(sampsByCompounds) = row.names(heat.class)	
	recnumcomps = rep(NA, nrow(sampsByCompounds))
	for(x in 1:length(recnumcomps)){
		recnumcomps[x] = length(which(sampsByCompounds[x,] > 0))
	}
	sampsByCompounds = sampsByCompounds[which(recnumcomps > 1),]
	comm.mat = as.data.frame(matrix(1, nrow=2, ncol = ncol(sampsByCompounds)))
	names(comm.mat) = names(sampsByCompounds)
	# tree = prune.sample(samp = comm.mat, phylo = qtree.itol)
	tree = prune.sample(samp = comm.mat, phylo = chemtree)
	maxdist = max(cophenetic(tree))
	pairwise.comps = 1-(cophenetic(tree)/maxdist)
	net.comps = nrow(pairwise.comps)
	nspp = nrow(sampsByCompounds)
	ncomps = ncol(sampsByCompounds)
	pairwise.spp = as.data.frame(matrix(0,nrow = nspp, ncol = nspp))
	names(pairwise.spp) = row.names(sampsByCompounds)
	row.names(pairwise.spp) = row.names(sampsByCompounds)
	sampsCompsStand = sampsByCompounds
	for(i in 1:nrow(sampsByCompounds)){	
		sampsCompsStand[i,] = sampsByCompounds[i,]/sum(sampsByCompounds[i,])
	}
	diags = pairwise.spp
	for (k in 1:nspp){
		sppX = as.character(row.names(sampsCompsStand)[k])
		cat("Comparing ", sppX, " to itself", "\n", sep = "")
		sppXonly = sampsCompsStand[k,which(sampsCompsStand[k,]>0)]
		ncomps = length(sppXonly)
		pairwise.comps.temp = pairwise.comps[names(sppXonly),names(sppXonly)]
		diags[k,k] = sum(((outer(as.numeric(sppXonly), as.numeric(sppXonly)))*pairwise.comps.temp),na.rm = T)	
	}
	for (i in 1:nspp){
		spp1 = as.character(row.names(sampsCompsStand)[i])
		for (j in i:nspp){
			spp2 = as.character(row.names(sampsCompsStand)[j])
			cat("Comparing ", spp1, " to ", spp2, "\n", sep = "")
			#identify which compounds are in each species
			spp1comps = sampsCompsStand[spp1,]
			spp2comps = sampsCompsStand[spp2,]
			spp_pair = rbind(spp1comps,spp2comps)
			paircomps = spp_pair[,which(colSums(spp_pair)>0)]
			#make a pairwise.comps matrix for only those compounds found in either species in the species pair
			ncomps = ncol(paircomps)
			pairwise.comps.temp = pairwise.comps[names(paircomps),names(paircomps)]
			pairwise.spp[i,j] = pairwise.spp[j,i] = sum(((outer(as.numeric(paircomps[1,]), as.numeric(paircomps[2,])))*pairwise.comps.temp), na.rm = T)/max(diags[i,i], diags[j,j])
		}
	}
	cscs = pairwise.spp
  save(class, date, heat.class, sampsByCompounds, cscs, file = outfile)
  return(cscs)
  cat("Completed cacluation of CSCS for all sample pairs","\n")
}

## Calculate CSCS between all pairs of species with respect to metabolites in the upper or lower quartile of five chemical properties
calcCSCS(spp = F, class = "nAtomP_75", date = 20251104, proj = "madidi_cscs_")
cscs.nAtomP.75 = cscs

calcCSCS(spp = F, class = "ALogP_25", date = 20251104, proj = "madidi_cscs_")
cscs.ALogP.25 = cscs

calcCSCS(spp = F, class = "TopoPSA_75", date = 20251104, proj = "madidi_cscs_")
cscs.TopoPSA.75 = cscs

calcCSCS(spp = F, class = "TopoPSA_25", date = 20251104, proj = "madidi_cscs_")
cscs.TopoPSA.25 = cscs

calcCSCS(spp = F, class = "Fsp3_75", date = 20251104, proj = "madidi_cscs_")
cscs.Fsp3.75 = cscs

## Create empty matrices to record community-weighted mean CSCS for each chemical property
metapropstraits = metaprops
madplotsprops = metaplotsprops
row.names(madplotsprops) = madplotsprops$Plot

nAtomP.75.plots = ALogP.25.plots = TopoPSA.75.plots = TopoPSA.25.plots = Fsp3.75.plots = rep(NA, nrow(madplotsprops))
names(nAtomP.75.plots) = names(ALogP.25.plots) = names(TopoPSA.75.plots) = names(TopoPSA.25.plots) = names(Fsp3.75.plots) = row.names(madplotsprops)
metapropstraits$nAtomP.75 = NA
metapropstraits$ALogP.25 = NA
metapropstraits$TopoPSA.75 = NA
metapropstraits$TopoPSA.25 = NA
metapropstraits$Fsp3.75 = NA
for(i in 1:nrow(madplotsprops)){
	plot = row.names(madplotsprops)[i]
	sprec.nAtomP.75 = rep(NA, length(which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.nAtomP.75)))))
	names(sprec.nAtomP.75) = metaprops$sampleCode[which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.nAtomP.75)))]
	sprec.ALogP.25 = rep(NA, length(which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.ALogP.25)))))
	names(sprec.ALogP.25) = metaprops$sampleCode[which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.ALogP.25)))]
	sprec.TopoPSA.75 = rep(NA, length(which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.TopoPSA.75)))))
	names(sprec.TopoPSA.75) = metaprops$sampleCode[which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.TopoPSA.75)))]
	sprec.TopoPSA.25 = rep(NA, length(which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.TopoPSA.25)))))
	names(sprec.TopoPSA.25) = metaprops$sampleCode[which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.TopoPSA.25)))]
	sprec.Fsp3.75 = rep(NA, length(which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.Fsp3.75)))))
	names(sprec.Fsp3.75) = metaprops$sampleCode[which((metaprops$Plot == plot)&(metaprops$sampleCode %in% names(cscs.Fsp3.75)))]	
	cscs.nAtomP.75.plot = cscs.nAtomP.75[names(sprec.nAtomP.75),names(sprec.nAtomP.75)]
	cscs.ALogP.25.plot = cscs.ALogP.25[names(sprec.ALogP.25),names(sprec.ALogP.25)]
	cscs.TopoPSA.75.plot = cscs.TopoPSA.75[names(sprec.TopoPSA.75),names(sprec.TopoPSA.75)]
	cscs.TopoPSA.25.plot = cscs.TopoPSA.25[names(sprec.TopoPSA.25),names(sprec.TopoPSA.25)]
	cscs.Fsp3.75.plot = cscs.Fsp3.75[names(sprec.Fsp3.75),names(sprec.Fsp3.75)]
	abs.plot = metaprops$Abundance[which(metaprops$Plot == plot)]
	names(abs.plot) = metaprops$sampleCode[which(metaprops$Plot == plot)]
	abs.plot.nAtomP.75 = abs.plot[names(sprec.nAtomP.75)]
	abs.plot.ALogP.25 = abs.plot[names(sprec.ALogP.25)]
	abs.plot.TopoPSA.75 = abs.plot[names(sprec.TopoPSA.75)]
	abs.plot.TopoPSA.25 = abs.plot[names(sprec.TopoPSA.25)]
	abs.plot.Fsp3.75 = abs.plot[names(sprec.Fsp3.75)]
	for(j in 1:length(abs.plot.nAtomP.75)){
		sprec.nAtomP.75[j] = weighted.mean(cscs.nAtomP.75.plot[j,-j], abs.plot.nAtomP.75[-j])
		metapropstraits$nAtomP.75[which(metapropstraits$sampleCode == names(sprec.nAtomP.75)[j])] = sprec.nAtomP.75[j]
	}
	for(j in 1:length(abs.plot.ALogP.25)){
		sprec.ALogP.25[j] = weighted.mean(cscs.ALogP.25.plot[j,-j], abs.plot.ALogP.25[-j])
		metapropstraits$ALogP.25[which(metapropstraits$sampleCode == names(sprec.ALogP.25)[j])] = sprec.ALogP.25[j]
	}
	for(j in 1:length(abs.plot.TopoPSA.75)){
		sprec.TopoPSA.75[j] = weighted.mean(cscs.TopoPSA.75.plot[j,-j], abs.plot.TopoPSA.75[-j])
		metapropstraits$TopoPSA.75[which(metapropstraits$sampleCode == names(sprec.TopoPSA.75)[j])] = sprec.TopoPSA.75[j]
	}
	for(j in 1:length(abs.plot.TopoPSA.25)){
		sprec.TopoPSA.25[j] = weighted.mean(cscs.TopoPSA.25.plot[j,-j], abs.plot.TopoPSA.25[-j])
		metapropstraits$TopoPSA.25[which(metapropstraits$sampleCode == names(sprec.TopoPSA.25)[j])] = sprec.TopoPSA.25[j]
	}
	for(j in 1:length(abs.plot.Fsp3.75)){
		sprec.Fsp3.75[j] = weighted.mean(cscs.Fsp3.75.plot[j,-j], abs.plot.Fsp3.75[-j])
		metapropstraits$Fsp3.75[which(metapropstraits$sampleCode == names(sprec.Fsp3.75)[j])] = sprec.Fsp3.75[j]
	}
	nAtomP.75.plots[i] = weighted.mean(sprec.nAtomP.75, abs.plot.nAtomP.75)
	ALogP.25.plots[i] = weighted.mean(sprec.ALogP.25, abs.plot.ALogP.25)
	TopoPSA.75.plots[i] = weighted.mean(sprec.TopoPSA.75, abs.plot.TopoPSA.75)
	TopoPSA.25.plots[i] = weighted.mean(sprec.TopoPSA.25, abs.plot.TopoPSA.25)
	Fsp3.75.plots[i] = weighted.mean(sprec.Fsp3.75, abs.plot.Fsp3.75)
}

metapropstraits$Elevation = NA
for(i in 1:length(levels(as.factor(metapropstraits$Plot)))){
	plot = levels(as.factor(metapropstraits$Plot))[i]
	metapropstraits$Elevation[which(metapropstraits$Plot == plot)] = madplotsprops$Elevation[which(madplotsprops$Plot == plot)]
}

madplotsprops = cbind(madplotsprops, nAtomP.25.plots, nAtomP.75.plots, TopoPSA.25.plots, TopoPSA.75.plots, MW.25.plots, MW.75.plots, ALogP.25.plots, ALogP.75.plots, Fsp3.25.plots, Fsp3.75.plots, FMF.25.plots, FMF.75.plots)



## H3a: Abiotic stress (radiation, desiccation) select for convergence in radiation- (upper quartile nAtomP) and desiccation-medating metabolites (lower quartile ALogP, upper quartile TopoPSA) at high elevations


mod = lm(madplotsprops$nAtomP.75.plots ~ madplotsprops$Elevation) ##### trees in high-elevation plots are more similar with resp to high nAtomP (i.e. light-absorbing) compounds
summary(mod)
plot(madplotsprops$nAtomP.75.plots ~ madplotsprops$Elevation) 
abline(mod)
# Call:
# lm(formula = madplotsprops$nAtomP.75.plots ~ madplotsprops$elevation)

# Residuals:
      # Min        1Q    Median        3Q       Max 
# -0.106233 -0.021394  0.003177  0.029758  0.069209 

# Coefficients:
                         # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             5.781e-01  2.861e-02  20.206 9.36e-12 ***
# madplotsprops$elevation 6.952e-05  1.390e-05   5.001 0.000194 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.04561 on 14 degrees of freedom
# Multiple R-squared:  0.6411,	Adjusted R-squared:  0.6154 
# F-statistic: 25.01 on 1 and 14 DF,  p-value: 0.0001943

mod = lm(madplotsprops$nAtomP.75.plots ~ madplotsprops$Climate_PC1)  #### trees in colder and drier plots are more similar with resp to high nAtomP (i.e. light-absorbing) compounds
summary(mod)
plot(madplotsprops$nAtomP.75.plots ~ madplotsprops$Climate_PC1)
abline(mod)
# Call:
# lm(formula = madplotsprops$nAtomP.75.plots ~ madplotsprops$PCA1)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -0.14466 -0.02959  0.01289  0.03762  0.09756 

# Coefficients:
                    # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.709364   0.015031  47.193   <2e-16 ***
# madplotsprops$PCA1 -0.026740   0.009198  -2.907   0.0115 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.06012 on 14 degrees of freedom
# Multiple R-squared:  0.3764,	Adjusted R-squared:  0.3319 
# F-statistic: 8.451 on 1 and 14 DF,  p-value: 0.01148


mod = lm(madplotsprops$ALogP.25.plots ~ madplotsprops$Elevation) ##### trees in high-elevation plots are more similar with resp to low ALogP (i.e. less hydrophobic) compounds
summary(mod)
plot(madplotsprops$ALogP.25.plots ~ madplotsprops$Elevation) 
abline(mod)


mod = lm(madplotsprops$ALogP.25.plots ~ madplotsprops$Climate_PC1)  ### trees in colder and drier plots are more similar with resp to low ALogP (i.e. less hydrophobic) compounds
summary(mod)
plot(madplotsprops$ALogP.25.plots ~ madplotsprops$Climate_PC1)
abline(mod)


mod = lm(madplotsprops$TopoPSA.75.plots ~ madplotsprops$Elevation) ##### trees in high-elevation plots are more similar with resp to high TopoPSA (i.e. more polar) compounds
summary(mod)
plot(madplotsprops$TopoPSA.75.plots ~ madplotsprops$Elevation) 
abline(mod)


mod = lm(madplotsprops$TopoPSA.75.plots ~ madplotsprops$Climate_PC1) ### trees in colder and drier plots are more similar with resp to high TopoPSA (i.e. more polar) compounds
summary(mod)
plot(madplotsprops$TopoPSA.75.plots ~ madplotsprops$Climate_PC1)
abline(mod)


## H3b: Biotic interactions select for divergence in defensive metabolites (lower quartile TopoPSA, upper quartile Fsp3) at low elevations

mod = lm(madplotsprops$TopoPSA.25.plots ~ madplotsprops$Elevation) 
summary(mod)
plot(madplotsprops$TopoPSA.25.plots ~ madplotsprops$Elevation) 
abline(mod)


mod = lm(madplotsprops$TopoPSA.25.plots ~ madplotsprops$Climate_PC1)  
summary(mod)
plot(madplotsprops$TopoPSA.25.plots ~ madplotsprops$Climate_PC1)
abline(mod)

mod = lm(madplotsprops$Fsp3.75.plots ~ madplotsprops$Elevation) 
summary(mod)
plot(madplotsprops$Fsp3.75.plots ~ madplotsprops$Elevation) 
abline(mod)


mod = lm(madplotsprops$Fsp3.75.plots ~ madplotsprops$Climate_PC1) 
summary(mod)
plot(madplotsprops$Fsp3.75.plots ~ madplotsprops$Climate_PC1)
abline(mod)




###############
###############
###############
###############
############### 2025-02-12: H3a: Radiation (nAtomP.75), Drought (ALogP.25, TopoPSA.75), H3b: Biotic (TopoPSA.25, Fsp3.75)   ****Fig 4*****
############### 
###############
###############
###############

## Read DatasetS2, with environmental and community data for each of the 16 forest plots, if needed
metaplotsprops = read.csv("~/Documents/Madidi_Project/Chadwick_ProcB_10.1098:rspb.2025.1721_DatasetS2.csv")
row.names(metaplotsprops) = metaplotsprops$Plot

## Define the datapoints to be plotted as black dots (this could be modified to highlight particular forest plots)
shpvec = rep(16, 16)
# shpvec[which(metaplotsprops$ForestType == "dry forest")] = 1

## Create Figure 4
par(mfrow = c(5, 5), oma = c(2,2,2,2))
par(mar = c(1,0,0,0))
plot.new()

plot(metaplotsprops$nAtomP_75 ~ metaplotsprops$Elevation, pch = shpvec, cex = 2,  xaxt = 'n', ylab = "Upper Quartile Aromaticity (nAtomP)", ylim = c(0.40, 0.85))
segments(x0 = min(metaplotsprops$Elevation), y0 = min(metaplotsprops$Elevation)*(6.952e-05) + 5.781e-01, x1 = max(metaplotsprops$Elevation), y1 = max(metaplotsprops$Elevation)*(6.952e-05) + 5.781e-01, lwd = 2, lty = 1)
text(x = 724, y = 0.80, labels = "a")
text(x = 2700, y = 0.60, labels = "R2 = 0.62")
text(x = 2700, y = 0.55, labels = "p < 0.001")

plot(metaplotsprops$nAtomP_75 ~ metaplotsprops$Climate_PC1, pch = shpvec,, cex = 2, yaxt = 'n', xaxt = 'n', ylim = c(0.40, 0.85))
segments(x0 = min(metaplotsprops$Climate_PC1), y0 = min(metaplotsprops$Climate_PC1)*(-0.026740) + 0.709364, x1 = max(metaplotsprops$Climate_PC1), y1 = max(metaplotsprops$Climate_PC1)*(-0.026740) + 0.709364, lwd = 2, lty = 1)
# segments(x0 = min(plotdata$Climate_PC1), y0 = min(plotdata$Climate_PC1)*(0.012900) + 0.217415, x1 = max(plotdata$Climate_PC1), y1 = max(plotdata$Climate_PC1)*(0.012900) + 0.217415, lwd = 2, lty = 2)
text(x = -3.5, y = 0.55, labels = "b")
text(x = -1.9, y = 0.60, labels = "R2 = 0.33")
text(x = -1.9, y = 0.55, labels = "p = 0.011")

plot.new()
plot.new()
plot.new()

plot(metaplotsprops$ALogP_25 ~ metaplotsprops$Elevation, pch = shpvec, cex = 2,  xaxt = 'n', ylab = "Lower Quartile Hydrophobicity (ALogP)", ylim = c(0.40, 0.85))
segments(x0 = min(metaplotsprops$Elevation), y0 = min(metaplotsprops$Elevation)*(9.334e-05) + 4.477e-01, x1 = max(metaplotsprops$Elevation), y1 = max(metaplotsprops$Elevation)*(9.334e-05) + 4.477e-01, lwd = 2, lty = 1)
text(x = 724, y = 0.80, labels = "c")
text(x = 2700, y = 0.60, labels = "R2 = 0.67")
text(x = 2700, y = 0.55, labels = "p < 0.001")

plot(metaplotsprops$ALogP_25 ~ metaplotsprops$Climate_PC1, pch = shpvec,, cex = 2, yaxt = 'n', xaxt = 'n', ylim = c(0.40, 0.85))
segments(x0 = min(metaplotsprops$Climate_PC1), y0 = min(metaplotsprops$Climate_PC1)*(-0.04047) + 0.62378, x1 = max(metaplotsprops$Climate_PC1), y1 = max(metaplotsprops$Climate_PC1)*(-0.04047) + 0.62378, lwd = 2, lty = 1)
# segments(x0 = min(plotdata$Climate_PC1), y0 = min(plotdata$Climate_PC1)*(0.012900) + 0.217415, x1 = max(plotdata$Climate_PC1), y1 = max(plotdata$Climate_PC1)*(0.012900) + 0.217415, lwd = 2, lty = 2)
text(x = -3.5, y = 0.55, labels = "d")
text(x = -1.9, y = 0.60, labels = "R2 = 0.48")
text(x = -1.9, y = 0.55, labels = "p = 0.002")

plot.new()
plot.new()
plot.new()

plot(metaplotsprops$TopoPSA_75 ~ metaplotsprops$Elevation, pch = shpvec, cex = 2,  xaxt = 'n', ylab = "Upper Quartile Polar Surface Area (TopoPSA)", ylim = c(0.40, 0.85))
segments(x0 = min(metaplotsprops$Elevation), y0 = min(metaplotsprops$Elevation)*(5.144e-05) + 5.361e-01, x1 = max(metaplotsprops$Elevation), y1 = max(metaplotsprops$Elevation)*(5.144e-05) + 5.361e-01, lwd = 2, lty = 1)
text(x = 724, y = 0.80, labels = "d")
text(x = 2700, y = 0.60, labels = "R2 = 0.56")
text(x = 2700, y = 0.55, labels = "p < 0.001")

plot(metaplotsprops$TopoPSA_75 ~ metaplotsprops$Climate_PC1, pch = shpvec,, cex = 2, yaxt = 'n', xaxt = 'n', ylim = c(0.40, 0.85))
segments(x0 = min(metaplotsprops$Climate_PC1), y0 = min(metaplotsprops$Climate_PC1)*(-0.019396) + 0.633156, x1 = max(metaplotsprops$Climate_PC1), y1 = max(metaplotsprops$Climate_PC1)*(-0.019396) + 0.633156, lwd = 2, lty = 1)
# segments(x0 = min(plotdata$Climate_PC1), y0 = min(plotdata$Climate_PC1)*(0.012900) + 0.217415, x1 = max(plotdata$Climate_PC1), y1 = max(plotdata$Climate_PC1)*(0.012900) + 0.217415, lwd = 2, lty = 2)
text(x = -3.5, y = 0.55, labels = "e")
text(x = -1.9, y = 0.60, labels = "R2 = 0.28")
text(x = -1.9, y = 0.55, labels = "p = 0.019")

plot.new()
plot.new()
plot.new()

plot(metaplotsprops$TopoPSA_25 ~ metaplotsprops$Elevation, pch = shpvec, cex = 2,  xaxt = 'n', ylab = "Lower Quartile Polar Surface Area (TopoPSA)", ylim = c(0.40, 0.85))
# segments(x0 = min(metaplotsprops$Elevation), y0 = min(metaplotsprops$Elevation)*(6.952e-05) + 5.781e-01, x1 = max(metaplotsprops$Elevation), y1 = max(metaplotsprops$Elevation)*(6.952e-05) + 5.781e-01, lwd = 2, lty = 1)
text(x = 724, y = 0.90, labels = "f")
# text(x = 2700, y = 0.60, labels = "R2 = 0.62")
# text(x = 2700, y = 0.55, labels = "p < 0.001")

plot(metaplotsprops$TopoPSA_25 ~ metaplotsprops$Climate_PC1, pch = shpvec,, cex = 2, yaxt = 'n', xaxt = 'n', ylim = c(0.40, 0.85))
# segments(x0 = min(metaplotsprops$Climate_PC1), y0 = min(metaplotsprops$Climate_PC1)*(-0.026740) + 0.709364, x1 = max(metaplotsprops$Climate_PC1), y1 = max(metaplotsprops$Climate_PC1)*(-0.026740) + 0.709364, lwd = 2, lty = 1)
# segments(x0 = min(plotdata$Climate_PC1), y0 = min(plotdata$Climate_PC1)*(0.012900) + 0.217415, x1 = max(plotdata$Climate_PC1), y1 = max(plotdata$Climate_PC1)*(0.012900) + 0.217415, lwd = 2, lty = 2)
text(x = -3.5, y = 0.90, labels = "g")
# text(x = -1.9, y = 0.60, labels = "R2 = 0.33")
# text(x = -1.9, y = 0.55, labels = "p = 0.011")

plot.new()
plot.new()
plot.new()

plot(metaplotsprops$Fsp3_75 ~ metaplotsprops$Elevation, pch = shpvec, cex = 2,   xlab = "Elevation (m)", ylab = "Upper Quartile 3D Complexity (Fsp3)", ylim = c(0.40, 0.85))
# segments(x0 = min(metaplotsprops$Elevation), y0 = min(metaplotsprops$Elevation)*(5.144e-05) + 5.361e-01, x1 = max(metaplotsprops$Elevation), y1 = max(metaplotsprops$Elevation)*(5.144e-05) + 5.361e-01, lwd = 2, lty = 1)
text(x = 724, y = 0.90, labels = "h")
# text(x = 2700, y = 0.60, labels = "R2 = 0.56")
# text(x = 2700, y = 0.55, labels = "p < 0.001")

plot(metaplotsprops$Fsp3_75 ~ metaplotsprops$Climate_PC1, pch = shpvec,, cex = 2, yaxt = 'n', xlab = "Climate PC1", ylim = c(0.40, 0.85))
# segments(x0 = min(metaplotsprops$Climate_PC1), y0 = min(metaplotsprops$Climate_PC1)*(-0.019396) + 0.633156, x1 = max(metaplotsprops$Climate_PC1), y1 = max(metaplotsprops$Climate_PC1)*(-0.019396) + 0.633156, lwd = 2, lty = 1)
# segments(x0 = min(plotdata$Climate_PC1), y0 = min(plotdata$Climate_PC1)*(0.012900) + 0.217415, x1 = max(plotdata$Climate_PC1), y1 = max(plotdata$Climate_PC1)*(0.012900) + 0.217415, lwd = 2, lty = 2)
text(x = -3.5, y = 0.90, labels = "i")
# text(x = -1.9, y = 0.60, labels = "R2 = 0.28")
# text(x = -1.9, y = 0.55, labels = "p = 0.019")

## end Fig. 4






##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
#### For Figs. S1-S4: Interactive 3D plot of nAtomP_75 (aromaticity), ALogP_25 (hydrophilicity), and Elevation
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

library(plotly)
##### labels include plot, family, species, abundance, elevation and abundance-weighted CSCS for nAtomP_75 and ALogP_25 when hovered over


colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')
fig <- plot_ly(metaprops, x = ~nAtomP_75, y = ~ALogP_25, z = ~Elevation, color = ~Plot, size = ~Abundance, colors = colors,
             marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 50),
             text = ~paste('Plot:', Plot, '<br>Species:', Species_binomial, '<br>Family:', Family, '<br>Similarity Lower Quartile AlogP:', ALogP_25, '<br>Similarity Upper Quartile nAtomP:', nAtomP_75,
                           '<br>Abundance:', Abundance, '<br>Elevation', Elevation, 'm'))
fig <- fig %>% layout(title = 'Similarity for most light-absorbing (upper-quartile aromatic nAtomP) and hydrophilic (lower-quartile AlogP) metabolites increases with elevation',
         scene = list(xaxis = list(title = 'Similarity: most light-aborbing metabolites (upper-quartile nAtomP)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(0, 0.9),
                      # type = 'log',
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwidth = 2),
               yaxis = list(title = 'Similarity: most hydrophilic metabolites (lower-quartile AlogP)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(0, 0.9),
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwith = 2),
               zaxis = list(title = 'Elevation (m)',
                            gridcolor = 'rgb(255, 255, 255)',
                            # type = 'log',
                            zerolinewidth = 1,
                            ticklen = 5,
                            gridwith = 2)),
         paper_bgcolor = 'rgb(243, 243, 243)',
         plot_bgcolor = 'rgb(243, 243, 243)')

fig

htmlwidgets::saveWidget(fig, "FigS1_3Dplot_Elev_nAtomP_AlogP.html")




############ 
############ Figs. S2-S4: Prepare metapropstraits to facilitate plots that highlight certain families or generas
############ 


metapropstraits = metaprops

metapropstraits$HighGen5 = NA ## highlight 5 distantly related and abundance high-elevation genera
metapropstraits$LowGen5 = NA ##  highlight 5 species-rich lowland genera
metapropstraits$ExampleFam = NA ## highlight 5 species-rich families

for(i in 1:nrow(metapropstraits)){
	if(i != 669){
		if(metapropstraits$Family[i] == "Lauraceae"){
		metapropstraits$ExampleFam[i] = "Lauraceae"
		}
		if(metapropstraits$Family[i] == "Fabaceae"){
			metapropstraits$ExampleFam[i] = "Fabaceae"
		}
		if(metapropstraits$Family[i] == "Melastomataceae"){
			metapropstraits$ExampleFam[i] = "Melastomataceae"
		}
		if(metapropstraits$Family[i] == "Rubiaceae"){
			metapropstraits$ExampleFam[i] = "Rubiaceae"
		}
		if(metapropstraits$Family[i] == "Myrtaceae"){
			metapropstraits$ExampleFam[i] = "Myrtaceae"
		}
		# if(metapropstraits$Family[i] == "Cunoniaceae"){
			# metapropstraits$ExampleFam[i] = "Cunoniaceae"
		# }
		# if(metapropstraits$Family[i] == "Clethraceae"){
			# metapropstraits$ExampleFam[i] = "Clethraceae"
		# }
	
	
		if(metapropstraits$Genus[i] == "Weinmannia"){
			metapropstraits$HighGen5[i] = "Weinmannia"
		}
		if(metapropstraits$Genus[i] == "Prunus"){
			metapropstraits$HighGen5[i] = "Prunus"
		}
		if(metapropstraits$Genus[i] == "Persea"){
			metapropstraits$HighGen5[i] = "Persea"
		}
		if(metapropstraits$Genus[i] == "Clethra"){
			metapropstraits$HighGen5[i] = "Clethra"
		}
		if(metapropstraits$Genus[i] == "Hedyosmum"){
			metapropstraits$HighGen5[i] = "Hedyosmum"
		}
		# if(metapropstraits$Genus[i] == "Ilex"){
			# metapropstraits$HighGen5[i] = "Ilex"
		# }
	
	
		if(metapropstraits$Genus[i] == "Inga"){
			metapropstraits$LowGen5[i] = "Inga"
		}
		if(metapropstraits$Genus[i] %in% c("Ocotea", "Nectandra")){
			metapropstraits$LowGen5[i] = "Ocotea/Nectandra"
		}
		if(metapropstraits$Genus[i] == "Miconia"){
			metapropstraits$LowGen5[i] = "Miconia"
		}
		if(metapropstraits$Genus[i] == "Pouteria"){
			metapropstraits$LowGen5[i] = "Pouteria"
		}
		if(metapropstraits$Genus[i] == "Myrcia"){
			metapropstraits$LowGen5[i] = "Myrcia"
		}
	}
}


#####
#####
##### Fig. S2: Five Families: Fabaceae, Lauraceae, Melastomataceae, Myrtaceae, Rubiaceae
#####
#####


colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')
fig <- plot_ly(metapropstraits, x = ~nAtomP_75, y = ~ALogP_25, z = ~Elevation, color = ~ExampleFam, size = ~Abundance, colors = colors,
             marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 50),
             text = ~paste('Plot:', Plot, '<br>Species:', Species_binomial, '<br>Family:', Family, '<br>Similarity Lower Quartile AlogP:', ALogP_25, '<br>Similarity Upper Quartile nAtomP:', nAtomP_75,
                           '<br>Abundance:', Abundance, '<br>Elevation', Elevation, 'm'))
fig <- fig %>% layout(title = 'Similarity for most light-absorbing (upper-quartile aromatic nAtomP) and hydrophilic (lower-quartile AlogP) metabolites increases with elevation',
         scene = list(xaxis = list(title = 'Similarity: most light-aborbing metabolites (upper-quartile nAtomP)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(0, 0.9),
                      # type = 'log',
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwidth = 2),
               yaxis = list(title = 'Similarity: most hydrophilic metabolites (lower-quartile AlogP)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(0, 0.9),
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwith = 2),
               zaxis = list(title = 'Elevation (m)',
                            gridcolor = 'rgb(255, 255, 255)',
                            # type = 'log',
                            zerolinewidth = 1,
                            ticklen = 5,
                            gridwith = 2)),
         paper_bgcolor = 'rgb(243, 243, 243)',
         plot_bgcolor = 'rgb(243, 243, 243)')

fig

htmlwidgets::saveWidget(fig, "FigS2_3Dplot_Elev_nAtomP_AlogP_5families.html")


#####
#####
##### Fig. S3: Diverse Lowland Genera: Inga, Miconia, Myrcia, Ocotea/Nectandra, Pouteria
#####
#####


colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')
fig <- plot_ly(metapropstraits, x = ~nAtomP_75, y = ~ALogP_25, z = ~Elevation, color = ~LowGen5, size = ~Abundance, colors = colors,
             marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 50),
             text = ~paste('Plot:', Plot, '<br>Species:', Species_binomial, '<br>Family:', Family, '<br>Similarity Lower Quartile AlogP:', ALogP_25, '<br>Similarity Upper Quartile nAtomP:', nAtomP_75,
                           '<br>Abundance:', Abundance, '<br>Elevation', Elevation, 'm'))
fig <- fig %>% layout(title = 'Similarity for most light-absorbing (upper-quartile aromatic nAtomP) and hydrophilic (lower-quartile AlogP) metabolites increases with elevation',
         scene = list(xaxis = list(title = 'Similarity: most light-aborbing metabolites (upper-quartile nAtomP)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(0, 0.9),
                      # type = 'log',
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwidth = 2),
               yaxis = list(title = 'Similarity: most hydrophilic metabolites (lower-quartile AlogP)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(0, 0.9),
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwith = 2),
               zaxis = list(title = 'Elevation (m)',
                            gridcolor = 'rgb(255, 255, 255)',
                            # type = 'log',
                            zerolinewidth = 1,
                            ticklen = 5,
                            gridwith = 2)),
         paper_bgcolor = 'rgb(243, 243, 243)',
         plot_bgcolor = 'rgb(243, 243, 243)')

fig

htmlwidgets::saveWidget(fig, "FigS3_3Dplot_Elev_nAtomP_AlogP_5DivGen.html")


#####
#####
##### Fig. S4: Independent Highland Genera: Clethra, Hedyosmum, Persea, Prunus, Weinmannia 
#####
#####


colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')
fig <- plot_ly(metapropstraits, x = ~nAtomP_75, y = ~ALogP_25, z = ~Elevation, color = ~HighGen5, size = ~Abundance, colors = colors,
             marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 50),
             text = ~paste('Plot:', Plot, '<br>Species:', Species_binomial, '<br>Family:', Family, '<br>Similarity Lower Quartile AlogP:', ALogP_25, '<br>Similarity Upper Quartile nAtomP:', nAtomP_75,
                           '<br>Abundance:', Abundance, '<br>Elevation', Elevation, 'm'))
fig <- fig %>% layout(title = 'Similarity for most light-absorbing (upper-quartile aromatic nAtomP) and hydrophilic (lower-quartile AlogP) metabolites increases with elevation',
         scene = list(xaxis = list(title = 'Similarity: most light-aborbing metabolites (upper-quartile nAtomP)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(0, 0.9),
                      # type = 'log',
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwidth = 2),
               yaxis = list(title = 'Similarity: most hydrophilic metabolites (lower-quartile AlogP)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(0, 0.9),
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwith = 2),
               zaxis = list(title = 'Elevation (m)',
                            gridcolor = 'rgb(255, 255, 255)',
                            # type = 'log',
                            zerolinewidth = 1,
                            ticklen = 5,
                            gridwith = 2)),
         paper_bgcolor = 'rgb(243, 243, 243)',
         plot_bgcolor = 'rgb(243, 243, 243)')

fig

htmlwidgets::saveWidget(fig, "FigS4_3Dplot_Elev_nAtomP_AlogP_5HighGen.html")



