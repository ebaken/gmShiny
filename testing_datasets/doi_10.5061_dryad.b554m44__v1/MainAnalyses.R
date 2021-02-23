# Ancestral State Reconstruction

# Libraries #### 
library(ape)
library(corHMM) # ML ASR analyses
library(phytools) # Bayesian ASR analyses
library(geomorph) # Morphological analyses
library(geiger) # for pruning trees and datasets to each other (function datatree)

# Ancestral State Reconstruction, Transition Rates and Counts #### 

Tree <- read.tree("Bonett_and_Blair_Phylo.tre") # Full Tree
Microhabitats <- read.csv("Microhabitats.csv")  # Full Tree Microhabitat Data

# Identifying the root of family Plethodontidae and global root
Microhabitats.Pleth <- subset(Microhabitats$Species, Microhabitats$Pleth == "y") 
PlethNode <- getMRCA(Tree, tip = as.character(Microhabitats.Pleth))

# ML: model comparison and node estimates
Models <- list("ER", "SYM", "ARD")
Micros <- cbind(as.character(Microhabitats$Species), as.character(Microhabitats$M6))

AllModelOutput <- lapply(Models, function(x) {
  rayDISC(Tree, Micros, ntraits = 1, model = x, node.states = "marginal", verbose = F,  diagn = F)
  }) 
AICs <- lapply(1:length(AllModelOutput), function(x) {AllModelOutput[[x]]$AIC}) # AIC values
Nodes.ARD <- AllModelOutput[[3]]$states # all node estimates from best fit model (ARD)
write.csv(AllModelOutput[[3]]$solution, "Q.ML.ARD.csv", row.names = T)
      
# Bayesian stochastic mapping: node and transition estimates
Micros <- as.factor(Microhabitats$M6)
names(Micros) <- Microhabitats$Species
q <- read.csv("Q.ML.ARD.csv")
      row.names(q) <- q[,1]
      q <- q[,-1]
      q <- as.matrix(q)
      for (i in 1:6) {q[i,i] <- sum(q[i,-i]) }
      diag(q) <- diag(q)*-1
      
StochMap <- make.simmap(Tree, Micros, model = "ARD", nsim = 1000, Q = q)
StochMap.Summary <- summary(StochMap) 
StochMap.Summary$ace[1,] # ASR estimate at root of all Caudata
StochMap.Summary$ace[(PlethNode-length(Tree$tip.label)),] # ASR estimate at root of Plethodontidae 
colSums(StochMap.Summary$count)/1000 # estimated number of transitions for all Caudata

# Linear Measurement Analyses ####

Tree <- read.tree("Bonett_and_Blair_Phylo.tre") # Full Tree
Microhabitats <- read.csv("Microhabitats.csv")  # Full Tree Microhabitat Data
LinMeas <- read.csv("LinMeas.SpeciesSummaries.csv")

row.names(LinMeas) <- LinMeas$Species
TreeDataPruned <- treedata(Tree, LinMeas) 
TreePruned <- TreeDataPruned$phy
LinMeasOrdered <- LinMeas[match(TreePruned$tip.label, LinMeas$Species),]
LinMeasOrdered$Species <- as.character(LinMeasOrdered$Species)
MicrohabitatsPruned <- Microhabitats[match(TreePruned$tip.label, Microhabitats$Species),] # Getting linear measurements in the right order
MicrohabitatsPruned$Species <- as.character(MicrohabitatsPruned$Species)

Micros <- MicrohabitatsPruned$M6
names(Micros) <- MicrohabitatsPruned$Species

ShapeLogAdjusted <- log(LinMeasOrdered[,c("SE", "HL", "BWL", "FLL", "HLL", "TL")]/LinMeasOrdered$SVL)
row.names(ShapeLogAdjusted)<-LinMeasOrdered$Species
ShapeLogAdjusted <- as.matrix(ShapeLogAdjusted) # necessary for geomorph.data.frame

GDF <- geomorph.data.frame(shape = ShapeLogAdjusted, habitat = Micros, phy = TreePruned) 
Results <- procD.pgls(shape ~ habitat, iter = 999, data = GDF, phy = phy)
summary(Results) 
Pairwise <- pairwise(Results, groups = Micros)
summary(Pairwise) 
   
ShapeArray <- arrayspecs(ShapeLogAdjusted, 3, 2) # necessary for evol.rates function
Rates <- compare.evol.rates(ShapeArray, TreePruned, gp = Micros) 
Rates
Rates$pairwise.pvalue  
            
# Foot Shape Analyses ####

Tree <- read.tree("Bonett_and_Blair_Phylo.tre") # Full Tree
Microhabitats <- read.csv("Microhabitats.csv")  # Full Tree Microhabitat Data
FootLMs <- readland.tps("FootLandmarks.BB.SpeciesSummary.tps", specID = "ID")
n
dimnames(FootLMs)[[3]] <- stringr::str_replace_all(dimnames(FootLMs)[[3]], " ", "_")
FootCsize <- read.csv("FootLandmarks.SpeciesSummary.Csize.csv")
Csize <- FootCsize[,2]
names(Csize) <- stringr::str_replace_all(FootCsize[,1], " ", "_")
MicrohabitatsPruned <- Microhabitats[match(dimnames(FootLMs)[[3]],Microhabitats$Species),] 
row.names(MicrohabitatsPruned) <- as.character(MicrohabitatsPruned$Species)

PrunedPhyData <- treedata(Tree, MicrohabitatsPruned)
TreePruned <- PrunedPhyData$phy
MicrohabitatsPruned <- PrunedPhyData$data
MicrohabitatsPruned <- MicrohabitatsPruned[match(TreePruned$tip.label, MicrohabitatsPruned[,1]),]
FootLMs <- FootLMs[,,match(TreePruned$tip.label, dimnames(FootLMs)[[3]])]
Csize <- Csize[match(TreePruned$tip.label, names(Csize))]

Micros <- MicrohabitatsPruned[,2]
names(Micros) <- MicrohabitatsPruned[,1]

gdf <- geomorph.data.frame(shape = FootLMs, phy = TreePruned, grp = Micros, size=Csize)

Results <- procD.pgls(shape ~ grp, phy = phy, data = gdf)
summary(Results)
    
Pairwise <- pairwise(Results, groups = Micros) 
summary(Pairwise)

EvolRates <- compare.evol.rates(FootLMs, phy = TreePruned, gp = Micros, iter = 999) 
EvolRates
EvolRates$pairwise

# Test of Allometric convergence

plethAllometry <- procD.lm(shape ~ size*grp, data = gdf, print.progress = F)
NewAllometry <- plotAllometry(plethAllometry, size = gdf$size, logsz = FALSE, method = "PredLine")

minsize<-match(tapply(Csize, Micros, min),Csize)[1:2] # only selecting 1 and 2 to just be A and C
maxsize<-match(tapply(Csize, Micros, max),Csize)[1:2] 

yhat.full<-NewAllometry$PredLine #predict shapes
shape.min<-yhat.full[minsize]; shape.max<-yhat.full[maxsize]
Dist.diff.obs<-sum(dist(shape.min))-sum(dist(shape.max))

PDiff<-1
permute<-99 
Diff<-rep(NA, permute)

for(k in 1:99){
  Micros.rand <- sample(Micros)       # changing the microhabitat classifications
  names(Micros.rand) <- names(Micros) # keeping the names in the right order with random microhabitat classification
  GDF <- geomorph.data.frame(shape = FootLMs, grp = Micros.rand, size = Csize)
  X <- procD.lm(shape ~ size*grp, data = GDF, print.progress = F) 
  NewAllometry <- plotAllometry(X, size = GDF$size, logsz = FALSE, method = "PredLine")
  yhat.rand <- NewAllometry$PredLine
    
  minsize<-match(tapply(Csize, Micros.rand, min),Csize)[1:2]
  maxsize<-match(tapply(Csize, Micros.rand, max),Csize)[1:2] 
  shape.min.r<-yhat.rand[minsize]
  shape.max.r<-yhat.rand[maxsize]
  Dist.diff.rand<-sum(dist(shape.min.r))-sum(dist(shape.max.r)) 
  Diff<-rbind(Diff,Dist.diff.rand)
  PDiff<-ifelse(Dist.diff.rand>=Dist.diff.obs,PDiff+1,PDiff)

}  

PDiff<-PDiff/(100)
PDiff # p value 


# Convergence Analyses  
library(convevol)

ShapeLogAdjustedM <- matrix(FootLMs, nrow = 288, ncol = 42)
rownames(ShapeLogAdjustedM) <- dimnames(FootLMs)[[3]]
ASpecies <- subset(names(Micros), Micros == "A")

ObservedC5 <- convnum(TreePruned, ShapeLogAdjustedM, ASpecies, plot = F)
OutputC5 <- convnumsig(TreePruned, ShapeLogAdjustedM, ASpecies, nsim = 100,  
            ellipse = NULL, plot = F, plotellipse = NULL) 

ObservedC1 <- convrat(TreePruned, ShapeLogAdjustedM, ASpecies, plot = F)
OutputC1 <- convratsig(TreePruned, ShapeLogAdjustedM, ASpecies, nsim = 100, plot = F)
    