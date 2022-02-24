########     DIVERSIFICATION OF SOUTH AMERICAN VESPERTILIONIDAE (CHIROPTERA)
########          IS NOT CONSTRAINED BY EVOLUTIONARY PRIORITY EFFECTS

#####  Phylogenetic Handling

# Loading packages
require(ape)
require(phytools)
require(geiger)

# Loading list of species separated by group occurrence (according to IUCN 2020 range maps)
# species selected are the ones present in the phylogenetic trees described in the text
species.list<-read.table('species.list.txt', header = TRUE, stringsAsFactors = FALSE)

###  Loading the phylogeny of Shi & Rabosky (2015) - downloaded from https://doi.org/10.1111/evo.12681)
tree.full <- read.nexus("tree_Shi&Rabosky2015.nex")

# rename species according to IUCN 2020 maps
tree.full$tip.label[match('Pipistrellus_subflavus', tree.full$tip.label)]<-'Perimyotis_subflavus'
tree.full$tip.label[match('Rhogeessa_gracilis', tree.full$tip.label)]<-'Baeodon_gracilis'
tree.full$tip.label[match('Rhogeessa_alleni', tree.full$tip.label)]<-'Baeodon_alleni'

# remove absent species
comparison<-name.check(tree.full, species.list[,1], data.names = species.list[,1])
tree.untacted<-drop.tip(tree.full, comparison$tree_not_data)  # 110 species included in the phylogeny

# writing the backbone phylogeny used for TACT
write.tree(tree.untacted, "backbone.newick")

### Now run TACT with the 'backbone.newick' as the phylogeny and the 'poly.tree.csv'
# which contains missing species based on the topology of the trees mentioned in the text.
# A tutorial can be seen here --> https://tact.jonathanchang.org/tutorial/ 

### Loading the phylogeny stochastically resolved by TACT
tree<-read.tree('tree.tacted.newick.tre')   # contains 127 species tips

### Extracting vespertilionids from the phylogeny to estimate ancestral area
require(dplyr)
vesper.list<-filter(species.list, group.occurrence == 'vesper.na' | group.occurrence == 'vesper.sa' | group.occurrence == 'both')
comparison.vesper<-name.check(tree, vesper.list[,1], data.names = vesper.list[,1])
vesper.tree<-drop.tip(tree, comparison.vesper$tree_not_data)  # 71 species
write.tree(vesper.tree, "vesper.tree.newick")  # writing the phylogeny


#####  ANCESTRAL AREA ESTIMATION - BioGeoBEARS
library(BioGeoBEARS)

# Setting up the directory directory:
dir_data = getwd()

# Inserting the phylogeny
trfn = np(paste(addslash(dir_data), "vesper.tree.newick", sep=""))
moref(trfn)
tr<-read.tree(trfn)

# Loading geographic data
geogfn = np(paste(addslash(dir_data), "geodata.txt", sep=""))
moref(geogfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
max_range_size = 2

### Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Number of computer cores to use (WTF IS THIS, its not in the script, I saw only after trying to run the analyses and getting the error)
BioGeoBEARS_run_object$num_cores_to_use = TRUE

### RUN the DEC model (finally)
# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "Vesper_DEC_M0_unconstrained_v1.Rdata"
res.dec<-if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}
res.dec

### Plot ancestral states - DEC
analysis_titletxt ="BioGeoBEARS DEC on Vespertilionidae M0_unconstrained"
results_object = resDEC
scriptdir = getwd()   # make sure the function "plot_phylo3_nodecoords" is in the directory
# that function is found at '[...]R\win-library\3.6\BioGeoBEARS\extdata\a_scripts'
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.6, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


#####  DIVERSIFICATION ANALYSES - BAMM
require(BAMMtools)

### Setting up BAMM priors
setBAMMpriors(tree)   # automatically generates a file named 'myPriors' in the working directory

# Now run BAMM (with the control file 'bat.control.file.txt', which had the priors changed according
# to the 'setBAMMpriors' function) and return here afterwards to evaluate the output data (event.data.txt)
# A tutorial to run BAMM can be found here --> http://bamm-project.org/settingup.html#running

# Copy the 'event.data.txt' and the 'mcmc_out.txt' to the working directory
# Loading the event data
edata <- getEventData(tree, eventdata = 'event_data.txt', burnin=0.1)
edata

### Assessing MCMC convergence (ran with five million generations)
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)  # check graphical convergence

# Discard some as burnin
burnstart <- floor(0.1 * nrow(mcmcout))  # 0.1
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# Check the effective sample sizes of the log-likelihood and the number of shift events present in each sample
library(coda)

# Values should be above 200 to achieve convergence, according to http://bamm-project.org/postprocess.html#assessing-mcmc-convergence
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

# Posterior probability of the number of diversification rate shifts
shift_probs <- summary(edata)   # '0' (an homogeneous rate) has 86% of probability

### Prior distribution in BAMM
bayes.f <- computeBayesFactors(mcmcout, expectedNumberOfShifts = 1, burnin = 0.1)
bayes.f   # 0 scores higher than all others
plotPrior(mcmcout, expectedNumberOfShifts=1)

### Bayesian credible sets of shift configurations
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)   # just one credible set, with 0 shifts
plot.credibleshiftset(css)

### Plot diversification rates
plot.bamm <- plot.bammdata(edata, labels = T, cex = 0.25, lwd = 2)
addBAMMlegend(plot.bamm, location = c(5, 20, 46, 47))
addBAMMshifts(edata, cex=2)   # none

# Rate-through-time plot
st <- max(branching.times(tree))   # get start time
plotRateThroughTime(edata, avgCol="black", start.time=st, ylim=c(0,1), cex.axis=1, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)   # confidence interval denotes 5% through 95% Bayesian credible regions

### Macroevolutionary cohort analysis
cmat <- getCohortMatrix(edata)
cohorts(cmat, edata, lwd = 0.5)


#####  GEOGRAPHIC-DEPENDENT ANALYSES - GeoHiSSE
library(hisse)
library(diversitree)

# Phylogeny (only vesper bats were employed in GeoHiSSE)
vesper.tree

# Load Geographic data
geo.data<-read.table('geohisse.txt', header = TRUE)

###  Running the models

# Model 1. Geographic-independent diversification without hidden states
# 4 parameters: 1 speciation, 1 extinction, 2 dispersal
turnover <- c(1,1,1)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
mod1 <- GeoHiSSE(phy = vesper.tree, data = geo.data, f=c(0.66,0.66,0.66),
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate)

# Model 2. Geographic-independent diversification with hidden states
# 9 parameters: 1 speciation(+1 hidden), 1 extinction(+1 hidden), 2 dispersal(+2 hidden), 1 dispersal between hidden states
turnover <- c(1,1,1,2,2,2)
eps <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
mod2 <- GeoHiSSE(phy = vesper.tree, data = geo.data, f=c(0.66,0.66,0.66),
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate)

# Model 3. Geographic-independent speciation (unequal extinction) without hidden states
# 5 parameters: 1 speciation, 2 extinction, 2 dispersal
turnover <- c(1,1,1)
eps <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
mod3 <- GeoHiSSE(phy = vesper.tree, data = geo.data, f=c(0.66,0.66,0.66),
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate)

# Model 4. Geographic-independent speciation (unequal extinction) with hidden states
# 11 parameters: 1 speciation(+1 hidden), 2 extinction(+2 hidden), 2 dispersal(+2 hidden), 1 dispersal between hidden states
turnover <- c(1,1,1,2,2,2)
eps <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
mod4 <- GeoHiSSE(phy = vesper.tree, data = geo.data, f=c(0.66,0.66,0.66),
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate)

# Model 5. Geographic-dependent speciation (equal extinction) without hidden states
# 6 parameters: 3 speciation, 1 extinction, 2 dispersal
turnover <- c(1,2,3)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
mod5 <- GeoHiSSE(phy = vesper.tree, data = geo.data, f=c(0.66,0.66,0.66),
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate)

# Model 6. Geographic-dependent speciation (equal extinction) with hidden states
# 13 parameters: 3 speciation(+3 hidden), 1 extinction(+1 hidden), 2 dispersal(+2 hidden), 1 dispersal between hidden states
turnover <- c(1,2,3,4,5,6)
eps <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
mod6 <- GeoHiSSE(phy = vesper.tree, data = geo.data, f=c(0.66,0.66,0.66),
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate)

# Model 7. Geographic-dependent diversification without hidden states
# 7 parameters: 3 speciation, 2 extinction, 2 dispersal
turnover <- c(1,2,3)
eps <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
mod7 <- GeoHiSSE(phy = vesper.tree, data = geo.data, f=c(0.66,0.66,0.66),
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate)

# Model 8. Geographic-dependent diversification with hidden states
# 15 parameters: 3 speciation(+3 hidden), 2 extinction(+2 hidden), 2 dispersal(+2 hidden), 1 dispersal between hidden states
turnover <- c(1,2,3,4,5,6)
eps <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
mod8 <- GeoHiSSE(phy = vesper.tree, data = geo.data, f=c(0.66,0.66,0.66),
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate)

# Akaike weights
GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4,
                   model5 = mod5, model6 = mod6, model7 = mod7, model8 = mod8), criterion="AIC")
