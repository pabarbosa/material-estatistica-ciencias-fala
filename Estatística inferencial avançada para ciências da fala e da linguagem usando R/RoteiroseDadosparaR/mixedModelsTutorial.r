# ############################################################################################################# #
# THIS CODE BRIEFLY DEMONSTRATES HOW TO RUN A MIXED-EFFECTS MODEL IN R.  THE SAMPLE DATA ARE FROM THE PHOIBLE   #
# KNOWLEDGE BASE, AND COMPRISE 969 LANGUAGES WITH POPULATION, GENEALOGICAL, AND PHONOLOGICAL INFORMATION ABOUT  #
# EACH.  THE MODEL SHOWN HAS TWO (NESTED) RANDOM EFFECTS AND ONE FIXED EFFECT, SO IT IS RATHER DIFFERENT THAN   #
# MOST PSYCHOLINGUISTIC MODELS (WHICH TYPICALLY HAVE SEVERAL FIXED EFFECT PREDICTORS AND INTERACTIONS AMONG     #
# THEM).  THE MAIN FOCUS IS CHECKING THAT MODEL ASSUMPTIONS ARE NOT VIOLATED, EXPLAINING HOW TO UNDERSTAND      #
# MODEL SUMMARIES, AND COMPARISON OF COMPETING MODELS FOR THE SAME DATA.                                        #
#                                                                                                               #
# NOTE: THE REASONING IN THIS TUTORIAL ABOUT HOW TO INTERPRET THE MODELS ROUGHLY FOLLOWS THAT LAID OUT IN THE   #
# PRESENTATION Moran et al 2012 (see references below), and the related draft manuscript. Both the presentation #
# and the draft manuscript are available from: http://depts.washington.edu/phonlab/projects                     #
#                                                                                                               #
# TUTORIAL VERSION 0.2 (2012 01 10)                                                                             #
#                                                                                                               #
# AUTHOR: DANIEL MCCLOY (drmccloy@uw.edu)                                                                       #
# LICENSED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE: http://creativecommons.org/licenses/by/3.0/        #
# ############################################################################################################# #

# GETTING STARTED
# Packages required for this tutorial: "lme4" (for creating the mixed models)
# Optional: "languageR" (not used in this tutorial but often useful for mixed models when calculating p-values
# using MCMC simulation; see discussion below.)


# ############################ #
# LIBRARIES, SOURCE DATA, ETC. #
# ############################ #

# LOAD THE lme4 LIBRARY
library(lme4)

# TELL R WHERE TO FIND YOUR DATA
setwd("/home/dan/Desktop/demo/")

# LOAD IN THE DATA
phoibleData <- read.delim("phoibleData.tsv", header=T, sep="\t", quote="\"", dec=".", stringsAsFactors=F)

# TAKE A LOOK AT HOW THE DATA IS STRUCTURED
head(phoibleData)


# ######## #
# ANALYSIS #
# ######## #

# HERE WE DO A BASIC SCATTERPLOT OF THE DATA WITH ORDINARY LEAST SQUARES REGRESSION LINE, JUST TO SEE WHAT IT LOOKS LIKE.
plot(phoibleData$logPop, phoibleData$pho, pch=16, col=rgb(0,0,204,102,maxColorValue=255))
olsLine <- lm(phoibleData$pho ~ phoibleData$logPop)
abline(olsLine, col="red")

# THE SUMMARY OF THE BASIC O.L.S. REGRESSION SUGGESTS THAT THERE IS A STATISICALLY SIGNIFICANT CORRELATION BETWEEN LOG(POPULATION) AND NUMBER OF PHONEMES.
summary(olsLine)

# HOWEVER, WE HAVEN'T YET CHECKED THE DATA FOR ADHERANCE TO THE MODEL ASSUMPTIONS. FOR EXAMPLE, THE Q-Q PLOT OF RESIDUALS SHOWS THAT THE RESIDUALS ARE NOT
# NORMALLY DISTRIBUTED (INDICATED BY THEIR DEVIATION FROM THE Q-Q LINE), WHICH IS AN IMPORTANT ASSUMPTION OF LINEAR MODELLING.
qqnorm(residuals(olsLine))
qqline(residuals(olsLine))

# TO CORRECT THIS, WE CAN LOG-TRANSFORM THE OUTCOME VARIABLE TOO:
phoibleData$logPho <- log10(phoibleData$pho)
ols2 <- lm(phoibleData$logPho ~ phoibleData$logPop)
qqnorm(residuals(ols2))
qqline(residuals(ols2))

# AT THIS POINT THE DATA STILL VIOLATES THE MOST IMPORTANT ASSUMPTION OF LINEAR MODELLING: NAMELY, INDEPENDENCE OF THE DATA POINTS. VIOLATING THE ASSUMPTION OF
# INDEPENDENCE CAN INFLATE TYPE-I ERROR RATE DRAMATICALLY. LANGUAGES (OR THEIR PHONEME INVENTORY COUNTS) ARE "NESTED" WITHIN GENEALOGICAL GROUPS (WE HAVE GROUPINGS
# BY GENUS AND BY TOP-LEVEL FAMILY) AND AS SUCH ARE NON-INDEPENDENT. MOREOVER, OUR GROUPS ARE NOT ALL THE SAME SIZE, AS SEEN IN THIS TABLE. 
familyTable <- data.frame(table(phoibleData$fam))
colnames(familyTable) <- c("fam","numLangs")
familyTable

# FORTUNATELY, MIXED EFFECTS MODELS ARE DESIGNED TO HANDLE NESTED DATA AND UNEQUAL GROUP SIZES. HERE IS THE BASIC MIXED MODEL CALL.
mixMod.rirs <- lmer(pho~logPop+(1+logPop|genus)+(1+logPop|fam), data=phoibleData)

# EXPLANATION OF ABOVE: The name "mixMod.rirs" is a reminder that this is a mixed model with "random intercepts random slopes", and is merely a handy mnemonic.
# The syntax of the function call above goes like this:  lmer(outcomeVariable ~ fixedEffect + (1 + fixedEffect | groupingFactorA) + (1 + fixedEffect | groupingFactorB))
# The parts like (1 + fixedEffect | groupingFactorA) are the random effects terms: ...
# ...the "1" gets you a (potentially) different intercept for each group in A, while the "fixedEffect" part within the random effect gets you a (potentially) different
# slope for each group in A. Here we have two random effects: family and genus.  The lmer() function is smart enough to figure out that genus is nested within family,
# and treats it accordingly. If the grouping factors had been crossed, lmer() would also have known what to do.

# ALTHOUGH IT HANDLES NESTED DATA AND UNEQUAL GROUP SIZES WELL, MIXED EFFECTS MODELS STILL REQUIRE THAT THE DATA MEET THE "NORMALITY OF RESIDUALS" CRITERION.
qqnorm(residuals(mixMod.rirs))
qqline(residuals(mixMod.rirs))

# As you can see, the residuals are not normally distributed: languages with higher logPop values tend to have larger residuals. 
# Just like we did before with the O.L.S. model, we rerun the mixed model with a log-transformed outcome:
mixMod2.rirs <- lmer(logPho~logPop+(1+logPop|genus)+(1+logPop|fam), data=phoibleData)
qqnorm(residuals(mixMod2.rirs))
qqline(residuals(mixMod2.rirs))

# The residuals are still a little skewed, but much less than before, and close enough to normal that we can safely interpret the model results.
# So let's look at the model summary:
summary(mixMod2.rirs)

# Here is what we get:
# Linear mixed model fit by REML 
# Formula: logPho ~ logPop + (1 + logPop | genus) + (1 + logPop | fam) 
#    Data: phoibleData 
#    AIC   BIC logLik deviance REMLdev
#  -1377 -1333  697.6    -1411   -1395
# Random effects:
#  Groups   Name        Variance   Std.Dev.   Corr   
#  genus    (Intercept) 5.3211e-13 7.2946e-07        
#           logPop      7.7167e-05 8.7845e-03 0.000  
#  fam      (Intercept) 1.6178e-02 1.2719e-01        
#           logPop      1.2281e-04 1.1082e-02 -0.654 
#  Residual             1.0729e-02 1.0358e-01        
# Number of obs: 969, groups: genus, 321; fam, 100
# 
# Fixed effects:
#             Estimate Std. Error t value
# (Intercept) 1.442268   0.020359   70.84
# logPop      0.009275   0.004099    2.26
# 
# Correlation of Fixed Effects:
#        (Intr)
# logPop -0.774

# First look at the Fixed effects. We're mostly interested in the *relationship* between population and phonemes, not in any notion of *average* inventory
# size. As such, for the most part we ignore the intercepts and focus on the slope (which is on the "logPop" line).  The estimate of 0.009275 means that
# for every 1% increase in logPop, we get a predicted increase of 0.009275% in the phoneme inventory size (This "percent change" interpretation works because
# both the outcome and predictor are log-transformed).  We can make a small table just to get a sense of the range of predicted values:
pop <- c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
pho <- 10^1.442268 * pop^0.009275
predictions <- cbind(pop,pho)
predictions

# That shows us that the model predicts an increase of about 6 phonemes if we go from a language with 1 speaker to a language with 1 billion speakers.
# Not a very big predicted change!

# You may have noticed that lmer() summary output does not provide p-values for predictors. There are good reasons why (see, e.g., Douglas Bates's note here:
# https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html, or the article by Baayen et al 2008).  To assess the significance of population as a predictor
# in the absence of provided p-values, we can look at the t-value of the fixed effects.  Intercept is expected to have a high t-value since phoneme counts
# are always a positive number (so we expect the intercept to be different from zero).  What interests us is the t-value of the slope estimate, which is 2.26.
# Since our dataset is fairly large (969 observations), we expect the t-distribution to closely approximate the normal distribution, and as such we can estimate
# that a predictor is significant if its t-value is greater than 2 (or less than -2) (on this point see Baayen et al 2008).

# It is also possible to get p-values by running Markov-Chain Monte Carlo simulations (MCMC) and using the pvals.fnc() feature of Harald Baayen's languageR package.
# However, at present the MCMC function has not been fully implemented (i.e., it fails when the random effects have correlated coefficients; see the note from Douglas
# Bates at https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q3/001409.html).  In our model, the intercepts & slopes of the random effects (genera and families)
# are correlated (i.e., families with negative slopes have higher intercepts, families with positive slopes have lower intercepts). Such correlation is not
# necessarily a problem if we don't much care about interpreting the intercepts (Gelman & Hill 2007:288), except it keeps us from running the MCMC simulation.
# One way around this problem is to transform the predictor so that it is mean-centered at zero, so that the intercepts & slopes are uncorrelated.  The choice
# to mean-center a predictor is actually not trivial, as it changes the model specification slightly (see Snijders & Bosker 1999:80-81 for discussion).

# So, here, logPop is a marginally significant predictor. But recall that the predicted effect is quite small, and remember that larger and larger datasets
# make it easier and easier to get "significant" results. Thus "significant" does not necessarily equal "interesting" or "meaningful"!  Sometimes those
# significant results are just artefacts (especially when the predicted effect is so small)... that's likely the case here.

# Looking now at the random effects, we see that the slopes of lines for the grouping factor "family" vary with a standard deviation of 0.011082, 
# and that the genus lines have additional variability with a standard deviation of 0.0087845.  Compare those to the predicted "overall" slope of 0.009275 
# and you see that within-family slopes vary quite a bit compared to the overall slope, such that some families are modeled with positive slopes and others
# with negative slopes.  This is another clue that the "significant" fixed effect might be artefactual (if it were real, we would expect the pattern to hold
# within language families as well as across families).

# Since we suspect a spurious correlation, we can do a bit of model comparison.  Here we create a "null model" where there are no fixed effects and no random
# effects other than a unique intercept for each group (in other words, each genus / family is modeled as a horizontal line; this is a lot like a factorial ANOVA).
nullMod <- lmer(logPho ~ (1|genus) + (1|fam), data=phoibleData)

summary(nullMod)
# the summary of the null model:
# Linear mixed model fit by REML 
# Formula: logPho ~ (1 | genus) + (1 | fam) 
#    Data: phoibleData 
#    AIC   BIC logLik deviance REMLdev
#  -1385 -1365  696.3    -1399   -1393
# Random effects:
#  Groups   Name        Variance  Std.Dev.
#  genus    (Intercept) 0.0021626 0.046504
#  fam      (Intercept) 0.0112925 0.106266
#  Residual             0.0107254 0.103563
# Number of obs: 969, groups: genus, 321; fam, 100
# 
# Fixed effects:
#             Estimate Std. Error t value
# (Intercept)  1.47500    0.01318   111.9

# We can compare the null model to the earlier model to see which one better fits the data.  There are at least three ways to do this.  
# The first two are to compare AIC (Aikake Information Criterion) and BIC (Bayesian Information Criterion).
# The third way is to use the anova() function in lme4 to compare the log-likelihood of the two models.
anova(mixMod2.rirs, nullMod)

# Here is the output:
# Data: phoibleData
# Models:
# nullMod: logPho ~ (1 | genus) + (1 | fam)
# mixMod2.rirs: logPho ~ logPop + (1 + logPop | genus) + (1 + logPop | fam)
#              Df     AIC     BIC logLik  Chisq Chi Df Pr(>Chisq)  
# nullMod       4 -1391.5 -1372.0 699.73                           
# mixMod2.rirs  9 -1393.2 -1349.3 705.59 11.726      5    0.03874 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

# The first thing to note is that the values for AIC, BIC, and logLik have changed from the previous model summaries we looked at.  This is because the
# models typically use "Restricted Maximum Likelihood" to estimate the model parameters, but when comparing models it is necessary to use (Full) Maximum
# Likelihood estimation. The anova() function knows this and correctly uses ML instead of REML.  So, what we see is that the AIC is about the same in the
# two models, and the BIC is slightly lower in the full model (for both of these, a lower number indicates a better model).  For the log-likelihood, the
# full model will always (by definition) be higher than the null model, and the question is "is it significantly better?"  For this, we look at the last
# number, labeled "Pr(>Chisq)".  This is a p-value, so here we see that the full model is marginally significantly better than the null model.  Overall,
# the model comparison says that including population as a predictor is, if anything, only *slightly* better than omitting it entirely. Nonetheless, there
# does seem to be some variability in the data that is slightly correlated with population size. When interpreting such subtle effects it is important to 
# remember that it is quite possible that the correlation is artefactual (see, e.g., van der Laan & Rose 2010), or that population and phoneme inventory
# size covary due to a hidden factor affecting both of them (e.g., social network structure or something).


# POSTSCRIPT
# If you want to plot the prediction lines from the mixed model, you can do it like this:

plot(phoibleData$logPop, phoibleData$logPho, pch=16, col=rgb(0,0,204,102,maxColorValue=255)) #scatterplot
groupCoefficients <- data.frame(coef(mixMod2.rirs)$fam)
grandIntercept <- unlist(fixef(mixMod2.rirs)[1])
grandSlope <- unlist(fixef(mixMod2.rirs)[2])
for (i in 1:nrow(groupCoefficients)) {
  abline(groupCoefficients[i,1], groupCoefficients[i,2], col=rgb(0,0,0,102,maxColorValue=255))
}
abline(grandIntercept, grandSlope, col="red", lwd=1.5)


# REFERENCES
# Baayen, R. H., Davidson, D. J., & Bates, D. M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. Journal of Memory and Language, 59(4), 390-412. doi:10.1016/j.jml.2007.12.005
# Bates, D. M. (2006). [R] lmer, p-values and all that. https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html
# Bates, D. M. (2008). [R-sig-ME] mcmcsamp can't handle random effects syntax. https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q3/001409.html
# Gelman, A., & Hill, J. (2007). Data analysis using regression and multilevel/hierarchical models. New York: Cambridge University Press.
# van der Laan, M., & Rose, S. (2010). Statistics ready for a revolution. Amstat News. Retrieved October 7, 2011, from http://magazine.amstat.org/blog/2010/09/01/statrevolution/
# Moran, S., McCloy, D. R., & Wright, R. A. (2012). Revisiting the population vs phoneme-inventory correlation. Paper presented at the The 86th Meeting of the Linguistic Society of America, Portland, OR.
# Snijders, T. A. B., & Bosker, R. J. (1999). Multilevel analysis: An introduction to basic and advanced multilevel modeling. London: Sage.
 
