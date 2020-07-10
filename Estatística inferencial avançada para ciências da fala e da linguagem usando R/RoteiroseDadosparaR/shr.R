SHR <- function(r,firstfactor,secondfactor,withinteraction){
if (withinteraction == 1) lm1 <- lm(rank(r)~firstfactor*secondfactor)
else lm1 <- lm(rank(r)~firstfactor+secondfactor)
anolm1 <- anova(lm1)
if (withinteraction == 1)  ms <-  anolm1[1:4,1:3]
else ms <-  anolm1[1:3,1:3]
ms[,4] <- ms[,2]/(length(r)*(length(r)+1)/12)
ms[,5] <- (1-pchisq(ms[,4],ms[,1]))
colnames(ms)[4:5] <- c("H", "pvalue")
ms
}