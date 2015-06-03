
### Problem 1
util <- function(r) .62*log(.004*r+1)
exp.util.1a <- util(500)/3
exp.util.1b <- 2*util(-100)/3 + util(400)/3

### Problem 3
brier <-  function(p,act) (p-act)^2 + (-p+act)^2
lloss <-  function(p,act) -act*log(p) - (1-act)*log(1-p)
pred1 <-c(.59,.79,.54,.36,.77,.62,.37,.18,.39,.14)
pred2 <-c(.32,.57,.60,.53,.85,.21,.46,.07,.45,.29)
actual <- c(1,1,0,0,1,1,1,0,1,0)

brier1.3 <- sum(brier(pred1,actual))/10
brier2.3 <- sum(brier(pred2,actual))/10
lloss1.3 <- sum(lloss(pred1,actual))/10 
lloss2.3 <- sum(lloss(pred2,actual))/10 


