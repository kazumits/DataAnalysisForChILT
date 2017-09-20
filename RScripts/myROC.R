# ------
# My ROC
# ------

# prd: predictor (numeric vector)
# ans: answer label (logical vector)
tableTPFP <- function(prd,ans)
{
  ord <- order(prd,decreasing = TRUE)
  # skip ties
  noties <- !duplicated(prd[ord],fromLast=TRUE)
  TP <- cumsum(ans[ord])[noties]
  FP <- seq_along(ans)[noties] - TP
  TPR <- TP/tail(TP,1)
  FPR <- FP/tail(FP,1)
  data.frame(
    threshold = prd[ord][noties],
    FP, TP, FN = sum(ans)-TP, TN = sum(!ans)-FP,
    TPR, FPR, YJ = TPR-FPR, Jaccard = TP/(FP+sum(ans))
    #Dist = sqrt(FPR^2+(1-TPR)^2)
  )
}

labelTPFP <- function(prd,ans) {
  lv <- c("TN","FN","FP","TP")
  factor(lv,lv)[2*prd+ans+1]
}

# Partial apply
papp <- function(f,...) function(x) f(x,...)

# P: predictors (numeric matrix) 
plotROCs <- function(P,ans,pts=TRUE)
{
  require(ggplot2)
  require(dplyr)
  rocs <- apply(P,2,papp(tableTPFP,ans)) %>%
    bind_rows(.id="predictor") %>% 
    mutate(predictor=factor(predictor,colnames(P)))
  
  p <- ggplot(rocs,aes(FPR,TPR,colour=predictor)) +
    #geom_abline(slope=1,lty=1,size=0.3,col="#DDDDDD") +
    geom_line(size=0.4) + 
    scale_x_continuous(breaks = seq(0,1,0.5)) +
    scale_y_continuous(breaks = seq(0,1,0.5)) +
    theme_bw(base_size = 15) #+ coord_fixed()
  if(pts) p <- p + geom_point()
  return(p)
}
