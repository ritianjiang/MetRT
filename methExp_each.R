#this function returns the slope and p
#it judge whether site can impress the exp in one realtion
#using lme
#Wanghaotian
#20171201
methExp_each<-function(idx,relation,meth,exp,sampleinfo){
  sampleinfo$exp<-exp[as.character(relation[idx,]$transId),] %>% as.numeric
  sampleinfo$meth<-meth[as.character(relation[idx,]$loc),] %>% as.numeric
  null.model<-try(lme(exp~gender+age+site,
                      random = ~1|sub,
                      data = sampleinfo,method = "ML"))
  full.model<-try(lme(exp~gender+age+site+meth,
                      random = ~1|sub,
                      data = sampleinfo,method = "ML"))
  fit<-try(anova(null.model,full.model),silent = T)
  if(class(fit) == "try-error"){
    result<-c(NA,NA);
    return(result)} #avoid "$" error
  result<-c(fit$`p-value`[2],full.model$coefficients$fixed[4])
  return(result)
}