library(brms)
library(cmdstanr)

### I/ read data file ####

test_data = read.csv('./output/dataset mother-infant prox good one.csv')


head(test_data)

nrow(test_data) #2332 scan clusters


### I/ a) explications of the column variable names

colnames(test_data)

# "Mother_ID" : Identity of the mother
# "Infant_ID": Identity of the infant
# "Dyad_ID": Identity of the mother-infant dyad : Please note that the identity of individuals has been anonymised
# "Mother_sociality": the gregariousness of the mother
# "Num_dep_offspr" : Number of dependent offspring that the mother has. 
# "Mother_age" : Age of the mother      
# "Infant_age" : Age of the infant (0 indicates age <0, 1 age 1-2, 2: age 2-3 etc)       
# "Infant_sex" : sex of the infant
# "population" : name of the population (Budongo or Tai)
# "group" : name of the community
# "rdate" : date at which the scan cluster was recorded. 
# "Num_males" :  Number of males in the party    
# "Num_females" : Number of females in the party
# "Scan_clust_ID" : ID of the scan cluster
# "Num_scan"  : Number of scans within each scan clusters
# "Num_scan_in_prox" : number of scan in a scan cluster during which mother and offspring were in proximity
# "p.scan.prox"  : proportion of scans during which mother and offspring were in proximity. 
# "Oustrous_0_1" : presence of oestrous females (0 for NO and 1 for YES)



#### II/ z transform the continuous predictors ####


test_data$Mother_sociality_z = as.numeric(as.vector(scale(test_data$Mother_sociality)))
test_data$z.Num_dep_offspr = as.numeric(as.vector(scale(as.numeric(as.character(test_data$Num_dep_offspr)))))
test_data$z.Mother_age = as.numeric(as.vector(scale(as.numeric(as.character(test_data$Mother_age)))))
test_data$z.Num_males = as.numeric(as.vector(scale(as.numeric(as.character(test_data$Num_males)))))
test_data$z.Num_females = as.numeric(as.vector(scale(as.numeric(as.character(test_data$Num_females)))))
test_data$z.Infant_age = as.numeric(as.vector(scale(as.numeric(as.character(test_data$Infant_age)))))


#### create a scan cluster number ID for the random effect

test_data$Scan_clust_ID = seq(1,nrow(test_data),1)
test_data$population = as.factor(test_data$population)


#### III/ get the priors ####


mprior.full.b = get_prior(Num_scan_in_prox | trials(Num_scan) ~ z.Num_dep_offspr + 
                             z.Mother_age + 
                          z.Infant_age + Infant_sex + Oustrous_0_1 +
                          
                          Mother_sociality_z *z.Num_males*population+
                          Mother_sociality_z * z.Num_females *population+ 
                            (z.Num_males + z.Num_females||rdate) + 
                          ( z.Num_males + z.Num_females||Dyad_ID) +
                          ( z.Num_males + z.Num_females||group)+
                            (1|Scan_clust_ID), 
                        family = binomial, data = test_data)

mprior.full.b$prior


### define the prior to Normal(0,1) for the fixed effects

mprior.full.b$prior[2:17] <- "normal(0,1)"

## Check the model script

make_stancode(Num_scan_in_prox | trials(Num_scan) ~ z.Num_dep_offspr + 
                z.Mother_age + 
                z.Infant_age + Infant_sex + Oustrous_0_1 +
                
                Mother_sociality_z *z.Num_males*population+
                Mother_sociality_z * z.Num_females *population+ 
                (z.Num_males + z.Num_females||rdate) + 
                ( z.Num_males + z.Num_females||Dyad_ID) +
                ( z.Num_males + z.Num_females||group)+
                (1|Scan_clust_ID), 
              family = binomial, data = test_data, 
              prior = mprior.full.b)




### IV/ - run full model ####

full = brm(Num_scan_in_prox | trials(Num_scan) ~ z.Num_dep_offspr + 
             z.Mother_age + 
             z.Infant_age + Infant_sex + Oustrous_0_1 +
             
             Mother_sociality_z *z.Num_males*population+
             Mother_sociality_z * z.Num_females *population+ 
             (z.Num_males + z.Num_females||rdate) + 
             ( z.Num_males + z.Num_females||Dyad_ID) +
             ( z.Num_males + z.Num_females||group)+
             (1|Scan_clust_ID), 
           family = binomial(), data = test_data, 
           prior = mprior.full.b,
           chains = 16, cores = 16, iter = 6000, warmup = 4000,
           control = list(adapt_delta = 0.995), backend = "cmdstanr")



## check posterior predictive check. 

pp_check(full, type = 'hist', ndraw = 10)


### Model results

summary(full)


#get 89% CI

round(posterior_interval(full,prob=0.89),2)


mod.post.full <- posterior_samples(full)

# support for positive estimate

round(sum(mod.post.full[,"b_z.Num_dep_offspr"] < 0.00)/
        length(mod.post.full[,"b_z.Num_dep_offspr"]),3) #0.807

round(sum(mod.post.full[,"b_z.Mother_age"] < 0.00)/
        length(mod.post.full[,"b_z.Mother_age"]),3) #0.711


round(sum(mod.post.full[,"b_z.Infant_age"] < 0.00)/
        length(mod.post.full[,"b_z.Infant_age"]),3) #1


round(sum(mod.post.full[,"b_Infant_sex"] > 0.00)/
        length(mod.post.full[,"b_Infant_sex"]),3) #0.81


round(sum(mod.post.full[,"b_Oustrous_0_1"] > 0.00)/
        length(mod.post.full[,"b_Oustrous_0_1"]),3) #0.987


round(sum(mod.post.full[,"b_Mother_sociality_z"] < 0.00)/
        length(mod.post.full[,"b_Mother_sociality_z"]),3) #0.626

round(sum(mod.post.full[,"b_z.Num_males"] > 0.00)/
        length(mod.post.full[,"b_z.Num_males"]),3) #0.90

round(sum(mod.post.full[,"b_populationTai"] > 0.00)/
        length(mod.post.full[,"b_populationTai"]),3) #0.623

round(sum(mod.post.full[,"b_z.Num_females"] <0.00)/
        length(mod.post.full[,"b_z.Num_females"]),3) #0.968

round(sum(mod.post.full[,"b_Mother_sociality_z:z.Num_males"] > 0.00)/
        length(mod.post.full[,"b_Mother_sociality_z:z.Num_males"]),3) #0.914

round(sum(mod.post.full[,"b_Mother_sociality_z:populationTai"] > 0.00)/
        length(mod.post.full[,"b_Mother_sociality_z:populationTai"]),3) #0.765

round(sum(mod.post.full[,"b_z.Num_males:populationTai"] <0.00)/
        length(mod.post.full[,"b_z.Num_males:populationTai"]),3) #0.904

round(sum(mod.post.full[,"b_Mother_sociality_z:z.Num_females"] < 0.00)/
        length(mod.post.full[,"b_Mother_sociality_z:z.Num_females"]),3) # 0.992

round(sum(mod.post.full[,"b_populationTai:z.Num_females"] > 0.00)/
        length(mod.post.full[,"b_populationTai:z.Num_females"]),3) #0.902

round(sum(mod.post.full[,"b_Mother_sociality_z:z.Num_males:populationTai"] < 0.00)/
        length(mod.post.full[,"b_Mother_sociality_z:z.Num_males:populationTai"]),3) # 0.535

round(sum(mod.post.full[,"b_Mother_sociality_z:populationTai:z.Num_females"] > 0.00)/
        length(mod.post.full[,"b_Mother_sociality_z:populationTai:z.Num_females"]),3) #0.942


### Effect size

#cond.r2
bayes_R2(full, re_formula = NULL, summary = T)

#marginal.r2
bayes_R2(full, re_formula = NA, summary = T)



###########################################

#### V/ Plot the model ####

### rerun the model with centered everything


test_data$Oustrous_0_1.centered = (as.numeric(as.factor(test_data$Oustrous_0_1)))-mean(as.numeric(as.factor(test_data$Oustrous_0_1)))
summary(test_data$Oustrous_0_1.centered)

test_data$Infant_sex.centered = (as.numeric(as.factor(test_data$Infant_sex)))-mean(as.numeric(as.factor(test_data$Infant_sex)))
summary(test_data$Infant_sex.centered)



####### rerun the full model 


mprior.full.plot = get_prior(Num_scan_in_prox | trials(Num_scan) ~ z.Num_dep_offspr + 
                          z.Mother_age + 
                          z.Infant_age + Infant_sex.centered + Oustrous_0_1.centered +
                          
                            Mother_sociality_z *z.Num_males*population+
                            Mother_sociality_z * z.Num_females *population+ 
                          (z.Num_males + z.Num_females||rdate) + 
                          ( z.Num_males + z.Num_females||Dyad_ID) +
                          ( z.Num_males + z.Num_females||group)+
                          (1|Scan_clust_ID), 
                          family = binomial, data = test_data)

mprior.full.plot$prior

mprior.full.plot$prior[2:17] <- "normal(0,1)"





full.plot =  brm(Num_scan_in_prox | trials(Num_scan) ~ z.Num_dep_offspr + 
                     z.Mother_age + 
                     z.Infant_age + Infant_sex.centered + Oustrous_0_1.centered +
                     
                     Mother_sociality_z *z.Num_males*population+
                     Mother_sociality_z * z.Num_females *population+ 
                     (z.Num_males + z.Num_females||rdate) + 
                     ( z.Num_males + z.Num_females||Dyad_ID) +
                     ( z.Num_males + z.Num_females||group)+
                     (1|Scan_clust_ID),
                   
             family = binomial(), data = test_data, 
             prior = mprior.full.plot,
             chains = 16, cores = 16, iter = 6000, warmup = 4000,
             control = list(adapt_delta = 0.995))


summary(full.plot)




###########   V/ a) Define plot parameters N females ####

xx = seq(min(test_data$Mother_sociality_z),max(test_data$Mother_sociality_z), length.out=100)

### check the range of number of female in z transform space

range(test_data$z.Num_females)

# -1.532551  5.047718


##### generate predictive dataframe

pred.data=data.frame(expand.grid(
  z.Num_females = c(-1,0,1),
  Mother_sociality_z = seq(min(test_data$Mother_sociality_z),max(test_data$Mother_sociality_z), 
                           length.out=100),
  population = c("Budongo", "Tai")))


head(pred.data)


est.full.plot = summary(full.plot)$fixed[,"Estimate"]

names(est.full.plot) = rownames(summary(full.plot)$fixed)


yy=seq(0,0,length.out=1000)



## Create some predictive data frame

ests = posterior_samples(full.plot)

m.mat=model.matrix(object=~(Mother_sociality_z * z.Num_females *population), data=pred.data)

est.names = c("b_Intercept",                         
              "b_Mother_sociality_z","b_z.Num_females", "b_populationTai",
              "b_Mother_sociality_z:populationTai",
              "b_Mother_sociality_z:z.Num_females", "b_populationTai:z.Num_females",
              "b_Mother_sociality_z:populationTai:z.Num_females")

ests = ests[, est.names]


colnames(ests) = c("(Intercept)",                         
                   "Mother_sociality_z","z.Num_females", "populationTai",
                   "Mother_sociality_z:populationTai",
                   "Mother_sociality_z:z.Num_females", "z.Num_females:populationTai",
                   "Mother_sociality_z:z.Num_females:populationTai")


which(colnames(m.mat)%in%colnames(ests))

head(m.mat[, colnames(ests)])


ci=lapply(1:nrow(ests), function(x){
  as.vector(m.mat[, colnames(ests)]%*%unlist(ests[x, ]))
})

ci = matrix(unlist(ci), ncol=length(ci), byrow=F)


### put 89% CI
ci = t(apply(ci, 1, quantile, probs=c(0.055, 0.945)))

pred.data = data.frame(pred.data, ci)

head(pred.data)



##### extract the predicted values from the estimates

est.m = data.frame(t(summary(full.plot)$fixed[,"Estimate"]))


colnames(est.m) = rownames(summary(full.plot)$fixed)

est.names2 = c("Intercept", "Mother_sociality_z", "z.Num_females",
               "populationTai", "Mother_sociality_z:z.Num_females", "Mother_sociality_z:populationTai",
               "populationTai:z.Num_females", "Mother_sociality_z:populationTai:z.Num_females")

est.m = est.m[, est.names2]

colnames(est.m) = c("(Intercept)", "Mother_sociality_z", "z.Num_females",
                    "populationTai", "Mother_sociality_z:z.Num_females", "Mother_sociality_z:populationTai",
                    "z.Num_females:populationTai", "Mother_sociality_z:z.Num_females:populationTai")

###generate the predicted values

pred.v = as.vector(m.mat[, colnames(est.m)]%*%unlist(est.m))

pred.data = cbind(pred.data , pred.v)

head(pred.data)

pred.data$pred.v.lin = exp(pred.data$pred.v)/(1+exp(pred.data$pred.v))

pred.data$pred.X5.5 = exp(pred.data$X5.5)/(1+exp(pred.data$X5.5))
pred.data$pred.X94.5 = exp(pred.data$X94.5)/(1+exp(pred.data$X94.5))

#### allocate to each point whether it's low, average or high number of females

test_data$cat.num.f = NA

test_data$cat.num.f [which(test_data$z.Num_females < -0.5)] = 1
test_data$cat.num.f [which(test_data$z.Num_females >= -0.5 & test_data$z.Num_females  <= 0.5)] = 2
test_data$cat.num.f [which(test_data$z.Num_females > 0.5)] = 3

table(test_data$cat.num.f)

# 1   2   3 
# 905 658 769

test_data$cat.num.f = as.numeric(as.character(test_data$cat.num.f))

#######create mother gregariousness categories and allocate each female to the time categories

seq.soc = seq(min(test_data$Mother_sociality_z),max(test_data$Mother_sociality_z), 
              length.out=11)

m.soc = data.frame(matrix(0,ncol=2,nrow=10))

colnames(m.soc)=c("start.soc","end.soc")

m.soc$start.soc=seq.soc[-11]
m.soc$end.soc=seq.soc[-1]
m.soc$gp=(m.soc$start.soc+m.soc$end.soc)/2


#####allocate each datapoint to its gregariousness category

zz = rep(NA,nrow(test_data))

for (i in 1:nrow(test_data)){
  xx=which(m.soc$start.soc < test_data$Mother_sociality_z[i] &
             m.soc$end.soc >= test_data$Mother_sociality_z[i])
  
  if(length(xx)==0){xx=which(m.soc$start.soc <= test_data$Mother_sociality_z[i] &
                               m.soc$end.soc > test_data$Mother_sociality_z[i])}
  
  zz[i]=m.soc$gp[xx]
  
  
}


##calculate the proportion of scan in proximity

test_data$Prox_proportion = test_data$Num_scan_in_prox / test_data$Num_scan

test_data$Prox_proportion = as.numeric(as.character(test_data$Prox_proportion))

head(test_data$Prox_proportion, 10)

test_data$Prox_proportion = as.numeric(test_data$Num_scan_in_prox/ test_data$Num_scan,10)

plot.tab = aggregate(test_data$Prox_proportion, by=list(zz, test_data$population, 
                                                        test_data$Dyad_ID, test_data$cat.num.f),mean)



colnames(plot.tab)=c("z.soc","population", "dyad", "cat.num.fem", "prop.prox")

nrow(plot.tab) # 126

# had sample size
head(plot.tab)

xx = aggregate(test_data$Prox_proportion, by=list(zz, test_data$population, 
                                                  test_data$Dyad_ID, test_data$cat.num.f),length)

plot.tab$samp.size=xx$x


#### average it over all dyads

plot.tab.soc = aggregate(plot.tab$prop.prox, by= list(plot.tab$z.soc, 
                                                      plot.tab$cat.num.fem,plot.tab$population),mean)

colnames(plot.tab.soc) = c("z.soc", "cat.num.fem", "population",  "prop.prox")

head(plot.tab.soc)


#### add the sample size


xx = aggregate(plot.tab$prop.prox, by= list(plot.tab$z.soc, 
                                            plot.tab$cat.num.fem,plot.tab$population),length)

plot.tab.soc$samp.size = xx$x



#############  V/ b) Plot N. females ####

### define colors and arguments for the plot

col = c("lightsalmon", "darkorchid", "steelblue2")
pch = c(16,17,18)

lty = c(1,2,3)

###calculate the values of xx values to appear in the plot

xx= c(seq(0,0.5,0.1))

xx.at = (xx - mean (test_data$Mother_sociality)) / 
  sd(test_data$Mother_sociality)


############## actual plot
library(scales)

windows (25, 15)

par(mfrow = c(1,2), oma= c(1,2,0.5,0.5))


######First Budongo ####

plot.tab.bud = subset(plot.tab.soc, plot.tab.soc$population == "Budongo")

nrow(plot.tab.bud) 


plot(plot.tab.bud$prop.prox ~ 
       jitter(plot.tab.bud$z.soc, factor=1.2), 
     col = alpha(col[plot.tab.bud$cat.num.f], alpha=0.9), pch=pch[plot.tab.bud$cat.num.f], cex= sqrt(plot.tab.bud$samp.size)*1.4, 
     ylab = "Probability to be in proximity", xlab = "Female gregariousness", cex.lab = 1.5, 
     axes=F, ylim=c(0,1))

#axis(1)

#### add axes

axis(1, labels =xx, at = xx.at, cex.axis = 1.2)

axis(2, las=1, cex.axis = 1.2)

axis(1, labels = c("", ""), at = c(-100,100))
axis(2, labels = c("", ""), at = c(-100,100))
axis(3, labels = c("", ""), at = c(-100,100))
axis(4, labels = c("", ""), at = c(-100,100))




####subset to Budongo population to draw model lines and 89% CI

pred.data.b = subset(pred.data, pred.data$population == "Budongo")


for(i in 1:3){
  sel = pred.data.b$z.Num_females==sort(unique(pred.data.b$z.Num_females))[i]
  polygon(x=c(pred.data.b$Mother_sociality_z[sel], rev(pred.data.b$Mother_sociality_z[sel])),
          y=c(pred.data.b[sel, "pred.X5.5"], rev(pred.data.b[sel, "pred.X94.5"])),
          col=adjustcolor(col[i], alpha.f=0.15), border=NA)
  lines(x=pred.data.b$Mother_sociality_z[sel], y=pred.data.b$pred.v.lin[sel],col=col[i],lwd=3, lty=lty[i])
}



mtext("Budongo", side = 3, line = 2, cex= 1.6)



##### then Tai ####



plot.tab.tai = subset(plot.tab.soc, plot.tab.soc$population == "Tai")

nrow(plot.tab.tai)  


plot(plot.tab.tai$prop.prox ~ 
       jitter(plot.tab.tai$z.soc, factor=1.2), 
     col = alpha(col[plot.tab.tai$cat.num.f], alpha=0.9), pch=pch[plot.tab.tai$cat.num.f], cex= sqrt(plot.tab.tai$samp.size)*1.4, 
     ylab = "Probability to be in proximity", xlab = "Female gregariousness", cex.lab = 1.5, 
     axes=F, ylim=c(0,1))

#### add axes
axis(1, labels =xx, at = xx.at, cex.axis = 1.2)

axis(2, las=1, cex.axis = 1.2)

axis(1, labels = c("", ""), at = c(-100,100))
axis(2, labels = c("", ""), at = c(-100,100))
axis(3, labels = c("", ""), at = c(-100,100))
axis(4, labels = c("", ""), at = c(-100,100))


####subset to Tai population to draw model lines and 89% CI

pred.data.t = subset(pred.data, pred.data$population == "Tai")


for(i in 1:3){
  sel = pred.data.t$z.Num_females==sort(unique(pred.data.t$z.Num_females))[i]
  polygon(x=c(pred.data.t$Mother_sociality_z[sel], rev(pred.data.t$Mother_sociality_z[sel])),
          y=c(pred.data.t[sel, "pred.X5.5"], rev(pred.data.t[sel, "pred.X94.5"])),
          col=adjustcolor(col[i], alpha.f=0.15), border=NA)
  lines(x=pred.data.t$Mother_sociality_z[sel], y=pred.data.t$pred.v.lin[sel],col=col[i],lwd=3, lty=lty[i])
}



mtext("Tai", side = 3, line = 2, cex= 1.6)


legend('top', legend= c("Low N. females in the party", "Average N. females in the party", 
                        "High N. females in the party"), cex= 1.3, col= col, lty = lty, pch = pch, lwd =1.1, 
       pt.cex = 2)



###########################################################################################

#### IV/ a) PLOT the age effect ####

### IV/ a) Define plot parameters infant age ####

xx = seq(min(test_data$z.Infant_age),max(test_data$z.Infant_age), length.out=100)

### check the range of numbe rof female in z transfo space

range(test_data$z.Infant_age)

# -1.889133  1.581918


##### generate predictive dataframe

pred.data2=data.frame(expand.grid(
  z.Infant_age = seq(min(test_data$z.Infant_age),max(test_data$z.Infant_age), 
                     length.out=100),
  population = c("Budongo", "Tai")))


head(pred.data2)


est.full.plot = summary(full.plot)$fixed[,"Estimate"]

names(est.full.plot) = rownames(summary(full.plot)$fixed)


yy=seq(0,0,length.out=1000)



## Create some predictive data frame

ests = posterior_samples(full.plot)

m.mat=model.matrix(object=~(z.Infant_age + population), data=pred.data2)

est.names = c("b_Intercept",                         
              "b_z.Infant_age", "b_populationTai")

ests = ests[, est.names]

mean(ests$b_populationTai)

colnames(ests) = c("(Intercept)",                         
                   "z.Infant_age", "populationTai")


which(colnames(m.mat)%in%colnames(ests))

head(m.mat[, colnames(ests)])


ci=lapply(1:nrow(ests), function(x){
  as.vector(m.mat[, colnames(ests)]%*%unlist(ests[x, ]))
})

ci = matrix(unlist(ci), ncol=length(ci), byrow=F)


### put 89% CI
ci = t(apply(ci, 1, quantile, probs=c(0.055, 0.945)))

pred.data2 = data.frame(pred.data2, ci)

head(pred.data2,30)

pred.data2[100:120,]


##### extract the predicted values from the estimates

est.m = data.frame(t(summary(full.plot)$fixed[,"Estimate"]))


colnames(est.m) = rownames(summary(full.plot)$fixed)

est.names2 = c("Intercept", "z.Infant_age", "populationTai")

est.m = est.m[, est.names2]

colnames(est.m) = c("(Intercept)", "z.Infant_age", "populationTai")

###generate the predicted values

pred.v = as.vector(m.mat[, colnames(est.m)]%*%unlist(est.m))

pred.data2 = cbind(pred.data2 , pred.v)

head(pred.data2)

pred.data2$pred.v.lin = exp(pred.data2$pred.v)/(1+exp(pred.data2$pred.v))

pred.data2$pred.X5.5 = exp(pred.data2$X5.5)/(1+exp(pred.data2$X5.5))
pred.data2$pred.X94.5 = exp(pred.data2$X94.5)/(1+exp(pred.data2$X94.5))


### create an age table 

tab.age = aggregate(test_data$Prox_proportion, by =list(test_data$population, test_data$Dyad_ID, test_data$z.Infant_age), mean)

colnames(tab.age) = c("Pop", "Dyad_ID", "z.Infant_age", "prop.prox")

####sample size 

xx = aggregate(test_data$Prox_proportion, by =list(test_data$Dyad_ID, test_data$z.Infant_age), length)

tab.age$sample.size = xx$x

# add sample size
head(tab.age)

str(tab.age)


####    ###  IV b) Plot infant age ####

### define color for the plot
#col = c("blue", "black", "red")

col = c("tomato", "darkturquoise")
pch = c(16,17)

lty = c(1,3)

###calculate the values of xx values to appear in the pplot

xx= c(seq(0,5,1))

xx.at = (xx - mean (test_data$Infant_age)) / 
  sd(test_data$Infant_age)

xx.labels = c("<1", "1-2", "2-3", "3-4", "4-5", "5-6")

############## actual plot
library(scales)

windows (25, 15)

par(oma= c(1,1,0.5,0.5))


plot(tab.age$prop.prox ~ 
       jitter(tab.age$z.Infant_age, factor=0.3), 
     col = alpha(col[as.numeric(as.factor(tab.age$Pop))], alpha=0.9), pch=pch[as.numeric(as.factor(tab.age$Pop))], cex= sqrt(tab.age$sample.size)/4, 
     ylab = "Probability to be in proximity", xlab = "Infant age", cex.lab = 1.3, cex.axis=4, 
     axes=F, ylim=c(0,1))

#axis(1)

#### add nice axes

axis(1, labels =xx.labels, at = xx.at, cex.axis = 1.2)

axis(2, las=1, cex.axis = 1.2)

axis(1, labels = c("", ""), at = c(-100,100))
axis(2, labels = c("", ""), at = c(-100,100))
axis(3, labels = c("", ""), at = c(-100,100))
axis(4, labels = c("", ""), at = c(-100,100))




####subset to Budongo population

for(i in 1:2){
  sel = pred.data2$population==sort(unique(pred.data2$population))[i]
  polygon(x=c(pred.data2$z.Infant_age[sel], rev(pred.data2$z.Infant_age[sel])),
          y=c(pred.data2[sel, "pred.X5.5"], rev(pred.data2[sel, "pred.X94.5"])),
          col=adjustcolor(col[i], alpha.f=0.15), border=NA)
  lines(x=pred.data2$z.Infant_age[sel], y=pred.data2$pred.v.lin[sel],col=col[i],lwd=3, lty=lty[i])
}




legend('topright', legend= c("Budongo", "Tai"), cex= 1.5, col= col, lty = lty, pch = pch, lwd =2, 
       pt.cex = 2)







