library(hier.part)


data(urbanwq)
env <- urbanwq[,2:8]
all.regs(urbanwq$lec, env, fam = "gaussian", gof = "Rsqu",
         print.vars = TRUE)


hier.part(urbanwq$lec, env, fam = "gaussian", gof = "Rsqu")
data(amphipod)
env1 <- amphipod[,2:5]
all.regs(amphipod$australis, env1, fam = "binomial",
         gof = "logLik", print.vars = TRUE)


data(urbanwq)
env <- urbanwq[,2:8]
gofs <- all.regs(urbanwq$lec, env, fam = "gaussian",
                 gof = "Rsqu", print.vars = TRUE)
partition(gofs, pcan = 7, var.names = names(urbanwq[,2,8]))

library(relaimpo)
data(swiss)
calc.relimp(swiss,
            type = c("lmg", "last", "first", "betasq", "pratt", "genizi", "car") )
# calculation of all available relative importance metrics
# non-US version offers the additional metric "pmvd",
# i.e. call would be
# calc.relimp(cov(swiss),
# type = c("lmg", "pmvd", "last", "first", "betasq, "pratt"),
# rela = TRUE )


data(urbanwq)
env <- urbanwq[,2:8]
all.regs(urbanwq$lec, env, fam = "gaussian", gof = "Rsqu",
         print.vars = TRUE)
#logistic regression of an amphipod species occurrence in
#streams against four independent variables describing
#catchment characteristics (from Walsh et al. 2004)



data(amphipod)
env1 <- amphipod[,2:5]
gofs<-all.regs2(amphipod$australis, env1, fam = "binomial",
         gof = "Rsqu.glm", print.vars = TRUE)
gofs
partition(gofs, pcan = 4, var.names = names(amphipod[,2:5]))



data(amphipod)
env1 <- amphipod[,2:5]
gof<-all.regs(amphipod$australis, env1, fam = "binomial",
         gof = "logLik", print.vars = TRUE)

partition(gof, pcan = 4, var.names = names(amphipod[,2:5]))












1/(mean(amphipod$australis)*(1-mean(amphipod$australis)))


fit1<-glm(australis~fimp+fconn+densep+unseal, data = amphipod, family = binomial)
summary(fit1)
VarF <- var(plogis(as.vector(model.matrix(fit1) %*% coef(fit1))))

VarF <- var((as.vector(model.matrix(fit1) %*% coef(fit1))))

M=matrix(0, 100000)
for (i in 1:100000){
  M[i]=(sum(rbinom(length(amphipod$australis), 1,plogis(predict(fit1))))/length(amphipod$australis))
}
M*length(amphipod$australis)
var(M)

VarDS<-var(M)*length(amphipod$australis)


gf <- VarF/(VarF + VarDS)
gf

fit1<-glmer(australis~fimp+fconn+densep+unseal, data = amphipod, family = binomial)

summary(fit1)
r2beta(model = fit1, method = 'sgv', data = amphipod)

summary(fit1)


fit1<-glmer(australis~fimp+fconn+densep, data = amphipod, family = binomial)


p<-plogis(model.matrix(fit1)%*%coef(fit1))
v<-p*(1-p)
v.o<-1/v
sum(v.o)/length(v.o)


fit1<-glm(australis~fimp+fconn+densep+unseal, data = amphipod, family = binomial(link='probit'))
summary(fit1)
VarF <- var(as.vector(model.matrix(fit1) %*% coef(fit1)))
VarDS <- 4.539811
gf <- VarF/(VarF + 1)

fit1<-glm(australis~fimp+fconn+densep+unseal, data = amphipod, family = binomial(link='logit'))
summary(fit1)
VarF <- var(as.vector(model.matrix(fit1) %*% coef(fit1)))
VarDS <- pi^2/3
gf <- VarF/(VarF + VarDS)
gf





fit1<-glm(australis~fimp+unseal+densep+fconn, data = amphipod, family = binomial(link='logit'))
summary(fit1)
VarF <- var(as.vector(model.matrix(fit1) %*% coef(fit1)))
VarDS <- pi^2/3
gf <- VarF/(VarF + VarDS)
gf
fit1<-glm(australis~densep+fimp, data = amphipod, family = binomial(link='logit'))


mod.pql = glmPQL(fit1)

r2beta(model = mod.pql, method = 'svg', data = amphipod)



#get exact predictions also of random effects then to VarF/(VarFuR+SimulatedVar)

obs<-seq(nrow(amphipod))

fit1<-glmer(australis~densep+fimp+fconn+unseal+(1|obs), data = amphipod, family = binomial(link='logit'))
summary(fit1)


M=matrix(0, 10000)
for (i in 1:10000){
  M[i]=qlogis(sum(rbinom(length(amphipod$australis), 1,plogis(predict(fit1))))/length(amphipod$australis))
  #M[i]=sum(rbinom(length(y), 1,plogis(predict(mod))))/length(y)
}

VarDS<-var(M)*length(amphipod$australis)

VarF <- var(as.vector(model.matrix(fit1) %*% fixef(fit1)))

VarFF<-var(predict(fit1, type='response'))

VarDS <- 4.539811
gf <- VarF/(VarF + 1)

gf <- VarF/(VarF + VarDS+sum(VarCorr(fit1)$obs))
VarF <- var(plogis((model.matrix(fit1) %*% fixef(fit1))))
VarF



