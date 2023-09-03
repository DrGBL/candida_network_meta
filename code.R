library("netmeta")
library("gemtc")
library(tidyverse)

data<-read.delim("DataCandida.txt")
dataClasses<-subset(data, !(Study %in% c("Pappas 2007", "Thompson III 2020")))
dataCandidemia<-subset(data, Study != "Kujath 1993")
dataCandidemiaClasses<-subset(dataClasses, !(Study %in% c("Kujath 1993", "Pappas 2007", "Thompson III 2020")))

#sizes: 7x9 landscape

######### plots #########
#individual drugs
p1<-pairwise(treat=data$Drugs, event=data$Outcome, n=data$N, studlab=data$Study)
net1<-netmeta(p1)
sizes<-c(sum(data$N[which(data$Drugs=="Anidulafungin")]), sum(data$N[which(data$Drugs=="Caspofungin")]), sum(data$N[which(data$Drugs=="DAmphoB")]), sum(data$N[which(data$Drugs=="Fluconazole")]), sum(data$N[which(data$Drugs=="Isavuconazole")]), sum(data$N[which(data$Drugs=="LAmphoB")]), sum(data$N[which(data$Drugs=="Micafungin")]), sum(data$N[which(data$Drugs=="Rezafungin")]), sum(data$N[which(data$Drugs=="Voriconazole")]))
sizesRel<-sizes/sum(data$N)
netgraph(net1, number.of.studies = TRUE, plastic=FALSE, col="black", points = TRUE, cex.points = sizesRel*40, thickness="number.of.studies", offset=0.035)

#drug classes
p2<-pairwise(treat=dataClasses$Classes, event=dataClasses$Outcome, n=dataClasses$N, studlab=dataClasses$Study)
net2<-netmeta(p2)
sizes2<-c(sum(dataClasses$N[which(dataClasses$Classes=="Amphotericin")]), sum(dataClasses$N[which(dataClasses$Classes=="Azole")]), sum(dataClasses$N[which(dataClasses$Classes=="Echinocandin")]))
sizesRel2<-sizes2/sum(dataClasses$N)
netgraph(net2, number.of.studies = TRUE, plastic=FALSE, col="black", points = TRUE, cex.points = sizesRel2*40, thickness="number.of.studies", offset=0.07)

######## Survival #########
#individual drugs
net5<-mtc.network(data.ab=data.frame(study=data$Study, treatment=data$Drugs, sampleSize=data$N, responders=(data$N-data$Mortality)))
plot(net5)

mod5<-mtc.model(net5, likelihood="binom", link="logit")
run5<-mtc.run(mod5)
summary(run5)
forest(run5)
forest(relative.effect(run5, t1=c("DAmphoB")))

rank.probability(run5)
plot(rank.probability(run5))

round(exp(relative.effect.table(run5)),3)

#drug classes
net6<-mtc.network(data.ab=data.frame(study=dataClasses$Study, treatment=dataClasses$Classes, sampleSize=dataClasses$N, responders=(dataClasses$N-dataClasses$Mortality)))
plot(net6)

mod6<-mtc.model(net6, likelihood="binom", link="logit")
run6<-mtc.run(mod6)
summary(run6)
forest(run6)
forest(relative.effect(run6, t1=c("Amphotericin")))

rank.probability(run6)
plot(rank.probability(run6))

round(exp(relative.effect.table(run6)),3)

######## Primary Outcome ###########
net<-mtc.network(data.ab=data.frame(study=data$Study, treatment=data$Drugs, sampleSize=data$N, responders=data$Outcome))
plot(net)

mod<-mtc.model(net, likelihood="binom", link="logit")
run<-mtc.run(mod)
summary(run)
forest(run)
forest(relative.effect(run, t1=c("DAmphoB")))

rank.probability(run)
plot(rank.probability(run))

round(exp(relative.effect.table(run)),3)

#drug classes
net2<-mtc.network(data.ab=data.frame(study=dataClasses$Study, treatment=dataClasses$Classes, sampleSize=dataClasses$N, responders=dataClasses$Outcome))
plot(net2)

mod2<-mtc.model(net2, likelihood="binom", link="logit")
run2<-mtc.run(mod2)
summary(run2)
forest(run2)
forest(relative.effect(run2, t1=c("Amphotericin")))

rank.probability(run2)
plot(rank.probability(run2))

round(exp(relative.effect.table(run2)),3)

######## Candidemia primary Outcome #########
#individual drugs
net3<-mtc.network(data.ab=data.frame(study=dataCandidemia$Study, treatment=dataCandidemia$Drugs, sampleSize=dataCandidemia$Candidemia, responders=dataCandidemia$OutcomeCandidemia))
plot(net3)

mod3<-mtc.model(net3, likelihood="binom", link="logit")
run3<-mtc.run(mod3)
summary(run3)
forest(run3)
forest(relative.effect(run3, t1=c("DAmphoB")))

rank.probability(run3)
plot(rank.probability(run3))

round(exp(relative.effect.table(run)),3)

#drug classes
net4<-mtc.network(data.ab=data.frame(study=dataCandidemiaClasses$Study, treatment=dataCandidemiaClasses$Classes, sampleSize=dataCandidemiaClasses$Candidemia, responders=dataCandidemiaClasses$OutcomeCandidemia))
plot(net4)

mod4<-mtc.model(net4, likelihood="binom", link="logit")
run4<-mtc.run(mod4)
summary(run4)
forest(run)
forest(relative.effect(run4, t1=c("Amphotericin")))

rank.probability(run4)
plot(rank.probability(run4))

round(exp(relative.effect.table(run4)),3)



