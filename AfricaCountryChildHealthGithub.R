## Africa country child health analysis
## Corey J.A. Bradshaw
## September 2019
## 
## Accompanies the paper: Bradshaw, CJA, SP Otto, Z Mehrabi, AA Annamalay, S Heft-Neal, Z Wagner, PN Le Souëf. 2019.
## Testing the socio-economic and environmental determinants of better child-health outcomes in Africa: a cross-sectional
## study among nations. BMJ Open In press

## remove everything
rm(list = ls())

## libraries
library(boot)
library(biomod2)
library(sem)
library(semPlot)
library(semGOF)
library(dismo)
library(gbm)
library(ape)
library(geosphere)
library(ncf)
library(ModelMetrics)
library(pgirmess)
library(lme4)



# background data for sorting
contreg <- read.table("continent.country.csv", sep=",", header=T)
contnum <- read.table("cntry.num.code.csv", sep=",", header=T)

## regions of Africa
regions.AFR <- read.table("Africa.regions.csv", sep=",", header=T)

## remove small-island nations (except Madagascar): Cabo Verde, Seychelles, Réunion, Comoros, Mauritius, São Tomé & Príncipe, Mayotte, Saint Helena-Ascension-Tristan da Cunha
smisnat <- c("CPV", "SYC", "REU", "COM", "MUS", "STP", "MYT", "SHN") 

## country code data for WHO datasets
WHO.cntry <- read.table("WHO.cntry.csv", sep=",", header=T)
colnames(WHO.cntry)[1] <- "country"

## country codes for FAO data
FAO.cntry <- read.table("faocountrycode.csv", sep=",", header=T)

## PPP
ppp <- read.table("PPP11.15.csv", sep=",", header=T)

## ENVIRONMENT
# Ecological Footprint consumption per capita
EFcons.PC <- read.table("EFconsPC12.csv", sep=",", header=T)
EFcons.reg <- merge(EFcons.PC, contreg, by="cntry.code")
AFR.EFcons.reg <- subset(EFcons.reg, cont=="AFR")
AFR.EFcons.sort <- AFR.EFcons.reg[order(AFR.EFcons.reg[,8], decreasing=F),]

# Megafauna Conservation Index
MCI <- read.table("megafaunaconserv.csv", sep=",", header=T)
mci.reg <- merge(MCI, contreg, by="cntry.code")
AFR.mci.reg <- subset(mci.reg, cont=="AFR")
AFR.mci.sort <- AFR.mci.reg[order(AFR.mci.reg[,8], decreasing=T),]

# proportion of land area protected (Population Reference Bureau)
landprot <- read.table("landprot.csv", sep=",", header=T)
landprot.reg <- merge(landprot, contreg, by="cntry.code")
AFR.landprot.reg <- subset(landprot.reg, cont=="AFR")
AFR.landprot.sort <- AFR.landprot.reg[order(AFR.landprot.reg[,4], decreasing=F),]

# Threatened species (IUCN 2016 data ([CR+EN+VU+NT]/[CR+EN+VU+NT+LC]))
AFR.threat <- read.table("AFR.threat.csv", sep=",", header=T)
AFR.threat$prop.threat <- apply(AFR.threat[,c(5:8)], 1, sum, na.rm=T)/apply(AFR.threat[,c(5:8,11)], 1, sum, na.rm=T)
AFR.threat.reg <- merge(AFR.threat, contreg, by="cntry.code")
AFR.threat.sort <- AFR.threat.reg[order(AFR.threat.reg[,13], decreasing=F),]

# freshwater removals (% of internal resources)
freshwrem <- read.table("freshwrem.csv", sep=",", header=T)
freshwrem.reg <- merge(freshwrem, contreg, by="cntry.code")
AFR.freshwrem.reg <- subset(freshwrem.reg, cont=="AFR")
AFR.freshwrem.sort <- AFR.freshwrem.reg[order(AFR.freshwrem.reg[,2], decreasing=F),]

# forest loss 2000 - 2012 (Hansen Science 2013)
forestloss <- read.table("forestloss.csv", sep=",", header=T)
forestloss.reg <- merge(forestloss, contreg, by="cntry.code")
landarea <- read.table("land.area.csv", sep=",", header=T)
forestlossla.reg <- merge(forestloss.reg, landarea, by="cntry.code")
AFR.forestloss.reg <- subset(forestlossla.reg, cont=="AFR")
AFR.forestloss.reg$netloss.area <- (AFR.forestloss.reg$tot.loss - AFR.forestloss.reg$tot.gain)/AFR.forestloss.reg$land.area15
AFR.forestloss.sort <- AFR.forestloss.reg[order(AFR.forestloss.reg[,18], decreasing=F),]

# livestock (cattle, pigs, buffaloes, sheep, goats) per ha of arable land
livestock <- read.table("livestock.csv", sep=",", header=T)
livestock.tot <- xtabs(livestock$no.ha ~ livestock$country)
livestock.tot.dat <- data.frame(names(livestock.tot), as.numeric(livestock.tot))
colnames(livestock.tot.dat) <- c("country", "livestock.tot")
cntry.unique <- as.data.frame(unique(livestock$country))
colnames(cntry.unique) <- "country"
cntry.codes <- (livestock[, c(1,2)])
cntry.code.unique <- unique(merge(cntry.unique, cntry.codes, by="country", all=F))
livestock.tot.out <- merge(livestock.tot.dat, cntry.code.unique, by="country")
livestock.tot.reg <- merge(livestock.tot.out, contreg, by="cntry.code")
AFR.livestock <- subset(livestock.tot.reg, cont=="AFR")
AFR.livestock.sort <- AFR.livestock[order(AFR.livestock[,3], decreasing=F),]

## arable land (to estimate total livestock numbers)
## ha arable land in 2014
arable <- read.table("arableland.csv", sep=",", header=T)
arable.reg <- merge(arable, contreg, by="cntry.code")
AFR.arable <- subset(arable.reg, cont="AFR")
AFR.livestock.num <- merge(AFR.livestock.sort, AFR.arable)
AFR.livestock.num$livestock.num <- round(AFR.livestock.num$arable2014 * AFR.livestock.num$livestock.tot, 0)
AFR.livestock.final <- AFR.livestock.num[!(AFR.livestock.num$cntry.code %in% smisnat),]

# permanent cropland (% of total land area; World Bank)
cropland <- read.table("cropland.csv", sep=",", header=T)
cropland.reg <- merge(cropland, contreg, by="cntry.code")
AFR.cropland.reg <- subset(cropland.reg, cont=="AFR")
AFR.cropland.sort <- AFR.cropland.reg[order(AFR.cropland.reg[,2], decreasing=F),]

# greenhouse gas emissions (CO2e 2013 per capita) World Bank ()
emissions <- read.table("pcemiss.csv", sep=",", header=T)
emissions.reg <- merge(emissions, contreg, by="cntry.code")
AFR.emissions.reg <- subset(emissions.reg, cont=="AFR")
AFR.emissions.sort <- AFR.emissions.reg[order(AFR.emissions.reg[,2], decreasing=F),]


# rank, sort, combine, remove small-island nations
# ecological footprint
AFR.EFcons.list <- AFR.EFcons.reg[!(AFR.EFcons.reg$cntry.code %in% smisnat),]
AFR.EFcons.list$EFrank <- rank(AFR.EFcons.list$EFconsPC12, na.last="keep", ties.method="average")
EFrank.dat <- AFR.EFcons.list[,c(1,13)]
# megafauna conservation index
AFR.mci.list <- AFR.mci.reg[!(AFR.mci.reg$cntry.code %in% smisnat),]
AFR.mci.list$MCIrank <- rank(-AFR.mci.list$MCI, na.last="keep", ties.method="average")
MCIrank.dat <- AFR.mci.list[,c(1,15)]
# threatened species
AFR.thr.list <- AFR.threat.reg[!(AFR.threat.reg$cntry.code %in% smisnat),]
AFR.thr.list$THRrank <- rank(AFR.thr.list$prop.threat, na.last="keep", ties.method="average")
THRrank.dat <- AFR.thr.list[,c(1,18)]
# freshwater removal
AFR.fwr.list <- AFR.freshwrem.reg[!(AFR.freshwrem.reg$cntry.code %in% smisnat),]
AFR.fwr.list$FWRrank <- rank(AFR.fwr.list$freshw.rem14, na.last="keep", ties.method="average")
FWRrank.dat <- AFR.fwr.list[,c(1,7)]
# forest loss
AFR.fsl.list <- AFR.forestloss.reg[!(AFR.forestloss.reg$cntry.code %in% smisnat),]
AFR.fsl.list$FSLrank <- rank(AFR.fsl.list$netloss.area, na.last="keep", ties.method="average")
FSLrank.dat <- AFR.fsl.list[,c(1,20)]
# livestock
AFR.lvs.list <- AFR.livestock.sort[!(AFR.livestock.sort$cntry.code %in% smisnat),]
AFR.lvs.list$LVSrank <- rank(AFR.lvs.list$livestock.tot, na.last="keep", ties.method="average")
LVSrank.dat <- AFR.lvs.list[,c(1,8)]
# cropland
AFR.crp.list <- AFR.cropland.sort[!(AFR.cropland.sort$cntry.code %in% smisnat),]
AFR.crp.list$CRPrank <- rank(AFR.crp.list$cropland14, na.last="keep", ties.method="average")
CRPrank.dat <- AFR.crp.list[,c(1,7)]
# emissions
AFR.emi.list <- AFR.emissions.sort[!(AFR.emissions.sort$cntry.code %in% smisnat),]
AFR.emi.list$EMIrank <- rank(AFR.emi.list$pcemiss13, na.last="keep", ties.method="average")
EMIrank.dat <- AFR.emi.list[,c(1,7)]

# merge
ranks.dat1 <- merge(EFrank.dat, MCIrank.dat, by="cntry.code", all.x=T, all.y=T)
ranks.dat2 <- merge(ranks.dat1, THRrank.dat, by="cntry.code", all.x=T, all.y=T)
ranks.dat3 <- merge(ranks.dat2, FWRrank.dat, by="cntry.code", all.x=T, all.y=T)
ranks.dat4 <- merge(ranks.dat3, FSLrank.dat, by="cntry.code", all.x=T, all.y=T)
ranks.dat5 <- merge(ranks.dat4, LVSrank.dat, by="cntry.code", all.x=T, all.y=T)
ranks.dat6 <- merge(ranks.dat5, CRPrank.dat, by="cntry.code", all.x=T, all.y=T)
ranks.dat <- merge(ranks.dat6, EMIrank.dat, by="cntry.code", all.x=T, all.y=T)

# remove countries with < 6 metrics
cntry.vec <- names(table(ranks.dat$cntry.code))
lcntr <- length(cntry.vec)
flag.kp <- rep(0,lcntr)
for (i in 1:lcntr) {
  cntry.dat <- subset(ranks.dat, cntry.code == cntry.vec[i])
  flag.kp[i] <- ifelse(length(which(is.na(cntry.dat[,c(2:9)])==F)) < 6, 0, 1)
}
flag.cntry <- data.frame(cntry.vec, flag.kp)
colnames(flag.cntry)[1] <- "cntry.code"
ranks.dat.mrg <- merge(ranks.dat, flag.cntry, by="cntry.code", all.x=T, all.y=T)
ranks.dat.enough <- subset(ranks.dat.mrg, flag.kp==1)

# median rank
ranks.dat.enough$median.rank <- apply(ranks.dat.enough[,2:9], 1, median, na.rm=T)

# geometric mean rank
ranks.dat.enough$geom.rank <- scale(10^(apply(log10(ranks.dat.enough[,2:9]), 1, mean, na.rm=T)), center=T, scale=T)
plot(ranks.dat.enough$median.rank, ranks.dat.enough$geom.rank, pch=19)

# sort
ranks.dat.sort <- ranks.dat.enough[order(ranks.dat.enough[,12], decreasing=F),]
ranks.dat.sort
env.dat <- ranks.dat.sort[,c(1,12)]
colnames(env.dat)[2] <- c("env.grank")

## corelation matrix of component parts
env.cor.mat <- cor(na.omit(ranks.dat.sort[,c(2:9)]), method="kendall")
env.cor.mat


## CHILD HEALTH METRICS
child.health <- read.table("childhealthmetrics.csv", sep=",", header=T)
table(child.health$year)
child.health15 <- subset(child.health, year==2015)
child.health15.u5 <- subset(child.health15, agegrp=="months1.59")
child.health15.reg <- merge(child.health15.u5, contreg, by="cntry.code")
child.health15.AFR <- subset(child.health15.reg, cont=="AFR")

# stunting
# % children under 5 considered underweight
pcu5.stunt <- read.table("pcu5stunt.csv", sep=",", header=T)
pcu5stunt <- merge(pcu5.stunt, WHO.cntry, by="country")
# most up-to-date year
cntr.vec <- as.vector(unique(names(table(pcu5stunt$cntry.code))))
lcntr <- length(cntr.vec)
pcu5stunt.max <- pcu5stunt[1,]
for (i in 1:lcntr) {
  cntry.dat <- subset(pcu5stunt, cntry.code == cntr.vec[i])
  cntry.max <- subset(cntry.dat, year == max(cntry.dat$year))
  pcu5stunt.max <- rbind(pcu5stunt.max,cntry.max)
}
pcu5stunt.max <- pcu5stunt.max[-1, ]
pcu5stunt.reg <- merge(pcu5stunt.max, contreg, by="cntry.code")
pcu5stunt.sort <- pcu5stunt.reg[order(pcu5stunt.reg[,4], decreasing=F),]
pcu5stunt.AFR <- subset(pcu5stunt.sort, cont=="AFR")

stunt <- data.frame(pcu5stunt.AFR$cntry.code, pcu5stunt.AFR$pcu5stunt.both)
colnames(stunt) <- c("cntry.code", "stuntingu5")

# deaths from respiratory disease (deaths/1000 live births for acute lower respiratory infection; WHO)
resp.dis <- read.table("respiratorydisease.csv", sep=",", header=T)
resp.dis15 <- subset(resp.dis, year==2015)
respdisp15 <- merge(resp.dis15, WHO.cntry, by="country")
respdisp.reg <- merge(respdisp15, contreg, by ="cntry.code")
respdisp.AFR <- subset(respdisp.reg, cont=="AFR")

# diarrhoea (under 5s)
diarrh.dat <- subset(child.health15.AFR, cause=="Diarrhoea")
diarrh.mort <- data.frame(diarrh.dat$cntry.code, diarrh.dat$mort)
colnames(diarrh.mort) <- c("cntry.code", "diarrhoea.mort")

# malaria/AIDS/pneum/measles/pert/TB
  # malaria
  malaria.dat <- subset(child.health15.AFR, cause=="Malaria")
  malaria.deaths <- data.frame(malaria.dat$cntry.code, malaria.dat$ndeaths, malaria.dat$livebirths)
  colnames(malaria.deaths) <- c("cntry.code", "mal.deaths", "livebirths")
  
  # AIDS
  AIDS.dat <- subset(child.health15.AFR, cause=="AIDS")
  AIDS.deaths <- data.frame(AIDS.dat$cntry.code, AIDS.dat$ndeaths)
  colnames(AIDS.deaths) <- c("cntry.code", "AIDS.deaths")

  # pneumonia (remove? already included in RESPIRATORY INFECTION)
  pneum.dat <- subset(child.health15.AFR, cause=="Pneumonia")
  pneum.deaths <- data.frame(pneum.dat$cntry.code, pneum.dat$ndeaths)
  colnames(pneum.deaths) <- c("cntry.code", "pneum.deaths")
  
  # pertussis
  pert.dat <- subset(child.health15.AFR, cause=="Pertussis")
  pert.deaths <- data.frame(pert.dat$cntry.code, pert.dat$ndeaths)
  colnames(pert.deaths) <- c("cntry.code", "pert.deaths")
  
#infectious total
infect1 <- merge(malaria.deaths, AIDS.deaths, by="cntry.code")
infect.dis <- merge(infect1, pert.deaths, by="cntry.code")
rm(infect1)
infect.dis$dismort <- 100 * apply(infect.dis[,c(2,4,5)], 1, sum, na.rm=T) / infect.dis$livebirths

# injuries (1-59 months)
injury.dat <- subset(child.health15.AFR, cause=="Injuries")
injury.mort <- data.frame(injury.dat$cntry.code, injury.dat$mort)
colnames(injury.mort) <- c("cntry.code", "injury.mort")


## combined health rank
chealth.comb1 <- merge(pcu5stunt.AFR[,c(1,4)], respdisp.AFR[,c(1,4)], by = "cntry.code")
chealth.comb2 <- merge(chealth.comb1, diarrh.mort, by = "cntry.code")
chealth.comb3 <- merge(chealth.comb2, infect.dis[,c(1,6)], by = "cntry.code")
chealth.comb4 <- merge(chealth.comb3, injury.mort, by = "cntry.code")
chealth.comb <- chealth.comb4[!(chealth.comb4$cntry.code %in% smisnat),]
chealth.comb$stunt.rnk <- rank(-(chealth.comb$pcu5stunt.both), na.last="keep", ties.method="average")
chealth.comb$resp.rnk <- rank(-chealth.comb$respdismonths1.59, na.last="keep", ties.method="average")
chealth.comb$diar.rnk <- rank(-chealth.comb$diarrhoea.mort, na.last="keep", ties.method="average")
chealth.comb$inf.rnk <- rank(-chealth.comb$dismort, na.last="keep", ties.method="average")
chealth.comb$inj.rnk <- rank(-chealth.comb$injury.mort, na.last="keep", ties.method="average")
chealth.comb$chealth.grank <- 10^(apply(log10(chealth.comb[,7:11]), 1, mean, na.rm=T))
chealth.comb$st.stunt <- scale(chealth.comb$pcu5stunt.both, center=T, scale=T)
chealth.comb$st.resp <- scale(chealth.comb$respdismonths1.59, center=T, scale=T)
chealth.comb$st.diar <- scale(chealth.comb$diarrhoea.mort, center=T, scale=T)
chealth.comb$st.inf <- scale(chealth.comb$dismort, center=T, scale=T)
chealth.comb$st.inj <- scale(chealth.comb$injury.mort, center=T, scale=T)
chealth.comb$chealth.mstanmn <- scale(apply(chealth.comb[,c(13:17)], 1, mean, na.rm=T), center=T, scale=T)
chealth.comb$chealth.mstan <- scale(chealth.comb$chealth.grank, center=T, scale=T)
chealth.out <- chealth.comb[,c(1,19)]
pairs(chealth.comb[,c(13:17,19)], pch=19, cex=0.7)

# mean vs. geometric mean rank
plot(chealth.comb$chealth.mstanmn, chealth.comb$chealth.mstan, pch=19, xlab="mean", ylab="geometric mean rank")
cor(x=chealth.comb$chealth.mstanmn, y=chealth.comb$chealth.mstan, method="spearman")

## corelation matrix of component parts
chlth.cor.mat <- cor(na.omit(chealth.comb[,c(7:11)]), method="kendall")
chlth.cor.mat

# DALY5
  # life expectancy
  lifeexp.f.dat <- read.table("life.exp.birth.fem.csv", sep=",", header=T)
  lifeexp.m.dat <- read.table("life.exp.birth.mal.csv", sep=",", header=T)
  lifeexp.dat <- merge(lifeexp.f.dat, lifeexp.m.dat, by="cntry.code")
  lifeexp.dat$lifexp.mn <- apply(lifeexp.dat[,c(3,6)], 1, mean, na.rm=T)
  lifeexp.reg <- merge(lifeexp.dat, contreg, by="cntry.code")
  lifeexp.AFR <- subset(lifeexp.reg, cont=="AFR")
  DALY5 <- lifeexp.AFR[,c(1,8)]
  DALY5$DALY5 <- DALY5$lifexp.mn - 5
  DALY5.AFR <- DALY5[!(DALY5$cntry.code %in% smisnat),]

  
    
## CORRELATES
# population density
popN <- read.table("pop.yr.csv", sep=",", header=T)
popN.land <- merge(popN, landarea, by="cntry.code")
popN.reg <- merge(popN.land, contreg, by="cntry.code")
popN.reg$popD15 <- popN.reg$X2015/popN.reg$land.area15
AFR.popN.reg <- subset(popN.reg, cont=="AFR")
AFR.popN.reg$popDrnk <- rank(AFR.popN.reg$popD15, na.last="keep", ties.method="average")
AFR.popN.reg$lPOPD <- scale(log10(AFR.popN.reg$popD15), center=T, scale=T)

# household size instead (https://population.un.org/Household/index.html#/countries)
house.size <- read.table("household size Africa.csv", sep=",", header=T)
AFR.popN.reg2 <- merge(house.size,AFR.popN.reg,by="cntry.code", all.x = T, all.y = T)
AFR.popN.reg2$hsize.sc <- scale(log10(AFR.popN.reg2$hsize.avg), center=T, scale=T)
popD.AFR <- AFR.popN.reg2[,c(1,67,69)] # choose household size instead
colnames(popD.AFR) <- c("cntry.code","popDrnk","lPOPD")

## examine relationship between household size and population density
hist(AFR.popN.reg2$lPOPD)
hist(AFR.popN.reg2$hsize.sc)
plot(AFR.popN.reg2$lPOPD, AFR.popN.reg2$hsize.sc, pch=19, ylab="household size", xlab="population density")
fit.HS.PD <- lm(AFR.popN.reg2$hsize.sc ~ AFR.popN.reg2$lPOPD)
abline(fit.HS.PD,col="red",lty=2)
cor(na.omit(AFR.popN.reg2[,c(67,68)]))

# GDP
gdppcppp <- read.table("gdppcppp.csv", sep=",", header=T) # income share of lowest 10%
gdppcppp.reg <- merge(gdppcppp, contreg, by="cntry.code")
AFRgdppcppp <- subset(gdppcppp.reg, cont=="AFR")
AFR.gdppcppp <- AFRgdppcppp[!(AFRgdppcppp$cntry.code %in% smisnat),]
AFR.gdppcppp$GDP.rnk <- rank(-AFR.gdppcppp$gdppcppp.1115, na.last="keep", ties.method="average")
AFR.gdppcppp$lGDP <- scale(log10(AFR.gdppcppp$gdppcppp.1115), center=T, scale=T)
GDP.AFR <- AFR.gdppcppp[,c(1,7,8)]

# poverty gap (Poverty gap at national poverty lines is the mean shortfall 
# from the poverty lines (counting the nonpoor as having zero shortfall) as 
# a percentage of the poverty lines. This measure reflects the depth of poverty as well as its incidence)
povgap <- read.table("povgap.csv", sep=",", header=T)
povgap.reg <- merge(povgap, contreg, by="cntry.code")
AFRpovgap <- subset(povgap.reg, cont=="AFR")
AFR.povgap <- AFRpovgap[!(AFRpovgap$cntry.code %in% smisnat),]
AFR.povgap$POVGAP.rnk <- rank(AFR.povgap$povgap0615, na.last="keep", ties.method="average")
povgap.AFR <- AFR.povgap[,c(1,7)]

## Ibrahim Index of African Governance
## Money, Politics & Transparency
# set working directory
iiag.dat <- read.table("overallgovernance.csv", sep=",", header=T)
iiag.2015 <- subset(iiag.dat, year==2015)
iiag.2015$cntr.code2 <- as.character(iiag.2015$cntr.code2)
iiag.2015$cntr.code2[which(iiag.2015$country == "Namibia")] <- c("NM") # fix 2-code Namibia from NA to NM
cc2to3 <- read.table("cc2to3.csv", sep=",", header=T)
iiag <- merge(iiag.2015, cc2to3, by="cntr.code2", all.x=T, all.y=T)
iiag.list <- iiag[!(iiag$cntry.code %in% smisnat),]
iiagsort <- iiag.list[order(iiag.list[,4], decreasing=T),]
iiag.sort <- iiagsort[!(iiagsort$cntry.code %in% smisnat),]
iiag.sort$lGOV <- scale((iiag.sort$governance/100), center=T, scale=T)
iiag.sort$GOV.rnk <- rank(-iiag.sort$governance, na.last="keep", ties.method="average")

# Health investment
  # aid (Institute for Health Metrics and Evaluation IHME DAH Database)
  aid.dat <- read.table("IHME_DAH_DATABASE_1990_2015_Y2016M04D25.csv", sep=",", header=T)
  colnames(aid.dat)[4] <- c("cntry.code")
  aid.dec <- subset(aid.dat, year > 2005) # health aid data from last decade
  aid.mrg <- merge(aid.dec, contreg, by="cntry.code")
  aid.AFR <- subset(aid.mrg, cont=="AFR")
  aid.tot.AFR.tmp <- aid.AFR[,c(1,2,4,13,15)]
  aid.tot.AFR <- subset(aid.tot.AFR.tmp, elim_ch == 0) 
  aid.tot.AFR$cntry.code <- as.factor(aid.tot.AFR$cntry.code)
  aid.sum.AFR.cntry <- names(xtabs(dah_15 ~ cntry.code, data=aid.tot.AFR))
  aid.sum.AFR.sum <- as.numeric(xtabs(dah_15 ~ cntry.code, data=aid.tot.AFR))
  aid.sum.AFR <- data.frame(aid.sum.AFR.cntry, aid.sum.AFR.sum)
  colnames(aid.sum.AFR) <- c("cntry.code", "aid.sum")
  aid.sum.AFR.mrg <- merge(aid.sum.AFR, contreg, by="cntry.code")  
  aid.out.AFR <- subset(aid.sum.AFR.mrg, cont=="AFR")
  popN <- read.table("pop.yr.csv", sep=",", header=T)
  aid.pc.AFR <- merge(aid.out.AFR, popN, by="cntry.code") 
  aid.pc.AFR$aid.pc <- aid.pc.AFR$aid.sum / aid.pc.AFR$X2015
  aidpc.AFR <- aid.pc.AFR[,c(1:5,62,63)]

  # total health expenditure (US$ per capita)
  hcexp.USDpc <- read.table("hc.exp.USDpc.csv", sep=",", header=T)
  hcexp.USDpc.reg <- merge(hcexp.USDpc, contreg, by="cntry.code")
  hcexp.USDpc.sort <- hcexp.USDpc.reg[order(hcexp.USDpc.reg[,2], decreasing=T),]
  hcexp.USDpc.AFR <- subset(hcexp.USDpc.sort, cont=="AFR")
  hinvest.tot1 <- merge(aidpc.AFR, hcexp.USDpc.AFR, by="cntry.code")
  hinvest.tot1$totpc.inv <- scale(log10(apply(hinvest.tot1[,c(7,8)], 1, sum, na.rm=T)), center=T, scale=T)
  hinvest.tot1$totpc.inv.rnk <- rank(hinvest.tot1$totpc.inv, na.last="keep", ties.method="average")
  hinvest <- hinvest.tot1[,c(1,13)]
  
# air quality
## PM2.5 from Sam Heft-Neal (Stanford)
PM25 <- as.data.frame(readRDS("country_level_pm25.rds", refhook = NULL))
colnames(PM25)[1] <- "cntry.code"
PM25AFR <- merge(PM25, contreg, by="cntry.code")
PM25AFR$lPM25 <- scale((PM25AFR$pm25_mean_pop/100), center=T, scale=T)
PM25.AFR <- PM25AFR[,c(1,9)]

# % population with access to improved water sources
impwater.pcaccess <- read.table("imp.water.pcaccess.csv", sep=",", header=T)
impwater.reg <- merge(impwater.pcaccess, contreg, by="cntry.code")  
impwater.sort <- impwater.reg[order(impwater.reg[,2], decreasing=T),]
impwater.AFR <- subset(impwater.sort, cont=="AFR")

# % population with access to improved sanitation facilities
impsan.pcaccess <- read.table("imp.san.pcaccess.csv", sep=",", header=T)
impsan.reg <- merge(impsan.pcaccess, contreg, by="cntry.code")  
impsan.sort <- impsan.reg[order(impsan.reg[,2], decreasing=T),]
impsan.AFR <- subset(impsan.sort, cont=="AFR")

waterhealthmrg <- merge(impwater.AFR, impsan.AFR, by="cntry.code",all.x=T, all.y=T)
waterhealth <- waterhealthmrg[,c(1,2,7)]
waterhealth.mrg <- waterhealth[!(waterhealth$cntry.code %in% smisnat),]
waterhealth.mrg$impwat.rnk <- rank(-waterhealth.mrg$imp.wat.pcaccess, na.last="keep", ties.method="average")  
waterhealth.mrg$impsan.rnk <- rank(-waterhealth.mrg$imp.san.pcaccess2011.15, na.last="keep", ties.method="average")  
waterhealth.mrg$waterhealth.grank <- 10^(apply(log10(waterhealth.mrg[,4:5]), 1, mean, na.rm=F))
waterhealth.mrg$lH20 <- scale((apply(waterhealth.mrg[,c(2,3)],1,"mean",na.rm=T)/100), center=T, scale=T)
waterhealthrnk <- waterhealth.mrg[,c(1,7)]

## food supply
foodsupply <- read.table("foodsupply.csv", sep=",", header=T)
foodsupply.cntry <- merge(foodsupply, FAO.cntry, by="countrycode")
foodsupply.reg <- merge(foodsupply.cntry, contreg, by="cntry.code")  
foodsupply.AFR <- subset(foodsupply.reg, cont=="AFR")
foodsupply13.AFR <- subset(foodsupply.AFR, year==2013)

  # kcalpcpd
  kcalpcpd.AFR <- subset(foodsupply13.AFR, element=="foodsupply")
  kcalpcpdAFR <- kcalpcpd.AFR[,c(1,6)]
  
  # protpcpd
  protpcpd.AFR <- subset(foodsupply13.AFR, element=="proteinsupply")
  colnames(protpcpd.AFR)[6] <- c("protpcpd")
  protpcpdAFR <- protpcpd.AFR[,c(1,6)]
  
foodsupplymrg <- merge(kcalpcpdAFR, protpcpdAFR, by="cntry.code")
foodsupply.mrg <- foodsupplymrg[!(foodsupplymrg$cntry.code %in% smisnat),]
foodsupply.mrg$kcal.rnk <- rank(-foodsupply.mrg$kcalpcpd, na.last="keep", ties.method="average")  
foodsupply.mrg$prot.rnk <- rank(-foodsupply.mrg$protpcpd, na.last="keep", ties.method="average")  
foodsupply.mrg$food.grank <- 10^(apply(log10(foodsupply.mrg[,4:5]), 1, mean, na.rm=F))
foodsupply.mrg$lkcal <- scale((foodsupply.mrg$kcalpcpd), center=T, scale=T)
foodsupplyrnk <- foodsupply.mrg[,c(1,7)]


# % infants exclusively breastfed for the first six months of life (%)
breastfed <- read.table("exclbreastfedAFR.csv", sep=",", header=T)
pcexclbf <- merge(breastfed, WHO.cntry, by="country")
colnames(pcexclbf)[2] <- "cntry.code"
pcexclbf$cntry.code <- factor(pcexclbf$cntry.code)
# most up-to-date year
cntr.vec <- as.vector(unique(names(table(pcexclbf$cntry.code))))
lcntr <- length(cntr.vec)
pcexclbf.max <- pcexclbf[1,]
for (i in 1:lcntr) {
  cntry.dat <- subset(pcexclbf, cntry.code == cntr.vec[i])
  cntry.max <- subset(cntry.dat, year == max(cntry.dat$year))
  pcexclbf.max <- rbind(pcexclbf.max,cntry.max)
}
pcexclbf.max <- pcexclbf.max[-1, ]
pcexclbf.reg <- merge(pcexclbf.max, contreg, by="cntry.code")
pcexclbf.sort <- pcexclbf.reg[order(pcexclbf.reg[,5], decreasing=F),]
pcexclbf.sort$breastfed.sc <- scale(logit(pcexclbf.sort$breastfed/100), center=T, scale=T)
breastfedrnk <- pcexclbf.sort[,c(1,11)]

# trade export value (USD)
export <- read.table("exportvalue.csv", sep=",", header=T)
export.cntry <- merge(export, FAO.cntry, by="countrycode")
export.reg <- merge(export.cntry, contreg, by="cntry.code")
export.AFR <- subset(export.reg, cont=="AFR")
export.pc.AFR <- merge(export.AFR, popN, by="cntry.code")  
export.pc.AFR$exportvaluepc <- export.pc.AFR$exportvalue / export.pc.AFR$X2013
export.pc.AFR$lTRAD <- log10(export.pc.AFR$exportvalue)
exportpc.AFR <- export.pc.AFR[,c(1,3,65,66,67)]
exportrnk <- exportpc.AFR[,c(1,5)]



## STRUCTURAL EQUATION MODELS (DALY5 or LIFE EXPECTANCY)

##################################################################
# CHOOSE WHICH RESPONSE
DALY5LE <- data.frame(DALY5.AFR$cntry.code, scale(sqrt(DALY5.AFR$DALY5), scale=T, center=T))
colnames(DALY5LE) <- c("cntry.code", "RESP")
##################################################################

# combine all data into one data.frame
dat.sem.1 <- merge(DALY5LE, env.dat, by="cntry.code", all.x=T, all.y=T)
dat.sem.2 <- merge(dat.sem.1, popD.AFR, by="cntry.code", all.x=T, all.y=T)
dat.sem.3 <- merge(dat.sem.2, GDP.AFR, by="cntry.code", all.x=T, all.y=T)
dat.sem.4 <- merge(dat.sem.3, iiag.sort[,c(5,6)], by="cntry.code", all.x=T, all.y=T)
dat.sem.5 <- merge(dat.sem.4, povgap.AFR, by="cntry.code", all.x=T, all.y=T)
dat.sem.6 <- merge(dat.sem.5, hinvest, by="cntry.code", all.x=T, all.y=T)
dat.sem.7 <- merge(dat.sem.6, waterhealthrnk, by="cntry.code", all.x=T, all.y=T)
dat.sem.8 <- merge(dat.sem.7, foodsupplyrnk, by="cntry.code", all.x=T, all.y=T)
#dat.sem.9 <- merge(dat.sem.8, exportrnk, by="cntry.code", all.x=T, all.y=T)
dat.sem.9 <- merge(dat.sem.8, breastfedrnk, by="cntry.code", all.x=T, all.y=T)
dat.sem.10 <- merge(dat.sem.9, PM25.AFR, by="cntry.code", all.x=T, all.y=T)
colnames(dat.sem.10) <- c("cntry.code", "RESP","ENVr","POPDr","lPOPD","GDPr","lGDP","lGOV","POVr","lHINV","lH2O","lFOOD","lBF","lPM25")
rm(dat.sem.1,dat.sem.2,dat.sem.3,dat.sem.4,dat.sem.5,dat.sem.6,dat.sem.7,dat.sem.8,dat.sem.9)
dat.sem <- dat.sem.10[!(dat.sem.10$cntry.code %in% smisnat),]

# remove NAs
datsem <- na.omit(dat.sem[,-c(4,6,9)])

# covariance matrix for path analysis
cov.mat <- cov(model.matrix(~ RESP + ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lBF + lPM25, data=datsem))[-1,-1]
cov.mat

# correlation matrix
cor.mat <- cor(na.omit(datsem[,-c(1,2)]), method="kendall")
cor.mat
cor.mat <- cor(na.omit(datsem[,-c(1,2)]), method="pearson")
cor.mat*cor.mat
pairs(datsem[,-c(1)], pch=19, cex=0.7)


# models
# saturated model
mod1 <- specifyModel(text="
  ENVr -> RESP, beta1, NA
  lPOPD -> RESP, beta2, NA
  lGDP -> RESP, beta3, NA
  lGOV -> RESP, beta4, NA
  lHINV -> RESP, beta5, NA
  lH2O -> RESP, beta6, NA
  lFOOD -> RESP, beta7, NA
  lBF -> RESP, beta8, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 2 ENV
mod2 <- specifyModel(text="
  ENVr -> RESP, beta1, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")


# model 3 POPD
mod3 <- specifyModel(text="
  lPOPD -> RESP, beta1, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 4 GDP
mod4 <- specifyModel(text="
  lGDP -> RESP, beta1, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 5 GDP + HINV
mod5 <- specifyModel(text="
  lGDP -> RESP, beta1, NA
  lHINV <- RESP, beta2, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 6 GDP + GOV
mod6 <- specifyModel(text="
  lGDP -> RESP, beta1, NA
  lGOV <- RESP, beta2, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 7 GOV
mod7 <- specifyModel(text="
  lGOV -> RESP, beta1, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 8 HINV
mod8 <- specifyModel(text="
  lHINV -> RESP, beta1, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 9 H2O
mod9 <- specifyModel(text="
  lH2O -> RESP, beta1, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 10 FOOD
mod10 <- specifyModel(text="
  lFOOD -> RESP, beta1, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 11 BREASTFEEDING
mod11 <- specifyModel(text="
  lBF -> RESP, beta1, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 12 FOOD & BREASTFEEDING
mod12 <- specifyModel(text="
  lFOOD -> RESP, beta1, NA
  lBF -> RESP, beta2, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 13 ENV & POPD
mod13 <- specifyModel(text="
  ENVr -> RESP, beta1, NA
  lPOPD -> RESP, beta2, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 14 ENV + H2O
mod14 <- specifyModel(text="
  ENVr -> RESP, beta1, NA
  lH2O -> RESP, beta2, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")

# model 15 ENV + POPD + H2O
mod15 <- specifyModel(text="
  ENVr -> RESP, beta1, NA
  lPOPD -> RESP, beta2, NA
  lH2O -> RESP, beta3, NA
  lPOPD -> lGOV, pi1, NA
  lGDP -> lHINV, pi2, NA 
  lPOPD -> ENVr, pi3, NA
  ENVr -> lFOOD, pi4, NA
  ENVr -> lH2O, pi5, NA
  lFOOD -> lBF, pi6, NA
  ENVr <-> ENVr, gam1, NA
  lPOPD <-> lPOPD, gam2, NA
  lGDP <-> lGDP, gam3, NA
  lGOV <-> lGOV, gam4, NA
  lHINV <-> lHINV, gam5, NA
  lH2O <-> lH2O, gam6, NA
  lFOOD <-> lFOOD, gam7, NA
  lBF <-> lBF, gam8
  RESP <-> RESP, gam9, NA
")


  model1 <- "ENV+POPD+GDP+GOV+HINV+H2O+FOOD+BF"
  model2 <- "ENV"
  model3 <- "POPD"
  model4 <- "GDP"
  model5 <- "GDP+HINV"
  model6 <- "GDP+GOV"
  model7 <- "GOV"
  model8 <- "HINV"
  model9 <- "H2O"
  model10 <- "FOOD"
  model11 <- "BF"
  model12 <- "FOOD+BF"
  model13 <- "ENV+POPD"
  model14 <- "ENV+H2O"
  model15 <- "ENV+POPD+H2O"
  
mod.lab.vec <- c(model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13,model14,model15)

  ## fit path models
  sem.mod1 <- sem(mod1, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod2 <- sem(mod2, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod3 <- sem(mod3, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod4 <- sem(mod4, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod5 <- sem(mod5, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod6 <- sem(mod6, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod7 <- sem(mod7, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod8 <- sem(mod8, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod9 <- sem(mod9, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod10 <- sem(mod10, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod11 <- sem(mod11, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod12 <- sem(mod12, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod13 <- sem(mod13, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod14 <- sem(mod14, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod15 <- sem(mod15, cov.mat, dim(datsem)[1], standardized=T)
  
  # plots
  semPaths(sem.mod1, what="col", whatLabels="stand", layout="spring", label.cex=1.3, edge.label.cex=1.3, nCharNodes=7)

  # summaries
  sum1 <- summary(sem.mod1)
  sum2 <- summary(sem.mod2)
  sum3 <- summary(sem.mod3)
  sum4 <- summary(sem.mod4)
  sum5 <- summary(sem.mod5)
  sum6 <- summary(sem.mod6)
  sum7 <- summary(sem.mod7)
  sum8 <- summary(sem.mod8)
  sum9 <- summary(sem.mod9)
  sum10 <- summary(sem.mod10)
  sum11 <- summary(sem.mod11)
  sum12 <- summary(sem.mod12)
  sum13 <- summary(sem.mod13)
  sum14 <- summary(sem.mod14)
  sum15 <- summary(sem.mod15)
  
  # standardised coefficients
  stdCoef(sem.mod1)
  stdCoef(sem.mod2)
  stdCoef(sem.mod3)
  stdCoef(sem.mod4)
  stdCoef(sem.mod5)
  stdCoef(sem.mod6)
  stdCoef(sem.mod7)
  stdCoef(sem.mod8)
  stdCoef(sem.mod9)
  stdCoef(sem.mod10)
  stdCoef(sem.mod11)
  stdCoef(sem.mod12)
  stdCoef(sem.mod13)
  stdCoef(sem.mod14)
  stdCoef(sem.mod15)
  
  # SEM goodness-of-fit
  # a cutoff value close to .95 for TLI, BL89, CFI, RNI, and Gamma Hat;
  # a cutoff value close to .90 for Mc; a cutoff value close to .08 for SRMR;
  # and a cutoff value close to .06 for RMSEA
  # incremental fit index, also known as Bollen's IFI, is also relatively insensitive to sample size.
  # Values that exceed .90 are regarded as acceptable, although this index can exceed 1.
  # To compute the IFI, first the difference between the chi square of the independence model--
  # in which variables are uncorrelated--and the chi-square of the target model is calculated.
  # Next, the difference between the chi-square of the target model and the df for the target model is calculated.
  # The ratio of these values represents the IFI.
  # McDonald's Non-Centrality Index
  Mc.vec <- c(summaryGOF(sem.mod1, digits=4)$Mc,summaryGOF(sem.mod2, digits=4)$Mc,summaryGOF(sem.mod3, digits=4)$Mc,summaryGOF(sem.mod4, digits=4)$Mc,summaryGOF(sem.mod5, digits=4)$Mc,summaryGOF(sem.mod6, digits=4)$Mc,summaryGOF(sem.mod7, digits=4)$Mc,summaryGOF(sem.mod8, digits=4)$Mc,summaryGOF(sem.mod9, digits=4)$Mc,summaryGOF(sem.mod10, digits=4)$Mc,summaryGOF(sem.mod11, digits=4)$Mc,summaryGOF(sem.mod12, digits=4)$Mc,summaryGOF(sem.mod13, digits=4)$Mc,summaryGOF(sem.mod14, digits=4)$Mc,summaryGOF(sem.mod15, digits=4)$Mc)
  # Bollen's Incremental Fit Index
  IFI.vec <- c(summaryGOF(sem.mod1, digits=4)$IFI,summaryGOF(sem.mod2, digits=4)$IFI,summaryGOF(sem.mod3, digits=4)$IFI,summaryGOF(sem.mod4, digits=4)$IFI,summaryGOF(sem.mod5, digits=4)$IFI,summaryGOF(sem.mod6, digits=4)$IFI,summaryGOF(sem.mod7, digits=4)$IFI,summaryGOF(sem.mod8, digits=4)$IFI,summaryGOF(sem.mod9, digits=4)$IFI,summaryGOF(sem.mod10, digits=4)$IFI,summaryGOF(sem.mod11, digits=4)$IFI,summaryGOF(sem.mod12, digits=4)$IFI,summaryGOF(sem.mod13, digits=4)$IFI,summaryGOF(sem.mod14, digits=4)$IFI,summaryGOF(sem.mod15, digits=4)$IFI)
  # Expected cross-validation index (ECVI) (Browne & Cudeck 1992)
  ECVI.vec <- c(summaryGOF(sem.mod1, digits=4)$ECVI,summaryGOF(sem.mod2, digits=4)$ECVI,summaryGOF(sem.mod3, digits=4)$ECVI,summaryGOF(sem.mod4, digits=4)$ECVI,summaryGOF(sem.mod5, digits=4)$ECVI,summaryGOF(sem.mod6, digits=4)$ECVI,summaryGOF(sem.mod7, digits=4)$ECVI,summaryGOF(sem.mod8, digits=4)$ECVI,summaryGOF(sem.mod9, digits=4)$ECVI,summaryGOF(sem.mod10, digits=4)$ECVI,summaryGOF(sem.mod11, digits=4)$ECVI,summaryGOF(sem.mod12, digits=4)$ECVI,summaryGOF(sem.mod13, digits=4)$ECVI,summaryGOF(sem.mod14, digits=4)$ECVI,summaryGOF(sem.mod15, digits=4)$ECVI)
  # chi-square/df ratio (2 to 5 acceptable)
  CSDFR.vec <- c(summaryGOF(sem.mod1, digits=4)$chisq.df,summaryGOF(sem.mod2, digits=4)$chisq.df,summaryGOF(sem.mod3, digits=4)$chisq.df,summaryGOF(sem.mod4, digits=4)$chisq.df,summaryGOF(sem.mod5, digits=4)$chisq.df,summaryGOF(sem.mod6, digits=4)$chisq.df,summaryGOF(sem.mod7, digits=4)$chisq.df,summaryGOF(sem.mod8, digits=4)$chisq.df,summaryGOF(sem.mod9, digits=4)$chisq.df,summaryGOF(sem.mod10, digits=4)$chisq.df,summaryGOF(sem.mod11, digits=4)$chisq.df,summaryGOF(sem.mod12, digits=4)$chisq.df,summaryGOF(sem.mod13, digits=4)$chisq.df,summaryGOF(sem.mod14, digits=4)$chisq.df,summaryGOF(sem.mod15, digits=4)$chisq.df)

  # summary table
  mod.lab <- 1:15
  df.vec <- c(sum1$df, sum2$df, sum3$df, sum4$df, sum5$df, sum6$df, sum7$df, sum8$df, sum9$df, sum10$df, sum11$df, sum12$df, sum13$df, sum14$df, sum15$df)
  chisq.vec <- c(sum1$chisq, sum2$chisq, sum3$chisq, sum4$chisq, sum5$chisq, sum6$chisq, sum7$chisq, sum8$chisq, sum9$chisq, sum10$chisq, sum11$chisq, sum12$chisq, sum13$chisq, sum14$chisq, sum15$chisq)

  # BIC ranks
  delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
  weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
  
  BIC.vec <- c(sum1$BIC, sum2$BIC, sum3$BIC, sum4$BIC, sum5$BIC, sum6$BIC, sum7$BIC, sum8$BIC, sum9$BIC, sum10$BIC, sum11$BIC, sum12$BIC, sum13$BIC, sum14$BIC, sum15$BIC)
  dBIC.vec <- delta.IC(BIC.vec)
  wBIC.vec <- weight.IC(dBIC.vec)

  AIC.vec <- c(sum1$AIC, sum2$AIC, sum3$AIC, sum4$AIC, sum5$AIC, sum6$AIC, sum7$AIC, sum8$AIC, sum9$AIC, sum10$AIC, sum11$AIC, sum12$AIC, sum13$AIC, sum14$AIC, sum15$AIC)
  dAIC.vec <- delta.IC(AIC.vec)
  wAIC.vec <- weight.IC(dAIC.vec)

  # BIC results dataframe
  table<-cbind(mod.lab,df.vec,chisq.vec,BIC.vec,dBIC.vec,wBIC.vec,Mc.vec,IFI.vec,CSDFR.vec)
  colnames(table)<-c("model","df","chisq","BIC","dBIC","wBIC","Mc","IFI","CSDFR")
  rownames(table)<- mod.lab.vec
  
  # table sorted by wBIC (Mc & IFI > 0.95 considered 'good' fit); chi-square/df ratio 2-5 ok
  summary.table<-table[order(table[,6],decreasing=TRUE),1:9]
  summary.table

  # AIC results dataframe
  table<-cbind(mod.lab,df.vec,chisq.vec,AIC.vec,dAIC.vec,wAIC.vec,Mc.vec,IFI.vec,CSDFR.vec)
  colnames(table)<-c("model","df","chisq","AIC","dAIC","wAIC","Mc","IFI","CSDFR")
  rownames(table)<- mod.lab.vec
  
  # table sorted by wAIC (Mc & IFI > 0.95 considered 'good' fit)
  summary.table<-table[order(table[,6],decreasing=TRUE),1:9]
  summary.table
  
  ########################################################
  ## boosted regression trees
  ########################################################
  
  ## BRT
  brt.fit <- gbm.step(datsem, gbm.x = attr(datsem, "names")[c(3:9)], gbm.y = attr(datsem, "names")[2], family="gaussian", tolerance = 0.001, learning.rate = 0.00001, bag.fraction=0.8, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
  summary(brt.fit)
  gbm.plot(brt.fit)
  tmp <- gbm.plot.fits(brt.fit, v=0)
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
  D2 # % deviance explained
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  print(c(CV.cor, CV.cor.se))
  
 
  
  ################################# 
  ## GLMMs
  #################################
  ## functions
  delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
  weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
  
  # Set functions
  AICc <- function(...) {
    models <- list(...)
    num.mod <- length(models)
    AICcs <- numeric(num.mod)
    ns <- numeric(num.mod)
    ks <- numeric(num.mod)
    AICc.vec <- rep(0,num.mod)
    for (i in 1:num.mod) {
      if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
      if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
      AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
      ns[i] <- n
      ks[i] <- k
      AICc.vec[i] <- AICcs[i]
    }
    return(AICc.vec)
  }
  
  linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
    fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
    AIC.vec <- c(AICc(fit.full),AICc(fit.null))
    dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
    ER <- wAIC.vec[1]/wAIC.vec[2]
    r.sq.adj <- as.numeric(summary(fit.full)[9])
    return(ER)
  }
  
  ## add regional categories for mixed effects
  dat.GLMM <- merge(datsem, regions.AFR, by="cntry.code")
  str(dat.GLMM)
  # UN, AU, or WHO
  dat.GLMM$region <- factor(dat.GLMM$UN.region)
  table(dat.GLMM$region)
  
  ## define model set
  mod1 <- "RESP ~ ENVr+ lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + (1|region)"
  mod2 <- "RESP ~ ENVr + (1|region)"
  mod3 <- "RESP ~ lPOPD + (1|region)"
  mod4 <- "RESP ~ lGDP + (1|region)"
  mod5 <- "RESP ~ lGDP + lHINV + (1|region)"
  mod6 <- "RESP ~ lGDP + lGOV + (1|region)"
  mod7 <- "RESP ~ lGOV + (1|region)"
  mod8 <- "RESP ~ lHINV + (1|region)"
  mod9 <- "RESP ~ lH2O + (1|region)"
  mod10 <- "RESP ~ lFOOD + (1|region)"
  mod11 <- "RESP ~ ENVr + lPOPD + (1|region)"
  mod12 <- "RESP ~ ENVr + lH2O + (1|region)"
  mod13 <- "RESP ~ ENVr + lPOPD + lH2O + (1|region)"
  mod14 <- "RESP ~ 1 + (1|region)"
  
  ## Make model vector
  mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14)
  
  ## Define n.mod
  n.mod <- length(mod.vec)
  
  # Model fitting and logLik output loop
  Modnum <- length(mod.vec)
  LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
  mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
  mod.num <- seq(1,Modnum,1)

  for(i in 1:Modnum) {
    fit <- lmer(as.formula(mod.vec[i]),data=dat.GLMM, na.action=na.omit)
    assign(paste("fit",i,sep=""), fit)
    mod.list[[i]] <- fit
    LL.vec[i] <- as.numeric(logLik(fit))
    k.vec[i] <- attr(logLik(fit),"df")
    AICc.vec[i] <- r.squared(fit)$AIC
    
    Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
    Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
    
    print(i)
  }
  dAICc <- delta.IC(AICc.vec)
  wAICc <- weight.IC(dAICc)
  
  sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),round(Rm,4),round(Rc,4))
  colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
  row.names(sumtable) <- as.character(mod.vec)
  summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
  summary.table

  ## saturated residual diagnostic
  mod.saturated.glm <- "RESP ~ ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD"
  fit <- glm(as.formula(mod.saturated.glm),family=gaussian(link="identity"), data=dat.GLMM, na.action=na.omit)
  plot(fit)
  
  top.model <- "RESP ~ ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD"
  topfit <- glm(as.formula(top.model),family=gaussian(link="identity"), data=dat.GLMM, na.action=na.omit)
  par(mfrow=c(3,3))
  plot1 <- termplot(topfit,terms="lGDP",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="GDP",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  plot1 <- termplot(topfit,terms="lH2O",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="H2O",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  plot1 <- termplot(topfit,terms="lPOPD",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="POPD",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  plot1 <- termplot(topfit,terms="lGOV",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="GOV",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  plot1 <- termplot(topfit,terms="ENVr",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="ENV",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  plot1 <- termplot(topfit,terms="lFOOD",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="FOOD",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  plot1 <- termplot(topfit,terms="lHINV",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="HINV",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  par(mfrow=c(1,1))
  
  
  
  
####################################################################
  ## STRUCTURAL EQUATION MODELS (COMBINED CHILD HEALTH INDICATOR)
  # combine all data into one data.frame
  dat.sem.1 <- merge(chealth.out, env.dat, by="cntry.code", all.x=T, all.y=T)
  dat.sem.2 <- merge(dat.sem.1, popD.AFR, by="cntry.code", all.x=T, all.y=T)
  dat.sem.3 <- merge(dat.sem.2, GDP.AFR, by="cntry.code", all.x=T, all.y=T)
  dat.sem.4 <- merge(dat.sem.3, iiag.sort[,c(5,6)], by="cntry.code", all.x=T, all.y=T)
  dat.sem.5 <- merge(dat.sem.4, povgap.AFR, by="cntry.code", all.x=T, all.y=T)
  dat.sem.6 <- merge(dat.sem.5, hinvest, by="cntry.code", all.x=T, all.y=T)
  dat.sem.7 <- merge(dat.sem.6, waterhealthrnk, by="cntry.code", all.x=T, all.y=T)
  dat.sem.8 <- merge(dat.sem.7, foodsupplyrnk, by="cntry.code", all.x=T, all.y=T)
  dat.sem.9 <- merge(dat.sem.8, breastfedrnk, by="cntry.code", all.x=T, all.y=T)
  dat.sem.10 <- merge(dat.sem.9, PM25.AFR, by="cntry.code", all.x=T, all.y=T)
  colnames(dat.sem.10) <- c("cntry.code", "stCHLTH","ENVr","POPDr","lPOPD","GDPr","lGDP","lGOV","POVr","lHINV","lH2O","lFOOD","lBF","lPM25")
  rm(dat.sem.1,dat.sem.2,dat.sem.3,dat.sem.4,dat.sem.5,dat.sem.6,dat.sem.7,dat.sem.8,dat.sem.9)
  dat.sem <- dat.sem.10[!(dat.sem.10$cntry.code %in% smisnat),]
  
  # remove NAs
  datsem <- na.omit(dat.sem[,-c(4,6,9)])

  # covariance matrix for path analysis
  cov.mat <- cov(model.matrix(~ stCHLTH + ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lBF + lPM25, data=datsem))[-1,-1]
  cov.mat
  
  # correlation matrix
  cor.mat <- cor(na.omit(datsem[,-c(1:2)]), method="kendall")
  cor.mat
  cor.mat <- cor(na.omit(datsem[,-c(1:2)]), method="pearson")
  cor.mat*cor.mat
  pairs(datsem[,-c(1:2)], pch=19, cex=0.7)
  
# models
# saturated model
  mod1 <- specifyModel(text="
                       ENVr -> stCHLTH, beta1, NA
                       lPOPD -> stCHLTH, beta2, NA
                       lGDP -> stCHLTH, beta3, NA
                       lGOV -> stCHLTH, beta4, NA
                       lHINV -> stCHLTH, beta5, NA
                       lH2O -> stCHLTH, beta6, NA
                       lFOOD -> stCHLTH, beta7, NA
                       lBF -> stCHLTH, beta8, NA
                       lPM25 -> stCHLTH, beta9, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")

# model 2 ENV + PM25
  mod2 <- specifyModel(text="
                       ENVr -> stCHLTH, beta1, NA
                       lPM25 -> stCHLTH, beta2, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 3 ENV
  mod3 <- specifyModel(text="
                       ENVr -> stCHLTH, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 4 POPD
  mod4 <- specifyModel(text="
                       lPOPD -> stCHLTH, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 5 GDP
  mod5 <- specifyModel(text="
                       lGDP -> stCHLTH, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 6 GDP + HINV
  mod6 <- specifyModel(text="
                       lGDP -> stCHLTH, beta1, NA
                       lHINV <- stCHLTH, beta2, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 7 GDP + GOV
  mod7 <- specifyModel(text="
                       lGDP -> stCHLTH, beta1, NA
                       lGOV <- stCHLTH, beta2, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 8 GOV
  mod8 <- specifyModel(text="
                       lGOV -> stCHLTH, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 9 HINV
  mod9 <- specifyModel(text="
                       lHINV -> stCHLTH, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 10 H2O
  mod10 <- specifyModel(text="
                       lH2O -> stCHLTH, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       stCHLTH <-> stCHLTH, gam10, NA
                       ")
  
  # model 11 FOOD
  mod11 <- specifyModel(text="
                        lFOOD -> stCHLTH, beta1, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        stCHLTH <-> stCHLTH, gam10, NA
                        ")
  
  # model 12 BREASTFEEDING
  mod12 <- specifyModel(text="
                        lBF -> stCHLTH, beta1, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        stCHLTH <-> stCHLTH, gam10, NA
                        ")
  
  # model 13 FOOD & BREASTFEEDING
  mod13 <- specifyModel(text="
                        lFOOD -> stCHLTH, beta1, NA
                        lBF -> stCHLTH, beta2, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        stCHLTH <-> stCHLTH, gam10, NA
                        ")
  
  # model 14 ENV & POPD
  mod14 <- specifyModel(text="
                        ENVr -> stCHLTH, beta1, NA
                        lPOPD -> stCHLTH, beta2, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        stCHLTH <-> stCHLTH, gam10, NA
                        ")
  
  # model 15 ENV + H2O
  mod15 <- specifyModel(text="
                        ENVr -> stCHLTH, beta1, NA
                        lH2O -> stCHLTH, beta2, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        stCHLTH <-> stCHLTH, gam10, NA
                        ")
  
  # model 16 ENV + POPD + H2O
  mod16 <- specifyModel(text="
                        ENVr -> stCHLTH, beta1, NA
                        lPOPD -> stCHLTH, beta2, NA
                        lH2O -> stCHLTH, beta3, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        stCHLTH <-> stCHLTH, gam10, NA
                        ")
  
  # model 17 PM25
  mod17 <- specifyModel(text="
                        lPM25 -> stCHLTH, beta1, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        stCHLTH <-> stCHLTH, gam10, NA
                        ")
  
  model1 <- "ENV+POPD+GDP+GOV+HINV+H2O+FOOD+BF+PM25"
  model2 <- "ENV+PM25"
  model3 <- "ENV"
  model4 <- "POPD"
  model5 <- "GDP"
  model6 <- "GDP+HINV"
  model7 <- "GDP+GOV"
  model8 <- "GOV"
  model9 <- "HINV"
  model10 <- "H2O"
  model11 <- "FOOD"
  model12 <- "BF"
  model13 <- "FOOD+BF"
  model14 <- "ENV+POPD"
  model15 <- "ENV+H2O"
  model16 <- "ENV+POPD+H2O"
  model17 <- "PM25"
  
  mod.lab.vec <- c(model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13,model14,model15,model16,model17)
  
  ## fit path models
  sem.mod1 <- sem(mod1, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod2 <- sem(mod2, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod3 <- sem(mod3, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod4 <- sem(mod4, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod5 <- sem(mod5, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod6 <- sem(mod6, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod7 <- sem(mod7, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod8 <- sem(mod8, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod9 <- sem(mod9, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod10 <- sem(mod10, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod11 <- sem(mod11, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod12 <- sem(mod12, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod13 <- sem(mod13, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod14 <- sem(mod14, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod15 <- sem(mod15, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod16 <- sem(mod16, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod17 <- sem(mod17, cov.mat, dim(datsem)[1], standardized=T)
  
  # plots
  semPaths(sem.mod1, what="col", whatLabels="stand", layout="spring", label.cex=1.3, edge.label.cex=1.3, nCharNodes=7)
  
  # summaries
  sum1 <- summary(sem.mod1)
  sum2 <- summary(sem.mod2)
  sum3 <- summary(sem.mod3)
  sum4 <- summary(sem.mod4)
  sum5 <- summary(sem.mod5)
  sum6 <- summary(sem.mod6)
  sum7 <- summary(sem.mod7)
  sum8 <- summary(sem.mod8)
  sum9 <- summary(sem.mod9)
  sum10 <- summary(sem.mod10)
  sum11 <- summary(sem.mod11)
  sum12 <- summary(sem.mod12)
  sum13 <- summary(sem.mod13)
  sum14 <- summary(sem.mod14)
  sum15 <- summary(sem.mod15)
  sum16 <- summary(sem.mod16)
  sum17 <- summary(sem.mod17)
  
  # standardised coefficients
  stdCoef(sem.mod1)
  stdCoef(sem.mod2)
  stdCoef(sem.mod3)
  stdCoef(sem.mod4)
  stdCoef(sem.mod5)
  stdCoef(sem.mod6)
  stdCoef(sem.mod7)
  stdCoef(sem.mod8)
  stdCoef(sem.mod9)
  stdCoef(sem.mod10)
  stdCoef(sem.mod11)
  stdCoef(sem.mod12)
  stdCoef(sem.mod13)
  stdCoef(sem.mod14)
  stdCoef(sem.mod15)
  stdCoef(sem.mod16)
  stdCoef(sem.mod17)
  
  # SEM goodness-of-fit
  # a cutoff value close to .95 for TLI, BL89, CFI, RNI, and Gamma Hat;
  # a cutoff value close to .90 for Mc; a cutoff value close to .08 for SRMR;
  # and a cutoff value close to .06 for RMSEA
  # incremental fit index, also known as Bollen's IFI, is also relatively insensitive to sample size.
  # Values that exceed .90 are regarded as acceptable, although this index can exceed 1.
  # To compute the IFI, first the difference between the chi square of the independence model--
  # in which variables are uncorrelated--and the chi-square of the target model is calculated.
  # Next, the difference between the chi-square of the target model and the df for the target model is calculated.
  # The ratio of these values represents the IFI.
  # McDonald's Non-Centrality Index
  Mc.vec <- c(summaryGOF(sem.mod1, digits=4)$Mc,summaryGOF(sem.mod2, digits=4)$Mc,summaryGOF(sem.mod3, digits=4)$Mc,summaryGOF(sem.mod4, digits=4)$Mc,summaryGOF(sem.mod5, digits=4)$Mc,summaryGOF(sem.mod6, digits=4)$Mc,summaryGOF(sem.mod7, digits=4)$Mc,summaryGOF(sem.mod8, digits=4)$Mc,summaryGOF(sem.mod9, digits=4)$Mc,summaryGOF(sem.mod10, digits=4)$Mc,summaryGOF(sem.mod11, digits=4)$Mc,summaryGOF(sem.mod12, digits=4)$Mc,summaryGOF(sem.mod13, digits=4)$Mc,summaryGOF(sem.mod14, digits=4)$Mc,summaryGOF(sem.mod15, digits=4)$Mc,summaryGOF(sem.mod16, digits=4)$Mc,summaryGOF(sem.mod17, digits=4)$Mc)
  # Bollen's Incremental Fit Index
  IFI.vec <- c(summaryGOF(sem.mod1, digits=4)$IFI,summaryGOF(sem.mod2, digits=4)$IFI,summaryGOF(sem.mod3, digits=4)$IFI,summaryGOF(sem.mod4, digits=4)$IFI,summaryGOF(sem.mod5, digits=4)$IFI,summaryGOF(sem.mod6, digits=4)$IFI,summaryGOF(sem.mod7, digits=4)$IFI,summaryGOF(sem.mod8, digits=4)$IFI,summaryGOF(sem.mod9, digits=4)$IFI,summaryGOF(sem.mod10, digits=4)$IFI,summaryGOF(sem.mod11, digits=4)$IFI,summaryGOF(sem.mod12, digits=4)$IFI,summaryGOF(sem.mod13, digits=4)$IFI,summaryGOF(sem.mod14, digits=4)$IFI,summaryGOF(sem.mod15, digits=4)$IFI,summaryGOF(sem.mod16, digits=4)$IFI,summaryGOF(sem.mod17, digits=4)$IFI)
  # Expected cross-validation index (ECVI) (Browne & Cudeck 1992)
  ECVI.vec <- c(summaryGOF(sem.mod1, digits=4)$ECVI,summaryGOF(sem.mod2, digits=4)$ECVI,summaryGOF(sem.mod3, digits=4)$ECVI,summaryGOF(sem.mod4, digits=4)$ECVI,summaryGOF(sem.mod5, digits=4)$ECVI,summaryGOF(sem.mod6, digits=4)$ECVI,summaryGOF(sem.mod7, digits=4)$ECVI,summaryGOF(sem.mod8, digits=4)$ECVI,summaryGOF(sem.mod9, digits=4)$ECVI,summaryGOF(sem.mod10, digits=4)$ECVI,summaryGOF(sem.mod11, digits=4)$ECVI,summaryGOF(sem.mod12, digits=4)$ECVI,summaryGOF(sem.mod13, digits=4)$ECVI,summaryGOF(sem.mod14, digits=4)$ECVI,summaryGOF(sem.mod15, digits=4)$ECVI,summaryGOF(sem.mod16, digits=4)$ECVI,summaryGOF(sem.mod17, digits=4)$ECVI)
  # chi-square/df ratio (2 to 5 acceptable)
  CSDFR.vec <- c(summaryGOF(sem.mod1, digits=4)$chisq.df,summaryGOF(sem.mod2, digits=4)$chisq.df,summaryGOF(sem.mod3, digits=4)$chisq.df,summaryGOF(sem.mod4, digits=4)$chisq.df,summaryGOF(sem.mod5, digits=4)$chisq.df,summaryGOF(sem.mod6, digits=4)$chisq.df,summaryGOF(sem.mod7, digits=4)$chisq.df,summaryGOF(sem.mod8, digits=4)$chisq.df,summaryGOF(sem.mod9, digits=4)$chisq.df,summaryGOF(sem.mod10, digits=4)$chisq.df,summaryGOF(sem.mod11, digits=4)$chisq.df,summaryGOF(sem.mod12, digits=4)$chisq.df,summaryGOF(sem.mod13, digits=4)$chisq.df,summaryGOF(sem.mod14, digits=4)$chisq.df,summaryGOF(sem.mod15, digits=4)$chisq.df,summaryGOF(sem.mod16, digits=4)$chisq.df,summaryGOF(sem.mod17, digits=4)$chisq.df)
  
  # summary table
  mod.lab <- 1:17
  df.vec <- c(sum1$df, sum2$df, sum3$df, sum4$df, sum5$df, sum6$df, sum7$df, sum8$df, sum9$df, sum10$df, sum11$df, sum12$df, sum13$df, sum14$df, sum15$df, sum16$df, sum17$df)
  chisq.vec <- c(sum1$chisq, sum2$chisq, sum3$chisq, sum4$chisq, sum5$chisq, sum6$chisq, sum7$chisq, sum8$chisq, sum9$chisq, sum10$chisq, sum11$chisq, sum12$chisq, sum13$chisq, sum14$chisq, sum15$chisq, sum16$chisq, sum17$chisq)
  
  # BIC ranks
  delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
  weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
  
  BIC.vec <- c(sum1$BIC, sum2$BIC, sum3$BIC, sum4$BIC, sum5$BIC, sum6$BIC, sum7$BIC, sum8$BIC, sum9$BIC, sum10$BIC, sum11$BIC, sum12$BIC, sum13$BIC, sum14$BIC, sum15$BIC, sum16$BIC, sum17$BIC)
  dBIC.vec <- delta.IC(BIC.vec)
  wBIC.vec <- weight.IC(dBIC.vec)

  AIC.vec <- c(sum1$AIC, sum2$AIC, sum3$AIC, sum4$AIC, sum5$AIC, sum6$AIC, sum7$AIC, sum8$AIC, sum9$AIC, sum10$AIC, sum11$AIC, sum12$AIC, sum13$AIC, sum14$AIC, sum15$AIC, sum16$AIC, sum17$AIC)
  dAIC.vec <- delta.IC(AIC.vec)
  wAIC.vec <- weight.IC(dAIC.vec)

  # BIC results dataframe
  table<-cbind(mod.lab,df.vec,chisq.vec,BIC.vec,dBIC.vec,wBIC.vec,Mc.vec,IFI.vec,CSDFR.vec)
  colnames(table)<-c("model","df","chisq","BIC","dBIC","wBIC","Mc","IFI","CSDFR")
  rownames(table)<- mod.lab.vec
  
  # table sorted by wBIC (Mc & IFI > 0.95 considered 'good' fit); chi-square/df ratio 2-5 ok
  summary.table<-table[order(table[,6],decreasing=TRUE),1:9]
  summary.table

  # AIC results dataframe
  table<-cbind(mod.lab,df.vec,chisq.vec,AIC.vec,dAIC.vec,wAIC.vec,Mc.vec,IFI.vec,CSDFR.vec)
  colnames(table)<-c("model","df","chisq","AIC","dAIC","wAIC","Mc","IFI","CSDFR")
  rownames(table)<- mod.lab.vec
  
  # table sorted by wAIC (Mc & IFI > 0.95 considered 'good' fit)
  summary.table<-table[order(table[,6],decreasing=TRUE),1:9]
  summary.table
  
  
  ########################################################
  ## boosted regression trees
  ########################################################
  
  ## BRT
  # full dataset first
  brt.fit <- gbm.step(datsem, gbm.x = attr(datsem, "names")[c(3:11)], gbm.y = attr(datsem, "names")[2], family="gaussian", tolerance = 0.0001, learning.rate = 0.0001, bag.fraction=0.8, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
  summary(brt.fit)
  gbm.plot(brt.fit)
  gbm.plot.fits(brt.fit)
  
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
  D2 # % deviance explained
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se
  print(c(CV.cor, CV.cor.se))
  eq.sp.points <- 100
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=9)
  ## output average predictions
  for (p in 1:9) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names

  # spatial autocorrelation test of residuals
  centroids <- read.table("country.centroids.csv", sep=",", header=T) # import country centroids (https://worldmap.harvard.edu/data/geonode:country_centroids_az8)
  centroids.reg <- merge(centroids, contreg, by="cntry.code")
  centroids.AFR1 <- subset(centroids.reg, cont=="AFR")
  centroids.AFR <- centroids.AFR1[!(centroids.AFR1$cntry.code %in% smisnat),]
  centroids.datsem <- centroids.AFR[(centroids.AFR$cntry.code %in% as.vector(datsem$cntry.code)),]
                             
  centr.List <- setNames(split(centroids.datsem[,c(2,3)], seq_len(nrow(centroids.datsem))), centroids.datsem$cntry.code)
  distMat <- outer(centr.List, centr.List, Vectorize(distVincentyEllipsoid), list(a=6378137, b=6356752.3142, f=1/298.257223563))
  distMat.inv <- 1/distMat
  diag(distMat.inv) <- 0
  
  # Moran's I
  Moran.I(as.vector(brt.fit$residuals), distMat.inv)
  
  # NCF
  detach(package:pgirmess)
  detach(package:ncf)
  rm(centr.NCF)
  rm(mean.dist)
  mean.dist <- mean(distMat[lower.tri(distMat, diag=F)]) / 111111 / 10
  centr.NCF <- correlog(x=as.vector(centroids.datsem$centr.long), y=as.vector(centroids.datsem$centr.lat), z=as.vector(brt.fit$residuals), increment=mean.dist, resamp=1000, latlon=F, quiet=F)
  centr.NCF
  par(mfrow=c(2,2))
  plot(as.numeric(names(centr.NCF$mean.of.class)),centr.NCF$correlation, xlab="distance (increments of 334 km)", ylab="Moran's I", pch=19, type="b")
  abline(h=0,lty=2)
  abline(v=7,lty=2)

  plot(as.numeric(names(centr.NCF$mean.of.class)), centr.NCF$p, pch=19, type="b",xlab="distance (increments of 334 km", ylab="p")
  abline(h=0.05, lty=2)
  
  # PGI
  centr.PGI <- correlog(coords=centroids.datsem[,c(2,3)], z=as.vector(brt.fit$residuals), method="Moran")
  centr.PGI
  plot(centr.PGI)
  centr.PGI <- correlog(coords=centroids.datsem[,c(2,3)], z=as.vector(brt.fit$residuals), method="Geary")
  centr.PGI
  plot(centr.PGI)
  par(mfrow=c(1,1))
  
  
  ##############################################################
  # randomise to test probablity of non-randomness per variable
  ##############################################################
  
  single.col.rnd = function(datfram, col.choose) { # function to randomise a single column of 'datfram' at a time (according to col 'col.choose')
    ext.col <- datfram[,which(colnames(datfram) == col.choose)]
    datrnd <- datfram
    datrnd[,which(colnames(datfram) == col.choose)] <- sample(ext.col, replace=F)
    return(datrnd)
  }
  
  iter <- 1000; itdiv <- iter/100
  var.vec <- colnames(datsem)[-c(1,2)]
  lvarvec <- length(var.vec)
  RMSE.obs <- rmse(datsem$stCHLTH,brt.fit$fitted) # root mean-squared error
  NRMSE.obs <- RMSE.obs/sd(brt.fit$fitted) # normalised root mean-squared error
  
  RMSE.pr <- NRMSE.pr <- rep(0,lvarvec)
  for (v in 1:lvarvec) {
    RMSE.pr.vec <- NRMSE.pr.vec <- rep(0,iter)
    for (i in 1:iter) {
      ran.it <- single.col.rnd(datsem, var.vec[v]) # randomised-column data frame
      brt.fit.ran <- gbm.step(ran.it, gbm.x = attr(ran.it, "names")[c(3:11)], gbm.y = attr(ran.it, "names")[2], family="gaussian", tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.8, tree.complexity = 2, tolerance.method = "auto", plot.main=F, plot.folds=F, max.trees=100000, verbose=F, silent=T)
      #RMSE.ran <- sqrt(sum((brt.fit.ran$fitted - brt.fit.ran$data$y)^2)/brt.fit.ran$nTrain)
      RMSE.ran <- rmse(datsem$stCHLTH, brt.fit.ran$fitted)
      NRMSE.ran <- RMSE.ran/sd(brt.fit.ran$fitted)
      RMSE.pr.vec[i] <- ifelse(RMSE.ran <= RMSE.obs, 1, 0)
      NRMSE.pr.vec[i] <- ifelse(NRMSE.ran <= NRMSE.obs, 1, 0)
      if (i %% itdiv==0) print(i) 
    }
    RMSE.pr[v] <- sum(RMSE.pr.vec)/iter
    NRMSE.pr[v] <- sum(NRMSE.pr.vec)/iter
    print("--------")
    print(paste("variable (",var.vec[v], ") = ", v, sep=""))
    print("--------")
  }
  out.pr <- data.frame(var.vec,RMSE.pr,NRMSE.pr)
  colnames(out.pr) <- c("var","rmsePr","nrmsePr")
  out.mrg <- merge(out.pr,brt.fit$contributions,by="var")
  out.sort <- out.mrg[order(out.mrg[,4],decreasing=T),1:4]
  out.sort
  
  
  
  ## re-do, but randomise all other variables except focal per iteration
  iter <- 1000; itdiv <- iter/100
  var.vec <- colnames(datsem)[-c(1,2)]
  lvarvec <- length(var.vec)
  RMSE.obs <- rmse(datsem$stCHLTH,brt.fit$fitted) # root mean-squared error
  NRMSE.obs <- RMSE.obs/sd(brt.fit$fitted) # normalised root mean-squared error
  
  RMSE.pr <- NRMSE.pr <- rep(0,lvarvec)
  for (v in 1:lvarvec) {
    RMSE.pr.vec <- NRMSE.pr.vec <- rep(0,iter)
    for (i in 1:iter) {
      ran.var.vec <- var.vec[-which(var.vec == var.vec[v])]
      ran.it <- randomise_data(datsem, variable=ran.var.vec, method="full_rand")
      brt.fit.ran <- gbm.step(ran.it, gbm.x = attr(ran.it, "names")[c(3:11)], gbm.y = attr(ran.it, "names")[2], family="gaussian", tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.8, tree.complexity = 2, tolerance.method = "auto", plot.main=F, plot.folds=F, max.trees=100000, verbose=F, silent=T)
      if(is.null(brt.fit.ran$fitted)==TRUE) {
        RMSE.pr.vec[i] <- NRMSE.pr.vec[i] <- NA
      }
      if(is.null(brt.fit.ran$fitted)==FALSE) {
        RMSE.ran <- rmse(datsem$stCHLTH, brt.fit.ran$fitted)
        NRMSE.ran <- RMSE.ran/sd(brt.fit.ran$fitted)
        RMSE.pr.vec[i] <- ifelse(RMSE.ran <= RMSE.obs, 1, 0)
        NRMSE.pr.vec[i] <- ifelse(NRMSE.ran <= NRMSE.obs, 1, 0)
        }
      if (i %% itdiv==0) print(i) 
    }
    RMSE.pr[v] <- sum(RMSE.pr.vec,na.rm=T)/(iter - length(which(is.na(RMSE.pr.vec)==TRUE)))
    NRMSE.pr[v] <- sum(NRMSE.pr.vec,na.rm=T)/(iter - length(which(is.na(NRMSE.pr.vec)==TRUE)))
    print("--------")
    print(paste("variable (",var.vec[v], ") = ", v, sep=""))
    print("--------")
  }
  out.pr <- data.frame(var.vec,RMSE.pr,NRMSE.pr)
  colnames(out.pr) <- c("var","rmsePr","nrmsePr")
  out.mrg <- merge(out.pr,brt.fit$contributions,by="var")
  out.sort <- out.mrg[order(out.mrg[,4],decreasing=T),1:4]
  out.sort
  
  
  
  ## bootstrap iteration loop to estimate prediction confidence interval
  iter <- 1000
  eq.sp.points <- 100
  
  # create storage arrays
  val.arr <- pred.arr <- array(data = 0, dim = c(eq.sp.points, 9, iter), dimnames=list(paste("x",1:eq.sp.points,sep=""), attr(datsem.boot, "names")[c(3:11)], paste("b",1:iter,sep="")))
  
  # create storage vectors
  D2.vec <- CV.cor.vec <- CV.cor.se.vec <- ENV.ri <- HS.ri <- GDP.ri <- GOV.ri <- HINV.ri <- H2O.ri <- FOOD.ri <- BF.ri <- PM25.ri <- rep(0,iter)

  # iterate
  for (b in 1:iter) {  # start b loop

    # boostrap data
    boot.sub <- sort(sample(x = 1:dim(datsem)[1], size = dim(datsem)[1], replace=TRUE))
    datsem.boot <- datsem[boot.sub,]
    
    brt.fit <- gbm.step(datsem.boot, gbm.x = attr(datsem.boot, "names")[c(3:11)], gbm.y = attr(datsem.boot, "names")[2], family="gaussian", tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.8, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000, silent=TRUE)
    summ.fit <- summary(brt.fit)

    ENV.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][1])]
    HS.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][2])]
    GDP.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][3])]
    GOV.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][4])]
    HINV.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][5])]
    H2O.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][6])]
    FOOD.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][7])]
    BF.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][8])]
    PM25.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(datsem.boot, "names")[c(3:11)][9])]
    
    gbm.plot(brt.fit)
    D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
    #D2 # % deviance explained
    D2.vec[b] <- D2
    CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
    CV.cor.vec[b] <- CV.cor
    CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
    CV.cor.se.vec[b] <- CV.cor.se

    RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=9)
    ## output average predictions
    for (p in 1:9) {
      RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
      RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
    }
    RESP.val.dat <- as.data.frame(RESP.val)
    colnames(RESP.val.dat) <- brt.fit$var.names
    RESP.pred.dat <- as.data.frame(RESP.pred)
    colnames(RESP.pred.dat) <- brt.fit$var.names
    
    val.arr[, , b] <- as.matrix(RESP.val.dat)
    pred.arr[, , b] <- as.matrix(RESP.pred.dat)
    
    print(b)
  }  # end b loop
  
  # kappa method to reduce effects of outliers on bootstrap estimates
  kappa <- 2
  kappa.n <- 5
  pred.update <- pred.arr
  
  for (k in 1:kappa.n) {
    boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
    boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
    
    for (z in 1:iter) {
      pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
    }

    print(k)
  }
  
  pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
  pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
  pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
  
  val.med <- apply(val.arr, MARGIN=c(1,2), median)
 
  par(mfrow=c(3,3)) 
  plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←worse) evironment (better →)" )
  lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
  lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

  plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←low) household size (high →)" )
  lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
  lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

  plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←poorer) GDP (higher" )
  lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
  lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red")

  plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←worse) governance (better→)" )
  lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
  lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

  plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←less) health investment (more→)" )
  lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
  lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

  plot(val.med[,6],pred.med[,6],type="l",ylim=c(min(pred.lo[,6]),max(pred.up[,6])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←less) improved water/sanitation (more→)" )
  lines(val.med[,6], pred.lo[,6], type="l", lty=2, col="red")
  lines(val.med[,6], pred.up[,6], type="l", lty=2, col="red")

  plot(val.med[,7],pred.med[,7],type="l",ylim=c(min(pred.lo[,7]),max(pred.up[,7])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←less) food supply (more→)" )
  lines(val.med[,7], pred.lo[,7], type="l", lty=2, col="red")
  lines(val.med[,7], pred.up[,7], type="l", lty=2, col="red")

  plot(val.med[,8],pred.med[,8],type="l",ylim=c(min(pred.lo[,8]),max(pred.up[,8])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←less) % breastfed (more→)" )
  lines(val.med[,8], pred.lo[,8], type="l", lty=2, col="red")
  lines(val.med[,8], pred.up[,8], type="l", lty=2, col="red")

  plot(val.med[,9],pred.med[,9],type="l",ylim=c(min(pred.lo[,9]),max(pred.up[,9])), lwd=2, ylab="(←worse) child health (better →)", xlab="(←lower) PM2.5 (higher→)" )
  lines(val.med[,9], pred.lo[,9], type="l", lty=2, col="red")
  lines(val.med[,9], pred.up[,9], type="l", lty=2, col="red")

  par(mfrow=c(1,1)) 
  
  # kappa method for output vectors
  D2.update <- D2.vec
  CV.cor.update <- CV.cor.vec
  CV.cor.se.update <- CV.cor.se.vec
  ENV.ri.update <- ENV.ri
  HS.ri.update <- HS.ri
  GDP.ri.update <- GDP.ri
  GOV.ri.update <- GOV.ri
  HINV.ri.update <- HINV.ri
  H2O.ri.update <- H2O.ri
  FOOD.ri.update <- FOOD.ri
  BF.ri.update <- BF.ri
  PM25.ri.update <- PM25.ri
  
  for (k in 1:kappa.n) {
    D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
    CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
    CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)

    ENV.mean <- mean(ENV.ri.update, na.rm=T); ENV.sd <- sd(ENV.ri.update, na.rm=T)
    HS.mean <- mean(HS.ri.update, na.rm=T); HS.sd <- sd(HS.ri.update, na.rm=T)
    GDP.mean <- mean(GDP.ri.update, na.rm=T); GDP.sd <- sd(GDP.ri.update, na.rm=T)
    GOV.mean <- mean(GOV.ri.update, na.rm=T); GOV.sd <- sd(GOV.ri.update, na.rm=T)
    HINV.mean <- mean(HINV.ri.update, na.rm=T); HINV.sd <- sd(HINV.ri.update, na.rm=T)
    H2O.mean <- mean(H2O.ri.update, na.rm=T); H2O.sd <- sd(H2O.ri.update, na.rm=T)
    FOOD.mean <- mean(FOOD.ri.update, na.rm=T); FOOD.sd <- sd(FOOD.ri.update, na.rm=T)
    BF.mean <- mean(BF.ri.update, na.rm=T); BF.sd <- sd(BF.ri.update, na.rm=T)
    PM25.mean <- mean(PM25.ri.update, na.rm=T); PM25.sd <- sd(PM25.ri.update, na.rm=T)
    
    for (u in 1:iter) {
      D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
      CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
      CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
      
      ENV.ri.update[u] <- ifelse((ENV.ri.update[u] < (ENV.mean-kappa*ENV.sd) | ENV.ri.update[u] > (ENV.mean+kappa*ENV.sd)), NA, ENV.ri.update[u])
      HS.ri.update[u] <- ifelse((HS.ri.update[u] < (HS.mean-kappa*HS.sd) | HS.ri.update[u] > (HS.mean+kappa*HS.sd)), NA, HS.ri.update[u])
      GDP.ri.update[u] <- ifelse((GDP.ri.update[u] < (GDP.mean-kappa*GDP.sd) | GDP.ri.update[u] > (GDP.mean+kappa*GDP.sd)), NA, GDP.ri.update[u])
      GOV.ri.update[u] <- ifelse((GOV.ri.update[u] < (GOV.mean-kappa*GOV.sd) | GOV.ri.update[u] > (GOV.mean+kappa*GOV.sd)), NA, GOV.ri.update[u])
      HINV.ri.update[u] <- ifelse((HINV.ri.update[u] < (HINV.mean-kappa*HINV.sd) | HINV.ri.update[u] > (HINV.mean+kappa*HINV.sd)), NA, HINV.ri.update[u])
      H2O.ri.update[u] <- ifelse((H2O.ri.update[u] < (H2O.mean-kappa*H2O.sd) | H2O.ri.update[u] > (H2O.mean+kappa*H2O.sd)), NA, H2O.ri.update[u])
      FOOD.ri.update[u] <- ifelse((FOOD.ri.update[u] < (FOOD.mean-kappa*FOOD.sd) | FOOD.ri.update[u] > (FOOD.mean+kappa*FOOD.sd)), NA, FOOD.ri.update[u])
      BF.ri.update[u] <- ifelse((BF.ri.update[u] < (BF.mean-kappa*BF.sd) | BF.ri.update[u] > (BF.mean+kappa*BF.sd)), NA, BF.ri.update[u])
      PM25.ri.update[u] <- ifelse((PM25.ri.update[u] < (PM25.mean-kappa*PM25.sd) | PM25.ri.update[u] > (PM25.mean+kappa*PM25.sd)), NA, PM25.ri.update[u])
    }
    
    print(k)
  }
  
  
  D2.med <- median(D2.update, na.rm=TRUE)
  D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
  D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
  print(c(D2.lo,D2.med,D2.up))
  
  CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
  CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
  CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
  print(c(CV.cor.lo,CV.cor.med,CV.cor.up))
  
  ENV.ri.lo <- quantile(ENV.ri.update, probs=0.025, na.rm=TRUE)
  ENV.ri.med <- median(ENV.ri.update, na.rm=TRUE)
  ENV.ri.up <- quantile(ENV.ri.update, probs=0.975, na.rm=TRUE)
  
  HS.ri.lo <- quantile(HS.ri.update, probs=0.025, na.rm=TRUE)
  HS.ri.med <- median(HS.ri.update, na.rm=TRUE)
  HS.ri.up <- quantile(HS.ri.update, probs=0.975, na.rm=TRUE)
  
  GDP.ri.lo <- quantile(GDP.ri.update, probs=0.025, na.rm=TRUE)
  GDP.ri.med <- median(GDP.ri.update, na.rm=TRUE)
  GDP.ri.up <- quantile(GDP.ri.update, probs=0.975, na.rm=TRUE)

  GOV.ri.lo <- quantile(GOV.ri.update, probs=0.025, na.rm=TRUE)
  GOV.ri.med <- median(GOV.ri.update, na.rm=TRUE)
  GOV.ri.up <- quantile(GOV.ri.update, probs=0.975, na.rm=TRUE)
  
  HINV.ri.lo <- quantile(HINV.ri.update, probs=0.025, na.rm=TRUE)
  HINV.ri.med <- median(HINV.ri.update, na.rm=TRUE)
  HINV.ri.up <- quantile(HINV.ri.update, probs=0.975, na.rm=TRUE)
  
  H2O.ri.lo <- quantile(H2O.ri.update, probs=0.025, na.rm=TRUE)
  H2O.ri.med <- median(H2O.ri.update, na.rm=TRUE)
  H2O.ri.up <- quantile(H2O.ri.update, probs=0.975, na.rm=TRUE)
  
  FOOD.ri.lo <- quantile(FOOD.ri.update, probs=0.025, na.rm=TRUE)
  FOOD.ri.med <- median(FOOD.ri.update, na.rm=TRUE)
  FOOD.ri.up <- quantile(FOOD.ri.update, probs=0.975, na.rm=TRUE)
  
  BF.ri.lo <- quantile(BF.ri.update, probs=0.025, na.rm=TRUE)
  BF.ri.med <- median(BF.ri.update, na.rm=TRUE)
  BF.ri.up <- quantile(BF.ri.update, probs=0.975, na.rm=TRUE)
  
  PM25.ri.lo <- quantile(PM25.ri.update, probs=0.025, na.rm=TRUE)
  PM25.ri.med <- median(PM25.ri.update, na.rm=TRUE)
  PM25.ri.up <- quantile(PM25.ri.update, probs=0.975, na.rm=TRUE)
  
  ri.lo <- c(ENV.ri.lo,HS.ri.lo,GDP.ri.lo,GOV.ri.lo,HINV.ri.lo,H2O.ri.lo,FOOD.ri.lo,BF.ri.lo,PM25.ri.lo)
  ri.med <- c(ENV.ri.med,HS.ri.med,GDP.ri.med,GOV.ri.med,HINV.ri.med,H2O.ri.med,FOOD.ri.med,BF.ri.med,PM25.ri.med)
  ri.up <- c(ENV.ri.up,HS.ri.up,GDP.ri.up,GOV.ri.up,HINV.ri.up,H2O.ri.up,FOOD.ri.up,BF.ri.up,PM25.ri.up)

  ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
  colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
  rownames(ri.out) <- attr(datsem.boot, "names")[c(3:11)]
  ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
  ri.sort
  

################################# 
## GLMMs
#################################
  ## functions
  delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
  weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
  
  # Set functions
  AICc <- function(...) {
    models <- list(...)
    num.mod <- length(models)
    AICcs <- numeric(num.mod)
    ns <- numeric(num.mod)
    ks <- numeric(num.mod)
    AICc.vec <- rep(0,num.mod)
    for (i in 1:num.mod) {
      if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
      if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
      AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
      ns[i] <- n
      ks[i] <- k
      AICc.vec[i] <- AICcs[i]
    }
    return(AICc.vec)
  }
  
  linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
    fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
    AIC.vec <- c(AICc(fit.full),AICc(fit.null))
    dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
    ER <- wAIC.vec[1]/wAIC.vec[2]
    r.sq.adj <- as.numeric(summary(fit.full)[9])
    return(ER)
  }
  
  ## add regional categories for mixed effects
  dat.GLMM <- merge(datsem, regions.AFR, by="cntry.code")
  str(dat.GLMM)
  # UN, AU, or WHO
  dat.GLMM$region <- factor(dat.GLMM$AU.region)
  table(dat.GLMM$region)
  
  ## define model set
  mod1 <- "stCHLTH ~ ENVr+ lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lBF + lPM25 + (1|region)"
  mod2 <- "stCHLTH ~ ENVr + lPM25 + (1|region)"
  mod3 <- "stCHLTH ~ ENVr + (1|region)"
  mod4 <- "stCHLTH ~ lPOPD + (1|region)"
  mod5 <- "stCHLTH ~ lGDP + (1|region)"
  mod6 <- "stCHLTH ~ lGDP + lHINV + (1|region)"
  mod7 <- "stCHLTH ~ lGDP + lGOV + (1|region)"
  mod8 <- "stCHLTH ~ lGOV + (1|region)"
  mod9 <- "stCHLTH ~ lHINV + (1|region)"
  mod10 <- "stCHLTH ~ lH2O + (1|region)"
  mod11 <- "stCHLTH ~ lFOOD + (1|region)"
  mod12 <- "stCHLTH ~ lBF + (1|region)"
  mod13 <- "stCHLTH ~ lFOOD + lBF + (1|region)"
  mod14 <- "stCHLTH ~ ENVr + lPOPD + (1|region)"
  mod15 <- "stCHLTH ~ ENVr + lH2O + (1|region)"
  mod16 <- "stCHLTH ~ ENVr + lPOPD + lH2O + (1|region)"
  mod17 <- "stCHLTH ~ lPM25 + (1|region)"
  mod18 <- "stCHLTH ~ 1 + (1|region)"
  
  ## Make model vector
  mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18)
  
  ## Define n.mod
  n.mod <- length(mod.vec)
  
  # Model fitting and logLik output loop
  Modnum <- length(mod.vec)
  LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
  mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
  mod.num <- seq(1,Modnum,1)

  for(i in 1:Modnum) {
    fit <- lmer(as.formula(mod.vec[i]),data=dat.GLMM, na.action=na.omit)
    assign(paste("fit",i,sep=""), fit)
    mod.list[[i]] <- fit
    LL.vec[i] <- as.numeric(logLik(fit))
    k.vec[i] <- attr(logLik(fit),"df")
    #AICc.vec[i] <- -2*LL.vec[i] + ((2*k.vec[i]*Modnum)/(Modnum-k.vec[i]-1))
    AICc.vec[i] <- r.squared(fit)$AIC
    
    Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
    Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
    
    print(i)
  }
  dAICc <- delta.IC(AICc.vec)
  wAICc <- weight.IC(dAICc)
  
  sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),round(Rm,4),round(Rc,4))
  colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
  row.names(sumtable) <- as.character(mod.vec)
  summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
  summary.table

  ## saturated residual diagnostic
  mod.saturated.glm <- "stCHLTH ~ ENVr+ lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lPM25"
  fit <- glm(as.formula(mod.saturated.glm),family=gaussian(link="identity"), data=dat.GLMM, na.action=na.omit)
  plot(fit)
  
  top.model <- "stCHLTH ~ lGDP + lGOV"
  topfit <- glm(as.formula(top.model),family=gaussian(link="identity"), data=dat.GLMM, na.action=na.omit)
  par(mfrow=c(1,2))
  plot1 <- termplot(topfit,terms="lGDP",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="GDP",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  plot1 <- termplot(topfit,terms="lGOV",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="GOV",ylab="CHILD HEALTH partial resid",cex.axis=1.2,cex.lab=1.2)
  par(mfrow=c(1,1))
  
 
   
  #####################################
  ## single environmental variable SEMs
  #####################################
  
  # STUNTING
  SNT.sc <- chealth.comb$st.stunt
  SNT.dat <- data.frame(chealth.comb$cntry.code, SNT.sc)
  colnames(SNT.dat) <- c("cntry.code", "RESP")
  
  # RESPIRATORY
  RSP.sc <- chealth.comb$st.resp
  RSP.dat <- data.frame(chealth.comb$cntry.code, RSP.sc)
  colnames(RSP.dat) <- c("cntry.code", "RESP")
  
  # DIARRHOEA
  DIA.sc <- chealth.comb$st.diar
  DIA.dat <- data.frame(chealth.comb$cntry.code, DIA.sc)
  colnames(DIA.dat) <- c("cntry.code", "RESP")

  # INFECTIOUS DISEASE
  INF.sc <- chealth.comb$st.inf
  INF.dat <- data.frame(chealth.comb$cntry.code, INF.sc)
  colnames(INF.dat) <- c("cntry.code", "RESP")

  # INJURY
  INJ.sc <- chealth.comb$st.inj
  INJ.dat <- data.frame(chealth.comb$cntry.code, INJ.sc)
  colnames(INJ.dat) <- c("cntry.code", "RESP")

  
  # combine with dat.sem
  # STUNTING
  SNT.comb <- merge(dat.sem, SNT.dat, by="cntry.code", all.x=T, all.y=T)
  SNT.all <- na.omit(SNT.comb[,-c(2,4,6,9)])
  SNT.cov <- cov(model.matrix(~ RESP + ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lBF + lPM25, data=SNT.all))[-1,-1]
  
  # combine with dat.sem
  # RESPIRATORY
  RSP.comb <- merge(dat.sem, RSP.dat, by="cntry.code", all.x=T, all.y=T)
  RSP.all <- na.omit(RSP.comb[,-c(2,4,6,9)])
  RSP.cov <- cov(model.matrix(~ RESP + ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lBF + lPM25, data=RSP.all))[-1,-1]
  
  # combine with dat.sem
  # DIARRHOEA
  DIA.comb <- merge(dat.sem, DIA.dat, by="cntry.code", all.x=T, all.y=T)
  DIA.all <- na.omit(DIA.comb[,-c(2,4,6,9)])
  DIA.cov <- cov(model.matrix(~ RESP + ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lBF + lPM25, data=DIA.all))[-1,-1]

  # combine with dat.sem
  # INFECTIOUS DISEASE
  INF.comb <- merge(dat.sem, INF.dat, by="cntry.code", all.x=T, all.y=T)
  INF.all <- na.omit(INF.comb[,-c(2,4,6,9)])
  INF.cov <- cov(model.matrix(~ RESP + ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lBF + lPM25, data=INF.all))[-1,-1]

  # combine with dat.sem
  # INJURY
  INJ.comb <- merge(dat.sem, INJ.dat, by="cntry.code", all.x=T, all.y=T)
  INJ.all <- na.omit(INJ.comb[,-c(2,4,6,9)])
  INJ.cov <- cov(model.matrix(~ RESP + ENVr + lPOPD + lGDP + lGOV + lHINV + lH2O + lFOOD + lBF + lPM25, data=INJ.all))[-1,-1]

  #################################################################
  ## Choose which variable (SNT, RSP, DIA, INF, INJ)
  cov.mat <- SNT.cov; var.dat <- SNT.all
  #cov.mat <- RSP.cov; var.dat <- RSP.all
  #cov.mat <- DIA.cov; var.dat <- DIA.all
  #cov.mat <- INF.cov; var.dat <- INF.all
  #cov.mat <- INJ.cov; var.dat <- INJ.all
  #################################################################

  
  # models
  # saturated model
  mod1 <- specifyModel(text="
                       ENVr -> RESP, beta1, NA
                       lPOPD -> RESP, beta2, NA
                       lGDP -> RESP, beta3, NA
                       lGOV -> RESP, beta4, NA
                       lHINV -> RESP, beta5, NA
                       lH2O -> RESP, beta6, NA
                       lFOOD -> RESP, beta7, NA
                       lBF -> RESP, beta8, NA
                       lPM25 <- RESP, beta9, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 2 ENV + PM25
  mod2 <- specifyModel(text="
                       ENVr -> RESP, beta1, NA
                       lPM25 -> RESP, beta2, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 3 ENV
  mod3 <- specifyModel(text="
                       ENVr -> RESP, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  
  # model 4 POPD
  mod4 <- specifyModel(text="
                       lPOPD -> RESP, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 5 GDP
  mod5 <- specifyModel(text="
                       lGDP -> RESP, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 6 GDP + HINV
  mod6 <- specifyModel(text="
                       lGDP -> RESP, beta1, NA
                       lHINV <- RESP, beta2, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 7 GDP + GOV
  mod7 <- specifyModel(text="
                       lGDP -> RESP, beta1, NA
                       lGOV <- RESP, beta2, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 8 GOV
  mod8 <- specifyModel(text="
                       lGOV -> RESP, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 9 HINV
  mod9 <- specifyModel(text="
                       lHINV -> RESP, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 10 H2O
  mod10 <- specifyModel(text="
                       lH2O -> RESP, beta1, NA
                       lPOPD -> lGOV, pi1, NA
                       lGDP -> lHINV, pi2, NA 
                       lPOPD -> ENVr, pi3, NA
                       ENVr -> lFOOD, pi4, NA
                       ENVr -> lH2O, pi5, NA
                       lFOOD -> lBF, pi6, NA
                       ENVr <-> ENVr, gam1, NA
                       lPOPD <-> lPOPD, gam2, NA
                       lGDP <-> lGDP, gam3, NA
                       lGOV <-> lGOV, gam4, NA
                       lHINV <-> lHINV, gam5, NA
                       lH2O <-> lH2O, gam6, NA
                       lFOOD <-> lFOOD, gam7, NA
                       lBF <-> lBF, gam8
                       lPM25 <-> lPM25, gam9
                       RESP <-> RESP, gam10, NA
                       ")
  
  # model 11 FOOD
  mod11 <- specifyModel(text="
                        lFOOD -> RESP, beta1, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        RESP <-> RESP, gam10, NA
                        ")
  
  # model 12 BREASTFEEDING
  mod12 <- specifyModel(text="
                        lBF -> RESP, beta1, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        RESP <-> RESP, gam10, NA
                        ")
  
  # model 13 FOOD & BREASTFEEDING
  mod13 <- specifyModel(text="
                        lFOOD -> RESP, beta1, NA
                        lBF -> RESP, beta2, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        RESP <-> RESP, gam10, NA
                        ")
  
  # model 14 ENV & POPD
  mod14 <- specifyModel(text="
                        ENVr -> RESP, beta1, NA
                        lPOPD -> RESP, beta2, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        RESP <-> RESP, gam10, NA
                        ")
  
  # model 15 ENV + H2O
  mod15 <- specifyModel(text="
                        ENVr -> RESP, beta1, NA
                        lH2O -> RESP, beta2, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        RESP <-> RESP, gam10, NA
                        ")
  
  # model 16 ENV + POPD + H2O
  mod16 <- specifyModel(text="
                        ENVr -> RESP, beta1, NA
                        lPOPD -> RESP, beta2, NA
                        lH2O -> RESP, beta3, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        RESP <-> RESP, gam10, NA
                        ")
  
  # model 17 PM25
  mod17 <- specifyModel(text="
                        lPM25 -> RESP, beta1, NA
                        lPOPD -> lGOV, pi1, NA
                        lGDP -> lHINV, pi2, NA 
                        lPOPD -> ENVr, pi3, NA
                        ENVr -> lFOOD, pi4, NA
                        ENVr -> lH2O, pi5, NA
                        lFOOD -> lBF, pi6, NA
                        ENVr <-> ENVr, gam1, NA
                        lPOPD <-> lPOPD, gam2, NA
                        lGDP <-> lGDP, gam3, NA
                        lGOV <-> lGOV, gam4, NA
                        lHINV <-> lHINV, gam5, NA
                        lH2O <-> lH2O, gam6, NA
                        lFOOD <-> lFOOD, gam7, NA
                        lBF <-> lBF, gam8
                        lPM25 <-> lPM25, gam9
                        RESP <-> RESP, gam10, NA
                        ")
  
  model1 <- "ENV+HS+GDP+GOV+HINV+H2O+FOOD+BF+PM25"
  model2 <- "ENV+PM25"
  model3 <- "ENV"
  model4 <- "HS"
  model5 <- "GDP"
  model6 <- "GDP+HINV"
  model7 <- "GDP+GOV"
  model8 <- "GOV"
  model9 <- "HINV"
  model10 <- "H2O"
  model11 <- "FOOD"
  model12 <- "BF"
  model13 <- "FOOD+BF"
  model14 <- "ENV+HS"
  model15 <- "ENV+H2O"
  model16 <- "ENV+HS+H2O"
  model17 <- "PM25"
  
  mod.lab.vec <- c(model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13,model14,model15,model16,model17)
  
  ## fit path models
  sem.mod1 <- sem(mod1, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod2 <- sem(mod2, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod3 <- sem(mod3, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod4 <- sem(mod4, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod5 <- sem(mod5, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod6 <- sem(mod6, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod7 <- sem(mod7, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod8 <- sem(mod8, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod9 <- sem(mod9, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod10 <- sem(mod10, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod11 <- sem(mod11, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod12 <- sem(mod12, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod13 <- sem(mod13, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod14 <- sem(mod14, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod15 <- sem(mod15, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod16 <- sem(mod16, cov.mat, dim(datsem)[1], standardized=T)
  sem.mod17 <- sem(mod17, cov.mat, dim(datsem)[1], standardized=T)
  
  # plots
  semPaths(sem.mod1, what="col", whatLabels="stand", layout="spring", label.cex=1.3, edge.label.cex=1.3, nCharNodes=7)
  
  # summaries
  sum1 <- summary(sem.mod1)
  sum2 <- summary(sem.mod2)
  sum3 <- summary(sem.mod3)
  sum4 <- summary(sem.mod4)
  sum5 <- summary(sem.mod5)
  sum6 <- summary(sem.mod6)
  sum7 <- summary(sem.mod7)
  sum8 <- summary(sem.mod8)
  sum9 <- summary(sem.mod9)
  sum10 <- summary(sem.mod10)
  sum11 <- summary(sem.mod11)
  sum12 <- summary(sem.mod12)
  sum13 <- summary(sem.mod13)
  sum14 <- summary(sem.mod14)
  sum15 <- summary(sem.mod15)
  sum16 <- summary(sem.mod16)
  sum17 <- summary(sem.mod17)
  
  # standardised coefficients
  stdCoef(sem.mod1)
  stdCoef(sem.mod2)
  stdCoef(sem.mod3)
  stdCoef(sem.mod4)
  stdCoef(sem.mod5)
  stdCoef(sem.mod6)
  stdCoef(sem.mod7)
  stdCoef(sem.mod8)
  stdCoef(sem.mod9)
  stdCoef(sem.mod10)
  stdCoef(sem.mod11)
  stdCoef(sem.mod12)
  stdCoef(sem.mod13)
  stdCoef(sem.mod14)
  stdCoef(sem.mod15)
  stdCoef(sem.mod16)
  stdCoef(sem.mod17)
  
  # SEM goodness-of-fit
  # a cutoff value close to .95 for TLI, BL89, CFI, RNI, and Gamma Hat;
  # a cutoff value close to .90 for Mc; a cutoff value close to .08 for SRMR;
  # and a cutoff value close to .06 for RMSEA
  # incremental fit index, also known as Bollen's IFI, is also relatively insensitive to sample size.
  # Values that exceed .90 are regarded as acceptable, although this index can exceed 1.
  # To compute the IFI, first the difference between the chi square of the independence model--
  # in which variables are uncorrelated--and the chi-square of the target model is calculated.
  # Next, the difference between the chi-square of the target model and the df for the target model is calculated.
  # The ratio of these values represents the IFI.
  # McDonald's Non-Centrality Index
  Mc.vec <- c(summaryGOF(sem.mod1, digits=4)$Mc,summaryGOF(sem.mod2, digits=4)$Mc,summaryGOF(sem.mod3, digits=4)$Mc,summaryGOF(sem.mod4, digits=4)$Mc,summaryGOF(sem.mod5, digits=4)$Mc,summaryGOF(sem.mod6, digits=4)$Mc,summaryGOF(sem.mod7, digits=4)$Mc,summaryGOF(sem.mod8, digits=4)$Mc,summaryGOF(sem.mod9, digits=4)$Mc,summaryGOF(sem.mod10, digits=4)$Mc,summaryGOF(sem.mod11, digits=4)$Mc,summaryGOF(sem.mod12, digits=4)$Mc,summaryGOF(sem.mod13, digits=4)$Mc,summaryGOF(sem.mod14, digits=4)$Mc,summaryGOF(sem.mod15, digits=4)$Mc,summaryGOF(sem.mod16, digits=4)$Mc,summaryGOF(sem.mod17, digits=4)$Mc)
  # Bollen's Incremental Fit Index
  IFI.vec <- c(summaryGOF(sem.mod1, digits=4)$IFI,summaryGOF(sem.mod2, digits=4)$IFI,summaryGOF(sem.mod3, digits=4)$IFI,summaryGOF(sem.mod4, digits=4)$IFI,summaryGOF(sem.mod5, digits=4)$IFI,summaryGOF(sem.mod6, digits=4)$IFI,summaryGOF(sem.mod7, digits=4)$IFI,summaryGOF(sem.mod8, digits=4)$IFI,summaryGOF(sem.mod9, digits=4)$IFI,summaryGOF(sem.mod10, digits=4)$IFI,summaryGOF(sem.mod11, digits=4)$IFI,summaryGOF(sem.mod12, digits=4)$IFI,summaryGOF(sem.mod13, digits=4)$IFI,summaryGOF(sem.mod14, digits=4)$IFI,summaryGOF(sem.mod15, digits=4)$IFI,summaryGOF(sem.mod16, digits=4)$IFI,summaryGOF(sem.mod17, digits=4)$IFI)
  # Expected cross-validation index (ECVI) (Browne & Cudeck 1992)
  ECVI.vec <- c(summaryGOF(sem.mod1, digits=4)$ECVI,summaryGOF(sem.mod2, digits=4)$ECVI,summaryGOF(sem.mod3, digits=4)$ECVI,summaryGOF(sem.mod4, digits=4)$ECVI,summaryGOF(sem.mod5, digits=4)$ECVI,summaryGOF(sem.mod6, digits=4)$ECVI,summaryGOF(sem.mod7, digits=4)$ECVI,summaryGOF(sem.mod8, digits=4)$ECVI,summaryGOF(sem.mod9, digits=4)$ECVI,summaryGOF(sem.mod10, digits=4)$ECVI,summaryGOF(sem.mod11, digits=4)$ECVI,summaryGOF(sem.mod12, digits=4)$ECVI,summaryGOF(sem.mod13, digits=4)$ECVI,summaryGOF(sem.mod14, digits=4)$ECVI,summaryGOF(sem.mod15, digits=4)$ECVI,summaryGOF(sem.mod16, digits=4)$ECVI,summaryGOF(sem.mod17, digits=4)$ECVI)
  # chi-square/df ratio (2 to 5 acceptable)
  CSDFR.vec <- c(summaryGOF(sem.mod1, digits=4)$chisq.df,summaryGOF(sem.mod2, digits=4)$chisq.df,summaryGOF(sem.mod3, digits=4)$chisq.df,summaryGOF(sem.mod4, digits=4)$chisq.df,summaryGOF(sem.mod5, digits=4)$chisq.df,summaryGOF(sem.mod6, digits=4)$chisq.df,summaryGOF(sem.mod7, digits=4)$chisq.df,summaryGOF(sem.mod8, digits=4)$chisq.df,summaryGOF(sem.mod9, digits=4)$chisq.df,summaryGOF(sem.mod10, digits=4)$chisq.df,summaryGOF(sem.mod11, digits=4)$chisq.df,summaryGOF(sem.mod12, digits=4)$chisq.df,summaryGOF(sem.mod13, digits=4)$chisq.df,summaryGOF(sem.mod14, digits=4)$chisq.df,summaryGOF(sem.mod15, digits=4)$chisq.df,summaryGOF(sem.mod16, digits=4)$chisq.df,summaryGOF(sem.mod17, digits=4)$chisq.df)
  
  # summary table
  mod.lab <- 1:17
  df.vec <- c(sum1$df, sum2$df, sum3$df, sum4$df, sum5$df, sum6$df, sum7$df, sum8$df, sum9$df, sum10$df, sum11$df, sum12$df, sum13$df, sum14$df, sum15$df, sum16$df, sum17$df)
  chisq.vec <- c(sum1$chisq, sum2$chisq, sum3$chisq, sum4$chisq, sum5$chisq, sum6$chisq, sum7$chisq, sum8$chisq, sum9$chisq, sum10$chisq, sum11$chisq, sum12$chisq, sum13$chisq, sum14$chisq, sum15$chisq, sum16$chisq, sum17$chisq)
  
  # BIC ranks
  delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
  weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
  
  BIC.vec <- c(sum1$BIC, sum2$BIC, sum3$BIC, sum4$BIC, sum5$BIC, sum6$BIC, sum7$BIC, sum8$BIC, sum9$BIC, sum10$BIC, sum11$BIC, sum12$BIC, sum13$BIC, sum14$BIC, sum15$BIC, sum16$BIC, sum17$BIC)
  dBIC.vec <- delta.IC(BIC.vec)
  wBIC.vec <- weight.IC(dBIC.vec)

  AIC.vec <- c(sum1$AIC, sum2$AIC, sum3$AIC, sum4$AIC, sum5$AIC, sum6$AIC, sum7$AIC, sum8$AIC, sum9$AIC, sum10$AIC, sum11$AIC, sum12$AIC, sum13$AIC, sum14$AIC, sum15$AIC, sum16$AIC, sum17$AIC)
  dAIC.vec <- delta.IC(AIC.vec)
  wAIC.vec <- weight.IC(dAIC.vec)

  # BIC results dataframe
  table<-cbind(mod.lab,df.vec,chisq.vec,BIC.vec,dBIC.vec,wBIC.vec,Mc.vec,IFI.vec,CSDFR.vec)
  colnames(table)<-c("model","df","chisq","BIC","dBIC","wBIC","Mc","IFI","CSDFR")
  rownames(table)<- mod.lab.vec
  
  # table sorted by wBIC (Mc & IFI > 0.95 considered 'good' fit); chi-square/df ratio 2-5 ok
  summary.table<-table[order(table[,6],decreasing=TRUE),1:9]
  summary.table

  # AIC results dataframe
  table<-cbind(mod.lab,df.vec,chisq.vec,AIC.vec,dAIC.vec,wAIC.vec,Mc.vec,IFI.vec,CSDFR.vec)
  colnames(table)<-c("model","df","chisq","AIC","dAIC","wAIC","Mc","IFI","CSDFR")
  rownames(table)<- mod.lab.vec
  
  # table sorted by wBIC (Mc & IFI > 0.95 considered 'good' fit)
  summary.table<-table[order(table[,6],decreasing=TRUE),1:9]
  summary.table

  

  ########################################################
  ## boosted regression trees
  ########################################################
  
  ## BRT
  brt.fit <- gbm.step(var.dat, gbm.x = attr(var.dat, "names")[c(2:10)], gbm.y = attr(var.dat, "names")[c(11)], family="gaussian", tolerance = 0.0001, learning.rate = 0.0001, bag.fraction=0.8, tree.complexity = 2, tolerance.method = "auto", plot.main=T, plot.folds=F, max.trees=100000)
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
  D2 # % deviance explained
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.se <- 100 * brt.fit$cv.statistics$correlation.se
  print(c(CV.cor, CV.cor.se))
  summary(brt.fit)
  gbm.plot(brt.fit, smooth=T)
  gbm.plot.fits(brt.fit)
  
  # raw values
  brt.fit$gbm.call$dataframe
  
  eq.sp.points <- 100
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=9)
  ## output average predictions
  for (p in 1:9) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names


## unit change estimation
# positive
p.lower <- -1.5977374; p.upper <- 2.49193458
r.lower <- 0.09266582; r.upper <- 0.1542678

xs <- c(p.lower,p.upper)
ys <- c(r.lower,r.upper)
plot(xs,ys,pch=3)
fit.lin <- lm(ys~xs)
abline(fit.lin,lty=2,col+"red")

# unit change in predictor
unit.x <- 1 + p.lower
unit.y <- as.numeric(coef(fit.lin)[1] + coef(fit.lin)[2]*unit.x)
diff.y <- unit.y - r.lower
round(diff.y * 100, 1)


# negative
p.lower <- -1.3420259; p.upper <- 2.7783756
r.lower <- 0.25611351; r.upper <- -0.044847

xs <- c(p.lower,p.upper)
ys <- c(r.lower,r.upper)
plot(xs,ys,pch=3)
fit.lin <- lm(ys~xs)
abline(fit.lin,lty=2,col+"red")

# unit change in predictor
unit.x <- 1 + p.lower
unit.y <- as.numeric(coef(fit.lin)[1] + coef(fit.lin)[2]*unit.x)
diff.y <- r.lower - unit.y
round(diff.y * 100, 1)
