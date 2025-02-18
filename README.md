# Africa Child Health
<img align="right" src="www/AfricaChildHealth.png" alt="Child Health" width="200" style="margin-top: 20px">

Data and R code to reproduce analyses examining socio-economic and environmental determinants of child-health outcomes among African nations

Accompanies the paper: <a href="https://globalecologyflinders.com/people/#DIRECTOR">Bradshaw, CJA</a>, <a href="https://www.zoology.ubc.ca/~otto/">SP Otto</a>, <a href="https://www.linkedin.com/in/alicia-annamalay-64876a43">AA Annamalay</a>, <a href="http://stanford.edu/~samhn/">S Heft-Neal</a>, <a href="https://www.rand.org/about/people/w/wagner_zachary.html">Z Wagner</a>, <a href="https://research-repository.uwa.edu.au/en/persons/peter-le-souef">PN Le Souëf</a>. 2019. Testing the socio-economic and environmental determinants of better child-health outcomes in Africa: a cross-sectional study among nations. <em>BMJ Open</em> 9: e029968. doi:<a href="http://doi.org/10.1136/bmjopen-2019-029968">10.1136/bmjopen-2019-029968</a>

## ABSTRACT
### Objective
We sought to test hypotheses regarding the principal correlates of child-health performance among African nations based on previous evidence collected at finer spatial scales.
### Design
Retrospective, cross-sectional study.
### Setting
All countries in Africa, excluding small-island nations.
### Primary and secondary outcome measures
We defined a composite child-health indicator for each country comprising the incidence of stunting, deaths from respiratory disease, deaths from diarrhoeal disease, deaths from other infectious disease and deaths from injuries for children aged under 5 years. We also compiled national-level data for Africa to test the effects of countrylevel water quality, air pollution, food supply, breast feeding, environmental performance, per capita wealth, healthcare investment, population density and governance quality on the child-health indicator.
### Results
Across nations, child health was lowest when water quality, improved sanitation, air quality and environmental performance were lowest. There was also
an important decline in child health as household size (a proxy for population density) increased. The remaining variables had only weak effects, but in the directions we hypothesised.
### Conclusions
These results emphasise the importance of continued investment in clean water and sanitation services, measures to improve air quality and efforts to
restrict further environmental degradation, to promote the UN’s Sustainable Development Goal 3 target to ‘... end preventable deaths of newborns and children under 5’ and Goal 6 to ‘... ensure access to water and sanitation for all’ by 2030.

<br>
<strong>Contact</strong>: Professor Corey J. A. Bradshaw <br>
<a href="http://globalecologyflinders.com">Global Ecology</a> <br>
<a href="mailto:corey.bradshaw@flinders.edu.au">e-mail</a>

03 September 2019

## Main R script
<code>AfricaCountryChildHealthGithub.R</code>

## Data files
The following files can be found in the <a href="https://github.com/cjabradshaw/AfricaChildHealth/tree/master/data">data</a> subdirectory:

- 'AFR.threat.csv'
- 'Africa.regions.csv'
- 'arableland.csv'
- 'cc2to3.csv'
- 'childhealthmetrics.csv'
- 'cntry.num.code.csv'
- 'continent.country.csv'
- 'country.centroids.csv'
- 'country_level_pm25.rds'
- 'cropland.csv'
- 'EFconsPC12.csv'
- 'exclbreastfedAFR.csv'
- 'exportvalue.csv'
- 'faocountrycode.csv'
- 'foodsupply.csv'
- 'forestloss.csv'
- 'freshwrem.csv'
- 'gdppcppp.csv'
- 'hc.exp.USDpc.csv'
- 'HME_DAH_DATABASE_1990_2015_Y2016M04D25.CSV.zip'
- 'household size Africa.csv'
- 'imp.san.pcaccess.csv'
- 'imp.water.pcaccess.csv'
- 'land.area.csv'
- 'landprot.csv'
- 'life.exp.birth.fem.csv'
- 'life.exp.birth.mal.csv'
- 'livestock.csv'
- 'megafaunaconserv.csv'
- 'overallgovernance.csv'
- 'pcemiss.csv'
- 'pcu5stunt.csv'
- 'pop.yr.csv'
- 'povgap.csv'
- 'PPP11.15.csv'
- 'respiratorydisease.csv'
- 'WHO.cntry.csv'

Note, you must unzip the file 'IHME_DAH_DATABASE_1990_2015_Y2016M04D25.CSV.zip' within the same directory to access the .csv file with the R script

## Required R libraries
- <code>boot</code>
- <code>biomod2</code>
- <code>sem</code>
- <code>semPlot</code>
- <code>semGOF</code>
- <code>dismo</code>
- <code>gbm</code>
- <code>ape</code>
- <code>geosphere</code>
- <code>ncf</code>
- <code>ModelMetrics</code>
- <code>pgirmess</code>
- <code>lme4</code>

## Source files
- <code>new_lmer_AIC_tables3.r</code>
- <code>r.squared.R</code>

<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" width="200" style="margin-top: 20px"></a>
<a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="200" style="margin-top: 20px"></a> 
<a href="https://github.com/FutureChildHealth"><img align="bottom-left" src="https://github.com/cjabradshaw/AfricaChildHealth/blob/master/www/FCHlogoFinaltransp.png" alt="Future Child Health logo" width="220" style="margin-top: 20px"></a>
