###############################################
#data preparation
###############################################

#set working directory
setwd("...")

#load data
Uni0 = read.csv('Unified clutch.csv', header=TRUE, sep=",", 
              na.strings="NA")

#load packages
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(brms)
library(rstan)
library(shinystan)
library(numDeriv) 
library(reshape2)
library(ggplot2)
library(tidybayes)
library(latex2exp)
library(RColorBrewer)
library(cowplot)
library(cmdstanr)
library(posterior)

#organize data & transform variables
{
  #exclude problematic cases
  Uni0.1 = subset(Uni0, Exclude != 'yes') #Deletes KY cases of experiments
  Uni0.1$Exclude = ifelse(Uni0.1$Clutch == 0, "yes", NA) #Checks for clutch = 0
  Uni0.2 = subset(Uni0.1, is.na(Uni0.1$Exclude)) #Deletes Lundy clutches = 0
  
  #Adjust variables to proper format and add some new ones
  Uni0.2$Clutch = as.numeric(Uni0.2$Clutch)
  Uni0.2$MinHatch = as.numeric(Uni0.2$MinHatch)
  Uni0.2$Prhatch = Uni0.2$MinHatch/Uni0.2$Clutch
  Uni0.2$Exclude = ifelse(Uni0.2$Prhatch > 1, "yes", NA)
  Uni0.3 = subset(Uni0.2, is.na(Uni0.2$Exclude))
  Uni0.3$Prhatch = as.numeric(Uni0.3$Prhatch)
  
  #Make Year a factor
  Uni0.3$Year = as.factor(Uni0.3$Year)
  Uni0.3$Attempt = as.numeric(Uni0.3$Attempt)
  Uni0.3$JulianFED = as.numeric(Uni0.3$JulianFED)
  Uni0.3$Banded = as.numeric(Uni0.3$Banded)
  Uni0.3$Nage = as.numeric(Uni0.3$Nage)
  Uni0.3$Mass = as.numeric(Uni0.3$Mass)
  
  #Mean-center Julian DFE within individual
  Uni0.3$UI = paste(Uni0.3$Fband)   # Create unique identifier (this is useful for all kind of things!)
  m = tapply(Uni0.3$JulianFED, Uni0.3$UI, mean, na.rm=TRUE)  # calculate mean of time per UI
  n = as.data.frame(m)  # transform into dataframe
  n$UI= names(m)			# add column UI to dataframe
  n$UI = as.factor(n$UI)	# tranfrom into categorical variable
  nl = match(Uni0.3$UI, n$UI)	
  m = vector(length=length(Uni0.3$Year)) 	# specifies number of rows (doesn't matter which column name you are using)
  x = n$m
  for (i in 1:length(Uni0.3$Year)) {
    m[i] = x[nl[i]] 
    Uni0.3$m = m}
  Uni0.3$DateMCW = Uni0.3$JulianFED - Uni0.3$m   # calculate difference between mean and raw score
  
  #Setting attempt so 1 = 0
  Uni0.3$AttemptAdj =Uni0.3$Attempt - 1
  #Setting Nestling age so 10 = 0
  Uni0.3$NageCent =Uni0.3$Nage - 10
  #Creating unique identifier for each brood
  Uni0.3$BroodID = paste(Uni0.3$Population, Uni0.3$NestID)
  #Creating indicator variable for first nestling in each attempt
  Uni0.3$Norder = 0 #Creates our indicator variable
  #dplyr code that retains row numbers in the data frame, groups by nest ID, and then adds a counter to Norder. 
  Uni1 = Uni0.3 %>% group_by(NestID) %>% dplyr::mutate(Norder = row_number())
  
  #remove missing female IDs in Lundy dataset
  Uni2 = subset(Uni1, Fband != "NA")
  
  ##Creating dataset with 1 observation of clutch and hatch per brood
  Uni2$PClutch = ifelse(Uni2$Norder == 1, Uni2$Clutch, NA) #Creating 1 clutch size value per attempt, for first nestling listed
  Uni2$PHatch = ifelse(Uni2$Norder == 1, Uni2$MinHatch, NA) #Creating 1 hatch size value per attempt, for first nestling listed
  Uni2$PBand = ifelse(Uni2$Norder == 1, Uni2$Banded, NA) #Creating 1 number banded value per attempt, for first nestling listed

  #splitting dataset into two populations
  Uni2KY = subset(Uni2, Population =='KY')
  Uni2LY = subset(Uni2, Population =='LY')

  ##Standardizing
  normFunc = function(x, na.rm=TRUE){(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)} #standardizing function

  ##Standardize environments
  Uni2KY$DateWSD = normFunc(Uni2KY$DateMCW)  #Standardizes across individuals the within-individual metric of date
  Uni2KY$AttemptSD = normFunc(Uni2KY$AttemptAdj)  #Standardizes the adjusted attempt (prior attempts).
  Uni2LY$DateWSD = normFunc(Uni2LY$DateMCW)  #Standardizes across individuals the within-individual metric of date
  Uni2LY$AttemptSD = normFunc(Uni2LY$AttemptAdj)  #Standardizes the adjusted attempt (prior attempts).

  #add observation-level index for 
  #binomial overdispersion
  Uni2KY$ovd = seq(1:nrow(Uni2KY))
  Uni2LY$ovd = seq(1:nrow(Uni2LY))
  
  #calculate average hatchling body mass per female attempt
  library(dplyr)
  avg_massKY = aggregate(Mass ~ Fband*Attempt, data = Uni2KY, FUN = mean, na.rm = TRUE)
  avg_massKY$Mass = (avg_massKY$Mass) 
  avg_massLY = aggregate(Mass ~ Fband*Attempt, data = Uni2LY, FUN = mean, na.rm = TRUE)
  avg_massLY$Mass = (avg_massLY$Mass) 
  
  #remove empty rows
  Uni2KY$Exclude = ifelse(is.na(Uni2KY$Exclude),0,1)
  Uni2KY = na.omit(Uni2KY)
  Uni2LY$Exclude = ifelse(is.na(Uni2LY$Exclude),0,1)
  Uni2LY = na.omit(Uni2LY)
  
  #add in body mass
  Uni2KY = left_join(Uni2KY, avg_massKY, by = c("Fband" = "Fband", "Attempt" = "Attempt"))
  Uni2KY$Mass_avg = Uni2KY$Mass.y / mean(Uni2KY$Mass.y)
  Uni2LY = left_join(Uni2LY, avg_massLY, by = c("Fband" = "Fband", "Attempt" = "Attempt"))
  Uni2LY$Mass_avg = Uni2LY$Mass.y  / mean(Uni2LY$Mass.y) 
  
  #scale hatchling count for ordinal model
  Uni2KY$Hatch = Uni2KY$PHatch + 1 #set 0 = 1
  Uni2LY$Hatch = Uni2LY$PHatch + 1 #set 0 = 1
}

#directory for cmdstan installation
set_cmdstan_path("...")

#general settings for Stan
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#custom function for cmdstan posteriors
extract = function(fit_obj) {
  vars = fit_obj$metadata()$stan_variables
  draws = posterior::as_draws_rvars(fit_obj$draws())
  
  lapply(vars, \(var_name){  
    posterior::draws_of(draws[[var_name]], with_chains = FALSE)
  }) |> setNames(vars)
}

#######################################################################
#Laying attempt plot
#######################################################################

KYa = reshape2::melt(Uni2KY, id.var = c("Fband", "Attempt", "Population"), measure.vars = "Clutch")
LYa = reshape2::melt(Uni2LY, id.var = c("Fband", "Attempt", "Population"), measure.vars = "Clutch")
dfa = rbind(KYa,LYa)


pattempt = 
  ggplot(dfa, aes(x = factor(value), group = Population, fill = Population)) +
    geom_bar(position = "dodge") +
  scale_fill_manual(values = c("seagreen3","lightblue4"))+
  facet_wrap(.~Attempt, labeller = labeller(Attempt = function(x) paste("Attempt", x)))+
  labs(x = "\n Clutch size", y = "Number of clutches \n")+
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "top",
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(color = guide_legend(override.aes = list(alpha=1)))

ggsave("attempt_plot.png",pattempt, dpi = 600, width = 5.5, height = 4, units = "in")

#######################################################################
#Underdispersion
#######################################################################

#visualize underdispersion in raw clutch size data
df = data.frame(count = c(Uni2KY$Clutch,
                          rpois(length(Uni2KY$Clutch), mean(Uni2KY$Clutch,na.rm=TRUE))),
                type  = c(rep("Raw data (clutch size)", length(Uni2KY$Clutch)),
                          rep("Poisson(4)", length(Uni2KY$Clutch))))
barplot(table(df$type, df$count), legend = TRUE, beside = TRUE)

#fit linear model (KY data)
m_nrm = bf(Clutch ~ 1 + DateWSD*AttemptSD + (1 + DateWSD*AttemptSD|Fband) + (1|Year),
           sigma ~  1 + (1|c|Fband)) + gaussian()
ud_nrm = brm(m_nrm, data = Uni2KY,
                 prior = c(prior("normal(0,1)", class = "Intercept"),
                           prior("normal(0,1)", class = "Intercept", dpar = "sigma"),
                           prior("normal(0,1)", class = "b"),
                           prior("exponential(2)", class = "sd"),
                           prior("exponential(2)", class = "sd", dpar = "sigma"),
                           prior("lkj(2)", class = "cor")),
             warmup = 500, iter=1500, chains = 4,
                 init = 0, control=list(adapt_delta=0.9, max_treedepth=10) ) 
saveRDS(ud_nrm,"ud_nrm.RDS")

ud_nrmL = brm(m_nrm, data = Uni2LY,
             prior = c(prior("normal(0,1)", class = "Intercept"),
                       prior("normal(0,1)", class = "Intercept", dpar = "sigma"),
                       prior("normal(0,1)", class = "b"),
                       prior("exponential(2)", class = "sd"),
                       prior("exponential(2)", class = "sd", dpar = "sigma"),
                       prior("lkj(2)", class = "cor")),
             warmup = 500, iter=1500, chains = 4,
             init = 0, control=list(adapt_delta=0.9, max_treedepth=10) ) 
saveRDS(ud_nrmL,"ud_nrmL.RDS")

ud_nrm = readRDS("ud_nrm.RDS")
ud_nrmL = readRDS("ud_nrmL.RDS")

#plot model predictions
pg = pp_check(ud_nrm)

#fit standard Poisson model (KY data)
m_pois = bf(Clutch ~ 1 + DateWSD*AttemptSD +
                 (1 + DateWSD*AttemptSD|Fband) + (1|Year) + (1|ovd)) + poisson()
ud_pois = brm(m_pois, data = Uni2KY,
                 prior = c(prior("normal(0,1)", class = "Intercept"),
                           prior("normal(0,1)", class = "b"),
                           prior("cauchy(0,1)", class = "sd"),
                           prior("lkj(2)", class = "cor")),
              warmup = 500, iter=1500, chains = 4,
                 inits = 0, control=list(adapt_delta=0.84, max_treedepth=10) ) 
saveRDS(ud_pois,"ud_pois.RDS")

ud_poisL = brm(m_pois, data = Uni2LY,
              prior = c(prior("normal(0,1)", class = "Intercept"),
                        prior("normal(0,1)", class = "b"),
                        prior("cauchy(0,1)", class = "sd"),
                        prior("lkj(2)", class = "cor")),
              warmup = 500, iter=1500, chains = 4,
              inits = 0, control=list(adapt_delta=0.84, max_treedepth=10) ) 
saveRDS(ud_poisL,"ud_poisL.RDS")

ud_pois = readRDS("ud_pois.RDS")
ud_poisL = readRDS("ud_poisL.RDS")

#plot model predictions
pp = pp_check(ud_pois)

#fit ordinal model (KY data)
m_c = bf(Clutch ~ 1 + DateWSD*AttemptSD + (1 + DateWSD*AttemptSD|c|Fband) + (1|Year), 
            disc ~ 1 + (1|c|Fband)) + cumulative("logit")
ud_ord = brm(m_c,  data = Uni2KY,
                           prior = c(prior("normal(0,1)", class = "Intercept"),
                           prior("normal(0,1)", class = "Intercept", dpar = "disc"),
                           prior("normal(0,1)", class = "b"),
                           prior("exponential(2)", class = "sd"),
                           prior("exponential(2)", class = "sd", dpar = "disc"),
                           prior("lkj(2)", class = "cor")),
                 warmup = 500, iter=1500, chains = 4,
                 init = 0, control=list(adapt_delta=0.90, max_treedepth=10)) 
saveRDS(ud_ord,"ud_ord.RDS")

ud_ordL = brm(m_c,  data = Uni2LY,
             prior = c(prior("normal(0,1)", class = "Intercept"),
                       prior("normal(0,1)", class = "Intercept", dpar = "disc"),
                       prior("normal(0,1)", class = "b"),
                       prior("exponential(2)", class = "sd"),
                       prior("exponential(2)", class = "sd", dpar = "disc"),
                       prior("lkj(2)", class = "cor")),
             warmup = 500, iter=1500, chains = 4,
             init = 0, control=list(adapt_delta=0.90, max_treedepth=10)) 
saveRDS(ud_ordL,"ud_ordL.RDS")

ud_ord = readRDS("ud_ord.RDS")
ud_ordL = readRDS("ud_ordL.RDS")

#plot model predictions
pp_gaus = pp_check(ud_nrm, nsamples = 100)
pp_pois = pp_check(ud_pois, nsamples = 100)
pp_ord = pp_check(ud_ord, resp = "clutch", nsamples = 100)

pp_gausL = pp_check(ud_nrmL, nsamples = 100)
pp_poisL = pp_check(ud_poisL, nsamples = 100)
pp_ordL = pp_check(ud_ordL, resp = "clutch", nsamples = 100)

#organize results
pred = data.frame(pp_ord$data)
obs = pred[pred$is_y==TRUE,]
obs$mod = "Data"
com = pred[pred$is_y==FALSE,]
com$mod = "Ordinal (ordered logit)"
pred = data.frame(pp_gaus$data)
gaus = pred[pred$is_y==FALSE,]
gaus$mod = "Gaussian"
pred = data.frame(pp_pois$data)
pois = pred[pred$is_y==FALSE,]
pois$mod = "Poisson"

pred = data.frame(pp_ordL$data)
obsL = pred[pred$is_y==TRUE,]
obsL$mod = "Data"
comL = pred[pred$is_y==FALSE,]
comL$mod = "Ordinal (ordered logit)"
pred = data.frame(pp_gausL$data)
gausL = pred[pred$is_y==FALSE,]
gausL$mod = "Gaussian"
pred = data.frame(pp_poisL$data)
poisL = pred[pred$is_y==FALSE,]
poisL$mod = "Poisson"

pdf = rbind(com, gaus, pois)
pdf$mod = factor(pdf$mod, levels = c("Gaussian","Poisson","Ordinal (ordered logit)"))

pdfL = rbind(comL, gausL, poisL)
pdfL$mod = factor(pdfL$mod, levels = c("Gaussian","Poisson","Ordinal (ordered logit)"))

pred.comp = 
ggplot(pdf, aes(x = value, group = interaction(rep_id,mod))) +
  geom_line(aes(color = mod),stat="density", size=1, alpha=0.24)+
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  scale_color_manual(values = c("#8BA794","#90c0df","#67e0a4"))+
    geom_line(inherit.aes = FALSE, data = obs, aes(x = value), 
            stat="density", size = 0.74, alpha=1, linetype = "solid", color = "black")+
  coord_cartesian(xlim = c(0,10))+
  labs(x = "Clutch size", y = "Density")+
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "top",
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))+
  guides(color = guide_legend(override.aes = list(alpha=1)))

pred.compL = 
  ggplot(pdfL, aes(x = value, group = interaction(rep_id,mod))) +
  geom_line(aes(color = mod),stat="density", size=1, alpha=0.24)+
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  scale_color_manual(values = c("#8BA794","#90c0df","#67e0a4"))+
  geom_line(inherit.aes = FALSE, data = obsL, aes(x = value), 
            stat="density", size = 0.74, alpha=1, linetype = "solid", color = "black")+
  coord_cartesian(xlim = c(0,10))+
  labs(x = "Clutch size", y = "Density")+
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "top",
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(color = guide_legend(override.aes = list(alpha=1)))

ggsave("pred_plot.png",pred.comp, dpi = 600, width = 5, height = 4, units = "in")
ggsave("pred_plotL.png",pred.compL, dpi = 600, width = 5, height = 4, units = "in")

#######################################################################
#Reaction norm analysis
#######################################################################

#model formulas
#note that brms uses a distinct parameterization of the latent residual
#where sigma = 1 / exp(disc), with disc being the modelled parameter
#random effect covariances in this model are, therefore, inversely proportional to those
#in our Stan model below where predictions are made directly on log(sigma)
m_c = bf(Clutch ~ 1 + DateWSD*AttemptSD + (1 + DateWSD*AttemptSD|c|Fband) + (1|Year), 
            disc ~ 1 + (1|c|Fband)) + cumulative("logit")
m_p =   bf(PBand| trials(Clutch) ~ 1 + DateWSD*AttemptSD + (1|c|Fband) + (1|Year) + (1|ovd)) + binomial("logit")
m_h =   bf(Hatch ~ 1 + DateWSD*AttemptSD + (1|c|Fband) + (1|Year),
           disc ~ 1) + cumulative("logit")
m_m =   bf(Mass_avg ~ 1 + DateWSD*AttemptSD + (1|c|Fband) + (1|Year)) + gaussian()

#model priors
mv_prior = c(prior("normal(0,1)", class = "Intercept", resp = "Clutch"),
             prior("normal(0,1)", class = "Intercept", resp = "Clutch", dpar = "disc"),
             prior("normal(0,1)", class = "b", resp = "Clutch"),
             prior("exponential(2)", class = "sd", resp = "Clutch"),
             prior("exponential(2)", class = "sd", resp = "Clutch", dpar = "disc"),
             
             prior("normal(0,1)", class = "Intercept", resp = "PBand"),
             prior("normal(0,1)", class = "b", resp = "PBand"),
             prior("exponential(2)", class = "sd", resp = "PBand"),
             
             prior("normal(0,1)", class = "Intercept", resp = "Hatch"),
             prior("normal(0,1)", class = "Intercept", resp = "Hatch", dpar = "disc"),
             prior("normal(0,1)", class = "b", resp = "Hatch"),
             prior("exponential(2)", class = "sd", resp = "Hatch"),
             
             prior("normal(0,1)", class = "Intercept", resp = "Massavg"),
             prior("normal(0,1)", class = "b", resp = "Massavg"),
             prior("exponential(2)", class = "sd", resp = "Massavg"),
             prior("exponential(2)", class = "sigma", resp = "Massavg"),
             prior("lkj(2)", class = "cor"))

#######################################################################
#KY
#######################################################################

#use cmdstan to speed up estimation 
#remove backend line to use rstan
ls_ord_KY = brm(m_c + m_p + m_h + m_m + set_rescor(FALSE), data = Uni2KY,
                 prior = mv_prior,
                 warmup = 500, iter=2000, chains = 4,
                 backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                 init = 0, control=list(adapt_delta=0.90, max_treedepth=10)) 

saveRDS(ls_ord_KY,"ls_ord_KY.RDS")
ls_ord_KY = readRDS("ls_ord_KY.RDS")

#summarize fixed effects
summary(ls_ord_KY, robust = T, prob = 0.9)
hypothesis(ls_ord_KY, "Clutch_DateWSD<0")
hypothesis(ls_ord_KY, "Clutch_AttemptSD>0")
hypothesis(ls_ord_KY, "Clutch_DateWSD:AttemptSD<0")
hypothesis(ls_ord_KY, "Hatch_DateWSD<0")
hypothesis(ls_ord_KY, "Hatch_AttemptSD<0")
hypothesis(ls_ord_KY, "Hatch_DateWSD:AttemptSD<0")
hypothesis(ls_ord_KY, "PBand_DateWSD>0")
hypothesis(ls_ord_KY, "PBand_AttemptSD<0")
hypothesis(ls_ord_KY, "PBand_DateWSD:AttemptSD>0")
hypothesis(ls_ord_KY, "Massavg_DateWSD<0")
hypothesis(ls_ord_KY, "Massavg_AttemptSD>0")
hypothesis(ls_ord_KY, "Massavg_DateWSD:AttemptSD<0")

#plot and save posterior predictive checks
#(check fit of model to observed data)
pc = pp_check(ls_ord_KY, resp = "Clutch")
ph = pp_check(ls_ord_KY, resp = "Hatch")
ps = pp_check(ls_ord_KY, resp = "PBand")
pm = pp_check(ls_ord_KY, resp = "Massavg")
ppc = plot_grid(pc, ph, ps, pm, nrow = 1)
pdf("ord_ppcheck_KY.pdf", width = 10, height = 2.4)
ppc
dev.off()

#plot average reaction norms
colfunc = colorRampPalette(c("chartreuse2", "darkgreen"))
colfunc(3)
library(ggplot2)
newdata = list(AttemptSD = c(-1.016, -0.129, 0.758)) #corresponds to attempt 1,2,3
{
#clutch size
cp = conditional_effects(ls_ord_KY, categorical = F, int_conditions = newdata, effects = "DateWSD:AttemptSD",
                         prob = 0.9, resp = "Clutch", plot = F)
plot(cp)

cpp = 
  plot(cp, plot = F)[[1]] + 
  scale_x_continuous(expand=c(0,0), breaks = c(-2,0,2), labels = c("Early", "Average", "Late"))+
  coord_cartesian(ylim = c(3,7))+ 
  ylab("Clutch size\n")+
  xlab("\nDate within season")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = rev(colfunc(3)))+
  scale_color_manual(values = rev(colfunc(3)))+
  guides(color = 'none', fill = 'none')

#hatchling count
hp = conditional_effects(ls_ord_KY, categorical = F, int_conditions = newdata, effects = "DateWSD:AttemptSD",
                         prob = 0.9, resp = "Hatch", plot = F)
#subtract 1 added arbitrarily for ordinal model
hp$`Hatch.Hatch_DateWSD:AttemptSD`$estimate__ = hp$`Hatch.Hatch_DateWSD:AttemptSD`$estimate__ - 1
hp$`Hatch.Hatch_DateWSD:AttemptSD`$lower__ = hp$`Hatch.Hatch_DateWSD:AttemptSD`$lower__ - 1
hp$`Hatch.Hatch_DateWSD:AttemptSD`$upper__ = hp$`Hatch.Hatch_DateWSD:AttemptSD`$upper__ - 1
hpp = 
  plot(hp, plot = F)[[1]] + 
  scale_x_continuous(expand=c(0,0), breaks = c(-2,0,2), labels = c("Early", "Average", "Late"))+
  scale_y_continuous(breaks = c(3,4,5,6))+
  coord_cartesian(ylim = c(3,6))+ 
  ylab("Hatchling count\n")+
  xlab("\nDate within season")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = rev(colfunc(3)))+
  scale_color_manual(values = rev(colfunc(3)))+
  guides(color = 'none', fill = 'none')

#hatchling survival
sp = conditional_effects(ls_ord_KY, categorical = F, int_conditions = newdata, effects = "DateWSD:AttemptSD",
                         prob = 0.9, resp = "PBand", plot = F, conditions = data.frame(Clutch = 1))
spp = 
  plot(sp, plot = F)[[1]] + 
  scale_x_continuous(expand=c(0,0), breaks = c(-2,0,2), labels = c("Early", "Average", "Late"))+
  coord_cartesian(ylim = c(0.4,0.8))+ 
  ylab("Nestling survival\n")+
  xlab("\nDate within season")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = rev(colfunc(3)))+
  scale_color_manual(values = rev(colfunc(3)))+
  guides(color = 'none', fill = 'none')

#hatchling mass
mp = conditional_effects(ls_ord_KY, categorical = F, int_conditions = newdata, effects = "DateWSD:AttemptSD",
                         prob = 0.9, resp = "Massavg", plot = F)
mpp = 
  plot(mp, plot = F)[[1]] + 
  scale_x_continuous(expand=c(0,0), breaks = c(-2,0,2), labels = c("Early", "Average", "Late"))+
  coord_cartesian(ylim = c(0.90, 1.10))+ 
  ylab("Nestling mass\n")+
  xlab("\nDate within season")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = rev(colfunc(3)))+
  scale_color_manual(values = rev(colfunc(3)))+
  guides(color = 'none', fill = 'none')
}
pp_ky = plot_grid(cpp, hpp, spp, mpp, ncol = 1, align = "v")
saveRDS(pp_ky,"pp_ky.RDS")

#variance
library(reshape2)
vc = VarCorr(ls_ord_KY, summary = F)
var_id = vc$Fband$sd^2
ky_var = reshape2::melt(var_id)
ky_var[ky_var$variable=="Massavg_Intercept","value"] = 
  ky_var[ky_var$variable=="Massavg_Intercept","value"] / vc$residual__$sd
ky_var$site = "KY"

median(var_id[,"Clutch_Intercept"])
quantile(var_id[,"Clutch_Intercept"], c(0.05,0.95))
median(var_id[,"Clutch_DateWSD"])
quantile(var_id[,"Clutch_DateWSD"], c(0.05,0.95))
median(var_id[,"Clutch_AttemptSD"])
quantile(var_id[,"Clutch_AttemptSD"], c(0.05,0.95))
median(var_id[,"Clutch_DateWSD:AttemptSD"])
quantile(var_id[,"Clutch_DateWSD:AttemptSD"], c(0.05,0.95))
median(var_id[,"disc_Clutch_Intercept"])
quantile(var_id[,"disc_Clutch_Intercept"], c(0.05,0.95))

median(var_id[,"Hatch_Intercept"])
quantile(var_id[,"Hatch_Intercept"], c(0.05,0.95))
median(var_id[,"PBand_Intercept"])
quantile(var_id[,"PBand_Intercept"], c(0.05,0.95))
median(var_id[,"Massavg_Intercept"] / vc$residual__$sd)
quantile(var_id[,"Massavg_Intercept"] / vc$residual__$sd, c(0.05,0.95))

#######################################################################
#LY
#######################################################################

#use cmdstan to speed up estimation 
#remove backend line to use rstan
ls_ord_LY = brm(m_c + m_p + m_h + m_m + set_rescor(FALSE), data = Uni2LY,
                prior = mv_prior,
                warmup = 500, iter = 2000, chains = 4,
                backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                init = 0, control=list(adapt_delta=0.90, max_treedepth=10)) 

saveRDS(ls_ord_LY,"ls_ord_lY.RDS")
ls_ord_LY = readRDS("ls_ord_LY.RDS")

#summarize fixed effects
summary(ls_ord_LY, robust = T, prob = 0.9)
hypothesis(ls_ord_LY, "Clutch_DateWSD<0")
hypothesis(ls_ord_LY, "Clutch_AttemptSD>0")
hypothesis(ls_ord_LY, "Clutch_DateWSD:AttemptSD<0")
hypothesis(ls_ord_LY, "Hatch_DateWSD<0")
hypothesis(ls_ord_LY, "Hatch_AttemptSD>0")
hypothesis(ls_ord_LY, "Hatch_DateWSD:AttemptSD<0")
hypothesis(ls_ord_LY, "PBand_DateWSD>0")
hypothesis(ls_ord_LY, "PBand_AttemptSD>0")
hypothesis(ls_ord_LY, "PBand_DateWSD:AttemptSD<0")
hypothesis(ls_ord_LY, "Massavg_DateWSD>0")
hypothesis(ls_ord_LY, "Massavg_AttemptSD>0")
hypothesis(ls_ord_LY, "Massavg_DateWSD:AttemptSD<0")

#plot and save posterior predictive checks
#(check fit of model to observed data)
pc = pp_check(ls_ord_LY, resp = "Clutch")
ph = pp_check(ls_ord_LY, resp = "Hatch")
ps = pp_check(ls_ord_LY, resp = "PBand")
pm = pp_check(ls_ord_LY, resp = "Massavg")
ppc = plot_grid(pc, ph, ps, pm, nrow = 1)
pdf("ord_ppcheck_LY.pdf", width = 10, height = 2.4)
ppc
dev.off()

#plot average reaction norms
colfunc = colorRampPalette(c("chartreuse2", "darkgreen"))
colfunc(3)
library(ggplot2)
{
#clutch size
cp = conditional_effects(ls_ord_LY, categorical = F, int_conditions = newdata, effects = "DateWSD:AttemptSD",
                         prob = 0.9, resp = "Clutch", plot = F)
cpp = 
  plot(cp, plot = F)[[1]] + 
  scale_x_continuous(expand=c(0,0), breaks = c(-2,0,2), labels = c("Early", "Average", "Late"))+
  coord_cartesian(ylim = c(3,7))+ 
  ylab("Clutch size\n")+
  xlab("\nDate within season")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = rev(colfunc(3)))+
  scale_color_manual(values = rev(colfunc(3)))+
  guides(color = 'none', fill = 'none')

#hatchling count
hp = conditional_effects(ls_ord_LY, categorical = F, int_conditions = newdata, effects = "DateWSD:AttemptSD",
                         prob = 0.9, resp = "Hatch", plot = F)
#subtract 1 added arbitrarily for ordinal model
hp$`Hatch.Hatch_DateWSD:AttemptSD`$estimate__ = hp$`Hatch.Hatch_DateWSD:AttemptSD`$estimate__ - 1
hp$`Hatch.Hatch_DateWSD:AttemptSD`$lower__ = hp$`Hatch.Hatch_DateWSD:AttemptSD`$lower__ - 1
hp$`Hatch.Hatch_DateWSD:AttemptSD`$upper__ = hp$`Hatch.Hatch_DateWSD:AttemptSD`$upper__ - 1
hpp = 
  plot(hp, plot = F)[[1]] + 
  scale_x_continuous(expand=c(0,0), breaks = c(-2,0,2), labels = c("Early", "Average", "Late"))+
  scale_y_continuous(breaks = c(3,4,5,6))+
  coord_cartesian(ylim = c(3,6))+ 
  ylab("Hatchling count\n")+
  xlab("\nDate within season")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = rev(colfunc(3)))+
  scale_color_manual(values = rev(colfunc(3)))+
  guides(color = 'none', fill = 'none')

#hatchling survival
sp = conditional_effects(ls_ord_LY, categorical = F, int_conditions = newdata, effects = "DateWSD:AttemptSD",
                         prob = 0.9, resp = "PBand", plot = F, conditions = data.frame(Clutch = 1))
spp = 
  plot(sp, plot = F)[[1]] + 
  scale_x_continuous(expand=c(0,0), breaks = c(-2,0,2), labels = c("Early", "Average", "Late"))+
  coord_cartesian(ylim = c(0.4,0.8))+ 
  ylab("Nestling survival\n")+
  xlab("\nDate within season")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = rev(colfunc(3)))+
  scale_color_manual(values = rev(colfunc(3)))+
  guides(color = 'none', fill = 'none')

#hatchling mass
mp = conditional_effects(ls_ord_LY, categorical = F, int_conditions = newdata, effects = "DateWSD:AttemptSD",
                         prob = 0.9, resp = "Massavg", plot = F)
mpp = 
  plot(mp, plot = F)[[1]] + 
  scale_x_continuous(expand=c(0,0), breaks = c(-2,0,2), labels = c("Early", "Average", "Late"))+
  coord_cartesian(ylim = c(0.9, 1.10))+ 
  ylab("Nestling mass\n")+
  xlab("\nDate within season")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = rev(colfunc(3)))+
  scale_color_manual(values = rev(colfunc(3)))+
  guides(color = 'none', fill = 'none')
}
pp_ly = plot_grid(cpp, hpp, spp, mpp, ncol = 1, align = "v")
saveRDS(pp_ly,"pp_ly.RDS")

#variance
library(reshape2)
vc = VarCorr(ls_ord_LY, summary = F)
var_id = vc$Fband$sd^2
ly_var = reshape2::melt(var_id)
ly_var[ly_var$variable=="Massavg_Intercept","value"] = 
  ly_var[ly_var$variable=="Massavg_Intercept","value"] / vc$residual__$sd
ly_var$site = "LY"

median(var_id[,"Clutch_Intercept"])
quantile(var_id[,"Clutch_Intercept"], c(0.05,0.95))
median(var_id[,"Clutch_DateWSD"])
quantile(var_id[,"Clutch_DateWSD"], c(0.05,0.95))
median(var_id[,"Clutch_AttemptSD"])
quantile(var_id[,"Clutch_AttemptSD"], c(0.05,0.95))
median(var_id[,"Clutch_DateWSD:AttemptSD"])
quantile(var_id[,"Clutch_DateWSD:AttemptSD"], c(0.05,0.95))
median(var_id[,"disc_Clutch_Intercept"])
quantile(var_id[,"disc_Clutch_Intercept"], c(0.05,0.95))

median(var_id[,"Hatch_Intercept"])
quantile(var_id[,"Hatch_Intercept"], c(0.05,0.95))
median(var_id[,"PBand_Intercept"])
quantile(var_id[,"PBand_Intercept"], c(0.05,0.95))
median(var_id[,"Massavg_Intercept"] / vc$residual__$sd)
quantile(var_id[,"Massavg_Intercept"] / vc$residual__$sd, c(0.05,0.95))


#######################################################################
#compare, combine and plot
#######################################################################

#compare mean clutch size
pKY = rowSums(fitted(ls_ord_KY, newdata = data.frame(DateWSD = 0, AttemptSD = 0), 
              resp = "Clutch", re_formula = .~ NULL, summary = F)[,1,] * 1:9)
pLY = rowSums(fitted(ls_ord_LY, newdata = data.frame(DateWSD = 0, AttemptSD = 0), 
              resp = "Clutch", re_formula = .~ NULL, summary = F)[,1,] * 1:7)
diffm = pKY - pLY
mean(diffm); quantile(diffm, c(0.05,0.95))

#compare slopes
postKY = posterior_samples(ls_ord_KY)
postLY = posterior_samples(ls_ord_LY)
diff1 = postLY$b_Clutch_DateWSD - postKY$b_Clutch_DateWSD
diff2 = postLY$b_Clutch_AttemptSD - postKY$b_Clutch_AttemptSD 
diff3 = postLY$`b_Clutch_DateWSD:AttemptSD` - postKY$`b_Clutch_DateWSD:AttemptSD`
mean(diff1); quantile(diff1, c(0.05,0.95)); sum(diff1<0)/length(diff1)
mean(diff2); quantile(diff2, c(0.05,0.95)); sum(diff2>0)/length(diff2)
mean(diff3); quantile(diff3, c(0.05,0.95)); sum(diff3<0)/length(diff3)

#combine plots
pp_ky = readRDS("pp_ky.RDS")
pp_ly = readRDS("pp_ly.RDS")
pp_both = plot_grid(pp_ky, pp_ly, ncol = 2, align = "h")
saveRDS(pp_both, "pp_both.RDS")
png("avg_rn.png", width = 6, height = 8, unit = "in", res = 600)
pp_both
dev.off()

#plot
library(tidybayes)
both_var = rbind(ky_var, ly_var)
both_var$variable = factor(both_var$variable, 
                           levels = c("Massavg_Intercept", "PBand_Intercept", "Hatch_Intercept",
                                     "disc_Clutch_Intercept", "Clutch_DateWSD:AttemptSD",
                                     "Clutch_DateWSD", "Clutch_AttemptSD", "Clutch_Intercept"))

var.p =
  ggplot(both_var, aes(x = value, y = variable, group = site, color = site, fill = site))+
  stat_pointinterval(position = position_dodge(0.24), .width=c(0.9), size = 4)+
  xlab("\nAmong-individual variance")+
  ylab(" ")+
  scale_y_discrete(labels = c("Nestling mass", "Nestling survival", "Hatchling count",
                              "Clutch: residuals", "Clutch: attempt x date",
                              "Clutch: date", "Clutch: attempt", "Clutch: intercept"))+
  scale_fill_manual(values = c("deeppink1","purple"))+
  scale_color_manual(values = c("deeppink1","purple"))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = 'none', color = 'none')

ggsave("var_p.png", var.p, dpi = 600, height = 7, width = 5.5, units = "in")

library(grid)
v_b = plot_grid(nullGrob(), pp_both, nullGrob(), ncol = 1, rel_heights = c(0.04,0.9,0.04))
v_p = plot_grid(nullGrob(),var.p, nullGrob(), ncol = 1, rel_heights=c(0.04,0.9,0.04))
pp_all = plot_grid(v_b, v_p, ncol = 2, rel_widths = c(0.6,0.4))
saveRDS(pp_all, "pp_all.RDS")
png("figure 1_cs.png", width = 10, height = 8, unit = "in", res = 600)
pp_all
dev.off()

#########################################################
#Nonlinear selection analysis
#########################################################

#data preparation for Stan
#Kentucky########################################################

#prepare data for Stan
df = Uni2KY[,c("Clutch","Hatch", "PBand", "Mass_avg", "Year", 
               "DateWSD", "AttemptSD", "Fband")]
df = na.omit(df)
df$id = as.integer(factor(df$Fband, levels=unique(df$Fband)))
df$year_id = as.integer(factor(df$Year, levels=unique(df$Year)))

stan_data_KY = list(I = length(unique(df$Fband)), # subjects
                    N = nrow(df), # data frame rows
                    nyear = length(unique(df$year_id)), # years
                    id = df$id, # index for subjects
                    year_id = df$year_id, #index for years
                    
                    nthres_clutch = max(df$Clutch) - 1,
                    nthres_hatch = max(df$Hatch) - 1,
                    X_clutch = model.matrix(~ DateWSD*AttemptSD, df)[,-1],
                    k_clutch = ncol(model.matrix(~ DateWSD*AttemptSD, df)[,-1]),
                    X_hatch = model.matrix(~ DateWSD*AttemptSD, df)[,-1],
                    k_hatch = ncol(model.matrix(~ DateWSD*AttemptSD, df)[,-1]),
                    X_survive = model.matrix(~ DateWSD*AttemptSD, df)[,-1],
                    k_survive = ncol(model.matrix(~ DateWSD*AttemptSD, df)[,-1]),
                    X_mass = model.matrix(~ DateWSD*AttemptSD, df)[,-1],
                    k_mass = ncol(model.matrix(~ DateWSD*AttemptSD, df)[,-1]),
                    
                    #random slope predictors
                    dateWSD = df$DateWSD,
                    attemptSD = df$AttemptSD,
                    date_x_attempt = df$DateWSD * df$AttemptSD,

                    #response variables
                    clutch = df$Clutch,
                    hatch = df$Hatch,
                    survive = df$PBand,
                    mass = df$Mass_avg)

#Lundy###################################################

df = Uni2LY[,c("Clutch","Hatch", "PBand", "Mass_avg", "Year", 
               "DateWSD", "AttemptSD", "Fband")]
df = na.omit(df)
df$id = as.integer(factor(df$Fband, levels=unique(df$Fband)))
df$year_id = as.integer(factor(df$Year, levels=unique(df$Year)))

stan_data_LY = list(I = length(unique(df$Fband)), # subjects
                    N = nrow(df), # data frame rows
                    nyear = length(unique(df$year_id)), # years
                    id = df$id, # index for subjects
                    year_id = df$year_id, #index for years
                    
                    #without year
                    nthres_clutch = max(df$Clutch) - 1,
                    nthres_hatch = max(df$Hatch) - 1,
                    X_clutch = model.matrix(~ DateWSD*AttemptSD, df)[,-1],
                    k_clutch = ncol(model.matrix(~ DateWSD*AttemptSD, df)[,-1]),
                    X_hatch = model.matrix(~ DateWSD*AttemptSD, df)[,-1],
                    k_hatch = ncol(model.matrix(~ DateWSD*AttemptSD, df)[,-1]),
                    X_survive = model.matrix(~ DateWSD*AttemptSD, df)[,-1],
                    k_survive = ncol(model.matrix(~ DateWSD*AttemptSD, df)[,-1]),
                    X_mass = model.matrix(~ DateWSD*AttemptSD, df)[,-1],
                    k_mass = ncol(model.matrix(~ DateWSD*AttemptSD, df)[,-1]),
                    
                    #random slope predictors
                    dateWSD = df$DateWSD,
                    attemptSD = df$AttemptSD, #centered on most common attempt number
                    date_x_attempt = df$DateWSD * df$AttemptSD,

                    #response variables
                    clutch = df$Clutch,
                    hatch = df$Hatch,
                    survive = df$PBand,
                    mass = df$Mass_avg)

#########################################################
#Estimate NL model
#########################################################

#final model
m_NLS = cmdstan_model(stan_file = "NLS clutch_cumord.stan", 
                      stanc_options = list("O1"))

#full model
m_NLSf = cmdstan_model(stan_file = "NLS clutch_cumord_full.stan", 
                      stanc_options = list("O1"))

#unadjusted model
m_NLSua = cmdstan_model(stan_file = "NLS clutch_cumord_uadj.stan", 
                      stanc_options = list("O1"))

#date selection effect model
m_NLSds = cmdstan_model(stan_file = "NLS clutch_cumord_ds.stan", 
                        stanc_options = list("O1"))

#Kentucky#################################################
fitKY = m_NLS$sample(
        data = stan_data_KY,
        iter_sampling = 500,
        iter_warmup = 1000,
        init = 1e-4,
        chains = 4,  
        parallel_chains = 4,
        adapt_delta = 0.80,
        refresh = 10) 

postKY = extract(fitKY)
saveRDS(postKY, "NLS_KY_m.RDS")
launch_shinystan(fitKY)

#robustness checks
{

#check if inclusion of RN slopes is warranted
fitKYf = m_NLSf$sample(
  data = stan_data_KY,
  iter_sampling = 500,
  iter_warmup = 1000,
  init = 1e-4,
  chains = 4,  
  parallel_chains = 4,
  adapt_delta = 0.80,
  refresh = 10) 

NLS_KYf = as_draws_rvars(fitKYf$draws())
saveRDS(NLS_KYf, "NLS_KYf_m.RDS")

sum_fun = function(x) {return(c(paste0(round(median(x),2), 
                                       ",", round(quantile(x, c(0.05,0.95))[1],2), 
                                       ",", round(quantile(x, c(0.05,0.95))[2],2))))}
b_ky = data.frame(
           c(apply(draws_of(NLS_KYf$b_h), 2, sum_fun),
           apply(draws_of(NLS_KYf$b_s), 2, sum_fun),
           apply(draws_of(NLS_KYf$b_m), 2, sum_fun)))
write.csv(b_ky, "fullcoef_bky.csv")

q_ky = data.frame(
  c(apply(draws_of(NLS_KYf$q_h), 2, sum_fun),
    apply(draws_of(NLS_KYf$q_s), 2, sum_fun),
    apply(draws_of(NLS_KYf$q_m), 2, sum_fun)))
write.csv(q_ky, "fullcoef_qky.csv")

#check if exclusion of adjusted effects is warranted
fitKYua = m_NLSua$sample(
  data = stan_data_KY,
  iter_sampling = 500,
  iter_warmup = 1000,
  init = 1e-4,
  chains = 4,  
  parallel_chains = 4,
  adapt_delta = 0.80,
  refresh = 10) 

NLS_KYua = as_draws_rvars(fitKYua$draws())
saveRDS(NLS_KYua, "NLS_KYua_m.RDS")

sum_fun = function(x) {return(c(paste0(round(median(x),2), 
                                       ",", round(quantile(x, c(0.05,0.95))[1],2), 
                                       ",", round(quantile(x, c(0.05,0.95))[2],2))))}
b_ky = data.frame(
  c(apply(postKY$b_h - draws_of(NLS_KYua$b_h), 2, sum_fun),
    apply(postKY$b_s - draws_of(NLS_KYua$b_s), 2, sum_fun),
    apply(postKY$b_m - draws_of(NLS_KYua$b_m), 2, sum_fun)))
write.csv(b_ky, "coefdiff_bky.csv")

q_ky = data.frame(
  c(apply(postKY$q_h - draws_of(NLS_KYua$q_h), 2, sum_fun),
    apply(postKY$q_s - draws_of(NLS_KYua$q_s), 2, sum_fun),
    apply(postKY$q_m - draws_of(NLS_KYua$q_m), 2, sum_fun)))
write.csv(q_ky, "coefdiff_qky.csv")

#check for date effects on selection
fitKYd = m_NLSds$sample(
  data = stan_data_KY,
  iter_sampling = 500,
  iter_warmup = 1000,
  init = 1e-4,
  chains = 4,  
  parallel_chains = 4,
  adapt_delta = 0.80,
  refresh = 10) 

NLS_KYd = as_draws_rvars(fitKYd$draws())
saveRDS(NLS_KYd, "NLS_KYd_m.RDS")

sum_fun = function(x) {return(c(paste0(round(median(x),2), 
                                       ",", round(quantile(x, c(0.05,0.95))[1],2), 
                                       ",", round(quantile(x, c(0.05,0.95))[2],2))))}
b_ky = data.frame(
  c(apply(draws_of(NLS_KYd$db_h), 2, sum_fun),
    apply(draws_of(NLS_KYd$db_s), 2, sum_fun),
    apply(draws_of(NLS_KYd$db_m), 2, sum_fun)))
write.csv(b_ky, "coefd_bky.csv")

q_ky = data.frame(
  c(apply(draws_of(NLS_KYd$dq_h), 2, sum_fun),
    apply(draws_of(NLS_KYd$dq_s), 2, sum_fun),
    apply(draws_of(NLS_KYd$dq_m), 2, sum_fun)))
write.csv(q_ky, "coefd_qky.csv")
}

#Lundy#################################################
fitLY = m_NLS$sample(
  data = stan_data_LY,
  iter_sampling = 500,
  iter_warmup = 1000,
  init = 1e-4,
  chains = 4,  
  parallel_chains = 4,
  adapt_delta = 0.80,
  refresh = 10) 

postLY = extract(fitLY)
saveRDS(postLY, "NLS_LY_m.RDS")
launch_shinystan(fitLY)

#robustness checks
{
#check if inclusion of RN slopes is warranted
fitLYf = m_NLSf$sample(
  data = stan_data_LY,
  iter_sampling = 500,
  iter_warmup = 1000,
  init = 1e-4,
  chains = 4,  
  parallel_chains = 4,
  adapt_delta = 0.80,
  refresh = 10) 

NLS_LYf = as_draws_rvars(fitLYf$draws())
saveRDS(NLS_LYf, "NLS_LYf_m.RDS")

sum_fun = function(x) {return(c(paste0(round(median(x),2), 
                                       ",", round(quantile(x, c(0.05,0.95))[1],2), 
                                       ",", round(quantile(x, c(0.05,0.95))[2],2))))}
b_ly = data.frame(
  c(apply(draws_of(NLS_LYf$b_h), 2, sum_fun),
    apply(draws_of(NLS_LYf$b_s), 2, sum_fun),
    apply(draws_of(NLS_LYf$b_m), 2, sum_fun)))
write.csv(b_ly, "fullcoef_bly.csv")

q_ly = data.frame(
  c(apply(draws_of(NLS_LYf$q_h), 2, sum_fun),
    apply(draws_of(NLS_LYf$q_s), 2, sum_fun),
    apply(draws_of(NLS_LYf$q_m), 2, sum_fun)))
write.csv(q_ly, "fullcoef_qly.csv")

#check if exclusion of adjusted effects is warranted
fitLYua = m_NLSua$sample(
  data = stan_data_LY,
  iter_sampling = 500,
  iter_warmup = 1000,
  init = 1e-4,
  chains = 4,  
  parallel_chains = 4,
  adapt_delta = 0.80,
  refresh = 10) 

NLS_LYua = as_draws_rvars(fitLYua$draws())
saveRDS(NLS_LYua, "NLS_LYua_m.RDS")

sum_fun = function(x) {return(c(paste0(round(median(x),2), 
                                       ",", round(quantile(x, c(0.05,0.95))[1],2), 
                                       ",", round(quantile(x, c(0.05,0.95))[2],2))))}
b_ly = data.frame(
  c(apply(postLY$b_h - draws_of(NLS_LYua$b_h), 2, sum_fun),
    apply(postLY$b_s - draws_of(NLS_LYua$b_s), 2, sum_fun),
    apply(postLY$b_m - draws_of(NLS_LYua$b_m), 2, sum_fun)))
write.csv(b_ly, "coefdiff_bly.csv")

q_ly = data.frame(
  c(apply(postLY$q_h - draws_of(NLS_LYua$q_h), 2, sum_fun),
    apply(postLY$q_s - draws_of(NLS_LYua$q_s), 2, sum_fun),
    apply(postLY$q_m - draws_of(NLS_LYua$q_m), 2, sum_fun)))
write.csv(q_ly, "coefdiff_qly.csv")

#check for date effects on selection
fitLYd = m_NLSds$sample(
  data = stan_data_LY,
  iter_sampling = 500,
  iter_warmup = 1000,
  init = 1e-4,
  chains = 4,  
  parallel_chains = 4,
  adapt_delta = 0.80,
  refresh = 10) 

NLS_LYd = as_draws_rvars(fitLYd$draws())
saveRDS(NLS_LYd, "NLS_LYd_m.RDS")

sum_fun = function(x) {return(c(paste0(round(median(x),2), 
                                       ",", round(quantile(x, c(0.05,0.95))[1],2), 
                                       ",", round(quantile(x, c(0.05,0.95))[2],2))))}
b_ly = data.frame(
  c(apply(draws_of(NLS_LYd$db_h), 2, sum_fun),
    apply(draws_of(NLS_LYd$db_s), 2, sum_fun),
    apply(draws_of(NLS_LYd$db_m), 2, sum_fun)))
write.csv(b_ly, "coefd_bly.csv")

q_ly = data.frame(
  c(apply(draws_of(NLS_LYd$dq_h), 2, sum_fun),
    apply(draws_of(NLS_LYd$dq_s), 2, sum_fun),
    apply(draws_of(NLS_LYd$dq_m), 2, sum_fun)))
write.csv(q_ly, "coefd_qly.csv")
}

###########################################################
#Distribution of individual RNs
###########################################################

#extract posteriors
samplesKY = readRDS("NLS_KY_m.RDS")
samplesLY = readRDS("NLS_LY_m.RDS")

#visualize the distributions of random intercepts and residuals
muKY = apply(samplesKY$I_c[,,1],2, median)
sigKY = apply(samplesKY$I_c[,,2],2, median)

muLY = apply(samplesLY$I_c[,,1], 2, median)
sigLY = apply(samplesLY$I_c[,,2], 2, median)

par(mfrow=c(1,2))
hist(sigKY)
hist(sigLY)

#visualize female BLUPs on data scale
thresKY = stan_data_KY$nthres_clutch
mean_muKY = apply(samplesKY$mu_c_0, 2, median)
dKY = median(samplesKY$int_c)

thresLY = stan_data_LY$nthres_clutch
mean_muLY = apply(samplesLY$mu_c_0, 2, median)
dLY = median(samplesLY$int_c)

logistic = function(x){
  p = 1/(1 + exp(-x))
  p = ifelse(x == Inf, 1, p)
  p
}
predf = function(mu, sig, mu0, sig0, thres){
  p = matrix(NA, nrow = length(mu), ncol = thres + 1)
  for(j in 1:thres){
    if(j == 1){p[,j] = logistic( (mu0[j] - mu) / (exp(sig0 + sig)) ) }
    else p[,j] = logistic( (mu0[j] - mu) / (exp(sig0 + sig)) ) - 
        logistic( (mu0[j-1] - mu) / (exp(sig0 + sig)) )
  }
  p[,thres+1] = 1 - rowSums(p[,1:thres])
  return(p)
}
predKY = predf(mu = muKY, sig = sigKY, mu0 = mean_muKY, sig0 = dKY, thres = thresKY)
predLY = predf(mu = muLY, sig = sigLY, mu0 = mean_muLY, sig0 = dLY, thres = thresLY)

pKY = reshape2::melt(predKY)
pKY$site = "KY"
pKY$col = "col"
for(i in 1:max(pKY$Var1)){
  cond = pKY[pKY$Var1==i & pKY$Var2==5, "value"] >
          sum(pKY[pKY$Var1==i & pKY$Var2>5, "value"]) &
         pKY[pKY$Var1==i & pKY$Var2==5, "value"] >
          sum(pKY[pKY$Var1==i & pKY$Var2<5, "value"])
  condb = pKY[pKY$Var1==i & pKY$Var2==5, "value"] <
           sum(pKY[pKY$Var1==i & pKY$Var2>5, "value"])
  conds = pKY[pKY$Var1==i & pKY$Var2==5, "value"] <
    sum(pKY[pKY$Var1==i & pKY$Var2<5, "value"])
  
  pKY[pKY$Var1==i, "col"] =
    ifelse(cond == T, "b: p(median) > p(< median <)", 
           ifelse(condb == T,"c: p(median) < p(median <)",  
                  "a: p(median) < p(< median)"))
}

pLY = reshape2::melt(predLY)
pLY$site = "LY"
pLY$col = "col"
for(i in 1:max(pLY$Var1)){
  cond = pLY[pLY$Var1==i & pLY$Var2==4, "value"] >
    sum(pLY[pLY$Var1==i & pLY$Var2>4, "value"]) &
    pLY[pLY$Var1==i & pLY$Var2==4, "value"] >
    sum(pLY[pLY$Var1==i & pLY$Var2<4, "value"])
  condb = pLY[pLY$Var1==i & pLY$Var2==4, "value"] <
    sum(pLY[pLY$Var1==i & pLY$Var2>4, "value"])
  conds = pLY[pLY$Var1==i & pLY$Var2==4, "value"] <
    sum(pLY[pLY$Var1==i & pLY$Var2<4, "value"])
  
  pLY[pLY$Var1==i, "col"] =
    ifelse(cond == T, "b: p(median) > p(< median <)", 
         ifelse(condb == T,"c: p(median) < p(median <)",  
                "a: p(median) < p(< median)"))
  }

ppred = rbind(pKY,pLY)
ppred$Var2 = as.factor(ppred$Var2)

ggplot(ppred, aes(x = Var2, y = value, group = Var1))+
  geom_line(alpha = 0.3, col = "darkgrey")+
  xlab("Clutch size")+
  ylab("Probability")+
  facet_wrap(.~ site*col, nrow = 2, ncol = 3, drop = F)+
  theme(panel.border=element_rect(
      fill=NA,color="black", linewidth=1, linetype="solid"),
    strip.background = element_blank(),
    panel.background= element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  guides(col = "none")


###########################################################
#Selection gradients
###########################################################

#extract posteriors
samplesKY = readRDS("NLS_KY_m.RDS")
samplesLY = readRDS("NLS_LY_m.RDS")


########################################################
#Kentucky
########################################################
#selection gradients for non-Gaussian fitness measures

#initialize lists
beta_hatchKY = list()
gamma_hatchKY = list()
beta_surviveKY = list()
gamma_surviveKY = list()
beta_massKY = list()
gamma_massKY = list()

#environmental covariates
X_env = as.matrix(cbind(stan_data_KY$dateWSD, stan_data_KY$attemptSD,
                        stan_data_KY$dateWSD * stan_data_KY$attemptSD), ncol = 3)
#environmental effects
env_hatch = samplesKY$coef_h
env_survive = samplesKY$coef_s
env_mass = samplesKY$coef_m

#individual effects (unexplained selection)
id_hatch = samplesKY$I_h_0
id_survive = samplesKY$I_s_0
id_mass = samplesKY$I_m_0

#year effects
year_hatch = samplesKY$Y_h_0
year_survive = samplesKY$Y_s_0
year_mass = samplesKY$Y_m_0

#overdispersion effects
ovd_survive = samplesKY$O_s_0


#inverse link function
logistic = function(x){
  p = 1/(1 + exp(-x))
  p = ifelse(x == Inf, 1, p)
  p
}

#hatchling count function
for(i in 1:nrow(env_hatch) ){
  #Wbar = returns mean value of fitness function based on model likelihood
  #z[1]-[2] are placeholders for females' reaction norm intercepts and residuals
  Wbar = function(z,i) {
    eta =     ( #latent mean (zero-centered, no intercept)
              X_env %*% env_hatch[i,] + #environmental effects
              
              id_hatch[i,stan_data_KY$id] + #unexplained selection
              year_hatch[i,stan_data_KY$year_id] + #yearly effects  
              
              z[1] * samplesKY$b_h[i,1] + #linear selection coef
              z[2] * samplesKY$b_h[i,2] +
              
              z[1]^2 * samplesKY$q_h[i,1] + #nonlin selection coef
              z[2]^2 * samplesKY$q_h[i,2] +
              z[1] * z[2] * samplesKY$q_h[i,3] 
              )  

   #calculate probabilities for each hatchling count
   nthres = stan_data_KY$nthres_hatch
   d = exp(samplesKY$int_h[i])
   p = matrix(NA, nrow = length(eta), ncol = nthres+1)
   for(j in 1:nthres){
     if(j == 1){p[,j] = logistic( (samplesKY$mu_h_0[i,j] - eta) / d ) }
     else p[,j] = logistic( (samplesKY$mu_h_0[i,j] - eta) / d ) - 
                  logistic(  (samplesKY$mu_h_0[i,j-1] - eta) / d )
   }
   p[,nthres+1] = 1 - rowSums(p[,1:nthres])
   
   #predict expected hatchling counts
   eggcount = seq(1:(nthres+1)) #total number of eggs possible
   W = apply(p, 1, FUN = function(x) sum(x * eggcount) - 1) #-1 for original scale
   return(mean(W)) #mean fitness 
     }
  
  #calculate derivatives with respect to Wbar
  z = rep(0,2)
  first.derivatives = grad(func = Wbar, x = z, i = i)
  second.derivatives = hessian(func = Wbar, x = z, i = i)
  denom = Wbar(z = z, i = i)
  beta_hatchKY[[i]] = first.derivatives/denom
  gamma_hatchKY[[i]] = second.derivatives/denom
  }

#organize results
betas_hatchKY = do.call(rbind.data.frame, beta_hatchKY)  
gammas_hatchKY = do.call(rbind.data.frame, lapply(gamma_hatchKY, reshape2::melt) )

#survival function
for(i in 1:nrow(env_survive) ){
  Wbar = function(z,i) {
    p = logistic(
        samplesKY$mu_s_0[i] +
        X_env %*% env_survive[i,] + #environmental effects
        
        id_survive[i,stan_data_KY$id] + #unexplained selection
        year_survive[i,stan_data_KY$year_id] + #yearly effects  
        
        z[1] * samplesKY$b_s[i,1] + #linear selection coef
        z[2] * samplesKY$b_s[i,2] +
        
        z[1]^2 * samplesKY$q_s[i,1] + #nonlin selection coef
        z[2]^2 * samplesKY$q_s[i,2] +
        z[1] * z[2] * samplesKY$q_s[i,3] +
        ovd_survive[i,])
    return(mean(p))  }
  z = rep(0,2)
  first.derivatives = grad(func = Wbar, x = z, i = i)
  second.derivatives = hessian(func = Wbar, x = z, i = i)
  denom = Wbar(z = z, i = i)
  beta_surviveKY[[i]] = first.derivatives/denom
  gamma_surviveKY[[i]] = second.derivatives/denom
  }

#organize results
betas_surviveKY = do.call(rbind.data.frame, beta_surviveKY)  
gammas_surviveKY = do.call(rbind.data.frame, lapply(gamma_surviveKY, reshape2::melt) )

#mass function
for(i in 1:nrow(env_survive) ){
  Wbar = function(z,i) {
    p = 
        samplesKY$mu_m_0[i] +
        X_env %*% env_mass[i,] + #environmental effects
        
        id_mass[i,stan_data_KY$id] + #unexplained selection
        year_mass[i,stan_data_KY$year_id] + #yearly effects  
        
        z[1] * samplesKY$b_m[i,1] + #linear selection coef
        z[2] * samplesKY$b_m[i,2] +
        
        z[1]^2 * samplesKY$q_m[i,1] + #nonlin selection coef
        z[2]^2 * samplesKY$q_m[i,2] +
        z[1] * z[2] * samplesKY$q_m[i,3]
    return(mean(p))  }
  z = rep(0,2)
  first.derivatives = grad(func = Wbar, x = z, i = i)
  second.derivatives = hessian(func = Wbar, x = z, i = i)
  denom = Wbar(z = z, i = i)
  beta_massKY[[i]] = first.derivatives/denom
  gamma_massKY[[i]] = second.derivatives/denom
}

#mass selection gradients
betas_massKY = samplesKY$b_m
gammas_massKY = data.frame(cbind(samplesKY$q_m[,1]*2, 
                                 samplesKY$q_m[,2]*2,
                                 samplesKY$q_m[,3]))

#organize results
betas_massKY = do.call(rbind.data.frame, beta_massKY)  
gammas_massKY = do.call(rbind.data.frame, lapply(gamma_massKY, reshape2::melt) )

#add correct names
colnames(betas_hatchKY) = colnames(betas_surviveKY) = colnames(betas_massKY) = c("m","v")

q_order = c("m^2","v^2","m*v")

gammas_hatchKY$par = ifelse(gammas_hatchKY$Var1==1 & gammas_hatchKY$Var2==1, "m^2",
                            ifelse(gammas_hatchKY$Var1==2 & gammas_hatchKY$Var2==2, "v^2",
                             ifelse(gammas_hatchKY$Var1==1 & gammas_hatchKY$Var2==2, "m*v",NA)))
gammas_hatchKY = na.omit(gammas_hatchKY) #remove redundant matrix elements
gammas_hatchKY = data.frame(data.table::dcast(setDT(gammas_hatchKY),
                                              rowid(par) ~ factor(gammas_hatchKY$par, levels = q_order ) ))[-1]

gammas_surviveKY$par = ifelse(gammas_surviveKY$Var1==1 & gammas_surviveKY$Var2==1, "m^2",
                            ifelse(gammas_surviveKY$Var1==2 & gammas_surviveKY$Var2==2, "v^2",
                                   ifelse(gammas_surviveKY$Var1==1 & gammas_surviveKY$Var2==2, "m*v",NA)))
gammas_surviveKY = na.omit(gammas_surviveKY) #remove redundant matrix elements
gammas_surviveKY = data.frame(data.table::dcast(setDT(gammas_surviveKY),
                                              rowid(par) ~ factor(gammas_surviveKY$par, levels = q_order ) ))[-1]

gammas_massKY$par = ifelse(gammas_massKY$Var1==1 & gammas_massKY$Var2==1, "m^2",
                              ifelse(gammas_massKY$Var1==2 & gammas_massKY$Var2==2, "v^2",
                                     ifelse(gammas_massKY$Var1==1 & gammas_massKY$Var2==2, "m*v",NA)))
gammas_massKY = na.omit(gammas_massKY) #remove redundant matrix elements
gammas_massKY = data.frame(data.table::dcast(setDT(gammas_massKY),
                                                rowid(par) ~ factor(gammas_massKY$par, levels = q_order ) ))[-1]


#convert selection gradients to standardized gradients
sd_selKY = data.frame(sd_m = apply(samplesKY$I_c[,,1],1, FUN=sd),
                      sd_v = apply(samplesKY$I_c[,,2],1, FUN=sd),
                      sd_m.m = apply(samplesKY$I_c[,,1]^2,1, FUN=sd),
                      sd_v.v = apply(samplesKY$I_c[,,2]^2,1, FUN=sd),
                      sd_m.v = apply(samplesKY$I_c[,,1]*samplesKY$I_c[,,2],1, FUN=sd))

#standardize
betas_sd_hatchKY = betas_hatchKY * sd_selKY[,1:2]
betas_sd_surviveKY = betas_surviveKY * sd_selKY[,1:2]
betas_sd_massKY = betas_massKY  * sd_selKY[,1:2] 
gammas_sd_hatchKY = gammas_hatchKY * sd_selKY[,3:5]
gammas_sd_surviveKY = gammas_surviveKY * sd_selKY[,3:5]
gammas_sd_massKY = gammas_massKY * sd_selKY[,3:5]

#total selection gradients

#linear selection
betas_totalKY = betas_hatchKY + betas_surviveKY + betas_massKY
betas_sd_totalKY = betas_totalKY * sd_selKY[,1:2]

#nonlinear selection

#make Gamma matrices
Gam_h = list() #initialize list of matrices
for(i in 1:nrow(betas_totalKY)) {
  #diagonal with quad gradients
  gamma = diag(c(gammas_hatchKY[i,1], gammas_hatchKY[i,2]))
  #add in off-diagonal elements
  gamma[1,2] = gammas_hatchKY$m.v[i]
  #make symmetric
  gamma[lower.tri(gamma)] = t(gamma)[lower.tri(gamma)]
  #add to list
  Gam_h[[i]] = gamma
}

Gam_s = list()
for(i in 1:nrow(betas_totalKY)) {
  gamma = diag(c(gammas_surviveKY[i,1], gammas_surviveKY[i,2]))
  gamma[1,2] = gammas_surviveKY$m.v[i]
  gamma[lower.tri(gamma)] = t(gamma)[lower.tri(gamma)]
  Gam_s[[i]] = gamma
}

Gam_m = list()
for(i in 1:nrow(betas_totalKY)) {
  gamma = diag(c(gammas_massKY[i,1], gammas_massKY[i,2]))
  gamma[1,2] = gammas_massKY$m.v[i]
  gamma[lower.tri(gamma)] = t(gamma)[lower.tri(gamma)]
  Gam_m[[i]] = gamma
}

gammas_totalKY = gammas_hatchKY #initialize df of correct size
for(i in 1:nrow(gammas_hatchKY)){
  Gam_tot = Gam_h[[i]] + Gam_s[[i]] + Gam_m[[i]] +
    t(betas_totalKY[i,]) %*% t(t(betas_totalKY[i,]))
  gammas_totalKY[i,] =
    c(Gam_tot[1,1],Gam_tot[2,2],Gam_tot[lower.tri(Gam_tot)][1])
}
gammas_sd_totalKY = gammas_totalKY * sd_selKY[,3:5]

#summarize
selection_KY =
  data.frame(par = c(" ", paste0("",colnames(betas_sd_hatchKY)),
                     paste0("",colnames(gammas_sd_massKY))),
             h_KY = c("median (90% CI)", 
                      paste0(round(apply(betas_sd_hatchKY,2,median),2)," (",
                             round(apply(betas_sd_hatchKY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(betas_sd_hatchKY,2,quantile,c(0.05,0.95))[2,],2),
                             ")"),
                      paste0(round(apply(gammas_sd_hatchKY,2,median),2)," (",
                             round(apply(gammas_sd_hatchKY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(gammas_sd_hatchKY,2,quantile,c(0.05,0.95))[2,],2),
                             ")")),
             s_KY = c("median (90% CI)", 
                      paste0(round(apply(betas_sd_surviveKY,2,median),2)," (",
                             round(apply(betas_sd_surviveKY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(betas_sd_surviveKY,2,quantile,c(0.05,0.95))[2,],2),
                             ")"),
                      paste0(round(apply(gammas_sd_surviveKY,2,median),2)," (",
                             round(apply(gammas_sd_surviveKY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(gammas_sd_surviveKY,2,quantile,c(0.05,0.95))[2,],2),
                             ")")),
             m_KY = c("median (90% CI)", 
                      paste0(round(apply(betas_sd_massKY,2,median),2)," (",
                             round(apply(betas_sd_massKY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(betas_sd_massKY,2,quantile,c(0.05,0.95))[2,],2),
                             ")"),
                      paste0(round(apply(gammas_sd_massKY,2,median),2)," (",
                             round(apply(gammas_sd_massKY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(gammas_sd_massKY,2,quantile,c(0.05,0.95))[2,],2),
                             ")")),
             tot_KY = c("median (90% CI)", 
                        paste0(round(apply(betas_sd_totalKY,2,median),2)," (",
                               round(apply(betas_sd_totalKY,2,quantile,c(0.05,0.95))[1,],2),
                               ",", round(apply(betas_sd_totalKY,2,quantile,c(0.05,0.95))[2,],2),
                               ")"),
                        paste0(round(apply(gammas_sd_totalKY,2,median),2)," (",
                               round(apply(gammas_sd_totalKY,2,quantile,c(0.05,0.95))[1,],2),
                               ",", round(apply(gammas_sd_totalKY,2,quantile,c(0.05,0.95))[2,],2),
                               ")")))

#save
write.csv(selection_KY,"selection_KY.csv")

#posterior probs for component-specific selection
apply(betas_sd_hatchKY, 2, FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(gammas_sd_hatchKY, 2,  FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(betas_sd_surviveKY, 2,  FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(gammas_sd_surviveKY, 2,  FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(betas_sd_massKY, 2,  FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(gammas_sd_massKY, 2,  FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))

#summarize total gradients
apply(betas_sd_totalKY, 2, median)
apply(betas_sd_totalKY, 2, quantile, c(0.05,0.95))
apply(gammas_sd_totalKY, 2, median)
apply(gammas_sd_totalKY, 2, quantile, c(0.05,0.95))

########################################################
#Lundy
########################################################
#selection gradients for non-Gaussian fitness measures

#initialize lists
beta_hatchLY = list()
gamma_hatchLY = list()
beta_surviveLY = list()
gamma_surviveLY = list()
beta_massLY = list()
gamma_massLY = list()

#environmental covariates
X_env = as.matrix(cbind(stan_data_LY$dateWSD, stan_data_LY$attemptSD,
                        stan_data_LY$dateWSD * stan_data_LY$attemptSD), ncol = 3)
#environmental effects
env_hatch = samplesLY$coef_h
env_survive = samplesLY$coef_s
env_mass = samplesLY$coef_m

#individual effects (unexplained selection)
id_hatch = samplesLY$I_h_0
id_survive = samplesLY$I_s_0
id_mass = samplesLY$I_m_0

#year effects
year_hatch = samplesLY$Y_h_0
year_survive = samplesLY$Y_s_0
year_mass = samplesLY$Y_m_0

#overdispersion effects
ovd_survive = samplesLY$O_s_0


#inverse link function
logistic = function(x){
  p = 1/(1 + exp(-x))
  p = ifelse(x == Inf, 1, p)
  p
}

#hatchling count function
for(i in 1:nrow(env_hatch) ){
  #Wbar = returns mean value of fitness function based on model likelihood
  #z[1]-[2] are placeholders for females' reaction norm intercepts and residuals
  Wbar = function(z,i) {
    eta =     ( #latent mean (zero-centered, no intercept)
      X_env %*% env_hatch[i,] + #environmental effects
        
        id_hatch[i,stan_data_LY$id] + #unexplained selection
        year_hatch[i,stan_data_LY$year_id] + #yearly effects  
        
        z[1] * samplesLY$b_h[i,1] + #linear selection coef
        z[2] * samplesLY$b_h[i,2] +
        
        z[1]^2 * samplesLY$q_h[i,1] + #nonlin selection coef
        z[2]^2 * samplesLY$q_h[i,2] +
        z[1] * z[2] * samplesLY$q_h[i,3] 
    )  
    
    #calculate probabilities for each hatchling count
    nthres = stan_data_LY$nthres_hatch
    d = exp(samplesLY$int_h[i])
    p = matrix(NA, nrow = length(eta), ncol = nthres+1)
    for(j in 1:nthres){
      if(j == 1){p[,j] = logistic( (samplesLY$mu_h_0[i,j] - eta) / d ) }
      else p[,j] = logistic( (samplesLY$mu_h_0[i,j] - eta) / d ) - 
          logistic( (samplesLY$mu_h_0[i,j-1] - eta) / d )
    }
    p[,nthres+1] = 1 - rowSums(p[,1:nthres])
    
    #predict expected hatchling counts
    eggcount = seq(1:(nthres+1)) #total number of eggs possible
    W = apply(p, 1, FUN = function(x) sum(x * eggcount) - 1) #-1 for original scale
    return(mean(W)) #mean fitness 
  }
  
  #calculate derivatives with respect to Wbar
  z = rep(0,2)
  first.derivatives = grad(func = Wbar, x = z, i = i)
  second.derivatives = hessian(func = Wbar, x = z, i = i)
  denom = Wbar(z = z, i = i)
  beta_hatchLY[[i]] = first.derivatives/denom
  gamma_hatchLY[[i]] = second.derivatives/denom
}

#organize results
betas_hatchLY = do.call(rbind.data.frame, beta_hatchLY)  
gammas_hatchLY = do.call(rbind.data.frame, lapply(gamma_hatchLY, reshape2::melt) )

#survival function
for(i in 1:nrow(env_survive) ){
  Wbar = function(z,i) {
    p = logistic(
      samplesLY$mu_s_0[i] +
        X_env %*% env_survive[i,] + #environmental effects
        
        id_survive[i,stan_data_LY$id] + #unexplained selection
        year_survive[i,stan_data_LY$year_id] + #yearly effects  
        
        z[1] * samplesLY$b_s[i,1] + #linear selection coef
        z[2] * samplesLY$b_s[i,2] +
        
        z[1]^2 * samplesLY$q_s[i,1] + #nonlin selection coef
        z[2]^2 * samplesLY$q_s[i,2] +
        z[1] * z[2] * samplesLY$q_s[i,3] +
        ovd_survive[i,])
    return(mean(p))  }
  z = rep(0,2)
  first.derivatives = grad(func = Wbar, x = z, i = i)
  second.derivatives = hessian(func = Wbar, x = z, i = i)
  denom = Wbar(z = z, i = i)
  beta_surviveLY[[i]] = first.derivatives/denom
  gamma_surviveLY[[i]] = second.derivatives/denom
}

#organize results
betas_surviveLY = do.call(rbind.data.frame, beta_surviveLY)  
gammas_surviveLY = do.call(rbind.data.frame, lapply(gamma_surviveLY, reshape2::melt) )

#mass function
for(i in 1:nrow(env_survive) ){
  Wbar = function(z,i) {
    p = 
      samplesLY$mu_m_0[i] +
      X_env %*% env_mass[i,] + #environmental effects
      
      id_mass[i,stan_data_LY$id] + #unexplained selection
      year_mass[i,stan_data_LY$year_id] + #yearly effects  
      
      z[1] * samplesLY$b_m[i,1] + #linear selection coef
      z[2] * samplesLY$b_m[i,2] +
      
      z[1]^2 * samplesLY$q_m[i,1] + #nonlin selection coef
      z[2]^2 * samplesLY$q_m[i,2] +
      z[1] * z[2] * samplesLY$q_m[i,3]
    return(mean(p))  }
  z = rep(0,2)
  first.derivatives = grad(func = Wbar, x = z, i = i)
  second.derivatives = hessian(func = Wbar, x = z, i = i)
  denom = Wbar(z = z, i = i)
  beta_massLY[[i]] = first.derivatives/denom
  gamma_massLY[[i]] = second.derivatives/denom
}

#organize results
betas_massLY = do.call(rbind.data.frame, beta_massLY)  
gammas_massLY = do.call(rbind.data.frame, lapply(gamma_massLY, reshape2::melt) )

#add correct names
colnames(betas_hatchLY) = colnames(betas_surviveLY) = colnames(betas_massLY) = c("m","v")

q_order = c("m^2","v^2","m*v")

gammas_hatchLY$par = ifelse(gammas_hatchLY$Var1==1 & gammas_hatchLY$Var2==1, "m^2",
                            ifelse(gammas_hatchLY$Var1==2 & gammas_hatchLY$Var2==2, "v^2",
                                   ifelse(gammas_hatchLY$Var1==1 & gammas_hatchLY$Var2==2, "m*v",NA)))
gammas_hatchLY = na.omit(gammas_hatchLY) #remove redundant matrix elements
gammas_hatchLY = data.frame(data.table::dcast(setDT(gammas_hatchLY),
                                              rowid(par) ~ factor(gammas_hatchLY$par, levels = q_order ) ))[-1]

gammas_surviveLY$par = ifelse(gammas_surviveLY$Var1==1 & gammas_surviveLY$Var2==1, "m^2",
                              ifelse(gammas_surviveLY$Var1==2 & gammas_surviveLY$Var2==2, "v^2",
                                     ifelse(gammas_surviveLY$Var1==1 & gammas_surviveLY$Var2==2, "m*v",NA)))
gammas_surviveLY = na.omit(gammas_surviveLY) #remove redundant matrix elements
gammas_surviveLY = data.frame(data.table::dcast(setDT(gammas_surviveLY),
                                                rowid(par) ~ factor(gammas_surviveLY$par, levels = q_order ) ))[-1]

gammas_massLY$par = ifelse(gammas_massLY$Var1==1 & gammas_massLY$Var2==1, "m^2",
                           ifelse(gammas_massLY$Var1==2 & gammas_massLY$Var2==2, "v^2",
                                  ifelse(gammas_massLY$Var1==1 & gammas_massLY$Var2==2, "m*v",NA)))
gammas_massLY = na.omit(gammas_massLY) #remove redundant matrix elements
gammas_massLY = data.frame(data.table::dcast(setDT(gammas_massLY),
                                             rowid(par) ~ factor(gammas_massLY$par, levels = q_order ) ))[-1]


#convert selection gradients to standardized gradients
sd_selLY = data.frame(sd_m = apply(samplesLY$I_c[,,1],1, FUN=sd),
                      sd_v = apply(samplesLY$I_c[,,2],1, FUN=sd),
                      sd_m.m = apply(samplesLY$I_c[,,1]^2,1, FUN=sd),
                      sd_v.v = apply(samplesLY$I_c[,,2]^2,1, FUN=sd),
                      sd_m.v = apply(samplesLY$I_c[,,1]*samplesLY$I_c[,,2],1, FUN=sd))

#standardize
betas_sd_hatchLY = betas_hatchLY * sd_selLY[,1:2]
betas_sd_surviveLY = betas_surviveLY * sd_selLY[,1:2]
betas_sd_massLY = betas_massLY  * sd_selLY[,1:2] 
gammas_sd_hatchLY = gammas_hatchLY * sd_selLY[,3:5]
gammas_sd_surviveLY = gammas_surviveLY * sd_selLY[,3:5]
gammas_sd_massLY = gammas_massLY * sd_selLY[,3:5]

#total selection gradients

#linear selection
betas_totalLY = betas_hatchLY + betas_surviveLY + betas_massLY
betas_sd_totalLY = betas_totalLY * sd_selLY[,1:2]

#nonlinear selection

#make Gamma matrices
Gam_h = list() #initialize list of matrices
for(i in 1:nrow(betas_totalLY)) {
  #diagonal with quad gradients
  gamma = diag(c(gammas_hatchLY[i,1], gammas_hatchLY[i,2]))
  #add in off-diagonal elements
  gamma[1,2] = gammas_hatchLY$m.v[i]
  #make symmetric
  gamma[lower.tri(gamma)] = t(gamma)[lower.tri(gamma)]
  #add to list
  Gam_h[[i]] = gamma
}

Gam_s = list()
for(i in 1:nrow(betas_totalLY)) {
  gamma = diag(c(gammas_surviveLY[i,1], gammas_surviveLY[i,2]))
  gamma[1,2] = gammas_surviveLY$m.v[i]
  gamma[lower.tri(gamma)] = t(gamma)[lower.tri(gamma)]
  Gam_s[[i]] = gamma
}

Gam_m = list()
for(i in 1:nrow(betas_totalLY)) {
  gamma = diag(c(gammas_massLY[i,1], gammas_massLY[i,2]))
  gamma[1,2] = gammas_massLY$m.v[i]
  gamma[lower.tri(gamma)] = t(gamma)[lower.tri(gamma)]
  Gam_m[[i]] = gamma
}

gammas_totalLY = gammas_hatchLY #initialize df of correct size
for(i in 1:nrow(gammas_hatchLY)){
  Gam_tot = Gam_h[[i]] + Gam_s[[i]] + Gam_m[[i]] +
    t(betas_totalLY[i,]) %*% t(t(betas_totalLY[i,]))
  gammas_totalLY[i,] =
    c(Gam_tot[1,1],Gam_tot[2,2],Gam_tot[lower.tri(Gam_tot)][1])
}
gammas_sd_totalLY = gammas_totalLY * sd_selLY[,3:5]

#summarize
selection_LY =
  data.frame(par = c(" ", paste0("",colnames(betas_sd_hatchLY)),
                     paste0("",colnames(gammas_sd_massLY))),
             h_LY = c("median (90% CI)", 
                      paste0(round(apply(betas_sd_hatchLY,2,median),2)," (",
                             round(apply(betas_sd_hatchLY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(betas_sd_hatchLY,2,quantile,c(0.05,0.95))[2,],2),
                             ")"),
                      paste0(round(apply(gammas_sd_hatchLY,2,median),2)," (",
                             round(apply(gammas_sd_hatchLY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(gammas_sd_hatchLY,2,quantile,c(0.05,0.95))[2,],2),
                             ")")),
             s_LY = c("median (90% CI)", 
                      paste0(round(apply(betas_sd_surviveLY,2,median),2)," (",
                             round(apply(betas_sd_surviveLY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(betas_sd_surviveLY,2,quantile,c(0.05,0.95))[2,],2),
                             ")"),
                      paste0(round(apply(gammas_sd_surviveLY,2,median),2)," (",
                             round(apply(gammas_sd_surviveLY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(gammas_sd_surviveLY,2,quantile,c(0.05,0.95))[2,],2),
                             ")")),
             m_LY = c("median (90% CI)", 
                      paste0(round(apply(betas_sd_massLY,2,median),2)," (",
                             round(apply(betas_sd_massLY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(betas_sd_massLY,2,quantile,c(0.05,0.95))[2,],2),
                             ")"),
                      paste0(round(apply(gammas_sd_massLY,2,median),2)," (",
                             round(apply(gammas_sd_massLY,2,quantile,c(0.05,0.95))[1,],2),
                             ",", round(apply(gammas_sd_massLY,2,quantile,c(0.05,0.95))[2,],2),
                             ")")),
             tot_LY = c("median (90% CI)", 
                        paste0(round(apply(betas_sd_totalLY,2,median),2)," (",
                               round(apply(betas_sd_totalLY,2,quantile,c(0.05,0.95))[1,],2),
                               ",", round(apply(betas_sd_totalLY,2,quantile,c(0.05,0.95))[2,],2),
                               ")"),
                        paste0(round(apply(gammas_sd_totalLY,2,median),2)," (",
                               round(apply(gammas_sd_totalLY,2,quantile,c(0.05,0.95))[1,],2),
                               ",", round(apply(gammas_sd_totalLY,2,quantile,c(0.05,0.95))[2,],2),
                               ")")))

#save
write.csv(selection_LY,"selection_LY.csv")

#posterior probs for component-specific selection
apply(betas_sd_hatchLY, 2, FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(gammas_sd_hatchLY, 2, FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(betas_sd_surviveLY, 2, FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(gammas_sd_surviveLY, 2, FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(betas_sd_massLY, 2, FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))
apply(gammas_sd_massLY, 2, FUN = function(x) sum(sign(x) == sign(median(x)))/length(x))

#summarize total gradients
apply(betas_sd_totalLY, 2, median)
apply(betas_sd_totalLY, 2, quantile, c(0.05,0.95))
apply(gammas_sd_totalLY, 2, median)
apply(gammas_sd_totalLY, 2, quantile, c(0.05,0.95))


##############################################################
#plot gradients
##############################################################

#organize in long format
ldf_hatchKY = rbind(reshape2::melt(betas_sd_hatchKY),
                    reshape2::melt(gammas_sd_hatchKY))
ldf_hatchKY$fit = "hatch"
ldf_surviveKY = rbind(reshape2::melt(betas_sd_surviveKY),
                      reshape2::melt(gammas_sd_surviveKY))
ldf_surviveKY$fit = "survive"
ldf_massKY = rbind(reshape2::melt(betas_sd_massKY),
                   reshape2::melt(gammas_sd_massKY))
ldf_massKY$fit = "mass"
ldf_totKY = rbind(reshape2::melt(betas_sd_totalKY),
                  reshape2::melt(gammas_sd_totalKY))
ldf_totKY$fit = "total"
ldf_KY = rbind(ldf_hatchKY, ldf_surviveKY, ldf_massKY, ldf_totKY)
ldf_KY$pop = "KY"

ldf_hatchLY = rbind(reshape2::melt(betas_sd_hatchLY),
                    reshape2::melt(gammas_sd_hatchLY))
ldf_hatchLY$fit = "hatch"
ldf_surviveLY = rbind(reshape2::melt(betas_sd_surviveLY),
                      reshape2::melt(gammas_sd_surviveLY))
ldf_surviveLY$fit = "survive"
ldf_massLY = rbind(reshape2::melt(betas_sd_massLY),
                   reshape2::melt(gammas_sd_massLY))
ldf_massLY$fit = "mass"
ldf_totLY = rbind(reshape2::melt(betas_sd_totalLY),
                  reshape2::melt(gammas_sd_totalLY))
ldf_totLY$fit = "total"
ldf_LY = rbind(ldf_hatchLY, ldf_surviveLY, ldf_massLY, ldf_totLY)
ldf_LY$pop = "LY"

#combine
ldf = rbind(ldf_KY, ldf_LY)

#reorder
ldf$fit2 = factor(ldf$fit, levels = c("hatch","survive","mass","total"))
ldf$variable2 = factor(ldf$variable, levels = 
                         c("m","v","m.2","v.2", "m.v"))

#plot
NLS_sdplot =
  ggplot(ldf, aes(x = value, y = variable2, color = pop, group = interaction(fit2,pop)) )+
  stat_pointinterval(.width=c(0.9), size = 5, position = position_dodge(-0.5))+
  facet_wrap(. ~ fit2, nrow = 1)+
  coord_cartesian(xlim=c(-0.2,0.2))+
  scale_y_discrete(limits = rev(levels(ldf$variable2)),
                   labels = c('m' = parse(text=TeX('$\\beta_{mu_0}$')),
                              'v' = parse(text=TeX('$\\beta_{sigma_0}$')),
                              'm.2' = parse(text=TeX('$\\gamma_{mu_0}$')),
                              'v.2' = parse(text=TeX('$\\gamma_{sigma_0}$')),
                              'm.v' = parse(text=TeX('$\\gamma_{mu_0 \\cdot \\sigma_0}$'))))+
  geom_vline(xintercept= 0, linetype = "dashed", size = 0.75)+
  xlab("\nStandardized selection gradients on clutch sizes (median +/- 90% CI)")+
  labs(color = "Population")+
  scale_fill_manual(values=c("deeppink1","purple"))+
  scale_color_manual(values=c("deeppink1","purple"))+
  theme(legend.title = element_text(face = "bold"),
        legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.margin = margin(0.1,0.1,0.1,0.3,unit = "in"))

ggsave("NLS_sdplot_tot.png", NLS_sdplot, height = 3, width = 8)

#############################################################
#Plot selection surfaces
#############################################################

#combo plot
library(plot3D)
png("surface_KY.png", units = "in", width = 8, height = 5, res = 600)
par(mar = c(1.2,1.2,1.2,1.2), mfrow = c(2,3))

#KY
{
  sd_m = median(samplesKY$sd_I_c[,1])
  sd_v = median(samplesKY$sd_I_c[,2])
  m_blup = apply(samplesKY$I_c[,,1], 2, median)
  m_blup = m_blup / sd_m
  v_blup = apply(samplesKY$I_c[,,2], 2, median)
  v_blup = v_blup / sd_v
  
  hb = apply(betas_sd_hatchKY,2,median)
  hq = apply(gammas_sd_hatchKY,2,median)
  Wh = function(m,v) {
    eta =     ( #latent mean (zero-centered, no intercept)
      mean(stan_data_KY$hatch - 1)+
        m * hb[1] + #linear selection coef
        v *  hb[2] + 
        m^2 * hq[1] + #nonlin selection coef
        v^2 * hq[2] +
        m * v * hq[3] 
    )  
    W = eta
    return(W) #mean fitness 
  }
  w_blup = Wh(m_blup,v_blup)
  m = seq(min(m_blup) - 1.5*sd_m , max(m_blup) + 1.5*sd_m, by = 0.05)
  v = seq(min(v_blup) - 1.5*sd_v, max(v_blup) + 1.5*sd_v, by = 0.1)
  w = outer(X = m, Y = v, FUN = Wh)
  pop = seq(min(w_blup),max(w_blup),0.05)
  
  colr = colorRampPalette(c("white","deeppink1"))
  persp3D(x = m, y = v, z = w, lwd = 0.5, border = NA, bty = "u", 
          col = colr(100), xlab = "mu", ylab = "sigma", zlab = "\nhatchling success",
          label.rotation = list(z = 90), phi = 30, theta = 500,
          colkey = FALSE, depth_mask = FALSE, contour = list(nlevels = 50, drawlabels = TRUE))
  
  hb = apply(betas_sd_massKY,2,median)
  hq = apply(gammas_sd_massKY,2,median)
  Wh = function(m,v) {
    eta =     ( #latent mean (zero-centered, no intercept)
      mean(stan_data_KY$mass)+
        m * hb[1] + #linear selection coef
        v * hb[2] + 
        m^2 * hq[1] + #nonlin selection coef
        v^2 * hq[2] +
        m * v * hq[3]
    )  
    W = eta
    return(W) #mean fitness 
  }
  w_blup = Wh(m_blup,v_blup)
  m = seq(min(m_blup) - 1.5*sd_m , max(m_blup) + 1.5*sd_m, by = 0.05)
  v = seq(min(v_blup) - 1.5*sd_v, max(v_blup) + 1.5*sd_v, by = 0.1)
  w = outer(X = m, Y = v, FUN = Wh)
  pop = seq(min(w_blup),max(w_blup),0.05)
  
  colr = colorRampPalette(c("white","deeppink1"))
  persp3D(x = m, y = v, z = w, lwd = 0.5, border = NA, bty = "u", 
          col = colr(100), xlab = "mu", ylab = "sigma", zlab = "\nnestling mass",
          label.rotation = list(z = 90), phi = 30, theta = 500,
          colkey = FALSE, depth_mask = FALSE, contour = list(nlevels = 50, drawlabels = TRUE))
  
  hb = apply(betas_sd_totalKY,2,median)
  hq = apply(gammas_sd_totalKY,2,median)
  Wh = function(m,v) {
    eta =     ( #latent mean (zero-centered, no intercept)
      1 +
        m * hb[1] + #linear selection coef
        v * -1 * hb[2] + 
        m^2 * hq[1] + #nonlin selection coef
        v^2 * hq[2] +
        m * v * hq[3]
    )  
    W = eta
    return(W) #mean fitness 
  }
  w_blup = Wh(m_blup,v_blup)
  m = seq(min(m_blup) - 1.5*sd_m , max(m_blup) + 1.5*sd_m, by = 0.05)
  v = seq(min(v_blup) - 1.5*sd_v, max(v_blup) + 1.5*sd_v, by = 0.1)
  w = outer(X = m, Y = v, FUN = Wh)
  pop = seq(min(w_blup),max(w_blup),0.05)
  
  colr = colorRampPalette(c("white","deeppink1"))
  persp3D(x = m, y = v, z = w, lwd = 0.5, border = NA, bty = "u", 
          col = colr(100), xlab = "mu", ylab = "sigma", zlab = "\nfitness",
          label.rotation = list(z = 90), phi = 30, theta = 500,
          colkey = FALSE, depth_mask = FALSE, contour = list(nlevels = 50, drawlabels = TRUE))
}

#LY
{
  sd_m = median(samplesLY$sd_I_c[,1])
  sd_v = median(samplesLY$sd_I_c[,2])
  m_blup = apply(samplesLY$I_c[,,1], 2, median)
  m_blup = m_blup / sd_m
  v_blup = apply(samplesLY$I_c[,,2], 2, median)
  v_blup = v_blup / sd_v
  
  hb = apply(betas_sd_hatchLY,2,median)
  hq = apply(gammas_sd_hatchLY,2,median)
  Wh = function(m,v) {
    eta =     ( #latent mean (zero-centered, no intercept)
      mean(stan_data_LY$hatch - 1)+
        m * hb[1] + #linear selection coef
        v * hb[2] + 
        m^2 * hq[1] + #nonlin selection coef
        v^2 * hq[2] +
        m * v * hq[3]
    )  
    W = eta
    return(W) #mean fitness 
  }
  w_blup = Wh(m_blup,v_blup)
  m = seq(min(m_blup) - 1.5*sd_m , max(m_blup) + 1.5*sd_m, by = 0.05)
  v = seq(min(v_blup) - 1.5*sd_v, max(v_blup) + 1.5*sd_v, by = 0.1)
  w = outer(X = m, Y = v, FUN = Wh)
  pop = seq(min(w_blup),max(w_blup),0.05)
  
  colr = colorRampPalette(c("white","purple"))
  persp3D(x = m, y = v, z = w, lwd = 0.5, border = NA, bty = "u", 
          col = colr(100), xlab = "mu", ylab = "sigma", zlab = "\nhatchling success",
          label.rotation = list(z = 90), phi = 30, theta = 500,
          colkey = FALSE, depth_mask = FALSE, contour = list(nlevels = 50, drawlabels = TRUE))
  
  hb = apply(betas_sd_massLY,2,median)
  hq = apply(gammas_sd_massLY,2,median)
  Wh = function(m,v) {
    eta =     ( #latent mean (zero-centered, no intercept)
      mean(stan_data_LY$mass)+
        m * hb[1] + #linear selection coef
        v * hb[2] + 
        m^2 * hq[1] + #nonlin selection coef
        v^2 * hq[2] +
        m * v * hq[3] 
    )  
    W = eta
    return(W) #mean fitness 
  }
  w_blup = Wh(m_blup,v_blup)
  m = seq(min(m_blup) - 1.5*sd_m , max(m_blup) + 1.5*sd_m, by = 0.05)
  v = seq(min(v_blup) - 1.5*sd_v, max(v_blup) + 1.5*sd_v, by = 0.1)
  w = outer(X = m, Y = v, FUN = Wh)
  pop = seq(min(w_blup),max(w_blup),0.05)
  
  colr = colorRampPalette(c("white","purple"))
  persp3D(x = m, y = v, z = w, lwd = 0.5, border = NA, bty = "u", 
          col = colr(100), xlab = "mu", ylab = "sigma", zlab = "\nnestling mass",
          label.rotation = list(z = 90), phi = 30, theta = 500,
          colkey = FALSE, depth_mask = FALSE, contour = list(nlevels = 50, drawlabels = TRUE))

  hb = apply(betas_sd_totalLY,2,median)
  hq = apply(gammas_sd_totalLY,2,median)
  Wh = function(m,v) {
    eta =     ( #latent mean (zero-centered, no intercept)
      1 +
        m * hb[1] + #linear selection coef
        v * hb[2] +
        m^2 * hq[1] + #nonlin selection coef
        v^2 * hq[2] +
        m * v *  hq[3]
    )  
    W = eta
    return(W) #mean fitness 
  }
  w_blup = Wh(m_blup,v_blup)
  m = seq(min(m_blup) - 1.5*sd_m , max(m_blup) + 1.5*sd_m, by = 0.05)
  v = seq(min(v_blup) - 1.5*sd_v, max(v_blup) + 1.5*sd_v, by = 0.1)
  w = outer(X = m, Y = v, FUN = Wh)
  pop = seq(min(w_blup),max(w_blup),0.05)
  
  colr = colorRampPalette(c("white","purple"))
  persp3D(x = m, y = v, z = w, lwd = 0.5, border = NA, bty = "u", 
          col = colr(100), xlab = "mu", ylab = "sigma", zlab = "\nfitness",
          label.rotation = list(z = 90), phi = 30, theta = 500,
          colkey = FALSE, depth_mask = FALSE, contour = list(nlevels = 50, drawlabels = TRUE))
}

dev.off()







