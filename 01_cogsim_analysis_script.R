#This script replicates the results from Guilbeault et al.'s (2023) paper in Management Science entitled: 
#Exposure to the Views of Opposing Others with Latent Cognitive Differences Results in Social Influence---
#But Only When Those Differences Remain Obscured
#Please direct correspondence to douglas.guilbeault@haas.berkeley.edu

###############
#Load packages#
###############
rm(list=ls());gc()
library(dplyr);library(ggplot2);library(tidyverse);library(tidyr)
library(aod);library(lme4);library(nlme);library(glmmTMB);library(sjPlot) 
library(clinfun);library(foreign);library(multcomp);library(sjPlot)
library(sjmisc);library(stats);library(lmtest);library(stringr);library(miceadds)
library(jtools);library(performance);library(lsa);library(mediation);library(multiwayvcov)

###########
#Functions#
###########

min_max_norm<-function(x){(x - min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T))}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

###############
#Load/org data#
###############
datapath<-"" #path to folder where your data is saved 
cogsim<-read.csv(paste(datapath, "exp_dt.csv", sep=""))
cogsim_long <- gather(cogsim, measure, measurement,  round1:round2_dissim, factor_key=TRUE)
cogsim_long$measure<-as.factor(cogsim_long$measure)
cogsim_long<-cogsim_long[complete.cases(cogsim_long),]
levels(cogsim_long$measure)<-c("Round1", "Round2","Round2","Round2","Round2")

cogsim_long_org<-cogsim_long %>% 
  group_by(subjID, showNUMcond, sim_cond, subj_SJT_choice, argument, Jaccard) %>% 
  dplyr::summarise(
    round1=measurement[measure=="Round1"],
    round2=measurement[measure=="Round2"],
    round1_neg=round1<0, 
    round2_neg=round2<0,
    cat_change=round1_neg!=round2_neg, 
    str_change=abs(round2) - abs(round1),
    total_change=abs(round2 - round1)
  )

cogsim_long_org$init_str<-abs(cogsim_long_org$round1)
cogsim_long_org_init_neg<-subset(cogsim_long_org, round1<1)
cogsim_long_org_init_neg$toward<-cogsim_long_org_init_neg$round2>cogsim_long_org_init_neg$round1
cogsim_long_org_init_pos<-subset(cogsim_long_org, round1>1)
cogsim_long_org_init_pos$toward<-cogsim_long_org_init_pos$round2<cogsim_long_org_init_pos$round1
cogsim_long_org<-rbind(cogsim_long_org_init_neg, cogsim_long_org_init_pos)
cogsim_long_org$revision_toward<-cogsim_long_org$total_change
cogsim_long_org[cogsim_long_org$toward==FALSE,]$revision_toward<-cogsim_long_org[cogsim_long_org$toward==FALSE,]$revision_toward * -1
cogsim_long_org$sim_cond<-as.factor(cogsim_long_org$sim_cond)
cogsim_long_org <- within(cogsim_long_org, sim_cond <- relevel(sim_cond, ref = 2))
cogsim_long_org$showNUMcond<-as.factor(cogsim_long_org$showNUMcond)
cogsim_long_org <- within(cogsim_long_org, showNUMcond <- relevel(showNUMcond, ref = 2))
cogsim_long_org$Dissimilar<-cogsim_long_org$sim_cond=="dissim"
cogsim_long_org$Hidden<-cogsim_long_org$showNUMcond=="No"
cogsim_long_org$Change_of_Stance<-cogsim_long_org$cat_change
cogsim_long_org$Initial_Conviction<-cogsim_long_org$init_str
cogsim_long_org$Magnitude_of_Revision_Toward_Opposing_Stance<-cogsim_long_org$revision_toward
cogsim_long_org$argument<-as.factor(cogsim_long_org$argument)

###############
#General Stats#
###############
cogsim_long_org %>% group_by(showNUMcond, sim_cond) %>% 
  dplyr::summarise(n = length(unique(subjID)),
                   probCHANGE=sum(cat_change)/length(cat_change), 
                   cilow=prop.test(sum(cat_change), length(cat_change), conf.level = 0.68)$conf.int[1],
                   cihi=prop.test(sum(cat_change), length(cat_change), conf.level = 0.68)$conf.int[2], 
                   round1 = mean(round1), 
                   round2 = mean(round2))
#########
#Figures#
#########

##########
#Figure 4#
##########
jaccard_plot<-cogsim_long_org
jaccard_plot$Jaccard_rev<-(1 - jaccard_plot$Jaccard)
jaccard_plot$sim_cond<-as.factor(jaccard_plot$sim_cond)
levels(jaccard_plot$sim_cond)<-c("Latent Cognitive\nSimilarity", 
                                 "Latent Cognitive\nDissimilarity")

ggplot(jaccard_plot, aes(x = Jaccard_rev, fill=sim_cond)) +
  geom_density(alpha=0.6, size=1.2) + theme_bw() + 
  scale_fill_manual(values=c("white", "grey")) + 
  ylab("Density") + xlab("Cognitive Distance\n(1 - Jaccard Index)") + 
  theme(legend.text=element_text(size=50),legend.position="top",
        legend.title = element_blank(),plot.title=element_blank(),
        axis.title.y=element_text(size = 50, hjust = 0.5),
        axis.title.x=element_text(size = 50, hjust = 0.5),
        axis.text.x=element_text(size = 60, hjust = 0.6), 
        axis.text.y=element_text(size = 50, hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_continuous(limits = c(0.15, 0.9)) 

###########
#Figure A1#
###########
cogsim$argument_num<-as.numeric(as.factor(cogsim$argument))

cogsimA_01<-subset(cogsim, 
                   subj_SJT_choice == 1 & 
                     option_A== "Point out to the employee how important his full commitment is to you. Openly communicate your criticism of the employee's current work ethics, but emphasize that you highly valued his performance on former projects.")
cogsimA_02<-subset(cogsim, 
                   subj_SJT_choice == 2 &
                     option_B== "Point out to the employee how important his full commitment is to you. Openly communicate your criticism of the employee's current work ethics, but emphasize that you highly valued his performance on former projects.")
cogsimA_simp<-rbind(cogsimA_01[,c("round1", "sim_cond", "showNUMcond")], cogsimA_02[,c("round1", "sim_cond", "showNUMcond")])
cogsimA_simp$Option<-1
cogsimA_simp$round1<- -1 * abs(cogsimA_simp$round1)

cogsimB_01<-subset(cogsim, 
                   subj_SJT_choice == 1 & 
                     option_A!= "Point out to the employee how important his full commitment is to you. Openly communicate your criticism of the employee's current work ethics, but emphasize that you highly valued his performance on former projects.")
cogsimB_02<-subset(cogsim, 
                   subj_SJT_choice == 2 & 
                     option_B!= "Point out to the employee how important his full commitment is to you. Openly communicate your criticism of the employee's current work ethics, but emphasize that you highly valued his performance on former projects.")

cogsimB_simp<-rbind(cogsimB_01[,c("round1", "sim_cond", "showNUMcond")], cogsimB_02[,c("round1", "sim_cond", "showNUMcond")])
cogsimB_simp$Option<-2
cogsimB_simp$round1<- abs(cogsimB_simp$round1)

option_distribution<-rbind(cogsimA_simp, cogsimB_simp)
option_distribution_sim<-subset(option_distribution, sim_cond=="sim")
option_distribution_dissim<-subset(option_distribution, sim_cond=="dissim")

ggplot(option_distribution_sim, aes(x=round1)) + 
  geom_histogram(colour="black", fill="grey") + 
  theme_bw() + theme(plot.title = element_text(size=50,  hjust=0.5), 
                     axis.text.x = element_text(size=50),
                     axis.text.y = element_text(size=50, hjust=0.5),
                     axis.title.x = element_text(size=50),
                     axis.title.y = element_text(size=50, angle=90, vjust=1),
                     strip.text.x = element_text(size=40),
                     legend.position = "top", 
                     legend.text=element_text(size=50),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     legend.title=element_blank(),
                     axis.line = element_line(colour = "black")) +
  ggtitle("Latent Cognitive Similarity") + 
  labs(y="Organizational Leaders", x="SJT Decision") + geom_vline(xintercept = 0, linetype="dotted", size=2) + 
  scale_x_continuous(limits=c(-52,52), breaks=c(-50,-25,0,25,50)) + 
  coord_cartesian(ylim=c(0,65))

sum(option_distribution_sim$round1<0)/nrow(option_distribution_sim)
sum(option_distribution_sim$round1>0)/nrow(option_distribution_sim)

ggplot(option_distribution_dissim, aes(x=round1)) + 
  geom_histogram(colour="black", fill="grey") + 
  theme_bw() + theme(plot.title = element_text(size=50,  hjust=0.5), 
                     axis.text.x = element_text(size=50),
                     axis.text.y = element_text(size=50, hjust=0.5),
                     axis.title.x = element_text(size=50),
                     axis.title.y = element_text(size=50, angle=90, vjust=1),
                     strip.text.x = element_text(size=40),
                     legend.position = "top", 
                     legend.text=element_text(size=50),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     legend.title=element_blank(),
                     axis.line = element_line(colour = "black")) +
  ggtitle("Latent Cognitive Dissimilarity") + 
  labs(y="Organizational Leaders",  x="SJT Decision") + geom_vline(xintercept = 0, linetype="dotted", size=2) + 
  scale_x_continuous(limits=c(-52,52), breaks=c(-50,-25,0,25,50)) + 
  coord_cartesian(ylim=c(0,65))

sum(option_distribution_dissim$round1<0)/nrow(option_distribution_dissim)
sum(option_distribution_dissim$round1>0)/nrow(option_distribution_dissim)
wilcox.test(option_distribution_sim$round1, option_distribution_dissim$round1)
wilcox.test(subset(option_distribution_sim, showNUMcond=="Yes")$round1, 
            subset(option_distribution_dissim, showNUMcond=="Yes")$round1)
wilcox.test(subset(option_distribution_sim, showNUMcond=="No")$round1, 
            subset(option_distribution_dissim, showNUMcond=="No")$round1)
wilcox.test(subset(option_distribution_dissim, showNUMcond=="Yes")$round1, 
            subset(option_distribution_dissim, showNUMcond=="No")$round1)
wilcox.test(subset(option_distribution_sim, showNUMcond=="Yes")$round1, 
            subset(option_distribution_sim, showNUMcond=="No")$round1)

option_distribution$combo_cond<-paste(option_distribution$sim_cond, option_distribution$showNUMcond, sep="_")
kruskal.test(round1 ~ combo_cond, data = option_distribution)

########
#Tables#
########
cogsim_long_org$argumentID<-as.character(as.numeric(as.factor(cogsim_long_org$argument)))
cogsim_long_org$Dissimilar<-as.numeric(cogsim_long_org$Dissimilar)
cogsim_long_org$Hidden<-as.numeric(cogsim_long_org$Hidden)
cogsim_long_org$Observable<-as.factor(cogsim_long_org$Hidden)
levels(cogsim_long_org$Observable)<-c(1,0) #note reverse coding of hidden condition
cogsim_long_org$Observable<-as.numeric(as.character(cogsim_long_org$Observable))
cogsim_long_org$Similar<-cogsim_long_org$sim_cond=="sim"
cogsim_long_org$Similar<-as.numeric(cogsim_long_org$Similar)

#########
#Table 2#
#########
table2 <- glm(data=subset(cogsim_long_org, Hidden==1), 
                   formula=Change_of_Stance ~ Dissimilar + argumentID, family="binomial")
tab_model(table2)
summary(table2)
exp(cbind("Odds ratio" = coef(table2), 
          confint.default(table2, level = 0.95)))
r2_nagelkerke(table2)

#########
#Table 3#
#########
table3<-lm(Magnitude_of_Revision_Toward_Opposing_Stance ~ 
                            Dissimilar + argumentID, subset(cogsim_long_org, Hidden==1))
summary(table3)
tab_model(table3, show.se = TRUE)

#########
#Table 4#
#########
table4 <- glm(data=cogsim_long_org, 
                       formula=Change_of_Stance ~ Dissimilar * Observable + argumentID,
                       family="binomial")
tab_model(table4, show.se = TRUE)
summary(table4)
exp(cbind("Odds ratio" = coef(table4), 
          confint.default(table4, level = 0.95)))
r2_nagelkerke(table4)

#########
#Table 5#
#########
table5_mod1<-lm(Magnitude_of_Revision_Toward_Opposing_Stance ~ Dissimilar * Observable + argumentID, cogsim_long_org)
summary(table5_mod1)
tab_model(table5_mod1, show.se = TRUE)

table5_mod2<-lm(Magnitude_of_Revision_Toward_Opposing_Stance ~ Initial_Conviction +  Dissimilar * Observable + argumentID, cogsim_long_org)
summary(table5_mod2)
tab_model(table5_mod2, show.se = TRUE)

####Robustness: Frequency of argument across conditions####
cogsim_long_org_arg_agg<-cogsim_long_org %>% 
  group_by(sim_cond, argument, showNUMcond) %>% 
  dplyr::summarise(num_subjects=length(unique(subjID)))

mean(cogsim_long_org_arg_agg$num_subjects)/250
sd(cogsim_long_org_arg_agg$num_subjects)/250

cogsim_long_org_arg_agg_observable<-subset(cogsim_long_org_arg_agg, showNUMcond=="No")

cogsim_long_org_arg_agg_observable_sim<-subset(cogsim_long_org_arg_agg_observable, sim_cond=="sim")
cogsim_long_org_arg_agg_observable_dissim<-subset(cogsim_long_org_arg_agg_observable, sim_cond=="dissim")

wilcox.test(cogsim_long_org_arg_agg_observable_sim$num_subjects, 
            cogsim_long_org_arg_agg_observable_dissim$num_subjects, paired=T)

cogsim_long_org_arg_agg_obscured<-subset(cogsim_long_org_arg_agg, showNUMcond=="Yes")
cogsim_long_org_arg_agg_obscured_sim<-subset(cogsim_long_org_arg_agg_obscured, sim_cond=="sim")
cogsim_long_org_arg_agg_obscured_dissim<-subset(cogsim_long_org_arg_agg_obscured, sim_cond=="dissim")

wilcox.test(cogsim_long_org_arg_agg_obscured_sim$num_subjects, 
            cogsim_long_org_arg_agg_obscured_dissim$num_subjects, paired=T)

###############################
#Supplementary Mechanism Study#
###############################
arg_rate<-read.csv(paste(datapath, "supp_dt.csv", sep=""))
arg_rate$ArgID<-as.factor(arg_rate$ArgID)
cor.test(arg_rate$Jaccard, arg_rate$novel)
cor.test(arg_rate$Jaccard, arg_rate$decide)
cor.test(arg_rate$novel, arg_rate$decide)

#Mediation Analysis

#Total effect
#all_data_clean$Argument<-as.factor(all_data_clean$arg_ID)
arg_rate$Cog.Dissim<-1-arg_rate$Jaccard
mod1<-lm(decide ~ Cog.Dissim + ArgID + SJT_dummy, data=arg_rate)
summary(mod1)

#effect of IV on mediator
mod2<-lm(novel ~ Cog.Dissim + ArgID + SJT_dummy, data=arg_rate)
summary(mod2)

#effect of mediator on the dependent variable
mod3<-lm(decide ~ Cog.Dissim + novel + ArgID + SJT_dummy, data=arg_rate)
summary(mod3)

#causal mediation analysis 
results = mediate(mod2, mod3, treat='Cog.Dissim', mediator='novel', boot=T)
summary(results)

##########
#Table C1#
##########
tablec1<-lm(decide ~ Jaccard + novel + convincing + interesting + informative + logical + emotional + SJT + ArgID, data=arg_rate)
summary(tablec1)
tablec1_vcov <- cluster.vcov(tablec1, arg_rate$SubjID)
coeftest(tablec1, tablec1_vcov)

###########
#Figure A2#
###########
subjects<-unique(arg_rate$SubjID)

compare_subj_df<-data.frame()

n<-0
for(subj_i in subjects){
  n<-n+1
  print(n)
  subj_i_df<-subset(arg_rate, SubjID == subj_i)
  
  subjects_j<-subjects[subjects != subj_i]
  for(subj_j in subjects_j){
    subj_j_df<-subset(arg_rate, SubjID == subj_j)
    subj_merge<-merge(subj_i_df, subj_j_df, by=c("Choice", "SJT", "Argument", "Arg.Content", "Arg.Bats"))
    compare_subj_df<-rbind(compare_subj_df, subj_merge)
  }
}


compare_subj_df$subj_jaccard<-sapply(1:nrow(compare_subj_df), function(x) jaccard(strsplit(compare_subj_df[x,]$Subj.BATs.x,",")[[1]], strsplit(compare_subj_df[x,]$Subj.BATs.y,",")[[1]]))
compare_subj_df$euc_dist<-sapply(1:nrow(compare_subj_df), function(x) as.numeric(dist(rbind(as.numeric(compare_subj_df[x,7:14]), as.numeric(compare_subj_df[x,21:28])))))
cor.test(compare_subj_df$subj_jaccard, compare_subj_df$euc_dist)

mod<-lm(euc_dist ~ subj_jaccard + ArgID.x, data=compare_subj_df)
summary(mod)

#Statistical tests 
compare_subj_df$comp_Jaccard_bin<-ntile(compare_subj_df$subj_jaccard,10)

compare_subj_df_agg<-compare_subj_df %>% group_by(comp_Jaccard_bin) %>%
  dplyr::summarise(euc_dist=mean(euc_dist, na.rm=T),subj_jaccard=mean(subj_jaccard))

ggplot(compare_subj_df_agg, aes(y = euc_dist, x=comp_Jaccard_bin)) + dougtheme_mod +  
  geom_point(colour="black", size = 8, position=position_dodge(0.8)) +
  geom_smooth(method='lm', formula= y~x, data=compare_subj_df, size=2, color="red") + 
  theme(axis.text.x = element_text(size=35),
        axis.text.y = element_text(size=35),
        axis.title.x = element_text(size=35, hjust=0.5, vjust=1),
        axis.title.y = element_text(size=35, hjust=0.5, vjust=1),
        strip.text.x = element_text(size =25), 
        legend.position=c(0.05,0.9),  
        legend.text = element_text(size=30),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
  ) + labs(x="Cognitive Similarity of Participants \n(Deciles)", 
           y = "Dissimilarity of Participants' Argument Ratings\n(Euclidean Distance)", linetype=NULL) + 
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10)) + 
  coord_cartesian(ylim=c(4.2, 4.31))
