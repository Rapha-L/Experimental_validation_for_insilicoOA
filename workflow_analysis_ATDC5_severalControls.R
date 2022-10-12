#  This file is modified from of the `BraDiPluS` R package (F. Educati, Website: https://github.com/saezlab/BraDiPluS)
#
#
#  File author(s): Raphaelle Lesage (raphaelle.lesage@kuleuven.be)
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://github.com/Rapha-L/Experimental_validation_for_insilicoOA 
#
# --------------------------------------------------------
#
#' Compute and plot z-score and p-value.

# in my ATDC5 dataset there are some conditions to be compared to a different control (DMSO) this script is made to compare
# the right conditions to the right control and display all result in the same boxplot graph

library(here)
here::i_am("workflow_analysis_ATDC5_severalControls.R")

  
library(corrplot)
library(ggplot2)
library(gridExtra)
library(corrplot) # note: this uses a slightly modified version of corrplot package (to handle NAs)
if (showLabels==T){
  library(ggrepel)
}
library(reshape2)
  
### Create dataset in List
  
#load  data in dataframe

allData.df = read.csv(here("ALP_data_FC_CMCM2_Fig5D.csv"),header = TRUE)

  
runs = c("run2","run3","run4","run5","run6","run7")
conditions = as.character(subset(allData.df,run== "run2")[,"conditions"])
ATDC5 = list()
  
# format dataset in List  
  for (r in runs){
    a =subset.data.frame(allData.df,run == r)
    MyRun = list()
    for (d in c(1:dim(a)[1])){
      temp = t(a[d,c('X','X.1','X.2')])
      rownames(temp) =NULL
      x =list(temp)
      x=as.data.frame(x)
      colnames(x) = "ALP"
      
      MyRun[[d]]= x
      
    } 
    MyRun = as.list(setNames(MyRun, conditions))
    rname = as.character(r)
    ATDC5[[rname]] = MyRun
  }
  
  as.list(setNames(ATDC5, runs))
  allData = list("ATDC5" = ATDC5)
  
  
  
##### manual computation : statistics & graphs

  myData = allData$ATDC5
  controlName = c("CM + CM","CM2 + CM2")
  nb_controls = length(controlName)
  subsample= F
  thSD = NA
  thPval = 0.05
  showLabels=T
  
  # minimum number of runs (for subsampling when computing combined p-value)
  minRuns<-min(sapply(allData, length))
  
  ## runs medians 
  runs_medians<-do.call(cbind, lapply(myData, function(myRun){
    apply(sapply(myRun, get, x="ALP"),2,median,na.rm = T)
  }))
  
  ## runs p-values
  
  runs_pvalue<-do.call(cbind, lapply(myData, function(myRun){
    
    idCM= list()
    ix_control= list()
    for (contr in 1:nb_controls){
      # store controls index
      ix_control[[contr]]<- which(names(myRun)==controlName[contr])
      
      # index of conditions related to each control
      myCM = levels(allData.df$Medium)[contr]
      idCM[[contr]] = as.vector(unlist(subset(allData.df,(Medium == myCM & run == "run2"),select = "Id")))
    }
    
    allControlRep = list()
    thePvals = list()
    for (contr in 1:nb_controls){
      ix = ix_control[[contr]]
      allControlRep_temp<-do.call(c, lapply(myRun[ix], function(x){
        x$ALP
      }))
      allControlRep[[contr]] = allControlRep_temp
      
      thePval_temp <-sapply(myRun[idCM[[contr]]], function(x){
        if (nrow(x)>=3 & !is.na(x)){ # needs at least 3 replicates
          wilcox.test(x$ALP, allControlRep_temp, alternative="less")$p.value
        }else{
          return(NA)
        }
      })
      thePvals[[contr]] <-p.adjust(thePval_temp, "BH")
    }
    
    thePval = NULL
    for (contr in 1:nb_controls) {
      thePval = c(thePval , thePvals[[contr]])
    }
    
    return(thePval)
  }))
  
  
  runs_pvalue_combined<-apply(runs_pvalue, 1, function(x){
    # NOTE: I could add as weight the sample size of each study
    # subsample runs when computing combined p-val (avoid differences in volcano just due to number of runs)
    if (length(x)>minRuns & subsample==T){x<-sample(x,2, replace = F)}
    run_number = dim(runs_pvalue)[2] - apply(apply(runs_pvalue,1,is.na),2,sum)
    survcomp::combine.test(x, weight = run_number, method = "fisher", na.rm = T)
    # mean(x, na.rm = T)
  })
  
  # compute corresponding z-score
  runs_zscore<-apply(runs_medians,2,scale)
  rownames(runs_zscore)<-rownames(runs_medians)
  # and the median z-score
  runs_zscore_median<-apply(runs_zscore, 1, function(x) median(x, na.rm = T))
  runs_zscore_sd<-apply(runs_zscore, 1, function(x) sd(x, na.rm = T))
  
  # find the control samples
  samplesNames<-rownames(runs_zscore)
  ix_control<- samplesNames %in% controlName

  #keep Zscore of controls for each run & put all runs in one single column
  controlData<-c(runs_zscore[ix_control,])
  controlData<-controlData[!is.na(controlData)]    
  
  
  # and remove those with high variability among replicates (put NA)
  if (!is.na(thSD)){
    ix_highSD<-which(runs_zscore_sd>thSD)
  }else{
    ix_highSD<-c()
  }
  
    # remove the control and the samples with high SD
    ix_remove<-c(ix_control, ix_highSD)
    runs_zscore_median_rm<-runs_zscore_median[-ix_remove]
    runs_pvalue_combined_rm<-runs_pvalue_combined[-ix_remove]
  
  ###### Results visualization 
    
  # prepare for boxplot
  is.significant<-as.factor(runs_zscore_median<0 & runs_pvalue_combined <= thPval)
  is.significant[ix_highSD]<-FALSE
  is.control<- as.factor(names(is.significant) %in% controlName)
  nRuns<-ncol(runs_zscore)
  dataBoxplot<-data.frame(ID=as.factor(rep(seq(1:nrow(runs_zscore)), nRuns)),
                          names=rep(names(is.significant),nRuns),
                          run=rep(colnames(runs_zscore), each=nrow(runs_zscore)),
                          values= do.call(rbind, lapply(split(runs_zscore, col(runs_zscore)), as.matrix)),
                          is.significant=rep(is.significant, nRuns),
                          is.control=rep(is.control, nRuns))
  
  library(extrafont)
  loadfonts(device = "win")
  
# define color of x labels according to in silico predictions
  Mylabel_colors2 = c(rep("grey40"),rep(c("blue","grey40"),5),"blue","blue", rep(c("grey40","blue"),3)) #depends on the order of conditions on the graph
  
  g_box <- ggplot(dataBoxplot, aes(x=ID, y=values, colour=is.significant, fill=is.control)) + geom_hline(yintercept=0, colour="grey70") +
    geom_boxplot(outlier.colour=NA) +
    geom_point(position = position_jitter(width = 0.2)) + theme_bw() +
    scale_y_continuous(sec.axis = dup_axis(name = waiver())) +
    ylab("z-score") + xlab("") +
    scale_x_discrete(breaks = factor(seq(1:nrow(runs_zscore))), labels = names(is.significant)) +
    theme(axis.title = element_text(size = 15, family = "Calibri"), axis.text.y =element_text(size = 15), axis.text.x=element_text(angle = -55, hjust = 0,color = Mylabel_colors2, size = 14), legend.position = "top", legend.text = element_text(size = 15)) +
    scale_colour_manual(values = c("grey81", "grey25"), name="", breaks=levels(is.significant), labels=c("not significant", "significant")) +
    scale_fill_brewer(palette="Accent", name="        ", breaks=levels(is.control), labels=c("drug response sample\n(single drug or pairwise combination)", "control sample\n(only Control medium, no drugs)"))
  
  g_box
  ggsave("ALP_invitro_validation_Fig5D.tiff",g_box, device= "tiff", dpi = 330)
  
