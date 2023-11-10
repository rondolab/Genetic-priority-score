library(plyr)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(DataCombine)
library(patchwork)
library(grid)
library(gridExtra)
library(gtable)
library(scales)
library(Hmisc)
library(corrplot)

#Create Fig 2 and Supl Fig. 1
#Create two plots, All drugs and drugs restricted to one gene target for Open Targets
lapply(c('allgenes', 'onegene'), function(file_type) {
    # Univariate association file
    Forestplotfile<- fread(paste0(Revision_folder_final, 'Univar_regression_Opentargets_',file_type,'.txt'), data.table=F)
    #Use filled circles for significant pvalues and open circles for non significant pvalues
    Forestplotfile$Pval_sig=ifelse(Forestplotfile$P.val>=0.05,0,1)
    Forestplotfile$Pval_sig<-as.factor(Forestplotfile$Pval_sig)
    if(length(unique(Forestplotfile$Pval_sig[!is.na(Forestplotfile$Pval_sig)]))==1){
    shape_vals=c(19)
    } else{
    shape_vals=c(1,19) }
    #Order of Y axis
    Forestplotfile$order=rep(length(unique(Forestplotfile$predictors)):1)
    Forestplotfile$predictors<-factor(Forestplotfile$predictors, levels=unique(Forestplotfile$predictors[order(Forestplotfile$order, decreasing=T)]))
    #For large CI add arrows rather than showing full bar. Add arrow at max OR + 5 
    maxor=max(Forestplotfile$OR)+5
    Forestplotfile$upperCI_cut<-ifelse(Forestplotfile$upperCI> maxor, maxor, NA)
    # include OR and CI text
    Forestplotfile$CI= paste0(round(Forestplotfile$OR,1),' (', round(Forestplotfile$lowerCI,1), ' - ', round(Forestplotfile$upperCI,1),')')
    Forestplotfile$CI<-factor(Forestplotfile$CI, levels=unique(Forestplotfile$CI[order(Forestplotfile$order, decreasing=T)]))
    #Create plot
    fp <- Forestplotfile %>% as_tibble() %>%
        ggplot(aes(y=predictors, x=OR))+
        geom_point(size=2, aes(color=Category, shape=Pval_sig),stroke =1.5) +
        scale_shape_manual(values=shape_vals) +guides(shape = "none") + 
        geom_linerange(aes(xmin=lowerCI, xmax=upperCI),color='black') +
        xlab("Odds ratio (95% CI)") + ylab('') +  coord_cartesian(c(0,maxor)) + #limit x axis from 0-maxOR values 
        scale_y_discrete(limits = rev(levels(Forestplotfile$predictors))) +
        geom_vline(xintercept=1, linetype='longdash', color='red') + 
        geom_text(aes(x=0, y=order+0.2), label = paste(round(Forestplotfile$No.genes,2)),size =2.5, color='#D41159')+ # add label for number of genes
        geom_text(aes(x=0, y=order-0.2), label = paste(round(Forestplotfile$No.parentterms,3)),size =2.5, color='#1A85FF') +# add label for number of predictors
        theme_classic() + 
        theme(axis.text=element_text(size=14,color = "black"), axis.title=element_text(size=15,color = "black"), legend.position="none", 
          axis.line = element_line(colour = 'black', size = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.5),
          axis.ticks.length = unit(0.25, "cm")) #plot formatting 

    #For large CI add arrows instead of full line
    fp1<-fp + scale_x_continuous(limits=c(0,maxor)) +
        geom_segment(aes(x = lowerCI , xend = upperCI_cut, y = predictors, yend = predictors ), arrow = arrow(length = unit(0.25, "cm")))
    # include OR and CI text
    p_right <- 
        Forestplotfile %>%
        ggplot() + geom_text(aes(x = 0, y = predictors, label = CI),size =3,hjust = 0 ) + 
        ggtitle('OR (95% CI)') + scale_y_discrete(limits = rev(levels(Forestplotfile$predictors))) +
        theme_void() + theme(plot.title = element_text(size =10, hjust = 0.75, face="bold"))
    #paste all plots together
    fp3 <- fp1 + plot_spacer() + p_right + plot_layout(widths = c(5, -1.1 ,3) ,guides = "collect")& theme(legend.position = "bottom",legend.justification='left') 
ggsave(fp3, file=paste0('Forestplot_Univar_Opentargets_',file_type,'.pdf'), width = 7, height=6)
})

### Fig 3, Extended data Fig 4 
lapply(c('Opentargets','SIDER'), function(datafile){ 
    binned_gps=fread(paste0(Revision_folder_final,'Binned_by_sum_binsize0.3_',datafile,'.txt'), data.table=F)
    if(length(unique(binned_gps$or_sig[!is.na(binned_gps$or_sig)]))==1){
    shape_vals=c(19)}
     else{
    shape_vals=c(1,19)}
    max_uci <- round(max(binned_gps$upperCI))
    #Create plot
    plot1=ggplot(data=binned_gps, aes(x=genescoresum, y=OR)) +
        geom_point(aes(color=OR,shape =or_sig),position=position_dodge(width=0.2),size=2.5, stroke=1.5)  + xlab('Genetic priority score bins') + ylab('Odds ratio (95% CI)') +
        scale_shape_manual(values=shape_vals) +guides(shape = "none") + 
        scale_x_discrete(limits = (levels(binned_gps$genescoresum))) +
        geom_linerange(aes(ymin=lowerCI, ymax=upperCI,color=OR),position=position_dodge(width=0.2)) +
        scale_color_gradient(low = ("blue"), high = ("red"), limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, max_value), labels = c(0,2, 4, 6, 8, max_value), oob = scales::squish  ) + #OR colour scale
        geom_hline(yintercept=1, color = "grey") +
        coord_cartesian(ylim=c(-0.5,max_uci), clip='off') +
        theme_classic() + 
        theme(axis.text.x=element_text(size=13,angle=45, vjust = 1, hjust=1,color = "black"),axis.text.y=element_text(color = "black",size=13),
          axis.title=element_text(size=15, color = "black"),legend.text=element_text(size=13), legend.title=element_text(size=13),
          axis.line = element_line(colour = 'black', size = 0.5), 
          axis.ticks = element_line(colour = "black", size = 0.5),
          axis.ticks.length = unit(0.25, "cm"),
          plot.margin=margin(10,20,10,10)) +
        scale_y_continuous(limits=c(0,max_uci)) 
    ggsave(plot1, file=paste0('Increase_GPS_',datafile,'_allgenes.pdf'), width = 7, height=7)

})


#Fig. 4
Simulation1=fread(paste0('Simulation_highgps_matchgenes_mi_sider.txt'), data.table=F)
#seperate data into Null and when a GPS threshold is used.
Nulldata=Simulation1[Simulation1$data=='Null',]
Abovethreshdata=Simulation1[Simulation1$data=='Above',]
yaxis=Nulldata %>%mutate(seround=round(MI_percent)) %>% select(threshold,Sampleon ,seround )  %>% group_by(seround,threshold,Sampleon )%>% tally()
#calculate empirical pvalue at each threshold
pvalcal1<-do.call(rbind,lapply(c(unique(Abovethreshdata$threshold)), function(thresh){
      Nulldata_threshold=Nulldata %>% filter(threshold==thresh)
      pvalcalc=data.frame(pval=length(Nulldata_threshold$MI_percent[Nulldata_threshold$MI_percent>Abovethreshdata$MI_percent[Abovethreshdata$threshold==thresh] ])/length(Nulldata_threshold$MI_percent))
      pvalcalc$threshold=thresh 
      return(pvalcalc)
}))

colnames(Abovethreshdata)[c(2,3,7)]=gsub('^*','Above_',colnames(Abovethreshdata)[c(2,3,7)])
Null_1= inner_join(Nulldata, Abovethreshdata[c(2,3,5,7)])
Null1=inner_join(Null_1, pvalcal1)

Null1$threshold=gsub('^', '>', Null1$threshold)
#position on x axis for label
Null1$xaxis=ifelse(Null1$Above_MI_percent<20,Null1$Above_MI_percent+2, Null1$Above_MI_percent-2 )
Null1$pval=round(Null1$pval,3)
histdata <- hist(Null1$MI_percent[Null1$threshold=='>0.3'],breaks=10, plot=FALSE, right=FALSE)
max.value <- max(histdata$counts)
max.value  # to position label on y axis. get maxium histogram count.  
#calculate mean percent of main indications %. To create null distribution 
mean_null_mi=Null1 %>% distinct(Iteration, MI_percent, threshold,xaxis)  %>% mutate(OR=(gsub('>','', threshold))) %>% group_by(threshold,xaxis,OR) %>% summarise(mean_mi=mean(MI_percent))  %>% arrange((OR))

xlim <- c(0, max(Null1$Above_MI_percent)) #limits for x axis
extra=max(Null1$Above_MI_percent)/4 #use for label positioning on x axis
Null1$pval= gsub('^','P-value = ', Null1$pval) 
Null1$pval[Null1$pval=='P-value = 0']<-'P-value < 0.001'
Null1$OR=(gsub('>','',Null1$threshold))
#To use >= sign. 
Null1$threshold2<-factor(Null1$threshold, labels=paste('GPS \u2265' ,unique(Null1$OR)))
mean_null_mi$threshold2<-factor(mean_null_mi$threshold, labels=paste('GPS \u2265' ,unique(mean_null_mi$OR)))
Null1<- inner_join(Null1, mean_null_mi[c('threshold2','mean_mi')])
### calculate fold ratio
Null1$ratio=round(Null1$Above_MI_percent/ Null1$mean_mi,2)

#Create plot
plot_null<-ggplot(data=Null1, aes(x=MI_percent)) + geom_histogram(binwidth = 0.5, color="black", fill='#619CFF') + theme_classic() +  facet_wrap(~threshold2)+ 
  geom_vline(aes(xintercept=Above_MI_percent), lwd=1, linetype=2, color="red") +
  scale_x_continuous(limits = xlim, oob = scales::squish) + #start at 0 
  geom_vline(data=mean_null_mi, aes(xintercept=mean_null_mi$mean_mi), lwd=1, linetype=2, color="blue") +
  xlab('Percent drug indications (%)') + ylab('Number of simulations') +
  geom_label(data=Null1,aes(x=max(Above_MI_percent)-extra, y=max.value-(max.value/2.5)), label = paste("Ratio =", Null1$ratio), col = "red",size =2.5) + #position ratio label on each facet
  geom_label(aes(x=max(Above_MI_percent)-extra, y=max.value - (max.value/7), fontface=3), label = paste(Null1$pval), col = "black",size =2.5)  +  #position pvalue label on each facet
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
    axis.ticks = element_line(colour = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    axis.text=element_text(color = "black",),
    axis.title=element_text(color = "black"),
    plot.margin=margin(10,20,10,10),
    legend.title = element_text(size=11), legend.text = element_text(size=11),
    legend.background = element_blank(),legend.box.background = element_rect(colour = "black")) 

ggsave(plot_null,file=paste0(plots_dir,'/Null_histogram_simulation_MIS_1000_threshold_parentterm_',file_type,'_',datafile,'_',run,'_binned_genescoresum_paper_',cancerd,'16_3.png'), width = 12, height=10)


#Fig 5 a and Extended Fig 8a
lapply(c('GPS','GPS-DOE'),function(score){
  clinicalphase<-fread(paste0(Revision_folder_final,'Stratified_clinicalphase_genescoresumregression_',score,'.txt'), data.table=F)
  #Format file
  clinicalphase$Pval_sig_raw=ifelse(clinicalphase$P.val_raw>=0.05,0,1)
  clinicalphase$Pval_sig_raw<-as.factor(clinicalphase$Pval_sig_raw)
  clinicalphase$clinicalphase=c('All','I','II','III','IV')
  clinicalphase$order=rep(length(unique(clinicalphase$clinicalphase)):1)
  clinicalphase$clinicalphase[-1] =gsub('^', 'Phase ', clinicalphase$clinicalphase[-1] )
  clinicalphase$clinicalphase<-factor(clinicalphase$clinicalphase, levels=unique(clinicalphase$clinicalphase))
#Create plot
  fp <- clinicalphase %>% as_tibble() %>%
    ggplot(aes(y=clinicalphase, x=OR_raw))+
    geom_point(size=3,stroke =1.5,aes(shape=Pval_sig_raw), color='#619CFF') +   
    scale_shape_manual(values=c(19,1)) +guides(shape = "none") + 
    geom_linerange(aes(xmin=lowerCI_raw, xmax=upperCI_raw)) +
    xlab("Odds ratio (95% CI)") + ylab('') + xlim(0, max(clinicalphase$upperCI_raw)) +
    scale_y_discrete(limits = rev(levels(clinicalphase$clinicalphase))) +
    geom_vline(xintercept=1, linetype='longdash', color='red') +
    geom_text(aes(x=0, y=order+0.2), label = paste(round(clinicalphase$No.MI,2)), size =3, color='#D41159')+
    theme_classic(base_size = 16) +
    theme(axis.text=element_text(size=13,color = "black"), plot.title = element_text(size = 13, face = "bold"), axis.title=element_text(size=15,color = "black"),
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(0.25, "cm"),
      plot.margin=margin(10,20,10,10)) 
    
  ggsave(fp, file=paste0(plots_dir,'Stratified_clinicalphase_genescoresumregression_',score,'.pdf'), width = 6, height=6)

})

#Fig 5 b and Extended Fig 8b
## Fold dif clinical phase
lapply(c('GPS','GPS-DOE'),function(score){

  folddif<-fread(paste0(,'Fold_enrichment_phasescomparedtophase1_',score,'.txt'), data.table=F)
  fold_dif2=reshape2::melt(folddif, variable.name='predictors')
  # fold_dif2$predictors=gsub('Fold_dif_','',fold_dif2$predictors)
  #restrict to thesholds 0.9, 1.5 and 2.1
  fold_dif1=fold_dif2 %>% filter( predictors%in% c(0.9, 1.5,2.1))
  fold_dif1$predictors=factor(fold_dif1$predictors, levels=unique(fold_dif1$predictors))
  #Create plot
  folddif_plot=ggplot(fold_dif1, aes(y=value, x=phase, fill=predictors)) + 
  geom_bar(position="dodge", stat="identity",color='black') + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_brewer(palette="Blues") +
  xlab('') + ylab('Fold difference')  + guides(fill=guide_legend(title=paste0(score, " bin\n cutoff"))) + 
  theme_classic(base_size = 16) +
  theme(axis.text=element_text(size=13,color = "black"),
    axis.text.x=element_text(color = "black", angle=45, vjust = 1, hjust=1),
    plot.title = element_text(size = 13, face = "bold"),
    axis.title=element_text(size=15,color = "black"),
    axis.line = element_line(colour = 'black'), 
    axis.ticks = element_line(colour = "black"),axis.ticks.length = unit(0.25, "cm"),
    plot.margin=margin(10,20,10,10)) 

  ggsave(folddif_plot, file=paste0(plots_dir, 'Clinicalphase_data_folddif_divisblephase1_Opentargets_score',score,'.pdf'), width = 6, height=6)
})

# Fig 6 and Extended Fig 7
lapply(c('Opentargets','SIDER'), function(datafile){
  binned_gps2=fread(paste0('Binned_by_sum_binsize0.3_',datafile,'_doe.txt'), data.table=F)
  binned_gps2$or_sig<-ifelse(binned_gps2$P.val>=0.05,0,1)
  binned_gps2$or_sig<-as.factor(binned_gps2$or_sig)
  binned_gps2$genescoresum=factor(binned_gps2$genescoresum)
  if(length(unique(binned_gps2$or_sig[!is.na(binned_gps2$or_sig)]))==1){
  shape_vals=c(19)
  } else{
    shape_vals=c(1,19)
  }
  #Format facet label titles
  rep_str = c('abs'='GPS-DOE','neg'='Gain-of-function\nonly','pos'='Loss-of-function\nonly')
  binned_gps2$DOE <- str_replace_all(binned_gps2$DOE, rep_str)
  binned_gps2$DOE<-factor(binned_gps2$DOE, levels=unique(binned_gps2$DOE))

  max_uci <- round(max(binned_gps2$upperCI))+1

  plot1=ggplot(data=binned_gps2, aes(x=genescoresum, y=OR)) + facet_wrap(~ DOE) + 
    geom_point(aes(color=OR,shape =or_sig),position=position_dodge(width=0.2),size=2, stroke=1)  + xlab('Genetic priority score bins') + ylab('Odds ratio (95% CI)') +
    scale_shape_manual(values=shape_vals) +guides(shape = "none") + 
    scale_x_discrete(limits = (levels(binned_gps2$genescoresum))) + scale_y_continuous(limits=c(0,max_uci)) +
    geom_linerange(aes(ymin=lowerCI, ymax=upperCI,color=OR),position=position_dodge(width=0.2)) +
    scale_color_gradient(low = ("blue"), high = ("red"), limits = c(0, 10),breaks = c(0, 2, 4, 6, 8, max_value),
      labels = c(0,2, 4, 6, 8, max_value),oob = scales::squish  ) +
    geom_hline(yintercept=1, color = "grey") +
    coord_cartesian(ylim=c(-0.5,max_uci), clip='off') +
    theme_classic() +
    theme(axis.text.x=element_text(size=13,angle=45, vjust = 1, hjust=1,color = "black"),axis.text.y=element_text(color = "black",size=13),
      axis.title=element_text(size=15,color = "black"), 
      legend.text=element_text(size=13), legend.title=element_text(size=13), 
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(0.25, "cm"),
      plot.margin=margin(10,20,10,10)) 

  ggsave(plot1, file=paste0('Increase_GPS_',datafile,'_allgenes_GPS_DOE.pdf'), width = 7, height=7)
})


### Extended Fig 1 
firthreg<- fread(paste0(Revision_folder_final, 'Firth_regression_Opentargets.txt'), data.table=F)
#Use filled circles for significant pvalues and open circles for non significant pvalues
if(length(unique(firthreg$beta_sig[!is.na(firthreg$beta_sig)]))==1){
shape_vals=c(19)
} else{
  shape_vals=c(1,19)}
# include plot facet labels
rep_str = c('allgenes'='All drugs','onegene'='Drugs with one gene target')
firthreg$file_type <- str_replace_all(firthreg$file_type, rep_str)
#Create plot                 
fp <- firthreg %>% as_tibble() %>%
    ggplot(aes(y=Predictor, x=beta,color=CV))+
    facet_wrap(~file_type) + # Create panels for plot for 'All drugs' and 'Drugs with one gene target'
    geom_point(size=1.4, stroke=0.5,aes(color=CV, shape =beta_sig), position =position_dodge(0.5))  +
    scale_shape_manual(values=shape_vals) +guides(shape = "none") + 
    geom_linerange(aes(xmin=lowerCI, xmax=upperCI),position = position_dodge(0.5)) +
    scale_y_discrete(limits = rev(levels(firthreg$Predictor))) +xlim(-0.6,max(firthreg$upperCI)) + xlab('Beta coefficients (95% CI)') + ylab('')+
    geom_vline(xintercept=0, linetype='longdash', color='red') +
    theme_classic()
    theme(axis.text=element_text(size=13,color = "black"), axis.title=element_text(size=14,color = "black"),
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(0.25, "cm"),
      plot.margin=margin(10,20,10,10))

ggsave(fp, file= 'Firth_regression_Opentargets.jpg', width = 7, height=6, dpi=300)

# Extended Fig 2
dataset<-fread('Drug_Dataset_Opentargets_allgenes_final_mi.txt.gz',data.table=F)
genescoredataset=fread(paste0('All_genescoresum_across_all_predictors_opentargets.txt.gz'), data.table=F) #GPS OT data 
Dataset_genescores1= genescoredataset %>%arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
#restrict dataset to phenotype predictor columns and add percent GPS
Dataset_genescores1_df=Dataset_genescores1[,grepl('phenotype', colnames(Dataset_genescores1))]
Dataset_genescores1_df=cbind(percent=Dataset_genescores1$percent, Dataset_genescores1_df) 
#format dataset
Dataset_genescores1_df_formatted=pivot_longer(data = Dataset_genescores1_df, cols = -c('percent'), names_to = "Predictor", values_to = "weights") #format dataset using percent name column and collapse phenotype columns to rows
Dataset_genescores1_df_formatted=as.data.frame(Dataset_genescores1_df_formatted[Dataset_genescores1_df_formatted$weights!=0,])  #restrict to percentiles with scores > 0
Dataset_genescores1_df_formatted$Predictor<-factor(Dataset_genescores1_df_formatted$Predictor, levels=unique(Dataset_genescores1_df_formatted$Predictor[order(Dataset_genescores1_df_formatted$weights, decreasing=TRUE)]))
weights=Dataset_genescores1_df_formatted %>% distinct(Predictor, weights,file_type) %>% group_by(Predictor,file_type) %>% dplyr::summarise(median.weights=median(weights))%>% arrange(desc(median.weights))  #get median weight for each predictor (median across the 5 CV)
weights$order=(nrow(weights):1)+0.2
Dataset_genescores1_df_formatted$Predictor<- factor(Dataset_genescores1_df_formatted$Predictor, levels = unique(weights$Predictor[order(weights$median.weights, decreasing=T)]))
## add sample size for each predictor
sample_size_allpred = Dataset_genescores1_df_formatted%>% group_by(Predictor) %>% tally() %>% dplyr::rename(num=n) 
Dataset_genescores1_df_formatted2=Dataset_genescores1_df_formatted %>% left_join(sample_size_allpred) %>%
  mutate(samplesize = paste0(Predictor, "\n", "(n=", num,")"))
weights1=weights %>% left_join(sample_size_allpred) %>%
  mutate(samplesize = paste0(Predictor, "\n", "(n=", num,")"))
Dataset_genescores1_df_formatted2$samplesize<- factor(Dataset_genescores1_df_formatted2$samplesize, levels = unique(weights1$samplesize[order(weights1$median.weights, decreasing=T)]))
#2. Make plot
plot4 <- ggplot(data = Dataset_genescores1_df_formatted2, mapping = aes(x =percent,y = samplesize, fill=Predictor)) +
     geom_violin(position = position_dodge(width = 0.2), scale = "width", width=0.3,adjust = 0.7, trim = T,size=0.2)+
     scale_x_continuous( limits=c( min(Dataset_genescores1_df_formatted2$percent), 100)) + 
     scale_y_discrete(limits = rev(levels(Dataset_genescores1_df_formatted2$samplesize))) +
     stat_summary(fun.data = "mean_cl_boot", geom = "point", shape=21, size=3.5, color='black') +
     guides(fill = F, color = F) + 
    geom_text(data= weights, aes(x=min(Dataset_genescores1_df_formatted2$percent), y=order) , label = paste(round(weights1$median.weights,3)),size=3) +
    xlab('Percentile of genetic priority score (%)') + ylab('')  + 
    scale_fill_manual(values=c('OMIM'='#DB72FB','HGMD'='#619CFF','ClinVar'=  '#00BA38', 'Gene Burden' = '#D39200',
      'Single Variant'=  '#00C19F', 'eQTL Phenotype'=  '#F8766D', 'Locus2gene'= '#00B9E3',  'pQTL Phenotype'=  '#FF61C3')) +
    theme_classic() +
    theme(axis.text=element_text(size=14,color = "black"),
      axis.title=element_text(size=15,color = "black"),
      legend.position="none",
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.ticks= element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(0.25, "cm"),
      plot.margin=margin(10,20,10,10)) 

ggsave(plot4 , file=paste0(plots_dir, 'Violin_plot_Opentargets_GPS.jpg'), width = 6, height=6, dpi=300)


# Extended Fig 3
Dataset_genescores1= genescoredataset %>%arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100) 
cols=names(Dataset_genescores1[grepl('phenotype', colnames(Dataset_genescores1))])
# calculate no predictors that contribute to each score.change each predictor value with score to 1 
Dataset_genescores1$no.pred=rowSums(Dataset_genescores1[grepl('phenotype|genescores_', colnames(Dataset_genescores1))] != 0)  
Dataset_genescores1[grepl('phenotype|genescores_', colnames(Dataset_genescores1))][Dataset_genescores1[grepl('phenotype|genescores_', colnames(Dataset_genescores1))] !=0 ] <- 1
#Add bin groupings for GPS from 0 - 2.1
binsize =0.3
breaks1=seq(0,2.1,binsize)
labels1=paste(breaks1,'-', breaks1+binsize)
breaks1=c(breaks1, Inf)
labels1=gsub(paste0('2.1 - ', 2.1+binsize),'2.1+', labels1)
Dataset_genescores1_group=Dataset_genescores1 %>% filter(genescoresum!=0)
Dataset_genescores1_group$group <- cut(Dataset_genescores1_group$genescoresum, breaks = breaks1, labels = labels1, right = FALSE, include.lowest=F)
Dataset_genescores1_df=Dataset_genescores1_group %>% select(-gene, -drugname, -parentterm, -category, -mi, -CV, -order,-genescoresum, -percent,-no.pred, -group)
Dataset_genescores1_df1=cbind(group=Dataset_genescores1_group$group, no.pred=Dataset_genescores1_group$no.pred, percent= Dataset_genescores1_group$percent, Dataset_genescores1_df)
#format dataset, collapse by phenotype columns.
Dataset_genescores1_df2=data.frame(pivot_longer(data = Dataset_genescores1_df1, cols = -c(1:3), names_to = "Predictor", values_to = "weights"))
#calculate counts for number of features within each 0.3 bin 
Dataset_genescores1_df2_count1=data.frame(Dataset_genescores1_df2 %>% group_by(group,no.pred,percent) %>% filter(weights!=0) %>%  summarise(Predictor=paste0(Predictor, collapse=','))) 
Dataset_genescores1_df2_count2=data.frame(Dataset_genescores1_df2_count1 %>% group_by(group,no.pred,Predictor) %>%  tally() %>%   mutate(numbering = row_number()) %>%
separate_rows(Predictor,sep=',') %>% rename(value=n))
Feature_bygroup <- tibble(Predictor = rep(unique(Dataset_genescores1_df2_count2$Predictor),length(unique(as.character(Dataset_genescores1_df2_count2$group)))),
  group =  rep(unique(as.character(Dataset_genescores1_df2_count2$group)),length(unique(as.character(Dataset_genescores1_df2_count2$Predictor)))))
Dataset_genescores1_df2_count2_pred<-full_join(Feature_bygroup, Dataset_genescores1_df2_count2) %>% arrange(group,Predictor )
Dataset_genescores1_df2_count2$group<-as.character(Dataset_genescores1_df2_count2$group)
#Create plot
bp <- ggplot(Dataset_genescores1_df2_count2, aes(fill=Predictor, y=value, x=no.pred)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_hline(yintercept=1, color = "grey") +
  facet_wrap(~group,scales="free_y") +
  xlab('Number of genetic features contributing to the GPS') + ylab('Count') +labs(fill = "Feature") +
  scale_fill_manual(values=c('OMIM'='#DB72FB','HGMD'='#619CFF','ClinVar'=  '#00BA38', 'Gene Burden' = '#D39200',
      'Single Variant'=  '#00C19F', 'eQTL Phenotype'=  '#F8766D', 'Locus2gene'= '#00B9E3',  'pQTL Phenotype'=  '#FF61C3'))+
  theme_classic() 
  
ggsave(bp, file=paste0(plots_dir,'/Barplot_counts_genescoresum_bins_Opentargets.jpg'), width=7, height=6, dpi=300)


###Extended data Fig 5
binned_gps=fread(paste0(Revision_folder_final,'Increase_GPS_onegene_OT_Sider.txt'), data.table=F)
if(length(unique(binned_gps$or_sig[!is.na(binned_gps$or_sig)]))==1){
shape_vals=c(19)
} else{
  shape_vals=c(1,19)
}
max_uci <- round(ifelse(max(binned_gps$upperCI[!is.na(binned_gps$upperCI)])>max_value,max_value , max(binned_gps$upperCI[!is.na(binned_gps$upperCI)])+1))
#Fix labels for plot
rep_str = c('Opentargets_independent'='Open Targets','Sider'='SIDER')
binned_gps$datafile <- str_replace_all(binned_gps$datafile, rep_str)
#Create plot
plot1=ggplot(data=binned_gps, aes(x=genescoresum, y=OR)) + facet_wrap(~datafile)+  # Create panels for plot for 'Open Target and 'SIDER'
      geom_point(aes(color=OR,shape =or_sig),position=position_dodge(width=0.2),size=2.5, stroke=1.5)  +
      xlab('Genetic priority score bins') + ylab('Odds ratio (95% CI)') +
      scale_shape_manual(values=shape_vals) +guides(shape = "none") + 
      scale_x_discrete(limits = (levels(binned_gps$genescoresum))) +
      geom_linerange(aes(ymin=lowerCI, ymax=upperCI,color=OR),position=position_dodge(width=0.2)) +
      scale_color_gradient(low = ("blue"), high = ("red"), limits = c(0, 10),
        breaks = c(0, 2, 4, 6, 8, max_value),labels = c(0,2, 4, 6, 8, max_value), oob = scales::squish) +
      geom_hline(yintercept=1, color = "grey") +
      coord_cartesian(ylim=c(-0.5,max_uci), clip='off') +
      theme_classic() + 
      theme(axis.text.x=element_text(size=13,angle=45, vjust = 1, hjust=1,color = "black"), axis.text.y=element_text(color = "black",size=13),
        axis.title=element_text(size=15,color = "black"),legend.text=element_text(size=13), legend.title=element_text(size=13),
        axis.line = element_line(colour = 'black', size = 0.5), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        plot.margin=margin(10,20,10,10)) +  
        scale_y_continuous(limits=c(0,max_uci)) 

ggsave(plot1, file=paste0('Increase_GPS_onegene.pdf'), width = 7, height=7)

