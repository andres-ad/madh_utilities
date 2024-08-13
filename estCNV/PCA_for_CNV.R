# Run est_CNV code through the Filter and reshape data section


# Modified function to get the fold changes for all targets
estCNV.forPCA <- function(data, sample.name="<sample>", verbose=F, k.gam=3,
                   plot.gam = F,
                   get.correction.factor = F,
                   target.info = NULL,
                   ...){
  
  # Make Pool and Reaction factors (Reaction no longer used as per formula below)
  data = data %>% 
    mutate(Pool.factor = as.factor(Pool),
           Reaction.factor = as.factor(Reaction))
  
  # Get a list of groups of targets to report fold changes on
  group.of.interest.list = unique(target.info$loci.of.interest$group)
  
  # Define formula (spline the same for all groups, with Groups and Pools as variables)
  form <- as.formula(paste("Reads ~ s(amplicon_length, k=k.gam) + Pool.factor +",
                           paste(group.of.interest.list,collapse=" + ")))
  
  # Get fit
  fit <- gam(form, data=data, family="poisson")
  
  
  if(get.correction.factor){  # This is run for controls (to get fold change in targets of interest to later normalize)
    
    # Make prediction frame to make plots with all data 
    pred.frame = expand_grid(amplicon_length=seq(floor(min(data$amplicon_length)),
                                                 ceiling(max(data$amplicon_length))),
                             Pool.factor = unique(data$Pool.factor))%>% 
      mutate(Reaction.factor = as.factor(substr(Pool.factor,1,1)))
    for (group.loci in group.of.interest.list) {
      pred.frame = pred.frame %>%
        mutate(!!group.loci := 0)
    }
    
    # Get predicted reads
    pred.frame$Reads = predict(fit,pred.frame,type="response")
    
    
    # Make a prediction frame for the targets we want residuals for
    targets.to.correct = data %>%
      select(Locus) %>% 
      mutate(Pool = sapply(str_split(Locus,"-"),tail,1),
             Reaction = substr(Pool,1,1)) %>% 
      mutate(Pool = ifelse(Pool == "1B","1Bor5",Pool)) %>% 
      left_join(amplicon.info %>% select(amplicon,amplicon_length),by = c("Locus"="amplicon")) 
    for (group.loci in unique(target.info$loci.of.interest$group)) {
      targets.to.correct = targets.to.correct %>%
        mutate(!!group.loci := 0) # When getting correction factor don't use Groups as variable 
    }
    
    
    targets.to.correct = targets.to.correct %>% 
      mutate(Pool.factor = as.factor(Pool),
             Reaction.factor = as.factor(Reaction))
    targets.to.correct$ExpectedReads =round(predict(fit, newdata = targets.to.correct,type="response"))
    
    # Make a data frame with the predicted values and get the observed fold change with respect to expected
    data.prediction = data %>%
      left_join(targets.to.correct %>% select(Locus,ExpectedReads),
                by = "Locus") %>% 
      mutate(FoldChange = Reads/ExpectedReads) %>% 
      select(SampleID,Group,Locus,Reads,ExpectedReads,FoldChange)
    
   
    
    return(data.prediction)
  }else{
    # If the function is used to get the fold change for each target group, return the coeficients
    coef <- summary(fit)$p.coef[group.of.interest.list]
    se <- summary(fit)$se[group.of.interest.list]
    
    return(list(Group = names(coef),FoldChangeGroup=exp(coef), Coef=coef, SE=se))
  }
}


fold.changes <- map_df(unique(paste(target.counts.filtered$SampleID, target.counts.filtered$Run, sep = ":")), 
                       function(sample_run) {
  parts <- strsplit(sample_run, ":")[[1]]
  sample.name <- parts[1]
  run.name <- parts[2]
  data <- filter(target.counts.filtered, SampleID == sample.name, Run == run.name)
  output <- estCNV.forPCA(data, get.correction.factor = T,target.info = target.info,plot.data=F)
  output %>% 
    mutate(Run = run.name)
})


# Filter diversity targets 
fc.1A = fold.changes %>% 
  filter(str_detect(Locus,"1A"))

# make wide data frame to prep for matrix
fc.1A.wide = fc.1A %>% select(SampleID,Locus,FoldChange) %>% 
  distinct() %>% 
  pivot_wider(names_from = Locus,values_from = FoldChange,values_fill=NA)

# make a matrix out of the wide data frame
fc.1A.matrix = as.matrix(fc.1A.wide[,2:ncol(fc.1A.wide)])
rownames(fc.1A.matrix)=fc.1A.wide[,1]$SampleID

# run pca
pc = prcomp(fc.1A.matrix,scale. = T)

# make pca results a dataframe
pc.df = data.frame(SampleID = rownames(pc$x), pc1 = pc$x[,1],pc2 = pc$x[,2]) %>% 
  left_join(samples.filtered)

ggplot(pc.df)+
  geom_point(aes(x = pc1,y=pc2,color=as.factor(paste(Run,Batch))))+
  facet_wrap(~SuperBatch,ncol=1)



## CONTROLS


# Filter diversity targets 
fc.1A = fold.changes %>% 
  left_join(samples.filtered,by=c("SampleID","Run")) %>% 
  filter(str_detect(Locus,"1A"),CNVControl)

# make wide data frame to prep for matrix
fc.1A.wide = fc.1A %>% select(SampleID,Locus,FoldChange) %>% 
  distinct() %>% 
  pivot_wider(names_from = Locus,values_from = FoldChange,values_fill=NA)

# make a matrix out of the wide data frame
fc.1A.matrix = as.matrix(fc.1A.wide[,2:ncol(fc.1A.wide)])
rownames(fc.1A.matrix)=fc.1A.wide[,1]$SampleID

# run pca
pc = prcomp(fc.1A.matrix,scale. = T)

# make pca results a dataframe
pc.df = data.frame(SampleID = rownames(pc$x), pc1 = pc$x[,1],pc2 = pc$x[,2]) %>% 
  left_join(samples.filtered)

ggplot(pc.df)+
  geom_point(aes(x = pc1,y=pc2,color=SuperBatch))


