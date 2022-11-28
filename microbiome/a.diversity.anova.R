alpha_div <- function(physeq,method){
  #==check for validity of selected methods
  method<- match.arg(method,c("richness", "fisher", "simpson", "invsimpson", "shannon", "evenness","pd"), several.ok = TRUE)
  
  abund_table <- otu_table(physeq)
  df <- NULL
  if("richness"%in%method){
    R<- vegan::rarefy(abund_table,min(rowSums(abund_table)))
    df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))
    if(is.null(df)){
      df<-df_R}
    else {
      df<-rbind(df,df_R)}
  }
  if("fisher"%in%method){
    alpha <- vegan::fisher.alpha(abund_table)
    df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))
    if(is.null(df)){
      df<-df_alpha}
    else {
      df<-rbind(df,df_alpha)}
  }
  if("simpson"%in%method){
    simp <- vegan::diversity(abund_table, "simpson")
    df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))
    if(is.null(df)){
      df<-df_simp}
    else {
      df<-rbind(df,df_simp)}
  }
  if("invsimpson"%in%method){
    simp <- vegan::diversity(abund_table, "invsimpson")
    df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("InvSimpson",length(simp)))
    if(is.null(df)){
      df<-df_simp}
    else {
      df<-rbind(df,df_simp)}
  }
  if("shannon"%in%method){
    H<- vegan::diversity(abund_table)
    df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))
    if(is.null(df)){
      df<-df_H}
    else {
      df<-rbind(df,df_H)}
  }
  if("evenness"%in%method){
    H<-vegan::diversity(abund_table)
    S <- specnumber(abund_table)
    J <- H/log(S)
    df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))
    if(is.null(df)){
      df<-df_J}
    else {
      df<-rbind(df,df_J)}
  }
  if("pd"%in%method){
    otu_tree <- phyloseq::phy_tree(physeq)
    PD <- pd(abund_table, otu_tree ,include.root = TRUE)
    df_PD<-data.frame(sample=names(PD),value=PD,measure=rep("PD",length(PD)))
    if(is.null(df)){
      df<-df_PD}
    else {
      df<-rbind(df,df_PD)}
  }
  return(df)
}

perform_anova <- function(df,meta_table,grouping_column,pValueCutoff){
  
  dt<-data.table::data.table(data.frame(df,.group.=meta_table[,grouping_column]))
  #specifying a p-value cutoff for the ggplot2 strips
  pval<-dt[, list(pvalue = sprintf("%.2g",
                                   tryCatch(summary(aov(value ~ .group.))[[1]][["Pr(>F)"]][1],error=function(e) NULL))),
           by=list(measure)]
  #Filter out pvals that we are not significant
  pval<-pval[!pval$pvalue=="",]
  pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]
  
  #using sapply to generate significances for pval$pvalue using the cut function.
  pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, pValueCutoff, Inf),label=c("***", "**", "*", "")))})
  
  #Update df$measure to change the measure names if the grouping_column has more than three classes
  if(length(unique(as.character(meta_table[,grouping_column])))>2){
    df$measure<-as.character(df$measure)
    if(dim(pval)[1]>0){
      for(i in seq(1:dim(pval)[1])){
        df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
      }
    }
    df$measure<-as.factor(df$measure)
  }
  #Get all possible pairwise combination of values in the grouping_column
  s<-combn(unique(as.character(df[,grouping_column])),2)
  
  #df_pw will store the pair-wise p-values
  df_pw<-NULL
  for(k in unique(as.character(df$measure))){
    #We need to calculate the coordinate to draw pair-wise significance lines
    #for this we calculate bas as the maximum value
    bas<-max(df[(df$measure==k),"value"])
    #Calculate increments as 10% of the maximum values
    inc<-0.1*bas
    #Give an initial increment
    bas<-bas+inc
    for(l in 1:dim(s)[2]){
      #Do a pair-wise anova
      tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
      #Ignore if anova fails
      if(!is.na(as.numeric(tmp[length(tmp)]))){
        #Only retain those pairs where the p-values are significant
        if(as.numeric(tmp[length(tmp)])<pValueCutoff){
          if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}
          #Generate the next position
          bas<-bas+inc
        }
      }
    }
  }
  if(!is.null(df_pw)){
    df_pw<-data.frame(row.names=NULL,df_pw)
    names(df_pw)<-c("measure","from","to","y","p")
  }
  out <- list("df_pw"=df_pw, "df"=df)
  return(out)
}

plot_anova_diversity <- function(physeq, method, grouping_column,pValueCutoff=0.05)
{
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- otu_table(physeq)
  meta_table <- sample_data(physeq)
  
  #get diversity measure using selected methods
  div.df <- alpha_div(physeq, method)
  
  #=add grouping information to alpha diversity measures
  df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
  
  #perform anova of diversity measure between groups
  anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  
  #Draw the boxplots
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  p<-p+theme_bw()
  p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Samples")
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Groups")
  
  #This loop will generate the lines and signficances
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    for(i in 1:dim(df_pw)[1]){
      p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
    }
  }
  return(p)
}