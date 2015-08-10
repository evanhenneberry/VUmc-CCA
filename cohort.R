#Evan Henneberry
setwd("~/Documents/VUmc CCA/data")
library("openxlsx")
library("dplyr")
library("ggplot2")
library('reshape2')
library('plyr')
###in the data, there are 5 cohorts, so make a list of drug types, and patients in each drug type
#want to break the patients up into their drug groups
drugs <- c("das","sor","sun","vem","erl") #names of drugs

#patients per drug type
das <- c("Pt12.","Pt14.","Pt16.","Pt19.","Pt20.")
sor <- c("Pt1.","Pt11.","Pt18.","Pt2.","Pt21.","Pt33.","Pt9.")
sun <- c("Pt3.","Pt6.","Pt7.","Pt8.","Pt15.")
vem <- c("Pt13.","Pt22.","Pt27.","Pt28.","Pt39.","Pt41.","Pt42.")
erl <- c("Pt5.","Pt23.","Pt24.","Pt26.","Pt30.","Pt32.","Pt40.","Pt43.")

#new pt1 contains amino acid, position, sequence, protein and ppid
new.pt <- read.xlsx("PrOEF-150421-OPL1013-RdH-pTyrIP.ICK.Study(Aggregated.Data.Five.Drugs)-tp150519.xlsx", sheet=7, startRow=1, rowNames=F)

#from the original df, drop unimportant columns, change colnames, remove extra data
prep.df <- function(df){
  #large df, so drop some unimportant columns
  ptx <- df[,c(1,2,5,6,11,81:218)]
  #treatment 2 and 3, too low protein concentrations, so drop these
  ptx <- ptx[,!grepl("\\.Tx2|\\.Tx3|\\.PreTx3|\\.PreTx2", names(ptx))]
  
  names(ptx) <- gsub("\\(", "", names(ptx)) #"\\([^)]+\\)"
  names(ptx) <- gsub("\\)", "", names(ptx))
  names(ptx) <- gsub("Spectral.count.","",names(ptx))
  names(ptx) <- gsub(".[A-z]*.pTyrIP","",names(ptx))
  names(ptx) <- gsub("Normalised.intensity.","",names(ptx))
  names(ptx) <- gsub("Spectral.count..pTyrIP.","",names(ptx))
  names(ptx) <- gsub("Normalised.intensity..pTyrIP.","",names(ptx))
  
  names(ptx)[5] <- "position"
  ptx$Proteins <- gsub(";.*","",ptx$Proteins)
  ptx$position <- gsub(";.*","",ptx$position) 
  ptx$Gene.Names <- gsub(";.*","",ptx$Gene.Names)
  ptx <- ptx[,order(names(ptx))]
  
}

ptx <- prep.df(new.pt)

#now, need information from the other df
new.pts <- read.xlsx("PrOEF-150421-OPL1013-RdH-pTyrIP.ICK.Study(Aggregated.Data.Five.Drugs)-tp150519.xlsx", sheet=12, startRow=1, rowNames=F)

prep.df1 <- function(df1){
  df1 <- df1[,c(1,3,7,8,9)]
  df1$Phospho.STY.probabilities <- gsub("\\([^)]+\\)","",df1$Phospho.STY.probabilities)
  names(df1)[3] <- "pp.sequence"
  df1$Proteins <- gsub(";.*","",df1$Proteins)
  df1$Positions.within.proteins <- gsub(";.*","",df1$Positions.within.proteins) 
  
}

new.pts <- prep.df1(new.pts)
#so, need to get the count of each protein and position to make sure it is taken into account if
#the same substrate and phosphorylated position are found multiple times...networkin doesnt account
#for it
full.df <- merge(ptx, new.pts,by.x=c("Sequence","Proteins","position"),by.y=c("pp.sequence","Proteins","ppSTY.ID"))

#networkin does not consider if the same site is uploaded two times, so count the number of
#times each phosphosite occurs, and will match it back later (double check method, because is
#taken care of in the merge).
ptx.count <- full.df %>% group_by(Proteins, Positions.within.proteins) %>% tally()

#for networkin
toNI <- full.df[,c("Proteins","Positions.within.proteins","Amino.acid")] 
write.table(toNI, file="ptx1.NI.csv",sep="\t",quote=F,col.names=F,row.names=F)

##############downloaded NI output
ptx.pt1 <- read.delim("~/Documents/VUmc CCA/data/ptx1/NI.ptx1.1-500.tsv",header=T)
ptx.pt2 <- read.delim("~/Documents/VUmc CCA/data/ptx1/NI.ptx1.501-1000.tsv",header=T)
ptx.pt3 <- read.delim("~/Documents/VUmc CCA/data/ptx1/NI.ptx1.1001-1500.tsv",header=T)
ptx.pt4 <- read.delim("~/Documents/VUmc CCA/data/ptx1/NI.ptx1.1501-1709.tsv",header=T)

ptx.full <- rbind(ptx.pt1,ptx.pt2,ptx.pt3,ptx.pt4)
ptx.full <- data.frame(lapply(ptx.full, as.character), stringsAsFactors=F)
#changed to ptx.counts from ptx1.counts
ptx.full <- merge(ptx.full, ptx.count, by.x=c("X.substrate","position"),by.y=c("Proteins","Positions.within.proteins"))
ptx.full <- ptx.full[,c("X.substrate","position","id","netphorest_group","n")]
##########################################################################################
##################here, we decide to do sorafenib cohort####################################
#it is set up now to run on all other cohorts as well, just need to get the right patients
##only some patients have pre and post treatment samples, so use these only
patientsBoth.sor <- c("Pt1.","Pt11.","Pt18.","Pt21.","Pt33.")
ptList.both.sor <- list()

#for every patient in the sor cohort, make a separate data.frame in a list with the sequence, 
#position, protein and peptide id plus the counts for each patient
for(i in patientsBoth.sor){
  newDF <- names(select(ptx, contains(i)))
  newDF1 <- ptx[newDF]
  newDF1 <- cbind(ptx[,c("Sequence","position","Proteins","ppModPeptide.ID")],newDF1)
  newDF1 <- merge(newDF1,new.pts,by.x=c("Sequence","Proteins","position"),by.y=c("pp.sequence","Proteins","ppSTY.ID"))
  ptList.both.sor[[i]] <- assign(paste0("df.both.", i), data.frame(newDF1))
  rm(newDF,newDF1,i)
}

#we need the "position" column to ensure the merge call in the above loop worked properly, now it
#is not needed, so delete it
ptList.both.sor <- lapply(ptList.both.sor, function(x) { x["position"] <- NULL; x })

#next, we want to know the actual gene name not the uniprot code of the substrate that is
#phosphorylated by each kinase. this data is found in the original networkin download, so make
#a data.frame of every protein and its actual gene name.
gene2protein <- ptx[,c("Proteins","Gene.Names")]
#this will be used in the for loop below

#here, we are getting the unique substrates, kinases and peptides that are found before and after
#the treatment grouped by their netphorest (kinase) group and writing them to a separate data frame 
#(named "info.pre.x" and "info.post.x"), where we then convert the substrate code to gene name 
#finally, we get the number of uniques before and after treatment (numeric) and write to data.frame
# called "full.counts.x"
for(i in patientsBoth.sor){
  #first, take one patient from the list at a time and unlist as a data.frame
  ul.df <- do.call(cbind.data.frame,ptList.both.sor[[i]])
  noSum <- c("Sequence","Proteins","ppModPeptide.ID","Amino.acid","Positions.within.proteins")
  
  #need to remove rows with only zeros in intensity and count
  ul.df <- subset(ul.df, rowSums(ul.df[,(-which(names(ul.df) %in% noSum))])>0)
  
  ##want to get the unique peptides pre and post treatment  |
  df <- merge(ptx.full,ul.df,by.x=c("X.substrate","position"),by.y=c("Proteins","Positions.within.proteins"))
  
  #need to change colnames for pre and post, so get everything before the first period, and
  #after the last
  names(df)[8] <- "PreTx"
  names(df)[9] <- "PreTx1"
  names(df)[10] <- "Stx"
  names(df)[11] <- "Stx1"
  #here, sort the last three names so the order is amino acid, pretx and tx
  df <- df[,c(names(df)[1:7], sort(names(df)[8:12]))]
  vars.pre <- names(df)[9]
  pre.df <- names(select(df, contains(vars.pre)))
  pre.df1 <- df[pre.df]
  pre.df <- cbind(df[,c(1,3,4,6)],pre.df1)
  #remove only 0 rows
  pre.df <- subset(pre.df,rowSums(pre.df[,c(5,6)])>0)
  ##now get uniques first before treatment
  pre.df <- pre.df[,c(1:4)] #remove the counts/intensity because no care
  pre.df <- data.frame(lapply(pre.df, as.character), stringsAsFactors=F)
  
  #now, replace the X.substrate with the Gene name
  pre.df$X.substrate <- gene2protein[match(pre.df$X.substrate, gene2protein$Proteins),2]
  #now aggregate and then count the number of diff subs/pps per NI group
  pre.count <- pre.df %>% group_by(X.substrate,netphorest_group,Sequence) %>% tally()
  #write this out to a data.frame to use later
  pre.df[[i]] <- assign(paste0("info.pre.", i), data.frame(pre.df))
  
  #pre.count <- pre.df %>% group_by(X.substrate,netphorest_group,Sequence) %>% tally()
  pre.count <- as.data.frame(pre.count)
  pre.count1 <- rle( sort( pre.count$netphorest_group ) )
  pre.count$pre.tx <- pre.count1[[1]][ match( pre.count$netphorest_group , pre.count1[[2]] ) ]
  pre.count <- pre.count[,c("netphorest_group","pre.tx")]
  pre.count <- pre.count[!duplicated(pre.count), ]
  
  #now post treatment
  vars.post <- colnames(df)[11] ##need to make sure pre always come before tx
  post.df <- names(select(df,contains(vars.post)))
  #post.df <- names(select(df,one_of(vars.post)))
  post.df <- df[post.df]
  post.df <- cbind(df[,c(1,3,4,6)],post.df)
  #remove only 0 rows
  post.df <- subset(post.df,rowSums(post.df[,c(5,6)])>0)
  ##now get uniques
  post.df <- post.df[,c(1:4)] #remove the counts/intensity because no care
  post.df <- data.frame(lapply(post.df, as.character), stringsAsFactors=F)
  #now, replace the X.substrate with the Gene name
  post.df$X.substrate <- gene2protein[match(post.df$X.substrate, gene2protein$Proteins),2]
  post.count <- post.df %>% group_by(X.substrate,netphorest_group,Sequence) %>% tally()
  #write this out to a data.frame
  post.df[[i]] <- assign(paste0("info.post.", i), data.frame(post.df))
  #post.count <- post.df %>% group_by(X.substrate,netphorest_group,Sequence) %>% tally()
  post.count <- as.data.frame(post.count)
  post.count1 <- rle( sort( post.count$netphorest_group ) )
  post.count$post.tx <- post.count1[[1]][ match( post.count$netphorest_group , post.count1[[2]] ) ]
  post.count <- post.count[,c("netphorest_group","post.tx")]
  post.count <- post.count[!duplicated(post.count), ]
  #merge the pre and post treatment counts into one df.
  full.count <- merge(pre.count,post.count,by="netphorest_group",all=T)
  #for each patient, write out a data.frame containing the counts per netphorest group pre and post tx
  full.count[[i]] <- assign(paste0("full.counts.", i), data.frame(full.count))
  
  rm(i,full.count,ul.df)
}

#done with that, now the next step is to match each original patient df to the downloaded output
merge.pt <- function(NI.df,ptX.df){
  merge(NI.df, ptX.df, by.y=c("X.substrate","position"),by.x=c("Proteins","Positions.within.proteins"),all.x=F,all.y=F)
}

sor.ptx <- lapply(ptList.both.sor, merge.pt, ptx.full) ##will use this for fold changes

#drop these unwanted columns now, we do not need them.
toDrop <- c("Sequence","position","X.substrate","id","ppModPeptide.ID","Amino.acid","Proteins",
            "Positions.within.proteins")
doei <- function(df){
  df <- df[ , -which(names(df) %in% toDrop)]
}
sor.ptx.list <- lapply(sor.ptx,doei)

#next, aggregate all the spectral counts per netphorest group
agg.by.group <- function(df){
  df <- split(df,f=df[,"netphorest_group"])
  df1 <- lapply(df, function(x) aggregate(.~netphorest_group,x,sum))
  df1 <- data.frame(matrix(unlist(df1),nrow=length(df1),byrow=T),stringsAsFactors=F)
}

sor.ptx.list <- lapply(sor.ptx.list,agg.by.group)

#just in case patients without treatment have snuck in, remove them here (if any)
sor.ptx.list <- sor.ptx.list[sapply(sor.ptx.list, function(x) ncol(x) > 4)]

#next, change the column names to be more explicit
sor.ptx.list <- lapply(sor.ptx.list, function(x) {
  colnames(x) <- c("Kinase.group","before.int","before.cnt","post.int","post.cnt","n")
  return(x)
}) 

#convert the values to numeric, and remove "_group" from the netphorest names
sor.ptx.list <- lapply(sor.ptx.list, function(df){
  df$before.int <- as.numeric(df$before.int)
  df$post.int <- as.numeric(df$post.int)
  df$before.cnt <- as.numeric(df$before.cnt)
  df$post.cnt <- as.numeric(df$post.cnt)
  df$n <- as.numeric(df$n)
  df$Kinase.group <- gsub("_group","",df$Kinase.group)
  return(df)
})

####here is where we choose either spectral counts OR intensities (using spectral counts)
sor.ptx.list <- lapply(sor.ptx.list, function(df){
  zelim <- df[,c(1,3,5)]
  zelim <- subset(zelim, rowMeans(zelim[,c("before.cnt","post.cnt")])>0)
  zelim <- melt(zelim,id.vars="Kinase.group")
  #convert to factor, then split by "_" so plot nicer
  zelim$Kinase.group <- as.factor(zelim$Kinase.group)
  zelim$Kinase.group <- gsub("_", "\n", levels(zelim$Kinase.group))
  return(zelim)
})

##this works, trying above to get before above
cntPlot <- function(k,drugs.i){
  ggplot(data=melt.List[[k]], aes(x=reorder(Kinase.group,-value),y=value, fill = variable)) +
    geom_bar(stat='identity', position='dodge',color="black",width=0.6)+
    labs(y="Spectral Count per Kinase Group")+
    labs(x="")+
    ggtitle(paste0("Spectral counts of phosphopeptides in ",k, " data\n before and after ", drugs.i," treatment"))+
    scale_fill_manual("Treatment\nTime",labels = c("Before", "Two weeks"), 
                      values = c("white","black"))+
    theme(axis.text.y=element_text(face="bold"))
}

#plots for each drug group, still one at a time per drug group
for(k in names(sor.ptx.list)){
  myplot <- cntPlot(k,drugs[[2]]) #change the drugs[[i]]
  ggsave(myplot,filename=paste("sor.",k,".pdf",sep=""),height=6, width=11,dpi=72)
}

#############################################################################END

##so right after the for loop, we have sor.ptx, all the five sor patients before any processing.
#now, take this list and find fold changes

reduced.ptL1 <- Reduce(function(x, y) merge(x, y, all=TRUE), sor.ptx)

prepFC <- function(full.df){
  full.df <- full.df[,c(1,3,7:9,11,13,15,17,19,21,23,25,27)]
  #remove if zero counts
  full.df <- subset(full.df, rowSums(full.df[,c(5:14)])>0)
  full.df[full.df==0] <- NA
  vars <- colnames(full.df[,c(5:14)])
  full.df[vars] <- lapply(full.df[vars], log10)
  #now, aggregate
  full.df1 <- aggregate(.~Proteins+Sequence+n+netphorest_group, data=full.df, FUN=mean,na.action=na.pass, na.rm=TRUE)
  full.df1 <- replace(full.df1, is.na(full.df1), 0)
  full.df1$pt1.f <- 10^(full.df1[,6]-full.df1[,5])
  full.df1$pt11.f <- 10^(full.df1[,8]-full.df1[,7])
  full.df1$pt18.f <- 10^(full.df1[,10]-full.df1[,9])
  full.df1$pt21.f <- 10^(full.df1[,12]-full.df1[,11])
  full.df1$pt33.f <- 10^(full.df1[,14]-full.df1[,13])
  
  full.df1 <- full.df1[,c(1,2,4,15:19)]
  specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
  folds <- colnames(full.df1[,c(4:8)])
  #now, replace all very large or very low with +/-10,000
  full.df1[folds] <- as.data.frame(lapply(full.df1[folds], function(x){replace(x, x <1,-1/x)}))
  full.df1[folds] <- as.data.frame(lapply(full.df1[folds], function(x){replace(x, x >10000,100000)}))
  full.df1[folds] <- as.data.frame(lapply(full.df1[folds], function(x){replace(x, x <(-10000),(-100000))}))
  return(full.df1)
}

#now, calculate the fold changes
calculateFC <- function(prepared.df){
  prepared.df <- prepared.df[,c(1:8,4:8)]
  folds1 <- colnames(prepared.df[,c(9:13)])
  fcup <- prepared.df
  fcup[folds1] <- as.data.frame(lapply(fcup[folds1], function(x){replace(x, x >=2,2)}))
  fcup[folds1] <- as.data.frame(lapply(fcup[folds1], function(x){replace(x, x <2,0)}))
  
  fcup <- subset(fcup, rowSums(fcup[,c(9:13)])>5)
  fcup$fc <- rep("up",nrow(fcup))
  
  fcdown <- prepared.df
  fcdown[folds1] <- as.data.frame(lapply(fcdown[folds1], function(x){replace(x, x >=(-1.5),0)}))
  fcdown[folds1] <- as.data.frame(lapply(fcdown[folds1], function(x){replace(x, x <(-1.5),2)}))
  
  fcdown <- subset(fcdown, rowSums(fcdown[,c(9:13)])>5)
  fcdown$fc <- rep("down",nrow(fcdown))
  fullfold <- rbind(fcup,fcdown)
  return(fullfold)
  
}

full.df1 <- prepFC(reduced.ptL1)
fullfold <- calculateFC(full.df1)

##########################################end of working
