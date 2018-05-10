#OGN
#########################
#preprocess: filter out proteins which are not in both ppi data and gene expression data;

library(pracma)
options(stringsAsFactors = FALSE);

#PPI network for biogrid 3.4.143
#only consider those proteins having gene expression data
dt <- read.delim('PPI_network.txt', sep='\t');
#names(dt)
intA=dt$interactorA;
intB=dt$interactorB;

pp <- strTrim(union(intA,intB));
#length(pp)

#gene expression data
gex <- read.table('geneexpressiondata.txt');
names(gex)
dim(gex)
length(unique(gex[,2]))

#overlap between proteins in network and in gene expression data
pp1 <- intersect(pp,strTrim(gex$V2));

#PPI network
net <- matrix(0,length(pp1),length(pp1));
for (i in 1:length(intA)){
  a=which(pp1==intA[i]);
  if (length(a)>0){
    b=which(pp1==intB[i]);
    if (length(b)>0){
      if (a!=b){
        net[a,b]=1;
        net[b,a]=1;
      }
    }	   
  }   
}

#test if there's protein which has no interaction
colNum2=vector();
for(i in 1:length(pp1)){
  colNum2[i]=sum(net[i,]);
}
min(colNum2)
max(colNum2)
a=which(colNum2>0);
netf=net[a,a];
ppf=pp1[a]; 
rownames(netf)=ppf;
colnames(netf)=ppf;

#find out gene expression for each protein in PPI network
gexd <- matrix(0,length(ppf),dim(gex)[2]-2);
N=0;
geneEXs <- strTrim(gex$V2);
for (i in 1:length(ppf)){
  a <- which(geneEXs==ppf[i]);
  #if (length(a)>0){
  N=N+1;
  if (length(a)>1){
    mx <- mean(as.matrix(gex[a[1],-c(1,2)]));
    ind <- 1;
    for (j in 2:length(a)){
      temp <- mean(as.matrix(gex[a[j],-c(1,2)]));
      if (mx < temp){
        mx <- temp;
        ind <- j;
      }
    }
    gexd[i,] <- as.numeric(gex[a[ind],-c(1,2)]);
  }else{
    gexd[i,] <- as.numeric(gex[a,-c(1,2)]);
  }
  # }
}

#essential proteins
ess <- read.table('combinedEssentialprotein.txt',header =TRUE);
names(ess)
essInd=vector();
for (i in 1:length(ppf)){
  a=which(ess$combinedEssentialProtein==ppf[i]);
  if (length(a)>0){
    essInd[i]=1;
  }else{
    essInd[i]=0;
  }
}

#correlation matrix
cordip <- matrix(0,length(ppf),length(ppf));
for (i in 1:length(ppf)-1){
  a=which(netf[i,]==1);
  if (length(a)>0){
    for (j in 1:length(a)){	    
      cordip[i,a[j]]=cor(gexd[i,],gexd[a[j],],use="p");
      cordip[a[j],i]=cordip[i,a[j]];
    } 
  }	     
}
###############
#calculate CoEWC
library('igraph')
gg <- simplify(graph_from_adjacency_matrix(netf,mode='undirected',weighted=NULL,diag=FALSE));

cco <- transitivity(gg,type='local',vids=V(gg),isolates='zero')
names(cco)=ppf

tpn=vector()
for(i in 1:length(ppf)){
  a=which(netf[i,]==1)
  if(length(a)>0){
    tpn[i]=sum(cordip[i,a]*cco[a])
  }else{
    tpn[i]=0
  }
}

#######################
#orthology
lm <- read.table(file='protein_allorthology.txt',header=TRUE)
names(lm)

orth=vector()
numOrgs=100
for (i in 1:length(ppf)){
  a=which(lm$protName==ppf[i])
  if(length(a)>0){
    orth[i]=lm$numOrthology[a]/numOrgs
  }else{
    orth[i]=0
  }
}
os=orth/max(orth)

#OGN
m1=seq(0,1,0.1)
ogn=matrix(0,length(ppf),length(m1))
for (i in 1:length(m1)){
  ogn[,i]=m1[i]*os+(1-m1[i])*(tpn/abs(max(tpn)))
}

#save to file
colnames(ogn) = paste(rep('OGN', length(m1)), m1, sep="");
info = data.frame(proteinName=ppf, ogn)
write.table(info,file="OGN.txt",row.names=FALSE,quote=FALSE)
