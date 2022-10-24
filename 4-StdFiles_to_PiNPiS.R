###################
#
# Input: files from 3-Make_gene_files_compute_stats.sh
# or reference fasta + gff and genotypes vcd
# output: PiN/PiS
# Kevin Sartori - 2022-09
#
###################

rm(list = ls())

library(vcfR)
library(RGenetics)
library(pegas)

fasta<-read.table(file = "fasta.fasta",sep = "\t",h=T,as.is = 1)
gff<-read.table("gff.gff")
vcf<-read.vcfR("vcf.recode.vcf")

# First step : make fasta in one row and then cut with the CDS sizes
names(fasta)<-unlist(strsplit(names(fasta),split = "\\."))[2]
for (i in 1:dim(fasta)[1]) { if (exists("pasta")) { pasta<-paste0(pasta,fasta[i,]) }else{pasta<-fasta[1,] } }
#nchar(pasta)

#dealing with strand sense
if (gff$V7[1]=="-") {
  seq_split<-unlist(strsplit(pasta,split = NULL))
  seq_comp<-rep(NA,length(seq_split))
  seq_comp[which(seq_split=="A")]<-"T"
  seq_comp[which(seq_split=="T")]<-"A"
  seq_comp[which(seq_split=="C")]<-"G"
  seq_comp[which(seq_split=="G")]<-"C"
  seq_comp[which(seq_split=="-")]<-"-"
  seq_comp_rev<-seq_comp[rev(1:length(seq_comp))]
  pasta<-paste0(seq_comp_rev,collapse = "")
  
  gff$V10<-gff$V4-min(gff$V4,gff$V5)+1
  gff$V11<-gff$V5-min(gff$V4,gff$V5)+1
  gff$length<-gff$V11-gff$V10
  
  new.fasta<-NULL
  for (i in 1:dim(gff)[1]) {
    new.fasta[i*2-1]<-paste0("CDS",i)
    new.fasta[i*2]<-substr(pasta,nchar(pasta)-gff$length[i],nchar(pasta))
    pasta<-substr(pasta,1,nchar(pasta)-gff$length[i]-1)
  
  }
}else{
  
  new.fasta<-NULL
  gff$V10<-gff$V4-min(gff$V4,gff$V5)+1
  gff$V11<-gff$V5-min(gff$V4,gff$V5)+1
  gff$length<-gff$V11-gff$V10
  
  for (i in 1:dim(gff)[1]) {
    new.fasta[i*2-1]<-paste0("CDS",i)
    new.fasta[i*2]<-substr(pasta,1,gff$length[i]+1)
    pasta<-substr(pasta,gff$length[i]+2,nchar(pasta))
    
  }
  
  }



#nchar(new.fasta)



##############



# gene CDS number 
CDS_nb<-dim(gff)[1]
# starting bp of the CDS
bp_start<-gff[,4]-1


# 1 - Pre-processing data
# initiating CDS list
CDS<-list(NA)
for(i in 1:CDS_nb){CDS[[i]]<-NA}
# initiating ALT list
ALT<-list(NA)
# table of alternative alleles
alt<-as.data.frame(vcf@fix[,c(2,4,5)],stringsAsFactors = F)
alt$POS<-as.numeric(alt$POS)
#table of genotypes
gt<-matrix(vcfR::extract.gt(vcf),nrow = dim(alt)[1])


for (k in 1:CDS_nb) {
  if (length(which(alt$POS %in% gff[k,4]:gff[k,5]))>0){
    # limits of CDS #k
    from<-which(alt$POS>gff[k,4])[1]
    to<-tail(which(alt$POS<gff[k,5]),n=1)
    # keep alt allele for this interval 
    alt_temp<-alt[from:to,]
    # Number of SNP (and indel) in this interval
    SNP<-dim(alt_temp)[1]
    # split alternative allele information in columns
    Y<-strsplit(alt_temp$ALT,split = ",")
    # initiate new columns in ALT
    X<-rep(NA,length(alt_temp$ALT))
    alt_temp<-data.frame(alt_temp,ALT1=X,ALT2=X,ALT3=X,ALT4=X,ALT5=X)
    # input new columns
    for (i in 1:length(alt_temp$ALT)) {
      for (j in 1:5) {
        alt_temp[i,3+j]<-Y[[i]][j]
      }
    }
    # remove the former "ALT" column 
    alt_temp<-alt_temp[,-3]
    # keep genotypes for the interval
    genotypes<-gt[from:to,]
    # Change missing values as reference values
    genotypes[which(is.na(genotypes))]<-"0/0"
    genotypes2<-genotypes
    # rotate matrix if only one SNP (weird R behavior)
    if (dim(alt_temp)[1]==1) { genotypes<-t(matrix(genotypes)) }
    if (dim(alt_temp)[1]==1) { genotypes2<-t(matrix(genotypes2)) }
    # Keep first parent's information 
    if(is.null(dim(genotypes)[2])==F){
      for (i in 1:dim(genotypes)[2]) {
        X<-unlist(strsplit(genotypes[,i],split = ""))
        X<-as.numeric(X[-seq(2,length(X),3)])
        X<-matrix(data = X,nrow = length(X)/2,ncol = 2,byrow = T,dimnames = NULL)
        for (j in 1:dim(genotypes)[1]) {
          genotypes[j,i]<-X[j,1]
        }
      }
    }else{
      for (i in 1:length(genotypes)) {
        X<-unlist(strsplit(genotypes[i],split = ""))
        X<-as.numeric(X[-seq(2,length(X),3)])
        X<-matrix(data = X,nrow = length(X)/2,ncol = 2,byrow = T,dimnames = NULL)
        genotypes[i]<-X[1]
        
      }
      
    }
    # Keep second parent's information 
    if(is.null(dim(genotypes2)[2])==F){
      for (i in 1:dim(genotypes2)[2]) {
        X<-unlist(strsplit(genotypes2[,i],split = ""))
        X<-as.numeric(X[-seq(2,length(X),3)])
        X<-matrix(data = X,nrow = length(X)/2,ncol = 2,byrow = T,dimnames = NULL)
        for (j in 1:dim(genotypes2)[1]) {
          genotypes2[j,i]<-X[j,2]
        }
      }
    }else{
      for (i in 1:length(genotypes2)) {
        X<-unlist(strsplit(genotypes2[i],split = ""))
        X<-as.numeric(X[-seq(2,length(X),3)])
        X<-matrix(data = X,nrow = length(X)/2,ncol = 2,byrow = T,dimnames = NULL)
        genotypes2[i]<-X[2]
        
      }
      
    }
    genotypes<-cbind(genotypes,genotypes2)
    # input list of genotype information for each CDS and each accession
    CDS[[k]]<-genotypes
    # input list of alternative alleles for each CDS
    ALT[[k]]<-alt_temp
  }
}

acc<-c(
  paste0(">",names(data.frame(vcfR::extract.gt(vcf))),"_A"),
  paste0(">",names(data.frame(vcfR::extract.gt(vcf))),"_B")
)

# create new empty fasta
last.fasta<-NULL  

# for each accession
for (j in 1:length(acc)) {
  # create new sequence from the reference that will receive the mutations
  alt_seq<-new.fasta
  # for each CDS
  for (i in 1:CDS_nb) {
    if(is.null(dim(CDS[[i]]))==F){
      # genotype of accession #j and CDS #i
      genotype<-as.numeric(CDS[[i]][,j])
      # first bp of CDS #i
      start<-bp_start[i]
      # for each SNP
      for (k in 1:length(genotype)) {
        # dealing with indel. Length of the reference locus
        l<-nchar(ALT[[i]]$REF[k])-1
        # The sequence of the CDS #i is located at line i*2 in the fasta
        # Next line change reference locus by accession locus according to the genotype file (REF, ALT1, ALT2 etc locus ?)
        if (ALT[[i]][k,2+genotype[k]]=="*") {
          substr(alt_seq[i*2],ALT[[i]]$POS[k]-start,ALT[[i]]$POS[k]-start+l)<-"-"
        }else{ substr(alt_seq[i*2],ALT[[i]]$POS[k]-start,ALT[[i]]$POS[k]-start+l) <- ALT[[i]][k,2+genotype[k]]}
        
      }
    }
    # dealing with strand sense. 
    if (gff$V7[i]=="-") {
      seq_split<-unlist(strsplit(alt_seq[i*2],split = NULL))
      seq_comp<-rep(NA,length(seq_split))
      seq_comp[which(seq_split=="A")]<-"T"
      seq_comp[which(seq_split=="T")]<-"A"
      seq_comp[which(seq_split=="C")]<-"G"
      seq_comp[which(seq_split=="G")]<-"C"
      seq_comp[which(seq_split=="-")]<-"-"
      seq_comp_rev<-seq_comp[rev(1:length(seq_comp))]
      alt_seq[i*2]<-paste0(seq_comp_rev,collapse = "")
      
    }
    
  }
  # save fasta for each
  last.fasta[j*2-1]<-acc[j]
  for (m in 1:CDS_nb) { if (exists("past.seq")) { past.seq<-paste0(past.seq,alt_seq[m*2]) }else{past.seq<-alt_seq[2] } }
  last.fasta[j*2]<-past.seq
  rm(past.seq)
}

#last.fasta<-as.data.frame(last.fasta)

write.table(x = last.fasta,file = "full_fasta.fasta",quote = F,row.names = F,col.names = F)


#########################
# compute pinpis
arg<-as.character(gff$V1[1])
rm(list = ls()[which(ls()!="arg")])

library(RGenetics)
gcode<-geneticCodeTable(DNA = TRUE)
gcode$fold_d<-NA
gcode$ST<-paste0(gcode$SecondPosition,gcode$ThirdPosition)
gcode$FT<-paste0(gcode$FirstPosition,gcode$ThirdPosition)
gcode$FS<-paste0(gcode$FirstPosition,gcode$SecondPosition)
conv<-t(rbind(c(1,2,3,4),c(4,"-","-",0)))

for (i in 1:dim(gcode)[1]) {
  nbaa1<-length(unique(gcode[which(gcode$ST==gcode$ST[i]),6]))
  nbaa2<-length(unique(gcode[which(gcode$FT==gcode$FT[i]),6]))
  nbaa3<-length(unique(gcode[which(gcode$FS==gcode$FS[i]),6]))
  nbaa1<-conv[conv[,1]==nbaa1,2]
  nbaa2<-conv[conv[,1]==nbaa2,2]
  nbaa3<-conv[conv[,1]==nbaa3,2]
  gcode$fold_d[i]<-paste0(nbaa1,nbaa2,nbaa3)
}
gcode<-gcode[,c(4,7)]
rm(conv,i,nbaa1,nbaa2,nbaa3)

# Open REF sequence
fasta<-read.table(file = "fasta.fasta",sep = "\t",h=T,as.is = 1)
names(fasta)<-unlist(strsplit(names(fasta),split = "\\."))[2]
for (i in 1:dim(fasta)[1]) { if (exists("pasta")) { pasta<-paste0(pasta,fasta[i,]) }else{pasta<-fasta[1,] } }
#nchar(pasta)

#split codons and convert to fold degenerate
pos<-seq(1,nchar(pasta),3)
for (i in 1:length(pos)) {
  a<-substr(pasta,pos[i],pos[i]+2)
  a<-gcode[gcode[,1]==a,2]
  if(!exists("fold")){fold<-a}else{fold<-paste0(fold,a)}
}
#make position lists of zero and four fold degenerate sites
zfold<-which(strsplit(fold,"",T)[[1]]=="0")
ffold<-which(strsplit(fold,"",T)[[1]]=="4")


# Open accessions sequences
zfasta<-read.FASTA("full_fasta.fasta")
ffasta<-zfasta
# split zero fold and four fold
for (i in 1:length(zfasta)) {
  zfasta[[i]]<-zfasta[[i]][zfold]
  ffasta[[i]]<-ffasta[[i]][ffold]
}

pin<-nuc.div(zfasta)
pis<-nuc.div(ffasta)
pinpis<-pin/pis

result<-read.table(paste0("../result_",arg,".txt"),h=T,sep=" ")
result[dim(result)[1]+1,]<-c(names(fasta),nchar(fasta),length(ffold),length(zfold),pis,pin,pinpis)
write.table(x = result,file = paste0("../result_",arg,".txt"),quote = F,row.names = F,col.names = T)

rm(list = ls())

