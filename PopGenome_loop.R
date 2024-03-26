library(getopt)
spec = matrix(c(
  'help' , 'h', 0, "logical",
  "vcf","v",1,"character",
  'chr',"c",1,"character",
  "group","g",1,"character",
  "num","n",2,"numeric",
  "width","w",2,"numeric",
  "jump","j",2,"numeric",
  "type","type",2,"numeric"
), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("
      Usage example: 
      Options: 
      --help            NULL                   help document
      --vcf             character              vcf.gz
      --chr             character              chrlist
      --group           character              group info
      --num             numeric                1000
      --width           numeric                100000  [30]
      --jump            numeric                10000  [10]
      --type            numeric                2  [1]

     Example:
     Rscript Pop_Genome.R --vcf Test.vcf.gz --chr chrlist --group groupinfo
\n")
  q(status=1);
}

if(!is.null(opt$help)|| length(opt)==0){Usage(spec)}
if(is.null(opt$vcf)){print_usage(spec)}
if(is.null(opt$chr)){print_usage(spec)}
if(is.null(opt$group)){print_usage(spec)}
if(is.null(opt$num)){opt$num = 1000}
if(is.null(opt$width)){opt$width = 100000}
if(is.null(opt$jump)){opt$jump = 10000}
if(is.null(opt$type)){opt$type = 2}

vcf = opt$vcf
chr = opt$chr
group = opt$group
num = opt$num
width = opt$width
jump = opt$jump
type = opt$type

popgenome <- function(
  vcf,
  chr,
  group,
  numcol = 1000,
  width = 100000,
  jump = 10000
){
  library(PopGenome)
  #read file
  group <- read.table(group, header = T, stringsAsFactors = F,sep = "\t")
  chrlist <- read.delim(chr,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
  for(i in 1:nrow(chrlist)){
    #read vcf 
    GENOME.class <- readVCF(filename = vcf, tid = chrlist[i,1], frompos = 1, topos = chrlist[i,2], numcols = numcol) 
    
    #Set the populations
    pop <- split(group,group[,2])
    for(j in 1:length(pop)){
      pop[[j]] <- pop[[j]][,1]
    }
    GENOME.class <- set.populations(object = GENOME.class, new.populations = pop, diploid = T)
    #Sliding window analyses
    GENOME.class <- sliding.window.transform(GENOME.class, width = width, jump = jump, type = type, whole.data=F)
  
    #get window information
    region <- as.data.frame(matrix(unlist(strsplit(GENOME.class@region.names, " ")), ncol = 7, byrow=T), stringsAsFactors=F)
    region$V1 <- chrlist[i,1]
    window <- region[,c(1,4,6)]
    names(window) <- c("Chr","Start","End")
    
    #neutrality test
    GENOME.class <- neutrality.stats(GENOME.class, detail = F,FAST = T)
    neu <- as.data.frame(get.neutrality(GENOME.class, stats = T)[1:length(pop)], check.names = F)
    names(neu) <- paste(names(neu), sort(rep(names(pop), 9)), sep = "_")
    
    #Fst and nuc diversity
    GENOME.class <- F_ST.stats(GENOME.class, mode = "nucleotide")
    Gpop <- as.data.frame(combn(names(pop),2),stringsAsFactors = F)
    Gpop <- paste0("Fst_",Gpop[1,],"_vs_",Gpop[2,])
    Fst <- as.data.frame(get.F_ST(object = GENOME.class,mode = "nucleotide",pairwise = T)[1])
    names(Fst) <- Gpop
    pi <- as.data.frame(GENOME.class@nuc.diversity.within)
    pi <- pi/width
    names(pi) <- paste0("Pi_",names(pop))
    
    #get window information
    region <- as.data.frame(GENOME.class@region.names)
    chr <- rep(chrlist[i,1],nrow(region))
    region <- cbind(chr,region)
    library(splitstackshape)
    window <- cSplit(indt = region,splitCols = names(region)[2],sep = " ")[,c(1,5,7)]
    names(window) <- c("Chr","Start","End")
    
    #combine all results
    segregating.sites <- neu[,paste0("n.segregating.sites_",names(pop))]
    TajimaD <- neu[,paste0("Tajima.D_",names(pop))]
    FuLiF <- neu[,paste0("Fu.Li.F_", names(pop))]
    FuLiD <- neu[,paste0("Fu.Li.D_", names(pop))]
    Result <- cbind(window,segregating.sites,pi,Fst,TajimaD,FuLiD,FuLiF)
    write.table(x = Result,file = paste0(chrlist[i,1],"_result.txt"),quote = F,sep = "\t",row.names = F,col.names = T)
  }
}

popgenome(vcf = vcf,chr = chr,group =  group,numcol = num,width = width,jump = jump)

