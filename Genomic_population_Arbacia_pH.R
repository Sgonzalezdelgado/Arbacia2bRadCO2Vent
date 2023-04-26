######## Arbacia 2bRad CO2vent analysis ########

### BENJAMINI AND YEKUTELY ("BY") CORRECTION FOR P-VALUES
library(readr)
library(FSA)

data <- read_delim("input.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
data

#order data by p-value, "P" is the column with p-values
data=data[order(data$P),]

#check if data is ordered
headtail(data) 

#Perform p-value adjustments and add to data frame
data$FDR<-p.adjust(data$P, method = "BY")
data$FDR
write.csv(data,file="Output_adjust.csv",quote=FALSE,sep = ";")


### PLOT FIGURES FOR THE SOFTWARE BAYESCAN IN R 
#function
plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}

plot_bayescan("Arbaciabayescan.txt",0,FDR=0.05)
results<-plot_bayescan("Arbaciabayescan.txt",0,FDR=0.05)

#show which loci are the outliers
results$outliers
#number of outliers
results$nb_outliers
#shows the different FST values
mydata<-read.table("Arbaciabayescan.sel",colClasses="numeric")
mydata


### RDA ANALYSIS

#we use the library vcfR to convert the vcf into the OutFLANK format
library(OutFLANK)
library(qvalue)
library(vcfR)
if (!("OutFLANK" %in% installed.packages())){install_github("whitlock/OutFLANK")}
library(ggplot2)
library(dartR)
library(vegan)

obj.vcfR <- read.vcfR("Arbacia.vcf")

#extract useful informations about snp id and position
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP

#gather this info in a dataframe
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) # save info about id, chr, position

#transform into a numeric
chr_pos$position<-as.numeric(as.character(chr_pos$position)) 

#we expect that it will be useful for subsequent analysis to have a file with snp id and position so let's write it in our folder 02_data
write.table(chr_pos, "ArbaciaSNPpos.txt", sep="\t", quote=F, row.names=TRUE)

#extract and format genotype matrix
geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes

#an empty matrix, (9 stands for missing data)
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno))

#that we fill with genotypes
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

#an overview of our data and its first 10 rows/10 columns
table(as.vector(G))
dim(G)
G[1:10,1:10]

#as it will be useful later, I suggest that we export it now
write.table(G, "geno_matrix.txt", sep="\t", col.names=TRUE, row.names=TRUE)

#import pop info
info_samples_mg<- read.table("AMBIENT.csv", header=T, fileEncoding="utf-16")
head(info_samples_mg)
pop_vector<- info_samples_mg$pop

# FST matrix with OutFLANK
my_fst <- utils.outflank.MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector)
??MakeDiploidFSTMat

geno<-read.table("geno_matrix.txt")

SNP_pos<-read.table("ArbaciaSNPpos.txt", header=T)

#transpose data and give meaningful colnames
gen<-t(geno)
colnames(gen)<-paste(SNP_pos$chromosome, SNP_pos$position, sep="_")
gen [1:10,1:10]

#replace 9 by NA
gen[which(gen=="9")]<- NA

#evaluate % of missing
sum(is.na(gen))/(dim(gen)[1]*dim(gen)[2]) # about 1.5% of missing data

#impute missing with the most common geno
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
gen.imp
sum(is.na(gen.imp)) # No NAs

info<-read.table("AMBIENT.csv", header=T, fileEncoding="utf-16")
head(info)

pH.rda <- rda(gen.imp ~ info$pH, scale=T)
pH.rda

RsquareAdj(pH.rda)

pH.signif.full <- anova.cca(pH.rda, parallel=getOption("mc.cores")) # default is permutation=999
pH.signif.full

##plots of RDA data
pdf("rda1_pH2.pdf")
plot(pH.rda, scaling=3)
info$pop=factor(info$pop, levels= c("C","T1","T2","V1", "V2"))
points(pH.rda, display="sites", pch=20, cex=1.3, col=info$pop, scaling=3) #colors from:https://myrbooksp.netlify.app/graph2.html
levels(info$pop) <- c("C","T1","T2","V1", "V2")
legend("bottomright", legend=levels(eco), bty="n", col = 1:5, pch=20, cex=1) #colors from:https://myrbooksp.netlify.app/graph2.html
dev.off()

load.pH.rda <- scores(pH.rda, choices=c(1), display="species") 

#load info about snps
SNP_pos<-read.table("ArbaciaSNPpos.txt", header=T)
load.pH.rda.pos<-cbind(SNP_pos,load.pH.rda )
head(load.pH.rda.pos)

#plot distribution
jpeg("06_rda/rda1_loading_hist.jpeg")
hist(load.pH.rda.pos$RDA1, main="Loadings on RDA1")
dev.off()

#chose the sd limit (Z values could be: 2-p=0.05; 2.25-p=0.025; 2.5-p=0.01)
z=2

lim_min<- mean(load.pH.rda.pos$RDA1) - ( z * sd(load.pH.rda.pos$RDA1) )
lim_max<- mean(load.pH.rda.pos$RDA1) + ( z * sd(load.pH.rda.pos$RDA1) )

#outliers
outlier_pH <- load.pH.rda.pos[load.pH.rda.pos$RDA1 >=lim_max | load.pH.rda.pos$RDA1 <= lim_min ,]
outlier_pH
nrow(outlier_pH)

#export them
write.table(outlier_pH, "outlier_pH_rda.txt", row.names=F, quote=F, sep="\t")


### GET OUTLIERS DETECTED BY TWO DIFFERENT METHODS (e.g., RADA and Bayescan)

library(dplyr)
library(ggVennDiagram)

#load outliers tables
outlier_ph_rda<-read.table("outlier_pH_rda.txt", header=T)
head(outlier_ph_rda)
nRDA<-dim(outlier_ph_rda)[1]
nRDA #how many outliers?

outlier_fst<-read.table("Fstoutliers.txt", header =T,fileEncoding="utf-16" )
head(outlier_fst)
nBP<-dim(outlier_fst)[1]
nBP #how many outliers?

#join outliers keeping positions present in either the 1st or the 2nd database (or both)
outlier_fulljoin<-full_join(outlier_ph_rda,outlier_fst)
head(outlier_fulljoin)
nALL<-dim(outlier_fulljoin)[1]
nALL # how many in total?

#export them
outlier_fulljoin2 <- select(outlier_fulljoin, id_snp,chromosome, position)
write.table(outlier_fulljoin2, "outlier_pH_all.txt", row.names=F, quote=F, sep="\t")

#join outliers keeping positions present in either the 1st or the 2nd database (or both)
outlier_ph_innerjoin<-inner_join(outlier_ph_rda,outlier_fst)
head(outlier_ph_innerjoin)
dim(outlier_ph_innerjoin)
nboth<-dim(outlier_ph_innerjoin)[1]
nboth #how many joint outliers?

#visualize
ggVennDiagram(list(rda = 1:nRDA, BP = (nRDA+1-nboth):(nRDA-nboth+nBP)))


#### SELECT THE OUTLIERS TO CREATE A NEW DATASET

#load libraries
library(genepop)
if (!require("devtools")) install.packages("devtools") # to install
#install the package from *Github*
devtools::install_github("rystanley/genepopedit") 
library(genepopedit) # load the library
library(data.table)

#Import data base and define the output directory
file.link<-"C:~Arbacia.txt"
GenePop <- data.table::fread(file.link,header=FALSE, sep="\t",stringsAsFactors=FALSE)
output_dir <- "C:~Newdataset"

#Define the variables
PopNames <- genepop_detective(GenePop, variable="Pops")
PopCounts <- genepop_detective(GenePop, variable="PopNum")

#inspect the names of the loci within the Genepop data
genepop_detective(GenePop,"Loci")

#subset the Genepop file and 'keep' the specified loci names, creating a new database 
#get outliers RDA
subloci <-c("chr1_102295","chr1_109465","chr1_11032","chr1_112372","chr1_114274","chr1_114797","chr1_117592","chr1_126500","chr1_126810","chr1_12761","chr1_13255","chr1_135593","chr1_1361","chr1_137186","chr1_139944","chr1_143161","chr1_144462","chr1_145020","chr1_147571","chr1_149447","chr1_150704","chr1_151032","chr1_152465","chr1_153244","chr1_154553","chr1_163485","chr1_167984","chr1_170872","chr1_174193","chr1_17703","chr1_179036","chr1_182829","chr1_189737","chr1_190812","chr1_19171","chr1_196458","chr1_198488","chr1_208744","chr1_211194","chr1_215193","chr1_218530","chr1_219662","chr1_223589","chr1_225191","chr1_22699","chr1_227280","chr1_241724","chr1_247360","chr1_252376","chr1_255944","chr1_262092","chr1_26741","chr1_269503","chr1_270079","chr1_276867","chr1_280769","chr1_28350","chr1_287234","chr1_296723","chr1_298624","chr1_299180","chr1_30035","chr1_303029","chr1_305585","chr1_315889","chr1_324019","chr1_3261","chr1_329812","chr1_33343","chr1_334593","chr1_33779","chr1_341408","chr1_348505","chr1_349470","chr1_356956","chr1_357102","chr1_358793","chr1_359759","chr1_370069","chr1_371942","chr1_382360","chr1_386688","chr1_38708","chr1_387198","chr1_390067","chr1_39087","chr1_394062","chr1_394283","chr1_394695","chr1_400662","chr1_400793","chr1_403751","chr1_405316","chr1_412784","chr1_4130","chr1_414207","chr1_415173","chr1_426943","chr1_429424","chr1_431810","chr1_433854","chr1_436376","chr1_45232","chr1_454876","chr1_456858","chr1_4625","chr1_463295","chr1_474073","chr1_474453","chr1_474596","chr1_475803","chr1_479617","chr1_48123","chr1_482338","chr1_482461","chr1_49079","chr1_503299","chr1_510238","chr1_511948","chr1_513674","chr1_516621","chr1_51886","chr1_524926","chr1_525401","chr1_529952","chr1_542719","chr1_54926","chr1_553150","chr1_556461","chr1_567056","chr1_569611","chr1_574536","chr1_578277","chr1_578639","chr1_58346","chr1_5916","chr1_592945","chr1_59308","chr1_59707","chr1_603833","chr1_6102","chr1_610470","chr1_613573","chr1_619245","chr1_619483","chr1_622636","chr1_624373","chr1_625111","chr1_630304","chr1_632296","chr1_636455","chr1_642158","chr1_643345","chr1_647326","chr1_654319","chr1_656795","chr1_65756","chr1_661364","chr1_669913","chr1_671992","chr1_673676","chr1_679482","chr1_680535","chr1_683579","chr1_685088","chr1_68684","chr1_687161","chr1_687769","chr1_690424","chr1_695898","chr1_69812","chr1_71668","chr1_73109","chr1_77678","chr1_78138","chr1_80265","chr1_82035","chr1_8220","chr1_85697","chr1_86201","chr1_88349","chr2_102921","chr2_109213","chr2_110063","chr2_110105","chr2_111502","chr2_111789","chr2_112195","chr2_112518","chr2_116377","chr2_11792","chr2_121413","chr2_122164","chr2_122994","chr2_125574","chr2_126064","chr2_126273","chr2_126704","chr2_12688","chr2_130373","chr2_130691","chr2_132029","chr2_13289","chr2_13808","chr2_139949","chr2_149736","chr2_150614","chr2_151585","chr2_152147","chr2_155224","chr2_156699","chr2_158167","chr2_165646","chr2_167309","chr2_17144","chr2_174187","chr2_174539","chr2_176298",
            "chr2_176987","chr2_183266","chr2_190529","chr2_19099","chr2_195175","chr2_196891","chr2_20136","chr2_205139","chr2_206119","chr2_213632","chr2_214606","chr2_230525","chr2_232318","chr2_241860","chr2_242624","chr2_253792","chr2_25421","chr2_261390","chr2_263023","chr2_275021","chr2_278084","chr2_278408","chr2_279905","chr2_280028","chr2_281572","chr2_281909","chr2_291338","chr2_293133","chr2_297838","chr2_299037","chr2_300368","chr2_305724","chr2_31224","chr2_318284","chr2_32286","chr2_328376","chr2_329643","chr2_332442","chr2_33346","chr2_337051","chr2_35577","chr2_363748","chr2_367947","chr2_371181","chr2_374750","chr2_383417","chr2_385288","chr2_39321","chr2_394373","chr2_396054","chr2_397028","chr2_40022","chr2_40219","chr2_407081","chr2_407320","chr2_40800","chr2_413265","chr2_415245","chr2_415861","chr2_419753","chr2_423345","chr2_42464","chr2_428272","chr2_429247","chr2_434999","chr2_445887","chr2_450989","chr2_451738","chr2_452388","chr2_455246","chr2_46024","chr2_469992","chr2_473525","chr2_47854","chr2_479899","chr2_48285","chr2_484663","chr2_493482","chr2_50167","chr2_508921","chr2_513497","chr2_520014","chr2_524295","chr2_52937","chr2_5321","chr2_547767","chr2_551402","chr2_570437","chr2_572086","chr2_573954","chr2_57734","chr2_579880","chr2_591239","chr2_59680","chr2_598347","chr2_599238","chr2_600749","chr2_601062","chr2_621055","chr2_625691","chr2_625799","chr2_626781","chr2_63091","chr2_63441","chr2_635427","chr2_643410","chr2_675728","chr2_677141","chr2_679115","chr2_680120","chr2_690499","chr2_73183","chr2_7544","chr2_81665","chr2_82325","chr2_84402","chr2_88257","chr2_91612","chr2_92865","chr2_94337","chr2_94843","chr2_97145","chr3_109577","chr3_115256","chr3_118782","chr3_136313","chr3_136781","chr3_141816","chr3_146864","chr3_169830","chr3_173744","chr3_175097","chr3_175890","chr3_208675","chr3_209424","chr3_220546","chr3_227784","chr3_236861","chr3_241153","chr3_245537","chr3_248340","chr3_254756","chr3_259941","chr3_275417","chr3_29072","chr3_29899","chr3_304697","chr3_320310","chr3_323585","chr3_326854","chr3_329016","chr3_347021","chr3_355949","chr3_372917","chr3_373188","chr3_41789","chr3_430124","chr3_446373","chr3_450668","chr3_470840","chr3_472098","chr3_489354","chr3_502788","chr3_505170","chr3_517193","chr3_518030","chr3_53312","chr3_547002","chr3_550433","chr3_55578","chr3_594556","chr3_606044","chr3_61067","chr3_612729","chr3_629507","chr3_658699","chr3_66786","chr3_675025","chr3_683661","chr3_686511","chr3_696184","chr3_83467","chr3_83577","chr3_89514","chr3_94554","chr3_99713","chr4_103222","chr4_119297","chr4_132199","chr4_151104","chr4_174788","chr4_186008","chr4_273572","chr4_31044","chr4_3400","chr4_380358","chr4_46035","chr4_540423","chr4_575431","chr4_647866","chr4_667566","chr4_79362","chr4_8444","chr4_94953","chr5_146440","chr5_251109","chr5_382543","chr5_543976","chr5_672831","chr5_83357","chr6_306273","chr7_120807","chr8_252191","chr8_286824")
subset_genepop(genepop= GenePop, keep = TRUE, subs = subloci, path = paste0(output_dir,"Genepop_RDA.txt"))

#remove all outliers and get the Neutral dataset
subloci <-c("chr1_102295","chr1_109465","chr1_11032","chr1_112372","chr1_114274","chr1_114797","chr1_117592","chr1_126500","chr1_126810","chr1_12761","chr1_13255","chr1_135593","chr1_1361","chr1_137186","chr1_139944","chr1_143161","chr1_144462","chr1_145020","chr1_147571","chr1_149447","chr1_150704","chr1_151032","chr1_152465","chr1_153244","chr1_154553","chr1_163485","chr1_167984","chr1_170872","chr1_174193","chr1_17703","chr1_179036","chr1_182829","chr1_189737","chr1_190812","chr1_19171","chr1_196458","chr1_198488","chr1_208744","chr1_211194","chr1_215193","chr1_218530","chr1_219662","chr1_223589","chr1_225191","chr1_22699","chr1_227280","chr1_241724","chr1_247360","chr1_252376","chr1_255944","chr1_262092","chr1_26741","chr1_269503","chr1_270079","chr1_276867","chr1_280769","chr1_28350","chr1_287234","chr1_296723","chr1_298624","chr1_299180","chr1_30035","chr1_303029","chr1_305585","chr1_315889","chr1_324019","chr1_3261","chr1_329812","chr1_33343","chr1_334593","chr1_33779","chr1_341408","chr1_348505","chr1_349470","chr1_356956","chr1_357102","chr1_358793","chr1_359759","chr1_370069","chr1_371942","chr1_382360","chr1_386688","chr1_38708","chr1_387198","chr1_390067","chr1_39087","chr1_394062","chr1_394283","chr1_394695","chr1_400662","chr1_400793","chr1_403751","chr1_405316","chr1_412784","chr1_4130","chr1_414207","chr1_415173","chr1_426943","chr1_429424","chr1_431810","chr1_433854","chr1_436376","chr1_45232","chr1_454876","chr1_456858","chr1_4625","chr1_463295","chr1_474073","chr1_474453","chr1_474596","chr1_475803","chr1_479617","chr1_48123","chr1_482338","chr1_482461","chr1_49079","chr1_503299","chr1_510238","chr1_511948","chr1_513674","chr1_516621","chr1_51886","chr1_524926","chr1_525401","chr1_529952","chr1_542719","chr1_54926","chr1_553150","chr1_556461","chr1_567056","chr1_569611","chr1_574536","chr1_578277","chr1_578639","chr1_58346","chr1_5916","chr1_592945","chr1_59308","chr1_59707","chr1_603833","chr1_6102","chr1_610470","chr1_613573","chr1_619245","chr1_619483","chr1_622636","chr1_624373","chr1_625111","chr1_630304","chr1_632296","chr1_636455","chr1_642158","chr1_643345","chr1_647326","chr1_654319","chr1_656795","chr1_65756","chr1_661364","chr1_669913","chr1_671992","chr1_673676","chr1_679482","chr1_680535","chr1_683579","chr1_685088","chr1_68684","chr1_687161","chr1_687769","chr1_690424","chr1_695898","chr1_69812","chr1_71668","chr1_73109","chr1_77678","chr1_78138","chr1_80265","chr1_82035","chr1_8220","chr1_85697","chr1_86201","chr1_88349","chr2_102921","chr2_109213","chr2_110063","chr2_110105","chr2_111502","chr2_111789","chr2_112195","chr2_112518","chr2_116377","chr2_11792","chr2_121413","chr2_122164","chr2_122994","chr2_125574","chr2_126064","chr2_126273","chr2_126704","chr2_12688","chr2_130373","chr2_130691","chr2_132029","chr2_13289","chr2_13808","chr2_139949","chr2_149736","chr2_150614","chr2_151585","chr2_152147","chr2_155224","chr2_156699","chr2_158167","chr2_165646","chr2_167309","chr2_17144","chr2_174187","chr2_174539","chr2_176298",
            "chr2_176987","chr2_183266","chr2_190529","chr2_19099","chr2_195175","chr2_196891","chr2_20136","chr2_205139","chr2_206119","chr2_213632","chr2_214606","chr2_230525","chr2_232318","chr2_241860","chr2_242624","chr2_253792","chr2_25421","chr2_261390","chr2_263023","chr2_275021","chr2_278084","chr2_278408","chr2_279905","chr2_280028","chr2_281572","chr2_281909","chr2_291338","chr2_293133","chr2_297838","chr2_299037","chr2_300368","chr2_305724","chr2_31224","chr2_318284","chr2_32286","chr2_328376","chr2_329643","chr2_332442","chr2_33346","chr2_337051","chr2_35577","chr2_363748","chr2_367947","chr2_371181","chr2_374750","chr2_383417","chr2_385288","chr2_39321","chr2_394373","chr2_396054","chr2_397028","chr2_40022","chr2_40219","chr2_407081","chr2_407320","chr2_40800","chr2_413265","chr2_415245","chr2_415861","chr2_419753","chr2_423345","chr2_42464","chr2_428272","chr2_429247","chr2_434999","chr2_445887","chr2_450989","chr2_451738","chr2_452388","chr2_455246","chr2_46024","chr2_469992","chr2_473525","chr2_47854","chr2_479899","chr2_48285","chr2_484663","chr2_493482","chr2_50167","chr2_508921","chr2_513497","chr2_520014","chr2_524295","chr2_52937","chr2_5321","chr2_547767","chr2_551402","chr2_570437","chr2_572086","chr2_573954","chr2_57734","chr2_579880","chr2_591239","chr2_59680","chr2_598347","chr2_599238","chr2_600749","chr2_601062","chr2_621055","chr2_625691","chr2_625799","chr2_626781","chr2_63091","chr2_63441","chr2_635427","chr2_643410","chr2_675728","chr2_677141","chr2_679115","chr2_680120","chr2_690499","chr2_73183","chr2_7544","chr2_81665","chr2_82325","chr2_84402","chr2_88257","chr2_91612","chr2_92865","chr2_94337","chr2_94843","chr2_97145","chr3_109577","chr3_115256","chr3_118782","chr3_136313","chr3_136781","chr3_141816","chr3_146864","chr3_169830","chr3_173744","chr3_175097","chr3_175890","chr3_208675","chr3_209424","chr3_220546","chr3_227784","chr3_236861","chr3_241153","chr3_245537","chr3_248340","chr3_254756","chr3_259941","chr3_275417","chr3_29072","chr3_29899","chr3_304697","chr3_320310","chr3_323585","chr3_326854","chr3_329016","chr3_347021","chr3_355949","chr3_372917","chr3_373188","chr3_41789","chr3_430124","chr3_446373","chr3_450668","chr3_470840","chr3_472098","chr3_489354","chr3_502788","chr3_505170","chr3_517193","chr3_518030","chr3_53312","chr3_547002","chr3_550433","chr3_55578","chr3_594556","chr3_606044","chr3_61067","chr3_612729","chr3_629507","chr3_658699","chr3_66786","chr3_675025","chr3_683661","chr3_686511","chr3_696184","chr3_83467","chr3_83577","chr3_89514","chr3_94554","chr3_99713","chr4_103222","chr4_119297","chr4_132199","chr4_151104","chr4_174788","chr4_186008","chr4_273572","chr4_31044","chr4_3400","chr4_380358","chr4_46035","chr4_540423","chr4_575431","chr4_647866","chr4_667566","chr4_79362","chr4_8444","chr4_94953","chr5_146440","chr5_251109","chr5_382543","chr5_543976","chr5_672831","chr5_83357","chr6_306273","chr7_120807","chr8_252191","chr8_286824",
            "chr2_545815","chr2_408805","chr1_437311","chr2_245946","chr2_131536","chr2_188364","chr2_131024","chr2_640898","chr2_112733","chr2_549053","chr2_630124","chr2_373965","chr1_679373")
subset_genepop(genepop= GenePop, keep = FALSE, subs = subloci, path = paste0(output_dir,"Genpop_neutral.txt"))
?subset_genepop
#subset the Genepop file and 'keep' the specified loci names.
#L<-genepop_detective(GenePop,"Loci")


### HEATMAP WITH GGPLOT2 FOR FST VALUES AND ALLELIC FREQUENCY

# Load required packages
library(ggplot2)
library(RColorBrewer)

data1 <-Matrix_FST_RDA2
data2 <-allec_freq_gene

#find your coloursç
colfunc<- colorRampPalette(c("white","blue"))
colfunc(20)
colfunc<- colorRampPalette(c("white","purple"))
colfunc(20)

#FST values
plot1<-ggplot(data, aes(Y, X, fill= FST_RDA), hc.order = TRUE) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "white", high = "blue", mid = "#A1A1FF", 
                       midpoint = 0.002, limit = c(0,0.004), space = "Lab", 
                       name="FST") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1), axis.text.y = element_text(
                                     size = 10)) +
  coord_fixed() 

plot1

#alleliq frequency
plot2 <-ggplot(data2, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "#AF43F2", mid = "#EBD0FB", 
                       midpoint = 0.8, limit = c(0.6,1), space = "Lab", 
                       name="Allelic freq") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

plot2

### DAPC ANALYSIS

library(genepop)
library(adegenet)
library(ade4)

data<-read.genepop("Newdataset_RDA.gen", ncode = 3L, quiet = FALSE)
data<-read.genepop("Newdataset_neutral.gen", ncode = 3L, quiet = FALSE)
data 
allelic.richness(data,min.n=NULL,diploid=TRUE)

# To know the optimal number of PCs retained
pdf("Arbacia_RDA_alfa_scores.pdf")
dapc1 <- dapc(data, n.da=1, n.pca=60) 
pH <- a.score(dapc1)
names(pH)
pH <- optim.a.score(dapc1)
dev.off()

# To know the number of clusters
grp <- find.clusters(data, max.n.clust=5)

names (grp)
head(grp$Kstat, 2) # it gives me the BIC of each K
head(grp$grp, 2)
table(pop(data), grp$grp)

pdf("Arbacia_RDA_K2.pdf")
table.value(table(pop(data), grp$grp), col.lab=paste("cluster", 1:6),row.lab=paste("pop", 1:11))
dev.off()

dapc<-dapc(data) 
dapc

pdf("DAPC_24_2_spatial_1_RDA.pdf")
scatter(dapc, posi.da="bottomright", posi.pca="topright", bg="white", pch=c(7, 10, 17, 16, 15, 0, 17,16, 5, 17, 16, 17, 16), cstar=1, clab=1, leg=TRUE, posi.leg = "topleft", cleg=0.80, scree.da=TRUE, cex.lab=1, scree.pca=TRUE, col=c("royalblue1", "lightgoldenrod1", "lightseagreen", "lightcoral", "lightsalmon"))
dev.off()

### STRUCTURE ANALYSIS

library(pophelper)
library(plotly)

#upload files

sfiles <- list.files(path = ".", pattern = "*_f") #from results carpet after structure
slist <- readQ(sfiles, filetype = "structure", indlabfromfile = T)
attributes(slist)

#tabulateQ and summarize####
slistab <- tabulateQ(slist)
slistsum <- summariseQ(slistab)

evannoMethodStructure(slistsum)

#align and merge####
slist_1 <- alignK(slist)
slist_2 <- mergeQ(slist_1)

#select colors
myCol <- c("deepskyblue3", "lightcyan2", "aquamarine3", "cyan4", "palegreen3", "turquoise2", "mediumpurple1","orchid4")
  
#graph RDA
p3 <- plotQ(as.qlist(slist_2), indlabsize = 3.5, indlabangle = 45, panelspacer = 0.3, barbordercolour="white", barbordersize = 0.05, showlegend=T, showyaxis=T,showticks=T, showindlab=T,showtitle=T,titlelab="Arbacia RDA Structure", exportplot=F,returnplot=T,imgoutput="join", clustercol = myCol)
p3

#graph neutral
p3 <- plotQ(as.qlist(slist_2), indlabsize = 3.5, indlabangle = 45, panelspacer = 0.3, barbordercolour="white", barbordersize = 0.05, showlegend=T, showyaxis=T,showticks=T, showindlab=T,showtitle=T,titlelab="Arbacia Neutrals Structure", exportplot=F,returnplot=T,imgoutput="join", clustercol = myCol)
p3


### MALTEN TEST
library(AFLP)
library(poppr)
library(dismo)
library(ape)
library(raster)
library(readr)

#make the matrix with geographical distances
Matrix_popcoord <- read_delim("Matrix_popcoord.csv", delim = ";", escape_double = FALSE, col_types = cols(V1 = col_number(), 
                                                                                   V2 = col_number(), T1 = col_number(), 
                                                                                   T2 = col_number(), C = col_number()), 
                              trim_ws = TRUE)
coord.dist<-Matrix_popcoord
as.matrix(coord.dist)

##import fst square matrix
#RDA
Matrix_FST_RDA <- read_delim("Matrix_FST_RDA.csv", delim = ";", escape_double = FALSE, col_types = cols(V1 = col_number(), 
                                                                                  V2 = col_number(), T1 = col_number(), 
                                                                                  T2 = col_number(), C = col_number()), 
                             trim_ws = TRUE)
fst.dist<-Matrix_FST_RDA
as.matrix(fst.dist)

#neutral
Matrix_FST_neutral <- read_delim("Matrix_FST_neutral.csv", 
                                 delim = ";", escape_double = FALSE, col_types = cols(V1 = col_number(), 
                                                                                      V2 = col_number(), T1 = col_number(), 
                                                                                      T2 = col_number(), C = col_number()), 
                                 trim_ws = TRUE)
fst.dist<-Matrix_FST_neutral
as.matrix(fst.dist)

##With the as.dist function we manage to accept the distance matrix without having to make it from 0.
data<-mantel.rtest(as.dist(coord.dist), as.dist(fst.dist), nrepet = 99999)
plot(data)
data


### ANNOTATION AND REVIGO REPRESENTATION
# Select annotation information of blast sequencing A.lixula
library(dplyr)
library(readr)

#load annotation information
Annotation <- read_delim("Annotation_assembly_CT+T22.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
Annotation

#load blast seq table
Blast <- read_delim("105_Blastresult_95_temp.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
Blast

merged.data <- merge(Annotation, Blast, by="Sequence Name")
write.table(merged.data, file="Blast_result_95_annotation.csv", sep= ",", row.names = TRUE,col.names = TRUE)

## treemap R script produced by the Revigo server at http://revigo.irb.hr/ (Anton Kratz)
library(treemap)
library(RColorBrewer)

#Biological process graph
revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006810","transport",18.3875755434523,1,1,-0,"transport"),
c("GO:0007154","cell communication",8.685360854268644,2,0.9859119590163058,0.01071393,"cell communication"),
c("GO:0007156","homophilic cell adhesion via plasma membrane adhesion molecules",0.21572517775130376,2,0.9903216956570806,0.00719301,"homophilic cell adhesion via plasma membrane adhesion molecules"),
c("GO:0007275","multicellular organism development",1.2483752134872164,3,1,-0,"multicellular organism development"),
c("GO:0008283","cell population proliferation",0.11222216744817197,1,0.9908477983887931,0.00679793,"cell population proliferation"),
c("GO:0009069","serine family amino acid metabolic process",0.606086131244822,4,0.897316715785535,0,"serine family amino acid metabolic process"),
c("GO:0006468","protein phosphorylation",4.356287697349047,4,0.780615755891554,0.19737155,"serine family amino acid metabolic process"),
c("GO:0006470","protein dephosphorylation",0.5032147030452196,1,0.8165504383777462,0.52180445,"serine family amino acid metabolic process"),
c("GO:0018108","peptidyl-tyrosine phosphorylation",0.28779866811301186,1,0.8195269049457224,0.5252144,"serine family amino acid metabolic process"),
c("GO:0040007","growth",0.09945091249154055,1,1,-0,"growth"),
c("GO:0044238","primary metabolic process",54.341480421137625,2,0.9347622549417456,0.05078471,"primary metabolic process"),
c("GO:0043170","macromolecule metabolic process",38.67120377521221,2,0.9200753466945987,0.18607852,"primary metabolic process"),
c("GO:0060070","canonical Wnt signaling pathway",0.02512367125513281,1,0.9334768210355101,0.00603837,"canonical Wnt signaling pathway"),
c("GO:0007165","signal transduction",8.104950197590462,1,0.9278061937307486,0.42409333,"canonical Wnt signaling pathway"),
c("GO:0090305","nucleic acid phosphodiester bond hydrolysis",1.0173790452958384,1,0.9119286197106595,0.07308809,"nucleic acid phosphodiester bond hydrolysis"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="revigo_treemap.pdf", width=16, height=9) # width and height are in inches

#treemap graph
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  palette = "Set3"
  inflate.labels = F,      # set this to TRUE for space-filling group labels
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "transparent",   # define background color of group labels
  position.legend = "none"
)
dev.off()
