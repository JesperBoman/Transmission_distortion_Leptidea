library(ggplot2)
library(reshape2)
library(ggsignif)
library(car)

#Reading and formatting data

setwd("~/Downloads")
dataF2_frq <- data.frame()
for(sample in c("BE", "DE", "DLDP", "F1ad", "F2ad_female", "F2ad_male", "F2RP")){
  tmp<- read.table(paste(sample, ".Swe.fixed.diff_noDel.mq20.af6.SweFreq_AncFreq.pol", sep=""))
  tmp$Sample <- paste(sample, "Swe", sep=".")
  dataF2_frq<-rbind(dataF2_frq, tmp)
}

for(sample in c("BE", "DE", "DLDP", "F1ad", "F2ad_female", "F2ad_male", "F2RP")){
  tmp<- read.table(paste(sample, ".Cat.fixed.diff_noDel.mq20.af6.SweFreq_AncFreq.pol", sep=""))
  tmp$Sample <- paste(sample, "Cat", sep=".")
  dataF2_frq<-rbind(dataF2_frq, tmp)
}

#colnames(dataF2_frq) <- c("Chromosome", "Position", "Coverage", "SWE_Frequency", "Swedish_allele", "Catalan_allele", "Sample")
colnames(dataF2_frq) <- c("Chromosome", "Position", "Coverage", "SWE_Frequency", "Swedish_allele", "Catalan_allele", "Ancestral_allele","ANC_Frequency", "Sample")

faidx <- read.table("Lsin_DToL.fasta.fai")
chr_file <- faidx[,1:2]
chr_file$Start <- 1 
colnames(chr_file) <- c("Chromosome", "End", "Start")
chr_file$Chr_num <- as.numeric(gsub("Chr_", "", chr_file$Chromosome))
chr_file$Chr_type <- ifelse(chr_file$Chr_num  == 2 | chr_file$Chr_num  == 3 | chr_file$Chr_num == 48, "Z", "A" )



dataF2_frq$region_ID <- paste(dataF2_frq$Chromosome, dataF2_frq$Position, sep="_")
#dataF2_frq_wide <- dcast(dataF2_frq, Chromosome+Position+region_ID~Sample, value.var = "SWE_Frequency")
dataF2_frq_wide <- dcast(dataF2_frq, Chromosome+Position+region_ID+Ancestral_allele~Sample, value.var = "SWE_Frequency")
dataF2_frq_wide_cov <- dcast(dataF2_frq, Chromosome+Position+region_ID~paste(Sample, ".Cov", sep=""), value.var = "Coverage")

dataF2_frq_wide <- cbind(dataF2_frq_wide, dataF2_frq_wide_cov[,4:17])
dataF2_frq_wide$Chr_num <- as.numeric(gsub("Chr_", "", dataF2_frq_wide$Chromosome))



## F2 BACKCROSS - MEIOTIC DRIVE ####

# Backcross - filtering #
cutoff=1
BE.filt <- dataF2_frq_wide[abs(dataF2_frq_wide$BE.Swe - dataF2_frq_wide$BE.Cat) < cutoff, c("Chromosome", "Position", "BE.Cat", "BE.Swe", "BE.Cat.Cov", "BE.Swe.Cov")]
#Here sites could be filtered

BE_data.Swe <- aggregate(SWE_Frequency~Chromosome,dataF2_frq[dataF2_frq$Sample == "BE.Swe" ,], mean)
BE_data.Cat <- aggregate(SWE_Frequency~Chromosome,dataF2_frq[dataF2_frq$Sample == "BE.Cat" ,], mean)
BE_data <- as.data.frame(cbind(Chromosome=BE_data.Swe[,1], Frequency=(BE_data.Swe[,2]+BE_data.Cat[,2])/2))
BE_data$Frequency <- as.numeric(BE_data$Frequency)




BE.filt$Frequency <- (BE.filt$BE.Cat+BE.filt$BE.Swe)/2
BE_data <- aggregate(Frequency~Chromosome,BE.filt, mean)

BE_data$Swe_counts <- BE_data$Frequency*67

fus_count <- (BE_data[BE_data$Chromosome == "Chr_7",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_24",]$Swe_counts)/2
fus_count <- fus_count+(BE_data[BE_data$Chromosome == "Chr_16",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_27",]$Swe_counts)/2
fus_count <- fus_count+(BE_data[BE_data$Chromosome == "Chr_20",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_31",]$Swe_counts)/2
fus_count <- fus_count+(BE_data[BE_data$Chromosome == "Chr_14",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_35",]$Swe_counts)/2
fus_count <- fus_count+(BE_data[BE_data$Chromosome == "Chr_25",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_37",]$Swe_counts)/2
fus_count <- fus_count+(BE_data[BE_data$Chromosome == "Chr_45",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_47",]$Swe_counts)/2

fus_freq_vector <- (BE_data[BE_data$Chromosome == "Chr_7",]$Frequency+BE_data[BE_data$Chromosome == "Chr_24",]$Frequency)/2
fus_freq_vector <- c(fus_freq_vector, (BE_data[BE_data$Chromosome == "Chr_16",]$Frequency+BE_data[BE_data$Chromosome == "Chr_27",]$Frequency)/2)
fus_freq_vector <- c(fus_freq_vector, (BE_data[BE_data$Chromosome == "Chr_20",]$Frequency+BE_data[BE_data$Chromosome == "Chr_31",]$Frequency)/2)
fus_freq_vector <- c(fus_freq_vector, (BE_data[BE_data$Chromosome == "Chr_14",]$Frequency+BE_data[BE_data$Chromosome == "Chr_35",]$Frequency)/2)
fus_freq_vector <- c(fus_freq_vector, (BE_data[BE_data$Chromosome == "Chr_25",]$Frequency+BE_data[BE_data$Chromosome == "Chr_37",]$Frequency)/2)
fus_freq_vector <- c(fus_freq_vector, (BE_data[BE_data$Chromosome == "Chr_45",]$Frequency+BE_data[BE_data$Chromosome == "Chr_47",]$Frequency)/2)





fis_count <- (BE_data[BE_data$Chromosome == "Chr_17",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_34",]$Swe_counts)/2
fis_count <- fis_count+(BE_data[BE_data$Chromosome == "Chr_40",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_15",]$Swe_counts)/2
fis_count <- fis_count+(BE_data[BE_data$Chromosome == "Chr_39",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_38",]$Swe_counts)/2
fis_count <- fis_count+(BE_data[BE_data$Chromosome == "Chr_10",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_18",]$Swe_counts)/2
fis_count <- fis_count+(BE_data[BE_data$Chromosome == "Chr_5",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_28",]$Swe_counts)/2

fis_freq_vector <- (BE_data[BE_data$Chromosome == "Chr_17",]$Frequency+BE_data[BE_data$Chromosome == "Chr_34",]$Frequency)/2
fis_freq_vector <- c(fis_freq_vector, (BE_data[BE_data$Chromosome == "Chr_40",]$Frequency+BE_data[BE_data$Chromosome == "Chr_15",]$Frequency)/2)
fis_freq_vector <- c(fis_freq_vector, (BE_data[BE_data$Chromosome == "Chr_39",]$Frequency+BE_data[BE_data$Chromosome == "Chr_38",]$Frequency)/2)
fis_freq_vector <- c(fis_freq_vector, (BE_data[BE_data$Chromosome == "Chr_10",]$Frequency+BE_data[BE_data$Chromosome == "Chr_18",]$Frequency)/2)
fis_freq_vector <- c(fis_freq_vector, (BE_data[BE_data$Chromosome == "Chr_5",]$Frequency+BE_data[BE_data$Chromosome == "Chr_28",]$Frequency)/2)





unpol_count <- (BE_data[BE_data$Chromosome == "Chr_33",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_8",]$Swe_counts)/2
unpol_count <- unpol_count+(BE_data[BE_data$Chromosome == "Chr_41",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_46",]$Swe_counts)/2
unpol_count <- unpol_count+BE_data[BE_data$Chromosome == "Chr_4",]$Swe_counts
unpol_count <- unpol_count+(BE_data[BE_data$Chromosome == "Chr_44",]$Swe_counts+BE_data[BE_data$Chromosome == "Chr_36",]$Swe_counts)/2

unpol_freq_vector <- (BE_data[BE_data$Chromosome == "Chr_33",]$Frequency+BE_data[BE_data$Chromosome == "Chr_8",]$Frequency)/2
unpol_freq_vector <- c(unpol_freq_vector, (BE_data[BE_data$Chromosome == "Chr_41",]$Frequency+BE_data[BE_data$Chromosome == "Chr_46",]$Frequency)/2)
unpol_freq_vector <- c(unpol_freq_vector, (BE_data[BE_data$Chromosome == "Chr_4",]$Frequency))
unpol_freq_vector <- c(unpol_freq_vector, (BE_data[BE_data$Chromosome == "Chr_44",]$Frequency+BE_data[BE_data$Chromosome == "Chr_36",]$Frequency)/2)



hom_count <- BE_data[BE_data$Chromosome == "Chr_23",]$Swe_counts
hom_count <- hom_count+BE_data[BE_data$Chromosome == "Chr_19",]$Swe_counts
hom_count <- hom_count+BE_data[BE_data$Chromosome == "Chr_11",]$Swe_counts
hom_count <- hom_count+BE_data[BE_data$Chromosome == "Chr_12",]$Swe_counts
hom_count <- hom_count+BE_data[BE_data$Chromosome == "Chr_1",]$Swe_counts

hom_freq_vector <- BE_data[BE_data$Chromosome == "Chr_23",]$Frequency
hom_freq_vector <- c(hom_freq_vector, BE_data[BE_data$Chromosome == "Chr_19",]$Frequency)
hom_freq_vector <- c(hom_freq_vector, BE_data[BE_data$Chromosome == "Chr_11",]$Frequency)
hom_freq_vector <- c(hom_freq_vector, BE_data[BE_data$Chromosome == "Chr_12",]$Frequency)
hom_freq_vector <- c(hom_freq_vector, BE_data[BE_data$Chromosome == "Chr_1",]$Frequency)

binom.test(round(fus_count), 67*6, p=3/4 )
binom.test(round(fis_count), 67*5, p=3/4 )
binom.test(round(unpol_count), 67*4, p=3/4 )
binom.test(round(hom_count), 67*5, p=3/4 )



freq_df <- as.data.frame(rbind(cbind(fus_freq_vector, rep("Fusion SWE")), cbind(fis_freq_vector, rep("Fission CAT")), cbind(unpol_freq_vector, rep("Unknown \n polarization")), cbind(hom_freq_vector, rep("Homologus"))))
colnames(freq_df) <- c("Frequency", "Type")
freq_df$Frequency <- as.numeric(freq_df$Frequency)

ggplot(freq_df, aes(y=Frequency, x=Type, fill=Type))+geom_point(size=3, shape=21, position=position_dodge2(width=0.2))+
  geom_hline(yintercept=0.75, lty=2)+
  theme_classic()+
  xlab("Chromosome type")+
  scale_fill_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40"))+
  geom_signif(y_position = 0.8, xmin = 1.75, tip_length = 0, textsize = 15,
              xmax = 2.25, annotation = c("*"), col="black")+
  theme(aspect.ratio=1, title= element_text(size=14), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=20), legend.text=element_text(size=18),  legend.title=element_text(size=18))






## F2 INTERCROSS - MEIOTIC DRIVE ####

# Inter - filtering ####

df <- dataF2_frq_wide
filt_df <- data.frame()
cutoff<-1
for(i in 1:length(df$Chromosome)){
  
  DE_sum <- ifelse(abs(df[i,7]-df[i,8])<cutoff, (df[i,7]+df[i,8])/2, NA)
  DLDP_sum <- ifelse(abs(df[i,9]-df[i,10])<cutoff, (df[i,9]+df[i,10])/2, NA)
  fem_sum <- ifelse(abs(df[i,13]-df[i,14])<cutoff, (df[i,13]+df[i,14])/2, NA)
  male_sum <- ifelse(abs(df[i,15]-df[i,16])<cutoff, (df[i,15]+df[i,16])/2, NA)
  rp_sum <- ifelse(abs(df[i,17]-df[i,18])<cutoff, (df[i,17]+df[i,18])/2, NA)
  
  n<-sum(!is.na(c(DE_sum, DLDP_sum, fem_sum, male_sum, rp_sum)))
  AVG <- mean(c(DE_sum, DLDP_sum, fem_sum, male_sum, rp_sum), na.rm=T)
  WAVG<- (DE_sum*298+DLDP_sum*72+fem_sum*80+male_sum*76+rp_sum*72)/598
  filt_df <- rbind(filt_df, cbind(df[i,1:2], AVG, WAVG, n))
}

filt_df<-filt_df[filt_df$n == 5,]

head(filt_df)


fus_inter_count <- ((mean(filt_df[filt_df$Chromosome == "Chr_7",]$WAVG)
                     +mean(filt_df[filt_df$Chromosome == "Chr_24",]$WAVG))/2)*598

fus_inter_count <- fus_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_16",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_27",]$WAVG))/2)*598

fus_inter_count <- fus_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_20",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_31",]$WAVG))/2)*598

fus_inter_count <- fus_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_14",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_35",]$WAVG))/2)*598

fus_inter_count <- fus_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_25",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_37",]$WAVG))/2)*598


fus_inter_count <- fus_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_45",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_47",]$WAVG))/2)*598

binom.test(round(fus_inter_count), 598*6, p=1/2 )



fus_inter_vector <- ((mean(filt_df[filt_df$Chromosome == "Chr_7",]$WAVG)
                      +mean(filt_df[filt_df$Chromosome == "Chr_24",]$WAVG))/2)
fus_inter_vector <- c(fus_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_16",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_27",]$WAVG))/2))
fus_inter_vector <- c(fus_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_20",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_31",]$WAVG))/2))
fus_inter_vector <- c(fus_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_14",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_35",]$WAVG))/2))
fus_inter_vector <- c(fus_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_25",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_37",]$WAVG))/2))
fus_inter_vector <- c(fus_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_45",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_47",]$WAVG))/2))






#Simple fissions inter

fis_inter_count <- ((mean(filt_df[filt_df$Chromosome == "Chr_17",]$WAVG)
                     +mean(filt_df[filt_df$Chromosome == "Chr_34",]$WAVG))/2)*598

fis_inter_count <- fis_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_40",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_15",]$WAVG))/2)*598

fis_inter_count <- fis_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_39",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_38",]$WAVG))/2)*598

fis_inter_count <- fis_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_10",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_18",]$WAVG))/2)*598

fis_inter_count <- fis_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_5",]$WAVG)
                                     +mean(filt_df[filt_df$Chromosome == "Chr_28",]$WAVG))/2)*598


binom.test(round(fis_inter_count), 598*5, p=1/2 )


fis_inter_vector <- ((mean(filt_df[filt_df$Chromosome == "Chr_17",]$WAVG)
                      +mean(filt_df[filt_df$Chromosome == "Chr_34",]$WAVG))/2)
fis_inter_vector <- c(fis_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_40",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_15",]$WAVG))/2))
fis_inter_vector <- c(fis_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_39",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_38",]$WAVG))/2))
fis_inter_vector <- c(fis_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_10",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_18",]$WAVG))/2))
fis_inter_vector <- c(fis_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_5",]$WAVG)
                                          +mean(filt_df[filt_df$Chromosome == "Chr_28",]$WAVG))/2))




#Unpol inter

unpol_inter_count <- ((mean(filt_df[filt_df$Chromosome == "Chr_33",]$WAVG)
                       +mean(filt_df[filt_df$Chromosome == "Chr_8",]$WAVG))/2)*598

unpol_inter_count <- unpol_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_41",]$WAVG)
                                         +mean(filt_df[filt_df$Chromosome == "Chr_46",]$WAVG))/2)*598

unpol_inter_count <- unpol_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_44",]$WAVG)
                                         +mean(filt_df[filt_df$Chromosome == "Chr_36",]$WAVG))/2)*598

unpol_inter_count <- unpol_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_4",]$WAVG)))*598

binom.test(round(unpol_inter_count), 598*4, p=1/2 )


unpol_inter_vector <- ((mean(filt_df[filt_df$Chromosome == "Chr_33",]$WAVG)
                        +mean(filt_df[filt_df$Chromosome == "Chr_8",]$WAVG))/2)
unpol_inter_vector <- c(unpol_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_41",]$WAVG)
                                              +mean(filt_df[filt_df$Chromosome == "Chr_46",]$WAVG))/2))
unpol_inter_vector <- c(unpol_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_44",]$WAVG)
                                              +mean(filt_df[filt_df$Chromosome == "Chr_36",]$WAVG))/2))
unpol_inter_vector <- c(unpol_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_4",]$WAVG))))





#Hom inter

hom_inter_count <- ((mean(filt_df[filt_df$Chromosome == "Chr_1",]$WAVG)))*598

hom_inter_count <- hom_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_11",]$WAVG)))*598

hom_inter_count <- hom_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_12",]$WAVG)))*598

hom_inter_count <- hom_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_19",]$WAVG)))*598

hom_inter_count <- hom_inter_count+((mean(filt_df[filt_df$Chromosome == "Chr_23",]$WAVG)))*598


binom.test(round(hom_inter_count), 598*5, p=1/2 )

hom_inter_vector <- ((mean(filt_df[filt_df$Chromosome == "Chr_1",]$WAVG)))
hom_inter_vector <- c(hom_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_11",]$WAVG))))
hom_inter_vector <- c(hom_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_12",]$WAVG))))
hom_inter_vector <- c(hom_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_19",]$WAVG))))
hom_inter_vector <- c(hom_inter_vector, ((mean(filt_df[filt_df$Chromosome == "Chr_23",]$WAVG))))



inter_freq_df <- as.data.frame(rbind(cbind(fus_inter_vector, rep("Fusion SWE")), cbind(fis_inter_vector, rep("Fission CAT")), cbind(unpol_inter_vector, rep("Unknown \n polarization")), cbind(hom_inter_vector, rep("Homologus"))))
colnames(inter_freq_df) <- c("Frequency", "Type")
inter_freq_df$Frequency <- as.numeric(inter_freq_df$Frequency)


ggplot(inter_freq_df, aes(y=Frequency, x=Type, fill=Type))+geom_point(size=3, shape=21, position=position_dodge2(width=0.2))+
  geom_hline(yintercept=0.5, lty=2)+
  theme_classic()+
  xlab("Chromosome type")+
  scale_fill_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40"))+
  geom_signif(y_position = c(0.6, 0.6), xmin = c(1.75, 3.75), xmax = c(2.25, 4.25), tip_length = 0, textsize = 15,
              annotation = c("*", "*"), col="black")+
  ylim(0.45,0.63)+
  theme(aspect.ratio=1, title= element_text(size=14), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=20), legend.text=element_text(size=18),  legend.title=element_text(size=18))




# Miscellaneous ####

# Fixed diffs. per chromosome ####

fd_per_chrom<-as.data.frame(table(dataF2_frq_wide$Chromosome))
mean(fd_per_chrom$Freq)


#Aneuploidy check ####
aggregate((BE.Swe.Cov+BE.Cat.Cov)/2~1, dataF2_frq_wide, mean)

aggregate((BE.Swe.Cov+BE.Cat.Cov)/2~Chromosome, dataF2_frq_wide, mean)

hist(aggregate((BE.Swe.Cov+BE.Cat.Cov)/2~Chromosome, dataF2_frq_wide, mean)[,2], breaks=30)
hist(aggregate(BE.Cat.Cov~Chromosome, dataF2_frq_wide, mean)[,2])


anDF <- aggregate((BE.Swe.Cov+BE.Cat.Cov)/2~Chromosome, dataF2_frq_wide, mean)
anDF <- aggregate(((DE.Swe.Cov+DE.Cat.Cov)*298+(DLDP.Swe.Cov+DLDP.Cat.Cov)*72+(F2ad_female.Swe.Cov+F2ad_female.Cat.Cov)*80+(F2ad_male.Swe.Cov+F2ad_male.Cat.Cov)*76+(F2RP.Swe.Cov+F2RP.Cat.Cov)*72)/(2*598)~Chromosome, dataF2_frq_wide, mean)

colnames(anDF) <- c("Chromosome", "Coverage")
anDF$Chr_type <- ifelse(anDF$Chromosome == "Chr_2" | anDF$Chromosome == "Chr_3" | anDF$Chromosome == "Chr_48", "Z", "A")
anDF$Chr_type <- ifelse(anDF$Chromosome == "Chr_7" | anDF$Chromosome == "Chr_24" | anDF$Chromosome == "Chr_16"
                        | anDF$Chromosome == "Chr_27" | anDF$Chromosome == "Chr_20" | anDF$Chromosome == "Chr_31"
                        | anDF$Chromosome == "Chr_14" | anDF$Chromosome == "Chr_35" | anDF$Chromosome == "Chr_25"
                        | anDF$Chromosome == "Chr_37" | anDF$Chromosome == "Chr_45" | anDF$Chromosome == "Chr_47" , "Fusion_SWE", anDF$Chr_type )

anDF$Chr_type <- ifelse(anDF$Chromosome == "Chr_17" | anDF$Chromosome == "Chr_34" | anDF$Chromosome == "Chr_40"
                        | anDF$Chromosome == "Chr_15" | anDF$Chromosome == "Chr_39" | anDF$Chromosome == "Chr_38"
                        | anDF$Chromosome == "Chr_10" | anDF$Chromosome == "Chr_18" | anDF$Chromosome == "Chr_5"
                        | anDF$Chromosome == "Chr_28" , "Fission_CAT", anDF$Chr_type )

anDF$Chr_type <- ifelse(anDF$Chromosome == "Chr_33" | anDF$Chromosome == "Chr_8" | anDF$Chromosome == "Chr_41"
                        | anDF$Chromosome == "Chr_46" | anDF$Chromosome == "Chr_44" | anDF$Chromosome == "Chr_36"
                        | anDF$Chromosome == "Chr_4", "Unpolarized", anDF$Chr_type )

anDF$Chr_type <- ifelse(anDF$Chromosome == "Chr_23" | anDF$Chromosome == "Chr_19" | anDF$Chromosome == "Chr_11"
                        | anDF$Chromosome == "Chr_12" | anDF$Chromosome == "Chr_1" , "Homologous", anDF$Chr_type )



Anova(lm(Coverage~Chr_type, type="II", anDF))

t<-TukeyHSD(aov(lm(Coverage~Chr_type, type="II", anDF[anDF$Chr_type != "A",])))
round(t$Chr_type, digits = 2)


ggplot(anDF[anDF$Chr_type != "A",], aes(x=Coverage, fill=Chr_type))+geom_histogram(col="black")+
  theme_classic()+
  xlab("Coverage")+
  ylab("Count")+
  scale_fill_manual(name="Chromosome type", values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40", "purple"))+
  geom_histogram()+
  theme(aspect.ratio=1, title= element_text(size=14), legend.position = "right", panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=20), legend.text=element_text(size=18),  legend.title=element_text(size=18))

