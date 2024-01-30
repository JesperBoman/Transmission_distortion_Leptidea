library(ggplot2)

#Use dataset from: 
#de Vos JM, Augustijnen H, Bätscher L, Lucek K (2020) Speciation through chromosomal fusion and fission in Lepidoptera. Philos Trans R Soc B Biol Sci 375:20190539. https://doi.org/10.1098/rstb.2019.0539

chrom_data <- read.table(file=file.choose(), sep = "\t", header=T)

chrom_data$Within_species_low <- ifelse(grepl("-", chrom_data$Chromosomes), sub("([0-9]+)-[0-9]+", "\\1", chrom_data$Chromosomes), chrom_data$Chromosomes )
chrom_data$Within_species_high <- ifelse(grepl("-", chrom_data$Chromosomes), sub("[0-9]+-([0-9]+)", "\\1", chrom_data$Chromosomes), chrom_data$Chromosomes )

max(as.double(chrom_data$Within_species_high), na.rm=T)
max(as.double(chrom_data$Within_species_low), na.rm=T)

chrom_data$Genus <- gsub(" " ,"", chrom_data$Genus)
chrom_data$Family <- gsub(" " ,"", chrom_data$Family)

chrom_data_min <- aggregate(as.double(Within_species_low)~Genus, chrom_data, min)
chrom_data_max <- aggregate(as.double(Within_species_high)~Genus, chrom_data, max)

genus_data<-cbind(chrom_data_min, chrom_data_max[,2] )
colnames(genus_data) <- c("Genus", "Min", "Max")

ggplot(genus_data, aes(y=reorder(Genus, Max),  xmin = Min, xmax=Max))+geom_linerange()+
  theme_classic()+
  geom_vline(xintercept=31, lty=2, alpha=0.5)+
  theme(axis.text.y=element_blank())

chrom_data$Within_species_high <- as.numeric(chrom_data$Within_species_high)
chrom_data$Within_species_low <- as.numeric(chrom_data$Within_species_low)

chrom_data$Group <- ifelse(chrom_data$Genus == "Leptidea", "Leptidea", "Other")
chrom_data$Butterfly <- ifelse(chrom_data$Family == "Lycaenidae" | chrom_data$Family == "Hesperidae" | chrom_data$Family == "Pieridae" | chrom_data$Family == "Riodinidae" | chrom_data$Family == "Nymphalidae" | chrom_data$Family == "Papilionidae", "Y", "N" )
chrom_data$Within_species_mean <- c(chrom_data$Within_species_high+chrom_data$Within_species_low)/2
chrom_data[chrom_data$Within_species_mean == min(chrom_data$Within_species_mean),]



ggplot(chrom_data, aes(x=Within_species_mean, y=reorder(Genus, Within_species_mean, FUN = max), col=Group, size=Group))+
  geom_point()+
  theme_classic()+
  scale_colour_manual(values = c("Red", "Black"))+
  scale_size_manual(values=c(2,0.2))+
  geom_vline(xintercept=31, lty=2, alpha=0.5)+
  ylab("Genus")+
  xlab("Haploid chromosome number")+
  theme( axis.text.y=element_blank(), axis.ticks.y = element_blank(), legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))


#log2 scale
ggplot(chrom_data, aes(x=log2(Within_species_mean), y=reorder(Genus, Within_species_mean, FUN = max), col=Group, size=Group))+
  geom_point()+
  theme_classic()+
  scale_colour_manual(values = c("Red", "Black"))+
  scale_size_manual(values=c(2,0.2))+
  geom_vline(xintercept=log2(31), lty=2, alpha=0.5)+
  ylab("Genus")+
  xlab("Haploid chromosome number")+
  scale_x_continuous(breaks = log2(seq(4.99, 225, 10)), labels = seq(5, 225, 10),
                     limits = log2(c(4.99, 225, 10)), guide = guide_axis(check.overlap = TRUE)) +
  
  theme(aspect.ratio = 1, axis.text.y=element_blank(), axis.ticks.y = element_blank(), legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

install.packages("ggridges")
library(ggridges)

#chrom_data
ggplot(chrom_data[chrom_data$Butterfly == "N",], aes(x=Within_species_mean, y=reorder(Family, Within_species_mean, FUN = sd), fill=after_stat(x)))+
  geom_density_ridges_gradient(scale=3)+
  theme_classic()+
  #scale_colour_manual(values = c("Red", "Black"))+
  scale_fill_viridis_c(option="magma", limits=c(0,40), oob = scales::squish)+
  geom_vline(xintercept=31, lty=2, alpha=0.5)+
  ylab("Family")+
  xlab("Haploid chromosome number")+
  theme( axis.ticks.y = element_blank(), legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))




ggplot(chrom_data, aes(x=Within_species_mean, fill=Butterfly))+geom_density(alpha=0.3)+
  theme_classic()+
  xlab("Chromosome number")
  theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

  aggregate(Within_species_high~Family+Butterfly, chrom_data, sd)

