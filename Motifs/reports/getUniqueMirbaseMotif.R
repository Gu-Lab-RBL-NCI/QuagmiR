
## PURPOSE: generate unique motifs from miRBase data
## INPUT: miRBase data  miRNA.xls
## OUTPUT: fasta file motif_list_species_mirbase_22.fa
##         table containing motifs  species-miRBase22.tsv

require(gdata)
require(plyr)
require(data.table)


base_path<-'C:/Users/Susanna/Documents/NCI/miRBase/'
setwd(base_path)
miRBase_file <- 'miRNA.xls' #file from miRBase v.22
perl <- 'C:/Strawberry/perl/bin/perl5.26.2.exe'
miRBase_data <- read.xls(miRBase_file, sheet = 1, header = TRUE, perl=perl)

motif_len <- 13 #length of motif
include_pos <- F #include motif position in fasta file?
motif_start <- 1 #set start position of motif

#formats file from miRBase
getMirnaData <- function(miRBase_file){
  miRBase_data <- read.xls(miRBase_file, sheet = 1, header = TRUE, perl=perl)
  miRBase_data <- miRBase_data[, -grep("X", colnames(miRBase_data), fixed=T)]
  write.table(miRBase_data, 'miRNA.tsv', sep='\t', row.names=F)
  miRBase_data <- read.table('miRNA.tsv', sep='\t', header=T, stringsAsFactors=F)
  miRBase_data <- na.omit(miRBase_data)
  miRBase_data <- Filter(function(x)!all(is.na(x)), miRBase_data)
  
  names(miRBase_data)[names(miRBase_data) == 'Accession'] <- 'PRI.ACCESSION'
  names(miRBase_data)[names(miRBase_data) == 'Sequence'] <- 'PRI.SEQUENCE'
  names(miRBase_data)[names(miRBase_data) == 'ID'] <- 'PRIMIRNA'
  
  mature2_names <- c('PRI.ACCESSION', 'PRIMIRNA', 'PRI.SEQUENCE', 'Mature2_Acc', 'Mature2_ID', 'Mature2_Seq')
  mature2 <- miRBase_data[, mature2_names]
  miRBase_data <- miRBase_data[,!names(miRBase_data) %in% c('Mature2_Acc', 'Mature2_ID', 'Mature2_Seq')]
  names(mature2)[names(mature2) == 'Mature2_Acc'] <- 'ACCESSION'
  names(mature2)[names(mature2) == 'Mature2_ID'] <- 'MIRNA'
  names(mature2)[names(mature2) == 'Mature2_Seq'] <- 'SEQUENCE'
  names(miRBase_data)[names(miRBase_data) == 'Mature1_Acc'] <- 'ACCESSION'
  names(miRBase_data)[names(miRBase_data) == 'Mature1_ID'] <- 'MIRNA'
  names(miRBase_data)[names(miRBase_data) == 'Mature1_Seq'] <- 'SEQUENCE'
  miRBase_data <- rbind.fill(miRBase_data, mature2)
  
  miRBase_data$SEQUENCE<-gsub('U', 'T', miRBase_data$SEQUENCE)
  miRBase_data$PRI.SEQUENCE<-gsub('U', 'T', miRBase_data$PRI.SEQUENCE)
  
  #miRBase_data <- miRBase_data[, -grep("X", colnames(miRBase_data), fixed=T)]
  miRBase_data <- miRBase_data[!(is.na(as.vector(miRBase_data$SEQUENCE)) | nchar(miRBase_data$SEQUENCE)==0 | as.character(miRBase_data$SEQUENCE)==''), ]
  return(miRBase_data)
}

#get miRNA that are indistinguishable from given miRNA
getDupliMirna<-function(mirna, mirna_data){
  sequence<-mirna_data[which(mirna_data$MIRNA==mirna), 'SEQUENCE'][1]
  unique_mirna<-unique(mirna_data$MIRNA)
  dupli_id<-c()
  for(m in unique_mirna){
    if(!identical(toString(m), toString(mirna))){
      curr_seq<-mirna_data[which(mirna_data$MIRNA==m), 'SEQUENCE']
      for(c in curr_seq){
        if(gregexpr(sequence, c, fixed=TRUE)[[1]]!=-1 || gregexpr(c, sequence, fixed=TRUE)[[1]]!=-1){
          dupli_id<-c(dupli_id, m)
        }
      }
    }
  }
  return(unique(dupli_id))
}

#generate motif at center of sequence
getUniqueMotif<-function(length, seq, primary_sequences, start){
  #mid<-trunc(nchar(seq)/2)
  #start<-mid-(trunc(length/2))
  #end<-mid+(trunc(length/2))
  motif<-substr(seq, start, start + length - 1)
  if(!inPriSeq(motif, primary_sequences)){
    return(motif)
  }else{
    return(NA)
  }
  
}

#check if sequence may be found in list of sequences
inPriSeq<-function(seq, primary_sequences){
  for(i in 1:length(primary_sequences)){
    p<-primary_sequences[i]
    found<-gregexpr(seq, p, fixed=TRUE)[[1]]
    if(found!=-1){
      return(TRUE)
    }
  }
  return(FALSE)
}

#get position of motif within sequence
getMotifPosition<-function(motif, primary_seq){
  
  start<-gregexpr(motif, primary_seq, fixed=TRUE)[[1]]
  end<-start+nchar(motif)
  return(c(start, end))
  
}

#get miRNA sequences excluding sequences associated with given miRNA
getPrimarySequences<-function(mirna){
  accession<-mirna_data[which(mirna_data$MIRNA==mirna),'ACCESSION']
  pri_accession<-mirna_data[which(mirna_data$MIRNA==mirna),'PRI.ACCESSION']
  dupli_id<-as.vector(strsplit(mirna_data[which(mirna_data$MIRNA==mirna), 'DUPLI.ID'], ',', fixed=TRUE)[[1]])
  primary_seq<-c()
  #get one primary sequence per accession
  for(i in 1:nrow(mirna_data)){
    curr_acc<-mirna_data[i,'ACCESSION']
    curr_pri_acc<-mirna_data[i,'PRI.ACCESSION']
    curr_mirna<-mirna_data[i, 'MIRNA']
    curr_seq<-mirna_data[i,'SEQUENCE']
    if(!(curr_acc%in%accession) && !(curr_pri_acc%in%pri_accession) && !(curr_mirna%in%dupli_id)){
      primary_seq<-c(primary_seq, curr_seq)
    }
  }
  primary_seq<-unique(primary_seq)
  primary_seq<-primary_seq[!is.na(primary_seq)]
  return(primary_seq)
}


miRBase_data <- getMirnaData(miRBase_file)
species_list <- unique(sapply(strsplit(as.character(miRBase_data$PRIMIRNA), "-", fixed=T), `[`, 1))
new_species_list <- species_list

#generate motifs
for(species in species_list[1:2]){
  
  if(file.exists(paste0(species, '-miRBase22.tsv'))){ next }
  
  mirna_data <- miRBase_data[miRBase_data$PRIMIRNA %like% species, ]
  
  if(length(unique(mirna_data$PRI.ACCESSION)) == 1 || length(unique(mirna_data$SEQUENCE)) == 1){  #skip if species only contains one primiRNA sequence
    new_species_list <- new_species_list[new_species_list!=species]
    next 
  }
  
  
  #rename miRNA which have paralogs
  for(a in unique(mirna_data$ACCESSION)){
    
    mirna_id <- toString(unique(mirna_data[which(mirna_data$ACCESSION==a), 'MIRNA']))
    
    num_paralogs <- nrow(mirna_data[which(mirna_data$ACCESSION==a),])
    
    if(num_paralogs > 1){
      for(i in 1:num_paralogs){
        mirna_id <- paste0(mirna_id, '-', i)
      }
    }
    
    mirna_data[which(mirna_data$ACCESSION==a), 'MIRNA'] <- toString(mirna_id)
  }
  
  
  
  unique_mirna<-unique(mirna_data$MIRNA)
  unique_mirna<-unique_mirna[!is.na(unique_mirna)]
  for(m in unique_mirna){
    #find miRNA each miRNA is indistinguishable from 
    mirna_data[which(mirna_data$MIRNA==m),'DUPLI.ID']<-paste(getDupliMirna(m, mirna_data), collapse=',')
    
    temp_data<-mirna_data[which(mirna_data$MIRNA==m),]
    curr_seq<-temp_data[1, 'SEQUENCE']
    pri_seq<-getPrimarySequences(m)
    
    #get unique motif for miRNA and generate id containing motif positions
    unique_motif<-getUniqueMotif(motif_len, curr_seq, pri_seq, motif_start)
    mirna_data[which(mirna_data$MIRNA==m), paste0('MOTIF.', motif_len)]<-unique_motif
    if(!is.na(unique_motif)){
      motif_position<-getMotifPosition(unique_motif, curr_seq)
      new_mirna_id<-paste(m, '[', motif_position[1], ',', motif_position[2], ']', sep='')
      mirna_data[which(mirna_data$MIRNA==m),'POS.MIRNA']<-new_mirna_id
    }
  }
  
  mirna_data$DUPLI.MIRNA<-NA
  
  #generate new id for combining miRNA which are indistinguishable from each other
  for(i in 1:nrow(mirna_data)){
    dupli_id<-mirna_data[i, 'DUPLI.ID']
    dupli_id<-strsplit(dupli_id,',', fixed=TRUE)[[1]]
    mirna<-mirna_data[i, 'MIRNA']
    motif<-mirna_data[i, paste0('MOTIF.', motif_len)]
    if(is.na(motif)){ next }
    if(length(dupli_id>0)){
      new_mirna_id<-c(toString(mirna))
      new_mirna_id<-paste(new_mirna_id, collapse=',')
      
      for(d in dupli_id){
        temp_dupli<-strsplit(d, '-', fixed=TRUE)[[1]]
        temp_dupli<-paste(temp_dupli[3:length(temp_dupli)], collapse='-')
        new_mirna_id<-paste(new_mirna_id, temp_dupli, sep='-', collapse='-')
      }
      
      for(d in dupli_id){
        if(is.na(mirna_data[which(mirna_data$MIRNA==d),'DUPLI.MIRNA'])){
          mirna_data[which(mirna_data$MIRNA==d),'DUPLI.MIRNA']<-new_mirna_id
          mirna_data[which(mirna_data$MIRNA==d),'MOTIF.13']<-motif
          curr_seq<-mirna_data[which(mirna_data$MIRNA==d),'SEQUENCE']
          motif_start<-gregexpr(motif, curr_seq, fixed=TRUE)[[1]][1]
          motif_end<-motif_start+nchar(motif)
          mirna_data[which(mirna_data$MIRNA==d),'POS.MIRNA']<-paste(new_mirna_id, '[', motif_start, ',', motif_end, ']', sep='')
        }
      }
      
      if(is.na(mirna_data[which(mirna_data$MIRNA==mirna), 'DUPLI.MIRNA'])){
        mirna_data[which(mirna_data$MIRNA==mirna), 'DUPLI.MIRNA']<-new_mirna_id
        curr_seq<-mirna_data[i,'SEQUENCE']
        motif_start<-gregexpr(motif, curr_seq, fixed=TRUE)[[1]][1]
        motif_end<-motif_start+nchar(motif)
        mirna_data[which(mirna_data$MIRNA==mirna),'POS.MIRNA']<-paste(new_mirna_id, '[', motif_start, ',', motif_end, ']', sep='')
      }
    }
  }
  
  
  
  write.table(mirna_data, paste0(species, '-miRBase22.tsv'), sep='\t', row.names=F)
  
}

species_list <- new_species_list #update list excluding species with only one primiRNA


#generate fasta file with mirna id w/ position and new unique motifs
for(species in species_list){
  
  mirna_file <- paste0(species, '-miRBase22.tsv')
  if(!file.exists(mirna_file)){ next }
  mirna_data <- read.table(mirna_file, sep='\t', header=T)
  data<-c()
  unique_mirna<-unique(mirna_data$POS.MIRNA)
  unique_mirna <- na.omit(unique_mirna)
  if(is.null(unique_mirna)){ next }
  for(m in unique_mirna){
    temp_data<-mirna_data[which(mirna_data$POS.MIRNA==m),]
    id <- m
    if(!include_pos){
      id <- unique(mirna_data[which(mirna_data$POS.MIRNA==toString(m)),'DUPLI.MIRNA'])
      if(is.na(id)){
        id <- toString(unique(mirna_data[which(mirna_data$POS.MIRNA==m),'MIRNA']))
      }
      
    }
    temp_data<-temp_data[1,]
    sequence<-temp_data[,'SEQUENCE']
    motif<-temp_data[,paste0('MOTIF.', motif_len)]
    
    first_line<-paste(paste('>', toString(id), sep=''), toString(motif), sep='\t')
    second_line<-toString(sequence)
    
    data<-c(data, first_line, '\n', second_line, '\n')
  }
  data <- data[-length(data)]
  motifs_file<-paste('motif_list_', species, '_mirbase_22.fa', sep='')
  file<-file(motifs_file)
  writeLines(data, file, sep='')
  close(file)
  
}


#calculate motif coverage
motif_coverage <- data.frame(SPECIES=as.character(species_list), COVERAGE=NA, stringsAsFactors = F)
for(species in species_list){
  
  mirna_file <- paste0(species, '-miRBase22.tsv')
  if(!file.exists(mirna_file)){ next }
  mirna_data <- read.table(mirna_file, sep='\t', header=T, stringsAsFactors = F)
  
  num_seq <- nrow(mirna_data)
  covered <- nrow(mirna_data[which(!is.na(as.vector(mirna_data$POS.MIRNA))),])
  
  motif_coverage[which(motif_coverage$SPECIES==species), 'COVERAGE'] <- as.numeric(covered/num_seq)
}
write.table(motif_coverage, 'motif_coverage.tsv', sep='\t', row.names=F)

