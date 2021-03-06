library(Biostrings)
# library(BSgenome)

# folder where fasta file is (there should be only one)
ID_Folder <- "Pythium_Brazil"

# finds fasta files
ID_fasta_files <- list.files(path = ID_Folder, pattern = "\\.fas$|\\.fasta$", recursive = FALSE)
#ID_fasta_files <- list.files(path = ID_Folder, pattern = "fasta\\.fasta$", recursive = FALSE)

# # makes a linsi command from this file with reorientation
# # cmd <- paste("/opt/bio/mafft/bin/linsi --adjustdirection --auto --reorder ",  ID_Folder, "/", ID_fasta_files, " > ", ID_Folder, "/", "query_aligned_fasta.fasta", sep = "")
# 
# # makes a linsi command from this file without  reorientation
# cmd <- paste("/opt/bio/mafft/bin/linsi --reorder ",  ID_Folder, "/", ID_fasta_files, " > ", ID_Folder, "/", "query_aligned_fasta.fasta", sep = "")
# 
# 
# # runs the command on the linux server
# system(cmd)
# 
# #####################
# #  Need a script to trim the alignment automatically.  Currently doing alignment trimming by hand
# 
# temp <- readDNAStringSet(paste(ID_Folder, "/", "query_aligned_fasta.fasta", sep =""))
# mat <- t(consensusMatrix(temp, baseOnly=TRUE))
# 
# #######################
# 
# 
# library("ape")
# align <- read.dna(paste(ID_Folder, "/", "query_aligned_fasta.fasta", sep =""), format="fasta")
# 
# #dm <- dist.dna(align, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
# 
# dm <- dist.dna(align, model = "RAW", pairwise.deletion = FALSE, as.matrix = FALSE)
# 
# tree <- njs(dm)
# 
# 
# pdf(file = paste(ID_Folder, "/", "NJ_tree_of_ID.pdf", sep =""), width = 8, height =11 )
# plot.phylo(tree, type = "phylogram")
# dev.off()
# 
# tree[[1]]
# 
# fit <- hclust(dm, method="average") 
# 
# # Number of groups based on NJ tree
# 
# num_clades <- 10
# 
# groups <-  cutree(fit, num_clades)
# group2 <- data.frame(names(groups), groups)
# 
# #####################
# # Get length  still needs work
# 
# seq_no_gaps <- as.DNAbin(del.gaps(align))
# 
# seq_length <- data.frame(names(seq_no_gaps), sapply(seq_no_gaps, length))
# 
# seq_gr_len <- merge(seq_length, group2, by.x = "names.seq_no_gaps.", by.y = "names.groups.")
# 
# 
# 
# unique_ID <- names(groups[!duplicated(groups)])
# 
# # unique_ID <- sub("NFIS4250_ITS_TB.seq", "NFIS4243_ITS_TB.seq",unique_ID, ignore.case = FALSE)
# # unique_ID <- sub("NFIS4248_ITS_TB.seq", "NFIS4249_ITS_TB.seq",unique_ID, ignore.case = FALSE)
# 
# #fasta_for_GenBank <- seq_no_gaps[[names=unique_ID]]
# 
# fasta_for_GenBank <- align[unique_ID,]
# 
# i <- 1
# 
# for(i in 1:nrow(fasta_for_GenBank)) {
#   write.dna(fasta_for_GenBank[i,], file=paste(ID_Folder, "/GenBank/",rownames(fasta_for_GenBank)[i], ".fasta", sep=""), format = "fasta")
# }


j <- 1

for(j in 1:length(ID_fasta_files)) {
cmd2 <- paste("/opt/bio/ncbi/bin/blastall -p blastn -d /isilon/biodiversity/reference/ncbi/blastdb/reference/nt/nt -i ", ID_Folder, "/",
              ID_fasta_files[j], " -e 10 -m 8 -v 250 -b 250 -F F -o ", ID_Folder, "/GenBank/", ID_fasta_files[j], ".out", sep="")
system(cmd2)
}
              
##################################################################################################################
####################################################################################################################



k <- 3


for(k in 1:nrow(ID_fasta_files)) {
  GB_Blast_table <- read.table(paste(ID_Folder, "/GenBank/", ID_fasta_files[k], ".out", sep=""),header=FALSE)


summary(GB_Blast_table)

colnames(GB_Blast_table) <- c("query id", "subject_ids", " %identity", "alignment length", "mismatches", "gap opens", "q.start", " q.end", "s.start", "s.end", "evalue", " bit_score")
GB_Blast_table$subject_ids <- sub("\\|$", "| ",GB_Blast_table$subject_ids, ignore.case = FALSE)
parsed_col <- data.frame(matrix(unlist(strsplit(as.character(GB_Blast_table$subject_ids), "\\|")), nrow=length(GB_Blast_table$subject_ids), byrow=T),stringsAsFactors = FALSE)
GB_Blast_table <- data.frame(parsed_col[,4],GB_Blast_table, stringsAsFactors = FALSE)
colnames(GB_Blast_table)[1] <- "GB_accession"


unique(GB_Blast_table$GB_accession)




sequences <- read.GenBank(unique(GB_Blast_table$GB_accession),  seq.names = unique(GB_Blast_table$GB_accession),
                          species.names = TRUE,gene.names = FALSE,  as.character = TRUE)

# bad <- read.GenBank("AB217667.1",  seq.names = "AB217667.1",
#                           species.names = TRUE,gene.names = FALSE,  as.character = TRUE)
GB_DNAstring[8]


sequences[[1]]
names(sequences)[1]

sequences_ord <- sequences[order(names(sequences))] 

names(sequences_ord)[1]
sequences_ord[[1]]

full_names <- data.frame(attr(sequences, "species"),names(sequences), stringsAsFactors = FALSE)

full_names$SpStrain <-  do.call(paste, c(full_names[1:2], sep="_"))

full_names_ord <- full_names[order(full_names$names.sequences.),] 

nrow(full_names_ord)

# sequences.bin <- as.DNAbin(sequences)

library(Biostrings)
GB_DNAstring <- DNAStringSet()
for(i in 1:length(sequences_ord)) {
  temp <- DNAStringSet(paste(sequences_ord[[i]], collapse = ""))
  names(temp) <- full_names_ord$SpStrain[i]
  GB_DNAstring <- c(GB_DNAstring,temp)
} 


# Creates a +/- column for orientation of sequences
GB_Blast_table$orientation <- sign(GB_Blast_table$s.end - GB_Blast_table$s.start)

GB_Blast_table$length <- abs(GB_Blast_table$s.end - GB_Blast_table$s.start)

#min_length <- 50

#GB_Blast_table <-  subset(GB_Blast_table, GB_Blast_table$length > min_length) 


GB_Blast_table$query.id <- as.character(GB_Blast_table$query.id)

# Put together data from multiple hits on single line
dat_agg <- aggregate(GB_Blast_table[c(2,4,5,12,13,14,10,11)], by=list(GB_Blast_table$GB_accession, GB_Blast_table$subject_ids),FUN=c)

# find min and max and sense of hits
for(i in 1:nrow(dat_agg)) {
  dat_agg$min[i] <- min(unlist(c(dat_agg$s.start[i],dat_agg$s.end[i])))
  dat_agg$max[i] <- max(unlist(c(dat_agg$s.start[i],dat_agg$s.end[i])))
  dat_agg$sense[i] <- mean(unlist(c(dat_agg$orientation[i]))) 
}


# add a column with bp (width) of sequences.
dat_agg$bp <- width(GB_DNAstring)
dat_agg$width <- dat_agg$max - dat_agg$min
dat_agg$diff <- dat_agg$width - dat_agg$max


dat_agg <- dat_agg[order(as.character(dat_agg$Group.1)), ] 

# to make sure that name order are the same
check_order <- data.frame(names(GB_DNAstring), dat_agg[,1])

# Very Important
# THE ORDER OF THESE TWO COLUMNS SHOULD BE THE SAME
# SEQUENCES WILL BE OFF IF NOT
write.table(check_order, file = paste(ID_Folder,"/CHECK ORDER OF FASTA file AND BLAST TABLE.csv",sep=""), append = FALSE, sep = ",", col.names = NA)


# from here http://stackoverflow.com/questions/13545547/how-to-write-a-data-frame-with-one-column-a-list-to-a-file
dataset2 <- dat_agg # make a copy just to be on the safe side
dataset2[sapply(dataset2, is.list)] <-
  sapply(dataset2[sapply(dataset2, is.list)], 
         function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

write.xlsx(dataset2, file = paste(ID_Folder, "/to check BLAST results table.xlsx", sep=""), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)
# 


attributes(GB_DNAstring)

i <- 8



x_trim <- GB_DNAstring
for(i in 1:length(x_trim)) {
  if (dat_agg$sense[i] == 1) {
    x_trim[i] <- DNAStringSet(x_trim[i], start=dat_agg$min[i], end=dat_agg$max[i], width=NA, use.names=TRUE)
    print(c(i,dat_agg$sense[i]))
  } else {
    if (dat_agg$sense[i] == -1) {
      x_trim[i] <- reverseComplement(DNAStringSet(x_trim[i], start=dat_agg$min[i], end=dat_agg$max[i], width=NA, use.names=TRUE))
      print(c(i,dat_agg$sense[i]))
    } else {
      x_trim[i] <- DNAStringSet("NNNN") 
      print(c(i,"NNNN"))
    } }
}




#show outliers
mean_length <- mean(width(x_trim))
#Calculate a confidence interval
Conf_interv <- 2*sd(width(x_trim))
#show the sequences that will be removed
x_trim[ (width(x_trim) > mean_length + 2*Conf_interv | width(x_trim) < mean_length - 0.5*Conf_interv) , ]
#create the file without the outliers
xtrim_no_outliers <- x_trim[!(width(x_trim) > mean_length + 2*Conf_interv | 
                                width(x_trim) < mean_length - 0.5*Conf_interv), ]

# write fasta file of trimmed sequences
writeXStringSet(xtrim_no_outliers, file=paste(ID_Folder, "/GB_csv_extracted_", ID_fasta_files[k], ".fasta", sep=""), append=FALSE, format="fasta") 



cmd2 <- paste("cat ",  ID_Folder, "/", ID_fasta_files[k], " ", ID_Folder, "/GB_csv_extracted_", ID_fasta_files[k], ".fasta > ",  ID_Folder, "/to_align.fasta", sep="")

system(cmd2)




cmd3 <- paste("/opt/bio/mafft/bin/mafft --auto --reorder ", ID_Folder, "/to_align.fasta > ", ID_Folder, "/", ID_fasta_files[k], "_aligned.fasta", sep="")

system(cmd3)


alignment_file2 <- paste(ID_Folder, "/", ID_fasta_files[k], "_aligned.fasta", sep="")
  align <- read.dna(alignment_file2, format="fasta")
  # redo dm because of name change
  
#   lines <- numeric()
#   for(jj in 1:length(selected_strains)) {   
#     line_num <- as.numeric(grep(selected_strains[jj],rownames(align)))
#     lines <- c(lines,line_num)
#   }
#   
#   align_subset <- align[lines,]
 

# rownames(align) <- sub("NFIS", ">>>>>>>>>> NFIS",rownames(align), ignore.case = FALSE)
# rownames(align) <- sub(".seq", ".seq <<<<<<<<<<",rownames(align), ignore.case = FALSE)
rownames(align) <- sub("_AY598", " REFERENCE STRAIN_AY598",rownames(align), ignore.case = FALSE)

rownames(align) <- sub("checksum.", ">>>>>>>>>> checksum.<<<<<<<<<<<<",rownames(align), ignore.case = FALSE)
# rownames(align) <- sub("^UWI", ">>>>>>>>>> UWI",rownames(align), ignore.case = FALSE)
# rownames(align) <- sub("PF$", "PF <<<<<<<<<<",rownames(align), ignore.case = FALSE)

  
  # raw distance is p-distance with substitution -> d: transition + transversion
  dm <- dist.dna(align, model = "K80", pairwise.deletion = FALSE, as.matrix = TRUE)
  
  MaxV <- max(rowSums(dm))

  my_root <- which(rowSums(dm) == MaxV)

  dm <- dist.dna(align, model = "raw", pairwise.deletion = FALSE, as.matrix = TRUE)

  tree <- njs(dm)

write.tree(tree, file = paste(ID_Folder, "/", ID_fasta_files[k], "_tree.newick", sep=""), append = FALSE, digits = 10, tree.names = FALSE)

  
  # plot(max(unmatrix(dm, byrow = TRUE)), type="h")
  max(dm, na.rm = TRUE)
  nrow(dm)
  # write.table(dm, file = "dm.csv", append = FALSE, sep = ",", col.names = NA)
  
  pdf(file = paste(ID_Folder, "/", ID_fasta_files[k],"_NJ_bionj_K80_tree_GenBank.pdf", sep=""), width = 8, height =14 )
  # "0.5-((nrow(dm)-50)/500)" is a rough equation to remove 0.1 to cex factor for every 50 taxa to keep font size small enough for larger data
  # if you want to have longer branches, reduce x.lim by 0.1 increments
  #plot.phylo(type = "phylogram", root(tree, root[i], node = NULL, resolve.root = TRUE), font=1, cex = 0.5,  x.lim = 0.7)
  plot.phylo(type = "phylogram", root(tree, my_root[1], node = NULL, resolve.root = TRUE), font=1, 
             cex = 0.52 - (sqrt(nrow(dm))/70),  x.lim = 0.05 + max(dm, na.rm = TRUE)/1.3, edge.width= 1.1 - (nrow(dm)/1200), no.margin=TRUE)
  title(main="NJ", outer=FALSE, cex.main=1, font.main=2)
  dev.off()   
 
}

#   #write fasta file
#   write.dna(align_subset, file=paste(genes[ii], "/alignment_rev_names_codes_subset.fasta", sep = ""), format = "fasta")
#   
#   
# }
# 
# 
















#I performed  tblastn with "Whole-genome shotgun contigs" and "nt"

list.dirs()

# remove any data from previous analyses (not necessary when you restart R)
rm(list=ls())

GB_fasta_list <- list.files(path = ".", pattern = "\\.txt", recursive = TRUE)

i <- 1

# Read a dump from BLAST search of full length fasta sequences (has some complete chromosomes)
#  the loop is to do this for multiple files from different genes and put them together
GB_fasta <- DNAStringSet()
for(i in 1:length(GB_fasta_list)) {
  temp <- readDNAStringSet(GB_fasta_list[i])
  GB_fasta <- c(GB_fasta,temp)
}

# Order the file
GB_fasta <- GB_fasta[order(names(GB_fasta)), ] 

attributes(GB_fasta)


# This is a function to make copies of sequences that identical and that were combined into one by NCBI
# This is to create a vector for these lines with multiple fasta files
# gregexpr creates a list and and cbind puts the multiple occurrences of ">" together,  
# "-1" means a single unique sequence
# This is to duplicate the lines with multiple identical sequences

mult_hits <- cbind(gregexpr(">", names(GB_fasta), ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE))

new_GB_fasta <- DNAStringSet()

for(i in 1:length(GB_fasta)) {
  if (unlist(mult_hits[i])[1] != -1) {   
    parsed_names <- strsplit(names(GB_fasta[i]), ">")[[1]]  # [[1]] to show first element
    temp_x <- rep(GB_fasta[i], length(unlist(mult_hits[i])) + 1)
    names(temp_x) <- parsed_names   
    new_GB_fasta <- append(new_GB_fasta, temp_x) 
  } else {
    new_GB_fasta <-  append(new_GB_fasta, GB_fasta[i]) } 
}

# Order the file
new_GB_fasta <- new_GB_fasta[order(names(new_GB_fasta)), ] 


# write fasta file if you want
#  writeXStringSet(new_x, file="seq_dump.fasta", append=FALSE, format="fasta") 

# This second part is to fix the GenBank names

# This is to get the csv files
GB_Blast_csv_list  <- list.files(path = ".", pattern = "\\.csv", recursive = TRUE)

i <-  1
GB_Blast_table <-read.table(GB_Blast_csv_list[i],sep=",",header=FALSE)
# # Read csv Hit Table output from GenBank (written to do multiple outputs, use above if only one)
# GB_Blast_table <- data.frame()
# for(i in 1:length(GB_Blast_csv_list)) {
#   temp <-read.table(GB_Blast_csv_list[i],sep=",",header=FALSE)
#   GB_Blast_table <- rbind(GB_Blast_table,temp)
# }

# Put proper column names on BLASTGenBank output
colnames(GB_Blast_table) <- c("query id", "subject_ids", " %identity", "alignment length", "mismatches", "gap opens", "q.start", " q.end", "s.start", "s.end", "evalue", " bit_score")
# Order Hit table like the sequences
GB_Blast_table <- GB_Blast_table[order(GB_Blast_table$subject_ids), ] 




##############################
##      IF THERE ARE TWO gi NUMBERS
# This is to duplicate the lines with multiple identical sequences

mult_gi <- cbind(gregexpr(";", GB_Blast_table$subject_ids, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE))

new_GB_Blast_table <- data.frame()

#i <- 2

for(i in 1:nrow(GB_Blast_table)) {
  if (unlist(mult_gi[i])[1] != -1) {   
    parsed_names <- strsplit(as.character(GB_Blast_table$subject_ids[i]), ";")[[1]]  # [[1]] to show first element  
    temp_table <- do.call("rbind", replicate((length(unlist(mult_gi[i])) + 1), GB_Blast_table[i,], simplify = FALSE))
    temp_table$subject_ids <- parsed_names   
    new_GB_Blast_table <- rbind(new_GB_Blast_table, temp_table) 
  } else {
    new_GB_Blast_table <-  rbind(new_GB_Blast_table, GB_Blast_table[i,]) } 
}


# Creates a +/- column for orientation of sequences
new_GB_Blast_table$orientation <- sign(new_GB_Blast_table$s.end - new_GB_Blast_table$s.start)

# Put together data from multiple hits on single line
dat_agg <- aggregate(new_GB_Blast_table[c(13,9,10)], by=list(new_GB_Blast_table$subject_ids),FUN=c)

# find min and max and sense of hits
for(i in 1:nrow(dat_agg)) {
  dat_agg$min[i] <- min(unlist(c(dat_agg$s.start[i],dat_agg$s.end[i])))
  dat_agg$max[i] <- max(unlist(c(dat_agg$s.start[i],dat_agg$s.end[i])))
  dat_agg$sense[i] <- mean(unlist(c(dat_agg$orientation[i]))) 
}


# add a column with bp (width) of sequences.
dat_agg$bp <- width(new_GB_fasta)
dat_agg$width <- dat_agg$max - dat_agg$min
dat_agg$diff <- dat_agg$width - dat_agg$max


dat_agg <- dat_agg[order(as.character(dat_agg$Group.1)), ] 

# to make sure that name order are the same
check_order <- data.frame(names(new_GB_fasta), dat_agg[,1])

# Very Important
# THE ORDER OF THESE TWO COLUMNS SHOULD BE THE SAME
# SEQUENCES WILL BE OFF IF NOT
write.table(check_order, file = "CHECK ORDER OF FASTA file AND BLAST TABLE.csv", append = FALSE, sep = ",", col.names = NA)


summary(dat_agg$s.start)
summary(dat_agg$s.end)

length(new_GB_fasta[1])

# code without reverse complement
# x_trim <- DNAStringSet()
# for(i in 1:length(GB_fasta)) {
#   x_trim[i] <- DNAStringSet(GB_fasta[i], start=dat_agg$min[i], end=dat_agg$max[i], width=NA, use.names=FALSE)
# }


x_trim <- new_GB_fasta
for(i in 1:length(x_trim)) {
  if (dat_agg$sense[i] == 1) {
    x_trim[i] <- DNAStringSet(x_trim[i], start=dat_agg$min[i], end=dat_agg$max[i], width=NA, use.names=TRUE)
  } else {
    if (dat_agg$sense[i] == -1) {
      x_trim[i] <- reverseComplement(DNAStringSet(x_trim[i], start=dat_agg$min[i], end=dat_agg$max[i], width=NA, use.names=TRUE))
    } else {
      x_trim[i] <- "inconsistent orientation of the BLAST matches for at least one accession number" } }
}


######  I have done much better below to rename
# temp <- data.frame(matrix(unlist(strsplit(as.character(names(x_trim)), "\\|")), nrow= length(x_trim), byrow=TRUE))
# temp$X6 <- sub(" ", "|",temp$X5, ignore.case = FALSE)
# temp$X6 <- sub(" ", "|",temp$X6, ignore.case = FALSE)
# temp$X6 <- sub(" ", "|",temp$X6, ignore.case = FALSE)
# temp2 <- data.frame(matrix(unlist(strsplit(as.character(temp$X6), "\\|")), nrow= length(temp$X6), byrow=TRUE))
# temp2$X3 <- sub("sp.", "sp",temp2$X3, ignore.case = FALSE)
# temp2$X2 <- sub("\\.", "-",temp2$X2, ignore.case = FALSE)
# temp2$X2 <- gsub("_", "-",temp2$X2, ignore.case = FALSE)
# temp2$X2 <- sub(",", "",temp2$X2, ignore.case = FALSE)
# temp2$X2 <- gsub(":", "-",temp2$X2, ignore.case = FALSE)
# 
# names(x_trim) <- paste(temp2$X2,"_",temp2$X3,"_gi",temp$X2,sep="")


#show outliers
mean_length <- mean(width(x_trim))
#Calculate a confidence interval
Conf_interv <- 3*sd(width(x_trim))
#show the sequences that will be removed
x_trim[ (width(x_trim) > mean_length + Conf_interv | width(x_trim) < mean_length - Conf_interv) , ]
#create the file without the outliers
xtrim_no_outliers <- x_trim[!(width(x_trim) > mean_length + Conf_interv | 
                                width(x_trim) < mean_length - Conf_interv), ]



# write fasta file of trimmed sequences
writeXStringSet(xtrim_no_outliers, file="GB_csv_extracted.fasta", append=FALSE, format="fasta") 

#pipe( "blastall -p blastp -i text.fasta -d data.fasta" )

#pipe("mafft GB_csv_extracted.fasta GB_csv_extracted_aligned.fasta")






###################################################################################################
# This module is to fix the GenBank output


GB_names <- names(x_trim)

library("stringr")

# This is a loop to make some taxonomy glitches replacement to make the parsing with space accurate
   for(i in 1:length(GB_names)) {
     # This is to replace abbreviations like "R.secalis" in GenBank
       regexp1 <- " ([[:upper:]]{1})(\\.)([[:lower:]]+) "
       replacement1 <- str_extract(GB_names[i],regexp1)
       replacement1 <- sub("\\.", " ", replacement1, ignore.case = FALSE)
       GB_names[i] <- sub(" ([[:upper:]]{1})(\\.)([[:lower:]]+) ", replacement1, GB_names[i], ignore.case = FALSE)
    # This is to replace " f. sp. " for "_f._sp._"   
    #   string <- "Marssonina brunnea f. sp. multigermtubi MB_m1 tubulin beta chain (MBM_05801), mRNA"
       regexp2 <- " ([[:lower:]]+) f. sp. (\\'?)([[:lower:]]+)(\\'?) "
       replacement2 <- str_extract(GB_names[i],regexp2)
       replacement2 <- sub(" f. sp. ", "_.f._sp._", replacement2, ignore.case = FALSE)
       GB_names[i] <- sub(" ([[:lower:]]+) f. sp. (\\'?)([[:lower:]]+)(\\'?) ", replacement2, GB_names[i], ignore.case = FALSE)
   }
  

GB_names <- sub(" ", "|", GB_names, ignore.case = FALSE)
GB_names <- sub(" ", "|", GB_names, ignore.case = FALSE)
GB_names <- sub(" ", "|", GB_names, ignore.case = FALSE)

parsed_col <- data.frame(matrix(unlist(strsplit(as.character(GB_names), "\\|")), nrow=length(GB_names), byrow=T),stringsAsFactors = FALSE)

# This is a loop to make extract strain numbers.  There are two approaches, one when there is "strain|isolate|voucher" and one 
# where the strain number is right after species of f.sp. name.  This may need to be tweaked for different strain codes
  for(i in 1:length(parsed_col$X8)) {
      if  (grepl("strain|isolate|voucher", parsed_col$X8[i])) {   
        regexp <- "(strain|isolate|voucher) ([[:alpha:]]*)( ?|\\-*|\\_*)([[:digit:]]*|[[:upper:]]*)(\\.*|\\-*)([[:digit:]]*|[[:upper:]]*)"
        parsed_col$strain[i] <- str_extract(parsed_col$X8[i],regexp)
      } else {
        # the trick for this was to leave a blank when there is no strain number
        regexp <- "([[:upper:]]*)( ?|\\-*|\\_*)([[:digit:]]*|[[:upper:]]*)(\\.*|\\-*)([[:digit:]]*|([[:alpha:]]+[[:digit:]]+))"
        parsed_col$strain[i] <- str_extract(parsed_col$X8[i],regexp) } 
    }

# remove strain|isolate|voucher in strain codes
parsed_col$strain2 <- sub("strain |isolate |voucher ", "", parsed_col$strain, ignore.case = FALSE)
parsed_col$strain2 <- sub("_", "", parsed_col$strain2, ignore.case = FALSE)
parsed_col$strain2 <- sub(" ", "", parsed_col$strain2, ignore.case = FALSE)
parsed_col$GB <- sub("\\.[[:digit:]]", "", parsed_col$X4, ignore.case = FALSE)
parsed_col$GB <- sub("_", "", parsed_col$GB, ignore.case = FALSE)
parsed_col$species <- gsub("_", "", parsed_col$X7, ignore.case = FALSE)
#for_names <- paste(parsed_col$X6,"_",parsed_col$species,"_",parsed_col$GB,"_strain(",parsed_col$strain2,")", sep = "")
names(x_trim) <- paste(parsed_col$X6,"_",parsed_col$species,"_",parsed_col$GB,"_strain(",parsed_col$strain2,")", sep = "")

writeXStringSet(x_trim, file="new_names.fasta", append=FALSE, format="fasta") 

