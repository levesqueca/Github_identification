library(ape)
library(openxlsx)

packageVersion("ape")

sharedPath   <- "/isilon/biodiversity/users/shared/Phytophthora_ID/"
ID_Folder <- paste(sharedPath,"Alejandro_CO1",sep="")
setwd(sharedPath)
getwd()

ID_fasta_files <- "Nucleotide_alignment_COI_24_oct_2016.fasta"

# read original names 
my_sequences <- read.dna(paste(ID_Folder, "/", ID_fasta_files, sep =""), format="fasta")

#ROWNAMES
temp1  <- sub("_","$",rownames(my_sequences))
parsed_col <- data.frame(matrix(unlist(strsplit(as.character(temp1), "\\$")), nrow=length(temp1), byrow=T),stringsAsFactors = FALSE)
parsed_col$X1 <- sub("LEV","",parsed_col$X1)
parsed_col$X1 <- sub("A","",parsed_col$X1)
rownames(my_sequences) <- parsed_col$X1
# Order the file
my_sequences <- my_sequences[order(rownames(my_sequences)), ] 


metadata_file2  <- "CTB_Specimen_Bulk_Edit_XSSF.xlsx"
metadata2 <- read.xlsx(metadata_file2, sheet = "Specimen Data", colNames = TRUE,
                       rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                       rows= c(1,3:62), check.names = FALSE, namedRegion = NULL)

colnames(metadata2) <- sub("\\/","_",colnames(metadata2)) 

# to pull out the metadata for those that have been sequenced
metadata3 <- metadata2[metadata2$Specimen.Identifier %in% rownames(my_sequences), ]

# sort on the specimen identifier
metadata3 <- metadata2[metadata2$Specimen.Identifier %in% rownames(my_sequences), ]
metadata3 <- metadata3[order(metadata3$Specimen.Identifier), ] 


#  check is the names of sequences are matching the order of the new mames
name_check <- data.frame(rownames(my_sequences), metadata3$Specimen.Identifier)


#NAMES
temp1  <- sub("_","$",names(my_sequences))
parsed_col <- data.frame(matrix(unlist(strsplit(as.character(temp1), "\\$")), nrow=length(temp1), byrow=T),stringsAsFactors = FALSE)
parsed_col$X1 <- sub("LEV","",parsed_col$X1)
names(my_sequences) <- parsed_col$X1
# Order the file
my_sequences <- my_sequences[order(names(my_sequences))] 


metadata_file2  <- "CTB_Specimen_Bulk_Edit_XSSF.xlsx"
metadata2 <- read.xlsx(metadata_file2, sheet = "Specimen Data", colNames = TRUE,
          rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
           rows= c(1,3:62), check.names = FALSE, namedRegion = NULL)

colnames(metadata2) <- sub("\\/","_",colnames(metadata2)) 

# to pull out the metadata for those that have been sequenced
metadata3 <- metadata2[metadata2$Specimen.Identifier %in% names(my_sequences), ]

# sort on the specimen identifier
metadata3 <- metadata2[metadata2$Specimen.Identifier %in% names(my_sequences), ]
metadata3 <- metadata3[order(metadata3$Specimen.Identifier), ] 


#  check is the names of sequences are matching the order of the new mames
name_check <- data.frame(names(my_sequences), metadata3$Specimen.Identifier)

##################
#start again


new_names <- paste(metadata3$Genus, "_", metadata3$Species, "|GB|",   metadata3$Collection.Code,	metadata3$Specimen.Identifier,
                   "|",   metadata3$Associated.Organism.Genus,"_",	metadata3$Associated.Organism.Species,"|",
                   metadata3$Country, "_", 
  #                 metadata3$Province_State, "_", metadata3$Region, "_",	
                   metadata3$City	, "_",
 #                  metadata3$Exact.Site,
                   sep="")

new_names <- gsub(" ","\\_",new_names)
new_names <- gsub("\\-","\\_",new_names)
new_names <- gsub("__","_",new_names)
new_names <- gsub("__","_",new_names)
new_names <- gsub("\\_$","",new_names)
new_names

#write the new name to the sequences
rownames(my_sequences) <- new_names
#names(my_sequences) <- new_names



# export fasta file
write.dna(my_sequences, paste(ID_Folder, "/my_sequences_with_metadata.fasta", sep=""), format = "fasta", append = FALSE)

