###################################
### Loading necessary libraries ###
###################################

library(tidyr)
library(plyr)
library(dplyr)

#################################
### Setting working directory ###
#################################

setwd("path/to/working/directory")

#########################################
### Loading the output from script 01 ###
#########################################

positive_exon_info <- read.table(
  "annotation_file.exons.targets.positive.2.3.4.gff3",
  sep="\t", header=FALSE)
negative_exon_info <- read.table(
  "annotation_file.exons.targets.negative.2.3.4.gff3",
  sep="\t", header=FALSE)
target_loci<- read.table("loci_file.2.3.4.gff3",sep="\t",header=FALSE)

###############################################
### Translating loci on the positive strand ###
###############################################

# Tidying up the data #
#######################

colnames(positive_exon_info)[colnames(positive_exon_info) == "V1"] <- "chr"
positive_exon_info <- subset(positive_exon_info, select = -c(V2,V3,V6,V7,V8) )
colnames(positive_exon_info)[colnames(positive_exon_info) == "V4"] <- "start"
colnames(positive_exon_info)[colnames(positive_exon_info) == "V5"] <- "end"
colnames(positive_exon_info)[colnames(positive_exon_info) == "V9"] <- "transcript_id"
colnames(positive_exon_info)[colnames(positive_exon_info) == "V10"] <- "id"
colnames(positive_exon_info)[colnames(positive_exon_info) == "V11"] <- "length"

# The code below is specific to TargetFinder output. Please edit this code so that transcript_loci is a dataframe with 5 columns:
# (1) TS, (2) TE, (3) locus_id, (4) Tlength, and (5) transcript_id. These correspond to 
# (1) Transcript-relative locus start site, (2) transcript-relative locus end site, (3) locus ID, (4) locus length, and (5) transctipt ID. 

transcript_loci <- subset(transcript_loci, select = -c(V1,V2,V3,V4,V7,V8,V9) )
colnames(transcript_loci)[colnames(transcript_loci) == "V12"] <- "transcript_id"
transcript_loci <- separate(transcript_loci,col=V10,into=c('col1','locus_id','col2','col3','col4'),sep='=')
transcript_loci <- subset(transcript_loci, select = -c(col1,col2,col3,col4) )
transcript_loci <- separate(transcript_loci,col=locus_id,into=c('locus_id','col1'),sep=';')
transcript_loci <- subset(transcript_loci, select = -c(col1) )
colnames(transcript_loci)[colnames(transcript_loci) == "V5"] <- "TS"
colnames(transcript_loci)[colnames(transcript_loci) == "V6"] <- "TE"
colnames(transcript_loci)[colnames(transcript_loci) == "V11"] <- "Tlength"

# Converting the exon information into wide format #
####################################################

positive_exon_info <- positive_exon_info %>%
  mutate(chr_transcript_id = paste0(chr, ",", transcript_id))
positive_exon_info_wide <- positive_exon_info %>%
  pivot_wider(
    id_cols = chr_transcript_id,
    names_from = id,
    values_from = c(start, end, length)
  )

positive_exon_info_wide <- separate(positive_exon_info_wide,col=chr_transcript_id,into=c('chr','transcript_id'),sep=',')

positive_exon_info_wide <- separate(positive_exon_info_wide,col=transcript_id,into=c('col1','transcript_id'),sep='=')
positive_exon_info_wide <- subset(positive_exon_info_wide, select = -c(col1) )
positive_exon_info_wide$transcript_id <- trimws(positive_exon_info_wide$transcript_id)

# Merging transcript and target info #
######################################

positive_targetexon_info <- merge(transcript_loci,positive_exon_info_wide,by="transcript_id", all.x = FALSE, all.y = FALSE)

# Translating the loci #
########################

positive_translatingloci <- function(row) {
  
  # Splitting the dataframe into rows, so each row is worked on independently
  row  <- as.data.frame(row)
  
  # Setting up the counter for unassigned target nucleotides. The number refers to the number of unassigned nucleotides remaining for the row. 
  unassigned_nt <- row$Tlength
  
  # Generating the numeric variables needed for the for the function
  exon_no_counter <- 1 # Variable for the exon number we are currently looking at
  next_exon_counter <- 2 # Variable for the next exon number from the one we are currently looking at
  region_no_counter <- 1 # Variable for the target region we are currently assigning
  current_transcript_loc_counter <- row$TS - row$length_exon_1 # Numeric variable for each TS-E1-E2 etc. value
  previous_transcript_loc_counter <- row$TS # Numeric variable for each TS-E1-E2 etc. value, but from the previous exon. I.e if we are on exon 3, this will be TS-E1-E2.
  
  while (unassigned_nt != 0) { # The script will run until the number of unassigned nucleotides is zero.  
    
    # Generating character variables I need for for loop 1
    region_start <- paste0("R", region_no_counter, "_start") # Generating the text to label the region start variable.
    region_end <- paste0("R", region_no_counter, "_end") # Generating the text to label the region end variable.
    exon_start_var <- paste0("start_exon_", exon_no_counter) # Generating the text to label the exon start variable.
    exon_end_var <- paste0("end_exon_", exon_no_counter) # Generating the text to label the exon end variable.
    exon_length_var <- paste0("length_exon_",exon_no_counter) # Generating the text to label the exon length variable.
    next_exon_start_var <- paste0("start_exon_", next_exon_counter) # Generating the text to label the exon start variable.
    next_exon_length_var <- paste0("length_exon_",next_exon_counter) # Generating the text to label the next exon length variable.
    
    if (current_transcript_loc_counter <0){ # Does the target start in this transcript? True if it does. Ifelse statement 1.
      
      ##OUTPUT## Calculating target start site in this exon.
      row[,ncol(row)+1] <- 0 # Make a new column
      names(row)[ncol(row)] <- region_start # Rename this new column
      row[,ncol(row)] <- ((row[[exon_start_var]] + previous_transcript_loc_counter) - 1) # Insert the value for this new column
      
      if (current_transcript_loc_counter > ((row$Tlength -1)*-1)){ # Is TS - E(n) > the negative miRNA length? This is true if the target site extends beyond this exon. Ifelse statement 2.
        
        ##OUTPUT## Calculating the end of this region (the end of the exon)
        row[,ncol(row)+1] <- 0 # Make a new column
        names(row)[ncol(row)] <- region_end # Rename this new column
        row[,ncol(row)] <- row[[exon_end_var]] # Insert the value for this new column (just the end of the exon).
        
        unassigned_nt <- (unassigned_nt - ( row[[region_end]] - row[[region_start]] +1)) # Updating the counter to show how many unassigned nucleotides are left.
        
        while (unassigned_nt != 0) { # While loop 2
          
          # Updating numeric variables for this for loop
          exon_no_counter <- exon_no_counter + 1 # Increasing the counter by one so the next loop will look at the next exon.
          region_no_counter <- region_no_counter + 1 # Increasing the counter by one so the next loop will label the region correctly. 
          
          # Update character variables for this for loop
          region_start <- paste0("R", region_no_counter, "_start") # Generating the text to label the region start variable.
          region_end <- paste0("R", region_no_counter, "_end") # Generating the text to label the region end variable.
          exon_start_var <- paste0("start_exon_", exon_no_counter) # Generating the text to label the exon start variable.
          exon_end_var <- paste0("end_exon_", exon_no_counter) # Generating the text to label the exon end variable.
          exon_length_var <- paste0("length_exon_",exon_no_counter) # Generating the text to label the exon length variable.
          
          ##OUTPUT## Calculating the start of this region (the start of the exon).
          row[,ncol(row)+1] <- 0 # Make a new column.
          names(row)[ncol(row)] <- region_start # Rename this new column.
          row[,ncol(row)] <- row[[exon_start_var]] # Insert the value for this new column (just the start of the exon).
          
          if (unassigned_nt > row[[exon_length_var]]) { # Are the remaining nulceotides to be assigned more than the length of this exon? This is true if the target site extends beyond this exon. Ifelse statement 3. 
            
            ##OUTPUT## Calculating the end of this region (the end of the exon).
            row[,ncol(row)+1] <- 0 # Make a new column
            names(row)[ncol(row)] <- region_end # Rename this new column
            row[,ncol(row)] <- row[[exon_end_var]] # Insert the value for this new column (just the end of the exon).
            
            unassigned_nt <- (unassigned_nt - ( row[[region_end]] - row[[region_start]] +1)) # Updating the counter to show how many unassigned nucleotides are left.
            
          } else { # Else condition for ifelse statement 3 - for if the region ends in this exon. 
            
            ##OUTPUT## Finding the target end site if it's within this exon.
            row[,ncol(row)+1] <- 0 # Make a new column
            names(row)[ncol(row)] <- region_end # Rename this new column
            row[,ncol(row)] <- ((row[[exon_start_var]] + unassigned_nt) -1) # Insert the value for this new column
            
            unassigned_nt <- (unassigned_nt - ( row[[region_end]] - row[[region_start]] +1)) # Updating the counter to show how many unassigned nucleotides are left.
            
          } # Closing bracket for ifelse statement 3
          
        } # Closing bracket for for while loop 2
        
      } else { # Else condition for ifelse statement 2
        
        ##OUTPUT## Finding the target end site if it's within this exon.
        row[,ncol(row)+1] <- 0 # Make a new column
        names(row)[ncol(row)] <- region_end # Rename this new column
        row[,ncol(row)] <- ((row[[region_start]] + row$Tlength) -1) # Insert the value for this new column
        
        unassigned_nt <- (unassigned_nt - ( row[[region_end]] - row[[region_start]] +1)) # Updating the counter to show how many unassigned nucleotides are left.
        
      } # Closing bracket for ifelse statement 2
      
    } else { # Else condition for ifelse statement 1
      
      exon_no_counter <- exon_no_counter + 1 # Variable for the exon number we are currently looking at
      next_exon_counter <- next_exon_counter + 1 # Variable for the next exon number from the one we are currently looking at
      exon_length_var <- paste0("length_exon_",exon_no_counter) # Generating the text to label the exon length variable.
      previous_transcript_loc_counter <- current_transcript_loc_counter
      current_transcript_loc_counter <- current_transcript_loc_counter - row[[exon_length_var]] # Numeric variable for each TS-E1-E2 etc. value
      
    } # Closing bracket for ifelse statement 1
    
  } # Closing bracket for while statement
  
  
  return(row)
  
} # Closing bracket for function

positive_list_total <- list()

for (i in 1:nrow(positive_targetexon_info)){
  result <- positive_translatingloci(positive_targetexon_info[i, ])
  positive_list_total <- append(positive_list_total, list(result))
}

positive_loci <- bind_rows(positive_list_total, .id = "source")
positive_loci$strand <- "+"

###############################################
### Translating loci on the negative strand ###
###############################################

# Tidying up the data #
#######################

colnames(negative_exon_info)[colnames(negative_exon_info) == "V1"] <- "chr"
negative_exon_info <- subset(negative_exon_info, select = -c(V2,V3,V6,V7,V8) )
colnames(negative_exon_info)[colnames(negative_exon_info) == "V4"] <- "start"
colnames(negative_exon_info)[colnames(negative_exon_info) == "V5"] <- "end"
colnames(negative_exon_info)[colnames(negative_exon_info) == "V9"] <- "transcript_id"
colnames(negative_exon_info)[colnames(negative_exon_info) == "V10"] <- "id"
colnames(negative_exon_info)[colnames(negative_exon_info) == "V11"] <- "length"

# Converting the exon information into wide format #
####################################################

negative_exon_info <- negative_exon_info %>%
  mutate(chr_transcript_id = paste0(chr, ",", transcript_id))
negative_exon_info_wide <- negative_exon_info %>%
  pivot_wider(
    id_cols = chr_transcript_id,
    names_from = id,
    values_from = c(start, end, length)
  )

negative_exon_info_wide <- separate(negative_exon_info_wide,col=chr_transcript_id,into=c('chr','transcript_id'),sep=',')

negative_exon_info_wide <- separate(negative_exon_info_wide,col=transcript_id,into=c('col1','transcript_id'),sep='=')
negative_exon_info_wide <- subset(negative_exon_info_wide, select = -c(col1) )
negative_exon_info_wide$transcript_id <- trimws(negative_exon_info_wide$transcript_id)

# Merging transcript and target info #
######################################

negative_targetexon_info <- merge(transcript_loci,negative_exon_info_wide,by="transcript_id", all.x = FALSE, all.y = FALSE)

# Translating the loci #
########################

negative_translatingloci <- function(row) {
  
  # Splitting the dataframe into rows, so each row is worked on independently
  row  <- as.data.frame(row)
  
  # Setting up the counter for unassigned target nucleotides. The number refers to the number of unassigned nucleotides remaining for the row. 
  unassigned_nt <- row$Tlength
  
  # Generating numeric variables I need for the for the function
  exon_no_counter <- 1 # Variable for the exon number we are currently looking at
  next_exon_counter <- 2 # Variable for the next exon number from the one we are currently looking at
  region_no_counter <- 1 # Variable for the target region we are currently assigning
  current_transcript_loc_counter <- row$TS - row$length_exon_1 # Numeric variable for each TS-E1-E2 etc. value
  previous_transcript_loc_counter <- row$TS # Numeric variable for each TS-E1-E2 etc. value, but from the previous exon. I.e if we are on exon 3, this will be TS-E1-E2.
  
  while (unassigned_nt != 0) { # The script will run until the number of unassigned nucleotides is zero.  
    
    # Generating character variables I need for for loop 1
    region_start <- paste0("R", region_no_counter, "_start") # Generating the text to label the region start variable.
    region_end <- paste0("R", region_no_counter, "_end") # Generating the text to label the region end variable.
    exon_start_var <- paste0("start_exon_", exon_no_counter) # Generating the text to label the exon start variable.
    exon_end_var <- paste0("end_exon_", exon_no_counter) # Generating the text to label the exon end variable.
    exon_length_var <- paste0("length_exon_",exon_no_counter) # Generating the text to label the exon length variable.
    next_exon_start_var <- paste0("start_exon_", next_exon_counter) # Generating the text to label the exon start variable.
    next_exon_length_var <- paste0("length_exon_",next_exon_counter) # Generating the text to label the next exon length variable.
    
    if (current_transcript_loc_counter <0){ # Does the target start in this transcript? True if it does. Ifelse statement 1.
      
      ##OUTPUT## Calculating target start site in this exon.
      row[,ncol(row)+1] <- 0 # Make a new column
      names(row)[ncol(row)] <- region_end # Rename this new column
      row[,ncol(row)] <- ((row[[exon_end_var]] - previous_transcript_loc_counter) + 1) # Insert the value for this new column
      
      if (current_transcript_loc_counter > ((row$Tlength -1)*-1)){ # Is TS - E(n) > the negative miRNA length? This is true if the target site extends beyond this exon. Ifelse statement 2.
        
        ##OUTPUT## Calculating the end of this region (the end of the exon)
        row[,ncol(row)+1] <- 0 # Make a new column
        names(row)[ncol(row)] <- region_start # Rename this new column
        row[,ncol(row)] <- row[[exon_start_var]] # Insert the value for this new column (just the end of the exon).
        
        unassigned_nt <- (unassigned_nt - ( row[[region_end]] - row[[region_start]] +1)) # Updating the counter to show how many unassigned nucleotides are left.
        
        while (unassigned_nt != 0) { # While loop 2
          
          # Updating numeric variables for this for loop
          exon_no_counter <- exon_no_counter + 1 # Increasing the counter by one so the next loop will look at the next exon.
          region_no_counter <- region_no_counter + 1 # Increasing the counter by one so the next loop will label the region correctly. 
          
          # Update character variables for this for loop
          region_start <- paste0("R", region_no_counter, "_start") # Generating the text to label the region start variable.
          region_end <- paste0("R", region_no_counter, "_end") # Generating the text to label the region end variable.
          exon_start_var <- paste0("start_exon_", exon_no_counter) # Generating the text to label the exon start variable.
          exon_end_var <- paste0("end_exon_", exon_no_counter) # Generating the text to label the exon end variable.
          exon_length_var <- paste0("length_exon_",exon_no_counter) # Generating the text to label the exon length variable.
          
          ##OUTPUT## Calculating the start of this region (the start of the exon).
          row[,ncol(row)+1] <- 0 # Make a new column.
          names(row)[ncol(row)] <- region_end # Rename this new column.
          row[,ncol(row)] <- row[[exon_end_var]] # Insert the value for this new column (just the start of the exon).
          
          if (unassigned_nt > row[[exon_length_var]]) { # Are the remaining nulceotides to be assigned more than the length of this exon? This is true if the target site extends beyond this exon. Ifelse statement 3. 
            
            ##OUTPUT## Calculating the end of this region (the end of the exon).
            row[,ncol(row)+1] <- 0 # Make a new column
            names(row)[ncol(row)] <- region_start # Rename this new column
            row[,ncol(row)] <- row[[exon_start_var]] # Insert the value for this new column (just the end of the exon).
            
            unassigned_nt <- (unassigned_nt - ( row[[region_end]] - row[[region_start]] +1)) # Updating the counter to show how many unassigned nucleotides are left.
            
          } else { # Else condition for ifelse statement 3 - for if the region ends in this exon. 
            
            ##OUTPUT## Finding the target end site if it's within this exon.
            row[,ncol(row)+1] <- 0 # Make a new column
            names(row)[ncol(row)] <- region_start # Rename this new column
            row[,ncol(row)] <- ((row[[exon_end_var]] - unassigned_nt) +1) # Insert the value for this new column
            
            unassigned_nt <- (unassigned_nt - ( row[[region_end]] - row[[region_start]] +1)) # Updating the counter to show how many unassigned nucleotides are left.
            
          } # Closing bracket for ifelse statement 3
          
        } # Closing bracket for for while loop 2
        
      } else { # Else condition for ifelse statement 2
        
        ##OUTPUT## Finding the target end site if it's within this exon.
        row[,ncol(row)+1] <- 0 # Make a new column
        names(row)[ncol(row)] <- region_start # Rename this new column
        row[,ncol(row)] <- ((row[[region_end]] - row$Tlength) +1) # Insert the value for this new column
        
        unassigned_nt <- (unassigned_nt - ( row[[region_end]] - row[[region_start]] +1)) # Updating the counter to show how many unassigned nucleotides are left.
        
      } # Closing bracket for ifelse statement 2
      
    } else { # Else condition for ifelse statement 1
      
      exon_no_counter <- exon_no_counter + 1 # Variable for the exon number we are currently looking at
      next_exon_counter <- next_exon_counter + 1 # Variable for the next exon number from the one we are currently looking at
      exon_length_var <- paste0("length_exon_",exon_no_counter) # Generating the text to label the exon length variable.
      previous_transcript_loc_counter <- current_transcript_loc_counter
      current_transcript_loc_counter <- current_transcript_loc_counter - row[[exon_length_var]] # Numeric variable for each TS-E1-E2 etc. value
      
    } # Closing bracket for ifelse statement 1
    
  } # Closing bracket for while statement
  
  
  return(row)
  
} # Closing bracket for function

negative_list_total <- list()

for (i in 1:nrow(negative_targetexon_info)){
  result <- negative_translatingloci(negative_targetexon_info[i, ])
  negative_list_total <- append(negative_list_total, list(result))
}

negative_loci <- bind_rows(negative_list_total, .id = "source")
negative_loci$strand <- "-"

#############################################################
### Combining the loci from positive and negative strands ###
#############################################################

genomic_loci <- rbind.fill(negative_loci, positive_loci)
df_for_export <- genomic_loci[, c("locus_id", "chr", "strand", grep("^R", names(genomic_loci), value = TRUE))]

write.table(df_for_export, file = "genomic_loci.txt", row.names = FALSE, col.names = TRUE,quote = FALSE)
