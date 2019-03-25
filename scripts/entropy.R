library(GrpString)
library(tibble)
library(ggplot2)

file_list <- list.files(path="../family_fasta", full.names = TRUE)

group_entropy_data <- c()
individual_entropy_data <- c()

for (i in 1:length(file_list)){
  family <- basename(file_list[i]) 
  lines <- readLines(file_list[i])
  seqs <- lines[seq(2, length(lines),2)] 
  group_entropy <- GrpString::TransEntro(seqs)
  n <- length(seqs)
  group_entropy_data <- rbind(group_entropy_data, (c(family, group_entropy, n)))
  
  all_entropy <- GrpString::TransEntropy(seqs)
  all_entropy <- cbind(family, all_entropy$Entropy)
  individual_entropy_data <- rbind(individual_entropy_data, all_entropy) 
}

group_entropy_data <- tibble("family" = group_entropy_data[,1], 
         "group_entropy" = group_entropy_data[,2],
         "n" = group_entropy_data[,3])
class(group_entropy_data$group_entropy) <- 'numeric'
class(group_entropy_data$n) <- 'numeric'

individual_entropy_data <- tibble("family" = individual_entropy_data[,1],
                                  "indiv_entropy" = individual_entropy_data[,2])
class(individual_entropy_data$indiv_entropy) <- 'numeric'
#plotting time