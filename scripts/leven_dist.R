library(GrpString)
library(tibble)
library(ggplot2)
library(dplyr)

file_list <- list.files(path="../family_fasta", full.names = TRUE)

distance_data <- c()

for (i in 1:length(file_list)){
  family <- basename(file_list[i]) 
  lines <- readLines(file_list[i])
  seqs <- lines[seq(2, length(lines),2)] 
  n <- length(seqs)
  
  if(n > 1){
  
    distances <- c()
    combinations <- combn(seqs, 2)
    for (j in 1:dim(combinations)[2]){
      dist <- adist(combinations[1,j], combinations[2,j])
      distances <- c(distances, dist)
    }
  }
  else {
    distances = c(0)
  }
  seq_dists <- cbind(family, distances, n)
  distance_data <- rbind(distance_data, seq_dists) 
}

distances <- tibble("family" = distance_data[,1], 
                             "ldist" = distance_data[,2],
                             "n" = distance_data[,3])
class(distances$ldist) <- "numeric"
 

max_dists <- distances %>% group_by(family) %>% summarise(max = max(ldist)) %>% arrange(max)
ggplot(max_dists, aes(x=max)) + 
  geom_density() + xlab('Edit Distances') + 
  ggtitle('Kernel Density Estimate of Max Edit Distances by Family')

no_zero <- max_dists %>% filter(max > 0)
tail(no_zero, 20) %>% ggplot(aes(y=max, x=reorder(family, -max))) + geom_bar(stat='identity') + coord_flip()

distances %>% ggplot(aes(x=reorder(family, ldist, FUN=median), y=ldist)) + geom_boxplot() + theme_void() 


max_family <- tail(no_zero, 20) %>% select(family)


distances %>% ggplot(aes(x=reorder(family, -ldist, FUN=median), y=ldist)) + 
  geom_boxplot() + 
  xlab("AMR Family") + ylab("Aitchinson Distance") + t
  theme_light() + 
    theme(axis.text.x = element_blank()) + 
    ggsave('distances.png')

  
medians <- distances %>% group_by(family) %>% summarise(med = median(ldist)) %>% arrange(-med)
top_median <- head(medians, 20) %>% select(family)
top_medians <- distances %>% filter(family %in% top_median$family) 
top_medians %>% ggplot(aes(x=reorder(family, ldist, FUN=median), y=ldist)) + geom_boxplot() + coord_flip()
