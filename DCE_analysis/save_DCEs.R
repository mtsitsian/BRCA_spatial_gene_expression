library(tidyverse) #data wrangling
library(magrittr) #pipeline control

dir.create("data")
dir.create("data/DCE")

#save DCEs

unify <- function(x){
  chromosome <- integer()
  start <- integer()
  end <- integer()
  for (chr in names(x)){
    chr_n <- chr %>%
      strsplit("chr") %>%
      unlist %>% 
      .[2]
    
    chromosome <- chr_n %>%
      rep(nrow(x[[chr]])) %>%
      append(x = chromosome)
    
    start %<>% c(x[[chr]][,"START"])
    end %<>% c(x[[chr]][,"END"])
    
  }
  data.frame(chromosome, start, end)
}

unify <- function(x){
  chromosome <- integer()
  start <- integer()
  end <- integer()
  for (chr in names(x)){
    
    chromosome <- chr %>%
      rep(nrow(x[[chr]])) %>%
      append(x = chromosome)
    
    start %<>% c(x[[chr]][,"START"])
    end %<>% c(x[[chr]][,"END"])
    
  }
  data.frame(chromosome, start, end)
}

#use the control group as an example (here, the control group is NP)
head(unify(c_CODs$pretreat))
cods <- lapply(c_CODs, unify)
summary(cods)
head(cods$cy6)




#for (group in names(cods)){

#  write.table(format(cods[[group]], 
#                     scientific = F,
#                     trim = T), file = paste0("DCEs_", group, ".tsv"),
#              sep = "\t", 
#              quote = F, row.names = F)
#}

find_mean_binsignal <- function(x){
  mean_binsignal <- list()
  for (group in names(x)){
    mean_binsignal[[group]] <- numeric()
    for (chr in names(x[[group]])){
      
      start <- x[[group]][[chr]][,"START"]
      end <- x[[group]][[chr]][,"END"]
      
      for (i in 1:length(start)){
        mean_cod_binsignal <- bin_signals[[group]][[chr]][start[i]:end[i], "binsignal"] %>% mean
        mean_binsignal[[group]] %<>% c(mean_cod_binsignal)
      }
    }
  }
  mean_binsignal
}

cod_binsignal <- find_mean_binsignal(CODs)

#i is equal to the conditions of the study (healthy, subtype1, subtype2 etc)
for (i in 1:6) {
  cods[[i]] %<>% cbind(cod_binsignal[[i]])
  colnames(cods[[i]])[4] <- "mean_binsignal"
}

for (group in names(cods)){
  
  write.table(format(cods[[group]], 
                     scientific = F,
                     trim = T), file = paste0("data/DCE/DCEs_", group), sep = "\t", 
              quote = F, dec = '.', row.names = F)
}
