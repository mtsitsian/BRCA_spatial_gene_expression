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


head(unify(c_CODs$healthy))
cods <- lapply(c_CODs, unify)
summary(cods)
head(cods$healthy)




for (group in names(cods)){
  
  write.table(format(cods[[group]], 
                     scientific = F,
                     trim = T), file = paste0("../bp-metric/DCEs/DCEs_", group, ".tsv"),
              sep = "\t", 
              quote = F, row.names = F)
}

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

for (i in 1:4) {
  cods[[i]] %<>% cbind(cod_binsignal[[i]])
  colnames(cods[[i]])[4] <- "mean_binsignal"
}

for (group in names(cods)){
  
  write.table(format(cods[[group]], 
                     scientific = F,
                     trim = T), file = paste0("RData/DCEs/DCEs_", group), sep = "\t", 
              quote = F, dec = '.', row.names = F)
}
