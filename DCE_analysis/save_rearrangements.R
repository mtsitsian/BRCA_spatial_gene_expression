#This is part of the former "meta_analysis.R" script
#save different categories of cods along with their mean binsignal for further analysis
#This script is intended to save the rearrangements files in your working directory
#IMPORTANT!!!! The ctrl_grp is the control group of your study! Set it at the beginning of the script right below

ctrl_grp <- "NP" 

library(tidyverse) #data wrangling
find_mean_binsignal <- function(group, chr, start, end){
  
  bin_signals[[group]][[chr]][start:end, "binsignal"] %>% mean
  
}


depleted <- list()
emerged <- list()
intact <- list()
split <- list()
merged <- list()
border_shift <- list()
expanded <- list()
contracted <- list()
shift <- list()

for (group in names(rearrangements)){
  depleted[[group]] <- data.frame(chromosome = character(), start = integer(), end  = integer(), mean_binsignal = numeric())
  emerged[[group]] <- data.frame(chromosome = character(), start = integer(), end  = integer(), mean_binsignal = numeric())
  intact[[group]] <- data.frame(chromosome = character(), start = integer(), end  = integer(), mean_binsignal = numeric())
  # split[[group]] <- data.frame(chromosome = character(), start = integer(), end  = integer(), mean_binsignal = numeric())
  # merged[[group]] <- data.frame(chromosome = character(), start = integer(), end  = integer(), mean_binsignal = numeric())
  split[[group]] <- data.frame(chromosome = character(), start_H = integer(), end_H  = integer(), mean_binsignal_H = numeric(),
                               start_D = integer(), end_D  = integer(), mean_binsignal_D = numeric(), Jaccard_index = numeric())
  merged[[group]] <- data.frame(chromosome = character(), start_H = integer(), end_H  = integer(), mean_binsignal_H = numeric(),
                                start_D = integer(), end_D  = integer(), mean_binsignal_D = numeric(), Jaccard_index = numeric())
  border_shift[[group]] <- data.frame(chromosome = character(), start_H = integer(), end_H  = integer(), mean_binsignal_H = numeric(),
                                      start_D = integer(), end_D  = integer(), mean_binsignal_D = numeric(), Jaccard_index = numeric(), 
                                      border_altered = character())
  expanded[[group]] <- data.frame(chromosome = character(), start_H = integer(), end_H  = integer(), mean_binsignal_H = numeric(),
                                  start_D = integer(), end_D  = integer(), mean_binsignal_D = numeric(), Jaccard_index = numeric())
  contracted[[group]] <- data.frame(chromosome = character(), start_H = integer(), end_H  = integer(), mean_binsignal_H = numeric(),
                                    start_D = integer(), end_D  = integer(), mean_binsignal_D = numeric(), Jaccard_index = numeric())
  shift[[group]] <- data.frame(chromosome = character(), start_H = integer(), end_H  = integer(), mean_binsignal_H = numeric(),
                               start_D = integer(), end_D  = integer(), mean_binsignal_D = numeric(), Jaccard_index = numeric())
  
  
  
  for (chr in names(rearrangements[[1]])){
    
    n_depleted <- rearrangements[[group]][[chr]][["depleted"]]
    n_emerged <- rearrangements[[group]][[chr]][["emerged"]]
    n_intact <- rearrangements[[group]][[chr]][["intact"]]
    # n_split <- unique(rearrangements[[group]][[chr]][["split"]][,"control"])
    # n_merged <- unique(rearrangements[[group]][[chr]][["merged"]][,"test"])
    n_split <- rearrangements[[group]][[chr]][["split"]]
    n_merged <- rearrangements[[group]][[chr]][["merged"]]
    border_alteration <- rearrangements[[group]][[chr]][["border_alteration"]]
    n_expanded <- rearrangements[[group]][[chr]][["expanded"]]
    n_contracted <- rearrangements[[group]][[chr]][["contracted"]]
    n_shift <- rearrangements[[group]][[chr]][["shift"]]
    
    if ( length(n_depleted) > 0 ){
      for (cod in n_depleted){
        depleted[[group]] <- 
          rbind(depleted[[group]], 
                data.frame(chromosome = chr,
                           start = c_CODs[[ctrl_grp]][[chr]][cod, "START"],
                           end  = c_CODs[[ctrl_grp]][[chr]][cod, "END"],
                           mean_binsignal = find_mean_binsignal(ctrl_grp, chr,
                                                                start = CODs[[ctrl_grp]][[chr]][cod, "START"],
                                                                end = CODs[[ctrl_grp]][[chr]][cod, "END"])))
      }
    }
    
    if ( length(n_emerged) > 0 ){
      for (cod in n_emerged){
        emerged[[group]] <- 
          rbind(emerged[[group]],
                data.frame(chromosome = chr,
                           start = c_CODs[[group]][[chr]][cod, "START"],
                           end  = c_CODs[[group]][[chr]][cod, "END"],
                           mean_binsignal = find_mean_binsignal(group, chr,
                                                                start = CODs[[group]][[chr]][cod, "START"],
                                                                end = CODs[[group]][[chr]][cod, "END"])))
      }
    }
    
    if ( length(n_intact) > 0 ){
      for (cod in n_intact){
        intact[[group]] <- 
          rbind(intact[[group]],
                data.frame(chromosome = chr,
                           start = c_CODs[[ctrl_grp]][[chr]][cod, "START"],
                           end  = c_CODs[[ctrl_grp]][[chr]][cod, "END"],
                           mean_binsignal = find_mean_binsignal(ctrl_grp, chr,
                                                                start = CODs[[ctrl_grp]][[chr]][cod, "START"],
                                                                end = CODs[[ctrl_grp]][[chr]][cod, "END"])))
      }
    }
    
    # if ( length(n_split) > 0 ){
    #   for (cod in n_split){
    #     split[[group]] <- 
    #       rbind(split[[group]],
    #             data.frame(chromsome = chr,
    #                        start = c_CODs[[ctrl_grp]][[chr]][cod, "START"],
    #                        end  = c_CODs[[ctrl_grp]][[chr]][cod, "END"],
    #                        mean_binsignal = find_mean_binsignal(ctrl_grp, chr,
    #                                                             start = CODs[[ctrl_grp]][[chr]][cod, "START"],
    #                                                             end = CODs[[ctrl_grp]][[chr]][cod, "END"])))
    
    if ( nrow(n_split) > 0 ){
      for (i in 1:nrow(n_split)){
        split[[group]] <-
          rbind(split[[group]],
                data.frame(chromosome = chr, 
                           start_H = c_CODs[[ctrl_grp]][[chr]][n_split[i, "control"], "START"], 
                           end_H  = c_CODs[[ctrl_grp]][[chr]][n_split[i, "control"], "END"], 
                           mean_binsignal_H = find_mean_binsignal(ctrl_grp, chr,
                                                                  start = CODs[[ctrl_grp]][[chr]][n_split[i, "control"], "START"],
                                                                  end = CODs[[ctrl_grp]][[chr]][n_split[i, "control"], "END"]),
                           start_D = c_CODs[[group]][[chr]][n_split[i, "test"], "START"], 
                           end_D  = c_CODs[[group]][[chr]][n_split[i, "test"], "END"], 
                           mean_binsignal_D = find_mean_binsignal(group, chr,
                                                                  start = CODs[[group]][[chr]][n_split[i, "test"], "START"],
                                                                  end = CODs[[group]][[chr]][n_split[i, "test"], "END"]), 
                           Jaccard_index = n_split[i, "overlap"]))
      }
    }
    
    # if ( length(n_merged) > 0 ){
    #   for (cod in n_merged){
    #     merged[[group]] <- 
    #       rbind(merged[[group]],
    #             data.frame(chromosome = chr,
    #                        start = c_CODs[[group]][[chr]][cod, "START"],
    #                        end  = c_CODs[[group]][[chr]][cod, "END"],
    #                        mean_binsignal = find_mean_binsignal(group, chr,
    #                                                             start = CODs[[group]][[chr]][cod, "START"],
    #                                                             end = CODs[[group]][[chr]][cod, "END"])))
    
    if ( nrow(n_merged) > 0 ){
      for (i in 1:nrow(n_merged)){
        merged[[group]] <-
          rbind(merged[[group]],
                data.frame(chromosome = chr, 
                           start_H = c_CODs[[ctrl_grp]][[chr]][n_merged[i, "control"], "START"], 
                           end_H  = c_CODs[[ctrl_grp]][[chr]][n_merged[i, "control"], "END"], 
                           mean_binsignal_H = find_mean_binsignal(ctrl_grp, chr,
                                                                  start = CODs[[ctrl_grp]][[chr]][n_merged[i, "control"], "START"],
                                                                  end = CODs[[ctrl_grp]][[chr]][n_merged[i, "control"], "END"]),
                           start_D = c_CODs[[group]][[chr]][n_merged[i, "test"], "START"], 
                           end_D  = c_CODs[[group]][[chr]][n_merged[i, "test"], "END"], 
                           mean_binsignal_D = find_mean_binsignal(group, chr,
                                                                  start = CODs[[group]][[chr]][n_merged[i, "test"], "START"],
                                                                  end = CODs[[group]][[chr]][n_merged[i, "test"], "END"]), 
                           Jaccard_index = n_merged[i, "overlap"]))
      }
    }
    
    if ( nrow(border_alteration) > 0){
      for (i in 1:nrow(border_alteration)){
        border_shift[[group]] <-
          rbind(border_shift[[group]],
                data.frame(chromosome = chr, 
                           start_H = c_CODs[[ctrl_grp]][[chr]][border_alteration[i, "control"], "START"], 
                           end_H  = c_CODs[[ctrl_grp]][[chr]][border_alteration[i, "control"], "END"], 
                           mean_binsignal_H = find_mean_binsignal(ctrl_grp, chr,
                                                                  start = CODs[[ctrl_grp]][[chr]][border_alteration[i, "control"], "START"],
                                                                  end = CODs[[ctrl_grp]][[chr]][border_alteration[i, "control"], "END"]),
                           start_D = c_CODs[[group]][[chr]][border_alteration[i, "test"], "START"], 
                           end_D  = c_CODs[[group]][[chr]][border_alteration[i, "test"], "END"], 
                           mean_binsignal_D = find_mean_binsignal(group, chr,
                                                                  start = CODs[[group]][[chr]][border_alteration[i, "test"], "START"],
                                                                  end = CODs[[group]][[chr]][border_alteration[i, "test"], "END"]), 
                           Jaccard_index = border_alteration[i, "overlap"], 
                           border_altered = border_alteration[i, "border_altered"]))
      }
    }
    
    if ( nrow(n_expanded) > 0 ){
      for (i in 1:nrow(n_expanded)){
        expanded[[group]] <-
          rbind(expanded[[group]],
                data.frame(chromosome = chr, 
                           start_H = c_CODs[[ctrl_grp]][[chr]][n_expanded[i, "control"], "START"], 
                           end_H  = c_CODs[[ctrl_grp]][[chr]][n_expanded[i, "control"], "END"], 
                           mean_binsignal_H = find_mean_binsignal(ctrl_grp, chr,
                                                                  start = CODs[[ctrl_grp]][[chr]][n_expanded[i, "control"], "START"],
                                                                  end = CODs[[ctrl_grp]][[chr]][n_expanded[i, "control"], "END"]),
                           start_D = c_CODs[[group]][[chr]][n_expanded[i, "test"], "START"], 
                           end_D  = c_CODs[[group]][[chr]][n_expanded[i, "test"], "END"], 
                           mean_binsignal_D = find_mean_binsignal(group, chr,
                                                                  start = CODs[[group]][[chr]][n_expanded[i, "test"], "START"],
                                                                  end = CODs[[group]][[chr]][n_expanded[i, "test"], "END"]), 
                           Jaccard_index = n_expanded[i, "overlap"]))
      }
    }
    
    if ( nrow(n_contracted) > 0 ){
      for (i in 1:nrow(n_contracted)){
        contracted[[group]] <-
          rbind(contracted[[group]],
                data.frame(chromosome = chr, 
                           start_H = c_CODs[[ctrl_grp]][[chr]][n_contracted[i, "control"], "START"], 
                           end_H  = c_CODs[[ctrl_grp]][[chr]][n_contracted[i, "control"], "END"], 
                           mean_binsignal_H = find_mean_binsignal(ctrl_grp, chr,
                                                                  start = CODs[[ctrl_grp]][[chr]][n_contracted[i, "control"], "START"],
                                                                  end = CODs[[ctrl_grp]][[chr]][n_contracted[i, "control"], "END"]),
                           start_D = c_CODs[[group]][[chr]][n_contracted[i, "test"], "START"], 
                           end_D  = c_CODs[[group]][[chr]][n_contracted[i, "test"], "END"], 
                           mean_binsignal_D = find_mean_binsignal(group, chr,
                                                                  start = CODs[[group]][[chr]][n_contracted[i, "test"], "START"],
                                                                  end = CODs[[group]][[chr]][n_contracted[i, "test"], "END"]), 
                           Jaccard_index = n_contracted[i, "overlap"]))
      }
    }
    
    if ( nrow(n_shift) > 0 ){
      for (i in 1:nrow(n_shift)){
        shift[[group]] <-
          rbind(shift[[group]],
                data.frame(chromosome = chr, 
                           start_H = c_CODs[[ctrl_grp]][[chr]][n_shift[i, "control"], "START"], 
                           end_H  = c_CODs[[ctrl_grp]][[chr]][n_shift[i, "control"], "END"], 
                           mean_binsignal_H = find_mean_binsignal(ctrl_grp, chr,
                                                                  start = CODs[[ctrl_grp]][[chr]][n_shift[i, "control"], "START"],
                                                                  end = CODs[[ctrl_grp]][[chr]][n_shift[i, "control"], "END"]),
                           start_D = c_CODs[[group]][[chr]][n_shift[i, "test"], "START"], 
                           end_D  = c_CODs[[group]][[chr]][n_shift[i, "test"], "END"], 
                           mean_binsignal_D = find_mean_binsignal(group, chr,
                                                                  start = CODs[[group]][[chr]][n_shift[i, "test"], "START"],
                                                                  end = CODs[[group]][[chr]][n_shift[i, "test"], "END"]), 
                           Jaccard_index = n_shift[i, "overlap"]))
      }
    }
  }
}

for (group in names(depleted)){
  write.table(format(depleted[[group]], trim = T, scientific = F), paste0("depleted_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(format(emerged[[group]], trim = T, scientific = F), paste0("emerged_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(format(intact[[group]], trim = T, scientific = F), paste0("intact_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(format(split[[group]], trim = T, scientific = F), paste0("split_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(format(merged[[group]], trim = T, scientific = F), paste0("merged_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(format(border_shift[[group]], trim = T, scientific = F), paste0("border_shift_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(format(expanded[[group]], trim = T, scientific = F), paste0("expanded_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(format(contracted[[group]], trim = T, scientific = F), paste0("contracted_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(format(shift[[group]], trim = T, scientific = F), paste0("shift_DCEs_", group),
              quote = F, sep = "\t", row.names = F, col.names = T)
}
