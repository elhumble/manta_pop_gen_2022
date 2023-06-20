#Script: newick.split
#License: GPLv3 or later
#Written by: Marco Milanesi
#Contact: marco.milanesi.mm@gmail.com
#Description: Split newick branches and give the bootstrap value
#Modification date: 2017-04-11

newick.split = function(
  phylip.file
){
  # Inport file
  tree.boot <- read.table(phylip.file, header = F, stringsAsFactors = F)
  
  # Ouput file
  vertices_boot.df <- NULL
  
  # Subset the tree
  for (pos.start in 2:nchar(tree.boot)){
    
    # a new branch is found
    if (substr(tree.boot, start = pos.start, stop = pos.start) == "("){
      count.open <- 0
      count.close <- 0
      pos.open <- 0
      boot.value <- NULL
      
      # Analyse the branch 
      for (pos in pos.start:nchar(tree.boot)){
        # count the open bracket
        if (substr(tree.boot, start = pos, stop = pos) == "("){
          if (count.open == 0){
            pos.open = pos
          }
          count.open = count.open + 1
        }
        
        # count the close bracket
        if (substr(tree.boot, start = pos, stop = pos) == ")"){
          count.close = count.close + 1
        }
        
        # When the open and close braket number is teh same
        if (count.open == count.close & count.open !=0 & count.close !=0){
          # Grep the population names
          tmp <- gsub(pattern = "\\)", replacement = "", 
                      x = gsub(pattern = "\\(",replacement = "", 
                               x = substr(tree.boot, start = pos.open, stop = pos)))
          tmp2 <- unlist(strsplit(x = tmp, split = ","))
          tmp3 <- NULL
          for (a in 1:length(tmp2)){
            tmp3 <- c(tmp3,unlist(strsplit(x = tmp2[a], split = ":"))[1])
          }
          tmp <- tmp3
          rm(tmp2);rm(tmp3)

          # Find the bootstrap value
          for (pos.boot in (pos+2):nchar(tree.boot)){
            if (pos.boot == nchar(tree.boot)){
              boot.value <- substr(tree.boot, start = pos+2, stop = pos.boot)
              break()
            }
            if (substr(tree.boot, start = pos.boot, stop = pos.boot) == ","){
              boot.value <- substr(tree.boot, start = pos+2, stop = pos.boot-1)
              break()
            }
            if (substr(tree.boot, start = pos.boot, stop = pos.boot) == ")"){
              boot.value <- substr(tree.boot, start = pos+2, stop = pos.boot-1)
              break()
            }
          }
          boot.value <- as.numeric(boot.value)
          
          # Write output 
          if (length(vertices_boot.df) == 0){
            vertices_boot.df <- data.frame(cbind(boot.value, paste(tmp, collapse = " ")), stringsAsFactors = FALSE)
          }else{
            vertices_boot.df <- rbind(vertices_boot.df, 
                                      data.frame(cbind(boot.value, paste(tmp, collapse = " ")), stringsAsFactors = FALSE))
          }
          break()
        }
      }
    }
  }
  
  colnames(vertices_boot.df) <- c("VALUE", "POP")
  vertices_boot.df$VALUE <- as.numeric(vertices_boot.df$VALUE)
  return(vertices_boot.df)
}

