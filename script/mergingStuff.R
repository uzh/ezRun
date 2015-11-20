
# examples do not work anymore, as there are no "Read Count" or "Read1 [File]" columns in them

# overlapping columns and rows + new cols in each
df1 = data.frame(nam=1:10, a=1:10,b=11:20,c=21:30, row.names = T)
df2 = data.frame(nam=6:13, a=101:108,b=71:78,d=111:118, row.names = T)

merge.dfs(df1, df1)
merge.dfs(df1, df2)
merge.dfs(df2, df1)
merge.dfs(df2, df2)

# overlapping columns and rows + new cols in each
df1 = data.frame(nam=1:10, a=1:10,b=11:20,c=21:30,e=61:70, row.names = T)
df2 = data.frame(nam=6:13, a=101:108,b=71:78,d=111:118, row.names = T)

merge.dfs(df1, df1)
merge.dfs(df1, df2)
merge.dfs(df2, df1)
merge.dfs(df2, df2)


merge.dfs = function(df1, df2, newDsDir=NULL){
  
  # these are used a lot
  rowDf1 = rownames(df1)
  rowDf2 = rownames(df2)
  colDf1 = colnames(df1)
  colDf2 = colnames(df2)
  
  cols = intersect(colDf1, colDf2)
  rowdiff1 = setdiff(rowDf1, rowDf2)
  rowdiff2 = setdiff(rowDf2, rowDf1)
  
  # create new df starting with the rows in df1 and adding those only found in df2 and set the right rownames
  dfNew = rbind(df1[, cols], df2[rowdiff2, cols])
  rownames(dfNew) = c(rowDf1, rowdiff2)
  rowDfNew = rownames(dfNew)
    
  # if colnames are not equal, we're not finished with this part:
  # this code should be able to merge columns correctly while filling in NA's no matter how column names occur in df1 and df2
  if (!setequal(colDf1, colDf2)){
    
    # add df1 columns
    cdiff1 = setdiff(colDf1, colDf2)
    for (i in 1:length(cdiff1)){
      dfNew = cbind(dfNew, c(df1[, cdiff1[i]], rep(NA, length(rowdiff2))))
    }
    
    # add df2 columns while making sure to add everything in the right place
    cdiff2 = setdiff(colDf2, colDf1)
    for (i in 1:length(cdiff2)){
      for (j in 1:nrow(dfNew)){
        colToAdd[j] = ifelse(rowDfNew[j] %in% rowDf2,
                             df2[rowDf2==rowDfNew[j], cdiff2[i]],
                             NA)
      }
      dfNew = cbind(dfNew, colToAdd)
    }
    
    # set the colnames for dfNew
    colnames(dfNew) = c(cols, cdiff1, cdiff2)
  }
  
  # adjust the read count and set the directory
  dfNew$"Read Count" = NA
  cwd = getwd()
  setwdNew(newDsDir)
  
  # loop through rows of dfNew to apply the merging
  for (nm in rowDfNew){
    if (nm %in% rowDf1 && !(nm %in% rowDf2)){ # nm is in df1, but not in df2
      dfNew[nm, "Read Count"] = df1[nm, "Read Count"]
      file = file.path("/srv/gstore/projects", df1[nm, "Read1 [File]"])
      ezSystem(paste("cp", file, "."))
      ds[nm, "Read1 [File]"] = file.path(newDsDir, basename(file))
    } else if (nm %in% rowDf2 && !(nm %in% rowDf1)){ # nm is in df2, but not in df1
      dfNew[nm, "Read Count"] = df2[nm, "Read Count"]
      file = file.path("/srv/gstore/projects", df2[nm, "Read1 [File]"])
      ezSystem(paste("cp", file, "."))
      ds[nm, "Read1 [File]"] = file.path(newDsDir, basename(file))
    } else { # nm is in df1 and df2, thus they need to be merged. there should be no other case.
      dfNew[nm, "Read Count"] = df1[nm, "Read Count"] + df2[nm, "Read Count"]
      file1 = file.path("/srv/gstore/projects", df1[nm, "Read1 [File]"])
      file2 = file.path("/srv/gstore/projects", df2[nm, "Read1 [File]"])
      fileMerged = paste0("combined-", nm, "_R1.fastq.gz")
      cmd = paste("gunzip -c", file1, file2, "|", "pigz -p4 --best >", fileMerged)
      cat(cmd, "\n")
      ezSystem(cmd)
      ds[nm, "Read1 [File]"] = file.path(newDsDir, fileMerged)
    }
  }
  setwd(cwd)
  return(dfNew)
}




