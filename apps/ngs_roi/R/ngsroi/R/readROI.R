readROI <-function(fpname){
   fi <- function(x,i){
      x[i]
   }
   fni <- function(x,i){
      as.numeric(x[[i]])
   }
   getVec <- function(x){
      mylen=as.numeric(x[5])
      veclen=length(x)
      t1=veclen - mylen + 1
      mylist = unlist(x)
      as.list(as.numeric(mylist[c(t1:veclen)]))
   }
   getVals <- function(y, x, columnNames){
      mylen=as.numeric(x[[1]][5])
      veclen=length(x[[1]])
      t1=veclen - mylen
      for ( colN in c(7:t1)) {
   	y[,columnNames[colN]]=unlist(lapply(x,fni,colN))
      }
      return(y)
   }

   
   rLines = readLines(gzfile(fpname))
   values = rLines[substr(rLines,1,1)!="#"]
   values = strsplit(values,"\t")
   if (length(rLines[substr(rLines,1,1)=="##"])==0){
      columnNames=unlist(list("# 1.  column: reference name",
         "# 2.  column: starting position on reference genome",
         "# 3.  column: end position on reference genome" ,
         "# 4.  column: strand +/-(. if empty)" ,
         "# 5.  column: length of region" ,
         "# 6.  column: annotation" ,
         "# 7.  column: max_Count" ,
         "# 8+. column: individual counts" ))
   }else{
      columnNames = unlist(rLines[substr(rLines,1,2)=="##"])
   }



   df=data.frame(id=unlist(lapply(values,fi,1)),
                 start=unlist(lapply(values,fni,2)),
                 end=unlist(lapply(values,fni,3)),
                 strand=unlist(lapply(values,fi,4)),
                 length=unlist(lapply(values,fni,5)),
                 annot=unlist(lapply(values,fi,6)))
   df=getVals(df, values, columnNames)
   df$vec = lapply(values,getVec)
   df$vec=lapply(df$vec,unlist)
   roiNames=names(df)
   names(df) <- sub("# .*. column: ","",roiNames)
   return (df)
}

