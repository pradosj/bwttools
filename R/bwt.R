

strrotate <- function(str) {
  m <- matrix(0L,nchar(str),nchar(str))
  m <- ((col(m) + row(m) - 2L) %% nchar(str)) + 1L
  m <- matrix(substring(str,m,m),nchar(str),nchar(str))
  apply(m,1,paste0,collapse="")
}

bwt <- function(str) {
  x <- strrotate(str)
  x <- sort(x)
  cat(x,"\n",sep="\n")
  x <- substr(x,nchar(x),nchar(x))
  paste0(x,collapse="")
}

ibwt <- function(str) {
  x <- substring(str,1:nchar(str),1:nchar(str))
  X <- x
  for(i in 2:nchar(str)) {
    X <- paste0(x,sort(X))
  }
  cat(paste(X,collapse="\n"),"\n")
  grep("\\$$",X,value=TRUE)
}



x <- bwt("ANANAS$BANANA$")
x
ibwt(x)



bwt("ANA$NAS$ANA$")
