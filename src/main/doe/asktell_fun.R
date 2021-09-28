## These functions are emulating the ask/tell hack. If you put ask.Y function as argument instead of a true function (like sin), it will wait as long as you put the Y values using tell.Y() in the same directory. So a 2nd R session will be used to get the X values asked (using ask.X) and then call tell.Y() which will unlock first session. This finally allows an asynchronized IO between many R/matlab/... sessions.
#' @author Y. Richet, from an idea by D. Sinoquet. Async IO principle was defined by G. Pujol.
#' @test x=matrix(runif(10)); f=sin; print(f(x)); parallel::mcparallel({tell.Y(f(ask.X()))}); Sys.sleep(1); print(ask.Y(x));
#' @test x=matrix(runif(10),ncol=2); f=function(X)rowSums(sin(X)); print(f(x)); parallel::mcparallel({tell.Y(f(ask.X()))}); Sys.sleep(1); print(ask.Y(x));
#' @test x=matrix(runif(10),ncol=2); f=function(X)rowSums(sin(X)); print(f(x));library(future); plan("multisession",workers=2); future(lazy=FALSE,{tell.Y(f(ask.X()))}); Sys.sleep(1); print(ask.Y(x));
#' @example f=sin; future::future(evaluator=plan("multisession"),{while(TRUE) {tell.Y(f(ask.X()))}}); optim(par=.5, fn=ask.Y, lower=0, upper=pi, method="L-BFGS-B")
#' @example f=sin; parallel::mcparallel({while(TRUE) {tell.Y(f(ask.X()))}}); optim(par=.5, fn=ask.Y, lower=0, upper=pi, method="L-BFGS-B")
#' @example optim(par=c(.5,.5), fn=function(x) ask.Y(matrix(x,ncol=2)), lower=c(0,0), upper=c(pi,pi), method="L-BFGS-B")

.default_dev.path="/tmp"
.default_dev.X  =  "X.todo"
.default_dev.dX = "dX.todo"
.default_dev.Y  =  "Y.done"
.default_dev.dY = "dY.done"
#dX = 0.001
.trace = TRUE
.sleep.step=1
.timeout=400000 # about 1 week

ask.Y <- function(x, id=0, dev.X=.default_dev.X, dev.Y=.default_dev.Y, dev.path=NULL, sleep.step=.sleep.step, sleep.init=0, trace=ifelse(exists(".trace"),.trace,TRUE), timeout=.timeout) {
    
    if (trace) cat("?Y(",paste0(collapse = ",",x),") ")
    
    write.io(x,file = X.file(id,dev.X,dev.path))
    
    Sys.sleep(sleep.init)
    t=0
    lock=paste0("ask.Y_",id)
    file.create(lock)
    while (!file.exists(file = Y.file(id,dev.Y,dev.path)) & (timeout>0 & t<timeout)) {
        Sys.sleep(sleep.step)
        t=t+sleep.step
        if (!file.exists(lock)) stop("ask.Y break !")
        if (trace) cat(".")
    }
    file.remove(lock)
    if (timeout>0 & t>=timeout) stop("ask.Y timeout !")
    Sys.sleep(sleep.step)
    if (trace) cat(",")
    
    y = read.io(file = Y.file(id,dev.Y,dev.path))
    
    if (trace) cat("(",paste0(collapse = ",",y),")")

    return(y)
}

ask.dY <- function(x, dX=0.001, id=0, dev.dX=.default_dev.dX, dev.dY=.default_dev.dY, dev.path=NULL, sleep.step=.sleep.step, sleep.init=0, trace=ifelse(exists(".trace"),.trace,TRUE), timeout=.timeout) {
    
    if (trace) cat("?dY(",paste0(collapse = ",",x),") ")

    d = ncol(x)
    # build the finite difference X matrix
    xdx = t(matrix(x,d,d+1))
    for (i in 1:d) {
        dx = dX
        if (xdx[1+i,i] + dx > 1) {dx = -dx}
        xdx[1+i,i] = xdx[1+i,i] + dx
    }
    
    write.io(xdx,file = dX.file(id,dev.dX,dev.path))
    
    Sys.sleep(sleep.init)
    t=0
    lock=paste0("ask.dY_",id)
    file.create(lock)
    while (!file.exists(file = dY.file(id,dev.dY,dev.path)) & (timeout>0 & t<timeout)) {
        Sys.sleep(sleep.step)
        t=t+sleep.step
        if (!file.exists(lock)) stop("ask.dY break !")
        if (trace) cat(".")
    }
    file.remove(lock)
    if (timeout>0 & t>=timeout) stop("ask.dY timeout !")
    Sys.sleep(sleep.step)
    if (trace) cat(",")
    
    ydy = read.io(file = dY.file(id,dev.dY,dev.path))
    
    # extract the finite differences for Y vector
    dy = array(-1,d)
    for (i in 1:d) {
        dy[i] = (ydy[i+1] - ydy[1]) / (xdx[i+1,i] - xdx[1,i])
    }
    
    if (trace) cat("(",paste0(collapse = ",",dy),")")
    
    return(dy)
}


#' @example while(x <- ask.X()) tell.Y(-sin(x))
#' @example while(x <- ask.X()) tell.Y(-sin(x[,1]*cos(x[,2])))
ask.X <- function(id=0, dev.X=.default_dev.X, dev.path=NULL, sleep.step=.sleep.step, sleep.init=0, trace=ifelse(exists(".trace"),.trace,TRUE), timeout=.timeout,clean = T) {
    
    if (trace) cat("?X ")
    
    Sys.sleep(sleep.init)
    t=0
    lock=paste0("ask.X_",id)
    file.create(lock)
    while (!file.exists(file = X.file(id,dev.X,dev.path)) & (timeout>0 & t<timeout)) {
        Sys.sleep(sleep.step)
        t=t+sleep.step
        if (!file.exists(lock)) stop("ask.X break !")
        if (trace) cat(":")
    }
    file.remove(lock)
    if (timeout>0 & t>=timeout) stop("ask.X timeout !")
    Sys.sleep(sleep.step)
    if (trace) cat(";")
    
    x = read.io(file = X.file(id,dev.X,dev.path),clean = clean)
    
    if (trace) cat("(",paste0(collapse = ",",x),")")
    
    return(x)
}

ask.dX <- function(id=0, dev.dX=.default_dev.dX, dev.path=NULL, sleep.step=.sleep.step, sleep.init=0, trace=ifelse(exists(".trace"),.trace,TRUE), timeout=.timeout, clean=T) {
    
    if (trace) cat("?dX ")
    
    Sys.sleep(sleep.init)
    t=0
    lock=paste0("ask.dX_",id)
    file.create(lock)
    while (!file.exists(file = dX.file(id,dev.dX,dev.path)) & (timeout>0 & t<timeout)) {
        Sys.sleep(sleep.step)
        t=t+sleep.step
        if (!file.exists(lock)) stop("ask.dX break !")
        if (trace) cat(":")
    }
    file.remove(lock)
    if (timeout>0 & t>=timeout) stop("ask.dX timeout !")
    Sys.sleep(sleep.step)
    if (trace) cat(";")
    
    dx = read.io(file = dX.file(id,dev.dX,dev.path),clean = clean)
    
    if (trace) cat(paste0(collapse = ",",dx))
    
    return(dx) #as.matrix(dx[,2:ncol(dx)]))
}

tell.Y <- function(y, id=0, dev.Y=.default_dev.Y, dev.path=NULL, trace=ifelse(exists(".trace"),.trace,TRUE)) {
    if (trace) cat("!Y=",y)
    
    write.io(y,file = Y.file(id,dev.Y,dev.path))
}

tell.dY <- function(y, id=0, dev.dY=.default_dev.dY, dev.path=NULL, trace=ifelse(exists(".trace"),.trace,TRUE)) {
    if (trace) cat("!dY=",y)
    
    write.io(y,file = dY.file(id,dev.dY,dev.path))
}

#' @test x=123;write.io(x,"x.dat");read.io("x.dat")
#' @test x=matrix(c(123,456),ncol=2);write.io(x,"x.dat");read.io("x.dat")
#' @test x=matrix(c(123,456,789,101),ncol=2);write.io(x,"x.dat");read.io("x.dat")
write.io <- function(data,file, trace=ifelse(exists(".trace"),.trace,TRUE)) {
    i=0
    while(file.exists(file) & i<100) {Sys.sleep(0.05); i=i+1; if (trace) cat(" ")}
    if (i>=100) stop("file ",file, " already exists !")
    if (trace) cat(">")
    write.table(data,file=file,row.names=FALSE)
}

read.io <- function(file,clean=TRUE, trace=ifelse(exists(".trace"),.trace,TRUE)) {
    t = NULL;
    i=0
    try(t <- as.matrix(read.table(file=file,header=TRUE)),silent=TRUE)
    while(is.null(t) & i<10) {
        Sys.sleep(0.05)
        i=i+1
        if (trace) cat(" ")
        try(t <- as.matrix(read.table(file=file,header=TRUE)),silent=TRUE)
    }
    if (is.null(t) & i<10 & trace) cat("\n:)\n")
    #t <- as.matrix(read.table(file=file,header=TRUE))
    if (clean) {file.remove(file); if (trace) cat("-")}
    if (trace) cat("<")
    return(t)
}

X.file <- function(id=0,dev.X=.default_dev.X,dev.path=.default_dev.path) {
    if (!is.null(dev.path)) {
        dev.X = file.path(dev.path,dev.X)
    }
    X.file = paste(sep="_",dev.X,id)
    #if(isTRUE(.trace)) print(X.file)
    
    return(X.file)
}

dX.file <- function(id=0,dev.dX=.default_dev.dX,dev.path=.default_dev.path) {
    if (!is.null(dev.path)) {
        dev.dX = file.path(dev.path,dev.dX)
    }
    dX.file = paste(sep="_",dev.dX,id)
    
    return(dX.file)
}

Y.file <- function(id=0,dev.Y=.default_dev.Y,dev.path=.default_dev.path) {
    if (!is.null(dev.path)) {
        dev.Y = file.path(dev.path,dev.Y)
    }
    Y.file = paste(sep="_",dev.Y,id)
    #if(isTRUE(.trace)) print(Y.file)
    
    return(Y.file)
}

dY.file <- function(id=0,dev.dY=.default_dev.dY,dev.path=NULL) {
    if (!is.null(dev.path)) {
        dev.dY = file.path(dev.path,dev.dY)
    }
    dY.file = paste(sep="_",dev.dY,id)
    
    return(dY.file)
}
