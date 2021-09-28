#help: Particle Swarm Optimization
#tags: Optimization
#authors: Yann Richet <yann.richet@irsn.fr>, Claus Bendtsen <papyrus.bendtsen@gmail.com>
#references: Clerc, M. et al. (2010) http://www.particleswarm.info/standard\_pso\_2007.c
#require: future
#options: maxit='10';seed='123'
#options.help: maxit=maximum number of iterations

PSO <- function(options) {
    algorithm = new.env()
    algorithm$maxit <- as.integer(options$maxit)
    algorithm$seed <- as.integer(options$seed)

    algorithm$id <- floor(1000*runif(1))
    
    return(algorithm)
}

getInitialDesign <- function(algorithm, input, output) {
    set.seed(algorithm$seed)
    algorithm$input <- input
    algorithm$output <- output
    d = length(input)
    library(future)
    plan("multisession",workers=2)
    algorithm$job = future(lazy = FALSE,{
        sink(paste0('PSO_',algorithm$id,'.out'),type='output')
        set.seed(algorithm$seed)
        o = psoptim(par=rep(.5,d),fn=function(x) {ask.Y(id=algorithm$id,matrix(x,ncol=d))}, lower=rep(0,d), upper=rep(1,d),control=list(maxit=algorithm$maxit,vectorize=T))
        print(o)
        sink(type='output')
        ask.Y(id=algorithm$id,NULL) # End of loop
    })
    
    algorithm$i = 0
    
    Sys.sleep(.1)
    
    Xn = ask.X(id=algorithm$id)
    algorithm$s = nrow(Xn)

    names(Xn) <- names(algorithm$input)
    return(from01(Xn,algorithm$input))
}

getNextDesign <- function(algorithm, X, Y) {
    #if (algorithm$maxit < algorithm$i) return(NULL)
    algorithm$i = algorithm$i + 1
    
    y = Y[(nrow(Y)-algorithm$s+1):nrow(Y),]
    if (is.na(y)) return(NULL)
    
    tell.Y(id=algorithm$id,y)
    
    Sys.sleep(.1)
    
    Xn = ask.X(id=algorithm$id)
    if (is.null(Xn)) return(NULL)
    algorithm$s = nrow(Xn)

    names(Xn) <- names(algorithm$input)
    return(from01(Xn,algorithm$input))
}

displayResults <- function(algorithm, X, Y) {
    X = X[-nrow(X),]
    Y = Y[-nrow(Y),]
    
    algorithm$files <- "plot.png"
    png(file = algorithm$files, height = 600, width = 600)
    red=(Y-min(Y))/(max(Y)-min(Y))
    #pairs(X,col=rgb(r=red,g=0,b=1-red),Y=Y,d=nrow(X),panel=panel.vec)
    pairs(X,col=rgb(r=red,g=0,b=1-red))
    #pairs(cbind(X,Y))
    dev.off()
    
    ## return HTML string containing plot image
    return(paste(sep = "",
                 "<HTML name='points'>",
                 "<img src='",algorithm$files,"' width='600' height='600'/>",
                 "<br/></HTML>",
                 "<min>",min(Y),"</min>",
                 "<argmin>",paste0(collapse = ",",X[which.min(Y),]),"</argmin>")
          )
}

displayResultsTmp <- displayResults

#' @ref https://cran.r-project.org/web/packages/pso/ , with slight mods for vectorized objective fun
psoptim <- function (par, fn, gr = NULL, ..., lower=-1, upper=1,
                     control = list()) {
    
    #apply.vectorized <<- TRUE
    
    fn1 <- function(par) fn(par, ...)/p.fnscale
    mrunif <- function(n,m,lower,upper) {
        return(matrix(runif(n*m,0,1),nrow=n,ncol=m)*(upper-lower)+lower)
    }
    norm <- function(x) sqrt(sum(x*x))
    rsphere.unif <- function(n,r) {
        temp <- runif(n)
        return((runif(1,min=0,max=r)/norm(temp))*temp)
    }
    svect <- function(a,b,n,k) {
        temp <- rep(a,n)
        temp[k] <- b
        return(temp)
    }
    mrsphere.unif <- function(n,r) {
        m <- length(r)
        temp <- matrix(runif(n*m),n,m)
        return(temp%*%diag(runif(m,min=0,max=r)/apply(temp,2,norm)))
    }
    npar <- length(par)
    lower <- as.double(rep(lower, ,npar))
    upper <- as.double(rep(upper, ,npar))
    con <- list(trace = 0, fnscale = 1, maxit = 1000L, maxf = Inf,
                abstol = -Inf, reltol = 0, REPORT = 10,
                s = NA, k = 3, p = NA, w = 1/(2*log(2)),
                c.p = .5+log(2), c.g = .5+log(2), d = NA,
                v.max = NA, rand.order = TRUE, max.restart=Inf,
                maxit.stagnate = Inf,
                vectorize=FALSE, hybrid = FALSE, hybrid.control = NULL,
                trace.stats = FALSE, type = "SPSO2007")
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    ## Argument error checks
    if (any(upper==Inf | lower==-Inf))
        stop("fixed bounds must be provided")
    
    p.type <- pmatch(con[["type"]],c("SPSO2007","SPSO2011"))-1
    if (is.na(p.type)) stop("type should be one of \"SPSO2007\", \"SPSO2011\"")
    
    p.trace <- con[["trace"]]>0L # provide output on progress?
    p.fnscale <- con[["fnscale"]] # scale funcion by 1/fnscale
    p.maxit <- con[["maxit"]] # maximal number of iterations
    p.maxf <- con[["maxf"]] # maximal number of function evaluations
    p.abstol <- con[["abstol"]] # absolute tolerance for convergence
    p.reltol <- con[["reltol"]] # relative minimal tolerance for restarting
    p.report <- as.integer(con[["REPORT"]]) # output every REPORT iterations
    p.s <- ifelse(is.na(con[["s"]]),ifelse(p.type==0,floor(10+2*sqrt(npar)),40),
                  con[["s"]]) # swarm size
    p.p <- ifelse(is.na(con[["p"]]),1-(1-1/p.s)^con[["k"]],con[["p"]]) # average % of informants
    p.w0 <- con[["w"]] # exploitation constant
    if (length(p.w0)>1) {
        p.w1 <- p.w0[2]
        p.w0 <- p.w0[1]
    } else {
        p.w1 <- p.w0
    }
    p.c.p <- con[["c.p"]] # local exploration constant
    p.c.g <- con[["c.g"]] # global exploration constant
    p.d <- ifelse(is.na(con[["d"]]),norm(upper-lower),con[["d"]]) # domain diameter
    p.vmax <- con[["v.max"]]*p.d # maximal velocity
    p.randorder <- as.logical(con[["rand.order"]]) # process particles in random order?
    p.maxrestart <- con[["max.restart"]] # maximal number of restarts
    p.maxstagnate <- con[["maxit.stagnate"]] # maximal number of iterations without improvement
    p.vectorize <- as.logical(con[["vectorize"]]) # vectorize?
    if (is.character(con[["hybrid"]])) {
        p.hybrid <- pmatch(con[["hybrid"]],c("off","on","improved"))-1
        if (is.na(p.hybrid)) stop("hybrid should be one of \"off\", \"on\", \"improved\"")
    } else {
        p.hybrid <- as.integer(as.logical(con[["hybrid"]])) # use local BFGS search
    }
    p.hcontrol <- con[["hybrid.control"]] # control parameters for hybrid optim
    if ("fnscale" %in% names(p.hcontrol))
        p.hcontrol["fnscale"] <- p.hcontrol["fnscale"]*p.fnscale
    else
        p.hcontrol["fnscale"] <- p.fnscale
    p.trace.stats <- as.logical(con[["trace.stats"]]) # collect detailed stats?
    
    if (p.trace) {
        message("S=",p.s,", K=",con[["k"]],", p=",signif(p.p,4),", w0=",
                signif(p.w0,4),", w1=",
                signif(p.w1,4),", c.p=",signif(p.c.p,4),
                ", c.g=",signif(p.c.g,4))
        message("v.max=",signif(con[["v.max"]],4),
                ", d=",signif(p.d,4),", vectorize=",p.vectorize,
                ", hybrid=",c("off","on","improved")[p.hybrid+1])
        if (p.trace.stats) {
            stats.trace.it <- c()
            stats.trace.error <- c()
            stats.trace.f <- NULL
            stats.trace.x <- NULL
        }
    }
    ## Initialization
    if (p.reltol!=0) p.reltol <- p.reltol*p.d
    if (p.vectorize) {
        lowerM <- matrix(lower,nrow=npar,ncol=p.s)
        upperM <- matrix(upper,nrow=npar,ncol=p.s)
    }
    X <- mrunif(npar,p.s,lower,upper)
    if (!any(is.na(par)) && all(par>=lower) && all(par<=upper)) X[,1] <- par
    if (p.type==0) {
        V <- (mrunif(npar,p.s,lower,upper)-X)/2
    } else { ## p.type==1
        V <- matrix(runif(npar*p.s,min=as.vector(lower-X),max=as.vector(upper-X)),npar,p.s)
        p.c.p2 <- p.c.p/2 # precompute constants
        p.c.p3 <- p.c.p/3
        p.c.g3 <- p.c.g/3
        p.c.pg3 <- p.c.p3+p.c.g3
    }
    if (!is.na(p.vmax)) { # scale to maximal velocity
        temp <- apply(V,2,norm)
        temp <- pmin.int(temp,p.vmax)/temp
        V <- V%*%diag(temp)
    }
    f.x <- fn1(t(X)) #apply(X,2,fn1)#,mode="vectorized") # first evaluations
    stats.feval <- p.s
    P <- X
    f.p <- f.x
    P.improved <- rep(FALSE,p.s)
    i.best <- which.min(f.p)
    error <- f.p[i.best]
    init.links <- TRUE
    if (p.trace && p.report==1) {
        message("It 1: fitness=",signif(error,4))
        if (p.trace.stats) {
            stats.trace.it <- c(stats.trace.it,1)
            stats.trace.error <- c(stats.trace.error,error)
            stats.trace.f <- c(stats.trace.f,list(f.x))
            stats.trace.x <- c(stats.trace.x,list(X))
        }
    }
    ## Iterations
    stats.iter <- 1
    stats.restart <- 0
    stats.stagnate <- 0
    while (stats.iter<p.maxit && stats.feval<p.maxf && error>p.abstol &&
           stats.restart<p.maxrestart && stats.stagnate<p.maxstagnate) {
        stats.iter <- stats.iter+1
        if (p.p!=1 && init.links) {
            links <- matrix(runif(p.s*p.s,0,1)<=p.p,p.s,p.s)
            diag(links) <- TRUE
        }
        ## The swarm moves
        if (!p.vectorize) {
            if (p.randorder) {
                index <- sample(p.s)
            } else {
                index <- 1:p.s
            }
            for (i in index) {
                if (p.p==1)
                    j <- i.best
                else
                    j <- which(links[,i])[which.min(f.p[links[,i]])] # best informant
                temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
                V[,i] <- temp*V[,i] # exploration tendency
                if (p.type==0) {
                    V[,i] <- V[,i]+runif(npar,0,p.c.p)*(P[,i]-X[,i]) # exploitation
                    if (i!=j) V[,i] <- V[,i]+runif(npar,0,p.c.g)*(P[,j]-X[,i])
                } else { # SPSO 2011
                    if (i!=j)
                        temp <- p.c.p3*P[,i]+p.c.g3*P[,j]-p.c.pg3*X[,i] # Gi-Xi
                    else
                        temp <- p.c.p2*P[,i]-p.c.p2*X[,i] # Gi-Xi for local=best
                    V[,i] <- V[,i]+temp+rsphere.unif(npar,norm(temp))
                }
                if (!is.na(p.vmax)) {
                    temp <- norm(V[,i])
                    if (temp>p.vmax) V[,i] <- (p.vmax/temp)*V[,i]
                }
                X[,i] <- X[,i]+V[,i]
                ## Check bounds
                temp <- X[,i]<lower
                if (any(temp)) {
                    X[temp,i] <- lower[temp]
                    V[temp,i] <- 0
                }
                temp <- X[,i]>upper
                if (any(temp)) {
                    X[temp,i] <- upper[temp]
                    V[temp,i] <- 0
                }
                ## Evaluate function
                if (p.hybrid==1) {
                    temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                                  upper=upper,control=p.hcontrol)
                    V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
                    X[,i] <- temp$par
                    f.x[i] <- temp$value
                    stats.feval <- stats.feval+as.integer(temp$counts[1])
                } else {
                    f.x[i] <- fn1(X[,i])
                    stats.feval <- stats.feval+1
                }
                if (f.x[i]<f.p[i]) { # improvement
                    P[,i] <- X[,i]
                    f.p[i] <- f.x[i]
                    if (f.p[i]<f.p[i.best]) {
                        i.best <- i
                        if (p.hybrid==2) {
                            temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                                          upper=upper,control=p.hcontrol)
                            V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
                            X[,i] <- temp$par
                            P[,i] <- temp$par
                            f.x[i] <- temp$value
                            f.p[i] <- temp$value
                            stats.feval <- stats.feval+as.integer(temp$counts[1])
                        }
                    }
                }
                if (stats.feval>=p.maxf) break
            }
        } else {
            if (p.p==1)
                j <- rep(i.best,p.s)
            else # best informant
                j <- sapply(1:p.s,function(i)
                    which(links[,i])[which.min(f.p[links[,i]])]) 
            temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
            V <- temp*V # exploration tendency
            if (p.type==0) {
                V <- V+mrunif(npar,p.s,0,p.c.p)*(P-X) # exploitation
                temp <- j!=(1:p.s)
                V[,temp] <- V[,temp]+mrunif(npar,sum(temp),0,p.c.p)*(P[,j[temp]]-X[,temp])
            } else { # SPSO 2011
                temp <- j==(1:p.s)
                temp <- P%*%diag(svect(p.c.p3,p.c.p2,p.s,temp))+
                    P[,j]%*%diag(svect(p.c.g3,0,p.s,temp))-
                    X%*%diag(svect(p.c.pg3,p.c.p2,p.s,temp)) # G-X
                V <- V+temp+mrsphere.unif(npar,apply(temp,2,norm))
            }
            if (!is.na(p.vmax)) {
                temp <- apply(V,2,norm)
                temp <- pmin.int(temp,p.vmax)/temp
                V <- V%*%diag(temp)
            }
            X <- X+V
            ## Check bounds
            temp <- X<lowerM
            if (any(temp)) {
                X[temp] <- lowerM[temp] 
                V[temp] <- 0
            }
            temp <- X>upperM
            if (any(temp)) {
                X[temp] <- upperM[temp]
                V[temp] <- 0
            }
            ## Evaluate function
            if (p.hybrid==1) { # not really vectorizing
                for (i in 1:p.s) {
                    temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                                  upper=upper,control=p.hcontrol)
                    V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
                    X[,i] <- temp$par
                    f.x[i] <- temp$value
                    stats.feval <- stats.feval+as.integer(temp$counts[1])
                }
            } else {
                f.x <- fn1(t(X))
                #f.x <- apply(X,2,fn1)#,mode="vectorized")
                stats.feval <- stats.feval+p.s
            }
            temp <- sapply(isTRUE,X=as.numeric(f.x)<f.p)
            if (any(temp)) { # improvement
                P[,temp] <- X[,temp]
                f.p[temp] <- f.x[temp]
                i.best <- which.min(f.p)
                if (temp[i.best] && p.hybrid==2) { # overall improvement
                    temp <- optim(X[,i.best],fn,gr,...,method="L-BFGS-B",lower=lower,
                                  upper=upper,control=p.hcontrol)
                    V[,i.best] <- V[,i.best]+temp$par-X[,i.best] # disregards any v.max imposed
                    X[,i.best] <- temp$par
                    P[,i.best] <- temp$par
                    f.x[i.best] <- temp$value
                    f.p[i.best] <- temp$value
                    stats.feval <- stats.feval+as.integer(temp$counts[1])
                }
            }
            if (stats.feval>=p.maxf) break
        }
        if (p.reltol!=0) {
            d <- X-P[,i.best]
            d <- sqrt(max(colSums(d*d)))
            if (d<p.reltol) {
                X <- mrunif(npar,p.s,lower,upper)
                V <- (mrunif(npar,p.s,lower,upper)-X)/2
                if (!is.na(p.vmax)) {
                    temp <- apply(V,2,norm)
                    temp <- pmin.int(temp,p.vmax)/temp
                    V <- V%*%diag(temp)
                }
                stats.restart <- stats.restart+1
                if (p.trace) message("It ",stats.iter,": restarting")
            }
        }
        init.links <- f.p[i.best]==error # if no overall improvement
        stats.stagnate <- ifelse(init.links,stats.stagnate+1,0)
        error <- f.p[i.best]
        if (p.trace && stats.iter%%p.report==0) {
            if (p.reltol!=0) 
                message("It ",stats.iter,": fitness=",signif(error,4),
                        ", swarm diam.=",signif(d,4))
            else
                message("It ",stats.iter,": fitness=",signif(error,4))
            if (p.trace.stats) {
                stats.trace.it <- c(stats.trace.it,stats.iter)
                stats.trace.error <- c(stats.trace.error,error)
                stats.trace.f <- c(stats.trace.f,list(f.x))
                stats.trace.x <- c(stats.trace.x,list(X))
            }
        }
    }
    if (error<=p.abstol) {
        msg <- "Converged"
        msgcode <- 0
    } else if (stats.feval>=p.maxf) {
        msg <- "Maximal number of function evaluations reached"
        msgcode <- 1
    } else if (stats.iter>=p.maxit) {
        msg <- "Maximal number of iterations reached"
        msgcode <- 2
    } else if (stats.restart>=p.maxrestart) {
        msg <- "Maximal number of restarts reached"
        msgcode <- 3
    } else {
        msg <- "Maximal number of iterations without improvement reached"
        msgcode <- 4
    }
    if (p.trace) message(msg)
    o <- list(par=P[,i.best],value=f.p[i.best],
              counts=c("function"=stats.feval,"iteration"=stats.iter,
                       "restarts"=stats.restart),
              convergence=msgcode,message=msg)
    if (p.trace && p.trace.stats) o <- c(o,list(stats=list(it=stats.trace.it,
                                                           error=stats.trace.error,
                                                           f=stats.trace.f,
                                                           x=stats.trace.x)))
    #apply.vectorized <<- FALSE
    
    return(o)
}

######################################## [min,max] <-> [0,1] ########################################

from01 = function(X, inp) {
  nX = names(X)
  for (i in 1:ncol(X)) {
    namei = nX[i]
    X[,i] = X[,i] * (inp[[ namei ]]$max-inp[[ namei ]]$min) + inp[[ namei ]]$min
  }
  return(X)
}

to01 = function(X, inp) {
  nX = names(X)
  for (i in 1:ncol(X)) {
    namei = nX[i]
    X[,i] = (X[,i] - inp[[ namei ]]$min) / (inp[[ namei ]]$max-inp[[ namei ]]$min)
  }
  return(X)
}

######################################## asktell_fun.R ########################################

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
.trace = FALSE
.sleep.step = 1
.timeout = 3600*24 # 1 day before timeout

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

