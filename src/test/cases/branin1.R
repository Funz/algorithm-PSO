## This file should provide following objects, when loaded:
# f : function
# input.f : list of input dimensions, contains list of properties like lower & upper bounds of each dimensions
# output.f : list of output dimensions
# *.f : list of math properties. To be compared with algorithm results
# [print.f] : method to print/plot the function for information

f <- function(x) {
	x1 <- x[1]*15-5   
	x2 <- 0.75 #x[2]*15     
	#matrix(
        (x2 - 5/(4*pi^2)*(x1^2) + 5/pi*x1 - 6)^2 + 10*(1 - 1/(8*pi))*cos(x1) + 10
    #,ncol=1)
}
input.f = list(
    x1=list(min=0,max=1)#,x2=list(min=0,max=1)
)
output.f = "branin1"
#optimise(f,interval = c(0,1))
argmin.f = 0.5566755
min.f = 2.400338

library(testthat)
if (!isTRUE(test_that("f(armgin.f) == f.min",{expect_equal(f(argmin.f),min.f,tolerance = .0001)}))) quit(status=1)

test = function(algorithm_file) {
    results = run.algorithm(algorithm_file, options=list(maxit='30'),fun=list(input=input.f,output=output.f,fun=f))
    if (!isTRUE(test_that("branin min",{expect_equal(as.numeric(results$min),min.f,tolerance = .01)}))) quit(status=1)
}
