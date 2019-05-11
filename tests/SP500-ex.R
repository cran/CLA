require(CLA)

source(system.file("xtraR", "platform-sessionInfo.R", # <<- ../inst/xtraR/platform-sessionInfo.R
		   package = "CLA", mustWork=TRUE)) # moreSessionInfo(), withAutoprint() ..
mS <- moreSessionInfo(print. = TRUE)

data(muS.sp500)

if(requireNamespace("FRAPO")) {
    data(SP500, package = "FRAPO")
    print(dim(SP500)) # 265 * 476
}

(n <- length(muS.sp500$mu)) # 476
system.time(# ~ 9 sec (64-bit); 13.8 sec (32-b florence); seen 27.44 sec on Winb.32
    CLs5c.0.120 <- CLA(muS.sp500$mu, muS.sp500$covar, lB=0, uB=1/20)
)
CLs5c.0.120 # -> print() method

uncall <- function(x) `$<-`(x, call, NULL)

doExtras <- TRUE  # for experiments, not normally
doExtras <- FALSE

if(doExtras) system.time({
    tols <- 10^-c(1,3,5:9,11,14)
    names(tols) <- paste0("10^", round(log10(tols)))
    CLs5c.ls <- lapply(tols, function(tol)
        CLA(muS.sp500$mu, muS.sp500$covar, lB=0, uB=1/20, tol.lambda = tol))
}) #  78.101 elapsed [nb-mm4] ; 46.108 [lynne 2018-10]

if(doExtras) {
    identical(uncall(CLs5c.ls[["10^-7"]]), uncall(CLs5c.0.120))
    for(i in seq_along(tols)[-1]) {
        cat("--=--=--=--=--\n", (n1 <- names(tols[i-1])), " vs. ", (n2 <- names(tols[i])), ": ")
        ae <- all.equal(uncall(CLs5c.ls[[i-1]]),
                        uncall(CLs5c.ls[[ i ]]))
        if(isTRUE(ae)) cat(" are all.equal()\n")
        else {
            CLA.i. <- uncall(CLs5c.ls[[i-1]]) ; wgt <- CLA.i.$weights_set
            cat("are different [all.equal()]: dim(..[[",n1,"]]$weights_set) =",
                dim(wgt)[1],"x", dim(wgt)[2],"\n")
        }
    }
}
## 2018-10 lynne, 64b Fedora 28
##  10^-1  vs.  10^-3 : are different [all.equal()]: dim(..[[ 10^-1 ]]$weights_set) = 476 x 47
##  10^-3  vs.  10^-5 : are different [all.equal()]: dim(..[[ 10^-3 ]]$weights_set) = 476 x 156
##  10^-5  vs.  10^-6 :  are all.equal()
##   ................     "   "    "

op <- options(width = max(500, getOption("width"))) # then it actually fits

if(require(Matrix)) withAutoprint(local = FALSE, {
    ## visualize how weights change "along turning points"
    spWts <- Matrix(CLs5c.0.120$weights_set, sparse=TRUE)
    image(spWts, xlab = "turning point", ylab = "asset number")
    ##
    dim(spWts.non0 <- spWts[rowSums(spWts) > 0 , ])
    round(1000 * spWts.non0)
    ##
    image(spWts.non0, xlab = "turning point", ylab = "asset number")
    wts.non0 <- as(spWts.non0, "matrix")
}) else {
    warning("'Matrix' package not available -- should not happen!")
    wts.non0 <- CLs5c.0.120$weights_set[rowSums(CLs5c.0.120$weights_set) > 0 , ]
    if(is.null(colnames(wts.non0))) ## empty column names for nice printing:
        colnames(wts.non0) <- rep("", ncol(wts.non0))
    print.table(round(1000 * wts.non0), zero.print=".")
}
options(op)
stopifnot(nrow(wts.non0) == 79)
if(FALSE) # once, manually (into tests/ directory
    saveRDS(wts.non0, "wtsn0.rds")
file.info("wtsn0.rds")$size  # 2702049
wtsn0.ref <- readRDS("wtsn0.rds")

## see on all platforms what we get -- on OpenBLAS, the dim() differs !
all.equal(target = wtsn0.ref, current = wts.non0, tol=0)
                                        # expect TRUE only on 64bit (Lnx, R's BLAS)
                                        # 3.10416e-15 and 1.366427e-15 on other BLAS
differWts <- ncol(wtsn0.ref) != ncol(wts.non0)
if(differWts) {
    cat("Got",ncol(wts.non0), "weights from CLA() -- different than ref with",
        ncol(wtsn0.ref), "\n")
    strict <- FALSE # !
} else {
    strict <- mS$strictR
    stopifnot(all.equal(target = wtsn0.ref, current = wts.non0,
                        tol = 1e-13))
}
non.0.assets <- Filter(function(.) . > 0, apply(wts.non0, 1, function(c) sum(c > 0)))

b64.n0 <-
    c(AAPL = 135L, ADSK = 66L, AET = 147L, AMGN = 3L, ATI = 75L,
      AYE = 56L, AZO = 26L, BAX = 95L, BCR = 35L, BDX = 36L, BIIB = 118L,
      BNI = 86L, BRL = 23L, BTU = 27L, BUD = 7L, CCE = 54L, CELG = 128L,
      CI = 69L, CL = 83L, CLX = 53L, CME = 140L, CNX = 16L, COST = 40L,
      CTL = 5L, CVS = 102L, DF = 36L, DGX = 33L, DVN = 14L, ED = 32L,
      EIX = 127L, ESRX = 48L, FCX = 54L, FE = 61L, GILD = 38L, HAL = 31L,
      HES = 40L, HST = 108L, HUM = 71L, INTU = 48L, JNJ = 34L, K = 61L,
      LH = 80L, LLL = 96L, LMT = 83L, LUK = 72L, MCD = 61L, MDT = 43L,
      MMC = 7L, MON = 53L, MRO = 136L, MTW = 66L, MUR = 97L, NEM = 45L,
      NOC = 74L, NUE = 30L, NVDA = 13L, PBG = 72L, PCP = 102L, PDCO = 71L,
      PEP = 69L, PG = 87L, RAI = 110L, RIG = 121L, RRC = 105L, RTN = 90L,
      SII = 27L, SSP = 14L, SYK = 19L, SYMC = 13L, TEX = 36L, TIE = 84L,
      TSO = 115L, TYC = 59L, UST = 127L, WAG = 17L, WFR = 5L, WMT = 6L,
      X = 43L, XTO = 102L)

## 32-bit Linux (Unfortunately, currently  the results are slighly *platform dependent*)
nn <- c("AZO", "BAX", "CLX", "COST", "DGX",  "DVN", "ESRX", "LMT", "MUR", "PEP",
        "RIG", "SYMC", "TYC", "UST")
b32.n0 <- b64.n0
b32.n0[nn]  <- b64.n0[nn] + 1L
nn <- c("AET", "BCR", "CI", "CL", "ED", "FE", "HAL", "MCD", "SII", "SYK")
b32.n0[nn]  <- b64.n0[nn] - 1L

## 64-bit Linux  -no-long-double:
b64nLD.n0 <- b64.n0
nn <- c("AZO", "BAX", "BNI", "CCE", "DVN", "GILD", "INTU", "JNJ", "LH",
        "MDT", "NOC", "PBG", "SYMC")
b64nLD.n0[nn] <- b64.n0[nn] + 1L
nn <- c("ADSK", "BCR", "BDX", "BUD", "CTL", "FE", "MCD", "SII", "SYK")
b64nLD.n0[nn] <- b64.n0[nn] - 1L
b64nLD.n0[["XTO"]] <- 99L # = b...  - 3L

non.0.TARG <- if(mS$ b64) { if(mS$ b64nLD) b64nLD.n0 else b64.n0
              } else b32.n0

## see on all platforms what we get;  typically no diff on 64bit *and* using R's BLAS/Lapack
if(all(non.0.assets == non.0.TARG)) { ## show differences:
    cat("Asset results == non.0.TARG;  showing differences b32 - b64 :\n")
    print(table(b32.n0 -b64.n0))
    dput(names(b64.n0)[b32.n0 -b64.n0 == +1])
    dput(names(b64.n0)[b32.n0 -b64.n0 == -1])
} else {
    cat("\n'non.0.assets' differing from non.0.TARG:\n")
    cat("+1:\n"); dput(names(b64.n0)[non.0.assets - non.0.TARG == +1])
    cat("-1:\n"); dput(names(b64.n0)[non.0.assets - non.0.TARG == -1])
    ## solaris :
      ## +1:
      ## c("AET", "BCR", "CI", "CL", "ED", "FE", "HAL", "MCD", "SII", "SYK")
      ## -1:
      ## c("AZO", "BAX", "CLX", "COST", "DGX", "DVN", "ESRX", "LMT", "MUR",
      ##   "PEP", "RIG", "SYMC", "TYC", "UST")
    ## ATLAS :
      ## +1:
      ## c("AZO", "BRL", "CCE", "CLX", "INTU", "JNJ", "K", "LMT", "MUR",
      ##   "PEP", "SSP", "TYC", "UST", "XTO")
      ## -1:
      ## c("AET", "BCR", "CI", "ED", "FE", "MCD", "NEM", "SII", "SYK", "WMT")
    ## MKL:
      ## +1:
      ## c("CLX", "INTU", "LH", "LLL", "LMT", "PBG", "SYMC", "TYC")
      ## -1:
      ## c("AMGN", "BCR", "BUD", "CL", "CTL", "ED", "HAL", "NEM", "XTO")
    ## OpenBLAS:
      ## +1:
      ## c("BAX", "COST", "HST", "JNJ", "MDT", "MUR", "NOC", "PDCO", "WAG")
      ## -1:
      ## c("BCR", "CI", "CTL", "ED", "EIX", "HAL", "MCD", "RAI", "SYK")

    if(any(isB <- abs(non.0.assets - non.0.TARG) > 1)) {
        cat("more different, showing differences:\n")
        dput((non.0.assets - non.0.TARG)[isB])
    }
    ## OpenBLAS (only!):
      ## c(XTO = -2L)
}

## They have the same names and only differ by  +/- 1:
stopifnot(exprs = {
    identical(names(b64.n0), names(b32.n0))

    if(mS$ b64) !strict || identical(non.0.assets, non.0.TARG) # identical(*) fails on ATLAS, MKL, OpenBLAS
    else if(.Platform$OS.type != "windows") identical(non.0.assets, b32.n0)
    else ## 32-bit Windows
        TRUE ## for now

    differWts || identical(head(CLs5c.0.120$free_indices, 12),
              list(c(295L, 453L), 453L, c(453L, 472L), c(19L, 453L, 472L),
                   c(19L, 453L), 453L, c(15L, 453L), 15L, c(15L, 320L),
                   c(15L, 105L, 320L), c(105L, 320L), c(105L, 320L, 472L)))
})

## Check some of the 'Env<n>' versions: ---------

##' Transform CLA() result to old style (= Env8 / Env9 results):
claStrip <- function(res) {
    class(res) <- NULL
    res$call <- NULL
    res$weights_set <- unname(res$weights_set)
    names(res)[[match("MS_weights", names(res))]] <- "MS_weight" # "MS_weights" w/ final "s"
    res
}

rCLA <- claStrip(CLs5c.0.120)

##' Drop "first turning point" from old, pre-0.95,  CLA() result:
claDrop1st <- function(res) {
    res$weights_set <- res$weights_set[, -1L , drop=FALSE] # drop 1st column
    if(is.matrix(res$MS_weight))
        res$MS_weight <- res$MS_weight[ -1L, , drop=FALSE] # drop 1st row
    for(nm in c("free_indices", "gammas", "lambdas"))
        res[[nm]] <- res[[nm]][-1L]
    res
}

## back compatibility to "old" Env8 results:
nsCLA <- asNamespace("CLA")
if(is.environment(e8 <- nsCLA$Env8)) local(withAutoprint({
    system.time(r8 <- e8$cla.solve(muS.sp500$mu, muS.sp500$covar,
                                   lB = rep(0,n), uB= rep(1/20, n)))
    ## lynne (2017): 9.6--9.8 sec; 2018: 6.1 sec
    if(ncol(claDrop1st(r8)$weights_set) == ncol(rCLA$weights_set))
        stopifnot(all.equal(claDrop1st(r8), rCLA, tol = 1e-14)) # they are the same!
    else cat("#{columns} differ in r8\n")
}))

if(is.environment(e9 <- nsCLA$Env9)) local(withAutoprint({
    system.time(r9 <- e9$cla.solve(muS.sp500$mu, muS.sp500$covar,
                                   lB = rep(0,n), uB= rep(1/20, n)))
    ## lynne(2017): 10.0 sec;  2018: 6.6 sec
    if(ncol(claDrop1st(r9)$weights_set) == ncol(rCLA$weights_set))
        stopifnot(all.equal(claDrop1st(r9), rCLA, tol = 1e-14)) # they are the same!
    else cat("#{columns} differ in r9\n")
}))
