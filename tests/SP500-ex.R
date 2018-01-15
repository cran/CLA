require(CLA)

b64 <- .Machine$sizeof.pointer == 8
cat(sprintf("%d bit platform type '%s'\n", if(b64) 64 else 32, .Platform$OS.type))
(nonWindows <- .Platform$OS.type != "windows")

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

un.call <- function(x) { x$call <- NULL ; x }

doExtras <- TRUE  # for experiments, not normally
doExtras <- FALSE

if(doExtras) system.time({
    tols <- 10^-c(1,3,5:9,11,14)
    names(tols) <- paste0("10^", round(log10(tols)))
    CLs5c.ls <- lapply(tols, function(tol)
        CLA(muS.sp500$mu, muS.sp500$covar, lB=0, uB=1/20, tol.lambda = tol))
}) #  78.101 elapsed [nb-mm4]

if(doExtras) {
    identical(un.call(CLs5c.ls[["10^-7"]]), un.call(CLs5c.0.120))
    for(i in seq_along(tols)[-1]) {
        cat("--=--=--=--=--\n", (n1 <- names(tols[i-1])), " vs. ", (n2 <- names(tols[i])), ": ")
        ae <- all.equal(un.call(CLs5c.ls[[i-1]]),
                        un.call(CLs5c.ls[[ i ]]))
        if(isTRUE(ae)) cat(" are all.equal()\n")
        else {
            CLA.i. <- un.call(CLs5c.ls[[i-1]]) ; wgt <- CLA.i.$weights_set
            cat("are different [all.equal()]: dim(..[[",n1,"]]$weights_set) =",
                dim(wgt)[1],"x", dim(wgt)[2],"\n")
        }
    }
}

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
file.info("wtsn0.rds")$size  # 27049
wtsn0.ref <- readRDS("wtsn0.rds")

## see on all platforms what we get
all.equal(target = wtsn0.ref, current = wts.non0, tol=0) # expect TRUE only on 64bit (Lnx)

stopifnot(all.equal(target = wtsn0.ref, current = wts.non0,
                    tol = 1e-13))

non.0.assets <- Filter(function(.) . > 0, apply(wts.non0, 1, function(c) sum(c > 0)))



b64.n0 <- c(AAPL = 136L, ADSK = 66L, AET = 147L, AMGN = 3L, ATI = 76L,
            AYE = 56L, AZO = 26L, BAX = 95L, BCR = 35L, BDX = 36L, BIIB = 118L,
            BNI = 86L, BRL = 23L, BTU = 28L, BUD = 7L, CCE = 54L, CELG = 129L,
            CI = 69L, CL = 83L, CLX = 53L, CME = 141L, CNX = 17L, COST = 40L,
            CTL = 5L, CVS = 102L, DF = 36L, DGX = 33L, DVN = 14L, ED = 32L,
            EIX = 127L, ESRX = 48L, FCX = 55L, FE = 61L, GILD = 38L, HAL = 31L,
            HES = 41L, HST = 108L, HUM = 71L, INTU = 48L, JNJ = 34L, K = 61L,
            LH = 80L, LLL = 96L, LMT = 83L, LUK = 72L, MCD = 61L, MDT = 43L,
            MMC = 7L, MON = 54L, MRO = 137L, MTW = 67L, MUR = 97L, NEM = 45L,
            NOC = 74L, NUE = 31L, NVDA = 14L, PBG = 72L, PCP = 103L, PDCO = 71L,
            PEP = 69L, PG = 87L, RAI = 110L, RIG = 121L, RRC = 106L, RTN = 90L,
            SII = 27L, SSP = 14L, SYK = 19L, SYMC = 13L, TEX = 37L, TIE = 85L,
            TSO = 116L, TYC = 59L, UST = 127L, WAG = 17L, WFR = 6L, WMT = 6L,
            X = 44L, XTO = 102L)

## 32-bit Linux (Unfortunately, currently  the results are slighly *platform dependent*)
nn <- c("AZO", "BAX", "CLX", "COST", "DGX",  "DVN", "ESRX", "LMT", "MUR", "PEP",
        "RIG", "SYMC", "TYC", "UST")
b32.n0 <- b64.n0
b32.n0[nn]  <- b64.n0[nn] + 1L
nn <- c("AET", "BCR", "CI", "CL", "ED", "FE", "HAL", "MCD", "SII", "SYK")
b32.n0[nn]  <- b64.n0[nn] - 1L

## see on all platforms what we get;  typically no diff on 64bit
if(all(non.0.assets ==  if(b64) b64.n0 else b32.n0)) { ## show differences:
    print(table(b32.n0 -b64.n0))
    dput(names(b64.n0)[b32.n0 -b64.n0 == +1])
    dput(names(b64.n0)[b32.n0 -b64.n0 == -1])
}

## They have the same names and only differ by  +/- 1:
stopifnot(
    identical(names(b64.n0), names(b32.n0))
    ##                              ______      ______
  , if(b64) identical(non.0.assets, b64.n0)
    else if(nonWindows) identical(non.0.assets, b32.n0)
    else ## 32-bit Windows
        TRUE ## for now
  , identical(head(CLs5c.0.120$free_indices, 12),
              list(295L, c(295L, 453L), 453L, c(453L, 472L), c(19L, 453L, 472L),
                   c(19L, 453L), 453L, c(15L, 453L), 15L, c(15L, 320L),
                   c(15L, 105L, 320L), c(105L, 320L)))
)

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

nsCLA <- asNamespace("CLA")
if(is.environment(e8 <- nsCLA$Env8)) local(withAutoprint({
    system.time(r8 <- e8$cla.solve(muS.sp500$mu, muS.sp500$covar,
                                   lB = rep(0,n), uB= rep(1/20, n)))
    ## lynne: 9.6--9.8 sec
    stopifnot(all.equal(r8, rCLA, tol = 1e-14)) # they are the same!
}))

if(is.environment(e9 <- nsCLA$Env9)) local(withAutoprint({
    system.time(r9 <- e9$cla.solve(muS.sp500$mu, muS.sp500$covar,
                                   lB = rep(0,n), uB= rep(1/20, n)))
    ## lynne: 10.0 sec
    stopifnot(all.equal(r9, rCLA, tol = 1e-14)) # they are the same!
}))
