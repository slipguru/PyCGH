require("MANOR")
data(flags)

manor_norm <- function(acgh){
    class(acgh) <- "arrayCGH"

    # Flags configuration -----------------------------------------------------
    local.spatial.flag$args <- list(var="LogRatio", by.var="Chromosome",
                                    nk=5, prop=0.25, thr=0.15, beta=1,
                                    family="gaussian")

    ref.snr.flag$args$var.FG <- "Ref_MedianSignal"
    ref.snr.flag$args$var.BG <- "Ref_BGMedianSignal"
    dapi.snr.flag$args$var.FG <- "Sample_MedianSignal"
    dapi.snr.flag$args$var.BG <- "Sample_BGMedianSignal"
    global.spatial.flag$args$by.var <- "Chromosome"         # Original- ChromosomeArm

    flag.list <- list(spatial=local.spatial.flag,           # S
                      ref.snr = ref.snr.flag,               # B
                      dapi.snr = dapi.snr.flag,             # D
                      global.spatial = global.spatial.flag,
                      replicate=rep.flag                    # E
                      )

    # Normalization -----------------------------------------------------------
    acgh.norm <- norm(acgh, flag.list=flag.list, FUN=median, na.rm=TRUE)
    acgh.norm <- sort(acgh.norm)

    # Quality Assessment ------------------------------------------------------
    data(qscores)
    profileCGH <- as.profileCGH(acgh.norm$cloneValues)
    profileCGH <- daglad(profileCGH, llambdabreak=6, lambdaclusterGen=20)
    acgh.norm$cloneValues <- as.data.frame(profileCGH)
    acgh.norm$cloneValues$ZoneGNL <- as.factor(acgh.norm$cloneValues$ZoneGNL)

    qscore.list <- list(smoothness=smoothness.qscore,
                        var.replicate=var.replicate.qscore,
                        dynamics = dynamics.qscore)
    acgh.norm$quality <- qscore.summary.arrayCGH(acgh.norm, qscore.list)

    return(acgh.norm)
}
