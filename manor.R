require("MANOR")

# Import ----------------------------------------------------------------------
dir.in <- "/home/sabba/Phd/Tonini_IST/Group1_extracted"

for (file.in in list.files(dir.in, pattern='*.txt')){
    spot.names <- c("ProbeID",
                    "Ref_MedianSignal", "Ref_BGMedianSignal",
                    "Sample_MedianSignal", "Sample_BGMedianSignal",
                    "LogRatio", "Chromosome")
    clone.names <- c("ProbeID", "Position", "Chromosome")

    acgh <- import(paste(dir.in, file.in, sep="/"),
                   type="default", sep="\t",
                   spot.names=spot.names,
                   clone.names=clone.names, add.lines=TRUE)

    # Flags configuration ---------------------------------------------------------
    data(flags)
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

    # Normalization ---------------------------------------------------------------
    acgh.norm <- norm(acgh, flag.list=flag.list, FUN=median, na.rm=TRUE)
    acgh.norm <- sort(acgh.norm)

    # Quality Assessment -----------------------------------------------------------
    data(qscores)
    profileCGH <- as.profileCGH(acgh.norm$cloneValues)
    profileCGH <- daglad(profileCGH, llambdabreak=6, lambdaclusterGen=20)
    acgh.norm$cloneValues <- as.data.frame(profileCGH)
    acgh.norm$cloneValues$ZoneGNL <- as.factor(acgh.norm$cloneValues$ZoneGNL)

    qscore.list <- list(smoothness=smoothness.qscore,
                        var.replicate=var.replicate.qscore,
                        dynamics = dynamics.qscore)
    acgh.norm$quality <- qscore.summary.arrayCGH(acgh.norm, qscore.list)

    html.report(acgh.norm, #acgh,
                dir.out=paste(dir.in, "/results/", sep=""),
                x="Position", y=c("LogRatioNorm", "LogRatio"),
                array.name=file.in, chrLim="LimitChr", light=FALSE,
                pch=20, zlim=c(-2, 2), file.name=file.in)

    #report.plot(acgh.norm, chrLim="LimitChr", zlim=c(-1, 1), cex=1)
}
