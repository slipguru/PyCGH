require("MANOR")
data(flags)

source('cghutils/manor_normalization.R')

# Import ----------------------------------------------------------------------
dir.in <- "/home/sabba/Phd/Tonini_IST/Group1_extracted"

#for (file.in in list.files(dir.in, pattern='*.txt')){
for (file.in in c("1912_G1_251495014704_1_2.txt")){
    spot.names <- c("ProbeID",
                    "Ref_MedianSignal", "Ref_BGMedianSignal",
                    "Sample_MedianSignal", "Sample_BGMedianSignal",
                    "LogRatio", "Chromosome")
    clone.names <- c("ProbeID", "Position", "Chromosome")

    acgh <- import(paste(dir.in, file.in, sep="/"),
                   type="default", sep="\t",
                   spot.names=spot.names,
                   clone.names=clone.names, add.lines=TRUE)

    acgh.norm <- manor_norm(acgh)

    html.report(acgh.norm, #acgh,
            dir.out=paste(dir.in, "/results_/", sep=""),
            x="Position", y=c("LogRatioNorm", "LogRatio"),
            array.name=file.in, chrLim="LimitChr", light=FALSE,
            pch=20, zlim=c(-2, 2), file.name=file.in)


    #data(cytoband)
    #genome.plot(acgh.norm, col.var="ZoneGNL", chrLim="LimitChr", cex=1)
}
