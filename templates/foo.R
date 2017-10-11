foo <- function(args) {
    print(names(args))
    for (key in names(args)) {
        print(args[[key]])
    }
}
df = read.csv("/Users/harry/AGRF/data/PCR1_B5HC6_CAACGACG-CGTGAT_L001/abundance/PCR1_L001_vjassoc.csv", stringsAsFactors=FALSE)
foo(list("PCR1"="PCR1_L001_ACA", "PCR2"="WHAT", df="PCR1"))
