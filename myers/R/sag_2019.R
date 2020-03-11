library(icesSAG)

## 1  Get stock keys
stocks <- getListStocks(2019)
stocks <- stocks[stocks$Purpose == "Advice",]

## 2  Get summary tables
sumtab <- sapply(stocks$AssessmentKey, getSummaryTable)
sumtab <- lapply(sumtab, type.convert, as.is=TRUE)  # convert char->num
names(sumtab) <- sapply(sumtab, function(x) x$fishstock[1])

## Alternative format
sumtab.df <- do.call(rbind, sumtab)
sumtab.main <- sumtab.df[c("fishstock", "Year", "SSB", "catches", "F")]
rownames(sumtab.main) <- NULL
all.na <- apply(is.na(sumtab.main[c("SSB","catches","F")]), 1, all)
sumtab.main <- sumtab.main[!all.na,]

## 3  Get reference points
refpts <- sapply(stocks$AssessmentKey, getFishStockReferencePoints)
names(refpts) <- sapply(refpts, function(x) x$StockKeyLabel)
refpts.df <- do.call(rbind, refpts)

save(list=ls(), file="sag_2019.RData")
