#### Counting the average number of samples in the GEO database
#### Written by: Amin Zollanvari Nov/24/2016
#### Download GEOmetadb.sqlite from https://gbnci-abcc.ncifcrf.gov/geo/
#### Change the working directory to where GEOmetadb.sqlite is saved

source(“http://bioconductor.org/biocLite.R”)
biocLite(“GEOmetadb”)
library(GEOmetadb)

# connect to database getSQLiteFile()
con = dbConnect(SQLite(), “GEOmetadb.sqlite”)
gds.count =dbGetQuery(con, “select gds,sample count from gds”)
gds.count = dbGetQuery(con, “select gds,sample count from gds”)

#To get the columns in con SQLite:
columnDescriptions(sqlite db name=’GEOmetadb.sqlite’)

#This part could be used to filter according to the year
gds.count.year = dbGetQuery(con, ”select gds,sample count, update date from gds”)
Year=strsplit(gds.count.year[,3],split=”-”)
Yearvec=numeric(0)

for (i in 1:length(Year))
    Yearvec=c(Yearvec,Year[[i]][1])

Yearvec=as.numeric(Yearvec)
gds.count.year[,3]=Yearvec
mean(gds.count.year[,2])
