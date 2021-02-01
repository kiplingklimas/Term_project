setwd("D:/Rissa/Data")
install.packages("rFIA")
library(rjags)
require(rgdal)
library(sf)
library(raster)

####None of this first section needs to be done again for now- layers are made####

#Bring in UFI perimeters from Megan
UFIPerimeters <- readOGR(dsn = "C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/UtahFireAtlas_20201020", layer = "UtahFireAtlas")
#Convert to SF
UFIP <- st_as_sf(UFIPerimeters)
#Restrict to fires > 2011 and get rid of Rx fires
UFIP <- UFIP[UFIP$FIRE_YEAR > 2011,]
UFIP <- UFIP[UFIP$FIRE_TYPE == "Wildfire",]
UFIP$Source <- "UFI"
write_sf(UFIP, "C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/UtahFireAtlas_20201020/UFIFires2012andAfter.shp")

#Downloaded MTBS wildfires only in Utah, 2012-2018
#Process files downloaded from MTBS: Unzip
setwd("D:/Rissa/Data/MTBSUtah2012_2018")
ziplist <- dir(getwd(), "*.zip", recursive = T, full.names = T)
for(x in ziplist){
  y <- gsub("\\..*","", x)
  y <- basename(y)
  unzip(x,exdir = paste0("D:/Rissa/Data/MTBSUtah2012_2018/Unzipped/", y))
}

#Continue processing. Open up all perimeter shapefiles, merge them into one giant shapefile
setwd("D:/Rissa/Data/MTBSUtah2012_2018/Unzipped")
MTBSlist <- dir(getwd(), "*burn_bndy.shp", recursive = T)
for(x in MTBSlist){
  y <- read_sf(x)
  unwanted <- "Mapping_ID"
  y <- y[,!(names(y) %in% unwanted)]
  z <- gsub("\\..*","", x)
  z <- dirname(z)
  write_sf(y, paste0("D:/Rissa/Data/MTBSUtah2012_2018/Unzipped/Temp/",z, ".shp"))
}
setwd("D:/Rissa/Data/MTBSUtah2012_2018/Unzipped/Temp")
firelist <- dir(getwd(), "*.shp")
firemerge <- do.call(rbind, lapply(firelist,readOGR))
name = "Utah_MTBS_WildFires_2012_2018"
writeOGR(firemerge, dsn = "D:/Rissa/Data/MTBSUtah2012_2018", 
         name, driver = "ESRI Shapefile", overwrite_layer = TRUE)
firemerge <- read_sf("D:/Rissa/Data/MTBSUtah2012_2018/Utah_MTBS_WildFires_2012_2018.shp")
MTBSP <- st_as_sf(firemerge)
MTBSP$Source <- "MTBS"

#Merge UFI and MTBS Perimeters together- first have to get columns & names identical
#Also have to make them have the same projection
names <- c("Fire_ID", "Fire_Name", "Year", "Source", "geometry")
UFIP2 <- UFIP[,c(1,4,10,44)]
colnames(UFIP2) <- names
MTBSP2 <- MTBSP[,c(5,6,7,13)]
MTBSP3 <- st_transform(MTBSP2, 32612)
colnames(MTBSP3) <- names

UtahPerims <- rbind(UFIP2, MTBSP3)
#Got rid of the MTBS "Sheep Creek" fire, it was duplicated in the UFI dataset
UtahPerims <- UtahPerims[!UtahPerims$Fire_ID %in% "UT4002311127520160801",]
#Got rid of two blank rows
UtahPerims <- UtahPerims[!is.na(UtahPerims$Fire_ID),]
UtahPerims$area <- st_area(UtahPerims)
write_sf(UtahPerims, "C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/UtahPerims.shp", overwrite_Layer = TRUE)

####This is all UFI and MTBS perimeters in one shapefile####
UtahPerims <- read_sf("C:/Users/A02248840/Box/Lab_Group/Kipling/UWRAPAnalysis/UtahPerims.shp")
UtahPerims <- read_sf("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/UtahPerims.shp")

#Put in random points
UTPoints <- st_sample(UtahPerims, size=100000)
st_write(UTPoints, dsn="C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis", layer="SamplingPoints", 
         driver = "ESRI Shapefile", overwrite_Layer = TRUE)
UTPoints <- read_sf("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/SamplingPoints.shp")
UTPts <- as_Spatial(UTPoints)

#Get MTBS raster clipped and put together

setwd("D:/Rissa/Data/MTBSUtah2012_2018/Unzipped")
folders <- dir(getwd(), full.names = T)
folders <- folders[-4]
folders <- folders[c(1:60,62:87,90:98,100:130)]#88,89,99 don't have a dnbr layer for some reason. Got rid of sheep creek (#61) because it's a repeat from UFI.
for(x in folders){
  dnbr <- raster(list.files(path = x, "*_dnbr.tif", recursive = T, full.names = T, include.dirs = T))
  perimeter <- read_sf(list.files(path=x, "*burn_bndy.shp", recursive = T, full.names = T, include.dirs = T))
  cropped <- mask(dnbr,perimeter)
  z <- basename(x)
  writeRaster(cropped, filename = paste0("D:/Rissa/Data/MTBSUtah2012_2018/CroppeddNBR/",z,".tif"), format="GTiff")
}

setwd("D:/Rissa/Data/MTBSUtah2012_2018/CroppeddNBR")
dnbrlist <- dir(getwd(), "*.tif$")
firerasters <- lapply(dnbrlist,raster)
firerasters$tolerance <- 1
firemerge <- do.call(merge, firerasters)
writeRaster(firemerge, filename = "D:/Rissa/Data/MTBSUtah2012_2018/MTBSRasters.tif",format="GTiff", overwrite=TRUE)

firemerge <- raster("D:/Rissa/Data/MTBSUtah2012_2018/MTBSRasters.tif")
firemerge <- raster("C:/Users/A02248840/Box/Lab_Group/Kipling/UWRAPAnalysis/MTBSUtah2012_2018/MTBSRasters.tif")
firemerge <- raster("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/MTBSUtah2012_2018/MTBSRasters.tif")
firemerge[firemerge < -500] <- NA

#Bring in UFI rasters
UFIrasters <- raster("C:/Users/A02248840/Box/Lab_Group/Kipling/UWRAPAnalysis/UtahFireAtlas_20201020/Utah Fire Atlas_dNBR_20201020/UtahFireAtlas_dNBR.tif")
UFIrasters <- raster("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/UtahFireAtlas_20201020/Utah Fire Atlas_dNBR_20201020/UtahFireAtlas_dNBR.tif")
UFIrasters[UFIrasters < -500] <- NA
UFIrast <- projectRaster(UFIrasters, crs = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m
+nadgrids=@null +no_defs")





#Bring in UWRAP layers (put them in a list so not all by hand)

setwd("C:/Users/A02248840/Box/Lab_Group/Kipling/UWRAPAnalysis/UWRAPtiffs")
setwd("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/UWRAPtiffs")
setwd("D:/Rissa/Data/UWRAPtiffs")
uwraplist <- dir(getwd(), "*.tif$")
ras <- lapply(uwraplist, raster)

ext <- lapply(ras, extract, UTPts, df=TRUE)

ext2 <- unlist(ext)

A <- matrix(ext2, nrow=100000, ncol=66)
A <- cbind(A[,2],A[,4],A[,6],A[,8],A[,10],A[,12],A[,14],A[,16],A[,18],A[,20],A[,22],A[,24]
           ,A[,26],A[,28],A[,30],A[,32],A[,34],A[,36],A[,38],A[,40],A[,42],A[,44],A[,46]
           ,A[,48],A[,50],A[,52],A[,54],A[,56],A[,58],A[,60],A[,62],A[,64],A[,66])
newcolnames <- sub("\\_.*", "", uwraplist)
colnames(A) <- newcolnames

B <- extract(UFIrast, UTPts, df=TRUE)
colnames(B)[2] <- "UFIdNBR"
C <- extract(firemerge, UTPts, df=TRUE)
colnames(C)[2] <- "MTBSdNBR"

AllData <- data.frame(cbind(A,B[,2],C[,2]))
colnames(AllData)[34:35] <- c("UFIdNBR", "MTBSdNBR")

AllData$severity <- do.call(pmax,c(AllData[34:35], na.rm=TRUE))
AllData2 <- AllData[!is.na(AllData$severity),]
AllData2 <- AllData2[!is.na(AllData2$asp),]
AllData2 <- AllData2[!is.na(AllData2$fri),]

library("readxl")
cctable <-read_excel("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/Tables/AllData2.xlsx", sheet = "CC")
AllData2$wwacc_ <- cctable$CC_PERCENT[match(AllData2$wwacc,cctable$VALUE)]
chtable <-read_excel("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/Tables/AllData2.xlsx", sheet = "CH")
AllData2$wwach_ <- chtable$CH_FT[match(AllData2$wwach,chtable$VALUE)]
fuel13table <- read_excel("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/Tables/AllData2.xlsx", sheet = "Fuels13")
AllData2$wwafuel13_ <- fuel13table$FBFM13[match(AllData2$wwafuel13,fuel13table$VALUE)]
fuel40table <- read_excel("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/Tables/AllData2.xlsx", sheet = "Fuels40")
AllData2$wwafuel40_ <- fuel40table$FBFM40[match(AllData2$wwafuel40,fuel40table$VALUE)]
cbhtable <- read_excel("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/Tables/AllData2.xlsx", sheet = "CBH")
AllData2$wwacbh_ <- cbhtable$CBH_F_X_10[match(AllData2$wwacbh,cbhtable$VALUE)]
vegtable <- read_excel("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/Tables/AllData2.xlsx", sheet = "VegType")
AllData2$EVT <- vegtable$EVT_NAME[match(AllData2$VegType.tif,vegtable$VALUE)]
AllData2$SystGrp <- vegtable$SYSTMGRPNA[match(AllData2$VegType.tif,vegtable$VALUE)]
AllData2$SAFSrm <- vegtable$SAF_SRM[match(AllData2$VegType.tif,vegtable$VALUE)]
AllData2$NVCSOrd <- vegtable$NVCSORDER[match(AllData2$VegType.tif,vegtable$VALUE)]

write.csv(AllData2, file="C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/Tables/FullTable.csv")

ForestData <- AllData2[AllData2$wwach > 0,]
write.csv(ForestData, file="C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/Tables/ForestTable.csv")

AllData2 <- read.csv("C:/Users/A02248840/Box/Lab_Group/Kipling/UWRAPAnalysis/SeverityData/SeverityCorrelations.csv")
AllData2 <- read.csv("C:/Users/ylari/Box/Lab_Group/Kipling/UWRAPAnalysis/SeverityData/SeverityCorrelations.csv")


# subset and re-bind
tree <- subset(data, NVCSOrd == "Tree-dominated")

trees <- as.data.frame(as.matrix(tree))

c <- subset(trees, SYSTMGRPPH == "Conifer")
o <- subset(trees, SYSTMGRPPH == "Hardwood")
s <- subset(trees, SYSTMGRPPH == "Conifer-Hardwood")
r <- subset(trees, SYSTMGRPPH == "Riparian")

library(corrplot)
L <- cor(AllData2, use = "complete.obs")
corrplot(L, method="number")

UFITable <- AllData[!is.na(AllData$UFIdNBR),]
UFITable <- UFITable[,-13]
write.csv(UFITable, file="C:/Users/A02248840/Box/Lab_Group/Kipling/UWRAPAnalysis/UFISeverityCorrelations.csv")
MTBSTable <- AllData[!is.na(AllData$MTBSdNBR),]
MTBSTable <- MTBSTable[,-12]
write.csv(MTBSTable, file="C:/Users/A02248840/Box/Lab_Group/Kipling/UWRAPAnalysis/MTBSSeverityCorrelations.csv")


M <- cor(UFITable, use = "complete.obs")
corrplot(M, is.corr = FALSE)

N <- cor(MTBSTable, use = "complete.obs")
corrplot(N, is.corr = FALSE)

#Then run correlations on the rasters:
library(spatialEco)
rasterCorrelation(UFIrasters, fti, s = 3, type = "pearson")
