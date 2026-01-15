library(sf)
library(terra)
library(beepr)
library(stringdist)
library(riem)
library(geosphere)
library(lubridate)
library(dplyr)

find_closest_matches <- function(AA, PO) {
  closest_matches <- matrix(nrow = length(AA), ncol = 2, dimnames = list(NULL, c("match1", "match2")))
  for (i in seq_along(AA)) {
    distances <- stringdist(AA[i], PO)
    closest_indices <- order(distances)[1:2]
    closest_matches[i,] <- PO[closest_indices]
  }
  return(closest_matches)
}
tutti<-read.csv('dets.csv')
posti<-read.csv('m3.csv')



delocs<-unique(tutti$location)
qualocs<-length(delocs)

alocs<-unique(posti$Site)

ref<-matrix(ncol=6,nrow=qualocs,dimnames=list(NULL,c("name_det","ind_det","name_ref","ind_ref","lat","long")))

a<-0
impossibles<-0

for(q in 1:qualocs) 
	{
	if(length(which(delocs[q]==alocs))==1) 
		{
			ref[q,1]<-delocs[q]
			ref[q,3]<-q
			ref[q,2]<-"em3"
			ref[q,4]<-which(alocs==delocs[q])
			impossibles<-c(impossibles,as.numeric(ref[q,4]))
		}
	else
		if(length(which(delocs[q]==alocs))==0) a<-c(a,q)
		else print("problem")
	}	
a<-a[-1]
orphans<-sort(delocs[-a])
hmo<-length(orphans)
impossibles<-impossibles[-1]
possibles<-sort(alocs[-impossibles])


em3<-posti

posti<-read.csv('2021DeploymentsSummary.csv')



alocs<-unique(posti$location)


impossibles<-0
ba<-0
for(q in 1:hmo) 
	{
	ref[a[q],1]<-orphans[q]
	if(length(which(orphans[q]==alocs))==1) 
		{
			ref[q,3]<-a[q]
			ref[a[q],2]<-"2021Deploy"
			ref[a[q],4]<-which(alocs==orphans[q])
			impossibles<-c(impossibles,as.numeric(ref[q,4]))
		}
	else
		if(length(which(orphans[q]==alocs))==0) 
			{
#			print(q)
			beep(1)
			ba<-c(ba,q)
			}
		else print("problem")
	}
		


eneis<-which(is.na(ref[,4]))
ref<-ref[-eneis,]

for(q in which(ref[,2]=="2021Deploy"))
	{
	wawa<-posti[which(posti$location==ref[q,1])[1],]
	ref[q,5]<-wawa[1,3]
	ref[q,6]<-wawa[1,4]
	}
	
for(q in which(ref[,2]=="em3"))
	{
	wawa<-em3[which(em3$Site==ref[q,1])[1],]
	ref[q,5]<-wawa[1,10]
	ref[q,6]<-wawa[1,11]
	}




lat <- as.numeric(ref[, 5])
lon <- as.numeric(ref[, 6])

coords <- data.frame(ref[,1],lon, lat)



points <- vect(coords, crs = "+proj=longlat")
	
choices<-array(list())

choices[[1]]<-c("MYOCIL","MYOEVO","MYOSEP","MYOVOL","MYOLUC")
choices[[2]]<-c("EPTFUS","LASNOC")
choices[[3]]<-"LASBOR"
choices[[4]]<-"LASCIN"

######## this will eventually become
# for(cuspe in choiches) or something 
#### but for now we work with one

cuspe<-choices[[1]]

questi<-tutti[which(tutti$species==cuspe),]



questi$recordingdate <- ymd_hms(questi$recordingdate)
questi$recordingdate <- if_else(hour(questi$recordingdate) < 12, questi$recordingdate - hours(12), questi$recordingdate)
questi$date_formatted <- questi$recordingdate %>% 
    as.Date() %>%
    format("%Y-%m-%d")
result <- questi %>% 
    group_by(location, date_formatted) %>% 
    summarise(num_rows = n()) %>% 
    ungroup()

result<-result[-which(is.na(result$date_formatted)),]

df_result <- as.data.frame(result)



stations <- riem_stations(network = "CA_AB_ASOS")

# Create empty columns for station_id, lat, and lon
coords$ns <- numeric(length(coords$lat))
coords$station_lon <- numeric(length(coords$lat))
coords$station_lat <- numeric(length(coords$lat))
coords$dist <- numeric(length(coords$lat))

for (q in 1:length(coords$lat)) {
  nearest_station <- stations[which.min(sqrt((stations$lat - coords$lat[q])^2 + (stations$lon - coords$lon[q])^2)),]
  
  # Assign station_id, lat, and lon to the corresponding columns
  coords$ns[q] <- nearest_station$id
  coords$station_lat[q] <- nearest_station$lat
  coords$station_lon[q] <- nearest_station$lon
  coords$dist[q] <- distHaversine(c(coords$lon[q], coords$lat[q]), c(coords$station_lon[q], coords$station_lat[q]))

}


this_date <- df_result$date_formatted[q]


#reduce df_result to the stations I have data for
ttM<-df_result
lis<-0
for(q in unique(df_result$location))	if(length(which(q==coords[,1]))<1) lis<-c(lis,q)
lis<-lis[-1]
for(q in 1:length(lis)) df_result<-df_result[-which(df_result$location==lis[q]),]



# add necessary columns to df_result
df_result$station_id <- ""
df_result$lat <- NA
df_result$lon <- NA
df_result$temp <- NA
df_result$dewpointemp <- NA
df_result$humidity <- NA
df_result$windspeed <- NA
df_result$ap <- NA
df_result$visib <- NA
df_result$precipi <- NA
lis<-0

#  function to get one data point per night
filter_midnight <- function(df) {
  return(df[strftime(df$valid, format="%H:%M") == "00:00",])
}


for (q in 1:nrow(df_result)) {
  current_site <- df_result$location[q]
  current_date <- df_result$date_formatted[q]
  nearest_station_id <- coords$ns[which(coords[, 1] == current_site)][1]
  if(is.na(nearest_station_id))
  	{
	beep(7)
	next
	}


  weather_data <- riem_measures(station = nearest_station_id, date_start = current_date, date_end = current_date)
  midnight_data <- filter_midnight(weather_data)
  if(length(midnight_data)<1) 
  	{
  	lis<-c(lis,q)
  	next
  	}
# convert from imperial to metric and add data points I want
  if (nrow(midnight_data) > 0) {
    df_result$station_id[q] <- nearest_station_id
    df_result$lat[q] <- coords$lat[which(coords[, 1] == current_site)][1]
    df_result$lon[q] <- coords$lon[which(coords[, 1] == current_site)][1]
    df_result$temp[q] <- (midnight_data$tmpf[1]-32)*5/9
    df_result$dewpointemp[q] <- (midnight_data$dwpf[1]-32)*5/9
    df_result$humidity[q] <- midnight_data$relh[1]
    df_result$windspeed[q] <- midnight_data$sknt[1]/0.53995680345572
    df_result$ap[q] <- midnight_data$mslp[1]
    df_result$visib[q] <- midnight_data$vsby[1]*0.62137119223733
    df_result$precipi[q] <- midnight_data$p01i[1]*25.4
  }
}


lis<-lis[-1]
df_result<-df_result[-lis,]

for(q in 4:length(df_result[1,])) print(paste(colnames(df_result)[q],"had",length(which(is.na(df_result[,q])==TRUE)),"NAs"))

hmcov<-7

trimming<-0

for(q in 7:length(df_result[2,])) trimming<-c(trimming,which(is.na(df_result[,q])))

trimming<-trimming[-1]
trimming<-sort(unique(trimming))

df_result<-df_result[-trimming,]

#plot points

oords <- data.frame(ref[,1],lon, lat)
oords<-oords[-which(oords$lon<(-125)),]

points <- vect(oords, crs = "+proj=longlat")

plot(points,main="Sampling Sites (black, full) & Weather Stations (red, hollow)")
lines(x=coords$station_lon,y=coords$station_lat,col=2,type='p')




### prepare mcmc list

list_a<-list(det=df_result$num_rows,lat=scale(df_result$lat),lon=scale(df_result$lon),temp=scale(df_result$temp),dewpointemp=scale(df_result$dewpointemp),rh=scale(df_result$humidity),ws=scale(df_result$windspeed),preci=scale(df_result$precipi),hmd=length(df_result$ap),n_cov=hmcov)


save.image('prepd.RData',version=2)
beep(3)
	


