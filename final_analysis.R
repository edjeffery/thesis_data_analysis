library(RMySQL)
library(anytime)
library(lubridate)
library(geosphere)
library(ggplot2)
library(geosphere)
library(dplyr)
library(ggmap)
library(robustbase)
library(robfilter)
library(zoo) # rollapply()

########################################################

# Get data from MySQL database
# 
# mydb = dbConnect(MySQL(),
#                  user='ej301',
#                  password='Aesoh8oWPu1xohMu',
#                  dbname='ej301',
#                  host='mysql5host.bath.ac.uk')
# # Route data
# rs = dbSendQuery(mydb, 'SELECT *
#                  FROM route
#                  ORDER BY time;')
# 
# routeData = fetch(rs, n=-1)
# dbClearResult(rs)
# 
# # Locus data
# rs = dbSendQuery(mydb, 'SELECT *
#                  FROM locus
#                  ORDER BY time;')
# 
# locusData = fetch(rs, n=-1)
# dbClearResult(rs)
# 
# # Event data
# rs = dbSendQuery(mydb, 'SELECT *
#                  FROM event
#                  ORDER BY time;')
# 
# eventData = fetch(rs, n=-1)
# dbClearResult(rs)
# 
# dbDisconnect(mydb)

routeData <- read.csv("R/routeData.csv") 
locusData <- read.csv("R/locusData.csv")
eventData <- read.csv("R/eventData.csv")

########################################################

# Add additional columns to each dataframe

AddColumns = function(data) {
  data$day <- weekdays(as.Date(anydate(data$time))) # Day of week column
  data$am <- (hour(anytime(data$time)) < 12) # AM or PM column (AM = TRUE)
  data$date <- anydate(data$time)
  data$period <- ifelse(data$am == TRUE, paste(data$day, "am", sep = "_"), paste(data$day, "pm", sep = "_"))
  data$weekend <- ifelse(data$day == "Saturday" | data$day == "Sunday", "weekend", "week")
  data$week_period <- ifelse(data$am == TRUE, paste(data$weekend, "am", sep = "_"), paste(data$weekend, "pm", sep = "_"))
  data$time_ <- anytime(data$time)
  data$mq131 <- ifelse(data[, c("date")] < anydate("2018-08-02"), data[, c("mq131")], data[, c("mq131")] * 1.1292) # Corrects slight error in calibration
  data$mq131[data$mq131 <= 1.5] <- NA
  data$temp[data$temp <= 0.1] <- NA
  data$pm_raw[data$pm_raw == 0] <- NA
  data$humidity[data$humidity <= 0.1 | data$humidity >= 99.9] <- NA
  data$o3 <- ((((((5 - data$mq131) / data$mq131) * 162.26) / 141.6) / 18.473) ^ (-1 / 0.917) ) * 1.9957 # 1.9957 is the ppb to ug/m^3 convertion factor (https://uk-air.defra.gov.uk/assets/documents/reports/cat06/0502160851_Conversion_Factors_Between_ppb_and.pdf)
  #data$correction_factor <- ifelse(is.na(data$humidity) | data$humidity < 60, 1, 1 + 0.25 * (data$humidity/100) * (data$humidity/100) / (1 - (data$humidity/100)))
  data$correction_factor <- ifelse(is.na(data$humidity), 1,
                              ifelse(data$humidity < 40, (data$humidity / 100) * 13,
                                   ifelse(data$humidity >= 40 & data$humidity < 50, (data$humidity / 100) * 8,
                                          ifelse(data$humidity >= 50 & data$humidity < 60, (data$humidity / 100) * 6,
                                                 ifelse(data$humidity >= 60 & data$humidity < 70, (data$humidity / 100) * 4,
                                                        ifelse(data$humidity >= 70 & data$humidity < 80, (data$humidity / 100) * 1.75,
                                                               ifelse(data$humidity >= 80 & data$humidity < 90, (data$humidity / 100) * 1.5, (data$humidity / 100) * 1)))))))
  data$pm <- data$pm_raw * data$correction_factor
  return(data)
}

dfList <- list(routeData, locusData, eventData)
dfList <- lapply(dfList, AddColumns)

routeData <- dfList[[1]]
locusData <- dfList[[2]]
eventData <- dfList[[3]]

########################################################

#Histograms

hist(routeData$pm_raw, breaks=50, freq=TRUE, main="",
     xlab=expression(paste("PM (", mu, g/m^3,")")), ylab="Frequency")

hist(routeData$o3, breaks=50, freq=TRUE, main="",
     xlab=expression(paste("Ozone (", mu, g/m^3,")")), ylab="Frequency")

hist(routeData$mq135, breaks=50, freq=TRUE, main="",
     xlab="Voltage", ylab="Frequency")

hist(routeData$temp, breaks=50, freq=TRUE, main="",
     xlab="Celsius", ylab="Frequency")

hist(routeData$humidity, breaks=50, freq=TRUE, main="",
     xlab="Percentage", ylab="Frequency")

# Boxplots

boxplot(routeData$pm_raw,
        xlab = "",
        ylab="")
title(xlab = "All PM measurements", ylab = expression(paste("PM (", mu, g/m^3,")")), line=2.5)

boxplot(routeData$o3,
        xlab="",
        ylab="")
title(xlab = "All MQ131 measurements", ylab = expression(paste("Ozone (", mu, g/m^3,")")), line=2.5)

boxplot(routeData$mq135,
        xlab="",
        ylab="")
title(xlab = "All MQ135 measurements", ylab ="MQ135 Voltage (V)", line=2.5)

boxplot(routeData$temp,
        xlab="",
        ylab="")
title(xlab = "All DHT22 temperature measurements", ylab ="Celsius", line=2.5)

boxplot(routeData$humidity,
        xlab="",
        ylab="")
title(xlab = "All DHT22 humidity measurements", ylab ="Percentage", line=2.5)

# Longitude and Latitude waveform plots

plot(eventData$time_, eventData$longitude, type = "l",
     xlab = "Time",
     ylab = "Longitude")
abline(h = -2.3908, col = "red")
abline(h = -2.3982, col = "red")

plot(eventData$time_, eventData$latitude, type = "l",
     xlab = "Time",
     ylab = "Latitude")
abline(h = 51.3814, col = "red")
abline(h = 51.3798, col = "red")

########################################################

# Outlier detection and removal

# (i) IQR outlier detection
iqr <- IQR(routeData$pm_raw, na.rm = TRUE)
quantiles <- quantile(routeData$pm_raw, na.rm = TRUE)
q1 <- quantiles[2]
q3 <- quantiles[4]
upper_bound <- q3 + 1.5 * iqr
lower_bound <- q1 - 1.5 * iqr
routeData$outlier_iqr_pm <- ifelse(!is.na(routeData$pm_raw), q3 + 1.5 * iqr < routeData$pm_raw, FALSE)
quantiles[3]
iqr
upper_bound
lower_bound
sum(routeData$outlier_iqr_pm)

iqr <- IQR(routeData$o3, na.rm = TRUE)
quantiles <- quantile(routeData$o3, na.rm = TRUE)
q1 <- quantiles[2]
q3 <- quantiles[4]
upper_bound <- q3 + 1.5 * iqr
lower_bound <- q1 - 1.5 * iqr
routeData$outlier_iqr_o3 <- ifelse(upper_bound < routeData$o3 | lower_bound > routeData$o3, TRUE, FALSE)
quantiles[3]
iqr
upper_bound
lower_bound
sum(routeData$outlier_iqr_o3, na.rm = TRUE)

iqr <- IQR(routeData$mq135, na.rm = TRUE)
quantiles <- quantile(routeData$mq135, na.rm = TRUE)
q1 <- quantiles[2]
q3 <- quantiles[4]
upper_bound <- q3 + 1.5 * iqr
lower_bound <- q1 - 1.5 * iqr
routeData$outlier_iqr_mq135 <- ifelse(q3 + 1.5 * iqr < routeData$mq135, TRUE, FALSE)
sum(routeData$outlier_iqr_mq135, na.rm = TRUE)

# (ii) MAD outlier detection
mad <- mad(routeData$pm_raw, na.rm = TRUE)
median <- median(routeData$pm_raw, na.rm = TRUE)
upper_bound <- median + 3 * mad
lower_bound <- median - 3 * mad
routeData$outlier_mad_pm <- ifelse(!is.na(routeData$pm_raw), upper_bound < routeData$pm_raw, FALSE)
median
mad
upper_bound
lower_bound
sum(routeData$outlier_mad_pm)

mad <- mad(routeData$o3, na.rm = TRUE)
median <- median(routeData$o3, na.rm = TRUE)
upper_bound <- median + 3 * mad
lower_bound <- median - 3 * mad
routeData$outlier_mad_o3 <- ifelse(upper_bound < routeData$o3, TRUE, FALSE)
sum(routeData$outlier_mad_o3, na.rm = TRUE)

mad <- mad(routeData$mq135, na.rm = TRUE)
median <- median(routeData$mq135, na.rm = TRUE)
upper_bound <- median + 3 * mad
lower_bound <- median - 3 * mad
routeData$outlier_mad_mq135 <- ifelse(upper_bound < routeData$mq135, TRUE, FALSE)
sum(routeData$outlier_mad_mq135, na.rm = TRUE)

# (iii) Sn
sn <- Sn(na.omit(routeData$pm_raw))
median <- median(routeData$pm_raw, na.rm = TRUE)
upper_bound <- median + 3 * sn
lower_bound <- median - 3 * sn
routeData$outlier_sn_pm <- ifelse(!is.na(routeData$pm_raw), upper_bound < routeData$pm_raw, FALSE)
median
sn
upper_bound
lower_bound
sum(routeData$outlier_sn_pm)

sn <- Sn(na.omit(routeData$o3))
median <- median(routeData$o3, na.rm = TRUE)
upper_bound <- median + 3 * sn
lower_bound <- median - 3 * sn
routeData$outlier_sn_o3 <- ifelse(upper_bound < routeData$o3, TRUE, FALSE)
sum(routeData$outlier_sn_o3, na.rm = TRUE)

sn <- Sn(na.omit(routeData$mq135))
median <- median(routeData$mq135, na.rm = TRUE)
upper_bound <- median + 3 * sn
lower_bound <- median - 3 * sn
routeData$outlier_sn_mq135 <- ifelse(upper_bound < routeData$mq135, TRUE, FALSE)
sum(routeData$outlier_sn_mq135, na.rm = TRUE)


# Time period outlier detection

routeData <- as.data.frame(routeData %>% group_by(week_period) %>%
   mutate(outlier_g_iqr_pm    = ifelse(is.na(pm_raw), FALSE, quantile(pm_raw, na.rm = TRUE)[4] + 1.5 * IQR(pm_raw, na.rm = TRUE) < pm_raw) | ifelse(is.na(pm_raw), FALSE, quantile(pm_raw, na.rm = TRUE)[2] - 1.5 * IQR(pm_raw, na.rm = TRUE) > pm_raw),
          outlier_g_iqr_o3    = ifelse(is.na(o3), FALSE, quantile(o3, na.rm = TRUE)[4] + 1.5 * IQR(o3, na.rm = TRUE) < o3 )            | ifelse(is.na(o3), FALSE, quantile(o3, na.rm = TRUE)[2] - 1.5 * IQR(o3, na.rm = TRUE) > o3),
          outlier_g_iqr_mq135    = quantile(mq135, na.rm = TRUE)[4] + 1.5 * IQR(mq135, na.rm = TRUE) < mq135 | quantile(mq135, na.rm = TRUE)[2] - 1.5 * IQR(mq135, na.rm = TRUE) > mq135,
          outlier_g_mad_pm    = ifelse(is.na(pm_raw), FALSE, quantile(pm_raw, na.rm = TRUE)[3] + 3 * mad(pm_raw, na.rm = TRUE) < pm_raw )       | ifelse(is.na(pm_raw), FALSE, quantile(pm_raw, na.rm = TRUE)[3] - 3 * mad(pm_raw, na.rm = TRUE) > pm_raw),
          outlier_g_mad_o3    = ifelse(is.na(o3), FALSE, quantile(o3, na.rm = TRUE)[3] + 3 * mad(o3, na.rm = TRUE) < o3 )                   | ifelse(is.na(o3), FALSE, quantile(o3, na.rm = TRUE)[3] - 3 * mad(o3, na.rm = TRUE) > o3),
          outlier_g_mad_mq135    = quantile(mq135, na.rm = TRUE)[3] + 3 * mad(mq135, na.rm = TRUE) < mq135        | quantile(mq135, na.rm = TRUE)[3] - 3 * mad(mq135, na.rm = TRUE) > mq135,
          outlier_g_sn_pm = ifelse(is.na(pm_raw), FALSE, quantile(pm_raw, na.rm = TRUE)[3] + 3 * Sn(na.omit(pm_raw)) < pm_raw )                 | ifelse(is.na(pm_raw), FALSE, quantile(pm_raw, na.rm = TRUE)[3] - 3 * Sn(na.omit(pm_raw)) > pm_raw),
          outlier_g_sn_o3 = ifelse(is.na(o3), FALSE, quantile(o3, na.rm = TRUE)[3] + 3 * Sn(na.omit(o3)) < o3 )                             | ifelse(is.na(o3), FALSE, quantile(o3, na.rm = TRUE)[3] - 3 * Sn(na.omit(o3)) > o3),
          outlier_g_sn_mq135 = quantile(mq135, na.rm = TRUE)[3] + 3 * Sn(na.omit(mq135)) < mq135                  | quantile(mq135, na.rm = TRUE)[3] - 3 * Sn(na.omit(mq135)) > mq135
   ))
sum(routeData$outlier_g_iqr_pm)
sum(routeData$outlier_g_iqr_o3)
sum(routeData$outlier_g_iqr_mq135)
sum(routeData$outlier_g_mad_pm)
sum(routeData$outlier_g_mad_o3)
sum(routeData$outlier_g_mad_mq135)
sum(routeData$outlier_g_sn_pm)
sum(routeData$outlier_g_sn_o3)
sum(routeData$outlier_g_sn_mq135)

# # (i) IQR
# iqrs_pm <- aggregate(routeData[, c("pm_raw")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000), m = median(x), iqr = IQR(x), ub = quantile(x)[4] + 1.5 * IQR(x), lb = quantile(x)[2] - 1.5 * IQR(x), outliers = sum(quantile(x)[4] + 1.5 * IQR(x) < x) ) )
# iqrs_o3 <- aggregate(routeData[, c("o3")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000, na.rm=TRUE), m = median(x, na.rm=TRUE), iqr = IQR(x, na.rm=TRUE), ub = quantile(x, na.rm=TRUE)[4] + 1.5 * IQR(x, na.rm=TRUE), lb = quantile(x, na.rm=TRUE)[2] - 1.5 * IQR(x, na.rm=TRUE), outliers = sum(quantile(x, na.rm=TRUE)[4] + 1.5 * IQR(x, na.rm=TRUE) < x, na.rm=TRUE) ) )
# iqrs_mq135 <- aggregate(routeData[, c("mq135")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000, na.rm=TRUE), m = median(x, na.rm=TRUE), iqr = IQR(x, na.rm=TRUE), ub = quantile(x, na.rm=TRUE)[4] + 1.5 * IQR(x, na.rm=TRUE), lb = quantile(x, na.rm=TRUE)[2] - 1.5 * IQR(x, na.rm=TRUE), outliers = sum(quantile(x, na.rm=TRUE)[4] + 1.5 * IQR(x, na.rm=TRUE) < x, na.rm=TRUE) ) )
# sum(iqrs_pm$x[,"outliers"])
# sum(iqrs_o3$x[,"outliers"])
# sum(iqrs_mq135$x[,"outliers"])
# 
# # (ii) MAD
# mads_pm <- aggregate(routeData[, c("pm_raw")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000), m = median(x), mad = mad(x), ub = quantile(x)[3] + 3 * mad(x), lb = quantile(x)[3] - 3 * mad(x), outliers = sum(quantile(x)[3] + 3 * mad(x) < x, TRUE, FALSE) ) )
# mads_o3 <- aggregate(routeData[, c("o3")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000, na.rm=TRUE), m = median(x, na.rm=TRUE), mad = mad(x, na.rm=TRUE), ub = quantile(x, na.rm=TRUE)[3] + 3 * mad(x, na.rm=TRUE), lb = quantile(x, na.rm=TRUE)[3] - 3 * mad(x, na.rm=TRUE), outliers = sum(quantile(x, na.rm=TRUE)[3] + 3 * mad(x, na.rm=TRUE) < x, na.rm=TRUE) ) )
# mads_mq135 <- aggregate(routeData[, c("mq135")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000, na.rm=TRUE), m = median(x, na.rm=TRUE), mad = mad(x, na.rm=TRUE), ub = quantile(x, na.rm=TRUE)[3] + 3 * mad(x, na.rm=TRUE), lb = quantile(x, na.rm=TRUE)[3] - 3 * mad(x, na.rm=TRUE), outliers = sum(quantile(x, na.rm=TRUE)[3] + 3 * mad(x, na.rm=TRUE) < x, na.rm=TRUE) ) )
# sum(mads_pm$x[,"outliers"])
# sum(mads_o3$x[,"outliers"])
# sum(mads_mq135$x[,"outliers"])
# 
# # (iii) Sn
# sns_pm <- aggregate(routeData[, c("pm_raw")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000), m = median(x), sn = Sn(x), ub = quantile(x)[3] + 3 * Sn(x), lb = quantile(x)[3] - 3 * Sn(x), outliers = sum(quantile(x)[3] + 3 * Sn(x) < x, TRUE, FALSE) ) )
# sns_o3 <- aggregate(routeData[, c("o3")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000, na.rm=TRUE), m = median(x, na.rm=TRUE), sn = Sn(na.omit(x)), ub = quantile(x, na.rm=TRUE)[3] + 3 * Sn(na.omit(x)), lb = quantile(x, na.rm=TRUE)[3] - 3 * Sn(na.omit(x)), outliers = sum(quantile(x, na.rm=TRUE)[3] + 3 * Sn(na.omit(x)) < x, na.rm=TRUE) ) )
# sns_mq135 <- aggregate(routeData[, c("mq135")], list(routeData$week_period), FUN = function(x) c(n = sum(x > -1000, na.rm=TRUE), m = median(x, na.rm=TRUE), sn = Sn(na.omit(x)), ub = quantile(x, na.rm=TRUE)[3] + 3 * Sn(na.omit(x)), lb = quantile(x, na.rm=TRUE)[3] - 3 * Sn(na.omit(x)), outliers = sum(quantile(x, na.rm=TRUE)[3] + 3 * Sn(na.omit(x)) < x, na.rm=TRUE) ) )
# sum(sns_pm$x[,"outliers"])
# sum(sns_o3$x[,"outliers"])
# sum(sns_mq135$x[,"outliers"])

# Moving average

window <- 29

routeData <- as.data.frame(routeData %>% group_by(period) %>%
  mutate(median = rollapplyr(pm_raw, window, median, fill = median(pm_raw, na.rm = TRUE), partial = TRUE),
         mad = rollapplyr(pm_raw, window, mad, fill = mad(pm_raw, na.rm = TRUE), partial = TRUE),
         iqr = rollapplyr(pm_raw, window, IQR, fill = IQR(pm_raw, na.rm = TRUE), na.rm = TRUE),
         sn = rollapplyr(pm_raw, window, function(x) Sn(na.omit(x)), fill = Sn(na.omit(pm_raw)), partial = TRUE),
         outlier_roll_iqr = median + 1.5 * iqr < pm_raw,
         outlier_roll_mad = median + 3 * mad < pm_raw,
         outlier_roll_sn = median + 3 * sn < pm_raw
        ))
sum(routeData$outlier_roll_iqr, na.rm = TRUE)
sum(routeData$outlier_roll_mad, na.rm = TRUE)
sum(routeData$outlier_roll_sn, na.rm = TRUE)

routeData <- as.data.frame(routeData %>% group_by(period) %>%
  mutate(median = rollapplyr(o3, window, median, fill = median(o3, na.rm = TRUE), partial = TRUE),
        mad = rollapplyr(o3, window, mad, fill = mad(o3, na.rm = TRUE), partial = TRUE),
        iqr = rollapplyr(o3, window, IQR, fill = IQR(o3, na.rm = TRUE), na.rm = TRUE),
        sn = rollapplyr(o3, window, function(x) Sn(na.omit(x)), fill = Sn(na.omit(o3)), partial = TRUE),
        outlier_roll_iqr = median + 1.5 * iqr < o3,
        outlier_roll_mad = median + 3 * mad < o3,
        outlier_roll_sn = median + 3 * sn < o3
      ))
sum(routeData$outlier_roll_iqr, na.rm = TRUE)
sum(routeData$outlier_roll_mad, na.rm = TRUE)
sum(routeData$outlier_roll_sn, na.rm = TRUE)

routeData <- as.data.frame(routeData %>% group_by(period) %>%
   mutate(median = rollapplyr(mq135, window, median, fill = median(mq135, na.rm = TRUE), partial = TRUE),
          mad = rollapplyr(mq135, window, mad, fill = mad(mq135, na.rm = TRUE), partial = TRUE),
          iqr = rollapplyr(mq135, window, IQR, fill = IQR(mq135, na.rm = TRUE), na.rm = TRUE),
          sn = rollapplyr(mq135, window, function(x) Sn(na.omit(x)), fill = Sn(na.omit(mq135)), partial = TRUE),
          outlier_roll_iqr = median + 1.5 * iqr < mq135,
          outlier_roll_mad = median + 3 * mad < mq135,
          outlier_roll_sn = median + 3 * sn < mq135
        ))
sum(routeData$outlier_roll_iqr, na.rm = TRUE)
sum(routeData$outlier_roll_mad, na.rm = TRUE)
sum(routeData$outlier_roll_sn, na.rm = TRUE)

########################################################

# Datasets with outliers removed: if any of outlier columns are TRUE then remove

pmData <- routeData[(routeData$outlier_g_iqr_pm == FALSE &
                      routeData$outlier_g_mad_pm == FALSE &
                      routeData$outlier_g_sn_pm == FALSE) &
                      routeData$pm_raw != 0,
                    c("time", "latitude", "longitude", "pm_raw", "altitude", "horizontal_speed", "temp", "humidity", "day", "am", "date", "period", "weekend", "week_period", "time_")]

o3Data <- routeData[routeData$outlier_g_iqr_o3 == FALSE &
                      routeData$outlier_g_mad_o3 == FALSE &
                      routeData$outlier_g_sn_o3 == FALSE,
                    c("time", "latitude", "longitude", "o3", "altitude", "horizontal_speed", "temp", "humidity", "day", "am", "date", "period", "weekend", "week_period", "time_")]

mq135Data <- routeData[routeData$outlier_g_iqr_mq135 == FALSE &
                      routeData$outlier_g_mad_mq135 == FALSE &
                      routeData$outlier_g_sn_mq135 == FALSE,
                    c("time", "latitude", "longitude", "mq135", "altitude", "horizontal_speed", "temp", "humidity", "day", "am", "date", "period", "weekend", "week_period", "time_")]

########################################################

# Locus: map of measurements
bbox <- make_bbox(longitude, latitude, locusData[,], f = 0.005)
m <- get_map(bbox, source = "stamen", zoom = 15)
ggmap(m) +
  geom_point(aes(x = longitude, y = latitude, color = mq135), data = locusData, size = 2)

# Locus: aggregation by time period
locusData$interval <- cut(locusData$time_, breaks = "15 min")
aggregate(locusData[, c("mq135")], list(locusData$interval), mean)
