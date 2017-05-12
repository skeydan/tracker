library(readr)
library(dplyr)
library(ggplot2)
library(GGally)
library(gridExtra)
library(dlm)

form <- function(theta){
  dlmModPoly(order = 1, dV = exp(theta[1]), dW = exp(theta[2]))
}


df <- read_csv('tracker.csv')
df

var <- ts(df$Latitude)
var <- ts(df$Longitude)
var <- Nile
model <- form(dlmMLE(var, parm = c(1, 1), form)$par)
model
filtered <- dlmFilter(var, model)

#model <- form(dlmMLE(Nile, parm = c(1, 1), form)$par)
#filtered <- dlmFilter(Nile, model)
#ggfortify:::autoplot.dlmFiltered(filtered, ts.linetype = 'dashed', fitted.colour = 'blue')
autoplot(filtered, ts.linetype = 'dashed', fitted.colour = 'cyan')

smoothed <- dlmSmooth(filtered)
class(smoothed)
autoplot(smoothed)

p <- autoplot(filtered, fitted.colour = 'cyan')
autoplot(smoothed, ts.colour = 'blue', p = p)



#df <- df %>% select(Stnr, Timestamp, Duration, Longitude, Latitude, LastLongitude, LastLatitude, Distance2D, Speed2D) 
df <- df %>% select(Stnr, Timestamp, Duration, Longitude, Latitude, Distance2D, Distance2DTot, Speed2D, Altitude, LastAltitude) %>%
  mutate(AltDiff = Altitude - LastAltitude, AltDiff2 = AltDiff^2, AltDiff3 = AltDiff^3, AltDiffExp = exp(AltDiff))
#df <- df[70:80, ]
df <- df %>% mutate(Stnr = factor(Stnr))

#distances <- c()
#for(i in 1:nrow(df)) {
#  distances[i] <- distm(c(df$LastLongitude[i], df$LastLatitude[i]), c(df$Longitude[i], df$Latitude[i]))
#}
#df <- df %>% mutate(d2 = distances)
#df

ggplot(df, aes(x = Timestamp, y = Altitude, color = Speed2D)) + geom_point() + facet_wrap(~ Stnr) + 
  scale_colour_gradient(low = "orange", high = "blue")
ggplot(df, aes(x = Distance2DTot, y = Altitude, color = Speed2D)) + geom_point() + facet_wrap(~ Stnr) +
  scale_colour_gradient(low = "orange", high = "blue")

hist(df$Duration[df$Duration>60])
hist(df$Duration[df$Duration<60])
#df <- df %>% filter(Duration < 30)

hist(df$Distance2D)
hist(df$AltDiff)

df69 <- df %>% filter(Stnr == 69)
df395 <- df %>% filter(Stnr == 395)
df428 <- df %>% filter(Stnr == 428)
df529 <- df %>% filter(Stnr == 529)

#df69$madist <- ma(df69$Distance2D, 6, centre = TRUE)

###########################################################################################

dfs <- list(df69, df395, df428, df529)
dfs

# Speed2D ~ Distance2DTot + AltDiff

fits <- Map(function (df) lm(Speed2D ~ Distance2DTot + AltDiff, data = df), dfs)
summaries <- Map(summary, fits)
summaries

for (df in dfs) print(ggpairs(df[c("AltDiff", "Distance2DTot", "Speed2D")]))

plots_speed_by_altdiff <- Map(function (df) ggplot(df, aes(x = AltDiff, y = Speed2D, color = Distance2DTot)) + geom_point() +
                                scale_colour_gradient(low = "orange", high = "blue"), dfs)
do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff, 'ncol' = 2, top = "Speed2D by AltDiff"))

predictions <- Map(function (fit) (predict(fit, interval = 'prediction')), fits)
#predictions <- Map(function (fit) (predict(fit, interval = 'confidence')), fits)

prediction_lists <- Map(function(df, prediction) cbind(AltDiff = df$AltDiff, Distance2DTot = df$Distance2DTot, Speed2D = df$Speed2D, prediction[ ,1:3]),
  dfs, predictions)
prediction_dfs <- Map(as.data.frame, prediction_lists)

prediction_plots <- Map(function(prediction_df) ggplot(prediction_df, aes(x = Distance2DTot)) +
               geom_point(aes(y = Speed2D, color = 'blue')) +
               geom_ribbon(aes(ymin = lwr, ymax = upr), fill = 'cyan', alpha = 0.2) ,
             prediction_dfs)

do.call('grid.arrange', list('grobs' = prediction_plots, 'ncol' = 2, top = "Speed2D by AltDiff and Distance2D: Prediction Intervals"))
#do.call('grid.arrange', list('grobs' = prediction_plots, 'ncol' = 2, top = "Speed2D by AltDiff and Distance2D: Confidence Intervals"))


#####################################################################################

#
plots_speed_by_altdiff <- Map(function (df) ggplot(df, aes(x = AltDiff, y = Speed2D)) + geom_point(color = 'green'), dfs)
do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff, 'ncol' = 2, top = "Speed2D by AltDiff"))

plots_speed_by_altdiff <- Map(function (df) ggplot(df, aes(x = AltDiff2, y = Speed2D)) + geom_point(color = 'blue'), dfs)
do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff2, 'ncol' = 2, top = "Speed2D by AltDiff2"))

plots_speed_by_altdiff <- Map(function (df) ggplot(df, aes(x = AltDiff3, y = Speed2D)) + geom_point(color = 'red'), dfs)
do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff3, 'ncol' = 2, top = "Speed2D by AltDiff3"))

plots_speed_by_altdiff <- Map(function (df) ggplot(df, aes(x = AltDiffExp, y = Speed2D)) + geom_point(color = 'cyan'), dfs)
do.call('grid.arrange', list('grobs' = plots_speed_by_altdiffexp, 'ncol' = 2, top = "Speed2D by AltDiffExp"))



df529 %>% filter(AltDiff > 30)





#####################################################################################

#

route <- read_csv('route.csv') 
route
route <- route %>% filter(RouteNumber == 1) %>% select(Longitude:Distance)
route


ggplot(route, aes(Distance, Altitude)) + geom_line() + 
  geom_line(mapping = aes(x = Distance2DTot, y = Altitude), data = df69, color = 'blue') +
  geom_line(mapping = aes(x = Distance2DTot, y = Altitude), data = df395, color = 'green') +
  geom_line(mapping = aes(x = Distance2DTot, y = Altitude), data = df428, color = 'red') +
  geom_line(mapping = aes(x = Distance2DTot, y = Altitude), data = df529, color = 'violet') 

ggplot(route, aes(Latitude, Longitude)) + geom_line() + 
  geom_line(mapping = aes(Latitude, Longitude), data = df69, color = 'blue') +
  geom_line(mapping = aes(Latitude, Longitude), data = df395, color = 'green') +
  geom_line(mapping = aes(Latitude, Longitude), data = df428, color = 'yellow') +
  geom_line(mapping = aes(Latitude, Longitude), data = df529, color = 'violet') 

ggplot(route, aes(Latitude, Longitude)) + geom_line() + 
  geom_line(mapping = aes(Latitude, Longitude), data = df69, color = 'blue') +
  geom_line(mapping = aes(Latitude, Longitude), data = df395, color = 'green') +
  geom_line(mapping = aes(Latitude, Longitude), data = df428, color = 'yellow') +
  geom_line(mapping = aes(Latitude, Longitude), data = df529, color = 'violet') +
  coord_cartesian(xlim = c(46.320, 46.330), ylim = c(7.200, 7.225))

ggplot(route, aes(Latitude, Longitude)) + geom_line() + 
  geom_line(mapping = aes(Latitude, Longitude), data = df69, color = 'blue') +
  geom_line(mapping = aes(Latitude, Longitude), data = df395, color = 'green') +
  geom_line(mapping = aes(Latitude, Longitude), data = df428, color = 'yellow') +
  geom_line(mapping = aes(Latitude, Longitude), data = df529, color = 'violet') +
  coord_cartesian(xlim = c(46.350, 46.355), ylim = c(7.225, 7.250))


