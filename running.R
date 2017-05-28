library(readr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(dlm)
library(lubridate)

df <- read_csv('tracker.csv')
df

# race starts at 10
df <- df %>% filter(Timestamp > ymd_hms('2016-08-06 10:00:00'))
df

df <- df %>% select(Stnr, Timestamp, Duration, Longitude,  LastLongitude, Latitude, LastLatitude,
                    Altitude, LastAltitude, Distance2D, Distance2DTot, Speed2D) 
df <- df %>% mutate(AltDiff = Altitude - LastAltitude)
df <- df %>% mutate(Stnr = factor(Stnr))

# remove all altdiffs < -10 or > +20
# see boxplots below
df <- df %>% filter(between(AltDiff, -10, 20))

mov <- function(x, n=5) {stats::filter(x,rep(1/n,n), sides=2)}
df <- df %>% group_by(Stnr) %>% mutate(maAltDiff = as.vector(unclass(mov(AltDiff))),
                                       maSpeed = as.vector(unclass(mov(Speed2D))))

# using moving average of altdiff here!
df <- df %>% mutate(AltDiffRatio = maAltDiff/Distance2D)
df <- df %>% mutate(ADR2 = AltDiffRatio^2, ADR3 = AltDiffRatio^3)

###

df69 <- df %>% filter(Stnr == 69)
df395 <- df %>% filter(Stnr == 395)
df428 <- df %>% filter(Stnr == 428)
df529 <- df %>% filter(Stnr == 529)

dfs <- list(df69, df395, df428, df529)






###########################################################################################
#                               Explore / Cleanup                                         #
###########################################################################################

ggplot(df, aes(x = Timestamp, y = Altitude, color = Speed2D)) + geom_point() + facet_wrap(~ Stnr) + 
  scale_colour_gradient(low = "orange", high = "blue")
ggplot(df, aes(x = Distance2DTot, y = Altitude, color = Speed2D)) + geom_point() + facet_wrap(~ Stnr) +
  scale_colour_gradient(low = "orange", high = "blue")

#distances <- c()
#for(i in 1:nrow(df)) {
#  distances[i] <- distm(c(df$LastLongitude[i], df$LastLatitude[i]), c(df$Longitude[i], df$Latitude[i]))
#}
#df <- df %>% mutate(d2 = distances)
#df


altdiff_hists <- Map(function (df) ggplot(df, aes(x = AltDiff)) + geom_histogram() +
                       coord_cartesian(xlim = c(-30,30)), dfs)
do.call('grid.arrange', list('grobs' = altdiff_hists, 'ncol' = 2, top = "Speed2D by AltDiff"))

altdiff_boxplot <- Map(function (df) ggplot(df, aes(x = 1, y = AltDiff)) + geom_boxplot() + 
                         ggtitle(quote(df)), dfs)
do.call('grid.arrange', list('grobs' = altdiff_boxplot, 'ncol' = 2, top = "Speed2D by AltDiff"))


# long time gap: duration doesn't fit! how is duration calculated?
df69strange <- df69 %>% filter(between(Timestamp, ymd_hms('2016-08-06 12:51:00'), ymd_hms('2016-08-06 12:59:00'))) %>%
  select(Stnr:LastAltitude)
df69strange

df529strange <- df529 %>% filter(between(Timestamp, ymd_hms('2016-08-06 12:16:00'), ymd_hms('2016-08-06 12:20:00'))) %>%
  select(- c(AltDiff2, AltDiff3, AltDiffExp))


###########################################################################################
#                               Kalman Filtering for lat/long/alt                         #
###########################################################################################

#var <- ts(df529$Latitude)
#var <- ts(df529$Longitude)
#var <- ts(df529$Altitude)
var <- ts(df529$Altitude)

####### local level (random walk plus noise)  - StructTS ######

#struct <- StructTS(var, type="level")
struct <- StructTS(var, type="trend")
if (struct$code != 0) stop("optimizer did not converge")

print(struct$coef)

cat("Transitional variance:", struct$coef["level"], "\n",
    "Slope variance:", struct$coef["slope"], "\n",
    "Observational variance:", struct$coef["epsilon"], "\n",
    "Initial level of mu:", struct$model0$a[1], "\n",
    "Initial level of lambda:", struct$model0$a[2], "\n",
    "Initial level:", struct$model0$a, "\n")
tsdiag(struct)

###### local level (random walk plus noise)  - DLM ######
buildModPoly1 <- function(v) {
  dV <- exp(v[1])
  dW <- exp(v[2])
  m0 <- v[3]
  dlmModPoly(1, dV=dV, dW=dW, m0=m0)
}

# local linear trend
buildModPoly2 <- function(v) {
  dV <- exp(v[1])
  dW <- exp(v[2:3])
  m0 <- v[4:5]
  dlmModPoly(order=2, dV=dV, dW=dW, m0=m0)
}

varGuess <- var(diff(var), na.rm=TRUE)
mu0Guess <- as.numeric(var[1])
lambda0Guess <- 0.0

#parm <- c(log(varGuess), log(varGuess), log(varGuess))
parm <- c(log(varGuess), log(varGuess), log(varGuess),
          mu0Guess, lambda0Guess)
#mle <- dlmMLE(var, parm=parm, buildModPoly1)
mle <- dlmMLE(var, parm=parm, buildModPoly2)

if (mle$convergence != 0) stop(mle$message)

#model <- buildModPoly1(mle$par)
model <- buildModPoly2(mle$par)
model

filtered <- dlmFilter(var, model)
tsdiag(filtered)
smoothed <- dlmSmooth(filtered)

p <- autoplot(filtered, ts.linetype = 'dashed', fitted.colour = 'blue', ts.colour = 'green')
p
autoplot(smoothed, ts.linetype = 'dotted', ts.colour = 'red', p = p)
plot(cbind(var, filtered$m[-1], smoothed$s[-1]), plot.type='s', col=c("black","red","blue"), ylab="Level", main="", lwd=c(1,2,2))




###########################################################################################
#                  maSpeed ~ AltDiffRatio + ADR2 + ADR3                                   #
###########################################################################################


plots_speed_by_altdiff <- Map(function (df) ggplot(df, aes(x = AltDiffRatio, y = maSpeed)) + geom_point(color = 'green'), dfs)
do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff, 'ncol' = 2, top = "maSpeed by AltDiff"))

plots_speed_by_altdiff2 <- Map(function (df) ggplot(df, aes(x = ADR2, y = maSpeed)) + geom_point(color = 'blue'), dfs)
do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff2, 'ncol' = 2, top = "maSpeed by AltDiff2"))

plots_speed_by_altdiff3 <- Map(function (df) ggplot(df, aes(x = ADR3, y = maSpeed)) + geom_point(color = 'red'), dfs)
do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff3, 'ncol' = 2, top = "maSpeed by AltDiff3"))

fits <- Map(function (df) lm(maSpeed ~ AltDiffRatio + ADR2 + ADR3, data = df), dfs)
#fits <- Map(function (df) lm(maSpeed ~ AltDiffRatio + ADR2 + ADR3 + Distance2DTot, data = df), dfs)

summaries <- Map(summary, fits)
summaries



###########################################################################################
#                                     prediction                                          #
###########################################################################################

pred_df <- df529

fit <- lm(maSpeed ~ AltDiffRatio + ADR2 + ADR3 , data = pred_df)
summary(fit)

r <- read_csv('route.csv') 
r
r <- r %>% filter(RouteNumber == 1) %>% select(AltitudeDifference:Distance)
r <- r %>% mutate(DistanceLag = Distance - lag(Distance))
r
r <- r %>% mutate(AltDiffRatio = AltitudeDifference/DistanceLag)
r

ggplot(pred_df, aes(x = AltDiffRatio)) + geom_histogram()
ggplot(r, aes(x = AltDiffRatio)) + geom_histogram()

r <- r %>% mutate(ADR2 = AltDiffRatio^2, ADR3 = AltDiffRatio^3, Distance2DTot = Distance)
r
r <- r %>% select(AltDiffRatio: Distance2DTot)
r

pred <- predict(fit, newdata = r, interval = 'prediction')

r <- r %>% cbind(pred = pred[,1], lower = pred[,2], upper=pred[,3])
r

ggplot(r, aes(Distance2DTot, pred)) + geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'cyan', alpha = 0.2) +
  coord_cartesian(ylim = c(0,20)) + geom_smooth()


# main problem are single out-of-range predictions that spoil it all
filter(r, pred > 20 | pred < 0)

#####################################################################################

#

route <- read_csv('route.csv') 
route
route <- route %>% filter(RouteNumber == 1) %>% select(Longitude:Distance)
route

ggplot(route, aes(Distance, Altitude)) + geom_line()

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

# zoom in further
ggplot(route, aes(Latitude, Longitude)) + geom_line() + 
  geom_line(mapping = aes(Latitude, Longitude), data = df69, color = 'blue') +
  geom_line(mapping = aes(Latitude, Longitude), data = df395, color = 'green') +
  geom_line(mapping = aes(Latitude, Longitude), data = df428, color = 'yellow') +
  geom_line(mapping = aes(Latitude, Longitude), data = df529, color = 'violet') +
  coord_cartesian(xlim = c(46.3250, 46.3255), ylim = c(7.200, 7.225))




