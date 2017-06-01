library(readr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(dlm)
library(lubridate)
library(e1071)
library(ggmap)



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

df <- df %>% mutate(AltDiffRatio = AltDiff/Distance2D)

# important: downhill does not lead to speedup as uphill leads to slowdown!!!
df <- df %>% mutate(AltDiffRatio = ifelse(AltDiffRatio < 0, 0 , AltDiffRatio))

# use moving average 
window_size <- 61
mov <- function(x, n=window_size) {stats::filter(x,rep(1/n,n), sides=2)}
df <- df %>% group_by(Stnr) %>% mutate(maAltDiffRatio = as.vector(unclass(mov(AltDiffRatio))),
                                       maSpeed = as.vector(unclass(mov(Speed2D))))

df <- df %>% mutate(maADR2 = ifelse(maAltDiffRatio < 0, -maAltDiffRatio^2, maAltDiffRatio^2), 
                    maADR3 = maAltDiffRatio^3)

###

df69 <- df %>% filter(Stnr == 69)
df395 <- df %>% filter(Stnr == 395)
df428 <- df %>% filter(Stnr == 428)
df529 <- df %>% filter(Stnr == 529)

dfs <- list(df69, df395, df428, df529)

##

normalize <- function(vec) {
  minval <- min(vec, na.rm = TRUE)
  maxval <- max(vec, na.rm = TRUE)
  (vec-minval) / (maxval-minval)
}
normalize(1:10)
ggplot(df529, aes(Distance2DTot, normalize(maAltDiffRatio))) +
  geom_point() + geom_line(aes(y=normalize(Altitude)))


ggplot(df529, aes(Distance2DTot, maAltDiffRatio)) +
  geom_point() +  
  geom_line(aes(y=maADR2))



###########################################################################################
#                               Explore / Cleanup                                         #
###########################################################################################

# ggplot(df, aes(x = Timestamp, y = Altitude, color = Speed2D)) + geom_point() + facet_wrap(~ Stnr) + 
#   scale_colour_gradient(low = "orange", high = "blue")
# ggplot(df, aes(x = Distance2DTot, y = Altitude, color = Speed2D)) + geom_point() + facet_wrap(~ Stnr) +
#   scale_colour_gradient(low = "orange", high = "blue")
# 
# #distances <- c()
# #for(i in 1:nrow(df)) {
# #  distances[i] <- distm(c(df$LastLongitude[i], df$LastLatitude[i]), c(df$Longitude[i], df$Latitude[i]))
# #}
# #df <- df %>% mutate(d2 = distances)
# #df
# 
# 
# altdiff_hists <- Map(function (df) ggplot(df, aes(x = AltDiff)) + geom_histogram() +
#                        coord_cartesian(xlim = c(-30,30)), dfs)
# do.call('grid.arrange', list('grobs' = altdiff_hists, 'ncol' = 2, top = "Speed2D by AltDiff"))
# 
# altdiff_boxplot <- Map(function (df) ggplot(df, aes(x = 1, y = AltDiff)) + geom_boxplot() + 
#                          ggtitle(quote(df)), dfs)
# do.call('grid.arrange', list('grobs' = altdiff_boxplot, 'ncol' = 2, top = "Speed2D by AltDiff"))
# 
# 
# # long time gap: duration doesn't fit! how is duration calculated?
# df69strange <- df69 %>% filter(between(Timestamp, ymd_hms('2016-08-06 12:51:00'), ymd_hms('2016-08-06 12:59:00'))) %>%
#   select(Stnr:LastAltitude)
# df69strange
# 
# df529strange <- df529 %>% filter(between(Timestamp, ymd_hms('2016-08-06 12:16:00'), ymd_hms('2016-08-06 12:20:00'))) %>%
#   select(- c(AltDiff2, AltDiff3, AltDiffExp))
# 
# 
# ###########################################################################################
# #                               Kalman Filtering for lat/long/alt                         #
# ###########################################################################################
# 
# #var <- ts(df529$Latitude)
# #var <- ts(df529$Longitude)
# #var <- ts(df529$Altitude)
# var <- ts(df529$Altitude)
# 
# ####### local level (random walk plus noise)  - StructTS ######
# 
# #struct <- StructTS(var, type="level")
# struct <- StructTS(var, type="trend")
# if (struct$code != 0) stop("optimizer did not converge")
# 
# print(struct$coef)
# 
# cat("Transitional variance:", struct$coef["level"], "\n",
#     "Slope variance:", struct$coef["slope"], "\n",
#     "Observational variance:", struct$coef["epsilon"], "\n",
#     "Initial level of mu:", struct$model0$a[1], "\n",
#     "Initial level of lambda:", struct$model0$a[2], "\n",
#     "Initial level:", struct$model0$a, "\n")
# tsdiag(struct)
# 
# ###### local level (random walk plus noise)  - DLM ######
# buildModPoly1 <- function(v) {
#   dV <- exp(v[1])
#   dW <- exp(v[2])
#   m0 <- v[3]
#   dlmModPoly(1, dV=dV, dW=dW, m0=m0)
# }
# 
# # local linear trend
# buildModPoly2 <- function(v) {
#   dV <- exp(v[1])
#   dW <- exp(v[2:3])
#   m0 <- v[4:5]
#   dlmModPoly(order=2, dV=dV, dW=dW, m0=m0)
# }
# 
# varGuess <- var(diff(var), na.rm=TRUE)
# mu0Guess <- as.numeric(var[1])
# lambda0Guess <- 0.0
# 
# #parm <- c(log(varGuess), log(varGuess), log(varGuess))
# parm <- c(log(varGuess), log(varGuess), log(varGuess),
#           mu0Guess, lambda0Guess)
# #mle <- dlmMLE(var, parm=parm, buildModPoly1)
# mle <- dlmMLE(var, parm=parm, buildModPoly2)
# 
# if (mle$convergence != 0) stop(mle$message)
# 
# #model <- buildModPoly1(mle$par)
# model <- buildModPoly2(mle$par)
# model
# 
# filtered <- dlmFilter(var, model)
# tsdiag(filtered)
# smoothed <- dlmSmooth(filtered)
# 
# p <- autoplot(filtered, ts.linetype = 'dashed', fitted.colour = 'blue', ts.colour = 'green')
# p
# autoplot(smoothed, ts.linetype = 'dotted', ts.colour = 'red', p = p)
# plot(cbind(var, filtered$m[-1], smoothed$s[-1]), plot.type='s', col=c("black","red","blue"), ylab="Level", main="", lwd=c(1,2,2))
# 
# 
# 
# 
# ###########################################################################################
# #                  maSpeed ~ AltDiffRatio + maADR2 + maADR3                                   #
# ###########################################################################################
# 
# 
# plots_speed_by_altdiff <- Map(function (df) ggplot(df, aes(x = maAltDiffRatio, y = maSpeed)) + geom_point(color = 'green'), dfs)
# do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff, 'ncol' = 2, top = "maSpeed by maAltDiffRatio"))
# 
# plots_speed_by_altdiff2 <- Map(function (df) ggplot(df, aes(x = maADR2, y = maSpeed)) + geom_point(color = 'blue'), dfs)
# do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff2, 'ncol' = 2, top = "maSpeed by AltDiff2"))
# 
# plots_speed_by_altdiff3 <- Map(function (df) ggplot(df, aes(x = maADR3, y = maSpeed)) + geom_point(color = 'red'), dfs)
# do.call('grid.arrange', list('grobs' = plots_speed_by_altdiff3, 'ncol' = 2, top = "maSpeed by AltDiff3"))
# 
fits <- Map(function (df) lm(maSpeed ~ maAltDiffRatio + maADR2, data = df), dfs)
#fits <- Map(function (df) lm(maSpeed ~ AltDiffRatio + maADR2 + maADR3 + Distance2DTot, data = df), dfs)

summaries <- Map(summary, fits)
summaries



###########################################################################################
#                                  prepare route dataset for prediction                   #
###########################################################################################
r <- read_csv('route.csv') 
r
r <- r %>% filter(RouteNumber == 1) 
r <- r %>% mutate(DistanceLag = Distance - lag(Distance))

# need to smooth Altitude here already because way too noisy
r <- r %>% mutate(maAltitude = as.vector(unclass(mov(Altitude))))
r <- r %>% mutate(maAltitudeLag = Altitude - lag(Altitude))

r <- r %>% mutate(AltDiffRatio = maAltitudeLag/DistanceLag)

# important: downhill does not lead to speedup as uphill leads to slowdown!!!
r <- r %>% mutate(AltDiffRatio = ifelse(AltDiffRatio < 0, 0 , AltDiffRatio))

r <- r %>% mutate(maAltDiffRatio = as.vector(unclass(mov(AltDiffRatio))))
#ggplot(pred_df, aes(x = maAltDiffRatio)) + geom_histogram()
#ggplot(r, aes(x = maAltDiffRatio)) + geom_histogram()

r <- r %>% mutate(maADR2 = maAltDiffRatio^2, maADR3 = maAltDiffRatio^3, Distance2DTot = Distance)
r

# remove NAs for prediction
na_indices <- which(is.na(r$maAltDiffRatio))
na_indices_start <- na_indices[1:ceiling(window_size/2)]
na_indices_end <- setdiff(na_indices, na_indices_start)
start <- head(na.omit(r),1)$maAltDiffRatio
start
end <- tail(na.omit(r),1)$maAltDiffRatio
end
r[na_indices_start, "maAltDiffRatio"] <- start
r[na_indices_end, "maAltDiffRatio"] <- end

r <- r %>% mutate(maADR2 = maAltDiffRatio^2, maADR3 = maAltDiffRatio^3, Distance2DTot = Distance)
r

###########################################################################################
#                                lm: predictions on same dataset                          #
###########################################################################################


pred_df <- df529
pred_df <- na.omit(pred_df)
fit <- lm(maSpeed ~ maAltDiffRatio + maADR2 + maADR3 , data = pred_df)
summary(fit)
preds <- predict(fit, interval = 'prediction')
preds

df_plus_pred <- pred_df %>% mutate(p = preds[,1])

ggplot(df_plus_pred, aes(x=Distance2DTot)) + geom_line(aes(y=p), color='red') + geom_line(aes(y=Speed2D), color='cyan') +
  geom_line(aes(y=maSpeed), color='violet')



###########################################################################################
#                                  lm: predictions for route                              #
###########################################################################################

#pred_df <- df428
pred_df <- df529
#pred_df <- df

fit <- lm(maSpeed ~ maAltDiffRatio + maADR2 + maADR3 , data = pred_df)
summary(fit)

pred <- predict(fit, newdata = r, interval = 'prediction')

r_lm <- r %>% cbind(pred = pred[,1], lower = pred[,2], upper=pred[,3])
r_lm

ggplot(r_lm, aes(Distance2DTot, pred)) + geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'cyan', alpha = 0.2) +
  coord_cartesian(ylim = c(0,20)) + geom_smooth()


###########################################################################################
#                                  smoothing spline                                       #
###########################################################################################

pred_df <- df529
pred_df <- na.omit(pred_df)
fit <- smooth.spline(pred_df$maAltDiffRatio, pred_df$maSpeed, df = 5) # all.knots = TRUE
fit

preds <- predict(fit, pred_df$maAltDiffRatio)
preds
df_plus_pred <- pred_df %>% mutate(expected_speed = preds$y)
ggplot(df_plus_pred, aes(x=Distance2DTot)) + geom_line(aes(y=expected_speed), color='red') + geom_line(aes(y=Speed2D), color='cyan') +
  geom_line(aes(y=maSpeed), color='violet') + 
  ggtitle("Actual speed measured(cyan), actual speed smoothed (violet), and predicted speed (red) for runner 529")

preds <- predict(fit, r$maAltDiffRatio)
preds
df_plus_pred <- r %>% mutate(expected_speed = preds$y)
ggplot(df_plus_pred, aes(x=Distance2DTot)) + #geom_line(aes(y=expected_speed), color='red') + 
  geom_smooth(aes(y = expected_speed)) + 
  ggtitle("Expected speed on given route data (runner 529)")

(lat <- r[[1, 'Latitude']])
(long <- r[[1, 'Longitude']])
m <- get_map(location = c(lon = 7.26, lat = 46.4), zoom=12)
ggmap(m) + geom_point(data= df_plus_pred, aes(x=Longitude, y = Latitude, color=expected_speed)) + 
  scale_color_continuous(low = 'violet', high = 'yellow')+ 
  ggtitle("Speed development as predicted for runner 529 (given route data)")


###########################################################################################
#                                  svm                                                    #
###########################################################################################

pred_df <- df529
pred_df <- na.omit(pred_df)

fit <- svm(maSpeed ~ maAltDiffRatio, data=pred_df, kernel = 'polynomial', degree=2, cost = 10)
summary(fit)
preds <- predict(fit, interval = 'prediction')
preds
df_plus_pred <- pred_df %>% mutate(p = preds)
ggplot(df_plus_pred, aes(x=Distance2DTot, y=p)) + geom_point()

###########################################################################################
#                                lat/long                                                 #
###########################################################################################

# route <- read_csv('route.csv') 
# route
# route <- route %>% filter(RouteNumber == 1) %>% select(Longitude:Distance)
# route
# 
# ggplot(route, aes(Distance, Altitude)) + geom_line()
# 
# ggplot(route, aes(Distance, Altitude)) + geom_line() + 
#   geom_line(mapping = aes(x = Distance2DTot, y = Altitude), data = df69, color = 'blue') +
#   geom_line(mapping = aes(x = Distance2DTot, y = Altitude), data = df395, color = 'green') +
#   geom_line(mapping = aes(x = Distance2DTot, y = Altitude), data = df428, color = 'red') +
#   geom_line(mapping = aes(x = Distance2DTot, y = Altitude), data = df529, color = 'violet') 
# 
# ggplot(route, aes(Latitude, Longitude)) + geom_line() + 
#   geom_line(mapping = aes(Latitude, Longitude), data = df69, color = 'blue') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df395, color = 'green') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df428, color = 'yellow') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df529, color = 'violet') 
# 
# ggplot(route, aes(Latitude, Longitude)) + geom_line() + 
#   geom_line(mapping = aes(Latitude, Longitude), data = df69, color = 'blue') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df395, color = 'green') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df428, color = 'yellow') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df529, color = 'violet') +
#   coord_cartesian(xlim = c(46.320, 46.330), ylim = c(7.200, 7.225))
# 
# # zoom in further
# ggplot(route, aes(Latitude, Longitude)) + geom_line() + 
#   geom_line(mapping = aes(Latitude, Longitude), data = df69, color = 'blue') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df395, color = 'green') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df428, color = 'yellow') +
#   geom_line(mapping = aes(Latitude, Longitude), data = df529, color = 'violet') +
#   coord_cartesian(xlim = c(46.3250, 46.3255), ylim = c(7.200, 7.225))
# 

