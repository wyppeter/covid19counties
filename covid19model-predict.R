# COVID-19 data analysis using county-level data from NYT
# Models the daily rate of case number growth for a number of states,
# based on an exponential decrease of the rate
# Note: This model is based on an empirical observation of trend shape, not a theoretically supported mathematical model.
# by Peter Wang 2020

library(RCurl)
library(tidyr)
library(dplyr)
library(nlstools)
library(ggplot2)
library(grid)
library(gridExtra)
library(lubridate)
library(zoo)
library(FME)
library(ggrepel)

# Get county-level data; could take some seconds
df.raw = read.csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")

# Choose states
STATES = c(
  "Massachusetts",
  "New York",
  "New Jersey",
  "Washington",
  "Florida",
  "Michigan",
  "California",
  "Louisiana",
  "Pennsylvania",
  "Connecticut",
  "Illinois",
  "Texas",
  "Maryland",
  "Tennessee",
  "Colorado"
)

# # Parameters
SAVEPLOTS = F  # save or view plots at the end


# Filters for N when calculating rates
CASEFILTER = 50  # rates calculated only when cases exceeed
DEATHFILTER = 10  # rates calculated only when deaths exceeed

# Plot components to set up
# General lin plot geom
geom_line_plot  = geom_line(size = 0.5, alpha = 0.7, color = "grey50")
# Growth rate y-axis
y.lim.rate = coord_cartesian(ylim = c(0.01,0.75))
y.ticks.pc = scale_y_log10(labels = scales::percent,
                           breaks = c(0.01,0.02,0.03,0.05,0.10,0.20,0.30,0.50,0.75),
                           minor_breaks = NULL)
colormodel = scale_colour_manual(values = c("black", "mediumblue", "firebrick"), name = "Prediction model")

# Set ggplot2 theme
theme1 = theme(
  panel.grid.major = element_line(colour = "gray85"),
  panel.grid.minor = element_line(colour = "gray90"),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  axis.text = element_text(color = "black", size = 11),
  axis.title = element_text(color = "black", size = 12),
  axis.title.y = element_blank(),
  plot.title = element_text(color = "black", size = 12),
  legend.background = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(color = "black", size = 10),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 11, colour = "black"),
  strip.text.y = element_text(size = 11, colour = "black", angle = 0)
)

# Date ranges
startDate = ymd("20200301")
maxDate = max(ymd(df.raw$date))
x.ticks.halfm = scale_x_date(breaks = c(ymd(
  "20200301","20200315",
  "20200401","20200415",
  "20200501","20200515",
  "20200601","20200615"
)),
date_labels = "%b %d",
limits = c(startDate,maxDate))


# Filter and process data
df = df.raw %>%
  filter(state %in% STATES) %>%  # filter states to look at
  arrange(state, county, date) %>%  # tidy things up
  mutate(date = ymd(date)) %>%
  select(state, county, date, cases, deaths) %>%
  group_by(state, date) %>%  # tally total of each state
  summarize_at(c("cases", "deaths"), sum, na.rm = T) %>%
  ungroup() %>%
  filter(date >= startDate) %>%  # filter data from too long ago
  group_by(state) %>%  # following couple lines arrange states by decreasing orders of latest case number
  mutate(cases.max = max(cases)) %>%
  arrange(-cases.max) %>%
  ungroup() %>%
  mutate(state.r = state,
         state = factor(state.r, levels = unique(state.r))) %>%
  select(-state.r, -cases.max) %>%
  arrange(state)

# Calculate deltas; NOTE: This is a different implementation from the main covid19.R script, which calculates delta with before, not ahead!!
df.trends = df %>%
  mutate(
    cases.growth = abs(lead(cases) - cases),  # abs is a hack to prevent -ve number math error in log for doubling time; will be cleaned up by date filter
    deaths.growth = abs(lead(deaths) - deaths),
    cases.growth.rate = cases.growth/cases,
    deaths.growth.rate = deaths.growth/deaths
  ) %>%
  # Calculate doubling times
  mutate(
    doublingtime = log(2)/log(1+cases.growth.rate),
    ddoublingtime = log(2)/log(1+deaths.growth.rate)
  ) %>%
  # Get rid of the last one that is a delta from the max of next state in the list
  group_by(state) %>%
  filter(date < max(date)) %>%
  ungroup()

################################################
# *!* Modeling section
# The daily rate of cases growth appear to be reducing -ve exponentially

# Dates for modeling purposes
modelStart = ymd("20200328")  # Model only when data starts to exhibit the observed trend
modelZero = maxDate - 2  # ODE origin is today, 2 days ago (allowing 5d-average)
daytln = as.numeric(modelZero - modelStart)

# Further work on data frame
df.rate = df.trends %>%
  mutate(cases.growth.rate.log = log10(cases.growth.rate)) %>%
  filter(cases >= CASEFILTER) %>%
  select(state, date, cases, cases.growth, cases.growth.rate, cases.growth.rate.log) %>%
  mutate(days.since = as.numeric(date - modelZero)) %>%  # use days instead data objects to fit as x
  filter(date >= modelStart)

# LM fitting
lmconf = 0.9  # 90% confidence interval
df.lm.o = df.rate %>%
  group_by(state) %>%
  summarise(
    rsq    = summary(lm(formula = cases.growth.rate.log ~ days.since))$r.squared,

    beta   = unname(coef(lm(formula = cases.growth.rate.log ~ days.since))[2]),
    yint   = unname(coef(lm(formula = cases.growth.rate.log ~ days.since))[1]),

    beta.b = unname(confint(lm(formula = cases.growth.rate.log ~ days.since), level = lmconf)[,1])[2],
    yint.b = unname(confint(lm(formula = cases.growth.rate.log ~ days.since), level = lmconf)[,1])[1],

    beta.w = unname(confint(lm(formula = cases.growth.rate.log ~ days.since), level = lmconf)[,2])[2],
    yint.w = unname(confint(lm(formula = cases.growth.rate.log ~ days.since), level = lmconf)[,2])[1]
  ) %>%
  arrange(state)

# Plot the fit
Glmfit = df.rate %>%
  ggplot(aes(x = days.since, y = cases.growth.rate.log)) +
  geom_line(color = "deeppink") +
  geom_smooth(method = "lm", formula = y ~ x, fullrange = T, se = T, color = "black") +

  geom_abline(data = df.lm.o, aes(slope = beta.b, intercept = yint.b), color = "blue") +
  geom_abline(data = df.lm.o, aes(slope = beta.w, intercept = yint.w), color = "red") +

  geom_text(inherit.aes = FALSE,
            data = df.lm.o,
            size = 3, x = -10, parse = F,
            aes(
              y = yint - 1.0,
              label = paste0(
                "log10(rate) = (",
                as.character(format(round(as.numeric(beta), 5), nsmall = 5)),
                ") * days + (",
                as.character(format(round(as.numeric(yint), 5), nsmall = 5)),
                ")"))) +
  geom_text(inherit.aes = FALSE,
            data = df.lm.o,
            size = 3, x = -10, parse = F,
            aes(
              y = yint - 1.1,
              label = paste0("R^2 = ", as.character(format(round(as.numeric(rsq), 3), nsmall = 3))))) +

  ggtitle("Fitting confirmed cases growth rate (log10)") +
  xlab(paste0("Days since ", modelZero)) +
  # xlim(c(-10, 25)) +
  facet_wrap(vars(state)) +
  theme1

Glmfit # view

#####
# Looks good? Now we model the ODE
# Parameters
backwarddays = as.numeric(modelZero - modelStart) + 7
modeldays = 100
day0estim.df = df %>%
  group_by(state) %>%
  filter(date <= modelZero + 2 & date >= modelZero - 2) %>%  # 5d-average in log scale
  mutate(caseslog10 = log10(cases))
df.zero = day0estim.df %>%
  summarize_at(vars(caseslog10), mean) %>%
  mutate(day0cases = 10^caseslog10) %>%
  select(state, day0cases) %>%
  arrange(state) %>%
  select(-state)

df.lm = cbind(df.lm.o, df.zero)

# ODE
DEmodel = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dcases = cases * (10^( b * days.since + a ))
    ddays.since = 1
    list(c(dcases, ddays.since))})}
DEmodel.b = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dcases = -cases * (10^( -b * (days.since) + a ))
    ddays.since = 1
    list(c(dcases, ddays.since))})}
# Days to model
t.f = seq(0, modeldays, by = 1)
t.b = seq(0, backwarddays, by = 1)
# ODE function
oderun = function(day0, yint, beta, mdl){

  p = c(a = yint, b = beta)

  # Forward
  tstate <- c(day0, 0)
  names(tstate) <- c("cases", "days.since")
  df.model.f = ode(y = tstate, times = t.f, func = DEmodel, parms = p) %>%
    as.data.frame() %>%
    mutate(cases.growth = lead(cases) - cases) %>%
    filter(!is.na(cases.growth)) %>%
    mutate(mdl = mdl)

  # Backward
  tstate <- c(day0, 0)
  names(tstate) <- c("cases", "days.since")
  df.model.b = ode(y = tstate, times = t.b, func = DEmodel.b, parms = p) %>%
    as.data.frame() %>%
    mutate(time = -time, days.since = -days.since) %>%
    arrange(time) %>%
    mutate(cases.growth = lead(cases) - cases) %>%
    filter(!is.na(cases.growth)) %>%
    mutate(mdl = mdl)

  # Combine
  rbind(df.model.b, df.model.f)
}
# Run the ODE function for each model scenario
df.model = data.frame(row.names = F)  # initialize to prevent cumulative rbind
for(i in 1:length(STATES)){
  s = STATES[i]
  dt.params = df.lm %>%
    filter(state == s)
  model.med = oderun(
    day0 = unname(unlist(dt.params$day0cases)),
    yint = unname(unlist(dt.params$yint)),
    beta = unname(unlist(dt.params$beta)),
    mdl = "best-fit")
  model.good = oderun(
    day0 = unname(unlist(dt.params$day0cases)),
    yint = unname(unlist(dt.params$yint.b)),
    beta = unname(unlist(dt.params$beta.b)),
    mdl = "optimistic")
  model.bad = oderun(
    day0 = unname(unlist(dt.params$day0cases)),
    yint = unname(unlist(dt.params$yint.w)),
    beta = unname(unlist(dt.params$beta.w)),
    mdl = "conservative")
  df.model.s = rbind(model.med,model.good,model.bad) %>%
    mutate(
      mdl = factor(mdl, levels = c("best-fit", "optimistic", "conservative")),
      state = s
    ) %>%
    ungroup() %>%
    select(state, mdl, time, days.since, cases, cases.growth)
  df.model = rbind(df.model, df.model.s)
}

# Pinpoint the peak day
df.peak = df.model %>%
  group_by(state, mdl) %>%
  mutate(dgrowth = lead(cases.growth) - cases.growth) %>%
  filter(dgrowth > 0, time != modeldays-1) %>%
  filter(dgrowth == min(dgrowth)) %>%
  ungroup()

# Plot the modeled curves only
GM1 = df.model %>%
  ggplot(aes(x = modelZero + time, y = cases, color = mdl)) +
  geom_line() +
  scale_y_log10(limits = c(1,NA)) +
  ggtitle("Modeled cases (log10)") +
  xlab("Date") +
  facet_wrap(vars(state), scales = "free") +
  colormodel +
  theme1

GM2 = df.model %>%
  ggplot(aes(x = modelZero + time, y = cases.growth, color = mdl)) +
  geom_line() +
  ggtitle("Modeled cases growth") +
  xlab("Date") +
  facet_wrap(vars(state), scales = "free") +
  colormodel +
  theme1

grid.arrange(GM1, GM2, nrow = 2) # view

######################

# Overlay data
# Axis labels
x.ticks.halfm.model = scale_x_date(
  breaks = c(ymd(
    "20200301",
    "20200401",
    "20200501",
    "20200601",
    "20200701",
    "20200801",
    "20200901"
  )),
  date_labels = "%b",
  limits = c(startDate, modelZero + modeldays))

# CASE TRENDS
GMDL = df.trends %>%
  ggplot(aes(x = date, y = cases.growth)) +
  geom_line_plot +
  x.ticks.halfm.model +
  geom_line(data = df.model,
            aes(x = time + modelZero, color = mdl)) +
  geom_text_repel(data = df.peak,
            inherit.aes = F,
            size = 3,
            alpha = 0.5,
            aes(color = mdl,
                y = cases.growth,
                x = modelZero + time,
                label = paste0(
                  format(time + modelZero, format="%m/%d")
                )
                )
            ) +
  ggtitle("Predicted daily case growth") +
  xlab("Date") +
  facet_wrap(vars(state), scales = "free", ncol = 5) +
  colormodel +
  theme1

#####################

# Plot time!
if(!SAVEPLOTS){
  GMDL  # for viewing
}else{
  ggsave(file = paste0("covid19_model_", maxDate, ".pdf"), GMDL, width = 18, height = 10, units = "in")
  ggsave(file = paste0("covid19_model_", maxDate, ".png"), GMDL, width = 18, height = 10, units = "in", dpi = 320)
}
