# COVID-19 data analysis using county-level data from NYT, for a given state
# by Peter Wang 2020

library(RCurl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(lubridate)
library(zoo)

# Get county-level data; could take some seconds
df.raw = read.csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")

# Choose state
# Recommended parameter values
STATE = "Massachusetts"  ;c=7;k=2
# STATE = "New York"       ;c=4;k=1
# STATE = "New Jersey"     ;c=5;k=3
# STATE = "Washington"     ;c=3;k=1
# STATE = "Florida"        ;c=5;k=1
# STATE = "Michigan"       ;c=5;k=1
# STATE = "California"     ;c=5;k=1
# STATE = "Louisiana"      ;c=5;k=2
# STATE = "Pennsylvania"   ;c=5;k=2
# STATE = "Connecticut"    ;c=3;k=1

# Parameters
SAVEPLOTS = F  # save or view plots at the end
SMOOTH = T  # 3-day smoothing of plot or raw trends
COLNUM = c  # number of counties to show color and label for
k = k  # divide state total count by k to make graph easier to see


# Filters for N when calculating rates
CASEFILTER = 50  # rates calculated only when cases exceeed
DEATHFILTER = 10  # rates calculated only when deaths exceeed

# MA policies (1. stay-at-home, 2. face masks) with delay-day delay for cases; use geom_blank() for other states
delay = 10
if (STATE == "Massachusetts"){
  MApolicy1 = geom_vline(xintercept = ymd(20200324)+delay, linetype = "dashed", color = "salmon", alpha = 0.5, size = 0.8)
  MApolicy2 = geom_vline(xintercept = ymd(20200410)+delay, linetype = "dashed", color = "cadetblue", alpha = 0.5, size = 0.8)
}else{
  MApolicy1 = geom_blank()
  MApolicy2 = geom_blank()
}

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


# Color palette
genCol = function(coln){
  fullcol = c(
    "darkgoldenrod",
    "firebrick",
    "darkgreen",
    "mediumblue",
    "deeppink3",
    "chartreuse3",
    "darkslateblue",
    "darkorange2",
    "cadetblue4",
    "darkorchid3"
  )
  greys = rep("grey80",100)
  c("black",fullcol[1:coln],greys)
}

# Date range
startDate = ymd("20200301")
maxDate = max(ymd(df.raw$date))
x.ticks.halfm = scale_x_date(breaks = c(ymd("20200301","20200315","20200401","20200415","20200501","20200515","20200601","20200615")),
                             date_labels = "%b %d",
                             limits = c(startDate,maxDate))
if(SMOOTH){date.axis.smooth = "Date (3-day average)"}else{date.axis.smooth = "Date"}

# Filter and process data
df = df.raw %>%
  arrange(state, county, date) %>%
  filter(state == STATE) %>%  # filter for this state
  mutate(date = ymd(date)) %>%
  filter(date >= startDate) %>%  # filter data from too long ago
  select(-state) %>%
  group_by(county) %>%  # following couple lines arrange counties by decreasing orders of latest case number
  mutate(cases.max = max(cases)) %>%
  arrange(-cases.max) %>%
  ungroup() %>%
  mutate(county.r = county,
         county = factor(county.r, levels = unique(county.r))) %>%
  select(-county.r, -cases.max) %>%
  arrange(county)

df.STATE.sum = df %>%  # we need state-total reference as well
  select(date, cases, deaths) %>%
  group_by(date) %>%
  summarize_all(sum) %>%
  mutate(county = paste0("TOTAL/", k), cases = cases/k, deaths = deaths/k) %>%  # divide by k to make data easier to look at
  ungroup()

df = df %>%
  filter(county != "Unknown") %>%
  select(date, cases, deaths, county) %>%
  rbind(df.STATE.sum, .) %>%  # putting everything together
  mutate(county = factor(county, c(unique(county),"Others"))) %>%  # add Others to factor for color legend display
  mutate(county.o = factor(county, unique(arrange(df,desc(county))$county)))  # use this to order the graphs so the higher ones cover the lower ones, not the other way round

df.trends = df %>%
  # Calculate deltas
  mutate(
    cases.growth = abs(cases - lag(cases)),  # abs is a hack to prevent -ve number math error in log for doubling time; will be cleaned up by date filter
    deaths.growth = abs(deaths - lag(deaths)),
    cases.growth.rate = cases.growth/lag(cases),
    deaths.growth.rate = deaths.growth/lag(deaths)
  ) %>%
  # Smoothing
  mutate(
    cases.growth.smooth = rollapply(
      cases.growth, mean, width = 3, by = 1, align = "center", fill = NA  # rolling window averaging
    ),
    deaths.growth.smooth = rollapply(
      deaths.growth, mean, width = 3, by = 1, align = "center", fill = NA
    ),
    cases.growth.rate.smooth = rollapply(
      cases.growth.rate, mean, width = 3, by = 1, align = "center", fill = NA
    ),
    deaths.growth.rate.smooth = rollapply(
      deaths.growth.rate, mean, width = 3, by = 1, align = "center", fill = NA
    )) %>%
  # Calculate doubling times
  mutate(
    doublingtime = log(2)/log(1+if(SMOOTH){cases.growth.rate.smooth}else{cases.growth.rate}),
    ddoublingtime = log(2)/log(1+if(SMOOTH){deaths.growth.rate.smooth}else{deaths.growth.rate})
  ) %>%
  # The first and last of 3-day averaging is not relevant; +2 gets rid of the first one that is a delta from the max of previous county in the list
  group_by(county) %>%
    filter(
      if(SMOOTH){
        date >= min(date)+2 & date <= max(date)-1
      }else{
        date >= min(date)+1 & date <= max(date)
      }
    ) %>%
  ungroup()

# We are only id-ing and giving colors to the top several counties
countiesToShow = c(as.character(unique(filter(df, as.numeric(county) < 2+COLNUM)$county)), "Others")
colorscale = scale_colour_manual(values = genCol(coln = COLNUM), breaks = countiesToShow, drop = FALSE)

# Plot components to set up
# General lin plot geom
geom_line_plot  = geom_line(size = 0.5, alpha = 0.7, aes(group = county.o))
# Growth rate y-axis
y.lim.rate = coord_cartesian(ylim = c(0.01,0.75))
y.ticks.pc = scale_y_log10(labels = scales::percent,
                           breaks = c(0.01,0.02,0.03,0.05,0.10,0.20,0.30,0.50,0.75),
                           minor_breaks = NULL)
# Doubling time y-axis
DOUBLINGWKS = ((unname(max(
  filter(
    df.trends,
    cases >= CASEFILTER,
    county == paste0("TOTAL/", k)
  )$doublingtime
))*1.2) %/% 7) + 1
y.lim.wks = coord_cartesian(ylim = c(0, 7*DOUBLINGWKS))
y.ticks.weeks = scale_y_continuous(breaks = seq(0,7*DOUBLINGWKS,7))

########################################
# We start graphing here:
########################################
# CASE TRENDS AND DAILY COUNT

# Graph #1: case trends
G1 = df %>%
  ggplot(aes(x = date, y = cases, color = county)) +
  geom_line_plot +
  scale_y_log10(limits = c(1,NA)) +
  x.ticks.halfm +
  ggtitle("Confirmed cases (log10)") +
  xlab("Date") +
  MApolicy1 + MApolicy2 +
  colorscale +
  theme1

# Graph #2: daily case counts, with 3-day smoothing
G2 = df.trends %>%
  ggplot(aes(x = date, y = if(SMOOTH){cases.growth.smooth}else{cases.growth}, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle("Daily case growth") +
  xlab(date.axis.smooth) +
  MApolicy1 + MApolicy2 +
  colorscale +
  theme1

# grid.arrange(G1, G2, nrow = 1)

########################################
# DEATH TRENDS AND DAILY COUNT

# Graph #3: death trends
G3 = df %>%
  ggplot(aes(x = date, y = deaths, color = county)) +
  geom_line_plot +
  scale_y_log10(limits = c(1,NA)) +
  x.ticks.halfm +
  ggtitle("Confirmed deaths (log10)") +
  xlab("Date") +
  colorscale +
  theme1

# Graph #4: daily case counts, with 3-day smoothing
G4 = df.trends %>%
  ggplot(aes(x = date, y = if(SMOOTH){deaths.growth.smooth}else{deaths.growth}, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle("Daily deaths growth") +
  xlab(date.axis.smooth) +
  colorscale +
  theme1

# grid.arrange(G3, G4, nrow = 1)

################################################
# DAILY GROWTH RATE

# Graph #5: daily case growth rate, with 3-day smoothing
G5 = df.trends %>%
  filter(cases >= CASEFILTER) %>%
  ggplot(aes(x = date, y = if(SMOOTH){cases.growth.rate.smooth}else{cases.growth.rate}, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle(paste0("Daily cases growth rate", " (N > ", CASEFILTER, ", log-scale)")) +
  xlab(date.axis.smooth) +
  y.lim.rate +
  y.ticks.pc +
  MApolicy1 + MApolicy2 +
  colorscale +
  theme1

# Graph #6: daily death growth rate, with 3-day smoothing
G6 = df.trends %>%
  filter(deaths >= DEATHFILTER) %>%
  ggplot(aes(x = date, y = if(SMOOTH){deaths.growth.rate.smooth}else{deaths.growth.rate}, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle(paste0("Daily deaths growth rate", " (N > ", DEATHFILTER, ", log-scale)")) +
  xlab(date.axis.smooth) +
  y.lim.rate +
  y.ticks.pc +
  colorscale +
  theme1

# grid.arrange(G5, G6, nrow = 1)

#############################
# DOUBLING TIME

# Graph #7: cases doubling time, with 3-day smoothing
G7 = df.trends %>%
  filter(cases >= CASEFILTER) %>%
  ggplot(aes(x = date, y = doublingtime, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle(paste0("Confirmed cases doubling time (days)", " (N > ", CASEFILTER, ")")) +
  xlab(date.axis.smooth) +
  y.lim.wks +
  y.ticks.weeks +
  MApolicy1 + MApolicy2 +
  colorscale +
  theme1

# Graph #8: deaths doubling time, with 3-day smoothing
G8 = df.trends %>%
  filter(deaths >= DEATHFILTER) %>%
  ggplot(aes(x = date, y = ddoublingtime, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle(paste0("Confirmed deaths doubling time (days)", " (N > ", DEATHFILTER, ")")) +
  xlab(date.axis.smooth) +
  y.lim.wks +
  y.ticks.weeks +
  colorscale +
  theme1

# grid.arrange(G7, G8, nrow = 1)

#####################

# Plot time!
titleobj = textGrob(paste0(STATE, " (up to ", maxDate, ")"), gp = gpar(fontface = "bold", fontsize = 24))

if(!SAVEPLOTS){
  grid.arrange(G1,G3,G2,G4,G5,G6,G7,G8, nrow = 4, top = titleobj)  # for viewing
}else{
  G = arrangeGrob(G1,G3,G2,G4,G5,G6,G7,G8, nrow = 4, top = titleobj)
  ggsave(file = paste0("covid19_", STATE, ".pdf"), G, width = 12, height = 12, units = "in")
  ggsave(file = paste0("covid19_", STATE, ".png"), G, width = 12, height = 12, units = "in", dpi = 320)
}
