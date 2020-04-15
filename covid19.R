# COVID-19 data analysis using county-level data from NYT, for a given state
# by Peter Wang 2020

library(RCurl)
library(tidyr)
library(dplyr)
library(nlstools)
library(ggplot2)
library(grid)
library(gridExtra)
library(lubridate)

# Get county-level data; could take some seconds
# df.raw = read.csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")

# Choose state
# Recommended parameter values
STATE = "Massachusetts"  ;c=7;k=2;d=2
# STATE = "New York"       ;c=4;k=1;d=4
# STATE = "New Jersey"     ;c=5;k=3;d=3
# STATE = "Washington"     ;c=3;k=1;d=8
# STATE = "Florida"        ;c=5;k=1;d=4
# STATE = "Michigan"       ;c=5;k=1;d=4
# STATE = "California"     ;c=5;k=1;d=5
# STATE = "Louisiana"      ;c=5;k=2;d=6
# STATE = "Pennsylvania"   ;c=5;k=2;d=4
# STATE = "Connecticut"    ;c=3;k=1;d=3

# Parameters
COLNUM = c  # number of counties to show color and label for
k = k  # divide total count by k to make graph easier to see
DOUBLINGWKS = d  # 2/3/4/5/6/etc. weeks as upper limit of doubling time graphs

# Filters for N when calculating rates
CASEFILTER = 50  # rates calculated only when cases exceeed
DEATHFILTER = 10  # rates calculated only when deaths exceeed

# MA policies (1. stay-at-home, 2. face masks) with delay-day delay for cases; use geom_blank() for other states
delay = 10
# MApolicy1 = geom_vline(xintercept = ymd(20200324)+delay, linetype = "dashed", color = "salmon", alpha = 0.5, size = 0.8)
# MApolicy2 = geom_vline(xintercept = ymd(20200410)+delay, linetype = "dashed", color = "cadetblue", alpha = 0.5, size = 0.8)
MApolicy1 = geom_blank()
MApolicy2 = geom_blank()


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
  select(date, cases, deaths, county) %>%
  rbind(df.STATE.sum, .) %>%  # putting everything together
  mutate(county = factor(county, c(unique(county),"Others"))) %>%  # add Others to factor for color legend display
  mutate(county.o = factor(county, unique(arrange(df,desc(county))$county)))  # use this to order the graphs so the higher ones cover the lower ones, not the other way round

# We are only id-ing and giving colors to the top several counties
countiesToShow = c(as.character(unique(filter(df, as.numeric(county) < 2+COLNUM)$county)), "Others")
colorscale = scale_colour_manual(values = genCol(coln = COLNUM), breaks = countiesToShow, drop = FALSE)

# This the key geom we are plotting every time
geom_line_plot  = geom_line(size = 0.5, alpha = 0.7, aes(group = county.o))


########################################
# We start graphing here:
########################################
# CASE TRENDS AND DAILY COUNT

# Graph #1: case trends
G1 = df %>%
  filter(county != "Unknown") %>%
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
df.growth = df %>%
  mutate(cases.growth = cases - lag(cases)) %>%
  filter(!is.na(cases.growth) & cases.growth >= 0) %>%
  mutate(cases.growth.smooth = zoo::rollapply(
    cases.growth, mean, width = 3, by = 1, align = "center", fill = NA
  )) %>%
  group_by(county) %>%
  filter(date > min(date) & date < max(date)) %>%
  ungroup()
G2 = df.growth %>%
  filter(county != "Unknown") %>%
  ggplot(aes(x = date, y = cases.growth.smooth, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle("Daily case growth") +
  xlab("Date (3-day average)") +
  MApolicy1 + MApolicy2 +
  colorscale +
  theme1

# grid.arrange(G1, G2, nrow = 1)

########################################
# DEATH TRENDS AND DAILY COUNT

# Graph #3: death trends
G3 = df %>%
  filter(county != "Unknown") %>%
  ggplot(aes(x = date, y = deaths, color = county)) +
  geom_line_plot +
  scale_y_log10(limits = c(1,NA)) +
  x.ticks.halfm +
  ggtitle("Confirmed deaths (log10)") +
  xlab("Date") +

  colorscale +
  theme1

# Graph #4: daily case counts, with 3-day smoothing
df.dgrowth = df %>%
  mutate(deaths.growth = deaths - lag(deaths)) %>%
  filter(!is.na(deaths.growth) & deaths.growth >= 0) %>%
  mutate(deaths.growth.smooth = zoo::rollapply(
    deaths.growth, mean, width = 3, by = 1, align = "center", fill = NA
  )) %>%
  group_by(county) %>%
  filter(date > min(date) & date < max(date)) %>%
  ungroup()
G4 = df.dgrowth %>%
  filter(county != "Unknown") %>%
  ggplot(aes(x = date, y = deaths.growth.smooth, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle("Daily deaths growth") +
  xlab("Date (3-day average)") +
  colorscale +
  theme1

# grid.arrange(G3, G4, nrow = 1)

################################################
# DAILY GROWTH RATE

y.lim.rate = coord_cartesian(ylim = c(0,0.8))
y.ticks.pc = scale_y_continuous(labels = scales::percent)

# Graph #5: daily case growth rate, with 3-day smoothing
df.growth.rate = df %>%
  mutate(cases.growth.rate = (cases-lag(cases))/lag(cases)) %>%
  filter(cases >= CASEFILTER) %>%
  filter(!is.na(cases.growth.rate) & cases.growth.rate >= 0) %>%
  mutate(cases.growth.rate.smooth = zoo::rollapply(
    cases.growth.rate, mean, width = 3, by = 1, align = "center", fill = NA
  )) %>%
  group_by(county) %>%
  filter(date > min(date) & date < max(date)) %>%
  ungroup() %>%
  mutate(doublingtime = log(2)/log(1+cases.growth.rate.smooth))
G5 = df.growth.rate %>%
  filter(county != "Unknown") %>%
  ggplot(aes(x = date, y = cases.growth.rate.smooth, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle(paste0("Daily cases growth rate", " (N > ", CASEFILTER, ")")) +
  xlab("Date (3-day average)") +
  y.lim.rate +
  y.ticks.pc +
  MApolicy1 + MApolicy2 +
  colorscale +
  theme1

# Graph #6: daily death growth rate, with 3-day smoothing
df.dgrowth.rate = df %>%
  mutate(deaths.growth.rate = (deaths-lag(deaths))/lag(deaths)) %>%
  filter(deaths >= DEATHFILTER) %>%
  filter(!is.na(deaths.growth.rate) & deaths.growth.rate >= 0) %>%
  mutate(deaths.growth.rate.smooth = zoo::rollapply(
    deaths.growth.rate, mean, width = 3, by = 1, align = "center", fill = NA
  )) %>%
  group_by(county) %>%
  filter(date > min(date) & date < max(date)) %>%
  ungroup() %>%
  mutate(ddoublingtime = log(2)/log(1+deaths.growth.rate.smooth))
G6 = df.dgrowth.rate %>%
  filter(county != "Unknown") %>%
  ggplot(aes(x = date, y = deaths.growth.rate.smooth, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle(paste0("Daily deaths growth rate", " (N > ", DEATHFILTER, ")")) +
  xlab("Date (3-day average)") +
  y.lim.rate +
  y.ticks.pc +
  colorscale +
  theme1

# grid.arrange(G5, G6, nrow = 1)

#############################
# DOUBLING TIME

y.ticks.weeks = scale_y_continuous(breaks = c(0,7,14,21,28,35,42,49,56))
y.lim.wks = coord_cartesian(ylim = c(0, 7*DOUBLINGWKS))

# Graph #7: cases doubling time, with 3-day smoothing
G7 = df.growth.rate %>%
  filter(county != "Unknown") %>%
  ggplot(aes(x = date, y = doublingtime, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle(paste0("Confirmed cases doubling time (days)", " (N > ", CASEFILTER, ")")) +
  xlab("Date (3-day average)") +
  y.lim.wks +
  y.ticks.weeks +
  MApolicy1 + MApolicy2 +
  colorscale +
  theme1

# Graph #8: deaths doubling time, with 3-day smoothing
G8 = df.dgrowth.rate %>%
  filter(county != "Unknown") %>%
  ggplot(aes(x = date, y = ddoublingtime, color = county)) +
  geom_line_plot +
  x.ticks.halfm +
  ggtitle(paste0("Confirmed deaths doubling time (days)", " (N > ", DEATHFILTER, ")")) +
  xlab("Date (3-day average)") +
  y.lim.wks +
  y.ticks.weeks +
  colorscale +
  theme1

# grid.arrange(G7, G8, nrow = 1)

#####################

# Plot time!
# grid.arrange(G1,G3,G2,G4,G5,G6,G7,G8, nrow = 4, top = textGrob(STATE, gp = gpar(fontface = "bold", fontsize = 24)))  # for viewing
G = arrangeGrob(G1,G3,G2,G4,G5,G6,G7,G8, nrow = 4, top = textGrob(STATE, gp = gpar(fontface = "bold", fontsize = 24)))
ggsave(file = paste0("covid19_", STATE, ".pdf"), G, width = 12, height = 12, units = "in")
ggsave(file = paste0("covid19_", STATE, ".png"), G, width = 12, height = 12, units = "in", dpi = 320)
