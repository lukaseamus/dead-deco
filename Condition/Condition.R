#### Dead deco: kelp decomposition without physiology ####
#### Luka Seamus Wright                               ####

# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

# 3rd February 2025 was the start date of the experiment
start <- dmy_hm("3.2.25 13:00")

sal <- here("Condition", "Salinity.csv") %>%
  read_csv() %>%
  mutate(Date = dmy(Date),
         Day = start %>% as_date() %--% Date / ddays()) %T>%
  print()

PAR <- here("Condition", "PAR.csv") %>%
  read_csv(col_types = list("f")) %>%
  mutate(Date = dmy(Date),
         Day = start %>% as_date() %--% Date / ddays()) %T>%
  print()

temp <- here("Condition", "Temperature") %>%
  list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Tank = Path %>% basename() %>% 
      str_extract("\\d+") %>% fct(),
    Data = Path %>% 
      map(
        ~.x %>% read_csv(skip = 1, col_select = 2:3) %>%
          rename(Datetime = 1, Temperature = 2) %>%
          mutate(Datetime = Datetime %>%
                   parse_date_time( # Varying date format
                     orders = c("mdy HM", "mdy IMS p")
                   ),
                 Day = start %--% Datetime / ddays())
      )
  ) %>%
  select(-Path) %>%
  unnest(Data) %T>%
  print()

# 1.2 Visualise data ####
# Date labelling function
date_label <- function(x) {
  scales::label_date("%e. %b")( start + ddays(x) )
}

sal %>%
  ggplot(aes(Day, Salinity)) +
    geom_point() +
    scale_x_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, 20),
      sec.axis = dup_axis(
        name = NULL,
        labels = date_label
      )
    ) +
    theme_minimal()

PAR %>%
  ggplot(aes(Day, PAR)) +
    geom_point() +
    scale_x_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, 20),
      sec.axis = dup_axis(
        name = NULL,
        labels = date_label
      )
    ) +
    facet_grid(~ Tank) +
    theme_minimal()

PAR %>%
  group_by(Day, Tank) %>%
  summarise(PAR_mean = mean(PAR),
            PAR_sd = sd(PAR)) %>%
  ggplot() +
    geom_pointrange(
      aes(Day, PAR_mean, ymin = PAR_mean - PAR_sd,
          ymax = PAR_mean + PAR_sd)
    ) +
    scale_x_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, 20),
      sec.axis = dup_axis(
        name = NULL,
        labels = date_label
      )
    ) +
    facet_grid(~ Tank) +
    theme_minimal()

temp %>%
  ggplot(aes(Day, Temperature)) +
    geom_line(linewidth = 0.3) +
    scale_x_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, 20),
      sec.axis = dup_axis(
        name = NULL,
        labels = date_label
      )
    ) +
    facet_grid(~ Tank) +
    theme_minimal()
# No need to analyse. Simply summarise and make 
# a figure with light and temperature.

# 2. Calculate summaries ####
# 2.1 Salinity ####
sal %>%
  summarise(mean = mean(Salinity),
            sd = sd(Salinity),
            median = median(Salinity),
            n = n())
# mean    sd median     n
# 35.2 0.142   35.2    19

# 2.2 PAR ####
# PAR across tanks
PAR %>%
  summarise(mean = mean(PAR),
            sd = sd(PAR),
            median = median(PAR),
            n = n())
# mean    sd median     n
# 14.5  2.95     15   660

# PAR per tank
PAR %>%
  group_by(Tank) %>%
  summarise(mean = mean(PAR),
            sd = sd(PAR),
            median = median(PAR),
            n = n())
# Tank   mean    sd median     n
# 1      14.1  3.05     14   220
# 2      13.9  2.94     14   220
# 3      15.4  2.62     16   220

# Pre-averaged PAR per tank
PAR %>%
  group_by(Tank, Date) %>%
  summarise(PAR = mean(PAR)) %>%
  group_by(Tank) %>%
  summarise(mean = mean(PAR),
            sd = sd(PAR),
            median = median(PAR),
            n = n())
# Tank   mean    sd median     n
# 1      14.1  2.54   13.4    22
# 2      13.9  2.46   13.6    22
# 3      15.4  2.16   16.2    22

# 2.3 Temperature ####
# Temperature across tanks
temp %>%
  summarise(mean = mean(Temperature),
            sd = sd(Temperature),
            median = median(Temperature),
            max = max(Temperature),
            n = n())
# mean    sd median   max     n
# 23.8 0.818   23.8  25.9 62916

# Temperature per tank
temp %>%
  group_by(Tank) %>%
  summarise(mean = mean(Temperature),
            sd = sd(Temperature),
            median = median(Temperature),
            max = max(Temperature),
            n = n())
# Tank   mean    sd median   max     n
# 1      23.8 0.818   23.8  25.9 20972
# 2      23.9 0.817   23.8  25.9 20972
# 3      23.8 0.820   23.8  25.8 20972

# Temperature at start and end
end <- dmy("17.4.25") # end of experiment
temp %>%
  mutate(Date = Datetime %>% as_date()) %>%
  # Temperature logging started after start
  # of experiment so first day is used instead
  filter(Date %in% c(min(Date), end)) %>%
  group_by(Date) %>%
  summarise(mean = mean(Temperature),
            sd = sd(Temperature),
            median = median(Temperature),
            max = max(Temperature),
            n = n())
# Date        mean     sd median   max     n
# 2025-02-04  25.5 0.0640   25.5  25.6   384
# 2025-04-17  22.6 0.128    22.6  22.8   864

temp %>%
  mutate(Date = Datetime %>% as_date()) %>%
  # Temperature logging started after start
  # of experiment so first day is used instead
  filter(Date %in% c(min(Date), end)) %>%
  group_by(Tank, Date) %>%
  summarise(mean = mean(Temperature),
            sd = sd(Temperature),
            median = median(Temperature),
            max = max(Temperature),
            n = n())
# Tank  Date        mean     sd median   max     n
# 1     2025-02-04  25.4 0.0524   25.5  25.5   128
# 1     2025-04-17  22.6 0.124    22.7  22.8   288
# 2     2025-02-04  25.5 0.0615   25.5  25.6   128
# 2     2025-04-17  22.7 0.123    22.7  22.8   288
# 3     2025-02-04  25.4 0.0532   25.5  25.5   128
# 3     2025-04-17  22.6 0.122    22.6  22.8   288

# 3. Figure S2 ####
# 3.1 Figure S2a ####
# Define custom theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing.x = unit(1, "cm"),
                 panel.spacing.y = unit(0.6, "cm"),
                 text = element_text(family = "Futura"))

Fig_S2a <- temp %>%
  ggplot(aes(Day, Temperature)) +
    geom_line(linewidth = 0.3) +
    scale_x_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, 20),
      sec.axis = dup_axis(
        name = NULL,
        breaks = seq(0, 80, 40),
        labels = date_label
      )
    ) +
    facet_grid(~ Tank) +
    labs(x = "Experimental day",
         y = "Temperature (°C)") +
    coord_cartesian(ylim = c(21, 26),
                    expand = F, clip = "off") +
    mytheme +
    theme(strip.text = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.x.bottom = element_blank(),
          # Vectorisation of hjust works but causes overlap which
          # needs to be removed afterwards
          axis.text.x.top = element_text(hjust = c(0, 0.5, 1)))

Fig_S2a

# 3.2 Figure S2b ####
Fig_S2b <- PAR %>%
  group_by(Day, Tank) %>%
  summarise(PAR_mean = mean(PAR),
            PAR_sd = sd(PAR)) %>%
  ggplot(aes(Day, PAR_mean)) +
     geom_pointrange(
      aes(ymin = PAR_mean - PAR_sd,
          ymax = PAR_mean + PAR_sd),
      shape = 16, size = 0.5
    ) +
    scale_x_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, 20),
      sec.axis = dup_axis(
        name = NULL,
        breaks = seq(0, 80, 40),
        labels = date_label
      )
    ) +
    labs(x = "Experimental day",
         y = expression("PAR (µmol photons m"^-2*" s"^-1*")")) +
    facet_grid(~ Tank) +
    coord_cartesian(ylim = c(5, 20),
                    expand = F, clip = "off") +
    mytheme +
    theme(strip.text = element_blank(),
          axis.text.x.top = element_blank(),
          axis.line.x.top = element_blank(),
          axis.ticks.x.top = element_blank())

Fig_S2b

# 3.3 Combine panels ####
require(patchwork)
Fig_S2 <- Fig_S2a / Fig_S2b

Fig_S2 %>%
  ggsave(filename = "Fig_S2.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 20, units = "cm")