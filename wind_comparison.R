# Wind data comparison 

library(worldmet)
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(readr)
library(tidyr)
library(stringr)
library(rlang)
library(patchwork)
library(RColorBrewer)
library(scales)

# Data --------------------

## Collection windows ----------
collections <- read_csv("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/metadata/collections.csv")

# Parse existing table of windows in UTC
windows <- collections |>
  rename(label = Experiment, start_str = Start, end_str = End) |>
  mutate(
    start_utc = dmy_hm(start_str, tz = "UTC"),
    end_utc   = dmy_hm(end_str,   tz = "UTC")
  ) |>
  select(label, start_utc, end_utc) |>
  arrange(start_utc)

## My Wind data --------------
wind_23 <- read_csv("../metadata/2023/10min_Oct22_Aug23_weather.csv")
wind_24 <- read_csv("../metadata/2024/10_min.csv")

# Standardise DateTime and select rows of interest 
wind_23 <- wind_23 %>%
  mutate(hour_utc = dmy_hm(`Timestamp (UTC+0)`, tz = "UTC")) %>%
  select(hour_utc, ws = Mean_WS)
wind_24 <- wind_24 %>%
  mutate(hour_utc = as.POSIXct(`Timestamp (UTC+0)`, tz = "UTC")) %>% 
  select(hour_utc, ws = Mean_WS)

# Bind with labels
your_all <- bind_rows(
  wind_23 %>% mutate(source = "Site", year = 2023),
  wind_24 %>% mutate(source = "Site", year = 2024)
)

# Aggregate to hourly
site_hourly <- your_all %>%
  mutate(hour = floor_date(hour_utc, "hour")) %>%
  group_by(hour, source, year) %>%
  summarise(
    ws = mean(ws, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(hour_utc = hour)


# Public wind data -----------
# define your site coordinates (Weather station = 52°37’52.29”, 1°10’23.57”)
site_lat <- 52.37   
site_lon <- 1.17  

# list nearby ISD stations
meta <- getMeta(lat = site_lat, lon = site_lon)

# Order by distance and show first few
meta <- meta %>%
  arrange(distance)
head(meta, 5)

# Pick TIBENHAM AIRFIELD which is 10km away (top of meta)
station_code <- paste0(meta$usaf[1], "-", meta$wban[1])

# Download hourly data for these years and this weather station 
years_needed <- c(2022, 2023, 2024)
met_data <- importNOAA(code = station_code, year = years_needed)

# Tidy: keep just time + wind speed ---
met_ws <- met_data %>%
  transmute(
    hour_utc = floor_date(as_datetime(date, tz = "UTC"), "hour"),
    ws   = ws) %>%       # wind speed in m/s
    mutate(
        source = "Public",
        year   = year(hour_utc)
      )


# Combine both datasets -----------
all_ws <- bind_rows(your_all, met_ws)
all_hourly_ws <- bind_rows(site_hourly, met_ws)

# function to filter all_ws by a given window
get_window_data <- function(data, lbl) {
  w <- windows %>% filter(label == lbl)
  data %>%
    filter(hour_utc >= w$start_utc, hour_utc <= w$end_utc) %>%
    mutate(window = lbl)
}

# Select the needed data for two collections 
sel_data <- bind_rows(
  get_window_data(all_ws, "all_2023"),
  get_window_data(all_ws, "all_2024")
)

sel_hrly_data <- bind_rows(
  get_window_data(all_hourly_ws, "all_2023"),
  get_window_data(all_hourly_ws, "all_2024")
)


# filter out extreme site values
all_ws_filt <- sel_data %>%
  filter(!(source == "Site" & ws > 100))

sel_hrly_filt <- sel_hrly_data %>%
  filter(!(source == "Site" & ws > 100))

# Plotting ----------
# Set plotting theme
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

my_cols <- brewer.pal(3, "Dark2")[1:2]

# Some clear too high spikes in 2023 site collected wind speed
p1 <- ggplot(sel_hrly_data, aes(x = hour_utc, y = ws, colour = source)) +
  geom_line(alpha = 0.5) +
  scale_colour_manual(values = my_cols) +
  facet_wrap(~year, scales = "free_x", ncol = 1) +
  labs(
    title = "Unfiltered Site Data",
    x = "Month", y = "Wind speed (m/s)", colour = "Source"
  )  +
  scale_x_datetime(
    date_labels = "%b\n%y",   # e.g. Jan 23, Feb 23
    date_breaks = "1 month"  # one tick per month
  ) + custom_theme


# Plot with these values removed
p2 <- ggplot(sel_hrly_filt, aes(x = hour_utc, y = ws, colour = source)) +
  geom_line(alpha = 0.5) +
  scale_colour_manual(values = my_cols) +
  facet_wrap(~year, scales = "free_x", ncol = 1) +
  labs(
    title = "Filtered Site Data (<100 m/s)",
    x = "Month", y = "Wind speed (m/s)", colour = "Source"
  )  +
  scale_x_datetime(
    date_labels = "%b\n%y",   # e.g. Jan 23, Feb 23
    date_breaks = "1 month"  # one tick per month
  ) + custom_theme

# Combine side by side
combined <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined
ggsave("../Graphs/wind_comparison.PDF", combined, height = 6, width = 10, dpi = 300)



# Pearson correlation  --------
# Just looking at the dates I reported
# correlation before filtering
cor_before <- sel_hrly_data %>%
  tidyr::pivot_wider(names_from = source, values_from = ws) %>%
  summarise(corr = cor(Site, Public, use = "complete.obs"))


# correlation after filtering
cor_after <- sel_hrly_filt %>%
  tidyr::pivot_wider(names_from = source, values_from = ws) %>%
  summarise(corr = cor(Site, Public, use = "complete.obs"))

# Look at public data when incorrect ---------
# find timestamps removed because ws > 100
removed_times <- sel_data %>%
  filter(source == "Site", ws > 100) %>%
  pull(hour_utc)

# look up public data for those times
removed_public <- sel_data %>%
  filter(source == "Public", hour_utc %in% removed_times)

removed_public

