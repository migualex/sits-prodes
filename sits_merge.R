#
# Merge test between MPC Sentinel-1 RTC and BDC Sentinel-2
#


#
# 1. Libraries
#


library(sits)
library(sf)


#
# 2. define the output directory
#


output_dir <- "data/merge_tile_bdc"


#
# 3. Creating BDC Sentinel-2 data cube
#


cube_s2 <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = c("B02", "B8A", "B11", "CLOUD"),
  tiles = c("012014", "012015", "013014", "013015"),
  start_date = "2024-07-01",
  end_date = "2025-07-31"
)

# Checking the timeline of the BDC cube
sits_timeline(cube_s2)


#
# 4. Creating MPC Sentinel-1 RTC data cube
#


# Part 1:


cube_s1_rtc_part_1 <-  sits_cube(
  source = "MPC",
  collection = "SENTINEL-1-RTC",
  bands = c("VV", "VH"),
  orbit = "descending",
  tiles = c("20MNS", "20MMS", "20LNR", "20LNQ", "20LMR", "20LMQ", "20LLR", "20LLQ", "20MLS"),
  start_date = "2024-07-01",
  end_date = "2025-03-01"
)


# Part 2:


cube_s1_rtc_part_2 <- sits_cube(
  source = "MPC",
  collection = "SENTINEL-1-RTC",
  bands = c("VV", "VH"),
  orbit = "descending",
  tiles = c("20MNS", "20MMS", "20LNR", "20LNQ", "20LMR", "20LMQ", "20LLR", "20LLQ", "20MLS"),
  start_date = "2025-03-15",
  end_date = "2025-07-31"
)


#
# 5. Merging Sentinel-1 RTC cubes parts 1 and 2
#


cube_s1_rtc <- sits_merge(cube_s1_rtc_part_1, cube_s1_rtc_part_2)


# Checking the timeline of the MPC cube
sits_timeline(cube_s1_rtc)


#
# 6. Creating an output directory for the regularized data
#


# Defining the output directory


output_dir_s1_s2 <- fs::path(output_dir) / "s1_s2_reg"


# Creating output dir


fs::dir_create(output_dir_s1_s2)


#
# 7. Merging Sentinel-1 RTC and Sentinel-2 cubes
#


cube_s1_s2 <- sits_merge(cube_s1_rtc, cube_s2)


#
# 8. Regularizing
#


cube_s1_s2_reg <-  sits_regularize(
  cube = cube_s1_s2,
  period = "P16D",
  timeline = sits_timeline(cube_s2),
  res = 40,
  tiles = c("012014", "012015", "013014", "013015"),
  grid_system = "BDC_SM_V2",
  memsize = 12,
  multicores =1,
  output_dir = output_dir_s1_s2
)


#
# 9. Loading samples
#


samples <- st_read("data/raw/samples_prodes_test.gpkg")


#
# 10. Creating an output directory for time series
#


# Defining the output directory
output_dir_ts <- fs::path(output_dir, "ts")

# Creating output dir
fs::dir_create(output_dir_ts)


#
# 11. Extracting samples from the regularized cube
#


ts <- sits_get_data(cube = cube_s1_s2_reg,
                    samples = samples,
                    multicores = 17
)

# Saving samples in .rds
saveRDS(ts, fs::path(output_dir_ts, "ts_s1_s2_reg_v1.rds"))


# 12. Visualizing the location of the samples


sits_view(ts)


# 13. Visualizing time series:


plot(ts)


# 14. Visualizing time series via GAM


# Sentinel-1-RTC:

s1 <- sits_select(ts, bands = c("VV", "VH"))

plot(sits_patterns(s1))

# Sentinel-2:

s2 <- sits_select(ts, bands = c("B02", "B8A", "B11"))

plot(sits_patterns(s2))
