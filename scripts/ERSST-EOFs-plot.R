## ---------- PLOT FIGURES as in CLUETT et al., 2025 ------------- ##
# Pan-Basin Warming Now Overshadows Robust Pacific Decadal Oscillation #

# Import libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(lubridate)
library(metR)
library(ncdf4)
library(purrr)
library(patchwork)
library(raster)
library(stringr)
library(tictoc)
library(tidyr)
library(zoo)

setwd('output')

# Plotting utilities ------------------------------------------------------

colfunc <- colorRampPalette(c("black", "gray80"))

# Read in Natural Earth coastline data (public domain)
world <-
  maps::map("world",
            fill = TRUE,
            col = "transparent",
            plot = FALSE)

# Center on North Pacific
worldSpP <-
  maptools::map2SpatialPolygons(world, world$names, sp::CRS("+proj=longlat +ellps=WGS84"))
worldSpP <- worldSpP[-grep("Antarctica", row.names(worldSpP)),]
worldSpP <- worldSpP[-grep("Ghana", row.names(worldSpP)),]
worldSpP <-
  worldSpP[-grep("UK:Great Britain", row.names(worldSpP)),]
worldSpPnr <- maptools::nowrapRecenter(worldSpP)
world_map <- fortify(worldSpPnr)

plot_eof_pattern <- function(file, pc = 1, title.text = "", rev = FALSE) {
  
  # Open NetCDF file
  open <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(open))  # ensure file closes if function exits early
  
  # Extract SST variable
  returnVar <- "sst"
  arrayVar  <- ncvar_get(open, returnVar)
  
  # Replace fill values and zeros with NA
  fillVal <- ncatt_get(open, returnVar, "_FillValue")$value
  arrayVar[arrayVar == fillVal | arrayVar == 0] <- NA
  
  sstArray <- arrayVar
  rm(arrayVar)
  
  # Extract dimensions
  returnDim <- c("time", "lat", "lon")
  dimList <- lapply(returnDim, function(dim) {
    dimVals <- ncvar_get(open, dim)
    units   <- ncatt_get(open, dim, "units")$value
    list(dimVals, units)
  })
  names(dimList) <- returnDim
  
  # Lat/Lon
  lat <- dimList$lat[[1]]
  lon <- dimList$lon[[1]]
  
  # Extract requested EOF pattern
  eof.pattern <- sstArray[, , pc]
  
  # Reshape for ggplot
  eof.melt <- reshape2::melt(eof.pattern, varnames = c("lon", "lat"))
  eof.melt$lat.val <- lat[eof.melt$lat]
  eof.melt$lon.val <- lon[eof.melt$lon]
  
  # Shift longitudes to [0, 360]
  eof.melt$lon.val[eof.melt$lon.val < 0] <- eof.melt$lon.val[eof.melt$lon.val < 0] + 360
  
  # Flip sign if requested
  if (rev) eof.melt$value <- -eof.melt$value
  
  # Plot
  eof.plot <- ggplot(eof.melt, aes(x = lon.val, y = lat.val)) +
    metR::geom_contour_fill(aes(z = value), binwidth = 0.2) +
    scale_fill_divergent(limits = c(-2.5, 2.5), oob = scales::squish) +
    geom_map(data = world_map, map = world_map,
             aes(x = long, y = lat, map_id = id),
             color = "black", fill = "black", size = 0.25) +
    theme_test() +
    xlab("Lon (°E)") + ylab("Lat (°N)") +
    scale_x_continuous(limits = c(0, 360), expand = c(0, 0), breaks = seq(0, 360, by = 45)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0, 0), breaks = seq(-90, 90, by = 45)) +
    coord_equal() +
    # Focus on North Pacific
    scale_x_continuous(limits = c(100, 260), expand = c(0, 0)) +
    scale_y_continuous(limits = c(20, 75), expand = c(0, 0)) +
    ggtitle(title.text) +
    labs(fill = "EOF Loading")
  
  return(eof.plot)
}


# --------------------------------------------------------------------
# Figure 1: Sliding window EOF plots (total & gmsstr, PC1 & PC2)
# --------------------------------------------------------------------

# Metadata: years and whether EOF sign should be flipped
windows <- tibble::tibble(
  end_year   = 2007:2023,
  start_year = end_year - 29,
  rev        = c(FALSE, FALSE, FALSE, FALSE, FALSE,
                 FALSE, FALSE, FALSE, TRUE, TRUE,
                 FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
) %>% 
  dplyr::filter((end_year %% 2) == 1) 

# Function to build windows with filenames/titles
make_windows <- function(base, prefix) {
  base %>%
    dplyr::mutate(
      file  = sprintf("eofs-%s-%d_%d.nc", prefix, start_year, end_year),
      title = sprintf("%d-%d", start_year, end_year)
    )
}

# Generate both total and gmsstr
windows_total  <- make_windows(windows, "total")
windows_gmsstr <- make_windows(windows, "gmsstr")

# Helper to generate plots
make_plots <- function(windows_tbl, pc, suffix) {
  purrr::pmap(
    windows_tbl,
    function(end_year, start_year, rev, file, title) {
      plot_eof_pattern(file = file, pc = pc, title.text = title, rev = rev)
    }
  ) %>%
    rlang::set_names(sprintf("eof%d.%d.%s", pc, windows_tbl$end_year, suffix))
}

# Generate all plots
plots <- list(
  pc1_total  = make_plots(windows_total,  1, "total"),
  pc2_total  = make_plots(windows_total,  2, "total"),
  pc1_gmsstr = make_plots(windows_gmsstr, 1, "gmsstr"),
  pc2_gmsstr = make_plots(windows_gmsstr, 2, "gmsstr")
)

# Assemble figure columns
fig1_col1 <- patchwork::wrap_plots(plots$pc1_total,  ncol = 1, guides = "collect") & 
  theme(legend.position = "bottom")

fig1_col2 <- patchwork::wrap_plots(plots$pc2_total,  ncol = 1, guides = "collect") & 
  theme(legend.position = "bottom")

fig1_col3 <- patchwork::wrap_plots(plots$pc1_gmsstr, ncol = 1, guides = "collect") & 
  theme(legend.position = "bottom")

fig1_col4 <- patchwork::wrap_plots(plots$pc2_gmsstr, ncol = 1, guides = "collect") & 
  theme(legend.position = "bottom")

# Combine four columns into a single figure
fig1 <- patchwork::wrap_plots(
  fig1_col1, fig1_col2, fig1_col3, fig1_col4,
  ncol = 4, guides = "collect"
) &
  theme(legend.position = "bottom")

# Save fig 1
ggplot2::ggsave(
  filename = "../figs/Figure1.pdf",   # output file name
  plot = fig1,                # patchwork figure object
  device = "pdf",             # ensure PDF output
  width = 16,                 # width in inches
  height = 12,                # height in inches
  units = "in",
  dpi = 300                   
)


# --------------------------------------------------------------------
# Figure 2: Sliding window % variance (total & gmsstr, PC1 & PC2)
# --------------------------------------------------------------------

read_variance <- function(files, pc = 1) {
  lst <- lapply(files, function(f) {
    read.csv(f, header = FALSE)[pc, , drop = TRUE]  # drop=TRUE returns a vector
  })
  
  df <- as.data.frame(do.call(rbind, lst))
  colnames(df) <- "value"
  df$filename <- files
  df
}

gmsstr_files <- list.files(pattern = "var-gmsstr.+\\.csv")
total_files  <- list.files(pattern = "var-total.+\\.csv")

gmsstr_var1 <- read_variance(gmsstr_files, 1)
gmsstr_var2 <- read_variance(gmsstr_files, 2)
total_var1  <- read_variance(total_files, 1)
total_var2  <- read_variance(total_files, 2)

plot_var <- function(df, x_vals = c(1994:2023)) {
  ggplot(df, aes(x = x_vals, y = value * 100)) +
    geom_col() +
    theme_minimal() +
    ylim(0, 50) +
    ylab("% Var Explained") +
    xlab("End of 30-yr Window")
}

# Generate plots
gmsstr_var1_plot <- plot_var(gmsstr_var1)
gmsstr_var2_plot <- plot_var(gmsstr_var2)
total_var1_plot  <- plot_var(total_var1)
total_var2_plot  <- plot_var(total_var2)

fig2 <- (total_var1_plot + total_var2_plot) / 
  (gmsstr_var1_plot + gmsstr_var2_plot)

# Save fig 2
ggplot2::ggsave(
  filename = "../figs/Figure2.pdf",   # output file name
  plot = fig2,                # patchwork figure object
  device = "pdf",             # ensure PDF output
  width = 6,                 # width in inches
  height = 4,                # height in inches
  units = "in",
  dpi = 300                   
)


# --------------------------------------------------------------------
# Figure 3: Global and regional time series
# --------------------------------------------------------------------

read_sst_nc <- function(file, varname = "sst") {
  nc <- nc_open(file)
  arrayvar <- ncvar_get(nc, varname)
  fillVal <- ncatt_get(nc, varname, "_FillValue")$value
  arrayvar[arrayvar == fillVal | arrayvar == 0] <- NA
  
  # Get dimensions and units
  dims <- lapply(c("time", "lat", "lon"), function(dn) {
    dim_vals <- ncvar_get(nc, dn)
    dim_units <- ncatt_get(nc, dn, "units")$value
    list(values = dim_vals, units = dim_units)
  })
  names(dims) <- c("time", "lat", "lon")
  nc_close(nc)
  list(data = arrayvar, dims = dims)
}

sst_total  <- read_sst_nc("total-anoms.nc")
sst_gmsstr <- read_sst_nc("gmsstr-anoms.nc")

# Access arrays and dimensions
sstArray       <- sst_total$data
sstArray.gmsstr <- sst_gmsstr$data
lat            <- sst_total$dims$lat$values
lon            <- sst_total$dims$lon$values
timesteps      <- sst_total$dims$time$values %>% as.Date(., origin = '1800-01-01')

read_pc <- function(file) {
  read.csv(file, header = FALSE)$V1 %>% as.vector()
}

gmsstr.pc1   <- read_pc("gmsstr-1994-2023-pc1.csv")
gmsstr.pc2   <- read_pc("gmsstr-1994-2023-pc2.csv")
total.pc1  <- read_pc("total-1994-2023-pc1.csv")
total.pc2  <- read_pc("total-1994-2023-pc2.csv")

official.pdo.dat <- read.delim(
  'https://psl.noaa.gov/pdo/data/pdo.timeseries.sstens.csv',
  header = TRUE,
  sep = ""
)

colnames(official.pdo.dat) <- c("year", as.character(1:12))

official.pdo.melt <- reshape2::melt(official.pdo.dat, id.vars = "year") %>%
  mutate(date = as.Date(paste0(variable, "-01-", year), format = "%m-%d-%Y"))

gmsst <- read.csv("ersst_globalmean.csv", header = FALSE)
gmAnom <- gmsst$V1 - mean(gmsst$V1)

# Calculate decomposed SST anoms

sst.regression <- function(pc = gmr.pc1) {
  
  int <-
    apply(
      sstArray,
      MARGIN = c(1, 2),
      FUN = function(x) {
        if (is.na(x[1])) {
          return(NA)
        } else{
          mod <- lm(x ~ pc)
          coef(mod) %>% .[1]
        }
      }
    )
  
  slope <-
    apply(
      sstArray,
      MARGIN = c(1, 2),
      FUN = function(x) {
        if (is.na(x[1])) {
          return(NA)
        } else{
          mod <- lm(x ~ pc)
          coef(mod) %>% .[2]
        }
      }
    )
  
  transfer.fun = function(x){
    x * slope + int}
  
  anom <- sapply(pc, FUN = transfer.fun, simplify = 'array')
  return(anom)
}

gmsstr.pc1.anoms <- gmsstr.pc1 %>% sst.regression()
total.pc1.anoms <- total.pc1 %>% sst.regression()
total.pc2.anoms <- total.pc2 %>% sst.regression()

# Calculate region means

# Read in masks
file = '../mapdata/LME_mask_ERSST.nc'

open <- file %>% ncdf4::nc_open()
returnVar = 'mask'
arrayVar <- ncvar_get(open, returnVar)
fillVal <- ncatt_get(open, returnVar, "_FillValue")
arrayVar[arrayVar == fillVal$value | arrayVar == 0] <-
  NA # replace fill val with NAs
maskArray <- arrayVar
remove(arrayVar)
nc_close(open)

region.ids <- c(1, 2, 3, 4, 10, 11, 13, 49, 51, 53, 54)
maskArray[which(!maskArray %in% region.ids)] <- 0

gw.file = '../mapdata/gridweights.nc'

open <- gw.file %>% ncdf4::nc_open()
returnVar = 'cell_weights'
arrayVar <- ncvar_get(open, returnVar)
fillVal <- ncatt_get(open, returnVar, "_FillValue")
arrayVar[arrayVar == fillVal$value | arrayVar == 0] <-
  NA # replace fill val with NAs
gw <- arrayVar
remove(arrayVar)
nc_close(open)

col.names = c('dates',
              'EBS',
              'GOA',
              'CCS',
              'GOC',
              'IPH',
              'PCAC',
              'HC',
              'KC',
              'OC',
              'WBS',
              'CS')

region_means <- function(anoms = gmsstr.pc1.anoms){
  
  region_ts_res <-
    data.frame(matrix(ncol = 12, nrow = length(timesteps)))
  
  colnames(region_ts_res) <- col.names
  
  for (j in 1:length(region.ids)) {
    idx <- region.ids[j]
    maskArray1 <- maskArray
    maskArray1[maskArray == idx] <- 1
    maskArray1[maskArray != idx] <- 0
    region.wts <- maskArray1 * gw
    
    tic()
    region.anom.res <-  apply(anoms,
                                     region.wts,
                                     MARGIN = c(3),
                                     FUN = '*')  %>% apply(.,
                                                           MARGIN = 2,
                                                           FUN = mean,
                                                           na.rm = T) / (region.wts %>% mean(na.rm = T))
    toc()
    
    region_ts_res[, j + 1] <- region.anom.res
    
  }
  
  return(region_ts_res)
}
  
region_ts_gmsstr1 <- region_means(gmsstr.pc1.anoms)
region_ts_total1 <- region_means(total.pc1.anoms)
region_ts_total2 <- region_means(total.pc2.anoms)

region_ts_anoms <- region_means(sstArray)
region_ts_anoms_gmsstr <- region_means(sstArray.gmsstr)

retrieve_pc_ts <-
  function(file,
           start.date = '1-01-1854',
           end.date = '07-24-2024') {
    ts <- read.csv(file, head = F)
    time <-
      seq(
        from = start.date %>% as.Date(., format = '%m-%d-%Y'),
        end.date %>% as.Date(., format = '%m-%d-%Y'),
        by = 'month'
      )
    data.frame(pc = ts %>% unlist(),
               time = time)
  }

pc1.total.2023 <- retrieve_pc_ts("total-1994-2023-pc1.csv")
pc2.total.2023 <- retrieve_pc_ts("total-1994-2023-pc2.csv")
pc1.gmsstr.2023 <- retrieve_pc_ts("gmsstr-1994-2023-pc1.csv")

pc1.gmsstr.2023$pc <- pc1.gmsstr.2023$pc * -1
pc2.total.2023$pc <- pc2.total.2023$pc

ggplot() + geom_point(aes(
  x = gmAnom,
  y = pc1.total.2023$pc / sd(pc1.total.2023$pc),
  col = pc1.total.2023$time
)) +
  scale_color_distiller(
    name = 'Date',
    type = "seq",
    direction = 1,
    palette = "Greys",
    trans = 'date'
  ) + xlab('GMSST (°C)') + ylab('PBP') + theme_minimal()

gmT.plot <-  ggplot() +
  geom_line(aes(x = pc1.total.2023$time,
                y = (gmsst$V1 - mean(gmsst$V1))),
            lwd = 0.8,
            alpha = 0.8) +
  scale_fill_manual(values = c('darkred', 'mediumblue')) + theme_minimal() + xlab('Time') + ylab('Global Mean SST (°C)') +
  theme(legend.position = 'none') + ggtitle('Global Mean SST')

pc1.gmsstr.plot <- ggplot() +
  ggh4x::stat_difference(aes(
    x = pc1.total.2023$time,
    ymin = 0,
    ymax = (pc1.gmsstr.2023$pc / sd(pc1.gmsstr.2023$pc))
  ),
  alpha = 0.7) +
  scale_fill_manual(values = c('darkred', 'mediumblue')) + theme_minimal() + xlab('Time') + ylab('PC1') +
  theme(legend.position = 'none') + ggtitle('GMSSTR PC1')

pc1.total.plot <- ggplot() +
  ggh4x::stat_difference(aes(
    x = pc1.total.2023$time,
    ymin = 0,
    ymax = (pc1.total.2023$pc / sd(pc1.total.2023$pc))
  ),
  alpha = 0.7) +
  scale_fill_manual(values = c('darkred', 'mediumblue')) + theme_minimal() + xlab('Time') + ylab('PC1') +
  theme(legend.position = 'none') + ggtitle('Total PC1')

pc2.total.plot <- ggplot() +
  ggh4x::stat_difference(aes(
    x = pc2.total.2023$time,
    ymin = 0,
    ymax = (pc2.total.2023$pc / sd(pc2.total.2023$pc))
  ),
  alpha = 0.7) +
  scale_fill_manual(values = c('darkred', 'mediumblue')) + theme_minimal() + xlab('Time') + ylab('PC1') +
  theme(legend.position = 'none') + ggtitle('Total PC2')

# Fig. 3 (left panel)
fig3ac <- gmT.plot + pc1.total.plot + pc2.total.plot  + pc1.gmsstr.plot + plot_layout(ncol = 1)

# Save fig 3A-C
ggplot2::ggsave(
  filename = "../figs/Figure3A-C.pdf",   # output file name
  plot = fig3ac,                # patchwork figure object
  device = "pdf",             # ensure PDF output
  width = 4,                 # width in inches
  height = 6,                # height in inches
  units = "in",
  dpi = 300                   
)

t.ax <- pc1.total.2023$time

plot_region_zoom <- function(region_name,
                             time_total = pc1.total.2023$time,
                             time_gmsstr = pc1.total.2023$time,
                             roll_width = 3,
                             alpha_total = 0.4,
                             alpha_gmsstr = 0.7,
                             ylims = c(-2.5, 5),
                             ybreaks = c(-2, 0, 2)) {
  
  ggplot() +
    ggh4x::stat_difference(
      aes(
        x = time_total,
        ymin = 0,
        ymax = zoo::rollmean(region_ts_total1[[region_name]], roll_width, na.pad = TRUE)
      ),
      alpha = alpha_total
    ) +
    ggh4x::stat_difference(
      aes(
        x = time_gmsstr,
        ymin = 0,
        ymax = zoo::rollmean(region_ts_gmsstr1[[region_name]], roll_width, na.pad = TRUE)
      ),
      alpha = alpha_gmsstr
    ) +
    scale_fill_manual(values = c('darkred', 'mediumblue')) +
    theme_minimal() +
    theme(legend.position = 'none') +
    # Line for anomalies
    geom_line(
      aes(
        x = time_gmsstr,
        y = zoo::rollmean(region_ts_anoms[[region_name]], roll_width, na.pad = TRUE)
      )
    ) +
    # Dashed line: total + gmsstr
    geom_line(
      aes(
        x = time_gmsstr,
        y = zoo::rollmean(region_ts_total1[[region_name]] + region_ts_gmsstr1[[region_name]], roll_width, na.pad = TRUE)
      ),
      lty = 'dashed'
    ) +
    xlim(as.Date('2000-01-01'), as.Date('2024-02-01')) +
    scale_y_continuous(limits = ylims, breaks = ybreaks) +
    xlab('Time') +
    ylab('SSTa (°C)')
}

# For multiple regions:
regions <- c("GOA", "CCS", "OC", "KC")
plots_zoom <- lapply(regions, plot_region_zoom)

# Combine into a vertical panel:
fig3dg <- wrap_plots(plots_zoom, ncol = 1)

# Save fig 3A-C
ggplot2::ggsave(
  filename = "../figs/Figure3D-G.pdf",   # output file name
  plot = fig3dg,                # patchwork figure object
  device = "pdf",             # ensure PDF output
  width = 3,                 # width in inches
  height = 4,                # height in inches
  units = "in",
  dpi = 300                   
)

# --------------------------------------------------------------------
# Figure 4: PDO-PBP biplots
# --------------------------------------------------------------------

# --- Create annual-mean data frames for each region ---
make_region_df <- function(region) {
  df <- data.frame(
    time = t.ax,
    PDO = region_ts_gmsstr1[[region]],    
    GM  = region_ts_total1[[region]],     
    obs = region_ts_anoms[[region]]
  )
  df$year <- year(df$time)
  
  list(
    PDO = aggregate(PDO ~ year, data = df, mean),
    GM  = aggregate(GM ~ year, data = df, mean),
    obs = aggregate(obs ~ year, data = df, mean)
  )
}

CCS <- make_region_df("CCS")
GOA <- make_region_df("GOA")
KC  <- make_region_df("KC")
OC  <- make_region_df("OC")

t.ax.year <- unique(year(t.ax))

# --- Triangle background data for XY plots ---
dt.triangle <- data.table(
  group = c(1, 1, 1),
  polygon.x = c(2, 4, 4),
  polygon.y = c(1, 1, 3)
)

# --- Function to create XY plot for each region ---
make_xy_plot <- function(region_data, xlim_vals = c(-2.5, 2.5), ylim_vals = c(-2.5, 2.5)) {
  ggplot() +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    
    # Triangle background (4 quadrants)
    geom_polygon(aes(x = c(2.5, 0, 2.5), y = c(-2.5, 0, 2.5)), fill = '#CF6F86', alpha = 0.7) +
    geom_polygon(aes(y = c(2.5, 0, 2.5), x = c(-2.5, 0, 2.5)), fill = '#8B172F', alpha = 0.7) +
    geom_polygon(aes(y = c(-2.5, 0, -2.5), x = c(-2.5, 0, 2.5)), fill = '#415AA9', alpha = 0.7) +
    geom_polygon(aes(x = c(-2.5, 0, -2.5), y = c(-2.5, 0, 2.5)), fill = '#6F9AAD', alpha = 0.7) +
    
    # Path for GM vs PDO
    geom_path(aes(x = region_data$GM$GM, y = region_data$PDO$PDO), alpha = 0.4, col = 'white') +
    
    # Points colored by observed anomalies
    geom_point(aes(x = region_data$GM$GM, y = region_data$PDO$PDO, fill = region_data$obs$obs),
               col = 'white', alpha = 0.85, size = 4, pch = 21) +
    scale_fill_divergent(limits = c(-3.5, 3.5)) +
    
    # Label years >= 2014
    geom_text(aes(
      x = region_data$GM$GM[which(t.ax.year >= 2014)],
      y = region_data$PDO$PDO[which(t.ax.year >= 2014)],
      label = t.ax.year[which(t.ax.year >= 2014)]
    ), hjust = 0, nudge_x = 0.1) +
    
    xlim(xlim_vals) + ylim(ylim_vals) +
    theme_classic() + coord_equal() +
    xlab('PBP Anom (°C)') + ylab('PDO Anom (°C)') +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    guides(size = "none")
}

# --- Generate XY plots ---
CCS.xy <- make_xy_plot(CCS)
GOA.xy <- make_xy_plot(GOA)
KC.xy  <- make_xy_plot(KC)
OC.xy  <- make_xy_plot(OC)

# --- Combine into figure ---
fig4 <- OC.xy + GOA.xy + KC.xy + CCS.xy +
  plot_layout(nrow = 2, guides = 'collect') &
  theme(legend.position = 'bottom')

# Save fig 4
ggplot2::ggsave(
  filename = "../figs/Figure4.pdf",   # output file name
  plot = fig4,                # patchwork figure object
  device = "pdf",             # ensure PDF output
  width = 9,                 # width in inches
  height = 9,                # height in inches
  units = "in",
  dpi = 300                   
)

# --------------------------------------------------------------------
# Figure 4: PDO-PBP biplots
# --------------------------------------------------------------------

# Example data
df <- data.frame(
  PBP = pc1.total.2023$pc/sd(pc1.total.2023$pc),
  time = pc1.total.2023$time
)

# Add year and month columns
df <- df %>%
  mutate(
    year = year(time),
    month = month(time)
  ) %>%
  dplyr::select(-time)  # drop original date if not needed

# Pivot to wide format: rows = year, columns = month
df_wide <- df %>%
  pivot_wider(
    names_from = month,
    values_from = PBP
  ) %>%
  rename_with(
    ~ month.abb[as.numeric(.)], 
    .cols = -year  # only apply to non-year columns
  )%>% mutate(across(-year, ~ round(.x, 2)))

# View
df_wide

# Write to CSV
write.csv(df_wide, "../PBP_Index_Values.csv", row.names = FALSE)
