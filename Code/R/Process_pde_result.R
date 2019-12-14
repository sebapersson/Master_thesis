library(tidyverse)
library(stringr)

# General plotting parameters 
my_theme <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                               plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=13)) 

# Colour-blind friendly palette
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Function that will read a t_end_data file and convert it to a grid-based format and plot the 
# concentration in form of a heat map. The created figure will finally be written to file 
# (stored in the result folder)
# Args:
#   path_data, the path to the data file
#   path_save, the path to where the file will be saved
#   plot_result, if true the result is plotted 
#   n_cores, the number of cores to use when doing the grid matching
#   is_circle, true if the geometry is a circle 
plot_t_end_data <- function(path_data, path_save, limit_grid, plot_result=T, n_cores=3, is_circle=F, h=0.05, seed=123)
{
  set.seed(seed)
  
  # Read the data and aggregate 
  data <- read_csv(path_data, col_types = cols()) %>% 
    select(-X1) %>%
    group_by(x, y) %>%
    summarise(u1_med = mean(u1), 
              u1_std = sd(u1), 
              u2_med = mean(u2), 
              u2_std = sd(u2))
  
  # Make a grid 
  x1s <- seq(-limit_grid, limit_grid, by = h)
  x2s <- seq(-limit_grid, limit_grid, by = h)
  X_pred <- expand.grid(x1s, x2s)
  # Find the closest point in the grid to the points in the data
  u_grid <- do.call(rbind, parallel::mclapply(1:dim(X_pred)[1], function(i){
    best_index <- which.min((data$x - X_pred$Var1[i])^2 + (data$y - X_pred$Var2[i])^2)[1]
    return(data$u1_med[best_index])}, mc.cores = n_cores))
  
  data_to_plot <- tibble(x = X_pred$Var1, 
                         y = X_pred$Var2, 
                         u1 = u_grid)
  
  
  # If circle set outside region to zero 
  if(is_circle){
    data_to_plot[data_to_plot$x^2 + data_to_plot$y^2 > 2.5^2, 3] <- 0
  }
  
  # Viridis will be used for plotting the result 
  p1 <- ggplot(data_to_plot, aes(x, y, fill = u1)) + 
    geom_raster(interpolate = T) + 
    scale_fill_viridis_c() + 
    labs(x = "x", y = "y") + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  # Write result to file 
  p1 
  ggsave(path_save, height = 6, width = 9)
  
  if(plot_result == T) print(p1)
}

# Function that will check if dir-exists, if a directory doesn't exist a 
# new directory will be created
# Args:
#   path_dir, the path to the directory 
check_if_dir_exists <- function(path_dir) if(!dir.exists(path_dir)) dir.create(path_dir)


# Function that for a certain geometry will plot the maximum concentration vs time, note that the 
# function works under the assumption of 0, 5 and twenty holes in the file-names. 
# Args:
#   dir_fil, the path to the directory where the files can be found 
#   geometry, the geometry considered 
#   path_save, the path to where the file should be saved
#   plot_data, if the created plot should be displayed or not, true by default
plot_max_conc_data <- function(dir_files, geometry, path_save, plot_data=T)
{
  # The files are saved as max_conc
  file_name <- "max_conc.csv"
  
  # Read the data
  path_data <- str_c(dir_files, "Zero_holes_files/", file_name)
  data_0_holes_rec <- read_csv(path_data, col_types = cols()) %>%
    mutate(geometry = geometry) %>%
    mutate(n_holes = 0)
  path_data <- str_c(dir_files, "Five_holes_files/", file_name)
  data_5_holes_rec <- read_csv(path_data, col_types = cols()) %>%
    mutate(geometry = geometry) %>%
    mutate(n_holes = 5)
  path_data <- str_c(dir_files, "Twenty_holes_files/", file_name)
  data_20_holes_rec <- read_csv(path_data, col_types = cols()) %>%
    mutate(geometry = geometry) %>%
    mutate(n_holes = 20)
  
  # Aggregate the result   
  data_tot_rec <- data_0_holes_rec %>%
    bind_rows(data_5_holes_rec) %>%
    bind_rows(data_20_holes_rec) %>% 
    mutate(geometry = as.factor(geometry)) %>%
    mutate(n_holes = as.factor(n_holes)) %>%
    select(-X1) %>%
    group_by(time, n_holes) %>%
    summarise(med = mean(Max_conc), 
              quant_low = quantile(Max_conc, 0.05),
              quant_high = quantile(Max_conc, 0.95)) 
  
  p1 <- ggplot(data_tot_rec, aes(time, med, color = n_holes, fill = n_holes)) + 
    geom_line(size = 1.2) + 
    geom_ribbon(aes(ymin = quant_low, ymax = quant_high), color = NA, alpha = 0.2) + 
    labs(x = "Time", y = "Maximum concentration") + 
    scale_color_manual(values = cbPalette[-1]) + 
    my_theme
  
  p1
  ggsave(path_save, height = 6, width = 9)
  
  if(plot_data == T) print(p1)
}

# Function that will plot 10 randomly selected end-time profiles, this with the goal to charactersise 
# what kind of patters are created when holes are inserted into the grid. 
# Args:
#   path_data, the path to the end-concentration file 
#   path_save, path to where the result will be saved 
#   limit_grid, the limits for the grid (differs between circles and rectangles)
#   h, the resolution of the grid (0.05 by default)
#   n_cores, the number of cores when doing grid search (3 by default)
#   is_circle, true of the underlaying geometry is a circle 
# Returns:
#   void 
plot_several_end_times <- function(path_data, path_save, limit_grid, h=0.05, n_cores=3, is_circle=F, seed=123)
{
  
  set.seed(seed)
  # Write function to plot a random sample of t-end results (like 10 of them)
  data <- read_csv(path_data, col_types = cols()) %>%
    select(-X1) 
  
  data <- data %>%
    mutate(id = as.factor(rep(1:20, each = length(data$x) / 20))) %>%
    select(id, everything())
  
  # Using 10 to (else the plot will be to messy)
  n_random_figs <- 10
  id_to_plot <- sample(1:20, n_random_figs, replace = F)
  
  # Create the grid 
  x1s <- seq(-limit_grid, limit_grid, by = h)
  x2s <- seq(-limit_grid, limit_grid, by = h)
  X_pred <- expand.grid(x1s, x2s)
  
  # Find the closest point in the grid to the points in the data
  plot_list <- lapply(1:length(id_to_plot), function(i){
    
    # Map grid to a rectagnel
    data_i <- data[data$id == id_to_plot[i], ]
    u_grid <- do.call(rbind, parallel::mclapply(1:dim(X_pred)[1], function(j){
      best_index <- which.min((data_i$x - X_pred$Var1[j])^2 + (data_i$y - X_pred$Var2[j])^2)[1]
      return(data_i$u1[best_index])}, mc.cores = n_cores))
    
    # Plot the final results using viridis palette
    data_to_plot <- tibble(x = X_pred$Var1, 
                           y = X_pred$Var2, 
                           u1 = u_grid)
    
    # If circle fix zeros 
    if(is_circle){
      data_to_plot[data_to_plot$x^2 + data_to_plot$y^2 > 2.5^2, 3] <- 0
    }
    
    p1 <- ggplot(data_to_plot, aes(x, y, fill = u1)) + 
      geom_raster(interpolate = T) +
      scale_fill_viridis_c() + 
      labs(x = "x", y = "y") + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 
    return(p1)})
  
  # Plot the list and write to file   
  p1 <- ggpubr::ggarrange(plotlist = plot_list, ncol = 5, nrow = 2, common.legend = T, legend = "bottom")
  p1
  ggsave(path_save, height = 6, width = 12)
  
}

# Check if the directory to save the figures in exists 
check_if_dir_exists(path_dir <- "../../Result/Schankenberg/Rectangle_figures/")
check_if_dir_exists(path_dir <- "../../Result/Schankenberg/Circle_figures/")
check_if_dir_exists(path_dir <- "../../Result/Gierer/Rectangle_figures/")
check_if_dir_exists(path_dir <- "../../Result/Gierer/Circle_figures/")

# ===================================================================================================
# Working with maximum concentrations rectangles and circles 
# ===================================================================================================
# Schankenberg-model 
dir_files <- "../../Intermediate/Schankenberg_files/Rectangles/"
geometry <- "rectangle"
path_save <- "../../Result/Schankenberg/Rectangle_figures/Max_conc.pdf"
plot_max_conc_data(dir_files, geometry, path_save)

dir_files <- "../../Intermediate/Schankenberg_files/Circles/"
geometry <- "circle"
path_save <- "../../Result/Schankenberg/Circle_figures/Max_conc.pdf"
plot_max_conc_data(dir_files, geometry, path_save)

# ---------------------------------------------------------------------------------------------------
# Last time-point data rectangles 
# ---------------------------------------------------------------------------------------------------
# Zero holes
limit_grid <- 2.2
path_data <- "../../Intermediate/Schankenberg_files/Rectangles/Zero_holes_files/t_end_data.csv"
path_save <- "../../Result/Schankenberg/Rectangle_figures/Zero_holes_end.pdf"
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F)
# Several end-points 
path_save <- "../../Result/Schankenberg/Rectangle_figures/Zero_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid)

# Five holes 
path_data <- "../../Intermediate/Schankenberg_files/Rectangles/Five_holes_files/t_end_data.csv"
path_save <- "../../Result/Schankenberg/Rectangle_figures/Five_holes_end.pdf"
limit_grid <- 2.2
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F)
# Several end-points 
path_save <- "../../Result/Schankenberg/Rectangle_figures/Five_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid)

# Twenty holes 
path_data <- "../../Intermediate/Schankenberg_files/Rectangles/Twenty_holes_files/t_end_data.csv"
path_save <- "../../Result/Schankenberg/Rectangle_figures/Twenty_holes_end.pdf"
limit_grid <- 2.2
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F)
# Several end-points 
path_save <- "../../Result/Schankenberg/Rectangle_figures/Twenty_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid)

# ---------------------------------------------------------------------------------------------------
# Last time-point data circles
# ---------------------------------------------------------------------------------------------------
# Zero holes
path_data <- "../../Intermediate/Schankenberg_files/Circles/Zero_holes_files/t_end_data.csv"
path_save <- "../../Result/Schankenberg/Circle_figures/Zero_holes_end.pdf"
limit_grid <- 2.6
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F, is_circle = T)
# Several end-points 
path_save <- "../../Result/Schankenberg/Circle_figures/Zero_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid, is_circle = T)

# Five holes 
path_data <- "../../Intermediate/Schankenberg_files/Circles/Five_holes_files/t_end_data.csv"
path_save <- "../../Result/Schankenberg/Circle_figures/Five_holes_end.pdf"
limit_grid <- 2.6
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F, is_circle = T)
# Several end-points 
path_save <- "../../Result/Schankenberg/Circle_figures/Five_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid, is_circle = T)

## Twenty holes 
path_data <- "../../Intermediate/Schankenberg_files/Circles/Twenty_holes_files/t_end_data.csv"
path_save <- "../../Result/Schankenberg/Circle_figures/Twenty_holes_end.pdf"
limit_grid <- 2.6
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F, is_circle = T)
# Several end-points 
path_save <- "../../Result/Schankenberg/Circle_figures/Twenty_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid, is_circle = T)
  

# ===================================================================================================
# Working with maximum concentrations rectangles and circles Gierer 
# ===================================================================================================
dir_files <- "../../Intermediate/Gierer_files/Rectangles/"
geometry <- "rectangle"
path_save <- "../../Result/Gierer/Rectangle_figures/Max_conc.pdf"
plot_max_conc_data(dir_files, geometry, path_save)

dir_files <- "../../Intermediate/Gierer_files/Circles/"
geometry <- "circle"
path_save <- "../../Result/Gierer/Circle_figures/Max_conc.pdf"
plot_max_conc_data(dir_files, geometry, path_save)

# ---------------------------------------------------------------------------------------------------
# Last time-point data rectangles Gierer 
# ---------------------------------------------------------------------------------------------------
# Zero holes
limit_grid <- 2.2
path_data <- "../../Intermediate/Gierer_files/Rectangles/Zero_holes_files/t_end_data.csv"
path_save <- "../../Result/Gierer/Rectangle_figures/Zero_holes_end.pdf"
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F)
# Several end-points 
path_save <- "../../Result/Gierer/Rectangle_figures/Zero_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid)

# Five holes 
path_data <- "../../Intermediate/Gierer_files/Rectangles/Five_holes_files/t_end_data.csv"
path_save <- "../../Result/Gierer/Rectangle_figures/Five_holes_end.pdf"
limit_grid <- 2.2
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F)
# Several end-points 
path_save <- "../../Result/Gierer/Rectangle_figures/Five_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid)

# Twenty holes 
path_data <- "../../Intermediate/Gierer_files/Rectangles/Twenty_holes_files/t_end_data.csv"
path_save <- "../../Result/Gierer/Rectangle_figures/Twenty_holes_end.pdf"
limit_grid <- 2.2
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F)
# Several end-points 
path_save <- "../../Result/Gierer/Rectangle_figures/Twenty_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid)

# ---------------------------------------------------------------------------------------------------
# Last time-point data circles Gierer 
# ---------------------------------------------------------------------------------------------------
# Zero holes
path_data <- "../../Intermediate/Gierer_files/Circles/Zero_holes_files/t_end_data.csv"
path_save <- "../../Result/Gierer/Circle_figures/Zero_holes_end.pdf"
limit_grid <- 2.6
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F, is_circle = T)
# Several end-points 
path_save <- "../../Result/Gierer/Circle_figures/Zero_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid, is_circle = T)

# Five holes 
path_data <- "../../Intermediate/Gierer_files/Circles/Five_holes_files/t_end_data.csv"
path_save <- "../../Result/Gierer/Circle_figures/Five_holes_end.pdf"
limit_grid <- 2.6
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F, is_circle = T)
# Several end-points 
path_save <- "../../Result/Gierer/Circle_figures/Five_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid, is_circle = T)

## Twenty holes 
path_data <- "../../Intermediate/Gierer_files/Circles/Twenty_holes_files/t_end_data.csv"
path_save <- "../../Result/Gierer/Circle_figures/Twenty_holes_end.pdf"
limit_grid <- 2.6
plot_t_end_data(path_data, path_save, limit_grid, plot_result = F, is_circle = T)
# Several end-points 
path_save <- "../../Result/Gierer/Circle_figures/Twenty_holes_end_several.pdf"
plot_several_end_times(path_data, path_save, limit_grid, is_circle = T)
