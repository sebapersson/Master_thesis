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
plot_t_end_data <- function(path_data, path_save, limit_grid, plot_result=T, n_cores=3, is_circle=F, h=0.04, seed=123)
{
  set.seed(seed)
  
  # Read the data and aggregate, only check u1 
  data <- read_csv(path_data, col_types = cols()) %>% 
    select(-X1) %>%
    filter(id_mol == 1)

  # Make a grid 
  x1s <- seq(-limit_grid, limit_grid, by = h)
  x2s <- seq(-limit_grid, limit_grid, by = h)
  X_pred <- expand.grid(x1s, x2s)
  # Find the closest point in the grid to the points in the data
  u_grid <- do.call(rbind, parallel::mclapply(1:dim(X_pred)[1], function(i){
    best_index <- which.min((data$x - X_pred$Var1[i])^2 + (data$y - X_pred$Var2[i])^2)[1]
    return(data$u1[best_index])}, mc.cores = n_cores))
  
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
#   file_list, a list of file paths (in order of 0, 5 and 20 holes)
#   geometry, the geometry considered 
#   path_save, the path to where the file should be saved
#   plot_data, if the created plot should be displayed or not, true by default
plot_max_conc_data <- function(file_list, geometry, path_save, plot_data=T)
{
  # Read the data
  data_tot_rec <- do.call(bind_rows, lapply(1:length(file_list), function(i){
    n_holes <- find_number_of_holes(file_list[i])
    data <- read_csv(file_list[i], col_types = cols()) %>%
      mutate(geometry = geometry) %>%
      mutate(n_holes = n_holes)}))
  
  # Aggregate the result   
  data_tot_rec <- data_tot_rec %>%
    mutate(geometry = as.factor(geometry)) %>%
    mutate(n_holes = as.factor(n_holes)) %>%
    select(-X1) %>%
    group_by(time, n_holes) %>%
    summarise(med = mean(Max_conc), 
              quant_low = quantile(Max_conc, 0.05),
              quant_high = quantile(Max_conc, 0.95)) 
  
  p1 <- ggplot(data_tot_rec, aes(time, med, color = n_holes, fill = n_holes, linetype = n_holes)) + 
    geom_line(size = 1.2) + 
    geom_ribbon(aes(ymin = quant_low, ymax = quant_high), color = NA, alpha = 0.2) + 
    labs(x = "Time", y = "Maximum concentration") + 
    scale_color_manual(values = cbPalette[-1], name = "Number of holes") + 
    scale_fill_manual(values = cbPalette[-1], name = "Number of holes") + 
    scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Number of holes") + 
    my_theme
  
  ggsave(path_save, plot=p1, height = 6, width = 9)
  
  if(plot_data == T) print(p1)
}

# Function that will extract and return the number of holes for a certain file. 
# Note that the number of holes are found using regular expressions 
# Args: 
#   file_path, the path to the file 
# Returns:
#   the number of holes the mesh has 
find_number_of_holes <- function(file_path)
{
  pattern <- "h([0-9]*)_"
  n_holes <- as.integer(str_match(file_path, pattern)[2])
  return(n_holes)
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
#   circle_data, contains information wheter or not a circle should be drawn 
# Returns:
#   void 
plot_several_end_times <- function(path_data, path_save, limit_grid, h=0.05, n_cores=3, is_circle=F, 
                                   circle_data=F, seed=123)
{
  
  set.seed(seed)
  # Write function to plot a random sample of t-end results (like 10 of them)
  data <- read_csv(path_data, col_types = cols()) %>%
    select(-X1) %>%
    filter(id_mol == 1)
  
  # Using 5 to (else the plot will be to messy)
  n_random_figs <- 5
  n_exp <- length(unique(data$id))
  id_to_plot <- sample(1:n_exp, n_random_figs, replace = F)
  
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
    
    # If a circle should be drawn 
    if(is.data.frame(circle_data)) p1 <- p1 + ggforce::geom_circle(data = circle_data, 
                                                                   mapping = aes(x0=circle_data$x, 
                                                                       y0=circle_data$y,
                                                                       r=circle_data$r), 
                                                                   inherit.aes = F, 
                                                                   size = 1.0)
    
    return(p1)})
  
  # Plot the list and write to file   
  p1 <- ggpubr::ggarrange(plotlist = plot_list, ncol = 5, nrow = 1, common.legend = T, legend = "bottom")
  ggsave(path_save, plot=p1, height = 6, width = 12)
  
}

# Function that will plot the lc-data (the data used to find the optimal 
# mesh size) for a given geometry and model 
# Args:
#   model, the model name
#   geom, the geometry the problem is solved over 
process_lc_data <- function(model, geom)
{
  # Select the directories with lc in the name
  path_files <- str_c("../../Intermediate/Schankenberg", "_files/", geom, "/")
  file_list <- list.files(path_files)[str_detect(list.files(path_files), "lc")]
  
  if(geom == "Rectangles"){
    is_circle <- F
    limit_grid <- 2.2
  }else{
    is_circle = T
    limit_grid <- 2.6
  }
  # Check that figures directory exists 
  fig_dir <- str_c("../../Result/", model, "/", geom, "/Figures")
  check_if_dir_exists(fig_dir)
  lc_dir <- str_c(fig_dir, "/lc_figs/")
  check_if_dir_exists(lc_dir)
  
  # Create each figure   
  for(i in 1:length(file_list)){
    # Check that the directory where the result is stored exists 
    file_path <- str_c(path_files, file_list[i], "/", "t_end_data.csv")
    dir_save <- str_c(lc_dir, "/", file_list[i])
    check_if_dir_exists(dir_save)
    
    # Save the t-end figure to disk 
    path_save <- str_c(dir_save, "/", file_list[i], "t_end_data.pdf")
    plot_t_end_data(file_path, path_save, limit_grid, plot_result = F, is_circle, n_cores = 1)
  }
  
}

# Function that for a certain experimenal condition will create the max-conc, t_end_several
# and t_mean plots for a provided model for all the geometries in the geom_list. 
# Args:
#   model, a string of the model 
#   geom_list, a list of the geometries to investigate 
#   sign, a signuare to get the experiment, note regex to avoid multiple matching
#   save_tag, save tag for the maximum plot (should reflect the experiment)
#   n_cores, the number of cores to use, by default 1
#   draw_circle, wheter or not a circle will be drawn, if it should be a list with circle information for each
#     file that is read. 
process_experiment <- function(model, geom_list, sign, save_tag, hole_list=c(0, 20, 5), n_cores=1, draw_circle=F)
{
  
  # Loop over the geometries 
  for(geom in geom_list){
    path_files <- str_c("../../Intermediate/", model, "_files/", geom, "/")
    file_list <- do.call(c, lapply(sign, function(i) {
      list.files(path_files)[str_detect(list.files(path_files), i)]}))
    
    # Save guard if regex fails
    if(length(file_list) > 3){
      print("File list is longer than three elements")
      print(file_list)
    }
    
    # Check that the directories to save in exist 
    dir_save <- str_c("../../Result/", model, "/", geom, "/Figures/")
    max_conc_dir <- str_c(dir_save, "Max_conc/")
    t_end_sev <- str_c(dir_save, "t_end_sev/")
    t_end_mean <- str_c(dir_save, "t_end_mean/")
    check_if_dir_exists(dir_save)
    check_if_dir_exists(max_conc_dir)
    check_if_dir_exists(t_end_mean)
    check_if_dir_exists(t_end_sev)
    
    # Maximum concentration 
    file_input_list <- str_c(path_files, file_list, "/max_conc.csv")
    path_save <- str_c(max_conc_dir, save_tag, "_max_cond.pdf")
    plot_max_conc_data(file_input_list, geom, path_save, plot_data = F)
    
    if(geom == "Rectangles"){
      is_circle <- F
      limit_grid <- 2.2
    }else{
      is_circle = T
      limit_grid <- 2.6
    }
    
    # The summarised max-conc data, plot for different number of holes 
    for(i in 1:length(file_list)){
      file <- file_list[i]
      
      path_file <- str_c(path_files, file, "/t_end_data.csv")
      path_save_mean <- str_c(t_end_mean, file, "_t_end_mean.pdf")
      path_save_several <- str_c(t_end_sev, file, "_t_end_sev.pdf")
      
      # In the case a circle should be drawn in the t-end several plot 
      if(is_tibble(draw_circle)){ 
        circle_data <- draw_circle[i, ]
      }else{ 
        circle_data <- F
      }
      
      plot_several_end_times(path_file, path_save_several, limit_grid, n_cores=n_cores, 
                             is_circle=is_circle, circle_data=circle_data)
      plot_t_end_data(path_file, path_save_mean, limit_grid, plot_result = F, n_cores=n_cores, 
                      is_circle = is_circle)
    }   
  }
}
  

# Function that will process and plot an illustration case, where the pde in questions 
# is solved in one dimension. 
# Args:
#   path_data, the path to the data 
#   path_save, the path to where the result is saved 
#   t_index, the time-indices to plot (should be five)
plot_illustration_case <- function(path_data, path_save, t_index)
{
  
  # Check that the directory where to save the result exist
  check_if_dir_exists("../../Result/Illustration")
  
  data_tot <- read_csv(path_data, col_types = cols(id_mol = col_factor())) %>%
    select(-X1) %>%
    select(-u2) %>% 
    rename("conc" = "u1") 
  
  t_points <- unique(data_tot$t)[t_index]
  
  plot_list <- lapply(t_points, function(i){
    data_t <- data_tot %>%
      filter(t == i)
    
    # For being able to compare the two cases 
    min_u1 <- min((data_tot %>% filter(id_mol == 1))$conc)
    min_u2 <- min((data_tot %>% filter(id_mol == 2))$conc)
    
    max_u1_tot <- max((data_tot %>% filter(id_mol == 1))$conc)
    min_u1_tot <- min((data_tot %>% filter(id_mol == 1))$conc)
    
    # Scale the data and make ready to plot 
    data_t <- data_t %>%
      mutate(conc_s = case_when(
        id_mol == 1 ~ conc - min_u1, 
        id_mol == 2 ~ conc - min_u2)) %>%
      mutate(id_mol = case_when(
        id_mol == 1 ~ "A", 
        id_mol == 2 ~ "I"))
    
    text_write = str_c("t = ", as.character(round(i, digits = 2)))
    
    p1 <- ggplot(data_t, aes(x, conc_s, color = id_mol, linetype = id_mol)) +
      geom_line(size = 1.1) + 
      annotate("text", x = Inf, y = Inf, label = text_write, hjust = 1.4, vjust = 2) + 
      scale_color_manual(values = cbPalette[-c(1)]) + 
      my_theme + 
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.title.y.left = element_blank(), 
            legend.title = element_blank(), 
            legend.text = element_text(size = 12)) + 
      ylim(0, max_u1_tot) 
    return(p1)})
  
  plot <- ggpubr::ggarrange(plotlist = plot_list, common.legend = T, ncol = 5, legend="bottom")
  ggsave(path_save, plot, height = 6, width = 12)
}


# ===============================================================================================
# Analyse the data 
# ===============================================================================================
# Process the lc-results,
# TODO: Processing of lc-data is currently broken 
#process_lc_data("Schankenberg", "Rectangles")
#process_lc_data("Schankenberg", "Circles")

# Process non-controlled disturbance case 
geom_list <- c("Rectangles", "Circles")
process_experiment(model="Schankenberg", geom_list=geom_list,sign="_d_k$", save_tag="d_k", n_cores=3)
process_experiment(model="Gierer", geom_list=geom_list,sign="_d_k$", save_tag="d_k", n_cores=3)
  
# Process the case when the initial value is disturbed (at a specific region)
geom_list <- c("Rectangles", "Circles")
circle_data <- tibble(x = c(0, -0.15, 0.55), y = c(0, 0.5, 1.05), r = 0.25)
sign_list <- c("h0_dr0D25x0y0_k$", "h5_dr0D25x-0D15y0D5_k$", "h20_dr0D25x0D55y1D05_k$")
process_experiment(model="Schankenberg", geom_list=geom_list, sign=sign_list, save_tag="d_rxy", 
                   n_cores=3, draw_circle = circle_data)
process_experiment(model="Gierer", geom_list=geom_list, sign=sign_list, save_tag="d_rxy", 
                   n_cores=3, draw_circle = circle_data)

# Information on where to draw the holes for the following plots 
hole_five <- tibble(x = 0, y = 0, r = 0.25)
hole_twenty <- tibble(x = 0.75, y = 0.85, r = 0.25)
hole_both <- hole_five %>% bind_rows(hole_twenty)

# Process when a subdomain has different parameters, note different models have diff param
geom_list <- c("Rectangles", "Circles")
sign_list <- c("h5_d_ka0D5b2D0ga10di100$", "h20_d_ka0D5b2D0ga10di100$")
process_experiment(model="Schankenberg", geom_list=geom_list, sign=sign_list, save_tag="d_k_sub", 
                   n_cores=3, draw_circle = hole_both)
sign_list <- c("h5_d_ka0D6b2D0ga20di50$", "h20_d_ka0D6b2D0ga20di50$")
process_experiment(model="Gierer", geom_list=geom_list, sign=sign_list, save_tag="d_k_sub", 
                   n_cores=3, draw_circle = hole_both)


# Process with subdomain parameters and stronger diffusion 
geom_list <- geom_list <- c("Rectangles", "Circles")
sign_list <- c("h5_d_ka0D6b2D0ga10di200_Kd200$")
process_experiment(model="Schankenberg", geom_list=geom_list, sign=sign_list, save_tag="test", 
                   n_cores=3, draw_circle = hole_five)
sign_list <- c("h5_d_ka0D65b2D0ga20di100_Kdi100$")
process_experiment(model="Gierer", geom_list=geom_list, sign=sign_list, save_tag="test", 
                   n_cores=3, draw_circle = hole_five)

# Process with very large activation 
geom_list <- geom_list <- c("Rectangles", "Circles")
sign_list <- c("h5_d_ka2D0b2D0ga10di100$")
process_experiment(model="Schankenberg", geom_list=geom_list, sign=sign_list, save_tag="test", 
                   n_cores=3, draw_circle = hole_five)
sign_list <- c("h5_d_ka2D0b2D0ga20di50$")
process_experiment(model="Gierer", geom_list=geom_list, sign=sign_list, save_tag="test", 
                   n_cores=3, draw_circle = hole_five)
 

# Process with very large activation and decreased break-down 
geom_list <- geom_list <- c("Rectangles", "Circles")
sign_list <- c("h5_d_ka2D5b0D2ga10di100_Kga5d800$")
process_experiment(model="Schankenberg", geom_list=geom_list, sign=sign_list, save_tag="test", 
                   n_cores=3, draw_circle = hole_five)
sign_list <- c("h5_d_ka2D0b0D5ga20di50$")
process_experiment(model="Gierer", geom_list=geom_list, sign=sign_list, save_tag="test", 
                   n_cores=3, draw_circle = hole_five)


# Process when both inital perturbation and parameters change 
geom_list <- geom_list <- c("Rectangles", "Circles")
sign_list <- c("h5_dr0D25x0D0y0D0_ka0D3b2D0ga15di50$", "h20_dr0D25x0D75y0D85_ka0D3b2D0ga15di50$")
process_experiment(model="Schankenberg", geom_list=geom_list, sign=sign_list, save_tag="test", 
                   n_cores=3, draw_circle = hole_both)
sign_list <- c("h5_dr0D25x0D0y0D0_ka0D6b2D0ga20di50", "h20_dr0D25x0D75y0D85_ka0D6b2D0ga20di50")
process_experiment(model="Gierer", geom_list=geom_list, sign=sign_list, save_tag="test", 
                   n_cores=3, draw_circle = hole_both)


# Process the illustration example 
# Schankenberg 
path_data <- "../../Intermediate/Illustration/Schankenberg_illustration.csv"
path_save <- "../../Result/Illustration/Schankenberg.pdf"
t_index <- c(1, 20, 40, 50, 74)
plot_illustration_case(path_data, path_save, t_index)

# Gierer 
path_data <- "../../Intermediate/Illustration/Gierer_illustration.csv"
path_save <- "../../Result/Illustration/Gierer.pdf"
t_index <- c(1, 80, 85, 90, 100)
plot_illustration_case(path_data, path_save, t_index)
