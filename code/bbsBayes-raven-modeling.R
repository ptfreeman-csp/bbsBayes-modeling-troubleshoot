#### Script written by: Patrick T. Freeman - patrick[at]csp-inc.org ####

#### LOAD LIBRARIES ####
library(bbsBayes)
library(tidyverse)
library(kableExtra)
library(sf)
library(spacetime)
library(ggpubr)

### Install package as necessary ####
#install.packages("bbsBayes")

#### YOU WILL NEED TO INSTALL JAGS IF YOU DO NOT CURRENTLY HAVE IT ####

### Reference tutorial: https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/Intro.html ####


##You can download BBS data by running fetch_bbs_data. This will save the data to a package-specific directory on your computer. You must agree to the terms and conditions of the data usage before the download will run [type yes at the prompt]. You only need run this function once for each annual update of the BBS database.###

#fetch_bbs_data()


#### Plot strata #### 
laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coordinate reference system for pretty mapping
stratifications <- c("bbs_cws","bbs_usgs","bcr","latlong","state")

strata_map = load_map(stratifications[4]) #selecting the map for "bbs_latlong" stratification
strata_map = st_transform(strata_map,crs = laea) #transforming the coordinate reference system

# mapping using geom_sf()
st_gg = ggplot(data = strata_map)+
  geom_sf()+
  labs(title = paste("stratify(by =",stratifications[4],")"))

### Load the cells to keep and remove (pre-determined by P.T. Freeman based on overlap with desert tortoise range)

target.cells <- read_rds("input/target-latlong-blocks.rds")
remove.cells <- read_rds("input/remove-latlong-blocks.rds")

#Download BBS data
#fetch_bbs_data()

# Stratify the data by 1 degree grid cell 
stratified_data_latlong <- stratify(by = "latlong")

# Prep the data for JAGS modelling
jags_data_1970_2019 <- prepare_data(strat_data = stratified_data_latlong,
                          min_year = 1970,
                          max_year = 2019,
                          min_n_routes = 1,
                          min_mean_route_years = 1,
                          species_to_run = "Common Raven",
                          model = "slope",
                          heavy_tailed = T,
                          strata_rem = remove.cells,
                          sampler = "jags")

# Run the model (dummy parameters for now)
jags_mod_1970_2019 <- run_model(jags_data = jags_data_1970_2019,
                      n_adapt = 10,
                      n_saved_steps = 10,
                      n_iter = 10,
                      n_burnin = 20,
                      n_chains = 3,
                      n_thin = 10,
                      parallel = FALSE,
                      parameters_to_save = c("n","beta","strata"))

jags_mod_1970_2019$Rhat # shows the Rhat values for each monitored parameter
jags_mod_1970_2019$n.eff # shows the Rhat values for each monitored parameter

### Save model object
#write_rds(jags_mod_1970_2019, "/Users/patrickfreeman-csp/Documents/GitHub/blm-pva/output/model_runs/jags_mod_1970_2019.rds")

#### ERROR ARISES HERE ####
# Generate indices at the continental, national, and stratum level
indices <- generate_indices(jags_mod = jags_mod_1970_2019,
                            jags_data = jags_data_1970_2019,
                            regions="stratum")


### Calculate decadal trends ####
## Population trends can be calculated from the output of generate_indices(). The trends are expressed as geometric mean rates of change (%/year) between two points in time.

trends_70_80 <- generate_trends(indices = indices,
                                    Min_year = 1970,
                                    Max_year = 1980) %>%
  dplyr::rename(ST_12 = Region)

trends_80_90 <- generate_trends(indices = indices,
                                Min_year = 1980,
                                Max_year = 1990) %>%
  dplyr::rename(ST_12 = Region)

trends_90_00 <- generate_trends(indices = indices,
                                Min_year = 1990,
                                Max_year = 2000) %>%
  dplyr::rename(ST_12 = Region)

trends_00_10 <- generate_trends(indices = indices,
                                Min_year = 2000,
                                Max_year = 2010) %>%
  dplyr::rename(ST_12 = Region)

trends_10_20 <- generate_trends(indices = indices,
                                Min_year = 2010,
                                Max_year = 2019) %>%
  dplyr::rename(ST_12 = Region)


### Create list of dataframes for decadal trends ###
trend_list <- list(trends_70_80, trends_80_90, trends_90_00, trends_00_10, trends_10_20)


### Perform spatial joins to trend information ####
trends_spatial <- list() 
for(i in 1:length(trend_list)){
  
  trend_df <- trend_list[[i]]
  
  ### Join the trend summary to the spatial object with grid cells 
  
  trend_spatial_df <- left_join(crop, trend_df, by="ST_12")
  
  ### Replace the missing start years
  trend_spatial_df$Start_year <- na.locf(na.locf(trend_spatial_df$Start_year),fromLast=TRUE)
  
  ### Replace the missing end years 
  trend_spatial_df$End_year <- na.locf(na.locf(trend_spatial_df$End_year),fromLast=TRUE)
  
  ### add to output list
  trends_spatial[[i]] <- trend_spatial_df
}

#### Bind all spatially joined trends ####
all_trends <- bind_rows(trends_spatial) %>%
  dplyr::mutate(span = paste0(Start_year,"_",End_year)) %>%
  dplyr::relocate(span, .after=End_year)

## refactor 
all_trends$span <- factor(all_trends$span, levels=c("1970_1980" ,"1980_1990", "1990_2000", "2000_2010","2010_2019"), ordered=T)

#### Plotting sampling effort, relative abundance, and population trends #### 

### Mean number of Routes surveyed in grid cell 
st_gg_routes <- ggplot(data = all_trends, aes(fill=Mean_Number_of_Routes)) +
  geom_sf() + 
  scale_fill_viridis_c() + 
  facet_wrap(~span, ncol=5) + 
  labs(title="Mean routes surveyed") + 
  theme(legend.position="top")

### Relative abundance
st_gg_relabund <- ggplot(data = all_trends, aes(fill=Relative_Abundance)) +
  geom_sf() + 
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA) + 
  facet_wrap(~span, ncol=5) + 
  labs(title="Relative abundance index") + 
  theme(legend.position="top")


### Geometric mean population trend
st_gg_trend = ggplot(data = all_trends, aes(fill=Trend)) +
  geom_sf() + 
  scale_fill_gradient(low = "green", high = "purple", na.value = NA) + 
  facet_wrap(~span, ncol=5) + 
  labs(title="Population Trend (geometric mean)") + 
  theme(legend.position="top")

ggarrange(st_gg_routes,
          st_gg_relabund,
          st_gg_trend,
          nrow=3, ncol=1,
          align="v")


 
