##### Building a map of the Eastern seaboard, before modifying on Adobe Illustrator to add species distribution shadings and making pretty ####
# Author: Daniel Wuitchik

# Libraries
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maptools)

#c lip gshhs_f.b to the broader area that you want to include in your map. gshhs_f.b is the high definition noaa coastline layer
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhs_f.b"

# Crop global map layer to desired extent
sf1 <- getRgshhsMap(gshhs.f.b,
                    xlim = c(-100, -55),
                    ylim = c(20, 55)) %>%
  fortify()

# Read in coordinates of sampling sites
a=read.csv('GPSCoordinates.csv')

#### Full Coast ####

full = ggplot() + 
  geom_polygon(data=sf1,
               aes(x=long, y=lat, group = group),
               fill = "grey70", color='black', lwd = 0.1) +
  geom_point(data=a[c(1:3)],
             aes(x=site_long,
                 y=site_lat,
                 shape = Site),
             size=5,
             col = 'red') +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio=1,
              xlim = c(-87, -66),
              ylim = c(25, 45)) +
  annotation_north_arrow(style = north_arrow_fancy_orienteering) +
  theme_bw()


#### Map of Woods Hole ####
as_map = ggplot() + 
  geom_polygon(data=sf1,
               aes(x=long, y=lat, group = group),
               fill = "grey70",
               color='black',
               lwd = 0.1)+
  geom_point(data=a[c(1:3)],
             aes(x=site_long,
                 y=site_lat,
                 shape = Site),
             size=5,
             col = 'red') +
  scale_shape_manual(values=c(17,17))+
  theme_bw()+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  coord_fixed(ratio=1,
              xlim = c(-71.5,-70),
              ylim = c(41,42.5))

#### Radio Island ####

os_map = ggplot() + 
  geom_polygon(data=sf1,
               aes(x=long,
                   y=lat,
                   group = group),
               fill = "grey70",
               color='black',
               lwd = 0.1) +
  geom_point(data=a[c(1:3)],
             aes(x=site_long,
                 y=site_lat,
                 shape = Site),
             size=5,
             col = 'red') +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio=1,
              xlim = c(-76.9,-76.42),
              ylim = c(34.5,35)) +
  theme_bw()
