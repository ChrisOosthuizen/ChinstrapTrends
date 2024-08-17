
## Specify MCMCglmm mixed model prior (Krüger 2023)

prior<- list(R = list(V = 1, nu = 0.002),
             G = list(G1 = list(V = diag(2), nu = 0.002,
                                alpha.mu = rep(0, 2),
                                alpha.V= diag(133, 2, 2))))

## Fit MCMCglmm mixed model (Krüger 2023)

mc1<-MCMCglmm(nests~season_starting, random=~us(1 + Lat):site_id, rcov=~units, 
              family="poisson", mev=NULL,
              data=nestM3, start=NULL, nodes="ALL", scale=TRUE, 
              nitt=13000, thin=10, burnin=3000, 
              pr=T, pl=FALSE, verbose=TRUE, DIC=TRUE, singular.ok=FALSE, 
              saveX=TRUE,prior=prior, saveZ=TRUE, saveXL=TRUE, slice=FALSE,
              ginverse=NULL, trunc=FALSE)


all_sites_list = list()
random_sites_list = list()


for (i in 1:9) {

# construct an hypothetical dataframe to generate the populations estimaters
years<-data.frame(season_starting=c(1960:2020)) 
pops<-data.frame(site_id=countsN$site_id[countsN$ncounts>1],
                 Lat=countsN$Lat[countsN$ncounts>1])
popy<-merge(pops,years)
popy$nests<-c(0) ### MCMCglmm needs a column with the response variable
popypred<-data.frame(predict(mc1,newdata=popy,type="response",
                             marginal=mc1$Random$formula,
                             interval="prediction", posterior="mean"))
popy$fit<-popypred$fit

# WCO: Add lower and upper prediction intervals to the data used for inference
popy$lwr<-popypred$lwr
popy$upr<-popypred$upr

random_sites = popy %>% 
             dplyr::filter(site_id =='PETE'| 
                           site_id =='PING'|
                           site_id =='PPPT'|
                           site_id =='PRAE'|
                           site_id =='RAYN'|
                           site_id =='RENI'|
                           site_id =='ROQU'|
                           site_id =='RUGG')

unique(random_sites$site_id)


random_sites_nestm3 = nestm3 %>% 
  dplyr::filter(site_id =='PETE'| 
                  site_id =='PING'|
                  site_id =='PPPT'|
                  site_id =='PRAE'|
                  site_id =='RAYN'|
                  site_id =='RENI'|
                  site_id =='ROQU'|
                  site_id =='RUGG')


random_sites_plot = ggplot(data = random_sites) + 
          geom_line(aes(x = season_starting, y = fit), col = "steelblue",linewidth=1.04) + 
          geom_line(aes(x = season_starting, y = lwr), col = "steelblue1", 
                    linetype="dotted", linewidth = 1.02) + 
          geom_line(aes(x = season_starting, y = upr), col = "steelblue1",
                    linetype="dotted",linewidth = 1.02) + 
          geom_point(data = random_sites_nestm3, aes(season_starting, y = nests), 
                     color = "red", cex = 2) + 
          geom_line(data = random_sites_nestm3, aes(season_starting, y = nests), 
                    color = "red",linewidth=0.8) +
          theme_bw() + 
          xlab("Year") +
          ylab("Predicted count") +  
          labs(subtitle = "Krüger (2023)") + 
          facet_wrap(~ site_id, ncol = 2, nrow = 4,
                              scales = 'free') 
  

random_sites_list[[i]] = random_sites_plot

# Blue solid lines are the predicted abundance (posterior mean) used by Krüger (2023)
# to predict regional declines. Light blue dots are the 95 % Highest Posterior Density
# interval for this prediction. Red points are the observed counts 
# (connected with a red line).

## Figure 3A (Krüger 2023)


p1v2<-ggplot(popy,aes(season_starting,fit/1000))+
  geom_smooth()+
  geom_point(alpha=0.15)+xlab("Year")+
  theme_bw()+th+ylab("Thousand nests")+
  ggtitle(label="a. Predicted count of nests")+
  scale_y_log10() # plot from the predicted fit

p1v2

all_sites_list[[i]] = p1v2

}

# Assume your list of ggplot2 figures is called 'plot_list'
combined_plot <- cowplot::plot_grid(plotlist = all_sites_list, ncol = 3)

# Assume your list of ggplot2 figures is called 'plot_list'
combined_plot_sites <- cowplot::plot_grid(plotlist = random_sites_list, ncol = 4)

library(patchwork)

# Save Plot 
pdf("./figure/repeat_predictions_sites.pdf", useDingbats = FALSE, width = 12, height = 12)
combined_plot_sites
dev.off()

ggsave("./figure/repeat_predictions.png", combined_plot, width = 12, height = 12, dpi = 1200)



# Assume your list of ggplot2 figures is called 'plot_list'
combined_plot_sites <- random_sites_list[[1]] +  random_sites_list[[2]] + 
                       random_sites_list[[3]] +  random_sites_list[[4]] + 
                       plot_layout(ncol = 2)

ggsave("./figure/repeat_predictions_sites.png", combined_plot_sites, width = 12, height = 12, dpi = 1200)



## Fit a better GLMM 

# Covariates should be standardized. 
nestM3$Zseason_starting = scale(nestM3$season_starting)
nestM3$ZLat = scale(nestM3$Lat)

mc2 <- MCMCglmm(nests ~ Zseason_starting * ZLat, 
                random=~us(1 + Zseason_starting):site_id, 
                rcov=~units,
                family="poisson", mev=NULL,
                data=nestM3,start=NULL, nodes="ALL", scale=TRUE, 
                nitt=23000, thin=10, burnin=13000, pr=T,
                pl=FALSE, verbose=TRUE, DIC=TRUE, singular.ok=FALSE, saveX=TRUE,
                prior=prior, saveZ=TRUE, saveXL=TRUE, slice=FALSE, 
                ginverse=NULL, trunc=FALSE)

## Predict using MCMCglmm mc2
# construct an hypothetical dataframe to predict to
# need to predict to z-standardized variables 
Z1 = dplyr::select(nestM3, season_starting, Lat)
Z2 <- scale(Z1)
attr(Z2,"scaled:center")
attr(Z2,"scaled:scale")

ave_ss = attr(Z2,"scaled:center")[[1]]
ave_lat = attr(Z2,"scaled:center")[[2]]

sd_ss = attr(Z2,"scaled:scale")[[1]]
sd_lat = attr(Z2,"scaled:scale")[[2]]

years<-data.frame(season_starting=c(1960:2020)) # extrapolate to 1960

pops<-data.frame(site_id=countsN$site_id[countsN$ncounts>1],
                 Lat=countsN$Lat[countsN$ncounts>1])
popy<-merge(pops,years)
popy$nests<-c(0) ### MCMCglmm needs a column with the response variable

popy$Zseason_starting = (popy$season_starting - ave_ss)/sd_ss
popy$ZLat = (popy$Lat - ave_lat)/sd_lat

# Don't extrapolate more than X years
first_last_season = nestM3 %>% 
  dplyr::group_by(site_id) %>%
  dplyr::summarise(minyear = min(season_starting),
                   maxyear = max(season_starting)) %>%
  dplyr::arrange(minyear)
first_last_season 

popy = merge(popy, first_last_season)

length(unique(popy$site_id))

popypred <- data.frame(predict(mc2, 
                               newdata=popy,
                               type="response",
                               marginal=NULL,      # crucial, and not default code.
                               interval="prediction",
                               posterior="all"))

popy$Zfit = popypred$fit
popy$Zlwr = popypred$lwr
popy$Zupr = popypred$upr

## How accurate are the predictions relative to observed data?

## Conditional model predictions

required_n_pages = round(133/16)+1

for(i in 1:required_n_pages){
  
  print(ggplot(data = popy) + 
          geom_line(aes(x = season_starting, y = Zfit), 
                    col = "steelblue", linewidth=1.04) + 
          geom_line(aes(x = season_starting, y = Zlwr), 
                    col = "steelblue1", linetype="dotted", linewidth = 1.02) + 
          geom_line(aes(x = season_starting, y = Zupr), 
                    col = "steelblue1", linetype="dotted", linewidth=1.02) + 
          geom_point(data = nestm3, aes(season_starting, y = nests), 
                     color = "red", cex = 2) + 
          geom_line(data = nestm3, aes(season_starting, y = nests), 
                    color = "red",linewidth=0.8) +
          theme_bw() + 
          xlab("Year") +
          ylab("Predicted count") + 
          #  theme(strip.text = element_text(size = 1.5)) +
          facet_wrap_paginate(~ site_id, ncol = 4, nrow = 4, 
                              page = i,
                              scales = 'free'))}


random_sites = popy %>% 
  dplyr::filter(site_id =='PETE'| 
                  site_id =='PING'|
                  site_id =='PPPT'|
                  site_id =='PRAE'|
                  site_id =='RAYN'|
                  site_id =='RENI'|
                  site_id =='ROQU'|
                  site_id =='RUGG')

random_sites_plot = ggplot(data = random_sites) + 
  geom_line(aes(x = season_starting, y = Zfit), col = "steelblue",linewidth=1.04) + 
  geom_line(aes(x = season_starting, y = Zlwr), col = "steelblue1", 
            linetype="dotted", linewidth = 1.02) + 
  geom_line(aes(x = season_starting, y = Zupr), col = "steelblue1",
            linetype="dotted",linewidth = 1.02) + 
  geom_point(data = random_sites_nestm3, aes(season_starting, y = nests), 
             color = "red", cex = 2) + 
  geom_line(data = random_sites_nestm3, aes(season_starting, y = nests), 
            color = "red",linewidth=0.8) +
  theme_bw() + 
  xlab("Year") +
  ylab("Predicted count") + 
  labs(subtitle = "Current (revised) analysis") + 
  facet_wrap(~ site_id, ncol = 2, nrow = 4,
             scales = 'free')

random_sites_plot 


print_plot = random_sites_list[[1]] +
     theme(panel.grid.major = element_blank(),
     panel.grid.minor = element_blank()) + 
  random_sites_plot + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  plot_layout(nrow = 1, widths = c(3, 3)) +
  plot_annotation(tag_levels = 'A')

print_plot 

ggsave("./figure/Figure 3.png", print_plot, width = 8, height = 6, dpi = 1200)


#  Save Plot 
pdf("./figure/Figure 3.pdf",
      useDingbats = FALSE, width = 8, height = 6)
print_plot
dev.off()
