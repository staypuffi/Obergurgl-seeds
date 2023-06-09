### plots for report -----

#source("exploration.R")
source("analysis.R", echo = TRUE) #, skip.echo = 5) #skip 5 lines

library(ggtext) # remotes::install_github("clauswilke/ggtext")
library(glue)
library(ggpubr) # to combine ggplots
library(colorRampPalette)  # for color gradients

## - general info -----
se <- function(x) {
  sqrt(sd(x))
}

unique(all_seeds$taxon)
sum(all_seeds$number_seedlings)
sum(all_seeds$number_seedlings)/(nrow(all_seeds)*20)

all_seeds %>%
  group_by(site, area, mowing) %>%
  summarise(n_seeds = sum(number_seedlings)) %>%
  group_by(site, mowing) %>%
  summarise(mean_seeds = mean(n_seeds),
            se = se(n_seeds))

all_seeds %>%
  group_by(site, area, mowing) %>%
  summarise(n_seeds = sum(number_seedlings)) %>%
  group_by(mowing) %>%
  summarise(mean_seeds = mean(n_seeds),
            se = se(n_seeds))


schran_niches %>%
  right_join(my_spec, by = c("taxon" = "species")) %>%
  filter(summer_niche == max(summer_niche))

schran_niches %>%
  filter(taxon %in% unique(all_seeds$taxon)) %>%
  filter(taxon != "Calluna vulgaris") %>%
  write.csv2(file = "data/seedling_species_chars.csv")


## -1- Regression germ_rate over temp & snow -----
area_seeds <- area_seeds %>%
  mutate(mowing = if_else(mowing == "C", "control", "mowed" ))
  

exp(coef(germ_temp_mod))

plot(germ_temp_mod)
summary(germ_temp_mod)
summary(aov(germ_temp_mod))

scale(summary(germ_temp_mod)$coefficients[,1])

ggplot(area_seeds, aes(x=mowing, y=tot_seedlings))+
  geom_boxplot()

write.csv2(summary(germ_temp_mod)$coefficients,
           file = "plots/seedlings over temperature model.csv")


# this would be the estimate at the end of the range
exp(coef(germ_temp_mod))[1]*(exp(coef(germ_temp_mod))[2])^max(area_seeds$summer_mean_from0)

exp(coef(germ_temp_mod)[1]+0*coef(germ_temp_mod)[2])


my_col_mowing <- c("deepskyblue3", "tan3")  


# temperature plot

## with points and ggplot lines
 (temp_plot <- area_seeds %>%
     ggplot(aes(x= summer_mean,
                y= tot_seedlings,
                col = mowing, 
                linetype = mowing,
                fill = mowing))+
     geom_point(size = 3,
                pch = 21,
                col = "black")+
    theme_classic()+
     # stat_smooth(method = "glm",
     #             method.args = list(family = "poisson"),
     #             formula = "y ~ x") +
     # stat_regline_equation(label.x = 13,
     #                       label.y = 150,
     #                       position = "stack",
     #    aes(label =  paste(..eq.label.., sep = "~~~~")),
     #    formula = "y ~ x"
     #  ) +

     geom_smooth(method = "glm",
                 method.args = list(family = "poisson"))+
    xlab("summer temperatures [°C]")+
    ylab("number of seedlings")+
    annotate(geom = "text", x= 15, y= 140,
             size = 4,
             col = "black",
             hjust = 0,
             label =  paste("p value < 0.001\nR² =",
                            round(my_varpart$Proportion[1], 3))
                      )+
     scale_color_manual(values = my_col_mowing)+
     scale_fill_manual(values = my_col_mowing)+
     theme(#axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
       legend.position = "none")
  )

# snowdays plot
(snow_plot <-  area_seeds %>%
    ggplot(aes(x= snow_days, y= tot_seedlings,
               col= mowing, linetype = mowing))+
    geom_point(size = 3)+
    geom_point(pch = 1,
               col = "black",
               size = 3)+
    annotate(geom = "text", 
             x= 80, y= 110,
             size = 4, 
             col = "black",
             hjust = 0,
             label = paste("p value < 0.001\nR² =",
                           round(my_varpart$Proportion[2], 3))
                      )+
    theme_classic()+
    geom_smooth(method = "glm",
                method.args = list(family = "poisson"))+
    xlab("duration of snow cover [days]")+
    ylab("number of seedlings")+
    scale_color_manual(values = my_col_mowing)+
    theme(legend.position = "none")
  
)

# snowmelt plot
max(area_seeds$last_snowday)
(snowmelt_plot <-  area_seeds %>%
    ggplot(aes(x= last_snowday, y= tot_seedlings,
               col= mowing, linetype = mowing))+
    geom_point(size = 3)+
    geom_point(pch = 1,
               col = "black",
               size = 3)+
    annotate(geom = "text", 
             x= date("2022-06-25"), y= 120,
             size = 4, 
             col = "black",
             hjust = 0,
             label = paste("p value < 0.001\nR² =",
                           round(my_varpart$Proportion[2], 3))
    )+
    theme_classic()+
    geom_smooth(method = "glm",
                method.args = list(family = "poisson"))+
    xlab("date of snowmelt")+
    ylab("number of seedlings")+
    scale_color_manual(values = my_col_mowing)+
    theme(legend.position = c(1.3,.5))
  
)
    

 png(filename = "plots/Seedlings perplot over summer temp.png", 
        height = 500, width = 600, units = "px")
    
ggarrange(snow_plot, temp_plot,
          snowmelt_plot,
          nrow = 2, ncol = 2,
          labels = c("a", "b", "c"), hjust = c(-1, +.6),
          widths = c(0.5, 0.5),
          common.legend = F)
  

dev.off()

## -2- Regression seedlings over dist from niche opt   ----
# to annotate plot with mean site temps
site_temp <- sum_temp %>%
  group_by(site) %>%
  summarise(site_summer_mean = round(mean(summer_mean),1)) %>% ungroup()
  

##
min_dist <- niche_seeds %>% 
  group_by(site) %>% 
  summarise(Min = min(dist_summer)) %>% pull()

summary(niche_distoptimum_glmer)

exp(coef(niche_distoptimum_glmer)$site[1,1]+
  coef(niche_distoptimum_glmer)$site[1,3])

## with points
png(filename = "plots/Distance from niche.png", 
    height = 400, width = 600, units = "px")

niche_seeds %>%  left_join(site_temp) %>%
  mutate(
      mowing = if_else(mowing == "C", "control", "mowed"),
      site_lab = as.character(
      glue("Site {site}<br><span style='font-size:12pt'> 
           Summer T.: {site_summer_mean} °C<span>") # size of "Temp."
    )) %>%
  ggplot(aes(x=dist_summer, y= number_seedlings,
             fill = mowing, 
             linetype = mowing, size = mowing))+
   geom_point(size = 2,
              pch = 21,
              alpha = 0.9)+
  geom_smooth(method = "glm", se = T, 
              size = 1.3,
              formula = "y ~ x",
              col = "gray10",
              method.args = list(family = "poisson"))+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(. ~ site_lab, scales = "free_x")+
  ylab("number of seedlings")+
  xlab("distance from specie's optimum temperature")+
theme_classic()+
  scale_fill_manual(values = my_col_mowing)+
  scale_linetype_manual(values = c(1, 3))+
  scale_size_manual(values = c(1.3, 1))+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
          legend.position = c(0.85, 0.25),
        strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_markdown(size = 16, # size of "Site"
                                    lineheight = 1.1))

dev.off()

## only with lines for sites
mycol_site <- c("goldenrod1",
                            "orange1",
                            "ivory4",
                            "turquoise3",
                            "plum3",
                            "gray20")
                            mycol_site <- colorRampPalette(c("goldenrod1", "turquoise3", "plum4"))
                            mycol_site(5)

png(filename = "plots/Distance from niche justlines.png", 
    height = 400, width = 600, units = "px")

niche_seeds %>%  left_join(site_temp) %>%
  mutate(
    mowing = if_else(mowing == "C", "control", "mowed"),
    site_lab = as.character(
      glue("Site {site}<br><span style='font-size:12pt'> 
           Summer T.: {site_summer_mean} °C<span>") # size of "Temp."
    )) %>%
  ggplot(aes(x=dist_summer, y= number_seedlings,
             fill = site_lab, 
             linetype = site_lab,# size = mowing,
             col = site_lab))+
  # geom_point(size = 2,
  #            pch = 21,
  #            alpha = 0.9)+
  geom_smooth(method = "glm", se = T, 
              size = 1.3,
              formula = "y ~ poly(x, 2)",
              method.args = list(family = "poisson"))+
  geom_vline(xintercept = 0, linetype = 2)+
  ylab("number of seedlings")+
  xlab("distance from specie's optimum temperature")+
  theme_classic()+
  scale_color_manual(values = c(mycol_site(5)))+
  scale_fill_manual(values = c(mycol_site(5)))+
  scale_linetype_manual(values = c(1,2,3,1,2,3))+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = c(0.8, 0.65),
        strip.background=element_rect(fill=NA, color=NA),
        legend.text =element_markdown(size = 12, # size of "Site"
                                    lineheight = 1.1),
        legend.title = element_blank())

dev.off()


## table
round(summary(niche_distoptimum_glmer)$coefficients, 3)

write.csv2(round(summary(niche_distoptimum_glmer)$coefficients, 3),
           file = "plots/mixed model niche distance table.csv")


## -3-  species germ_are against niche ------
p_values <- c(round(c(summary(niche_germrate_1)$coefficients[2,4],
                      summary(niche_germrate_2)$coefficients[2,4],
                      summary(niche_germrate_3)$coefficients[2,4],
                      summary(niche_germrate_4)$coefficients[2,4],
                      summary(niche_germrate_5)$coefficients[2,4],
                      summary(niche_germrate_all)$coefficients[2,4]), 5)) 

# to annotate plot with mean site temps
site_temp <- sum_temp %>%
  mutate(site = paste("Site", site)) %>%
  group_by(site) %>%
  summarise(site_summer_mean = round(mean(summer_mean),1)) %>% ungroup() %>%
  add_row(site = "All sites",
          site_summer_mean = mean(sum_temp$summer_mean, na.rm = T)) %>%
  mutate(site = factor(site),
         site_summer_mean = round(site_summer_mean, 1),
         site_lab = as.character(glue("{site} <br><span style='font-size:12pt'> Summer T.: {site_summer_mean} °C<span>")),
         site_lab2 = (glue("bquote(bold({site}) ~ Summer T.: {site_summer_mean} °C)"))
      ) %>%
  add_column(p_values) %>%
  mutate(p_values = if_else(p_values < 0.001, "< 0.001", ""))


## plot
png(filename = "plots/seedlings over niche.png", 
    height = 400, width = 600, units = "px")
niche_seeds_allsites %>%
  mutate(site = "All sites") %>%
  rbind(niche_seeds_overall) %>%
  left_join(site_temp) %>%
  ggplot(aes(x= summer_niche, y = number_seedlings))+
  geom_point(col = "gray70",
             size = 3)+
    geom_point(pch = 1,
               col = "gray20",
               size = 3)+
  #geom_text(aes(label = taxon))+
  geom_smooth(#col = "chartreuse3",
    col = "black",
    method = "glm",
    method.args = list(family = "poisson"),
    formula = "y ~ x")+
  theme_classic()+
  facet_wrap(~site_lab, scales = "free_y")+
  xlab("climatic niche of species [°C]")+
  ylab("number of seedlings")+
  #as all p-vals are < 0.001 unnecesarry
   #geom_text(data = site_temp,
   #          size = 2.5,
   #          aes(x = 8.5, y = 80,
   #              label = paste("p =", p_values)))+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
          strip.background = element_rect(fill=NA, color=NA),
          strip.text = element_markdown(size = 16, # size of "Site"
                                      lineheight = 1.1))
dev.off()

## plot only model lines
mycol_site <- c("goldenrod1",
                "orange1",
                "ivory4",
                "turquoise3",
                "plum3",
                "gray20")
mycol_site <- colorRampPalette(c("goldenrod1", "turquoise3", "plum4"))
mycol_site(5)

png(filename = "plots/seedlings over niche justlines.png", 
    height = 400, width = 600, units = "px")
niche_seeds_allsites %>%
    mutate(site = "All sites",
           number_seedlings = number_seedlings/5) %>%
    rbind(mutate(niche_seeds_overall,
                 site = paste("Site", site))) %>%
    left_join(site_temp) %>%
    ggplot( aes(x = summer_niche,
                y = number_seedlings,
                col = site_lab,
                fill = site_lab,
                #alpha = 0.5,
                linetype = site_lab
    )) +
    theme_classic() +
  xlab("climatic niche of species [°C]")+
  ylab("number of seedlings")+
    geom_smooth(method = "glm", 
                method.args = list(family = "poisson"),
                alpha = 0.3)+
    scale_y_continuous(name="number of seedlings for sites seperately", 
                       sec.axis = sec_axis(~ 5*., name="number of seedlings for all sites"))+
  scale_color_manual(values = c("gray40", mycol_site(5)))+
  scale_fill_manual(values = c("gray40", mycol_site(5)))+
  scale_linetype_manual(values = c(1,2,3,1,2,3))+
  theme(legend.position = c(.20,.7),
        legend.text = element_markdown(size = 15),
        legend.title = element_blank(),
        axis.title = element_text(size = 10))
dev.off()


## tables
summary(niche_germrate_1)

var_names <- c("Intercept", "Temperature niche", "Natural frequency")

nichovertemp_table <- as_tibble(rbind(
  "All sites", round(summary(niche_germrate_all)$coefficients, 3),
  "Site 1", round(summary(niche_germrate_1)$coefficients, 3),
  "Site 2", round(summary(niche_germrate_2)$coefficients, 3),
  "Site 3", round(summary(niche_germrate_3)$coefficients, 3),
  "Site 4", round(summary(niche_germrate_4)$coefficients, 3),
  "Site 5", round(summary(niche_germrate_5)$coefficients, 3))) %>%
  mutate(var_name = c("All sites", var_names,
                  "Site 1", var_names,
                  "Site 2", var_names,
                  "Site 3", var_names,
                  "Site 4", var_names,
                  "Site 5", var_names),
         `p value` = if_else(`Pr(>|z|)` == 0, "< 0.001",
                             `Pr(>|z|)`)) %>%
  select(-`Pr(>|z|)`) %>%
  relocate(var_name, .before = Estimate)
  
  
write.csv2(nichovertemp_table,
           file = "plots/model seedlings over niche.csv")

