#install.packages("tidyverse")
#install.packages("dplyr")
library(tidyverse)
library(dplyr)
library(readxl)

Journal_comparison_data <- read_excel("https:/vuw-my.sharepoint.com/my/Journal_comparison_data.xlsx")
View(Journal_comparison_data)

View(Journal_comparison_data)
data <- select(Journal_comparison_data, Source, Year, CiteScore)
View(data)

wider_data_citescore <- data %>% 
  pivot_wider(names_from = Year, values_from = CiteScore)
View(wider_data_citescore)

# Recreate a data table to see SJR, SNIP etc..
#data_SJR <- select(Journal_comparison_data, Source, Year, SJR)
#View(data_SJR)

#wider_data_SJR <- data_SJR %>% 
 # pivot_wider(names_from = Year, values_from = SJR)
#View(wider_data_SJR)

# Recreate a data table to see SJR, SNIP etc..
# data_SNIP <- select(Journal_comparison_data, Source, Year, SNIP)
# View(data_SNIP)

#wider_data_SNIP <- data_SJR %>% 
# pivot_wider(names_from = Year, values_from = SNIP)
#View(wider_data_SNIP)

# summary(Journal_comparison_data$Documents)

# Journal_comparison_data
# #drop_na(Source) %>% 
#   group_by(Journal_comparison_data$Source) %>% 
#   summarise(Lower = min(Journal_comparison_data$Citations),
#             Average = mean(Journal_comparison_data$Citations),
#             Upper = max(Journal_comparison_data$Citations),
#             Difference = max(Journal_comparison_data$Citations) - min(Journal_comparison_data$Citations)) %>% 
#   arrange(Average) %>% 
#   view()

# plot(pressure)
# starwars
# ggplot(data=starwars, mapping = aes(x=gender))+ geom_bar()  

ggplot(data = Journal_comparison_data, aes(x=Source))+geom_bar()

Journal_comparison_data %>% 
  ggplot(aes(Year, Citations, color = Source))+
  geom_point(size = 5, alpha = 0.5) + theme_bw() + labs(title = "Year and Citation by Source")

Journal_comparison_data %>% 
  ggplot(aes(Year, Citations, color = Source))+
  geom_point(size = 2, alpha = 0.8) + 
  geom_smooth() +
  facet_wrap(~ Source) +
  theme_bw() + 
  labs(title = "Year and Citation by Source")

# Linear model for year and citations
Journal_comparison_data %>% 
lm(Year ~ Citations, data = .) %>% 
summary()

#scatter plot for the data
ggplot(Journal_comparison_data, aes(Year, Citations))+
  geom_point() 
 # geom_smooth(method = "lm")


# mean of the each journal by citations
Journal_comparison_data %>% 
  group_by(Source) %>% 
summarise(mean(Citations), median(Citations))

#scatter plot with different color for journal names
ggplot(Journal_comparison_data, aes(x = Year, y = Citations, color = Source))+
  geom_point() 
       