# Load necessary library
library(ggplot2)

# Create the data frame
food_data <- data.frame(
  Food = c("Beef", "Pork", "Chicken", "Fish", "Vegetable"),
  Person1 = c(5, 4, 3, 2, 1),  # Preferences for Person 1
  Person2 = c(1, 2, 5, 4, 3)   # Preferences for Person 2
)

# Reshape the data for plotting
food_data_long <- reshape2::melt(food_data, id.vars = "Food", 
                                 variable.name = "Person", 
                                 value.name = "Preference")

# Create the line graph
ggplot(food_data_long, aes(x = Food, y = Preference, group = Person, color = Person)) +
  geom_line(size = 1.2) +                # Line graph
  geom_point(size = 4) +                # Points on the graph
  labs(title = "Food Preferences of Two People",
       x = "Food Type",
       y = "Preference Level",
       color = "Person") +
  theme_minimal()
