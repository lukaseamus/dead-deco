require(tidyverse)
require(magrittr)
require(here)

carbon <- here("Turnover", "Carbon.csv") %>% read_csv()














k <- here("Turnover", "k.csv") %>% 
  read_csv(col_types = list("f", "c", "c", "f", "c", "f", "f")) %T>%
  print()

k %>%
  filter(Kelp, !Wrack, !Dissolved, !Smallscale) %>%
  distinct(Reference) %>%
  print(n = 30)
