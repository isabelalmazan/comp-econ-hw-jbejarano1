library(tidyverse)
library(foreign)
library(plm)

df = as.tibble(read.dta('SeatBelts.dta'))
df

reg = lm(fatalityrate ~ sb_useage + speed65 + speed70 + ba08 + drinkage21 + log(income) + age, data = df)
summary(reg)

reg = plm(fatalityrate ~ sb_useage + speed65 + speed70 + ba08 + drinkage21 + log(income) + age, 
          data = df, index = c('state', 'year'), effect = 'individual')
summary(reg)

reg = plm(fatalityrate ~  sb_useage + speed65 + speed70 + ba08 + drinkage21 + log(income) + age, 
          data = df, index = c('state', 'year'), effect = 'twoway')
summary(reg)

