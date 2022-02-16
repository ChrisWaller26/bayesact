library(cowplot)
library(ggplot2)
library(dplyr)

dummy_data =
  data.frame(
    group = sample(c("A", "B"), 100, T),
    expo = rlnorm(100, 10, 2),
    profit = runif(100, -1, 1)
  )

dummy_data_grouped =
  dummy_data %>%
  group_by(
    group
  ) %>%
  summarise(
    expo        = sum(expo),
    profit_low  = quantile(profit, 0.025),
    profit_q1   = quantile(profit, 0.25),
    profit_q2   = quantile(profit, 0.5),
    profit_q3   = quantile(profit, 0.75),
    profit_high = quantile(profit, 0.975)
  ) %>%
  ungroup()

expo_trans =
  function(x){x / max(dummy_data_grouped$expo) *
      (max(dummy_data_grouped$profit_high) -
         min(dummy_data_grouped$profit_low))
  }

expo_trans_inv =
  function(x){x * max(dummy_data_grouped$expo) /
      (max(dummy_data_grouped$profit_high) -
         min(dummy_data_grouped$profit_low))
  }

dummy_data_grouped =
  dummy_data_grouped %>%
  mutate(
    expo_adj = expo_trans(expo),
    profit_low_adj  = profit_low - min(profit_low),
    profit_q1_adj   = profit_q1 - min(profit_low),
    profit_q2_adj   = profit_q2 - min(profit_low),
    profit_q3_adj   = profit_q3 - min(profit_low),
    profit_high_adj = profit_high - min(profit_low)
  )


profit_plot =
  ggplot() +
  geom_hline(
    aes(yintercept = -min(dummy_data_grouped$profit_low)),
    color = "red"
  ) +
  geom_bar(
    data = dummy_data_grouped,
    aes(x = group,
        y = expo_adj),
    alpha = 0.3,
    stat = "identity"
  ) +
  geom_boxplot(
    data = dummy_data_grouped,
    aes(x = group,
        ymin = profit_low_adj,
        lower = profit_q1_adj,
        middle = profit_q2_adj,
        upper = profit_q3_adj,
        ymax = profit_high_adj,
        color = group),
    stat = "identity",
    show.legend = F,
    width = 0.6
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5) - min(dummy_data_grouped$profit_low),
    labels =
      function(x){
        paste0(100 * round(x + min(dummy_data_grouped$profit_low), 2), "%")},
    sec.axis =
      sec_axis(expo_trans_inv,
               name = "Exposure")
  ) +
  labs(
    x = "Group",
    y = "Profit Margin"
  )

profit_plot

