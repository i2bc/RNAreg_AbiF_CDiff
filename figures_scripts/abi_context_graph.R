# Load
df <- read.table('abi_COG_up100_graph.txt', sep='\t', header=TRUE)
               
# Define factors
df$COG = factor(df$COG,
                       levels = c("C", "E", "F", "G", "J", "K", "L", "M", "O", "P", "T", "V", "X", "others", "R", "S", "empty"))
df$position = factor(df$position,
                     levels = c("-5","-4","-3","-2","-1","abi","+1","+2","+3","+4","+5"))

# Plot
abi_plot = ggplot2::ggplot(df, ggplot2::aes(x = position,
                        y = abundance,
                        fill = COG)) +
  ggplot2::geom_bar(stat = "identity",
                    position = ggplot2::position_fill()) +
  ggplot2::scale_fill_manual(values = c("red", "blue", "cyan", "green", "yellow", "salmon", "orange", "violet", "pink", "brown1", "coral2", "darkcyan", "black", "grey", "gray", "bisque", "white"),
                             breaks = c("C", "E", "F", "G", "J", "K", "L", "M", "O", "P", "T", "V", "X", "others", "R", "S", "empty")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0, NA))) +
  ggplot2::labs(title = "Distribution of the gene COG classification in the Abi context",
                x = "Position",
                y = "Abundance") +
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
ggplot2::ggsave("abi_context_graph.png", abi_plot)
