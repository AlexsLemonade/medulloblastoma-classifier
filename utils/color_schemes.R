# Okabe-Ito is a color-blind friendly color scheme
# Reference: Okabe, M., & Ito, K. (2008). Color universal design (CUD): How to make figures and presentations that are friendly to colorblind people. https://jfly.uni-koeln.de/color/#pallet (Original work published 2002)

# palette.colors returns a named list, which must be unnamed to not match color names to category labels
okabe_ito_colors <- unname(grDevices::palette.colors(palette = "Okabe-Ito"))

# set subgroup colors
subgroup_colors <- c(G3 = okabe_ito_colors[2],
                     G4 = okabe_ito_colors[3],
                     SHH = okabe_ito_colors[4],
                     WNT = okabe_ito_colors[8],
                     Unclassified = "#E0E0E0")

# set platform colors
platform_colors = c(Array = okabe_ito_colors[1],
                    `RNA-seq` = okabe_ito_colors[9])
