library(hexSticker)
library(dplyr)

set.seed(14)

df <- expand.grid(Age = 1:8, Year = 1999:2018) %>%
  mutate(plotresid = rnorm(8 * 20))

x <- df$Year
y <- df$Age
z <- df$plotresid
SAMxlim <- c(min(x) - 1, max(x) + 1)
SAMylim <- c(min(y) - 1, max(y) + 1)
zmax <- max(sqrt(abs(z)), na.rm = TRUE)
SAMcex <- sqrt(abs(z))/zmax*5
neg <- z<0

png("randplot.png", width=480, height=480, pointsize = 12)
par(mar=c(0,0,0,0))
plot(x, y, xlab="", ylab="", xlim=SAMxlim, ylim=SAMylim, type="n", axes=FALSE) 
points(x[neg], y[neg], cex=SAMcex[neg], col=rgb(1, 0, 0, alpha=0.5), pch=19)
points(x[!neg], y[!neg], cex=SAMcex[!neg], col=rgb(0, 0, 1, alpha=0.5), pch=19)
dev.off()

sticker("randplot.png", 
        package="R3", p_size = 48, p_color = "black", p_x=1, p_y=1,
        h_fill="green", h_color="black",
        s_x = 1, s_y = 1, s_width = 0.63, s_height = 0.63, 
        filename = "hexR3.png",
        url = "https://github.com/cmlegault/R3", u_size = 2.5)
