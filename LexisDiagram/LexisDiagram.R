library("tidyverse")
library("patchwork")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

# See this webpage for an example 
# https://maartenbijlsma.com/2017/11/25/homage-to-the-lexis-diagram-april-23-2013/
calendarTime <- c(0:6)
trialTime <- c(0:6)
m <- data.frame(calendarTime, trialTime)
labels <- data.frame(x = 1, y = 3:1, family = c("sans", "serif", "mono"))

lexis <- ggplot() + theme_bw() + 
  geom_blank(data = m, aes(x = .data$calendarTime, y = .data$trialTime)) + 
  # Trial 0 line
  geom_segment(lineend="square" ,size=1.25, aes(x = 0.03, y = 0.03, xend = 1.97, yend = 1.97, colour = "Z = 0")) +
  geom_point(aes(x=0, y=0), colour="red", shape=1, size=4, stroke=1.5) +
  geom_point(aes(x=2, y=2), colour="red", shape=2, size=4, stroke=2) +
  # Trial 1 line
  geom_segment(lineend="square" ,size=1.25, aes(x = 1.03, y = 0.03, xend = 1.97, yend = 0.97, colour = "Z = 0")) +
  geom_point(aes(x=1, y=0), colour="red", shape=1, size=4, stroke=1.5) +
  geom_point(aes(x=2, y=1), colour="red", shape=2, size=4, stroke=2) +
  # Trial 2 line
  geom_segment(lineend="square" ,size=1.25, aes(x = 2.03, y = 0.03, xend = 3.97,
                                                yend = 1.97, colour = "Z = 1")) +
  geom_point(aes(x=2, y=0), colour="blue", shape=1, size=4, stroke=1.5) +
  geom_point(aes(x=4, y=2), colour="blue", shape=4, size=4, stroke=2) +
  
  scale_x_continuous(breaks = c(0:4), limits=c(0,4.3), labels=c(0:4), 
                     minor_breaks = seq(0, 4, 1), 
                     sec.axis = dup_axis(breaks = c(2.286,3.614), labels=c("First\n vaccine dose",
                                                                           addline_format("COVID-19 hospitalization")), name = "")) +
  scale_y_continuous(breaks = c(0:2), limits=c(0,2.5), labels=c(0:2) ,
                     minor_breaks = seq(0, 2, 1)) +
  geom_vline(xintercept =2.286, size = 1, linetype=2) +
  geom_vline(xintercept =3.614, size = 1, linetype=2) + 
  scale_color_manual(labels=c(substitute(paste(Z, " = 0")),
                              substitute(paste(Z, " = 1"))), values=c("red", "blue")) +
  
  labs(x="Calendar time (l), in weeks", y="Trial time (k), in weeks") +
  theme(text = element_text(size = 14), legend.position=c(.14,.8853),
        axis.text.x.top= element_text(size = 14),
        axis.text.x.bottom= element_text(size = 14),
        # axis.text.y.top= element_text(size = 24),
        axis.text.y= element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(colour = "Treatment regimen") + 
  annotate("text",x=1,y=1.11, label="Trial 0", angle=45, size=6) +
  annotate("text",x=1.5,y=.61, label="Trial 1", angle=45, size=6) +
  annotate("text",x=3,y=1.11, label="Trial 2", angle=45, size=6) 

G_j <- c("1","1","1","0","0")
E_j <- c(0,0,0,0,1)
Z_j <- c("0","0","1","-","-")
j <- c(0:4)
DataDf <- data.frame(j,G_j, Z_j)
DataDf %>% 
  pivot_longer(c(G_j, Z_j), names_to = "layer", values_to = "label") 
p2 <- DataDf %>% 
  pivot_longer(c(G_j, Z_j), names_to = "layer", values_to = "label") %>%
  ggplot(aes(x = j)) +
  geom_text(aes(y = factor(layer, levels = c("Z_j", "G_j")), label = label),
            size = 4.5) +
  labs(y = "Data", x = NULL, size =5) +
  scale_x_continuous( limits=c(0,4.3)) +
  scale_y_discrete( labels=c(expression("Z"[j]), 
                             expression("G"[j])
                             )) +
  theme_gray() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y= element_text(size = 14, face="bold"),
        panel.grid = element_blank(), strip.text = element_blank(), 
        axis.title.y= element_text(size = 14),
        ) 

lexis