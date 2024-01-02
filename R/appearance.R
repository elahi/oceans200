# Constants to produce uniform look across figures and to allow changes as necessary
sc="grey60"   #Symbol colour
ac="grey60"   # Axis colour
lf="grey80"   # bars; light fill
lc="black"    #line color - use in single plots
ls=1.15       #line thickness for geom_line
ss=2          #symbol size
fs = 0.4      #starting point for bar fills
fe = 0.95     #end point for bar fills
pw = 7.5      #plot width
pww = 15      #plot width for wide plots
pwm = 10      #plot width for 1.5 column plots
ph = 7.5      #plot height
phm = 10      #plot height for 1.5 rows
phh = 15      #plot height for double height plots
sym2=c(0,1)   # symbol set for plot with 2 symbols with open symbols
sym3=c(0,1,2)
sym4=c(0,1,2,5)
sym2s=c(15,16)   #symbol sets when symbols are solid
sym3s=c(15,16,17)
sym4s=c(15,16,17,18)
theme_qk <- function(){ 
  font <- "sans"   #assign font family up front
  theme_classic() %+replace%    #replace elements we want to change
    theme(
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_line(color = ac), #grey axis ticks
      axis.line = element_line(color = ac),
      #text elements
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10),               #font size
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 9,                 #font size
        color=ac),                #font color       
    )
}
