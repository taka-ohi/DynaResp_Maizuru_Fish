####
#### R script for Ohigashi et al (2025)
#### 
#### 2025.05.19 written by Ohigashi
#### R 4.5.0

dir.create("00_SessionInfo")
dir.create("Data")
dir.create("FigCode")
dir.create("Function")


### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/00_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
