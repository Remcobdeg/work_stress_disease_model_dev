filename="comp2018-05-24.15h45"


#checks if a dedicated data output folder is available in the working directory and, if not, ask the user to define the directory
winDir <- "\\\\Cnas.ru.nl/u276198/Surfdrive/BS28A - Major project/simulations"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations"
if (file.exists(winDir)){setwd(winDir)}; if (file.exists(macDir)){setwd(macDir)}; 
if (!file.exists("simoutput")){setwd(tk_choose.dir(getwd(), "Choose folder for storing simulation data"))
  if (!file.exists("simoutput")){dir.create("simoutput")}}

render("makegraphscomp.Rmd",output_file = paste0("simoutput/",filename,".html"))
