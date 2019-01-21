sys.frame()
commandArgs()
getSrcDirectory()
sourceDir <- getSrcDirectory(function(dummy) {dummy})


library(rstudioapi)    
rstudioapi::getActiveDocumentContext()$path
