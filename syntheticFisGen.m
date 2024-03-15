% range_x = 


fis = mamfis;
fis = addInput(fis,[0 80],"Name","distance_x");
fis = addMF(fis,"distance_x","gaussmf",[2 5]);
fis = addOutput(fis,[0 100],"Name","output");
plotmf(fis,"input",1)