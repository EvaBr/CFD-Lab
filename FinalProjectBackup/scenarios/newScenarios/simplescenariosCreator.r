pipe <- function(){
  noslipwall <- matrix(rep(2, 82*52), nrow = 52)
  
  waterwall <-  rbind(rep(2,80), matrix(rep(0, 80*50), nrow = 50), rep(2,80))
  waterwall <- cbind(rep(2,52), waterwall, rep(4,52))
  
  mixed <- waterwall
  mixed[18:35, 1:30] <- 2
  mixed2 <- mixed
  mixed2[51, 10:30] <- 5
  mixed2[51, 9] <-2
  mixed2[51, 31] <- 2
  
  pipe2 <-mixed2
  pipe2[18, 1:30] <- 2
  pipe2[19:34, 1:30] <- 5
  pipe2[35, 1:30] <- 2
  
  pipe <- waterwall
  pipe[18, 1:30] <- 2
  pipe[19:34, 1:30] <- 5
  pipe[35, 1:30] <- 2
  
  array <- matrix(nrow=52*52, ncol=82)
  array[1:52, ] <- noslipwall
  array[53:(17*52), ] <- matrix( rep(t(waterwall), 16), ncol=82 , byrow = TRUE )
  array[(17*52+1):(18*52), ] <- mixed
  array[(18*52+1):(34*52), ] <- matrix( rep(t(pipe) , 16) , ncol=82 , byrow = TRUE )
  array[(34*52+1):(35*52), ] <- mixed
  array[(34*52+1):(35*52), ] <- mixed2  ##added for oral exam, not to be used
  array[(35*52+1):(51*52), ] <- matrix( rep(t(waterwall), 16), ncol=82 , byrow = TRUE )
  array[(51*52+1):(52*52), ] <- noslipwall
  
  
  
  arraybeg = matrix(c("P2", "82 2704", "6"))
  write.table(arraybeg, file="pipeinflow.pgm", sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(array, file="pipeinflow.pgm", append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}  

#########################################################
##pipe outflow, into air/half water:
pipeAir <- function(x,y,z, piplen, fill, fz, outfl){ #outlf=1 or0, whther we have outflow on the right or not.
  noslipwall <- matrix(rep(2, (x+2)*(y+2)), nrow = y+2)
  
  waterwall <-  rbind(rep(2,x), matrix(rep(0, x*y), nrow = y), rep(2,x))
  waterwall <- cbind(rep(2,y+2), waterwall, rep(2*(1+outfl),y+2))
  airwall <- waterwall
  airwall[2:(y+1), 2:(x+1)] <- 1 #make it air
    
  fy <- floor(y/3)
  mixed <- airwall
  mixed[(fy):(2*fy+1), 1:(piplen+1)] <- 2
  
  mixed2 <- mixed
  mixed2[(fy+1):(2*fy), 1:(piplen+1)] <- 5
  #mixed2[fy, 1:(piplen+1)] <-2
  #mixed2[2*fy+1, 1:(piplen+1)] <- 2
  #pipe2 <-mixed2
  #pipe2[18, 1:30] <- 2
  #pipe2[19:34, 1:30] <- 5
  #pipe2[35, 1:30] <- 2
  #pipe <- waterwall
  #pipe[18, 1:30] <- 2
  #pipe[19:34, 1:30] <- 5
  #pipe[35, 1:30] <- 2
  
  #fz = floor(z/3)  #fz=where the pipe starts, with regrard to z coordinate. has to be higher then fill(=water level in pipe)
  fz2 = floor(z/3) #how 'wide' the pipe is in z dimension
  array <- matrix(nrow=(y+2)*(z+2), ncol=x+2)
  array[1:(y+2), ] <- noslipwall
  array[(y+3):((fill+1)*(y+2)), ] <- matrix( rep(t(waterwall), fill), ncol=x+2 , byrow = TRUE )
  array[((fill+1)*(y+2)+1):((fz+1)*(y+2)), ] <- matrix( rep(t(airwall), fz-fill), ncol=x+2 , byrow = TRUE )
  array[((fz+1)*(y+2)+1):((fz+2)*(y+2)), ] <- mixed
  array[((fz+2)*(y+2)+1):((fz+fz2+2)*(y+2)), ] <- matrix( rep(t(mixed2) , fz2) , ncol=x+2 , byrow = TRUE )
  array[((fz+fz2+2)*(y+2)+1):((fz+fz2+3)*(y+2)), ] <- mixed
  
  array[((fz+fz2+3)*(y+2)+1):((z+1)*(y+2)), ] <- matrix( rep(t(airwall), z-2-fz-fz2), ncol=x+2 , byrow = TRUE )
  array[((z+1)*(y+2)+1):((z+2)*(y+2)), ] <- noslipwall
  
  
  
  arraybeg = matrix(c("P2", paste( x+2, (y+2)*(z+2),sep= " "), "6"))
  write.table(arraybeg, file=paste("pipeinflow2", x,y,z,"-", fill, ".pgm",sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(array, file=paste("pipeinflow2", x,y,z,"-", fill, ".pgm",sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}  

#pipeAir(x,y,z, piplen, fill, fz, outfl)
pipeAir(60,30,30, 20, 5, 15, 0)
#pipeAir(60,30,30, 20, 5, 10, 1)


###---------------------------------------------------###
buildCavity <- function(x, y, z){
  noslp <- matrix(rep(2, (x+2)*(y+2)), nrow=y+2);
  wall <- matrix(rep(6, (x+2)*(y+2)), nrow=y+2);
  middle <- noslp
  middle[2:(y+1), 2:(x+1)] <- 0
  simpleCavity <- noslp
  for (i in 2:(z+1)){
    simpleCavity <- rbind(simpleCavity, middle)
  }
  simpleCavity <- rbind(simpleCavity, wall)
  
  arraybeg2 = matrix(c("P2", paste(x+2, (y+2)*(z+2), sep=" "), "6"))
  write.table(arraybeg2, file=paste("cavity", x,y,z, ".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(simpleCavity, file=paste("cavity", x,y,z, ".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}
buildCavity2 <- function(x, y, z){
  noslp <- matrix(rep(2, (x+2)*(y+2)), nrow=y+2);
  wall <- matrix(rep(6, (x+2)*(y+2)), nrow=y+2);
  middle <- noslp
  middle[5:(y-2), 5:(x-2)] <- 0
  simpleCavity <- noslp
  for (i in 2:(z+1)){
    simpleCavity <- rbind(simpleCavity, middle)
  }
  simpleCavity <- rbind(simpleCavity, wall)
  
  arraybeg2 = matrix(c("P2", paste(x+2, (y+2)*(z+2), sep=" "), "6"))
  write.table(arraybeg2, file=paste("cavity", x,y,z, ".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(simpleCavity, file=paste("cavity", x,y,z, ".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}
#buildCavity(30, 30, 30)
#buildCavity(30, 20, 20)
#buildCavity(20, 30, 20)
#buildCavity(20, 20, 30)
#buildCavity2(20,20,2)

tum <- function(x,y,z, xdom, ydom, zdom, kjex, outpt){
  #create the TUM slice in yz, of size x*y*z. Whole domain pipe is of size xdom x ydom x zdom.
  #assumed, that y>=20 is multiple of 10, and z is multiple of 3.  
  hy <- y/10
  hz <- z/3
  ostZ <- (zdom-z)/2
  ostY <- floor((ydom-y)/2)+1
  
  noslp <- matrix(rep(2, (xdom+2)*(ydom+2)), nrow=ydom+2);
  bckgrd <- noslp
  bckgrd[2:(ydom+1), 2:(xdom+1)] <- 0
  bckgrd[2:(ydom+1), 1] <- 5
  bckgrd[2:(ydom+1), xdom+2] <- 4
  
  TUM1 <- bckgrd
  TUM1[(ostY+1):(ostY+hy), kjex:(kjex-1+x)] <- 2
  TUM1[(ostY+hy*2+1):(ostY+3*hy), kjex:(kjex-1+x)] <- 2
  TUM1[(ostY+hy*4+1):(ostY+hy*7), kjex:(kjex-1+x)] <- 2
  TUM1[(ostY+hy*8+1):(ostY+hy*9), kjex:(kjex-1+x)] <- 2
  
  TUM2 <- TUM1
  TUM2[(ostY+hy*5+1):(ostY+hy*6), kjex:(kjex-1+x)] <- 0 
  
  TUM3 <- bckgrd
  TUM3[(ostY+1):(ostY+hy*5), kjex:(kjex-1+x)] <- 2
  TUM3[(ostY+hy*6+1):(ostY+hy*10), kjex:(kjex-1+x)] <- 2
  
  TUM <- rbind(TUM1, TUM2, TUM3)
  matrika <- noslp
  if (ostZ>=1){
    for (i in 1:floor(ostZ)) {
      matrika <- rbind(matrika, bckgrd)
    }
  }
  for (i in 0:2) {
    for (j in 1:hz) {
      matrika <- rbind(matrika, TUM[((ydom+2)*i+1):((ydom+2)*(i+1)), ])
    }
  }
  if (ostZ>0){
    for (i in 1:ceiling(ostZ)) {
      matrika <- rbind(matrika, bckgrd)
    }
  }
  matrika <- rbind(matrika, noslp)

  if(outpt){
    return(matrika)
  } else {
    arraybg = matrix(c("P2", paste(xdom+2, (ydom+2)*(zdom+2), sep=" "), "6"))
    write.table(arraybg, file=paste("tum",y,z,":",kjex,"-",xdom,ydom,zdom,".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(matrika, file=paste("tum",y,z,":",kjex,"-",xdom,ydom,zdom,".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
  }
}


#tum(x,y,z, xdom,ydom,zdom, kjex)
#tum(6, 20, 9, 24, 24, 20, 6)
#tum(10, 20, 20, 80, 20, 20, 30)



#        //sprintf( szFileName, "/media/norbert/940CB6150CB5F27A/Documents/simulation/%s.%i.vtk", szProblem, timeStepNumber );

twopipes <- function(x,y,z){
  noslipwall <- matrix(rep(2, (y+2)*(x+2)), nrow = y+2)
  
  waterwall <-  rbind(rep(2,x), matrix(rep(0, x*y), nrow = y), rep(2,x))
  waterwall <- cbind(rep(2,(y+2)), waterwall, rep(2,(y+2)))
  
  mixed <- waterwall
  mixed[floor(y/3):floor(2*y/3), floor(3/5*x):floor(4/5*x)] <- 2

  mixedend <- mixed
  mixedend[(floor(y/3)+1):(floor(2*y/3)-1), (floor(3/5*x)+1):(floor(4/5*x)-1)] <- 4 #outflow
  
  
  end <- noslipwall
  end2 <- waterwall
  end2[1:floor(y/3), floor(1/5*x):floor(2/5*x)] <- 2
  end2[2:(floor(y/3)-1), (floor(1/5*x)+1):(floor(2/5*x)-1)] <-  5 #inflow
  
  mixed2 <- mixed
  mixed2[61, 10:30] <- 5
  mixed2[61, 9] <-2
  mixed2[61, 31] <- 2
  
  mixed3 <- mixed
  mixed3[61, 9:31] <- 2
  
  matr <- noslipwall
  for (i in 1:ceiling(z/3)){
    matr <- rbind(matr, waterwall)
  }
  matr <- rbind(matr, mixedend)
  for (i in 1:(z-ceiling(z/3)-2)) {
    if (i ==1 ){
      matr <- rbind(matr, mixed3)
    }
    else if(i==21){
      matr <- rbind(matr, mixed3)
    }
    else if (i < 21){
      matr<- rbind(matr, mixed2)
    }
    
    else {matr <- rbind(matr, mixed)}
  }
  #
  matr <- rbind(matr, end2)
  #
  matr <- rbind(matr, end)

  arraybg = matrix(c("P2", paste(x+2, (y+2)*(z+2), sep=" "), "6"))
  write.table(arraybg, file=paste("twopipes",x,y,z,".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(matr, file=paste("twopipes",x,y,z,".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}

#twopipes(50,30,30)
#twopipes(100,60,60)


####################################################################################################3
breakingDam <- function(x,y,z, damx, outpt) {
  noslipwall <- matrix(rep(2, (y+2)*(x+2)), nrow = y+2)
  
  #make inner cells air
  airwall <-  rbind(rep(2,x), matrix(rep(1, x*y), nrow = y), rep(2,x))
  airwall <- cbind(rep(2,(y+2)), airwall, rep(2,(y+2)))

  #make some of it water:
  airwall[2:(y+1), 2:(damx+1)] <- 0 
  
 matrik <- noslipwall
 for (i in 1:z) {
   matrik <- rbind(matrik, airwall)
 }
 matrik <- rbind(matrik, noslipwall)

 if (outpt==1) {
   return(matrik)
 } else {
   arraybg = matrix(c("P2", paste(x+2, (y+2)*(z+2), sep=" "), "6"))
   write.table(arraybg, file=paste("simpleDam",x,y,z,".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
   write.table(matrik, file=paste("simpleDam",x,y,z,".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
 }
}

#breakingDam(60,40,40,20)


######################################################################3
BeerDam <- function(x,y,z, damx, beer, beerheight, fill, size) {
  noslipwall <- matrix(rep(2, (y+2)*(x+2)), nrow = y+2)
  
  #make inner cells air
  airwall <-  rbind(rep(2,x), matrix(rep(1, x*y), nrow = y), rep(2,x))
  airwall <- cbind(rep(2,(y+2)), airwall, rep(2,(y+2)))
  
  #make some of it water:
  airwall[2:(y+1), 2:(damx+1)] <- 0 
  
  #beer. sirok bo size=14 cells
#  size <- 14
  beery = (y-size)/2             #y mora biti SOD!
  
  #beer bottom #once
  beer0 <- airwall
  beer0[(beery+1):(beery+size), (beer+4):(beer+size-5)] <- 2
  beer0[(beery+2):(beery+size-1), (beer+3):(beer+size-4)] <- 2
  beer0[(beery+3):(beery+size-2), (beer+2):(beer+size-3)] <- 2
  beer0[(beery+4):(beery+size-3), (beer+1):(beer+size-2)] <- 2
  beer0[(beery+5):(beery+size-4), beer:(beer+size-1)] <- 2
  #beer middle    #z/4-1 times lower, same upper
  beer1 <- beer0
  beer1[(beery+3):(beery+size-2), (beer+4):(beer+size-5)] <- fill
  beer1[(beery+4):(beery+size-3), (beer+3):(beer+size-4)] <- fill
  beer1[(beery+5):(beery+size-4), (beer+2):(beer+size-3)] <-fill
  

  #beer handle beginning #twice lower, twice upper
  beerH0 <- beer1
  beerH0[(beery+size/2):(beery+size/2+1), (beer+size):(beer+size-1+size/2)] <- 2
  
  #beer handle middle   #z/2-2
  beerH1 <- beer1
  beerH1[(beery+size/2):(beery+size/2+1), (beer+size-2+size/2):(beer+size-1+size/2)] <- 2
  
  #binding it together
  how <- beerheight/4-1
  matrik <- rbind(noslipwall, beer0)
  for (i in 1:(how-1)) {
    matrik <- rbind(matrik, beer1)
  }
  #handle:
  matrik <- rbind(matrik, beerH0, beerH0)
  for (i in 1:(how*2)) {
    matrik <- rbind(matrik, beerH1)
  }
  matrik <- rbind(matrik, beerH0, beerH0)
  #rest of the glass
  for (i in 1:how) {
    matrik <- rbind(matrik, beer1)
  }
  
  #emptyspace
  for (i in 1:(z-beerheight)) {
    matrik <- rbind(matrik, airwall)
  }
  #end wall
  matrik <- rbind(matrik, noslipwall)
  
  arraybg = matrix(c("P2", paste(x+2, (y+2)*(z+2), sep=" "), "6"))
  write.table(arraybg, file=paste("beerDam",x,y,z,"-", beer, beerheight, "-", fill, ".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(matrik,  file=paste("beerDam",x,y,z,"-", beer, beerheight, "-", fill, ".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}

#BeerDam(x,y,z, damx, beer, beerheight, fill, size)
#BeerDam(100,60,60, 40, 51, 24, 0, 14)
#BeerDam(100,60,60, 40, 51, 24, 1, 14)
#BeerDam(100,60,60, 40, 51, 24, 2, 14)
#BeerDam(60,40,40, 20, 31, 12, 1, 12)





#this one, as a change from previous one, writes the handle perpendicularly to x axis.
BeerDam2 <- function(x,y,z, damx, beer, beerheight, fill, size) {
  noslipwall <- matrix(rep(2, (y+2)*(x+2)), nrow = y+2)
  
  #make inner cells air
  airwall <-  rbind(rep(2,x), matrix(rep(1, x*y), nrow = y), rep(2,x))
  airwall <- cbind(rep(2,(y+2)), airwall, rep(2,(y+2)))
  
  #make some of it water:
  airwall[2:(y+1), 2:(damx+1)] <- 0 
  
  #beer. sirok bo size=14 cells
  #  size <- 14
  beery = (y-size)/2             #y mora biti SOD!
  
  #beer bottom #once
  beer0 <- airwall
  beer0[(beery+1):(beery+size), (beer+4):(beer+size-5)] <- 2
  beer0[(beery+2):(beery+size-1), (beer+3):(beer+size-4)] <- 2
  beer0[(beery+3):(beery+size-2), (beer+2):(beer+size-3)] <- 2
  beer0[(beery+4):(beery+size-3), (beer+1):(beer+size-2)] <- 2
  beer0[(beery+5):(beery+size-4), beer:(beer+size-1)] <- 2
  #beer middle    #z/4-1 times lower, same upper
  beer1 <- beer0
  beer1[(beery+3):(beery+size-2), (beer+4):(beer+size-5)] <- fill
  beer1[(beery+4):(beery+size-3), (beer+3):(beer+size-4)] <- fill
  beer1[(beery+5):(beery+size-4), (beer+2):(beer+size-3)] <-fill
  
  
  #beer handle beginning #twice lower, twice upper
  beerH0 <- beer1
  beerH0[(beery+size+1):(beery+size+size/2), (beer+size/2-1):(beer+size/2)] <- 2
  
  #beer handle middle   #z/2-2
  beerH1 <- beer1
  beerH1[(beery+size+size/2-1):(beery+size+size/2), (beer+size/2-1):(beer+size/2)] <- 2
  
  #binding it together
  how <- beerheight/4-1
  matrik <- rbind(noslipwall, beer0)
  for (i in 1:(how-1)) {
    matrik <- rbind(matrik, beer1)
  }
  #handle:
  matrik <- rbind(matrik, beerH0, beerH0)
  for (i in 1:(how*2)) {
    matrik <- rbind(matrik, beerH1)
  }
  matrik <- rbind(matrik, beerH0, beerH0)
  #rest of the glass
  for (i in 1:how) {
    matrik <- rbind(matrik, beer1)
  }
  
  #emptyspace
  for (i in 1:(z-beerheight)) {
    matrik <- rbind(matrik, airwall)
  }
  #end wall
  matrik <- rbind(matrik, noslipwall)
  
  arraybg = matrix(c("P2", paste(x+2, (y+2)*(z+2), sep=" "), "6"))
  write.table(arraybg, file=paste("beerDam2-",x,y,z,"-", beer, beerheight, "-", fill, ".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(matrik,  file=paste("beerDam2-",x,y,z,"-", beer, beerheight, "-", fill, ".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}


#BeerDam2(60,40,40, 20, 31, 12, 1, 12)
#BeerDam2(100,60,60, 40, 51, 24, 0, 14)
#BeerDam2(100,60,60, 40, 51, 24, 1, 14)
#BeerDam2(100,60,60, 40, 51, 24, 2, 14)



#the usual pipe scenario, but with the beer
BeerPipe2 <- function(x,y,z, beer, beerheight, fill, size) {
  noslipwall <- matrix(rep(2, (y+2)*(x+2)), nrow = y+2)
  
  #make inner cells water
  airwall <-  rbind(rep(2,x), matrix(rep(0, x*y), nrow = y), rep(2,x))
  airwall <- cbind(rep(5,(y+2)), airwall, rep(4,(y+2)))
  
  #beer. sirok bo size=14 cells
  #  size <- 14
  beery = (y-size)/2             #y mora biti SOD!
  
  #beer bottom #once
  beer0 <- airwall
  beer0[(beery+1):(beery+size), (beer+4):(beer+size-5)] <- 2
  beer0[(beery+2):(beery+size-1), (beer+3):(beer+size-4)] <- 2
  beer0[(beery+3):(beery+size-2), (beer+2):(beer+size-3)] <- 2
  beer0[(beery+4):(beery+size-3), (beer+1):(beer+size-2)] <- 2
  beer0[(beery+5):(beery+size-4), beer:(beer+size-1)] <- 2
  #beer middle    #z/4-1 times lower, same upper
  beer1 <- beer0
  beer1[(beery+3):(beery+size-2), (beer+4):(beer+size-5)] <- fill
  beer1[(beery+4):(beery+size-3), (beer+3):(beer+size-4)] <- fill
  beer1[(beery+5):(beery+size-4), (beer+2):(beer+size-3)] <-fill
  
  
  #beer handle beginning #twice lower, twice upper
  beerH0 <- beer1
  beerH0[(beery+size+1):(beery+size+size/2), (beer+size/2-1):(beer+size/2)] <- 2
  
  #beer handle middle   #z/2-2
  beerH1 <- beer1
  beerH1[(beery+size+size/2-1):(beery+size+size/2), (beer+size/2-1):(beer+size/2)] <- 2
  
  #binding it together
  how <- beerheight/4-1
  matrik <- rbind(noslipwall, beer0)
  for (i in 1:(how-1)) {
    matrik <- rbind(matrik, beer1)
  }
  #handle:
  matrik <- rbind(matrik, beerH0, beerH0)
  for (i in 1:(how*2)) {
    matrik <- rbind(matrik, beerH1)
  }
  matrik <- rbind(matrik, beerH0, beerH0)
  #rest of the glass
  for (i in 1:how) {
    matrik <- rbind(matrik, beer1)
  }
  
  #emptyspace
  for (i in 1:(z-beerheight)) {
    matrik <- rbind(matrik, airwall)
  }
  #end wall
  matrik <- rbind(matrik, noslipwall)
  
  arraybg = matrix(c("P2", paste(x+2, (y+2)*(z+2), sep=" "), "6"))
  write.table(arraybg, file=paste("beerPipe2-",x,y,z,"-", beer, beerheight, "-", fill, ".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(matrik,  file=paste("beerPipe2-",x,y,z,"-", beer, beerheight, "-", fill, ".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}

#BeerPipe2(100,60,60, 51, 24, 0, 14)

#xdom <-30; ydom <- 20; zdom <- 10
#damx <- 10; tumx <- 21; x<-5;y<-20;z<-9
#outflw<-0
tumDam <- function(xdom, ydom, zdom, damx, tumx, x,y,z, outflw){
  tumi <- tum(x,y,z, xdom,ydom,zdom, tumx, 1)
  dam <- breakingDam(xdom,ydom,zdom, damx, 1)

  #tumi[tumi==0] <- 1 #instead of water, have air
  sum <- tumi+dam
  sum[sum==4] <- 2
  sum[sum==7] <- 2
  sum[sum==6] <- 2*(1+outflw)
 
  arraybg = matrix(c("P2", paste(xdom+2, (ydom+2)*(zdom+2), sep=" "), "6"))
  write.table(arraybg, file=paste("tumDam-",xdom,ydom,zdom,"-", damx, tumx, "-", outflw, ".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sum,  file=paste("tumDam-",xdom,ydom,zdom,"-", damx, tumx, "-", outflw, ".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
}

#debug tryouts:
#tumDam(30,30,30, 10, 21, 10, 30, 30, 0)
#actual cases:
#tumDam(100,60,60, 40, 61, 10, 60, 60, 0)
#tumDam(100,60,60, 40, 61, 10, 60, 60, 1)

cityDam <- function(xd, yd, zd, damx, out){ #, obst){
  noslipwall <- matrix(rep(0, (xd+2)*(yd+2)), nrow = yd+2)
  #waterwall <-  rbind(rep(2,xd), matrix(rep(0, xd*yd), nrow = yd), rep(2,xd))
  #waterwall <- cbind(rep(2,yd+2), waterwall, rep(2*(1+out), yd+2))
  waterwall <- noslipwall; waterwall[1:(yd+2), xd+2] <- 2*out
  
  gemuse <- breakingDam(xd,yd,zd, damx, 1)
  
  #####house on a hill obstacle. assumed: all of it is in air.
  #make the hill
  hrib <- ceiling((xd-damx)/3)+damx
  apfel <- waterwall; apfel[2:(yd-1), hrib:(xd+2)] <- 1
  strudel <- apfel; strudel[2:(yd-1), hrib] <- 0
  ofen <- strudel; ofen[2:(yd-1), hrib+1] <- 0
  teller <- ofen; teller[2:(yd-1), hrib+2] <- 0
  mund <- teller; mund[2:(yd-1), hrib+3] <- 0
  
  #make the haus
  hswidth <- floor((xd-hrib-3)/3)
  hs <- hswidth + hrib + 3
  haus <- waterwall; haus[2:floor(3/2*hswidth), (hs+1):(hs+hswidth)] <- 1   ##haus walls are no slip too
  dach1 <- haus; dach1[2:floor(3/2*hswidth), (hs-1):(hs+hswidth+2)] <- 1
  dach2 <- haus; dach2[2:floor(3/2*hswidth), (hs):(hs+hswidth+1)] <- 1
  dach3 <- haus;
  dach4 <- waterwall; dach4[2:floor(3/2*hswidth), (hs+2):(hs+hswidth-1)] <- 1
  dach5 <- waterwall; dach5[2:floor(3/2*hswidth), (hs+4):(hs+hswidth-3)] <- 1
  
  haushei <- ceiling((zd-8-6)*3/7)
  obst <- rbind(noslipwall, apfel, apfel, apfel, strudel, strudel, ofen, teller, mund) #noslip + 8lines of hill
  for (i in 1:haushei){
    obst <- rbind(obst, haus)
  }
  obst <- rbind(obst, dach1, dach1, dach2, dach3, dach4, dach5)
  rem <- zd-haushei-6-8
  for (i in 1:rem){
    obst <- rbind(obst, waterwall)
  }
  obst <- rbind(obst, noslipwall)
  #######################################################
  ####Tree on same hill obstacle
  #hill ends and tree starts at z=10
  tree <- rbind(noslipwall, noslipwall, noslipwall); tree <- rbind(tree, tree, tree);
  #tree(till leaves) height:
  trehig <- haushei+4
  roots1 <- noslipwall; roots1[(yd-6):(yd-5),(xd-7):(xd-4)] <-1; roots1[(yd-7):(yd-4),(xd-6):(xd-5)] <-1
  for (i in 1:trehig){
    tree <- rbind(tree, roots1)
  }
  roots2 <- noslipwall; roots2[(yd-7):(yd-4),(xd-8):(xd-3)] <-1; roots2[(yd-8):(yd-3),(xd-7):(xd-4)] <-1
  roots3 <- roots2; roots3[(yd-6):(yd-5),(xd-9):(xd-2)] <-1; roots3[(yd-9):(yd-2),(xd-6):(xd-5)] <-1
  endtree <- noslipwall; endtree[(yd-6):(yd-5), (xd-6):(xd-5)] <- 1
  tree <- rbind(tree, roots2, roots2, roots3, roots3, roots3, roots3, roots2, roots1, endtree)
  
  rem <- zd-10-trehig-8
  for (i in 1:(rem+2)){
    tree <- rbind(tree, noslipwall)
  }
  ########################################################
  print(dim(tree)); print(dim(obst)); print(dim(gemuse))
  sum <- obst+gemuse+tree
  
  arraybg = matrix(c("P2", paste(xd+2, (yd+2)*(zd+2), sep=" "), "6"))
  write.table(arraybg, file=paste("cityDam-",xd,yd,zd,"-", damx,"-", out, ".pgm", sep=""), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sum,  file=paste("cityDam-",xd,yd,zd,"-", damx,"-", out, ".pgm", sep=""), append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
  
}
 
#cityDam(x,y,z, damx, outflow)
cityDam(100, 60, 60, 40, 0)