  #include "colors.inc" 
  
  
  
  background { color White }
  camera {
    location <175, 85, 175>
    look_at  <30, -40, 30>
  } 

#declare scale_factor = 1.17;
#declare Width = 125;
#declare Height0 = 40;

box {
    <-scale_factor*Width, -Height0/2, -scale_factor*Width>,  // Near lower left corner
    < scale_factor*Width, Height0/2,  scale_factor*Width>   // Far upper right corner
   texture {
          pigment { rgb<101/255, 32/255, 250/255>} 
        }
  }
  
  
#declare rad1 = 3;
#declare Height1 = 12*rad1;  
  
#for (X, -Width, Width+16, 14)
  #for (Z, -Width, Width+16, 14)   
      sphere {
        <X, Height0/2+Height1/2, Z>, rad1
        texture {
          pigment { color Copper }
        }    
        no_shadow
      }
  #end
#end   
       
box {
    <-scale_factor*Width, Height0/2, -scale_factor*Width>,  // Near lower left corner
    < scale_factor*Width, Height0/2+Height1,  scale_factor*Width>   // Far upper right corner
   texture {
          pigment { color rgbt<0, 0.5, 0, 0.7> }  
      finish {
      ambient .01
    }
        }
  }

 
#declare Height3 = 10;
box {
    <-scale_factor*Width, Height0/2+Height1, -scale_factor*Width>,  // Near lower left corner
    < scale_factor*Width, Height0/2+Height1+Height3,  scale_factor*Width>   // Far upper right corner
   texture {
          pigment { rgb<101/255, 32/255, 250/255> }
        }
  }     
 
#declare rad2 = 2; 

#for (X, -Width, Width+20, 14)
  #for (Z, -Width, Width+20, 14)   
      sphere {
        <X, Height0/2+Height1+Height3+rad2, Z>, rad2
        texture {
          pigment { color Silver }
        } 
        no_shadow 
      }
  #end
#end 

light_source {<175, 70, 175> color White} 
light_source {<0, 1000, 0> color White}
light_source {<800, 130, 0> color White}  
//light_source {<10, 300, 300> color White}          
