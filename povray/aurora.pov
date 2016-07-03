#version 3.7;

// Simple Auroral Emission field of view
// Michael Hirsch
global_settings { 
    assumed_gamma 1.0  
    ambient_light 0
    radiosity {
        media on
  }
}

#declare hist0=
 camera {
    location <0,0,0>
    look_at <0,100,0>
    angle 20 // degrees FOV
    up y
    right x
    rotate <0,0,0>
}

#declare hist1=
 camera {
    location <5,1,0>
    look_at <0,100,0>
    rotate <0,0,0>
    angle 10
    up <0,0,1>
    right <1,0,0>
    rotate <0,0,0>
}

camera{hist1}

#declare auroral_emissions= media {
      emission 0.5
      intervals 20
          density {
            gradient y 
            turbulence <5,.1,1> 
            scale <.02,1,.05>
            color_map {
                [0.0  color rgb <0,    0,   0>]
                [0.25 color rgb <0,    0,   0>] // sharp bottomside cutoff
                [0.3  color rgb <.15, .8,   .1>]  
                [0.375 color rgb <0.1,.6,  .1>]
                [0.45 color rgb <.025,.15,   .00>]
                [0.6  color rgb <.035, .1,   0>]
                [0.65 color rgb <0.04, .05,  0>]
                [0.85 color rgb <.05, .01,  0>]
                [1.0  color rgb <0,   0,    0>]
                }
          }
        }  


#declare arc1=
 box{
 <0, 0, -1> <.2,1,1>
 hollow // won't glow without "hollow" 
 }
 
#declare arc2=
 box{
 <0, 0, -1> <.2,1,1>
 hollow
 }

#declare group1=
    union{
    object{arc1 translate <1+clock*.2, 0, 0> }
    object{arc2 translate <1+clock*1,  0, 0> }
    }
    
#declare aurora1=
object {group1 
    interior {media{auroral_emissions}}
    pigment { rgbt 1}
} 

object {aurora1 
        scale     <1, 50,  20> 
        rotate    0 
        translate <0, 100, 0>
        }

