#version 3.7;
// Flaming Auroral Emission field of view
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
    angle 1.5 // degrees FOV
    up z
    right x
    //rotate <0,0,0>
}

#declare hist1=
 camera {
    location <5,1,0>
    look_at <0,100,0>
    angle 10
    up z
    right x
    //rotate <0,0,0>
}

camera{hist0}

#declare flame = media {
      emission 0.5
      intervals 20
      density {
            gradient y 
            turbulence <5,.1,1> 
            scale <.02,1,.05>
            color_map {
                [0.0              color rgb <0,    0,              0>]
                [0.25 +clock*0.02 color rgb <0,    0,              0>] // sharp bottomside cutoff
                [0.3  +clock*0.02 color rgb <.8*.7255-clock*.05*.15/.8, 
                                             .8 -clock*.05,   
                                             .1 -clock*.05*.1/.8>]  
                [0.375+clock*0.02 color rgb <.6*.7255 -clock*.05*.1/.6, 
                                            .6  -clock*.05,   
                                            .1  -clock*.05*.1/.6>]
                [0.45 +clock*0.02 color rgb <.15*.7255-clock*.01*.025/.15,
                                             .15 -clock*.01,  
                                             .0>]
                [0.6  +clock*0.02 color rgb <.035, 
                                             0.00294*.035  ,  
                                             0>]
                [0.65             color rgb <.05, 
                                             .00294*.05,
                                             0>]
                [0.85             color rgb <.05, .01,             0>]
                [1.0              color rgb <0,   0,               0>]
                }
          }
        }  
        
#declare static = media {
      emission 0.5
      intervals 20
      density {
            gradient y 
            turbulence <5,.1,1> 
            scale <.02,1,.05>
            color_map {
                [0.0              color rgb <0,    0,   0>]
                [0.25             color rgb <0,    0,   0>] // sharp bottomside cutoff
                [0.3              color rgb <.8*.7255, .8,   .1>]  
                [0.375            color rgb <.6*.7255,.6,  .1>]
                [0.45             color rgb <.15*.7255,.15,   .00>]
                [0.6              color rgb <.035, 0.00294*.035,   0>]
                [0.65             color rgb <.05,0.00294*.05,  0>]
                [0.85             color rgb <.01, .00294,  0>]
                [1.0              color rgb <0,   0,    0>]
                }
          }
        } 


#declare arc1=
 cylinder{ <-.6, 0, -.6>, <-.2,1,-.2>,.1
hollow // won't glow without "hollow" 
}
 
#declare arc2=
 cylinder{ <.6, 0, .6>, <.2,1,.2>,.1
hollow
}
 
#declare arc3=
 cylinder{ <.6, 0, -.6>, <.2,1,-.2>,.1
hollow // won't glow without "hollow" 
}

#declare arc4=
 cylinder{ <-.6, 0, .6>, <-.2,1,.2>,.1
hollow // won't glow without "hollow" 
}
#declare grp_static=
    union{
    object{arc2}
    object{arc3}
    object{arc4}
    }
    
#declare aurora_static=
object {grp_static 
    interior {media{static}}
    pigment { rgbt 1} // "t" transparency
} 

object {aurora_static 
        scale     <1, 50,  1> 
        rotate    0 
        translate <0, 100, 0>
        }
        
#declare aurora_flame=
object{ arc1
    interior {media{flame}}
    pigment { rgbt 1} // "t" transparency
} 

object {aurora_flame 
        scale     <1, 100,  1> 
        rotate    0 
        translate <0, 100, 0>
        }

