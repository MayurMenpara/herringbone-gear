----------------------------------------------------------------àª¶à«àª°à«€ àª—àª£à«‡àª¶àª¾àª¯ àª¨àª®àªƒ--------------------------------------------------------------

--> In this programe we will generate parametric helical gear with involute geomatery and simplified fillet radious.
--> The required input is Module, Number of teeth, pressure angle, fillet radious, Profile shift, helix angle, thickness of gear.
--------------------------------Notes-------------------------------------
--> Inputs of Pressure angle is in degree
--> There is Gear and Pinion defined at the ending
--> User need to find Appropriate Relation between Gear and Pinion Parameters
--> In Lua Function start with ("Capital Latter") are globle function you can use it only ones
--> In Lua Function start with ("Small Latter") are Local function.
--> The Credit of function Spiral_extrude goes to Prof. Dr. Ing. Stefan Scherbarth
--> The curves ( Left fillet, left Involute, Right Involute, Right fillet) are in this oder to maintain continuty
--> in this file we have added text with font ('Sansation-Regular.ttf') so plese add tis font if you find any.

--------------------------- Required Function To Generate Helical gear ---------------------

--F = font(Path .. 'Sansation_Regular.ttf')  -- Font for printing text

-----------------------------------------------------------------------------------------------
-------------------------------------------bearing --------------------------------------------
----------------------------------------------------------------------------------------------- 

-----------------------------------------------------------------------------------
---------------------------------------------inputs -------------------------------
-----------------------------------------------------------------------------------


m = 1
Beta_p = ui_number("Helix angle for Both Gear                 [Degree]",44,0,50);
m_t = m*ui_scalar("Module Of Both Gear                       [-]",3,0.5,10);
alpha_t = ui_scalar("Pressure Angle at Pitch circle    [Degree]",20,20,30);
Thickness =m* ui_number("Thickness Of Gear(Gear1)                  [mm]",15,1,50);
z1 = ui_number("Number Of Teeth(Gear1)                    [-]",30,3,50);
z2 = ui_number("Number Of Teeth(Gear2)                    [-]",20,3,50);
x_coef1 = ui_scalar("Profile Shift(Gear1)                      [-]",0,-0.75,1);
x_coef2 = ui_scalar("Profile Shift(Gear2)                      [-]",0,-0.75,1); -- Profile shift of the mashing gear
f_r1 = m*ui_scalar("Fillet Radius(Gear1)                      [mm]",1,0,5);
--f_r2 = m*ui_scalar("Fillet Radius(Gear2)                      [mm]",1,0,5);
-------------------------------------------------------------------------
-------------------------------------Calculation for center distance -----------------------
--------------------------------------------------------------------------

inver_inv = function(inv_Alpha)                                         -- inverse involute function 
    x = inv_Alpha
    Alpha = 3^(1/3)*x^(1/3) - 2/5*x + 9/175*3^(2/3)*x^(5/3) - 2/175*3^(1/3)*x^(7/3) 
            - 144/67375*x^(9/3) + 3258/3128125*3^(2/3)*x^(11/3) - 49711/153278125*3^(1/3)*x^(13/3) 
            - 1130112/9306171875*x^(15/3) + 5169659643/95304506171875*3^(2/3)*x^(17/3)
        return(Alpha)
    end


Involut_curve = function(base_radius, involute_angle)                   -- involute curve generating function
    return v(base_radius *
                 (math.sin(involute_angle) - involute_angle *
                 math.cos(involute_angle)), base_radius *
                 (math.cos(involute_angle) + involute_angle *
                 math.sin(involute_angle)))
end


Angle_involute = function(b, a)                                         -- finds angle of Involute at any given radious
    return (math.sqrt((a * a - b * b) / (b * b)))
end 


Mirror = function(coord)                                                --  for Mirroring of Contour
    return v(-coord.x, coord.y)
end                 


Rotation = function(rotate, coord)                                      -- Rotaion of any Contour
    return v(math.cos(rotate) * coord.x + math.sin(rotate) * coord.y,
             math.cos(rotate) * coord.y - math.sin(rotate) * coord.x,0)
end



Fillet_center = function(coontur,f_r,r_r)                               -- To fine center point of fillet radious
    Slop = (coontur[3].y - coontur[2].y) / (coontur[3].x - coontur[2].x)
    Slop_ang = math.atan(Slop)
    ---parallel point of involue tangent---
    x = coontur[2].x + f_r * math.cos(Slop_ang + math.pi / 2)
    y = coontur[2].y + f_r * math.sin(Slop_ang + math.pi / 2)
    d = (y - Slop * x) / math.sqrt(Slop * Slop+ 1)
    th1 = math.asin(d / (f_r + r_r)) + Slop_ang
    return(v((f_r + r_r) * math.cos(th1), (f_r + r_r) * math.sin(th1)))
end

Circle = function(a, b, r, th)                                           -- Finds any requrired point on Circle or Whole circle
    return (v(a + r * math.cos(th), b + r * math.sin(th)))
end


Slop_Contour = function(coontur)                                         -- slop of involute curve
    return ((coontur[2].y - coontur[1].y) / (coontur[2].x - coontur[1].x))
end


Spiral_extrude = function(Contour, angle, dir_v, scale_v, z_steps)      -- to Create spiral or helical extrude
    -- extrude a Contour to a shape by turning the contour to angle in z_steps
    -- extrude a Contour in a dircetion given by the vector dir_v
    -- extrude a Contour with a scaling factor given by vector scale_v 
    -- Contour: a table of vectors as a closed contour (start point and end point is the same)
    -- angle: roation angle of contour along z_steps in deg
    -- dir_v: vector(x,y,z) direction of extrusion
    -- sacle_v: vector(x,y,z) scaling factors for extrudion
    -- z_steps: number of steps for the shape, mostly z_steps=2 if angle is equal zero 
    local n_cont= #Contour
    local angle= angle/180*math.pi
 
    local Vertex= {}
    for j= 0,z_steps-1 do
       local phi= angle*j/(z_steps-1)
       local dir_vh= dir_v*j/(z_steps-1)
       local scale_vh= (scale_v - v(1,1,1))*(j/(z_steps-1)) + v(1,1,1)
       for i= 1,n_cont-1 do
           Vertex[i+j*n_cont]= v((dir_vh.x + scale_vh.x * (Contour[i].x*math.cos(phi) - Contour[i].y*math.sin(phi))),
                                 (dir_vh.y + scale_vh.y * (Contour[i].x*math.sin(phi) + Contour[i].y*math.cos(phi))),
                                 (dir_vh.z * scale_vh.z))
       end
       table.insert(Vertex,Vertex[1+j*n_cont])
    end
 
    local vertex_sum_1 = v(0,0,0)
    local vertex_sum_m = v(0,0,0)
    for i= 1,n_cont-1 do
       vertex_sum_1= vertex_sum_1 + Vertex[i]
       vertex_sum_m= vertex_sum_m + Vertex[i+n_cont*(z_steps-1)]
    end    
    table.insert(Vertex,vertex_sum_1/(n_cont-1)) --n_cont*m_cont + 1
    table.insert(Vertex,vertex_sum_m/(n_cont-1)) --n_cont*m_cont + 2
 
 
    Tri= {} -- !!! the index on table with Vertex starts with zero !!!
    local k= 1
    for j=0,z_steps-2 do
       for i= 0,n_cont-2 do
          Tri[k]=   v(i, i+1, i+n_cont) + v(1,1,1)*n_cont*j
          Tri[k+1]= v(i+1, i+n_cont+1, i+n_cont) + v(1,1,1)*n_cont*j
          k= k+2
       end
    end
    for i= 0,n_cont-2 do
       Tri[k]= v(i+1,i,n_cont*z_steps)
       k= k+1
    end
    for i= 0,n_cont-2 do
       Tri[k]= v(i+n_cont*(z_steps-1),i+1+n_cont*(z_steps-1),n_cont*z_steps+1)
       k= k+1
    end
    return(polyhedron(Vertex,Tri))
 end

helix_angle = function(thickness,z,m_t,Beta_p)                          -- It's give relation between Bottom Plane and Top plane
    r_p = z * m_t / 2
    beta = Beta_p*math.pi/180
    return(math.tan(beta)*thickness*180/(d_p*math.pi))

end

angle_gear = function(m_t,x_coef,alpha_t,z)                               -- for make the first tooth straight
    alpha_t = alpha_t
   return(((math.pi*m_t/2) + 2*m_t*x_coef*math.tan(alpha_t))/(z*m_t/2) +
                2*math.tan(alpha_t) - 2*alpha_t)
end



---------------------------------end of required function ------------------------

----------------------------------Function to Create one single tooth (Tooth)--------------
--> All Calculation is Done inside this Function.
--> Input of Pressure angle for function is in degree.
--> this involute curve is stariting from TIF Diameter.




tooth = function(m_t,z1,z2,alpha_t,x_coef1,x_coef2,f_r,a)
   
    Tooth = {}
    c = 0.167 * m_t                                              -- clearance of tooth
    ---------------------calculation for required Paeameters---------------------

    d_p = z1 * m_t                                               -- pitch Diamter
    r_p = d_p / 2                                                -- pitch radious

    d_b = d_p * math.cos(alpha_t)                                -- base diamter of gear
    r_b = d_b / 2                                                -- base raidpous

    h = (2.25 + c_coef - (x_coef1 + x_coef2))*m_t
    d_a = d_p + 2*(1 + c_coef - x_coef2 )*m_t                    -- addensum diameter
    r_a = d_a / 2                                                -- addendum radious
                                                 

    d_r = (d_a - 2*h)                                            -- root diamter
    r_r = d_r / 2                                                -- root_radius

    if r_b < r_r + f_r then r_r = r_b - f_r end
    h_r = m_t + c 
    h_a = m_t + c                                               -- root height


    S_0 = m_t*((math.pi/2)+ 2*x_coef1* math.tan(alpha_t));
    inv_a = math.tan(alpha_t) - alpha_t; 
    Q_fp = ((m_t*((math.pi/4)-(math.tan(alpha_t)))-c*math.tan(alpha_t))*(1+math.sin(alpha_t)))/(math.cos(alpha_t))
    d_TIF = math.sqrt(math.pow(d_p * math.sin(alpha_t) - 2 *(h_a - (m_t * x_coef1) - Q_fp *(1 - math.sin(alpha_t)))/(math.sin(alpha_t)), 2) +   d_b * d_b);                                -- True involute diameter
   
    r_TIF = d_TIF / 2

 
    tooth_ang = (((math.pi*m_t/2) + 2*m_t*x_coef1*math.tan(alpha_t))/r_p +
                2*math.tan(alpha_t) - 2*alpha_t)                -- Angle between two involute Curve.

    Start_Involute = Angle_involute(r_b,r_TIF);    
    End_Involute   = Angle_involute(r_b,r_a);                  -- ending of involute curve.

    Points = 10                                              -- points for Better Accuracy.

    -------------------------------------------------------------------------------------------

    -------------------------------Calculation For Fillet Radious------------------------------
    --> for Calculation we need Involute Curve to maintain the continuty of curve

    local involute = {}
   
    for i = 1,Points do
        involute[i] = Involut_curve(r_b,(Start_Involute + (End_Involute - Start_Involute)*i/Points))
    end
    center_f = {}
    center_f[1] = Fillet_center(involute,f_r,r_r)

    Slop_inv = Slop_Contour(involute)

    Slop_ang = math.atan(Slop_inv)


    Start_fillet = 2 * math.pi + math.atan(center_f[1].y / center_f[1].x)
    End_fillet = 3 * math.pi / 2 + Slop_ang
    -------------------------------------------------------------------------------------------

    if r_r ~= 0 then
        for i = 1,Points do                                                           -- Left Fillet
            Tooth[#Tooth + 1] = Circle(center_f[1].x,center_f[1].y, f_r,(Start_fillet + (End_fillet - Start_fillet) *i / Points))
        end
    end 
  
        for i = 1,Points do
            Tooth[#Tooth + 1] = Involut_curve(r_b,(Start_Involute+                    -- Left Involute
                                        (End_Involute - Start_Involute) * i / Points))
        end

        for i = Points,1,-1 do
            Tooth[#Tooth + 1] = Rotation(tooth_ang, Mirror(Involut_curve(r_b,         -- Right involut
                                        (Start_Involute+(End_Involute - Start_Involute) * i / Points))))      
        end


    if r_r ~= 0 then
        for i = Points,1,-1 do                                                        -- Right Fillet
        Tooth[#Tooth + 1] = Rotation(tooth_ang,Mirror(Circle(center_f[1].x,center_f[1].y,f_r,(Start_fillet + (End_fillet - Start_fillet)*i / Points))))
        end
    end 

    --Tooth[#Tooth + 1] = Tooth[1]
    return Tooth
end


gear = function(z,Tooth)           -- Whole gear
    local Gear = {}
    for i = 1,z do
    local  n_cont = #Tooth         -- total points on tooth 
    angle = 2*math.pi*i/z
        for j = 1,n_cont do
            Gear[#Gear + 1] = Rotation(angle,v(Tooth[j].x,Tooth[j].y))
        end
    end
    Gear[#Gear + 1] = Gear[1]
    return Gear
end

alpha_t = alpha_t*math.pi/180  
inv_Alpha = 2*math.tan(alpha_t)*(x_coef1 + x_coef2)/(z1 + z2) + math.tan(alpha_t) - alpha_t
alpha = inver_inv(inv_Alpha) 
c_coef = (z1 + z2)*(math.cos(alpha_t)/math.cos(alpha) - 1)/2
a = ((z1 + z2)/2 + c_coef)*m_t



-----------------------------------------------------------------------------------------------
-------------------------------------------bearing end --------------------------------------------
----------------------------------------------------------------------------------------------- 
---------------------------------------------Same Parameter for both Gear------------------------------------

----------------------------------------------------------
----------------------------Gear 1 -----------------------
----------------------------------------------------------


-------------------------------------------------------------------------------------------------
Tooth1 = tooth(m_t,z1,z2,alpha_t,x_coef1,x_coef2,f_r1)
Gear = gear(z1,Tooth1)
helix = helix_angle(Thickness,z1,m_t,Beta_p)
angle = angle_gear(m_t,x_coef1,alpha_t,z1) 
r= rotate(0,0, angle*90/math.pi)

---------------------- GEar details ----------------
o = 1
s = scale(o,o,o)



Helical_gear = Spiral_extrude(Gear, helix, v(0,0,-Thickness), v(1,1,1), 20)
Helical_gear2 = difference(r*Helical_gear, translate(0,0,-10)*ccylinder(8*m, Thickness*2))

emit(s*Helical_gear2,6)


-- ----------------------------------------------------------
-- ----------------------------Gear 2 -----------------------
-- ----------------------------------------------------------

Tooth1 = tooth(m_t,z1,z2,alpha_t,x_coef1,x_coef2,f_r1)
Gear = gear(z1,Tooth1)
helix = -helix_angle(Thickness,z1,m_t,Beta_p)
angle = angle_gear(m_t,x_coef1,alpha_t,z1) 
r= rotate(0,0, angle*90/math.pi)

---------------------- GEar details ----------------

r3 = rotate(0,180, 0)


Helical_gear = Spiral_extrude(Gear, helix, v(0,0,-Thickness), v(1,1,1), 20)
Helical_gear2 = difference(r*Helical_gear, translate(0,0,-10)*ccylinder(8*m, Thickness*2))

emit(r3*Helical_gear2,6)


-- -------------------------------------------------------------------------------------------------
-- Tooth1 = tooth(m_t,z2,z1,alpha_t,x_coef2,x_coef1,f_r2)
-- Gear = gear(z2,Tooth1)
-- helix = -helix_angle(Thickness,z2,m_t,Beta_p)
-- angle = angle_gear(m_t,x_coef2,alpha_t,z2) 
-- r1= rotate(0,0,- 180 - (360/z2 - angle*180/math.pi)/2)
-- r3 = rotate(0,0, -i*z1/z2)


-- ---------------------- geaar detail ----------------


-- Helical_gear = Spiral_extrude(Gear, helix, v(0,0,-Thickness), v(1,1,1), 100)
-- Helical_gear2 = r3*difference(r1*Helical_gear, translate(0,0,-10)*ccylinder(8, Thickness*2))

-- emit(translate(0,a,Thickness/2)*Helical_gear2,2)






-- --------base plate -----------

-- emit(translate(0,a/2,-Thickness - 0.1/2)*cube(2*a,3*a,0.1),2)


-- ---Shaft---

-- emit(translate(0,0,Thickness/2)*ccylinder(8,Thickness*3),0)
-- emit(translate(0,a,0)*ccylinder(8,Thickness*2),0)

-- ---------frem ------------

-- frem1 = (translate(0,0,-Thickness - 0.1/2)*cube(2*a,30,2*Thickness))
-- frem2 = (translate(0,0,-Thickness - 0.1/2)*cube(2*a-6,30,2*Thickness-6))
-- emit(difference(difference(frem1, frem2),ccylinder(14,Thickness*2)))
-- emit(translate(0,a,0)*difference(difference(frem1, frem2),ccylinder(14,Thickness*2)))

-- -----
-- Circle = function(a, b, r)  
--     n = 200
--     xy = {}
--     for i = 0,n do
--         xy[i] = v(a + r * math.cos(2*math.pi*i/n), b + r * math.sin(2*math.pi*i/n))
--     end
--                            -- Finds any requrired point on Circle or Whole circle
--     return xy
-- end
-- ---------inner body-----------

-- shaft_radious = 8
-- ball_radious = 2


-- cir = Circle(0, 0, shaft_radious + 2)
-- dir    = v(0,0,ball_radious*3)
-- Inner_1 = (linear_extrude(dir, cir) )

-- cir = Circle(0, 0, shaft_radious)
-- Inner_2 = (linear_extrude(dir, cir) )

-- inner = (difference(Inner_1, Inner_2))


-- ----------for bearing ----------------
-- cir = Circle(ball_radious + shaft_radious + 1, 0, ball_radious)
-- bearing_hole  = ( translate(0,0,3)*rotate_extrude( cir, 25 ) )

-- -----------outer body--------
-- cir = Circle(0, 0, shaft_radious + 3*ball_radious)

-- outer_1 = (linear_extrude(dir, cir)  )

-- cir = Circle(0, 0, shaft_radious + 2*ball_radious)
-- outer_2 = (linear_extrude(dir, cir)  )

-- outer = (difference(outer_1, outer_2))


-- body = (difference(union(outer,inner),bearing_hole))

-- n = 8
-- ball = (translate(0,11,3)*sphere(2))

-- for i = 1,n do
--     th = 360*i/n
--     ball = union{ ball,rotate(0,0,th)*ball}

-- end

-- ball = (ball)

-- emit(translate(0,a,Thickness-6)*r3*union(body,ball),6)
-- emit(r2*translate(0,0,Thickness-6)*union(body,ball),6)

-- --------------upper wheel ------------
-- cir = Circle(25, 0, ball_radious)
-- emit( translate(0,0,2*Thickness-2)*rotate_extrude( cir, 50 ) )

-- rim1 = cube(2,50, 2)
-- rim2 = rotate(0,0,120)*rim1
-- rim3  = rotate(0,0,120)*rim2
-- rim = union{rim1,rim2,rim3}
-- emit(r2*translate(0,0,2*Thickness-3)*rim,2)

