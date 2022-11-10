set cartoon_transparency, 0.5
set ray_opaque_background, 1
set ray_shadows, 0
set valence, 0
bg_color white

load PDB/Mos99_WT_NA_tetramer_fit.pdb, Mos99_NA
create NA_protomer_1, Mos99_NA and chain A and resi 82-465
create NA_protomer_2, Mos99_NA and chain B and resi 82-465
create NA_protomer_3, Mos99_NA and chain C and resi 82-465
create NA_protomer_4, Mos99_NA and chain D and resi 82-465

hide all
show surface, NA_protomer_1
show cartoon, NA_protomer_2
show cartoon, NA_protomer_3
show cartoon, NA_protomer_4
spectrum b, blue white, minimum=50, maximum=90
color grey50, NA_protomer_2 
color grey50, NA_protomer_3 
color grey50, NA_protomer_4

set_view (\
     0.485448331,   -0.073379867,    0.871174872,\
    -0.873858631,   -0.070882082,    0.480972916,\
     0.026454758,   -0.994769752,   -0.098537005,\
    -0.000178725,    0.000043211, -308.965179443,\
    -1.907482624,   -2.206866741,  -43.625488281,\
   252.302276611,  365.490966797,  -20.000000000 )
hide cartoon, NA_protomer_3
ray; png graph/Mos99_fit_structure_interface.png

set_view (\
    -0.934556365,    0.353700727,    0.038601313,\
     0.351653963,    0.934727311,   -0.051123798,\
    -0.054160401,   -0.034205634,   -0.997934103,\
    -0.000178725,    0.000043211, -341.060668945,\
    -1.907482624,   -2.206866741,  -43.625488281,\
   284.397705078,  397.586486816,  -20.000000000 )
show cartoon, NA_protomer_3
ray; png graph/Mos99_fit_structure_top.png

set_view (\
    -0.968248069,   -0.014799706,   -0.249534056,\
     0.249928832,   -0.075886443,   -0.965276778,\
    -0.004646836,   -0.996992946,    0.077176943,\
    -0.000218464,   -0.001332954, -308.965179443,\
    -3.466495037,   -2.194570780,  -48.400943756,\
   252.302276611,  365.490966797,  -20.000000000 )
ray; png graph/Mos99_fit_structure_side1.png

set_view (\
     0.628273726,    0.034595944,   -0.777216196,\
     0.776401460,    0.035761908,    0.629207313,\
     0.049560130,   -0.998746693,   -0.004389367,\
    -0.000218464,   -0.001332954, -308.965179443,\
    -3.466495037,   -2.194570780,  -48.400943756,\
   252.302276611,  365.490966797,  -20.000000000 )
ray; png graph/Mos99_fit_structure_side2.png

set_view (\
     0.932319820,    0.325315684,    0.157923386,\
    -0.336334646,    0.940498292,    0.048206195,\
    -0.132847831,   -0.098060042,    0.986258626,\
    -0.000025246,   -0.000044078, -308.941741943,\
    -0.617381573,   -1.537739754,  -48.865543365,\
   252.302276611,  365.490966797,  -20.000000000 )
ray; png graph/Mos99_fit_structure_bottom.png
