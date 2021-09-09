set ray_opaque_background, 1
set ray_shadows, 0
set valence, 0
bg_color white

load PDB/2aep_fit.pdb
symexp sym, 2aep_fit, 2aep_fit, 4
alter sym01000000 and chain A, chain='B'
alter sym02000000 and chain A, chain='C'
alter sym03000000 and chain A, chain='D'
create NA_protomer_1, 2aep_fit and chain A and resi 82-465
create NA_protomer_2, sym01000000 and chain B and resi 82-465
create NA_protomer_3, sym02000000 and chain C and resi 82-465
create NA_protomer_4, sym03000000 and chain D and resi 82-465

hide all
show surface, NA_protomer_1
show surface, NA_protomer_2
show surface, NA_protomer_3
show surface, NA_protomer_4
spectrum b, blue_white_red, minimum=55, maximum=100
color grey50, NA_protomer_2 
color grey50, NA_protomer_3 
color grey50, NA_protomer_4

set_view (\
    -0.802609622,    0.453452110,   -0.387552798,\
     0.595911503,    0.580507159,   -0.554889023,\
    -0.026639421,   -0.676305592,   -0.736138940,\
    -0.000000000,    0.000000000, -314.817321777,\
   155.332717896,   77.554649353,   46.808834076,\
   248.204376221,  381.430267334,  -20.000000000 )
ray; png graph/Mos99_fit_structure.png
