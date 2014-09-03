# add molecule
mol delete all
mol new
set molid [molinfo top get id]

# set fname "pore_W5_1_F0_20/2.run/2.out.xyz"
set fname "pore_W10_2_F0_40/4.run/4.out.xyz"

mol addfile $fname waitfor all
set fluid [atomselect top "name N"]
set wall  [atomselect top "name C"]

# set fname "0.out.dump"
# mol addfile $fname type lammpstrj waitfor all
# set fluid [atomselect top "type 1"]
# set wall  [atomselect top "type 2"]


# create representation
mol delrep 0 top

mol representation VDW 0.2 20.0
mol material Opaque
mol color colorID 23
mol selection "[$fluid text]"
mol addrep top

mol representation VDW 0.2 20.0
mol material Opaque
mol color colorID 2
mol selection "[$wall text]"
mol addrep top


# set display
axes location lowerleft
display resetview
display height 5.0
display projection orthographic
display antialias on
display shadows on
display update

# rotate x by -90
# rotate y by -30
# rotate x by  30
# rotate y by -20
# rotate z by -10
color Display Background white
scale by 1.0
animate style loop
animate forward
