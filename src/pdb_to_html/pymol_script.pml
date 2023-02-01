# Pymol script used by Python script generating HTML websites.
# It displays the structure with coloring of selected residues
#    and calls pymol2glmol script
# Parameters:
# 1. Path to folder with .pdb files
# 2. Name of protein
# 3. Number of X.Y strings where X is a position of modification and Y its color in #rrggbb format

from sys import argv
path = argv[1:][0]
name = argv[1:][1]

reinitialize
import pymol2glmol

cmd.load(path + name + ".pdb")
cmd.color("gray")
cmd.show(representation="ribbon" )

python

# se tof colors seen so far
colors = set()
for p in argv[1:][2:]:
    pos, color = p.split('.')
    color = color[1:]
    if color not in colors:
       cmd.set_color(color, [int(color[0:2],16), int(color[2:4],16), int(color[4:6],16)])
       colors.add(color)
    cmd.color(color ,"resi " + str(pos))
python end
cmd.do("pymol2glmol " +  name)
