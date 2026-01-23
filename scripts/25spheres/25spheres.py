import sac_de_billes as sdb
import merope

#NB merope is dimensioneless in principle, 
# but we choose to refer everything w.r.t. 
# µmeters for numerical precision 
dimensions_side:int = 8  # µm

# create the box for the RVE
L = [dimensions_side, dimensions_side, dimensions_side]

# set the pores features
radius = 0.85  # µm
pore_fraction = 0.11  # /
distMin = 0.3 # minimal distance between two seeds in the RSA algorithm

# you can change this to modify the random placement of the seeds
randomSeed = 1

# throw the spheres
theSpheres = sdb.throwSpheres_3D(sdb.TypeAlgo.RSA, sdb.NameShape.Tore, L, randomSeed, [
    [radius, pore_fraction]], [2], distMin)

# Now we have to collect the name of the entities we do not want to mesh.
# nb in merope every entity is a "phase"

tab_phase_not_to_mesh = []
for j in range(0, len(theSpheres)):
    sphere = theSpheres[j]
    sphere.phase = 2+j
    tab_phase_not_to_mesh.append(sphere.phase)

# we collect some information about the spheres 
# which will be needed for the mmm simulation
with open("bubbles.txt", 'w') as bbl:
    for sphere in theSpheres:
        bbl.write(str(sphere.phase) + " ")
        bbl.write(str(sphere.center[0]) + " " +
                  str(sphere.center[1]) + " " +
                  str(sphere.center[2]) + " ")
        bbl.write(str(sphere.radius) + "\n")


# we build now the µstructure
sphInc = merope.SphereInclusions_3D()
sphInc.setLength(L)
sphInc.setSpheres(theSpheres)

mi = merope.MultiInclusions_3D()
mi.setInclusions(sphInc)
mi.setMatrixPhase(1) #label the physical volume as 1 in the mesh, useful for gmsh

# mesh options and features
meshGenerator = merope.mesh.MeshGenerator()
meshGenerator.setMeshOrder(2)
meshGenerator.setMeshSize(0.25)
meshGenerator.setMultiInclusions(mi)
#crucial! in out case, the pores aren't meshed, we need to tell it to merope
meshGenerator.do_not_mesh(tab_phase_not_to_mesh)
test_name=str(dimensions_side)+"um_"+str(pore_fraction*100).split(".")[0]+"_pct"
meshGenerator.set_nameOutput([test_name+".vtk"])
meshGenerator.write("mesh_"+test_name+".geo") #generate the geometry file
