import copy
import os
import sys
import math
import time
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# This code assumes that the porin is aligned to the y-axis!


# Generic formula for ellipse
def myFunction(a, b, c, d, e, x, y):
    z = a + b*x + c*x*x + d*y + e*x*y - y*y
    return z

'''
def pore_area(PDB, VDW_radius, AngleSearchSize, YSearchSize, outputPrefix, YMin=-2, YMax=2, PDB_chain = "", User="mfa5147", exportScatter = True, exportPYMOL = True):
    # Inputs
    # PDB: PDB File Name
    # PDB_chain: Only necessary if multiple chains in "PDB"
    # Origin: Approximate center of ellipse
    # outputName: All output file names minus file extention
    # VDW_radius: Value from I-CAVER
    # AngleSearchSize: 360 degrees divided by this number number of searches
    # YSearchSize: number of steps between YMin and YMax
    # User: User Name
    # exportScatter: True if export scatter plot
    # exportPYMOL: True if export PYMOL figure

    # Code
    print("PORE AREA:", PDB)
    # Create the experiment dictionary
    experiment = {}
    experiment["User"] = User

    # Create the Molecule class object
    porin_models = list(PDBParser().get_structure(outputPrefix, PDB).get_models())
    if len(porin_models) == 1:
        porin_chains = porin_models[0].get_chains()
    else:
        raise ValueError("structure has {} (not 1) models".format(len(porin_models)))
    if PDB_chain != "":
        porin = porin_chains[PDB_chain]
    else:
        porin = list(porin_chains)[0]

    # Bin all of the atoms by angle and y-coordinate
    angleMin, angleMax = -180, 180
    angle_step = 360 / AngleSearchSize
    y_step = (YMax - YMin) / (YSearchSize - 1)
    bins = {(y, angle): list() for y in range(YSearchSize) for angle in range(AngleSearchSize)}
    for res in porin:
        for atom in res:
            x, y, z = atom.get_coord()
            angle = np.rad2deg(math.atan2(x, z))
            y_bin = int(round((y - YMin) / y_step))
            angle_bin = int(round((angle - angleMin) / angle_step))
            bin_ = y_bin, angle_bin
            if bin_ in bins:
                bins[bin_].append(atom)

    # Parameterize the Molecule
    VDW = {"C": 1.7, "H": 1.2, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.8, "S": 1.8}
    max_Radius = 1.8
    # Check pore size at each y coordinate.
    for y_bin in range(YSearchSize):
        y_coord = round(y_bin * y_step + YMin, 2)
        print("Y:", y_coord)
        # Store the results in these lists 
        Edge_Atoms = []
        Edge_Coordinates = []
        # Search through the angle divisions.
        for angle_bin in range(AngleSearchSize):
            bin_ = y_bin, angle_bin
            if len(bins[bin_]) == 0:
                continue
            theta = np.deg2rad(angle)
            # Get the distance between each atom and the y axis.
            atom_distances = [np.linalg.norm(atom.get_coord()[0: 3: 2]) for atom in bins[bin_]]
            # Find the atom closest to the y axis..
            closest_atom = bins[bin_][np.argmin(atom_distances)]
            min_Atom = str(closest_atom.get_parent().get_id()[1]) + "_" + closest_atom.get_name()
            Edge_Atoms.append(min_Atom)
            Edge_Coordinates.append(list(closest_atom.get_coord()))
        twoD = []
        x = []
        y = []
        for coor in Edge_Coordinates:
            twoD.append([coor[0], coor[2]])
            x.append(coor[0])
            y.append(coor[2])
        # Remove any outliers
        mean =  np.mean(x)
        std = np.std(x)
        tot_rem = []
        for i in range(len(Edge_Coordinates)):
            if x[i] < (mean-2*std):
                tot_rem.append(i)
            if x[i] > (mean+2*std):
                tot_rem.append(i)
        new2D = []
        newX = []
        newY = []
        for i in range(len(Edge_Coordinates)):
            if i not in tot_rem:
                new2D.append(twoD[i])
                newX.append(x[i])
                newY.append(y[i])
        twoD = new2D
        x = newX
        y = newY
        try:
            area, r2, model, coeffs = regression(twoD)
        except ValueError as e:
            print(e)
        text =  "Y:      {}\n".format(y_coord)
        text += "Area:   " + format(area, '.3f') + "\n"
        text += "R^2:    " + format(r2, '.3f') + "\n"
        text += "Removed " + str(len(tot_rem)) + " Atoms from Model\n"
        text += "Model:  " + model + "\n\nEdge Atoms:\n"
        for atom in Edge_Atoms:
            items = atom.split("_")
            res = "Position: " + items[0]
            ATOM = "Atom: " + items[1]
            text += res.ljust(17) + ATOM + "\n"
        text += "\n\nEdge Coordinates:\n"
        for coor in Edge_Coordinates:
            text += "(" + format(coor[0], '.3f') + ", " + format(coor[1], '.3f')
            text += ", " + format(coor[2], '.3f') + ")\n"
        outputName = "{}_{}".format(outputPrefix,  y_coord)
        fileName = outputName + ".txt"
        with open(fileName, "w") as file:
            file.write(text)
        if exportScatter:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(x, y, c='k')
            DX = max(x)-min(x)
            DY = max(y)-min(y)
    
            xlist = np.linspace(min(x)-.25*DX,max(x)+.25*DX,500)
            ylist = np.linspace(min(y)-.25*DY,max(y)+.25*DY,500)
    
            X,Y = np.meshgrid(xlist,ylist)
    
            Z = myFunction(coeffs[4], coeffs[1], coeffs[0], coeffs[3], coeffs[2], \
            X,Y)
            CS = ax.contour(X,Y,Z, [0])
            ax.set_xlim([min(x)-.25*DX, max(x)+.25*DX])
            ax.set_ylim([min(y)-.25*DY, max(y)+.25*DY]) 
            ax.annotate(r'$\mathregular{R^2=\/' + format(r2, '.3f') + '}$', \
            (max(x)-.17*DX, max(y)+.08*DY))
            ax.annotate(r'$\mathregular{Area\/=\/' + format(area, '.1f') + \
            '\/\AA^2}$', (max(x)-.17*DX, max(y)+.15*DY))
            ax.set_xlabel(r'$\mathregular{x\/(\AA)}$')
            ax.set_ylabel(r'$\mathregular{z\/(\AA)}$')
            ax.set_title('$\mathregular{' + outputName + '}$')
            fig.savefig(outputName + '_scatter.png')
            plt.close('all')
        if exportPYMOL:
            script = "cd ~/Downloads/OmpFMutants\nload " + PDB + "\ncreate obj1, i. "
            includeMe = []
            for atom in Edge_Atoms:
                pos = atom.split("_")[0]
                if pos not in includeMe:
                    includeMe.append(pos)
            for i in includeMe:
                script += i + "+"
            script = script[:-1] + """
util.cbay obj1
create obj2, """ + PDB.split(".")[0] + """ and not obj1
color green, obj2
delete """ + PDB.split(".")[0] + """
create obj4, obj1
show_as cartoon, obj2
show_as sticks, obj1
show_as surface, obj4
set transparency, 0.5, obj4 or obj2
color yellow, obj4
set_view (\
     0.992484033,    0.106027707,   -0.061102211,\
     0.078066587,   -0.164057046,    0.983357430,\
     0.094238915,   -0.980736375,   -0.171100885,\
     0.000000000,    0.000000000, -176.344696045,\
     1.304159164,    0.480228424,    6.366327286,\
   139.031494141,  213.657897949,  -20.000000000 )
zoom all
set ray_opaque_background, off
set ray_shadows, off
"""
        script += "png " + outputName + "_pyMOL.png, height=800, width=800, dpi=300, ray=1"
        with open(outputName + "_pymol_script.pml", "w") as file:
            file.write(script)
'''

def pore_area(PDB, VDW_radius, AngleSearchSize, RadiusSearchSize, YSearchSize, outputPrefix,
PDB_chain = "", Origin = [0.0, 0.0, 0.0], User = "Ratul Chowdhury",
exportScatter = True, exportPYMOL = True):
    # Inputs
    # PDB: PDB File Name
    # PDB_chain: Only necessary if multiple chains in "PDB"
    # Origin: Approximate center of ellipse
    # outputName: All output file names minus file extention
    # VDW_radius: Value from I-CAVER
    # AngleSearchSize: 360 degrees divided by this number number of searches
    # RadiusSearchSize: % increase of VDW radius after each step
    # User: User Name
    # exportScatter: True if export scatter plot
    # exportPYMOL: True if export PYMOL figure

    # Code
    print("PORE:", PDB)
    tm = time.time()
    
    # Create the experiment dictionary
    experiment = {}
    experiment["User"] = User

    # Create the Molecule class object
    porin_models = list(PDBParser().get_structure(outputPrefix, PDB).get_models())
    if len(porin_models) == 1:
        porin_chains = porin_models[0].get_chains()
    else:
        raise ValueError("structure has {} (not 1) models".format(len(porin_models)))
    if PDB_chain != "":
        porin = porin_chains[PDB_chain]
    else:
        porin = list(porin_chains)[0]
    YMin = -2
    YMax = 2
    y_step = (YMax - YMin) / (YSearchSize - 1)
    for y_bin in range(YSearchSize):
        y_coord = round(y_bin * y_step + YMin, 2)
        print("Y:", y_coord)
        Origin[1] = y_coord
        # Parameterize the Molecule
        VDW = {"C": 1.7, "H": 1.2, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.8, \
        "S": 1.8}
        max_Radius = 1.8
        # Store the results in these lists   
        Edge_Atoms = []
        Edge_Coordinates = []
        # Generate the search angle
        for angle in range(AngleSearchSize):
            theta = (angle*math.pi*2/AngleSearchSize)
            # Generate the search radius
            r = float(VDW_radius)
            Found = False
            # Continuously search until a clash is found
            while not Found:
                # Generate the new coordinates
                x = r*math.sin(theta)
                y = Origin[1]
                z = r*math.cos(theta)
                new = np.array([x, y, z])
                # Go through all of the atoms
                Min = 1000.0
                min_Atom = "NONE"
                for res in porin:
                    for atom in res:
                        diff = np.linalg.norm(atom.get_coord() - new) - VDW[atom.get_name()[0]] - \
                        VDW["H"]
                        if diff < Min:
                            Min = diff
                            min_Atom = atom.get_parent().get_resname() + "_" + atom.get_name()
                if Min < 0.0:
                    Edge_Atoms.append(min_Atom)
                    Edge_Coordinates.append([round(r*math.sin(theta),2), round(Origin[1],2),
                    round(r*math.cos(theta),2)])
                    Found = True
                    break
                r += (float(RadiusSearchSize)/100.0)*float(VDW_radius)
        
        twoD = []
        x = []
        y = []
        for coor in Edge_Coordinates:
            twoD.append([coor[0], coor[2]])
            x.append(coor[0])
            y.append(coor[2])
        # Remove any outliers
        mean =  np.mean(x)
        std = np.std(x)
        tot_rem = []
        for i in range(len(Edge_Coordinates)):
            if x[i] < (mean-2*std):
                tot_rem.append(i)
            if x[i] > (mean+2*std):
                tot_rem.append(i)
        new2D = []
        newX = []
        newY = []
        for i in range(len(Edge_Coordinates)):
            if i not in tot_rem:
                new2D.append(twoD[i])
                newX.append(x[i])
                newY.append(y[i])
        twoD = new2D
        x = newX
        y = newY
        try:
            area, r2, model, coeffs = regression(twoD)
        except ValueError as e:
            print(e)
            continue
        
        outputName = "{}_{}".format(outputPrefix,  y_coord)
        text = "Area:   " + format(area, '.3f') + "\n"
        text += "R^2:    " + format(r2, '.3f') + "\n"
        text += "Removed " + str(len(tot_rem)) + " Atoms from Model\n"
        text += "Model:  " + model + "\n\nEdge Atoms:\n"
        for atom in Edge_Atoms:
            items = atom.split("_")
            res = "Position: " + items[0]
            ATOM = "Atom: " + items[1]
            text += res.ljust(17) + ATOM + "\n"
        text += "\n\nEdge Coordinates:\n"
        for coor in Edge_Coordinates:
            text += "(" + format(coor[0], '.3f') + ", " + format(coor[1], '.3f')
            text += ", " + format(coor[2], '.3f') + ")\n"
        fileName = outputName + ".txt"
        with open(fileName, "w") as file:
            file.write(text)
        if exportScatter:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(x, y, c='k')
            DX = max(x)-min(x)
            DY = max(y)-min(y)
    
            xlist = np.linspace(min(x)-.25*DX,max(x)+.25*DX,500)
            ylist = np.linspace(min(y)-.25*DY,max(y)+.25*DY,500)
    
            X,Y = np.meshgrid(xlist,ylist)
    
            Z = myFunction(coeffs[4], coeffs[1], coeffs[0], coeffs[3], coeffs[2], \
            X,Y)
            CS = ax.contour(X,Y,Z, [0])
            ax.set_xlim([min(x)-.25*DX, max(x)+.25*DX])
            ax.set_ylim([min(y)-.25*DY, max(y)+.25*DY]) 
            ax.annotate(r'$\mathregular{R^2=\/' + format(r2, '.3f') + '}$', \
            (max(x)-.17*DX, max(y)+.08*DY))
            ax.annotate(r'$\mathregular{Area\/=\/' + format(area, '.1f') + \
            '\/\AA^2}$', (max(x)-.17*DX, max(y)+.15*DY))
            ax.set_xlabel(r'$\mathregular{x\/(\AA)}$')
            ax.set_ylabel(r'$\mathregular{z\/(\AA)}$')
            ax.set_title('$\mathregular{' + outputName + '}$')
            fig.savefig(outputName + '_scatter.png')
            plt.close('all')
        if exportPYMOL:
            script = "cd ~/Downloads/OmpFMutants\nload " + PDB + "\ncreate obj1, i. "
            includeMe = []
            for atom in Edge_Atoms:
                pos = atom.split("_")[0]
                if pos not in includeMe:
                    includeMe.append(pos)
            for i in includeMe:
                script += i + "+"
            script = script[:-1] + """
util.cbay obj1
create obj2, """ + PDB.split(".")[0] + """ and not obj1
color green, obj2
delete """ + PDB.split(".")[0] + """
create obj4, obj1
show_as cartoon, obj2
show_as sticks, obj1
show_as surface, obj4
set transparency, 0.5, obj4 or obj2
color yellow, obj4
set_view (\
     0.992484033,    0.106027707,   -0.061102211,\
     0.078066587,   -0.164057046,    0.983357430,\
     0.094238915,   -0.980736375,   -0.171100885,\
     0.000000000,    0.000000000, -176.344696045,\
     1.304159164,    0.480228424,    6.366327286,\
   139.031494141,  213.657897949,  -20.000000000 )
zoom all
set ray_opaque_background, off
set ray_shadows, off
"""
        script += "png " + outputName + "_pyMOL.png, height=800, width=800, dpi=300, ray=1"
        with open(outputName + "_pymol_script.txt", "w") as file:
            file.write(script)


def regression(points):
    xs = []
    ys = []
    combo = []
    for point in points:
        xs.append(float(point[0]))
        ys.append(float(point[1]))
        combo.append(2.0*float(point[1])*float(point[0]))

    C2 = np.array(xs)
    C1 = np.power(C2, 2)
    C3 = np.array(combo)
    C4 = np.array(ys)
    C5 = np.power(C4, 2)
    C6 = np.ones(C2.size)
    F = np.zeros((C2.size,5))
    F[:,0] = C1
    F[:,1] = C2
    F[:,2] = C3
    F[:,3] = C4
    F[:,4] = C6
    Y = np.ones((C2.size,1))
    Y[:,0] = C5
    M = np.mat(F.T) * np.mat(F)
    B = np.mat(F.T) * np.mat(Y)
    sol = np.linalg.solve(M, B)
    sol0 = -1.0*float(sol.item(0))
    sol1 = -1.0*float(sol.item(1))
    sol2 = -1.0*float(sol.item(2))
    sol3 = -1.0*float(sol.item(3))
    sol4 = -1.0*float(sol.item(4))
    """
    # Transformation of ellipse into "normal" form
    A = np.linalg.inv(np.matrix([[sol0, sol2],[sol2,1]]))
    b = np.matrix([[-1.0*sol1], [-1.0*sol3]])
    transformation = 0.5*A*b
    trans1 = float(transformation.item(0))
    trans2 = float(transformation.item(1))
    A = np.matrix([[sol0, sol2],[sol2,1]])
    c = sol0 * trans1 * trans1 + trans2 * trans2 - 2.0 * sol2 * trans1 * trans2
    c = c + sol4
    # Calculate the eigenvalues and eigenvectors
    val, vec = np.linalg.eig(A)
    # Divide the eigenvalues by the constant
    val1 = float(val.item(0)) / c
    val2 = float(val.item(1)) / c
    # Solve for a and b
    if val1 < 0:
        val1 = val1 * -1.0
    if val2 < 0:
        val2 = val2 * -1.0
    a = math.sqrt(1.0/val1)
    b = math.sqrt(1.0/val2)
    # Calculate the area of the ellipse
    """ 
    theta = 0.5*math.atan(2.0*sol2/(sol0-1.0))
    A = sol0*math.pow(math.cos(theta),2)
    A += 2.0*sol2*math.sin(theta)*math.cos(theta)
    A += math.pow(math.sin(theta),2)
    B = 0
    C = sol0*math.pow(math.sin(theta),2)
    C -= 2.0*sol2*math.sin(theta)*math.cos(theta)
    C += math.pow(math.cos(theta),2)
    D = sol1*math.cos(theta) + sol3*math.sin(theta)
    E = -sol1*math.sin(theta) + sol3*math.cos(theta)
    F = sol4
    a = math.sqrt((-4*F*A*C+C*D*D+A*E*E)/(4*A*C*C))
    b = math.sqrt((-4*F*A*C+C*D*D+A*E*E)/(4*A*A*C))
    print (a, b)
    area = math.pi * a * b

    # Report Error of Equation
    # Total number of values
    n = float(len(points))
    # Mean
    sum = 0.0
    for point in points:
        sum += math.pow(float(point[1]),2)
    mean = sum/n
    # SSTot
    sstot = 0.0
    for point in points:
        sstot += math.pow(math.pow(float(point[1]),2) - mean, 2)
    # SSRes
    ssres = 0.0
    for point in points:
        X = float(point[0])
        Y = float(point[1])
        yModel = 1.0*sol0*math.pow(X, 2) + 2.0*sol2*X*Y + 1.0*sol1*X
        yModel = yModel + 1.0*sol3*Y + 1.0*sol4 + Y*Y
        ssres += math.pow(yModel,2)
    # Calculate the R^2 value for the model
    r2 = 1.0 - ssres/sstot
    # Store the equation for the model
    model = "y^2 = " + format(-sol0, '.3f') + "x^2"
    if sol1 > 0:
        model += " - " + format(sol1, '.3f') + "x"
    else:
        model += " + " + format(-sol1, '.3f') + "x"
    if sol2 > 0:
        model += " - " + format(2*sol2, '.3f') + "xy"
    else:
        model += " + " + format(-2*sol2, '.3f') + "xy"
    if sol3 > 0:
        model += " - " + format(sol3, '.3f') + "y"
    else:
        model += " + " + format(-sol3, '.3f') + "y"
    if sol4 > 0:
        model += " - " + format(sol4, '.3f')
    else:
        model += " + " + format(-sol4, '.3f')
    coeffs = [-sol0, -sol1, -sol2, -sol3, -sol4]
    return area, r2, model, coeffs


try:
    results_dir, pdb_name, report_dir = sys.argv[1:]
except ValueError:
    raise IOError("usage: python analyze_pores.py /path/to/results/directory nameOfStructure.pdb /path/to/report/directory")

vdw_radius = 0.05
angle_nums = 20
radius_nums = 20
y_nums = 5

if not os.path.isdir(report_dir):
    os.mkdir(report_dir)

for result in os.listdir(results_dir):
    pdb_file = os.path.join(results_dir, result, pdb_name)
    pore_area(pdb_file, vdw_radius, angle_nums, radius_nums, y_nums, os.path.join(report_dir, result))