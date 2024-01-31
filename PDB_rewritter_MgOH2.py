import numpy as np
import ast as ast

# TIP4P flag (True or False)
TIP4P = True

Filename = '6mole_MgOH2'
Inputfile = Filename + '.pdb'
Outputfile = Filename + '_new.pdb'

# Read pdb file
box = np.genfromtxt(Inputfile, dtype=str, max_rows=1)
pdb = np.genfromtxt(Inputfile, dtype=str, skip_header=1, skip_footer=1)
#print(box)

# Create a new matrix with the same dimensions as the original matrix
box_new = []
pdb_new = []

# Convert the strings to the correct data type
converted_row = []
for item in box:
    try:
        converted_item = ast.literal_eval(item)
    except (ValueError, SyntaxError):
        converted_item = item
    converted_row.append(converted_item)
box_new.append(converted_row)

for row in pdb:
    converted_row = []
    for item in row:
        try:
            converted_item = ast.literal_eval(item)
        except (ValueError, SyntaxError):
            converted_item = item
        converted_row.append(converted_item)
    pdb_new.append(converted_row)

# MW insertion
def MW_insert(i, atom_counter, molecular_counter):     # i is the index of the oxygen atom of the water molecule
     # get the coordinates of the atoms
    OW= [float(pdb_new[i][5]), float(pdb_new[i][6]), float(pdb_new[i][7])]
    HW1= [float(pdb_new[i+1][5]), float(pdb_new[i+1][6]), float(pdb_new[i+1][7])]
    HW2= [float(pdb_new[i+2][5]), float(pdb_new[i+2][6]), float(pdb_new[i+2][7])]

    # standardlize the vector
    vector_x = (HW1[0]-OW[0]) + (HW2[0]-OW[0])
    vector_y = (HW1[1]-OW[1]) + (HW2[1]-OW[1])
    vector_z = (HW1[2]-OW[2]) + (HW2[2]-OW[2])
    Rc = (vector_x**2 + vector_y**2 + vector_z**2)**0.5 # angstrom

    # get the coordinates of the MW (dOM = 0.1577 angstrom)
    MW_x = OW[0] + vector_x/Rc*0.1577
    MW_y = OW[1] + vector_y/Rc*0.1577
    MW_z = OW[2] + vector_z/Rc*0.1577

    new_row = ['ATOM', (atom_counter+3), 'MW', 'SOL', molecular_counter, MW_x, MW_y, MW_z, 1.0, 0.0, ' ']
    pdb_new.insert(i+3,new_row)

# Write pdb file
atom_counter = 1
molecular_counter = 1
i=0
while i < (len(pdb_new)+1):
    pdb_new[i][1] = atom_counter
    
    if (i+2) > len(pdb_new)-1:   # Avoid exceptions when the last line is reached
        pdb_new[i+1][1] = atom_counter+1
        break
    else:
        # Change the residue number of Mg(OH)2
        desired_word_MG = 'MG'
        desired_word_O = 'O'
        desired_word_H = 'H'
        if desired_word_MG in pdb_new[i][2] and desired_word_MG in pdb_new[i+1][2] and desired_word_O in pdb_new[i+2][2] and desired_word_H in pdb_new[i+4][2]:
            for j in range(0, 2):
                pdb_new[i+j][2] = 'MG'
                pdb_new[i+j][3] = 'LIG'
                pdb_new[i+j][4] = molecular_counter
            
            for j in range(0, 2):
                pdb_new[i+2+j][2] = 'O'
                pdb_new[i+2+j][3] = 'LIG'
                pdb_new[i+2+j][4] = molecular_counter

            for j in range(0, 2):
                pdb_new[i+4+j][2] = 'H'
                pdb_new[i+4+j][3] = 'LIG'
                pdb_new[i+4+j][4] = molecular_counter

            for j in range(0, 2):
                pdb_new[i+6+j][2] = 'MG'
                pdb_new[i+6+j][3] = 'LIG'
                pdb_new[i+6+j][4] = molecular_counter
            
            for j in range(0, 6):
                pdb_new[i+8+j][2] = 'O'
                pdb_new[i+8+j][3] = 'LIG'
                pdb_new[i+8+j][4] = molecular_counter

            for j in range(0, 6):
                pdb_new[i+14+j][2] = 'H'
                pdb_new[i+14+j][3] = 'LIG'
                pdb_new[i+14+j][4] = molecular_counter

            molecular_counter += 1
    
        # Change the residue number of CO2
        desired_word_C = 'C'
        desired_word_O = 'O'
        if desired_word_C in pdb_new[i][2] and desired_word_O in pdb_new[i+1][2]:
            pdb_new[i][2] = 'C'
            pdb_new[i][3] = 'CO2'
            pdb_new[i][4] = molecular_counter
            pdb_new[i+1][2] = 'O1'
            pdb_new[i+1][3] = 'CO2'
            pdb_new[i+1][4] = molecular_counter
            pdb_new[i+2][2] = 'O2'
            pdb_new[i+2][3] = 'CO2'
            pdb_new[i+2][4] = molecular_counter
            molecular_counter += 1

        # Change the residue number of CH4
        desired_word_C = 'C'
        desired_word_H = 'H'
        if desired_word_C in pdb_new[i][2] and desired_word_H in pdb_new[i+1][2]:
            pdb_new[i][2] = 'C'
            pdb_new[i][3] = 'CH4'
            pdb_new[i][4] = molecular_counter
            pdb_new[i+1][2] = 'H1'
            pdb_new[i+1][3] = 'CH4'
            pdb_new[i+1][4] = molecular_counter
            pdb_new[i+2][2] = 'H2'
            pdb_new[i+2][3] = 'CH4'
            pdb_new[i+2][4] = molecular_counter
            pdb_new[i+3][2] = 'H3'
            pdb_new[i+3][3] = 'CH4'
            pdb_new[i+3][4] = molecular_counter
            pdb_new[i+4][2] = 'H4'
            pdb_new[i+4][3] = 'CH4'
            pdb_new[i+4][4] = molecular_counter
            molecular_counter += 1    

        
        # Change the residue number of H2O
        desired_word_O = 'O'
        desired_word_H = 'H'
        if desired_word_O in pdb_new[i][2] and desired_word_H in pdb_new[i+1][2] and desired_word_H in pdb_new[i+2][2]:
                # Check if Mg(OH)2 is present
                if desired_word_MG in pdb_new[i-2][2] or desired_word_MG in pdb_new[i-6][2]:
                    pass
                else:
                    pdb_new[i][2] = 'OW'
                    pdb_new[i][3] = 'SOL'
                    pdb_new[i][4] = molecular_counter
                    pdb_new[i+1][2] = 'HW1'
                    pdb_new[i+1][3] = 'SOL'
                    pdb_new[i+1][4] = molecular_counter
                    pdb_new[i+2][2] = 'HW2'
                    pdb_new[i+2][3] = 'SOL'
                    pdb_new[i+2][4] = molecular_counter
                    if TIP4P == True:
                        MW_insert(i, atom_counter, molecular_counter)
                    molecular_counter += 1
    
    # renew the counters
    atom_counter += 1
    i += 1

# Write the new pdb file
with open(Outputfile, 'w') as f:
    f.write("%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%3s\n"% (box_new[0][0],box_new[0][1],box_new[0][2],box_new[0][3],box_new[0][4],box_new[0][5],box_new[0][6],box_new[0][7]))
    for i in range(len(pdb_new)):
        f.write("%s%7d  %-4s%s%6d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n"% (pdb_new[i][0],pdb_new[i][1],pdb_new[i][2],pdb_new[i][3],pdb_new[i][4],pdb_new[i][5],pdb_new[i][6],pdb_new[i][7],pdb_new[i][8],pdb_new[i][9],pdb_new[i][10]))
    f.write("TER\n")


