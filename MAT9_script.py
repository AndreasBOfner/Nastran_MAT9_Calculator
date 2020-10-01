# NASTRAN MAT9 Stiffness matrix calculator / material card export tool for e-machine modeling
#
# 	V08: 2019-07-08
#	auth: Andreas Benjamin Ofner
	

# THIS PROGRAM is designed to calculate the stiffness matrix and export it as a MAT9 material card 
#			   for an e-machine stator to be used with NASTRAN solver.
# 	   
#		INPUTS:  Young's Modulus, Density of used materials
#				 Stack factor, Fill factor in your model 
#				 Filename for the Nastran MAT9 card
#
#		OUTPUTS: Calculated Stiffness Matrix for e-machine stator
#                BDF-File in Nastran MAT9 format                


# ----------------------------------------- INPUT SECTION -------------------------------------------


# ESSENTIAL PARAMETERS - basic material characteristics which characterize properties of laminates 

inp_ESteel = 210000 	# [MPa] Young's Modulus of Steel for Stator Lamination
inp_ECopper = 70000		# [MPa] Young's Modulus of Copper in Coils 
inp_EEpoxy01 = 2500 	# [MPa] Young's Modulus of Epoxy used in Stator Lamination
inp_EEpoxy02 = 1000		# [MPa] Young's Modulus of Epoxy used to fill Coil airgaps
inp_stack = 0.85		# [-]   Stack Factor - relative amount of Steel mass in total stator lamination mass
inp_fill = 0.55			# [-]   Fill Factor - relative amount of Copper mass in total coil mass

inp_fileName = 'New_Nastran_MAT9.bdf'	 
         
# ADDITIONAL PARAMETERS - optionally to be modified when using special alloys, epoxy mixtures, etc.

rho_Steel = 7.3		# [g/cm^3] Density of steel used in stator lamination
rho_Copper = 8.9	# [g/cm^3] Density of copper used in coils
rho_Epoxy01 = 1.2	# [g/cm^3] Density of epoxy in stator lamination
rho_Epoxy02 = 1.19	# [g/cm^3] Density of epoxy in coil airgaps

pr_Steel = 0.3		# [-] Poisson Ratio of steel used in stator lamination
pr_Copper = 0.17	# [-] Poisson Ratio of copper used in coils
pr_Epoxy01 = 0.4	# [-] Poisson Ratio of epoxy in stator lamination
pr_Epoxy02 = 0.37	# [-] Poisson Ratio of epoxy in coil airgaps

#
# ----------------------------------------- PROCESSING SECTION -------------------------------------------
#
# CALCULATION STEP 00

inp   [inp_ESteel, inp_ECopper, inp_EEpoxy01, inp_EEpoxy02, inp_stack, inp_fill]
pr = [pr_Steel, pr_Copper, pr_Epoxy01, pr_Epoxy02]
rho = [rho_Steel, rho_Copper, rho_Epoxy01, rho_Epoxy02]
Shear = []

for kk in range(0,4):
	Shear.append(inp[kk]/2/(1+pr[kk]))

#
# CALCULATION STEP 01 
# Index = 0 --> Lamination
# Index = 1 --> Coil
#
	
rho_C = []
pr_C  = []
pE_C  = []
pPR_C = []
pSH_C = []
aSH_C = []
E_C   = []


for xx in range(0,2):
	rho_C.append(round((rho[xx]*inp[xx+4])+(rho[xx+2]*(1-inp[xx+4])),2)*(10**(-9)))
	E_C.append((inp[xx]*inp[xx+4])+(inp[xx+2]*(1-inp[xx+4])))
	pr_C.append((pr[xx]*inp[xx+4])+(pr[xx+2]*(1-inp[xx+4])))
	pE_C.append(1/(inp[xx+4]/inp[xx] + (1-inp[xx+4])/inp[xx+2])*1.075)
	pPR_C.append(pr_C[xx]*pE_C[xx]/E_C[xx])
	pSH_C.append(E_C[xx]/2/(1+pr_C[xx]))
	aSH_C.append(Shear[xx]*Shear[xx+2]/(Shear[xx]*(1-inp[xx+4])+Shear[xx+2]*(inp[xx+4])))

#
# CALCULATION STEP 02
#

E1_Lam = E_C[0]
E2_Lam = E_C[0]
E3_Lam = pE_C[0]
GG12_Lam = pSH_C[0]
GG13_Lam = aSH_C[0]
GG23_Lam = aSH_C[0]
nu12_Lam = pr_C[0]
nu21_Lam = pr_C[0]
nu32_Lam = pPR_C[0]*pE_C[0]/E_C[0]
nu23_Lam = pPR_C[0]
nu31_Lam = pPR_C[0]*pE_C[0]/E_C[0]
nu13_Lam = pPR_C[0]


Q_Lam = (1-nu32_Lam*nu23_Lam-nu31_Lam*nu13_Lam-nu21_Lam*nu12_Lam-2*nu12_Lam*nu23_Lam*nu31_Lam)/(E1_Lam*E2_Lam*E3_Lam)

E1_Coil = pE_C[1]
E2_Coil = pE_C[1]
E3_Coil = E_C[1]
GG12_Coil = aSH_C[1]
GG13_Coil = pSH_C[1]
GG23_Coil = pSH_C[1]
nu12_Coil = pPR_C[1]
nu21_Coil = pPR_C[1]
nu32_Coil = pr_C[1]*pE_C[1]/E_C[1]
nu23_Coil = pr_C[1]
nu31_Coil = pr_C[1]*pE_C[1]/E_C[1]
nu13_Coil = pr_C[1]

Q_Coil = (1-nu32_Coil*nu23_Coil-nu31_Coil*nu13_Coil-nu21_Coil*nu12_Coil-2*nu12_Coil*nu23_Coil*nu31_Coil)/(E1_Coil*E2_Coil*E3_Coil)

#
# CALCULATION STEP 03
#

G11_Lam = int((1-nu32_Lam*nu23_Lam)/(E2_Lam*E3_Lam*Q_Lam))
G12_Lam = int((nu21_Lam+nu31_Lam*nu23_Lam)/(E2_Lam*E3_Lam*Q_Lam))
G13_Lam = int((nu31_Lam+nu21_Lam*nu32_Lam)/(E2_Lam*E3_Lam*Q_Lam))
G22_Lam = G11_Lam
G23_Lam = int((nu32_Lam+nu31_Lam*nu12_Lam)/(E3_Lam*E1_Lam*Q_Lam))
G33_Lam = int((1-nu12_Lam*nu21_Lam)/(E1_Lam*E2_Lam*Q_Lam))
G44_Lam = int(GG12_Lam)
G55_Lam = int(GG13_Lam)
G66_Lam = int(GG23_Lam)


Lam = [G11_Lam, G12_Lam, G13_Lam, G22_Lam, G23_Lam, G33_Lam, G44_Lam, G55_Lam, G66_Lam]
L_N = []

for ii in Lam:
	L_N.append('{:>8}'.format(str(ii)))

G11_Coil = int((1-nu32_Coil*nu23_Coil)/(E2_Coil*E3_Coil*Q_Coil))
G12_Coil = int((nu21_Coil+nu31_Coil*nu23_Coil)/(E2_Coil*E3_Coil*Q_Coil))
G13_Coil = int((nu31_Coil+nu21_Coil*nu32_Coil)/(E2_Coil*E3_Coil*Q_Coil))
G22_Coil = G11_Coil
G23_Coil = int((nu32_Coil+nu31_Coil*nu12_Coil)/(E3_Coil*E1_Coil*Q_Coil))
G33_Coil = int((1-nu12_Coil*nu21_Coil)/(E1_Coil*E2_Coil*Q_Coil))
G44_Coil = int(GG12_Coil)
G55_Coil = int(GG13_Coil)
G66_Coil = int(GG23_Coil)

Coil = [G11_Coil, G12_Coil, G13_Coil, G22_Coil, G23_Coil, G33_Coil, G44_Coil, G55_Coil, G66_Coil]
C_N = []


# ----------------------------------------- EXPORT SECTION -------------------------------------------


# Adapting calculated parameters to 8-digit Nastran format

for jj in Coil:
	C_N.append('{:>8}'.format(str(jj)))

# Create a placeholder variable for empty cells	

yy = '      0.'

# Writing out the MAT9 card. Stack factor and Fill factor are included in the material name for reasons of experience

f_out = open(inp_fileName, 'w+')
f_out.write('$ \n' + '$\n' + '$\n' + '$Lamination - Stack Factor ' + str(inp[4]) + '\n')	
f_out.write('MAT9    '+'       1' + L_N[0] + L_N[1] + L_N[2] + yy + yy + yy + L_N[3] + '+\n') 
f_out.write('        '+ L_N[4] + yy + yy + yy + L_N[5] + yy + yy + yy + '+\n')
f_out.write('        '+ L_N[6] + yy + yy + L_N[7] + yy + L_N[8] + '{:>8}'.format(str(rho_C[0])) + yy + '+\n')
f_out.write('        '+ yy + yy + yy + yy + yy + '     20.' + '    0.05\n')
f_out.write('$ \n' + '$\n' + '$\n' + '$Coil - Fill Factor ' + str(inp[5]) + '\n')
f_out.write('MAT9    '+'       2' + C_N[0] + C_N[1] + C_N[2] + yy + yy + yy + C_N[3] + '+\n') 
f_out.write('        '+ C_N[4] + yy + yy + yy + C_N[5] + yy + yy + yy + '+\n')
f_out.write('        '+ C_N[6] + yy + yy + C_N[7] + yy + C_N[8] + '{:>8}'.format(str(rho_C[1])) + yy + '+\n')
f_out.write('        '+ yy + yy + yy + yy + yy + '     20.' + '    0.09')
f_out.close()

#
#                       CREATION OF MAT9 CARD FINISHED - HAVE FUN WITH YOUR MODEL!
#

# ----------------------------------------- END OF PROGRAM -------------------------------------------





# SUGGESTED IMPROVEMENTS
#
# COMPLETED- Introduce clearly separated sections
# COMPLETED- Collect similar parameters in list variables
# COMPLETED- Introduce filename variable to avoid overwriting
# - Advanced Refactoring
# - Simplify list variables / create more functional variables
#

# NOTES
# 
# FIXED- Error in coil matrix calculation formula
#

