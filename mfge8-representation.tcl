mol new mfge8_trial5bfearlyf199.pdb
# mol new mfge8_trial5bflatef199.pdb

# mol addfile step7.1_production.dcd step 10
# mol addfile step7.2_production.dcd step 10
# mol addfile step7.3_production.dcd step 10
# mol addfile step7.4_production.dcd step 10
# mol addfile step7.5_production.dcd step 10
# mol addfile step7.6_production.dcd step 10
# mol addfile step7.7_production.dcd step 10
# mol addfile step7.8_production.dcd step 10
# mol addfile step7.9_production.dcd step 10
# mol addfile step7.10_production.dcd step 10
# mol addfile step7.11_production.dcd step 10

# protein
mol color Name
mol representation NewCartoon 0.3 10 4.1 0
mol selection protein
mol material Opaque
mol modrep 0 top

# POPS
mol color colorID 2
mol representation Bonds 0.3 12
mol selection resname POPS
mol material Opaque
mol addrep top

# select the three residues that define the binding groove
mol color colorID 1
mol representation Bonds 0.3 12
mol selection resname LYS and resid 45
mol material Opaque
mol addrep top

mol color colorID 4
mol representation Bonds 0.3 12
mol selection resname ARG and resid 146
mol material Opaque
mol addrep top

mol color colorID 4
mol representation Bonds 0.3 12
mol selection resname ARG and resid 79
mol material Opaque
mol addrep top
