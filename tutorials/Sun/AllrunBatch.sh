#!/bin/bash
#OAR -n granullarCollapseDense_Sun_perm_partial
#OAR -l /nodes=1/core=12,walltime=8:00:00
#OAR --stdout log.out
#OAR --stderr errors.err

source /etc/profile
module load openfoam/2312plus
export OMPI_MCA_plm_rsh_agent=oar-envsh

# Nombre de processus parallèles
NBCPUS=$(cat ${OAR_NODEFILE} | wc -l)

# create the mesh
#foamCleanPolyMesh
#blockMesh


# create the intial time folder
#cp -r 0_org 0
#module load openfoam/1906plus
#mapFields -sourceTime 202 ../phi580_refined64
# Decompose the case in order to run in parallel 
#funkySetFields -time 0

#decomposePar


# Lancement
mpirun -np ${NBCPUS} -machinefile  ${OAR_NODEFILE} sedFoam_rbgh -parallel 

