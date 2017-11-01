#PBS -l nodes=10:ppn=10
#PBS -m abe
#PBS -M 21516774@student.uwa.edu.au
source /etc/bash.bashrc
mpirun lattice -size 1000 -p 0.5 -s