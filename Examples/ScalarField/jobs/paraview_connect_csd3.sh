module load gcc/11.3.0
module use /rds/project/rds-YVo7YUJF2mk/shared/modules
module load paraview/5.11.0-osmesa
pvserver --no-mpi
# Then ssh -NL 11111:localhost:11111 fawcett.maths.cam.ac.uk locally