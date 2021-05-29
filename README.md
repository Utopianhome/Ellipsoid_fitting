# Ellipsoid_fitting
### API:

__int ellipsoid_fitting(struct Matrix data, struct Matrix D, struct ellipsoid character);__

__int ellipsoid_fitting_lasting(struct Matrix data, struct Matrix D, struct ellipsoid character, double appenddata[3]);__

input:  [x1, y1, z1, x2, y2, z2, .......] 

output: offset[3] and gain[3]
