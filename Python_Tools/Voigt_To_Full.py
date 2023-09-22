from matscipy.elasticity import measure_triclinic_elastic_constants, fit_elastic_constants
from matscipy.elasticity import full_3x3x3x3_to_Voigt_6x6
from matscipy.elasticity import Voigt_6x6_to_full_3x3x3x3

import numpy as np


A = np.matrix([[11,12,13,0,0,0],[21,22,23,0,0,0],[31,32,33,0,0,0],[0,0,0,44,0,0],[0,0,0,0,55,0],[0,0,0,0,0,66]]);

print("In Voigt Notation")
print("")
print(A)


B = Voigt_6x6_to_full_3x3x3x3(A);
print("Full3x3x3x3")
print("")
#print(B)

print("0000", B[0,0,0,0])
print("0010",B[0,0,1,0])
print("0020",B[0,0,2,0])
print("1000", B[1,0,0,0])
print("1010",B[1,0,1,0])
print("1020",B[1,0,2,0])
print("2000", B[2,0,0,0])
print("2010",B[2,0,1,0])
print("2020",B[2,0,2,0])

print("0001", B[0,0,0,1])
print("0011",B[0,0,1,1])
print("0021",B[0,0,2,1])
print("1001", B[1,0,0,1])
print("1011",B[1,0,1,1])
print("1021",B[1,0,2,1])
print("2001", B[2,0,0,1])
print("2011",B[2,0,1,1])
print("2021",B[2,0,2,1])

print("0002", B[0,0,0,2])
print("0012",B[0,0,1,2])
print("0022",B[0,0,2,2])
print("1002", B[1,0,0,2])
print("1012",B[1,0,1,2])
print("1022",B[1,0,2,2])
print("2002", B[2,0,0,2])
print("2012",B[2,0,1,2])
print("2022",B[2,0,2,2])



print("0100", B[0,1,0,0])
print("0110",B[0,1,1,0])
print("0120",B[0,1,2,0])
print("1100", B[1,1,0,0])
print("1110",B[1,1,1,0])
print("1120",B[1,1,2,0])
print("2100", B[2,1,0,0])
print("2110",B[2,1,1,0])
print("2120",B[2,1,2,0])

print("0101", B[0,1,0,1])
print("0111",B[0,1,1,1])
print("0121",B[0,1,2,1])
print("1101", B[1,1,0,1])
print("1111",B[1,1,1,1])
print("1121",B[1,1,2,1])
print("2101", B[2,1,0,1])
print("2111",B[2,1,1,1])
print("2121",B[2,1,2,1])

print("0102", B[0,1,0,2])
print("0112",B[0,1,1,2])
print("0122",B[0,1,2,2])
print("1102", B[1,1,0,2])
print("1112",B[1,1,1,2])
print("1122",B[1,1,2,2])
print("2102", B[2,1,0,2])
print("2112",B[2,1,1,2])
print("2122",B[2,1,2,2])



print("0200", B[0,2,0,0])
print("0210",B[0,2,1,0])
print("0220",B[0,2,2,0])
print("1200", B[1,2,0,0])
print("1210",B[1,2,1,0])
print("1220",B[1,2,2,0])
print("2200", B[2,2,0,0])
print("2210",B[2,2,1,0])
print("2220",B[2,2,2,0])

print("0201", B[0,2,0,1])
print("0211",B[0,2,1,1])
print("0221",B[0,2,2,1])
print("1201", B[1,2,0,1])
print("1211",B[1,2,1,1])
print("1221",B[1,2,2,1])
print("2201", B[2,2,0,1])
print("2211",B[2,2,1,1])
print("2221",B[2,2,2,1])

print("0202", B[0,2,0,2])
print("0212",B[0,2,1,2])
print("0222",B[0,2,2,2])
print("1202", B[1,2,0,2])
print("1212",B[1,2,1,2])
print("1222",B[1,2,2,2])
print("2202", B[2,2,0,2])
print("2212",B[2,2,1,2])
print("2222",B[2,2,2,2])

#C = np.zeros((3,3,3,3));
#C[0,0,0,0]=1111;
#C[1,1,1,1]=2222;
#E = np.zeros((3,3,3,3));


#ct = 0
#for i in range(3):
#    for j in range(3):
#        for k in range(3):
#            for l in range(3):
#                if C[i,j,k,l] == 0:
#                    C[i,j,k,l] = ct
#                    C[j,i,k,l] = ct
#                    C[i,j,l,k] = ct
#                    ct += 1

#print(C)

#D=full_3x3x3x3_to_Voigt_6x6(C)

#print(D)
