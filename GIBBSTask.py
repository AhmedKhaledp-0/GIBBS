import numpy as np
import math
r1 = [-294.32, 4265.1, 5986.7 ]
r2 = [-1365.5,3637.6, 6346.8 ]
r3 = [-2940.3 , 2473.7,6555.8]
r1Vector = np.array(r1)
r2Vector = np.array(r2)
r3Vector = np.array(r3)

r1Magnitude = np.linalg.norm(r1Vector)
r2Magnitude = np.linalg.norm(r2Vector)
r3Magnitude = np.linalg.norm(r3Vector)

print (r1Magnitude,r2Magnitude,r3Magnitude)

r12Coefficient = np.cross(r1Vector, r2Vector)
r23Coefficient = np.cross(r2Vector, r3Vector)
r31Coefficient = np.cross(r3Vector, r1Vector)

r23CoefficientMagitude = np.linalg.norm(r23Coefficient)

r1UnitVector = r1Vector / r1Magnitude
r23CoefficientUnitVector = r23Coefficient/r23CoefficientMagitude
print (r23CoefficientUnitVector)
print (r1UnitVector)

x =np.fix( r1UnitVector.dot(r23CoefficientUnitVector))
print(x)
if (x == 0):
    print(0)
else:
    print(1)

theN = r1Magnitude * r23Coefficient + r2Magnitude * r31Coefficient + r3Magnitude * r12Coefficient
print (theN)

theNMagnitude = np.linalg.norm(theN)
print (theNMagnitude)

theD = r12Coefficient + r23Coefficient + r31Coefficient

theDMagnitude = np.linalg.norm(theD)
print (theDMagnitude)

theS = r1Vector * (r2Magnitude - r3Magnitude) + r2Vector * (r3Magnitude - r1Magnitude) + r3Vector * (r1Magnitude -r2Magnitude)
print (theS)

v2Vector = ((math.sqrt(398600/(theNMagnitude * theDMagnitude)) ) * ( (np.cross(theD, r2Vector)/r2Magnitude)+theS ))
print(v2Vector)
