import numpy as np
import math
import os



r1_input = input("Enter values for r1 separated by commas: ")
r1 = np.array([float(val.strip()) for val in r1_input.split(",")])


r2_input = input("Enter values for r2 separated by commas: ")
r2 = np.array([float(val.strip()) for val in r2_input.split(",")])


r3_input = input("Enter values for r3 separated by commas: ")
r3 = np.array([float(val.strip()) for val in r3_input.split(",")])

r1Vector = np.array(r1)
r2Vector = np.array(r2)
r3Vector = np.array(r3)


def write_vectors_to_file(f, r1, r2, r3):
    f.write('\\\\')
    f.write("\\textbf{Vectors}:\n" + '\\vspace{1mm}')
    f.write("$\\vec{\\text{r}}_{1} = " + f"{r1[0]}i + {r1[1]}j + {r1[2]}k$ (\\text{{Km}}) \n" + '\n')
    f.write("$\\vec{\\text{r}}_{2} = " + f"{r2[0]}i + {r2[1]}j + {r2[2]}k$ (\\text{{Km}}) \n" + '\n')
    f.write("$\\vec{\\text{r}}_{3} = " + f"{r3[0]}i + {r3[1]}j + {r3[2]}k$ (\\text{{Km}}) \n " + '\n')


with open(f"equations/eq2.tex", "w") as f:
    write_vectors_to_file(f, r1, r2, r3)

# Call the function to write vectors to the file
r1Magnitude = np.linalg.norm(r1Vector)
r2Magnitude = np.linalg.norm(r2Vector)
r3Magnitude = np.linalg.norm(r3Vector)


def find_the_magnitude(f, r1, r2, r3, r1Magnitude, r2Magnitude, r3Magnitude):
    f.write('\\\\')
    f.write("\\textbf{Magnitude of the vector}:\n " +'\n \\vspace{1mm}')
    f.write("$\\text{r}_{1} = \\sqrt {" + f"{r1[0]}^2 +{r1[1]}^2 + {r1[2]}^2  $\n" + f" = {round(r1Magnitude, 3)}$ (\\text{{Km}}) \n" + '\n')
    f.write("$\\text{r}_{2} = \\sqrt {" + f"{r2[0]}^2 +{r2[1]}^2 + {r2[2]}^2  $\n" + f" = {round(r2Magnitude, 3)}$ (\\text{{Km}}) \n" + '\n')
    f.write("$\\text{r}_{3} = \\sqrt {" + f"{r3[0]}^2 +{r3[1]}^2 + {r3[2]}^2  $\n" + f" = {round(r3Magnitude, 3)}$ (\\text{{Km}}) \n " + '\n ')


with open(f"equations/magnitude.tex", "w") as f:
    find_the_magnitude(f, r1, r2, r3, r1Magnitude, r2Magnitude, r3Magnitude)


r12Coefficient = np.cross(r1Vector, r2Vector)
r23Coefficient = np.cross(r2Vector, r3Vector)
r31Coefficient = np.cross(r3Vector, r1Vector)


def find_the_coefficients(f, r1, r2, r3, r12Coefficient, r23Coefficient, r31Coefficient):
    f.write('\\\\')
    f.write("\\textbf{Step2: } Find \\hspace{1mm} \\vec{\\text{ C}}_{1}_{2}, \\vec{\\text{C}}_{2}_{3 } \\hspace{1mm} and \\hspace{1mm} \\vec{\\text{C}}_{3}_{1} \n"+'\n \\vspace{1mm}')
    f.write("$\\vec{\\text{C}}_{1}_{2} = $\n\n")
    f.write("\\begin{vmatrix}\n")
    f.write("\\begin{array}{|c @{\\hspace{0.5em}}c@{\\hspace{0.5em}}c|}\n")
    f.write("$i$ & $j$ & $k$ \\\\\n")
    f.write(f"${r1[0]}$ & ${r1[1]}$ & ${r1[2]}$ \\\\\n")
    f.write(f"${r2[0]}$ & ${r2[1]}$ & ${r2[2]}$ \\\\\n")
    f.write("\\end{array}\n")
    f.write("\\end{vmatrix}$\n")
    f.write(f"= {round(r12Coefficient[0], 2)}i + {round(r12Coefficient[1], 2)}j + {round(r12Coefficient[2], 2)}k$ (\\text{{Km}}^2) \n" + '\n')
    f.write("$\\vec{\\text{C}}_{2}_{3} = $\n\n")
    f.write("\\begin{vmatrix}\n")
    f.write("\\begin{array}{|c @{\\hspace{0.5em}}c@{\\hspace{0.5em}}c|}\n")
    f.write("$i$ & $j$ & $k$ \\\\\n")
    f.write(f"${r2[0]}$ & ${r2[1]}$ & ${r2[2]}$ \\\\\n")
    f.write(f"${r3[0]}$ & ${r3[1]}$ & ${r3[2]}$ \\\\\n")
    f.write("\\end{array}\n")
    f.write("\\end{vmatrix}$\n")
    f.write(f"= {round(r23Coefficient[0], 2)}i + {round(r23Coefficient[1], 2)}j + {round(r23Coefficient[2], 2)}k$ (\\text{{Km}}^2) \n" + '\n')
    f.write("$\\vec{\\text{C}}_{3}_{1} = $\n\n")

    f.write("\\begin{vmatrix}\n")
    f.write("\\begin{array}{|c @{\\hspace{0.5em}}c@{\\hspace{0.5em}}c|}\n")
    f.write("$i$ & $j$ & $k$ \\\\\n")
    f.write(f"${r3[0]}$ & ${r3[1]}$ & ${r3[2]}$ \\\\\n")
    f.write(f"${r3[0]}$ & ${r3[1]}$ & ${r3[2]}$ \\\\\n")
    f.write("\\end{array}\n")
    f.write("\\end{vmatrix}$\n")
    f.write(f"= {round(r31Coefficient[0], 2)}i + {round(r31Coefficient[1], 2)}j + {round(r31Coefficient[2], 2)}k$ (\\text{{Km}}^2) \n" + '\n \n')
    f.write("\n" + '\n')


with open(f"equations/coefficient.tex", "w") as f:
    find_the_coefficients(f, r1, r2, r3, r12Coefficient, r23Coefficient, r31Coefficient)

r23CoefficientMagitude = np.linalg.norm(r23Coefficient)
r1UnitVector = r1Vector / r1Magnitude
r23CoefficientUnitVector = r23Coefficient / r23CoefficientMagitude


def c23_magnitude(f, r23Coefficient):
    f.write('\\\\')
    f.write("\\textbf{Step3:} \\textbf{The magnitude of }  \\vec{\\text{C}}_{2}_{3 } :\n " +'\n \\vspace{3mm}')
    f.write(
        "$\\text{C}_{2}_{3} = \\sqrt {" + f"{r23Coefficient[0]}^2 +{r23Coefficient[1]}^2 + {r23Coefficient[2]}^2$\n")
    f.write(f" \\hspace{{9px}}= {round(r23CoefficientMagitude, 3)}$\n" +'\n')


with open(f"equations/c23CoefficientMagitude.tex", "w") as f:
    c23_magnitude(f, r23Coefficient)


def c23_unit_vector(f,r23Coefficient):
    r23CoefficientRound = [round(coord, 2) for coord in r23Coefficient]
    f.write('\\\\')
    f.write("\\textbf{The unit vector of } \\vec{\\text{C}}_{2}_{3 } :\n " +'\n \\vspace{3mm}')
    f.write(
        "$\\hat{\\text{C}}_{2}_{3}$ = $" + f"\\frac{{ {r23CoefficientRound[0]} i + {r23CoefficientRound[1]} j + {r23CoefficientRound[2]} k }}{{ {round(r23CoefficientMagitude, 3)} }} = {round(r23CoefficientUnitVector[0], 3)}i + {round(r23CoefficientUnitVector[1], 3)}j + {round(r23CoefficientUnitVector[2], 3)}k $\n")
    f.write("\\break")

with open(f"equations/c23UnitVector.tex", "w") as f:
    c23_unit_vector(f, r23Coefficient)

def r1_unit_vector(f,r1):
    f.write('\\\\')
    f.write("\\textbf{The unit vector of } \\vec{\\text{u}}_{r_{1} :\n " + '\n \\vspace{3mm}')
    f.write(
        "$\\hat{\\text{u}}_{r_{1}$ = $" + f"\\frac{{ {r1[0]} i + {r1[1]} j + {r1[2]} k }}{{ {round(r1Magnitude, 3)} }} = {round(r1UnitVector[0], 3)}i + {round(r1UnitVector[1], 3)}j + {round(r1UnitVector[2], 3)}k $\n"+'\n')


with open(f"equations/r1_unit_vector.tex", "w") as f:
    r1_unit_vector(f, r1)

def is_on_the_same_plane(r1UnitVector, r23CoefficientUnitVector):
    f.write('\\\\')
    f.write("\\textbf{Therefore,}:\n " + '\n \\vspace{3mm}')
    f.write("$\\hat{\\text{u}}_{r_{1}$ . \\hat{\\text{C}}_{2}_{3} = " + f"{(r1UnitVector.dot(r23CoefficientUnitVector))} \n"+'\n')
    if (np.fix(r1UnitVector.dot(r23CoefficientUnitVector)) == 0) :
        f.write("\\vspace{3mm}")
        f.write("\\begin{center}")
        f.write("This is close enough to or equal zero for our purposes. The three vectors \\text{r}_{1}, \\text{r}_{2}, and \\text{r}_{3}$ are coplanar.")
        f.write("\\end{center}")


with open(f"equations/is_on_the_same_plane.tex", "w") as f:
    is_on_the_same_plane(r1UnitVector, r23CoefficientUnitVector)
x = np.fix(r1UnitVector.dot(r23CoefficientUnitVector))

theN = r1Magnitude * r23Coefficient + r2Magnitude * r31Coefficient + r3Magnitude * r12Coefficient

theNMagnitude = np.linalg.norm(theN)

theD = r12Coefficient + r23Coefficient + r31Coefficient

theDMagnitude = np.linalg.norm(theD)

theS = r1Vector * (r2Magnitude - r3Magnitude) + r2Vector * (r3Magnitude - r1Magnitude) + r3Vector * (
        r1Magnitude - r2Magnitude)
def step_4():
    f.write('\\\\')
    f.write("\\textbf{Step4}:\n "+ '\n \\vspace{3mm}')
    f.write("\\vec{N} = $\\text{ r}_{1} \\vec{\\text{C}}_{2}_{3} + $\\text{r}_{2} \\vec{\\text{C}}_{3}_{1} + $\\text{r}_{3} \\vec{\\text{C}}_{1}_{2} \n" + '\n')
    f.write(f"\\hspace{{9px}} = {round(r1Magnitude,2)} ( {round(r23Coefficient[0], 2)}i + {round(r23Coefficient[1], 2)}j + {round(r23Coefficient[2], 2)}k) $ \n"+'\n')
    f.write(f"\\hspace{{9px}} + {round(r2Magnitude,2)} ( {round(r31Coefficient[0], 2)}i + {round(r31Coefficient[1], 2)}j + {round(r31Coefficient[2], 2)}k) $ \n"+'\n')
    f.write(f"\\hspace{{9px}} + {round(r3Magnitude,2)} ( {round(r12Coefficient[0], 2)}i + {round(r12Coefficient[1], 2)}j + {round(r12Coefficient[2], 2)}k) $ \n"+'\n')
    f.write("\\vspace{3mm}")
    f.write(f"\\vec{{N}} =( {round(theN[0], 2)}i + {round(theN[1], 2)}j + {round(theN[2], 2)}k) \n"+'\n' )
    f.write("\\vspace{3mm}")

    f.write("\\text{N} = \\sqrt {" + f"{round(theN[0],2)}^2 +{round(theN[1],2)}^2 + {round(theN[2],2)}^2$\n" +'\n')
    f.write("\\vspace{3mm}")

    f.write(f"\\hspace{{10px}}\\text = {round(theNMagnitude,3)} (\\text{{Km}}^2) \n" +'\n')
    f.write("\\vspace{5mm}")

    f.write(
        "\\vec{D} = $\\vec{\\text{ C}}_{1}_{2} + $ \\vec{\\text{C}}_{2}_{3} + $\\vec{\\text{C}}_{3}_{1} \n" + '\n')
    f.write(
        f"\\hspace{{9px}} =  ( {round(r12Coefficient[0], 2)}i + {round(r12Coefficient[1], 2)}j + {round(r12Coefficient[2], 2)}k) $ \n" + '\n')

    f.write(
        f"\\hspace{{9px}} +  ( {round(r31Coefficient[0], 2)}i + {round(r31Coefficient[1], 2)}j + {round(r31Coefficient[2], 2)}k) $ \n" + '\n')
    f.write(
        f"\\hspace{{9px}} +  ( {round(r23Coefficient[0], 2)}i + {round(r23Coefficient[1], 2)}j + {round(r23Coefficient[2], 2)}k) $ \n" + '\n')
    f.write("\\vspace{3mm}")
    f.write(f"\\vec{{D}} =( {round(theD[0], 2)}i + {round(theD[1], 2)}j + {round(theD[2], 2)}k) \n"+'\n' )
    f.write("\\vspace{3mm}")

    f.write("\\text{D} = \\sqrt {" + f"{round(theD[0], 2)}^2 +{round(theD[1], 2)}^2 + {round(theD[2], 2)}^2$\n" + '\n')
    f.write("\\vspace{3mm}")
    f.write(f"\\hspace{{10px}}\\text = {round(theDMagnitude,3)} (\\text{{Km}}^2) \n" +'\n')

    f.write("\\vspace{5mm}")

    f.write(
        "\\vec{S} = $(\\vec{\\text{ r}}_{1})(\\text{ r}}_{2} - \\text{ r}}_{3}) + $(\\vec{\\text{ r}}_{2})(\\text{ r}}_{3} - \\text{ r}}_{1}) + $(\\vec{\\text{ r}}_{3})(\\text{ r}}_{1} - \\text{ r}}_{2})  \n" + '\n')
    f.write(
    f"\\hspace{{9px}} =  ( {r1[0]}i + {r1[1]}j + {r1[2]}k$) ({round(r2Magnitude,3)} - {round(r3Magnitude,3)}) \n" + '\n')

    f.write(
        f"\\hspace{{9px}} =  ( {r2[0]}i + {r2[1]}j + {r2[2]}k$) ({round(r3Magnitude, 3)} - {round(r1Magnitude, 3)}) \n" + '\n')
    f.write(
        f"\\hspace{{9px}} =  ( {r3[0]}i + {r3[1]}j + {r3[2]}k$) ({round(r1Magnitude, 3)} - {round(r2Magnitude, 3)}) \n" + '\n')
    f.write("\\vspace{3mm}")

    f.write(f"\\vec{{S}} =( {round(theS[0], 2)}i + {round(theS[1], 2)}j + {round(theS[2], 2)}k) \n"+'\n' )
    f.write("\\vspace{3mm}")



with open(f"equations/step_4.tex", "w") as f:
    step_4()

v2Vector = ((math.sqrt(398600 / (theNMagnitude * theDMagnitude))) * ((np.cross(theD, r2Vector) / r2Magnitude) + theS))

def step_5():
    f.write('\\\\')
    f.write("\\textbf{Step5}:\n " + '\n \\vspace{3mm}')
    f.write("\\vec{\\text{v}}_{2} = \\sqrt{\\frac {\\mu}{\\text{N}{D}}} (\\frac{\\vec{D} \\times \\vec{\\text{r}}_{2}}  {\\text{r}}_{2}+ \\vec{S}) \n" +'\n')
    f.write(f"= \\sqrt{{\\frac{{398,600}}{{ ({round(theNMagnitude,2)} )({round(theDMagnitude,2)})}} }}\\times \n" +'\n \\vspace{5mm}')
    # f.write("\\frac{{b}{a}}")
    f.write(f"=\\left(\\frac{{\\begin{{vmatrix}}\n \\begin{{array}}{{|c @{{\\hspace{{0.5em}}}}c@{{\\hspace{{0.5em}}}}c|}}\n $i$ & $j$ & $k$ \\\\\n ${r1[0]}$ & ${r1[1]}$ & ${r1[2]}$ \\\\\n ${r2[0]}$ & ${r2[1]}$ & ${r2[2]}$ \\\\\n \\end{{array}}\n \\end{{vmatrix}}}}{{{r2Magnitude} }} +( {round(theS[0], 2)}i + {round(theS[1], 2)}j + {round(theS[2], 2)}k) \\right )\n\n")
    f.write("\\vspace{3mm}")
    f.write(f"=( {round(v2Vector[0], 2)}i + {round(v2Vector[1], 2)}j + {round(v2Vector[2], 2)}k) (\\text{{Km/s}}) ")

with open(f"equations/step_5.tex", "w") as f:
    step_5()
