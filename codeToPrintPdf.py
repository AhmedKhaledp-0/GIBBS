import numpy as np
import math

# Input the 3 spacecraft position vectors
# r1 = [-294.32, 4265.1, 5986.7 ]
# r2 = [-1365.5,3637.6, 6346.8 ]
# r3 = [-2940.3 , 2473.7,6555.8] 
r1_input = input("Enter values for r1 separated by commas: ")
r1 = np.array([float(val.strip()) for val in r1_input.split(",")])
r2_input = input("Enter values for r2 separated by commas: ")
r2 = np.array([float(val.strip()) for val in r2_input.split(",")])
r3_input = input("Enter values for r3 separated by commas: ")
r3 = np.array([float(val.strip()) for val in r3_input.split(",")])

def c_product(first, second):
    return(np.cross(first, second))

# Convert the vectors to arrays
r1_vector = np.array(r1)
r2_vector = np.array(r2)
r3_vector = np.array(r3)

# Get the 3 vectors magnitude
r1_magnitude = np.linalg.norm(r1_vector)
r2_magnitude = np.linalg.norm(r2_vector)
r3_magnitude = np.linalg.norm(r3_vector)

# Find the 3 coefficient
r12_coefficient = np.cross(r1_vector, r2_vector)
r23_coefficient = np.cross(r2_vector, r3_vector)
r31_coefficient = np.cross(r3_vector, r1_vector)

# Find the nuit vector for the C23 & r1
r23_coefficient_magitude = np.linalg.norm(r23_coefficient)
r1_unit_vector = r1_vector / r1_magnitude
r23_coefficientUnitVector = r23_coefficient / r23_coefficient_magitude

# Check if the 3 vectors lay on the same plane by dot product C23 & r1 uint vectors
x = np.fix(r1_unit_vector.dot(r23_coefficientUnitVector))

# Get the N, D, and S 
n_vector = r1_magnitude * r23_coefficient + r2_magnitude * r31_coefficient + r3_magnitude * r12_coefficient
n_magnitude = np.linalg.norm(n_vector)
d_vector = r12_coefficient + r23_coefficient + r31_coefficient
d_magnitude = np.linalg.norm(d_vector)
s_vector = r1_vector * (r2_magnitude - r3_magnitude) + r2_vector * (r3_magnitude - r1_magnitude) + r3_vector * (r1_magnitude - r2_magnitude)

# Find the V2 vector 
v2_vector = (math.sqrt(398600 / (n_magnitude * d_magnitude))) * ((c_product(d_vector,r2_vector) / r2_magnitude) + s_vector)

# Assign the r2 position vector & v2 to get the orbital 6 elements
r = r2_vector
v = v2_vector
r_magnitude = np.linalg.norm(r)
v_magnitude = np.linalg.norm(v)

# Get the rabidal velocity 
radial_velocity = (r.dot(v)) / r_magnitude

# Get the specific angular momentum vector & magnitude
specific_angular_momentum = np.cross(r, v)
specific_angular_momentum_magnitude = np.linalg.norm(specific_angular_momentum)

# Calculate the inclination 
i = math.acos(specific_angular_momentum[2]/specific_angular_momentum_magnitude) * (180.0 / math.pi)

# Calculate the vector that defines the node line & it's magnitude
node_line_vector = np.cross((0,0,1),specific_angular_momentum)
node_line_magnitude = np.linalg.norm(node_line_vector)

# Calculate the right ascension of the ascending node
if node_line_vector[1] >= 0:
    right_ascension_of_the_ascending_node = math.acos(node_line_vector[0]/node_line_magnitude) * (180.0 / math.pi)
else:
    right_ascension_of_the_ascending_node =360 - (math.acos(node_line_vector[0]/node_line_magnitude)) * (180.0 / math.pi)

# Calculate the eccentricity vector & magnitude
eccentricity_vector = (1/398600)*((v_magnitude**2 - 398600/r_magnitude)*r - (r_magnitude*radial_velocity)*v)
eccentricity_magnitude = np.linalg.norm(eccentricity_vector)

# Calculate the argument of perigee
if eccentricity_vector[2] >= 0:
     argument_of_perigee = math.acos(node_line_vector.dot(eccentricity_vector)/(node_line_magnitude*eccentricity_magnitude)) * (180.0 / math.pi)
else:
     argument_of_perigee =360 - math.acos(node_line_vector.dot(eccentricity_vector)/(node_line_magnitude*eccentricity_magnitude)) * (180.0 / math.pi)

# Calculate the true anomaly
if  radial_velocity >= 0:
    true_anomaly = math.acos ((eccentricity_vector.dot(r))/(eccentricity_magnitude * r_magnitude)) * (180 / math.pi)
else:
    true_anomaly = 360 - (math.acos ((eccentricity_vector.dot(r))/(eccentricity_magnitude * r_magnitude)) * (180 / math.pi))

# Calculate the perigee radius
perigee_radii = (specific_angular_momentum_magnitude**2/398600)*(1/(1+eccentricity_magnitude * math.cos(0*math.pi /180)))

# Calculate the apogee radius
apogee_radii = (specific_angular_momentum_magnitude**2/398600)*(1/(1+eccentricity_magnitude * math.cos(180*math.pi /180)))

# Calculate the semi major axis
semimajor_axis = (1/2) * (perigee_radii+apogee_radii)

# Print the algorithm by latex to PDF
def print_the_algorithm_by_latex():
    with open(f"the_output/main.tex", "w") as f:

        # Define the document
        f.write("\\documentclass{article} "+'\n')
        f.write("\\usepackage{amsmath}"+'\n')
        f.write("\\usepackage {geometry}"+'\n')        
        f.write("\\newcommand {\\cosinv}{\\cos^{-1}}"+'\n')
        f.write("\\geometry{a4paper}"+'\n')

        f.write("\\begin{document}"+'\n')
        f.write("\\section{GIBBS}"+'\n')
        f.write("\\begin{math}"+'\n')

        # Print the 3 posithion vectors 
        f.write("\\textbf{Vectors}:" + '$')
        f.write("$$\\vec{r}_1 = " + f"{r1[0]} \\hat{{i}} + {r1[1]} \\hat{{j}} + {r1[2]} \\hat{{k}} \\texttt{{ Km}}$$"+'\n')
        f.write("$$\\vec{r}_2 = " + f"{r2[0]} \\hat{{i}} + {r2[1]} \\hat{{j}} + {r2[2]} \\hat{{k}} \\texttt{{ Km}}$$"+'\n')
        f.write("$$\\vec{r}_3 = " + f"{r3[0]} \\hat{{i}} + {r3[1]} \\hat{{j}} + {r3[2]} \\hat{{k}} \\texttt{{ Km}}$$"+'\n')

        # Print the 3 posithion vectors magnitude 
        f.write("\\textbf{Magnitude of the vector}:\n"+'\n')
        f.write("$${r}_1 = \\sqrt {" + f"{r1[0]}^2 +{r1[1]}^2 + {r1[2]}^2 }} " + f" = {round(r1_magnitude, 3)} \\texttt{{ Km}}$$" +'\n')
        f.write("$${r}_2 = \\sqrt {" + f"{r2[0]}^2 +{r2[1]}^2 + {r2[2]}^2 }} " + f" = {round(r2_magnitude, 3)} \\texttt{{ Km}}$$" +'\n')
        f.write("$${r}_3 = \\sqrt {" + f"{r3[0]}^2 +{r3[1]}^2 + {r3[2]}^2 }} " + f" = {round(r3_magnitude, 3)} \\texttt{{ Km}}$$" +'\n')

        ## Find the 3 coefficient
        f.write("$\\textbf{Step2:} Find { } \\vec{C}_{12}, \\vec{C}_{23}  \\text{ and } \\vec{C}_{31} $\n"+'\n')
        f.write(f"$$\\vec{{C}}_{{12}}=  \\begin{{vmatrix}} \\hat{{i}} & \\hat{{j}} & \\hat{{k}} \\\\ {r1[0]} & {r1[1]} & {r1[2]} \\\\ {r2[0]} & {r2[1]} & {r2[2]} \\end{{vmatrix}}= {round(r12_coefficient[0], 2)}i + {round(r12_coefficient[1], 2)}j + {round(r12_coefficient[2], 2)}k (\\texttt{{Km}}^2)$$ \n" + '\n')
        f.write(f"$$\\vec{{C}}_{{23}}=  \\begin{{vmatrix}} \\hat{{i}} & \\hat{{j}} & \\hat{{k}} \\\\ {r2[0]} & {r2[1]} & {r2[2]} \\\\ {r3[0]} & {r3[1]} & {r3[2]} \\end{{vmatrix}}= {round(r23_coefficient[0], 2)}i + {round(r23_coefficient[1], 2)}j + {round(r23_coefficient[2], 2)}k (\\texttt{{Km}}^2)$$ \n" + '\n')
        f.write(f"$$\\vec{{C}}_{{31}}=  \\begin{{vmatrix}} \\hat{{i}} & \\hat{{j}} & \\hat{{k}} \\\\ {r3[0]} & {r3[1]} & {r3[2]} \\\\ {r1[0]} & {r1[1]} & {r1[2]} \\end{{vmatrix}}= {round(r31_coefficient[0], 2)}i + {round(r31_coefficient[1], 2)}j + {round(r31_coefficient[2], 2)}k (\\texttt{{Km}}^2)$$ \n")
        
        # Find C23 magnitude
        f.write("$\\textbf{Step3:The magnitude of } \\vec{C}_{23}:$\n" +'\n')
        f.write("$${C}_{23} = \\sqrt {" + f"{round(r23_coefficient[0],3)}^2 +{round(r23_coefficient[1],2)}^2 + {round(r23_coefficient[2],2)}^2 }} = {round(r23_coefficient_magitude, 3) }(\\texttt{{Km}}^2) $$\n")
        r23_coefficientRound = [round(coord, 2) for coord in r23_coefficient]

        # Find the nuit vector for the C23 & r1
        f.write("$\\textbf{The unit vector of } \\vec{C}_{23} : $\n " +'\n')
        f.write("$$\\hat{{C}_{23}} =" + f"\\left(\\frac{{ {r23_coefficientRound[0]} i + {r23_coefficientRound[1]} j + {r23_coefficientRound[2]} k }}{{ {round(r23_coefficient_magitude, 3)} }}\\right) = {round(r23_coefficientUnitVector[0], 3)}i + {round(r23_coefficientUnitVector[1], 3)}j + {round(r23_coefficientUnitVector[2], 3)}k $$\n")
        f.write("$\\textbf{The unit vector of } \\vec{r}_1 : $ \n" + '\n' +'\n')
        f.write("$$\\hat{{u}_{r_1}}= "+ f"\\left(\\frac{{ {r1[0]} i + {r1[1]} j + {r1[2]} k }}{{ {round(r1_magnitude, 3)} }}\\right) = {round(r1_unit_vector[0], 3)}i + {round(r1_unit_vector[1], 3)}j + {round(r1_unit_vector[2], 3)}k $$\n")

        # Check if the 3 vectors lay on the same plane by dot product C23 & r1 uint vectors
        f.write("$\\textbf{Therefore}: $\n " + '\n ')
        f.write("$$\\hat{{u}_{r_1}} . \\hat{{C}_{23}} = " + f"{(r1_unit_vector.dot(r23_coefficientUnitVector))} $$\n")
        if (np.fix(r1_unit_vector.dot(r23_coefficientUnitVector)) == 0) :
            f.write("\\begin{center}")
            f.write("\\text{This is close enough to or equal zero for our purposes. }\n"+'\n')
            f.write("\\text{The three vectors ${r}_{1}$, ${r}_{2}$, and ${r}_{3}$ are coplanar. }\n"+'\n')
            f.write("\\end{center}")

        # Get the N vector & magnitude 
        f.write("$\\textbf{Step4}:\n $"+'\n')
        f.write("\\begin{align*}\\vec{N} &= {r}_1 \\vec{C}_{23} + {r}_2 \\vec{C}_{31} + {r}_3 \\vec{C}_{12}  \\\\" + '\n')
        f.write(f"&= {round(r1_magnitude,2)} ( {round(r23_coefficient[0], 2)}i + {round(r23_coefficient[1], 2)}j + {round(r23_coefficient[2], 2)}k) \\\\"+'\n')
        f.write(f"&+ {round(r2_magnitude,2)} ( {round(r31_coefficient[0], 2)}i + {round(r31_coefficient[1], 2)}j + {round(r31_coefficient[2], 2)}k) \\\\"+'\n')
        f.write(f"&+ {round(r3_magnitude,2)} ( {round(r12_coefficient[0], 2)}i + {round(r12_coefficient[1], 2)}j + {round(r12_coefficient[2], 2)}k)  \\\\"+'\n')
        f.write(f"\\vec{{N}} &=({round(n_vector[0], 2)}i + {round(n_vector[1], 2)}j + {round(n_vector[2], 2)}k)  (\\texttt{{Km}}^3) \\\\"+'\n' )
        f.write("\\text{N} &= \\sqrt {" + f"{round(n_vector[0],2)}^2 +{round(n_vector[1],2)}^2 + {round(n_vector[2],2)}^2 }} \\\\ " +'\n')
        f.write(f"&= {round(n_magnitude,3)} (\\texttt{{Km}}^3) \\end{{align*}}\n" +'\n')

        # Get the D vector & magnitude 
        f.write("\\begin{align*}\\vec{D} &= \\vec{C}_{12} + \\vec{C}_{23} + \\vec{C}_{31} \\\\" + '\n')
        f.write(f"&=  ( {round(r12_coefficient[0], 2)}i + {round(r12_coefficient[1], 2)}j + {round(r12_coefficient[2], 2)}k)  \\\\" + '\n')
        f.write(f"&+  ( {round(r31_coefficient[0], 2)}i + {round(r31_coefficient[1], 2)}j + {round(r31_coefficient[2], 2)}k)  \\\\" + '\n')
        f.write(f"&+  ( {round(r23_coefficient[0], 2)}i + {round(r23_coefficient[1], 2)}j + {round(r23_coefficient[2], 2)}k)  \\\\" + '\n')
        f.write(f"\\vec{{D}} &=( {round(d_vector[0], 2)}i + {round(d_vector[1], 2)}j + {round(d_vector[2], 2)}k)  (\\texttt{{Km}}^2)\\\\"+'\n' )
        f.write("\\text{D} &= \\sqrt {" + f"{round(d_vector[0], 2)}^2 +{round(d_vector[1], 2)}^2 + {round(d_vector[2], 2)}^2 }} \\\\" + '\n')
        f.write(f"&= {round(d_magnitude,3)}  (\\texttt{{Km}}^2) \\end{{align*}}\n" +'\n')

        # Get the S vector 
        f.write("\\begin{align*}\\vec{S} &= (\\vec{r}_1)({r}_2 - {r}_3) + (\\vec{r}_2)({r}_3 - {r}_1) + (\\vec{r}_3)({r}_1 - {r}_2) \\\\ " + '\n')
        f.write(f"&=  ( {r1[0]}i + {r1[1]}j + {r1[2]}k) ({round(r2_magnitude,3)} - {round(r3_magnitude,3)}) \\\\" + '\n')
        f.write(f"&=  ( {r2[0]}i + {r2[1]}j + {r2[2]}k) ({round(r3_magnitude, 3)} - {round(r1_magnitude, 3)}) \\\\" + '\n')
        f.write(f"&=  ( {r3[0]}i + {r3[1]}j + {r3[2]}k) ({round(r1_magnitude, 3)} - {round(r2_magnitude, 3)}) \\\\" + '\n')
        f.write(f"\\vec{{S}} &=( {round(s_vector[0], 2)}i + {round(s_vector[1], 2)}j + {round(s_vector[2], 2)}k)  (\\texttt{{Km}}^2) \\end{{align*}} "+'\n' )
        
        # Find the V2 vector 
        f.write("$\\textbf{Step5}: $\n " + '\n')
        f.write("$$\\vec{v}_2 = \\sqrt{\\left(\\frac {\\mu}{{N}{D}}\\right) } \\left(\\frac{\\vec{D} \\times \\vec{\\text{r}}_{2}}  {{\\text{r}}_{2}} + \\vec{S}\\right)$$ \n" +'\n')
        f.write(f"\\begin{{align*}} \\vec{{v}}_2 &= \\sqrt{{\\frac{{398,600}}{{ ({round(n_magnitude,2)} )({round(d_magnitude,2)})}} }} \\\\" +'\n ')
        f.write(f" &\\times \\left[\\frac{{\\begin{{vmatrix}}\n $i$ & $j$ & $k$ \\\\\n ${r1[0]}$ & ${r1[1]}$ & ${r1[2]}$ \\\\\n ${r2[0]}$ & ${r2[1]}$ & ${r2[2]}$ \\\\\n  \\end{{vmatrix}}}}{{{round(r2_magnitude,2)} }} +( {round(s_vector[0], 2)}i + {round(s_vector[1], 2)}j + {round(s_vector[2], 2)}k) \\right ]\\\\")
        f.write(f"&=( {round(v2_vector[0], 2)}i + {round(v2_vector[1], 2)}j + {round(v2_vector[2], 2)}k) \\texttt{{ (Km/s)}} \\end{{align*}} \n" +'\n' )
        f.write("\\break \n"+'\n')

        # Assign the r2 position vector & v2 with them magnitude to get the orbital 6 elements
        f.write("$\\textbf{Calculate the orbital elements}$\n"+ '\n')
        f.write("$\\text{we have }\\vec{r} \\text{ and } \\vec{v}$ \n"+'\n')
        f.write(f"$$\\vec{{r}} = ({round(r[0],2)}) \\hat{{i}} + ({round(r[1],2)}) \\hat{{j}} + ({round(r[2],2)}) \\hat{{k}} \\texttt{{ (Km)}}$$ \n"+'\n')
        f.write(f"$$\\vec{{v}} = ({round(v[0],2)}) \\hat{{i}} + ({round(v[1],2)}) \\hat{{j}} + ({round(v[2],2)}) \\hat{{k}} \\texttt{{ (Km/s)}}$$ \n"+'\n')
        f.write("\\textbf{Step1} \n"+'\n')
        f.write(f"$${{r}} = \\sqrt{{\\vec{{r}} \\text{{ . }} \\vec{{r}}}} = \\sqrt{{(({round(r[0],2)}))^2 + (({round(r[1],2)}))^2 + (({round(r[2],2)}))^2}} = {round(r_magnitude,2)} \\texttt{{ (Km)}} $$ \n"+'\n')
        f.write("\\textbf{Step2} \n"+'\n')
        f.write(f"$${{v}} = \\sqrt{{\\vec{{v}} \\text{{ . }} \\vec{{v}}}} = \\sqrt{{(({round(v[0],2)}))^2 + (({round(v[1],2)}))^2 + (({round(v[2],2)}))^2}} = {round(v_magnitude,2)} \\texttt{{ (Km/s)}} $$ \n"+'\n')

        # Get the rabidal velocity & get information about spacecraft position 
        f.write("\\textbf{Step3} \n"+'\n')
        f.write(f"$${{v}}_r = \\frac{{\\vec{{v}}.\\vec{{r}}}}{{r}} = \\frac{{(({round(v[0],2)})).(({round(r[0],2)})) + (({round(v[1],2)})).(({round(r[1],2)})) + (({round(v[2],2)})).(({round(r[2],2)})) }}{{{round(r_magnitude,2)}}} = {round(radial_velocity,2)} \\texttt{{ (Km/s)}}$$ \n"+'\n')
        if radial_velocity > 0:
            f.write("$ \\text{Since }{v}_r >\\text{ 0 , the satellite is flying away from perigee.} $ \n"+'\n')
        else:
            f.write("$ \\text{Since }{v}_r <\\text{ 0 , the satellite is is flying toward perigee.} $ \n"+'\n')

        # Get the specific angular momentum vector & magnitude
        f.write("\\textbf{Step4} \n"+'\n')
        f.write(f"$$ \\vec{{h}} = \\vec{{r}} \\times \\vec{{v}} = \\begin{{vmatrix}} \\hat{{i}} & \\hat{{j}} & \\hat{{k}} \\\\ ({round(r[0],2)}) & ({round(r[1],2)}) & ({round(r[2],2)}) \\\\ ({round(v[0],2)}) & ({round(v[1],2)}) & ({round(v[2],2)}) \\end{{vmatrix}} = {round(specific_angular_momentum[0],2)} \\hat{{i}} + {round(specific_angular_momentum[1],2)} \\hat{{j}} + {round(specific_angular_momentum[2],2)} \\hat{{k}} \\texttt{{ (Km$^{2}$/s)}} $$\n"+'\n')
        f.write("\\textbf{Step5} \n"+'\n')
        f.write(f"$${{h}} = \\sqrt{{\\vec{{h}} \\text{{ . }} \\vec{{h}}}} = \\sqrt{{({round(specific_angular_momentum[0],2)})^2 + ({round(specific_angular_momentum[1],2)})^2 + ({round(specific_angular_momentum[2],2)})^2}} = {round(specific_angular_momentum_magnitude,2)} \\texttt{{ (Km$^{2}$/s)}} $$ \n"+'\n')

        # Calculate the inclination 
        f.write("\\textbf{Step6} \n"+'\n')
        f.write(f"$$ {{i}} = \\cosinv{{\\frac{{{{h}}_z}}{{h}}}} = \\cosinv{{\\left(\\frac{{{round(specific_angular_momentum[2],2)}}}{{{round(specific_angular_momentum_magnitude,2)}}}\\right)}} \\Rightarrow  {round(i,2)}^ \\circ $$ \n" +'\n')
        if i == 0:
            f.write("$ \\text{Since }{i} =\\text{ 0 , this is equatorial orbit.} $ \n"+'\n')
        elif i == 90:
            f.write("$ \\text{Since }{i} =\\text{ 90 , this is polar orbit.} $ \n"+'\n')
        elif i == 180:
            f.write("$ \\text{Since }{i} =\\text{ 180 , this is equatorial orbit.} $ \n"+'\n')
        elif i > 0 and i < 90:
            f.write("$ \\text{Since 0} > {i} < \\text{ 90 , this is prograde orbit.} $ \n"+'\n')
        else:
            f.write("$ \\text{Since 90} > {i} < \\text{ 180 , this is retrograde orbit.} $ \n"+'\n')

        # Calculate the vector that defines the node line 
        f.write("\\textbf{Step7} \n"+'\n')
        f.write(f"$$ \\vec{{N}} = \\vec{{\\hat{{K}}}} \\times \\vec{{h}} =  \\begin{{vmatrix}} \\hat{{i}} & \\hat{{j}} & \\hat{{k}} \\\\ {0} & {0} & {1} \\\\ {round(specific_angular_momentum[0],2)} & {round(specific_angular_momentum[1],2)} & {round(specific_angular_momentum[2],2)} \\end{{vmatrix}} = {round(node_line_vector[0],2)} \\hat{{i}} + {round(node_line_vector[1],2)} \\hat{{j}} + {round(node_line_vector[2],2)} \\hat{{k}} \\texttt{{ (Km$^{2}$/s)}} $$\n"+'\n')

        # Calculate the the node line magnitude
        f.write("\\textbf{Step8} \n"+'\n')
        f.write(f"$${{N}} = \\sqrt{{\\vec{{N}} \\text{{ . }} \\vec{{N}}}} = \\sqrt{{({round(node_line_vector[0],2)})^2 + ({round(node_line_vector[1],2)})^2 + ({round(node_line_vector[2],2)})^2}} = {round(node_line_magnitude,2)} \\texttt{{ (Km$^{2}$/s)}} $$ \n"+'\n')
        f.write("\\break")

        # Calculate the right ascension of the ascending node
        f.write("\\textbf{Step9} \n"+'\n')
        if node_line_vector[1] < 0:
            f.write("$$\\text{we know that }{N}_Y < \\text{ 0 ; therefore } \\Omega \\text{ must lie in the third or foruth quadrant}$$")
        else:
            f.write("$$\\text{we know that }{N}_Y \\geq  \\text{ 0 ; therefore } \\Omega \\text{ must lie in the first or second quadrant}$$")

        f.write(f"$$ {{\\Omega}} = \\cosinv{{\\frac{{{{N}}_x}}{{N}}}} = \\cosinv{{\\left(\\frac{{{round(node_line_vector[1],2)}}}{{{round(node_line_magnitude,2)}}}\\right)}} \\Rightarrow  {round(right_ascension_of_the_ascending_node,2)}^ \\circ $$ \n" +'\n')

        # Calculate the eccentricity vector & magnitude
        f.write("\\textbf{Step10} \n"+'\n')
        f.write(f"$$\\vec{{e}} = \\frac{{1}}{{\\mu}} \\left[ \\left( {{v}}^2 - \\frac{{\\mu}}{{r}} \\right) \\vec{{r}} - {{r}} {{v}}_r \\vec{{v}} \\right] $$\n"+'\n')       
        f.write(f"\\begin{{align*}}\\vec{{e}} &= \\frac{{1}}{{398600}} \\left[ \\left( {{{round(v_magnitude,2)}}}^2 - \\frac{{398600}}{{{round(r_magnitude,2)}}} \\right)\\right]  (({round(r[0],2)}) \\hat{{i}} + ({round(r[1],2)}) \\hat{{j}} + ({round(r[2],2)}) \\hat{{k}}) \\\\"+'\n')
        f.write(f"&- {round(r_magnitude,2)} * {round(radial_velocity,2)} \\text{{  }} (({round(v[0],2)}) \\hat{{i}} + ({round(v[1],2)}) \\hat{{j}} + ({round(v[2],2)}) \\hat{{k}}) \\\\"+'\n')
        f.write(f"\\vec{{e}} &= ({round(eccentricity_vector[0],2)} \\hat{{i}} + {round(eccentricity_vector[1],2)} \\hat{{j}} + {round(eccentricity_vector[2],2)} \\hat{{k}}) \\end{{align*}} \n"+'\n')

        # Calculate the eccentricity magnitude & check the orbit shape
        f.write("\\textbf{Step11} \n"+'\n')
        f.write(f"$${{e}} = \\sqrt{{\\vec{{e}} \\text{{ . }} \\vec{{e}}}} = \\sqrt{{({round(eccentricity_vector[0],2)})^2 + ({round(eccentricity_vector[1],2)})^2 + ({round(eccentricity_vector[2],2)})^2}} = {round(eccentricity_magnitude,2)} $$ \n"+'\n')
        if eccentricity_magnitude == 0:
            f.write("$\\text{Clearly, the orbit is a circle.} $")
        elif eccentricity_magnitude > 0 and eccentricity_magnitude <1 :
            f.write("$\\text{Clearly, the orbit is an ellipse.} $ \n"+'\n')
        elif eccentricity_magnitude == 1:
            f.write("$\\text{Clearly, the orbit is an parabola.} $ \n"+'\n')
        else:
            f.write("$\\text{Clearly, the orbit is an Hyperabola.} $ \n"+'\n')

        # Calculate the argument of perigee
        f.write("\\textbf{Step12} \n"+'\n')
        if eccentricity_vector[2] >= 0:
            f.write("$$\\text{We know that }{e}_z \\geq  \\text{ 0 ; therefore } \\omega \\text{ must lie in the first or second quadrant}$$")
        else:
            f.write("$$\\text{We know that }{e}_z <  \\text{ 0 ; therefore } \\omega \\text{ must lie in the third or fourth quadrant}$$")
        
        f.write(f"$$ {{\\omega}} = \\cosinv{{\\frac{{\\vec{{N}}.\\vec{{e}}}}{{N e}}}} = \\cosinv{{\\left[\\frac{{ ({round(node_line_vector[0],2)}) ({round(eccentricity_vector[0],2)}) + ({round(node_line_vector[1],2)}) ({round(eccentricity_vector[1],2)}) + ({round(node_line_vector[2],2)}) ({round(eccentricity_vector[2],2)}) }}{{({round(node_line_magnitude,2)})({round(eccentricity_magnitude,2)})}}\\right]}} \\Rightarrow  {round(argument_of_perigee,2)}^ \\circ $$ \n" +'\n')

        # Calculate the true anomaly
        f.write("\\textbf{Step13} \n"+'\n')
        if radial_velocity >= 0:
            f.write("$$\\text{We know that }{v}_r \\geq  \\text{ 0 , which means 0} ^ \\circ \\geq \\theta < \\text{ 180} ^ \\circ {. Therefore,}$$")
        else:
            f.write("$$\\text{We know that }{v}_r <  \\text{ 0 , which means 180} ^ \\circ \\geq \\theta < \\text{ 360} ^ \\circ {. Therefore,}$$")
        
        f.write(f"$$ {{\\theta}} = \\cosinv{{\\frac{{\\vec{{e}}.\\vec{{r}}}}{{e r}}}} = \\cosinv{{\\left[\\frac{{ ({round(eccentricity_vector[0],2)}) ({round(r[0],2)}) + ({round(eccentricity_vector[1],2)}) ({round(r[1],2)}) + ({round(eccentricity_vector[2],2)}) ({round(r[2],2)}) }}{{({round(eccentricity_magnitude,2)})({round(r_magnitude,2)})}}\\right]}} \\Rightarrow  {round(true_anomaly,2)}^ \\circ $$ \n" +'\n')
        
        # Calculate the semi major axis
        f.write("$$\\text{Having found the six orbital elements, we can go on to compute other parameters}$$\n"+'\n')
        f.write("$$\\text{The perigee and apogee radii are}$$")
        f.write(f"$$ {{r}}_p = \\frac{{{{h}}^2}}{{\\mu}} \\frac{{1}}{{1 + {{e}} \\cos (0^\\circ)}} = \\left( \\frac {{{round(specific_angular_momentum_magnitude,2)}^2}} {{398600}} \\right) \\left( \\frac {{1}} {{1 + {round(eccentricity_magnitude,2)} * {math.cos(0*math.pi/180)}}} \\right) = {round(perigee_radii,2)} \\texttt{{ (Km)}} $$")
        f.write(f"$$ {{r}}_a = \\frac{{{{h}}^2}}{{\\mu}} \\frac{{1}}{{1 + {{e}} \\cos (180^\\circ)}} = \\left( \\frac {{{round(specific_angular_momentum_magnitude,2)}^2}} {{398600}} \\right) \\left( \\frac {{1}} {{1 + {round(eccentricity_magnitude,2)} * {math.cos(180*math.pi/180)}}} \\right) = {round(apogee_radii,2)} \\texttt{{ (Km)}} $$")
        f.write("$$ \\text{From these it follows that the semimajor axis of the ellipse is}$$ \n"+'\n')
        f.write(f"$${{a}} = \\frac{{1}}{{2}} \\left( {{r}}_p + {{r}}_a \\right) = {round(semimajor_axis,2)} \\texttt{{ (Km)}} $$  \n"+'\n')
        
        # The end of the pdf
        f.write("$\\end{math}"+'\n')        
        f.write("\\end{document}"+'\n')       
print_the_algorithm_by_latex()