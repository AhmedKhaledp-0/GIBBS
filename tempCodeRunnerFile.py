 f.write("$\\vec{C}_{31} = $\n\n")
    f.write("$\\begin{vmatrix}$\n")
    f.write("$$i$ & $j$ & $k$ \\\\\n")
    f.write(f"${r3[0]}$ & ${r3[1]}$ & ${r3[2]}$ \\\\\n")
    f.write(f"${r3[0]}$ & ${r3[1]}$ & ${r3[2]}$ \\\\\n")
    f.write("$$\\end{vmatrix}\n")
    f.write(f"$= {round(r31Coefficient[0], 2)}i + {round(r31Coefficient[1], 2)}j + {round(r31Coefficient[2], 2)}k$ (\\text{{Km}}^2)$ \n" + '\n \n')