## Cooding the algorithm of GIBBS method using PYTHON
> supervisor : Dr . Ahmed badwy
>
> By : Ahmed khaled fathy

## Requirments 
- **[VS Code](https://code.visualstudio.com/)** or **[VSCodium](https://vscodium.com/)**
- **[LaTeX Workshop](https://marketplace.visualstudio.com/items?itemName=James-Yu.latex-workshop)** extension 
- TeX distribution such as **[MiKTeX](https://miktex.org/download)**

## instructure 
1. Debug the python code **use the example or any right vectors didn't make any error messages yet**
2. open the main.tex file -If u installed the requirments right you will find a TEX icon on the VS code-
3. Build the main.tex file and you will get a main.pdf wiht all the algorethm to get the v vector 

## Example 
> if you used the numpers on the refrenc algoritgm
>
> -294.32i + 4265.1j + 5956.7k
>
> -1365.5i + 3637.6j + 6346.8k
>
> -2940.3i + 2473.7j + 6555.8k
.
> you will get that pdf 


Vectors:
$$\vec{r}_1 = -294.32 \hat{i} + 4265.1 \hat{j} + 5986.7 \hat{k} \texttt{ Km}$$
$$\vec{r}_2 = -1365.5 \hat{i} + 3637.6 \hat{j} + 6346.8 \hat{k} \texttt{ Km}$$
$$\vec{r}_3 = -2940.3 \hat{i} + 2473.7 \hat{j} + 6555.8 \hat{k} \texttt{ Km}$$

Magnitude of the vector

$${r}_1 = \sqrt {-294.32^2 +4265.1^2 + 5986.7^2 }  = 7356.513 \texttt{ Km}$$

$${r}_2 = \sqrt {-1365.5^2 +3637.6^2 + 6346.8^2 }  = 7441.68 \texttt{ Km}$$

$${r}_3 = \sqrt {-2940.3^2 +2473.7^2 + 6555.8^2 }  = 7598.886 \texttt{ Km}$$


$\textbf{Step2:Find } \vec{C}_{12}$ 


$$\vec{C}_{12}=  \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ 
-294.32 & 4265.1 & 5986.7 \\
-1365.5 & 3637.6 & 6346.8 \end{vmatrix}= 5292516.76i + -6306848.67j + 4753375.62k (\texttt{Km}^2)$$ 

$$\vec{C}_{23}=  \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ 
-1365.5 & 3637.6 & 6346.8 \\
-2940.3 & 2473.7 & 6555.8 \end{vmatrix}= 8147298.92i + -9709551.14j + 7317797.93k (\texttt{Km}^2)$$ 

$$\vec{C}_{31}=  \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\
-2940.3 & 2473.7 & 6555.8 \\
-294.32 & 4265.1 & 5986.7 \end{vmatrix}= -13151842.79i + 15673190.95j + -11812614.15k (\texttt{Km}^2)$$ 

$\textbf{Step3:The magnitude of} \vec{C}_{23}:$

$${C}_{23} = \sqrt {8147298.92^2 +-9709551.14^2 + 7317797.93^2 } = 14635710.764(\texttt{Km}^2) $$

$\textbf{The unit vector of } \vec{C}_{23} : $
 
$$\hat{{C}_{23}} =\left(\frac{ 8147298.92 i + -9709551.14 j + 7317797.93 k }{ 14635710.764 }\right) = 0.557i + -0.663j + 0.5k $$

$\textbf{The unit vector of } \vec{r}_1 : $ 


$$\hat{{u}_{r_1}}= \left(\frac{ -294.32 i + 4265.1 j + 5986.7 k }{ 7356.513 }\right) = -0.04i + 0.58j + 0.814k $$

$\textbf{Therefore}:$
 
$$ \\texttt{the dot product of two unit vectors = } -6.118058187509767e-06 $$
 


$$ \text{This is close enough to or equal zero for our purposes. }$$

$$\text{The three vectors ${r}__{1}$, ${r}__{2}$, and ${r}_{3}$ are coplanar.}$$


$\textbf{Step4}:$


 
```math
\begin{align*}

\vec{N} &= r_1 \vec{C}_{23} + r_2 \vec{C}_{31} + r_3 \vec{C}_{12} \\

  &= 7356.51 ( 8147298.92i + -9709551.14j + 7317797.93k) \\
  
  &+ 7441.68 ( -13151842.79i + 15673190.95j + -11812614.15k) \\
  
  &+ 7598.89 ( 5292516.76i + -6306848.67j + 4753375.62k)  \\
  
  \vec{N} &=(2281140568.33i + -2718596493.52j + 2048144274.9k)  (\texttt{Km}^3) \\

  \text{N} &= \sqrt {2281140568.33^2 +-2718596493.52^2 + 2048144274.9^2 } \\
  
  &= 4097470458.447 (\texttt{Km}^3) 
  
\end{align*}
```
```math
\begin{align*}

  \vec{D} &= \vec{C}_{12} + \vec{C}_{23} + \vec{C}_{31} \\

&=  ( 5292516.76i + -6306848.67j + 4753375.62k)  \\

&+  ( -13151842.79i + 15673190.95j + -11812614.15k)  \\

&+  ( 8147298.92i + -9709551.14j + 7317797.93k)  \\

\vec{D} &=( 287972.89i + -343208.86j + 258559.4k)  (\texttt{Km}^2)\\

\text{D} &= \sqrt {287972.89^2 +-343208.86^2 + 258559.4^2 } \\

&= 517275.237  (\texttt{Km}^2)

\end{align*}
```
```math
\begin{align*}
\vec{S} &= (\vec{r}_1)({r}_2 - {r}_3) + (\vec{r}_2)({r}_3 - {r}_1) + (\vec{r}_3)({r}_1 - {r}_2) \\
 
&=  ( -294.32i + 4265.1j + 5986.7k) (7441.68 - 7598.886) \\

&=  ( -1365.5i + 3637.6j + 6346.8k) (7598.886 - 7356.513) \\

&=  ( -2940.3i + 2473.7j + 6555.8k) (7356.513 - 7441.68) \\

\vec{S} &=( -34275.77i + 478.57j + 38810.21k)  (\texttt{Km}^2)

 \end{align*}

```

$\textbf{Step5}: $

```math
 
$$\vec{v}_2 = \sqrt{\left(\frac {\mu}{{N}{D}}\right) } \left(\frac{\vec{D} \times \vec{\text{r}}_{2}}  {{\text{r}}_{2}} + \vec{S}\right)$$ 
```
```math
\begin{align*} \vec{v}_2 &= \sqrt{\frac{398,600}{ (4097470458.45 )(517275.24)} } \\
  &\times \left[\frac{\begin{vmatrix}
 $i$ & $j$ & $k$ \\
 $-294.32$ & $4265.1$ & $5986.7$ \\
 $-1365.5$ & $3637.6$ & $6346.8$ \\
  \end{vmatrix}}{7441.68 } +( -34275.77i + 478.57j + 38810.21k) \right ]\\&=( -6.22i + -4.01j + 1.6k) \texttt{ (Km/s)} \end{align*} 
```



$\textbf{Calculate the orbital elements}$

$\text{we have }\vec{r} \text{ and } \vec{v}$ 

$$\vec{r} = (-1365.5) \hat{i} + (3637.6) \hat{j} + (6346.8) \hat{k} \texttt{ (Km)}$$ 

$$\vec{v} = (-6.22) \hat{i} + (-4.01) \hat{j} + (1.6) \hat{k} \texttt{ (Km/s)}$$ 

$\textbf{Step1}$

$${r} = \sqrt{\vec{r} \text{ . } \vec{r}} = \sqrt{((-1365.5))^2 + ((3637.6))^2 + ((6346.8))^2} = 7441.68 \texttt{ (Km)} $$ 

$\textbf{Step2} $

$${v} = \sqrt{\vec{v} \text{ . } \vec{v}} = \sqrt{((-6.22))^2 + ((-4.01))^2 + ((1.6))^2} = 7.57 \texttt{ (Km/s)} $$ 

$\textbf{Step3} $

$${v}_r = \frac{\vec{v}.\vec{r}}{r} = \frac{((-6.22)).((-1365.5)) + ((-4.01)).((3637.6)) + ((1.6)).((6346.8)) }{7441.68} = 0.54 \texttt{ (Km/s)}$$ 

$$ \text{Since }{v}_r >\text{ 0 , the satellite is flying away from perigee.} $$

$\textbf{Step4} $

$$ \vec{h} = \vec{r} \times \vec{v} = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\
(-1365.5) & (3637.6) & (6346.8) \\
(-6.22) & (-4.01) & (1.6) \end{vmatrix} = 31280.88 \hat{i} + -37277.19 \hat{j} + 28095.03 \hat{k} \texttt{ (Km$^2$/s)} $$

$\textbf{Step5}$

$${h} = \sqrt{\vec{h} \text{ . } \vec{h}} = \sqrt{(31280.88)^2 + (-37277.19)^2 + (28095.03)^2} = 56190.86 \texttt{ (Km$^2$/s)} $$ 

$\textbf{Step6}$

$$ {i} = \cos^{-1}{\frac{{h}_z}{h}} = \cos^{-1}{\left(\frac{28095.03}{56190.86}\right)} \Rightarrow  60.0^ \circ $$ 

$$ \text{Since 0} > {i} < \text{ 90 , this is prograde orbit.} $$

$\textbf{Step7}$ 

$$ \vec{N} = \vec{\hat{K}} \times \vec{h} =  \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\
 0 & 0 & 1 \\
  31280.88 & -37277.19 & 28095.03 \end{vmatrix} = 37277.19 \hat{i} + 31280.88 \hat{j} + -0.0 \hat{k} \texttt{ (Km$^2$/s)} $$

$\textbf{Step8}$

$${N} = \sqrt{\vec{N} \text{ . } \vec{N}} = \sqrt{(37277.19)^2 + (31280.88)^2 + (-0.0)^2} = 48662.95 \texttt{ (Km$^2$/s)} $$ 

$\textbf{Step9}$

$$ {\Omega} = \cos^{-1}{\frac{{N}_x}{N}} = \cos^{-1}{\left(\frac{31280.88}{48662.95}\right)} \Rightarrow  40.0^ \circ $$ 

$\textbf{Step10}$ 

$$\vec{e} = \frac{1}{\mu} \left[ \left( {v}^2 - \frac{\mu}{r} \right) \vec{r} - {r} {v}_r \vec{v} \right] $$

```math
\begin{align*}\vec{e} &= \frac{1}{398600} \left[ \left( {7.57}^2 - \frac{398600}{7441.68} \right)\right]  ((-1365.5) \hat{i} + (3637.6) \hat{j} + (6346.8) \hat{k}) \\
&- 7441.68 * 0.54 \text{  } ((-6.22) \hat{i} + (-4.01) \hat{j} + (1.6) \hat{k}) \\
\vec{e} &= (0.05 \hat{i} + 0.07 \hat{j} + 0.04 \hat{k}) \end{align*}
```

$\textbf{Step11}$ 

$${e} = \sqrt{\vec{e} \text{ . } \vec{e}} = \sqrt{(0.05)^2 + (0.07)^2 + (0.04)^2} = 0.1 $$ 

$\text{Clearly, the orbit is an ellipse.} $ 

$\textbf{Step12}$ 

$$\text{We know that }{e}_z > \text{ 0 ; therefore } \omega \text{ must lie in the first quadrant}$$

$$ {\omega} = \cos^{-1}{\frac{\vec{N}.\vec{e}}{N e}} = \cos^{-1}{\left[\frac{ (37277.19) (0.05) + (31280.88) (0.07) + (-0.0) (0.04) }{(48662.95)(0.1)}\right]} \Rightarrow  30.07^ \circ $$ 

$\textbf{Step13}$

$$\text{We know that }{v}_r > \text{ 0 , which means 0° } \geq \theta < \text{ 180°. Therefore,}$$

$$ {\theta} = \cos^{-1}{\frac{\vec{e}.\vec{r}}{e r}} = \cos^{-1}{\left[\frac{ (0.05) (-1365.5) + (0.07) (3637.6) + (0.04) (6346.8) }{(0.1)(7441.68)}\right]} \Rightarrow  49.93^ \circ $$ 

---


$$\text{Having found the six orbital elements, we can go on to compute other parameters}$$

$$\text{The perigee and apogee radii are}$$

$$ {r}_p = \frac{{h}^2}{\mu} \frac{1}{1 + {e} \cos (0^\circ)} = \left( \frac {56190.86^2} {398600} \right) \left( \frac {1} {1 + 0.1 * 1.0} \right) = 7200.46 \texttt{ (Km)} $$

$$ {r}_a = \frac{{h}^2}{\mu} \frac{1}{1 + {e} \cos (180^\circ)} = \left( \frac {56190.86^2} {398600} \right) \left( \frac {1} {1 + 0.1 * -1.0} \right) = 8802.41 \texttt{ (Km)} $$

$$\text{From these it follows that the semimajor axis of the ellipse is}$$ 

$${a} = \frac{1}{2} \left( {r}_p + {r}_a \right) = 8001.44 \texttt{ (Km)} $$  

