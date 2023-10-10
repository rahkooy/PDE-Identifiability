#====================A Reaction-Diffusion Model of Cancer Invasion======

restart:

printf("\n================================================\n");
printf("=============A Reaction-Diffusion Model of Cancer Invasion======\n\n");
printf("================================================\n\n");

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):
with(VectorCalculus):
with(LinearAlgebra):

#Setting up the PDE:
printf("The PDE system is:\n\n");
sys := expand([
         diff(u(x,t),t) 
         - r1 * u(x,t) * (1 - u(x,t)/k1 - v(x,t)/k2 * a12) 
         + d1 * w(x,t) * u(x,t),

         diff(v(x,t),t) 
         - r2 * v(x,t) * (1 - v(x,t)/k2 - u(x,t)/k1 * a21) 
         - d2 * diff((1-u(x,t)) * diff(v(x,t),x),x),
 
        diff(w(x,t),t) 
        - d4 * diff(w(x,t),x$2) 
        - r3 * v(x,t) 
        + d3 * w(x,t),

         y1(x, t) - u(x, t),
         y2(x, t) - v(x, t),
         y3(x, t) - w(x, t)
]); 

#Setting up differential algebra:
R := DifferentialRing(blocks = [[u, v,w], [y1, y2,y3]], derivations = [x,t], arbitrary = [r1,k1,a12,k2,d1,r2,a21,d2,d3,r3,d4]):
RG := RosenfeldGroebner(sys,R):
IOeqs := Equations(RG[1])[-3..-1]:
for eq in IOeqs do print(expand(eq)); print() od:

printf("We consider each of the three equations separately  and compute the \\
Wronskian for the coeffs of each equation. \n\n");

printf("We divide the first equation by its coefficients so that, \\
We obtain an equation with a monomial that doesn't contain any parameter. \\
Then we collect the monomials:\n\n");
eq1:=expand(IOeqs[1]/(k1*k2));
Ms1:= coeffs(eq1,[k1,k2,d2,r2,a21], 'l1');
l1;

printf("The above line shows that we also have collected the coefficient of 1. \\
So we exclude that: \n\n");
Ms11:= Ms1[2..nops([Ms1])];

printf("Second equation contains a term without parameter coefficients \\
So we just collect the coefficients\n\n");
Ms2 := coeffs(expand(IOeqs[2]), [r1,k1,a12,k2,d1,r2,a21,d2,d3,r3,d4], 'l2');
l2;
Ms21:=Ms2[2..nops([Ms2])];

printf("We divide the third equation by its coefficients \\
So that we obtain an equation with a monomial that doesn't contain any parameter. \\
Then we collect the monomials.\n\n");
eq3:=expand(IOeqs[3]/(k1*k2));
Ms3 := coeffs( eq3, [k1,k2,r1,a12,d1], 'l3');
l3;
Ms31:=Ms3[2..nops([Ms3])];

#Computing Wronskians: 
Wr1 := Wronskian([Ms11], x):
Wr2 := Wronskian([Ms21], x):
Wr3 := Wronskian([Ms31], x):

#redefining Differential Ring:
R2 := DifferentialRing(blocks = [[y1,y2, y3,u,v,w]], derivations = [x,t], arbitrary = [r1,k1,a12,k2,d1,r2,a21,d2,d3,r3,d4]):
RG2 := RosenfeldGroebner(sys, R2):
Equations(RosenfeldGroebner(sys, R2))[1]:

#Computing normal form of the determinant of Wronskians:
NF1 := NormalForm(simplify(Determinant(Wr1)), RG2):
NF2 := NormalForm(simplify(Determinant(Wr2)), RG2):
NF3 := NormalForm(simplify(Determinant(Wr3)), RG2):

#Coeffs of  normal form of the first Wronskian:
cfs1 := coeffs(expand(numer(NF1))[1],[k1,k2,d2,r2,a21]):
#for cf in cfs1 do print(cf); print() od;

printf("Number of Coefficients of the normal form of the deteminat of the first Wronskian: \n\n");
nops([cfs1]);

printf("Rosenfeld-Groebner of the coefficients of the first equations: \n\n");
Equations(RosenfeldGroebner([cfs1[-1],cfs1[1..10]],R));

#Coeffs of  normal form of the second Wronskian:
cfs2 := coeffs(expand(numer(NF2))[1],[r1,k1,a12,k2,d1,r2,a21,d2,d3,r3,d4]):
#for cf in cfs2 do print(cf); print() od:


printf("Number of Coefficients of the normal form of the deteminat of the second Wronskian:\n\n");
nops([cfs2]);

#redefine differential ring:
Rnew := DifferentialRing(blocks = [[u, v,w,y1, y2, y3]], derivations = [x,t], arbitrary = [r1,k1,a12,k2,d1,r2,a21,d2,d3,r3,d4]):
eqs2 := Equations(RosenfeldGroebner([cfs2][3..3], Rnew)):
printf("Rosenfeld-Groebner of the coefficients of the second Wronskian\n\n");
for eq in eqs2 do print(eq) od;

printf("We take the first condition, solve using the specified initial condition. \\
We obtain only zero solution, hence this case is not possible: \n\n");
nops(eqs2):
subs(w(x,0) =0,eval(diff(w(x, t), t, x)*w(x, t) - diff(w(x, t), x)*diff(w(x, t), t), t=0));

printf("So, either derivative of w wrt t and x is 0 when t equal to 0, \\
Or derivative of w wrt t is zero when t equal to 0.\n\n");

pdsolve([diff(w(x, t), t)]); 
printf("This one with w equal to 0 when t being zero implies w is zero \n\n");

pdsolve([diff(w(x, t), x)]); 
printf("This contradicts the boundary condition \\
As the parameters are not allowed to be 0. \n\n");

#for i from 1 to 3 do Equations(RosenfeldGroebner(systfor2[i], Rnew)) od;


cfs3 := coeffs(expand(numer(NF3))[1],[k1,k2,r1,a12,d1,d2,d3,d4,r2,r3,a21]):
nops([cfs3]):
printf("Rosenfeld-Groebner of the coefficients for the third equation: \n\n");
Equations(RosenfeldGroebner([cfs3][1..10],Rnew));

printf("Adding the above conditions to the system:\n\n");
sysplusux := [op(sys), diff(u(x, t), x)]:
RGplusux := RosenfeldGroebner(sysplusux, R):
eqs := Equations(RGplusux):
for eq in eqs[1] do print(eq) od:

printf("The states are not constant w.r.t. x according to the conditions of the problem \\
(As it would imply several parameters are 0 - not permitted). \\
This is because the normal forms are: \n\n");
NormalForm(diff(w(x,t),x), RGplusux[1]);
NormalForm(diff(u(x,t),x), RGplusux[1]);
NormalForm(diff(v(x,t),x), RGplusux[1]);
