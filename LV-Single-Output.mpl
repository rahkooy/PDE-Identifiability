#====================Lotka Volterra Single Output======

restart:

printf("\n================================================\n");
printf("=============Lotka Volterra Single Output======\n\n");
printf("================================================\n\n");

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):
with(VectorCalculus):
with(LinearAlgebra):
with(PDETools):
with(RandomTools):
read "ComputeIdentifiableFunctions.mpl":


#Setting up the PDE:
printf("The PDE system is:\n\n");
sys := [
         diff(u(x, t), t) - d1 * diff(u(x, t), x$2) - u(x, t) * (a1 - b1 * u(x, t) - c1 * v(x, t)),
         diff(v(x, t), t) - d2 * diff(v(x, t), x$2) - v(x, t) * (a2 - b2 * u(x, t) - c2 * v(x, t)),
         y1(x, t) - u(x, t)         
       ];



#Setting up Differential Ring with single output y1
R := DifferentialRing(blocks = [[u, v], [y1]], derivations = [t,x], arbitrary = [a1, b1, c1, d1, a2, b2, c2, d2]):
RG := RosenfeldGroebner(sys, R):

\printf("The output equations are (Note: ignoring input Eq.s):\n\n");
IOeqs:= Equations(RG[1]):

\printf("Out of 3 equations, we are interested in the last one does not depend on u, v which are the state variables::\n\n");
IOeqs[3]:

\printf("since the equation has no monomial with coeff 1, we divide it by c1*d1, the coefficient of the first monomial:\n\n");
IOeqs3:=expand(IOeqs[3]/(c1*d1)):

\printf("Below we collect coefficients of the parameters in IOeqs3, that is the monomials: \n\n");
C3 := [coeffs(IOeqs3,[a1, b1, c1, d1, a2, b2, c2, d2], 'l1')]:l1:

#\printf("We exclude the first monomials whose coeff is 1:\n\n");
M3:=C3[2..nops(C3)]:

\printf("The Wronskians wrt x and t are the following (We use the Wronskian wrt x):\n\n");
W3 := Wronskian(M3, x):
W3t:=Wronskian(M3,t):

#\printf("Redefining differential ring:\n\n");
R2 := DifferentialRing(blocks = [[y1, u, v]], derivations = [t,x], arbitrary = [a1, b1, c1, d1, a2, b2, c2, d2]):
RG2 := RosenfeldGroebner(sys, R2):
Equations(RosenfeldGroebner(sys, R2));

\printf("Determinant of the Wronskian did not finish in reasonable time ");
#Determinant(W3):

\printf("So we try alternatively decomposing IQeqs3 into polynomials with less variables/coefficients, using ComputeIdentifiableFunctions.mpl as follows:\n
The output is two sets: \n
A: a set of polynomials in parameters, \n
B: a set of polynomials in monomials, such that,\n
IOeqs3=sum_i A[i]*B[i]. \n\n");
infolevel:=1:
IOeqs3_ders := indets(IOeqs[3])[11..-1]; nops(IOeqs3_ders); pars:={a1, b1, c1, d1, a2, b2, c2, d2}:
IOeq3_decomposed:=DecomposePolynomial(IOeqs[3], IOeqs3_ders, pars, infolevel); indets(IOeqs[3]);  [seq(el/c2, el in IOeq3_decomposed[1])]:
\printf("Decomposed version has 14 elements instead of 17 \n
We check identifiability of coefficients of these 14 elements using AllIdentifiable package:\n\n");
FilterGenerators(FieldToIdeal([seq(el/c2, el in IOeq3_decomposed[1])]));

\printf("In oder to be able to apply our Wronskian method, we need to normalise the decomposition, i.e., make the lead term one. \n
We divide by the 7th coeff in the decomposition, because the 7th polynomial in the decomposition is quite large.:");
IOeq3_decomposed[1][7]:
IOeq3_decomposed[2][7]:

\printf("So the 7th element has coeff 1, so we exclude it\n");
IOeq3normal:= subsop(7=NULL, IOeq3_decomposed[2]):

\printf("And then compute the Wronskian of the latter with 13 terms:\n\n");
W3d := Wronskian(IOeq3normal, x):

\printf("The determinant of the latter ran out of memory after 10GB and 2900Sec\n\n");
#NF3d := NormalForm(simplify(Determinant(W3d)), RG2):

\printf("So we first compute the normal forms of the entries of the Wronskian (It takes about 2-3 hours)\n
then in order to be able to apply Initial Conditions, we eliminate derivatives of y1 wrt to t. \n
For the latter, we use IOeq[2] to replace derivative of y1 wrt t with the rest of IOeq[2] and do recursively for the higher order derivatives of  y1:\n
IOeq[2] is:\n\n");
IOeqs[2]:

\printf("So pp:=derivative of y1 wrt t is replaced with below differntial polynomial dert:\n\n");
dert:=IOeqs[2] - diff(y1(x, t), t);
pp:=diff(y1(x,t),t):

\printf("The procedure subsdty1 takes as input a differential polynomail p then eliminates derivatives of y1 wrt t of any order using dert \n\n");
subsdty1 := proc(p) 
local dert_local, pp, pnew, ord, i;
dert_local := v(x, t)*y1(x, t)*c1 - diff(y1(x, t), x, x)*d1 + y1(x, t)^2*b1 - y1(x, t)*a1;
pp := diff(y1(x,t),t);
pnew :=p:
ord := difforder(p,t);
for i from 1 to ord do
    # print("i-th order", i);   
    pnew:= simplify(subs(pp=dert_local,pnew));
    # print("pnew",pnew);
od;
return simplify(pnew);
end proc:

\printf("The procedure subsdtu is exactly the same as subsdty1, except that it eliminates derivatives of u intead of y1 \n\n");
subsdtu := proc(p) 
local dert_local, pp, pnew, ord, i;
dert_local := v(x, t)*u(x, t)*c1 - diff(u(x, t), x, x)*d1 + u(x, t)^2*b1 - u(x, t)*a1;
pp := diff(u(x,t),t);
pnew :=p:
ord := difforder(p,t);
#print(ord);
for i from 1 to ord do
    #print("i-th order", i);   
    pnew:= simplify(subs(pp=dert_local,pnew));
    #print("pnew",pnew);
od;
return simplify(pnew);
end proc:

\printf("Assigning random values between 0.1 and 10 to parameters: \n\n");
randvec:=vector(8,Generate(rational(range=1/10..10))):
paramval:=[a1=randvec[1],a2=randvec[2],b1=randvec[3],b2=randvec[4],
c1=randvec[5],c2=randvec[6],d1=randvec[7],d2=randvec[8]];
#param:=[a1,a2,b1,b2,c1,c2,d1,d2];

\printf("evaluate random values of par in Wronskian\n\n");
W3drand:=eval(W3d,paramval):
printf("entry 1,1 of wronskian w random par value",W3drand[1,1]);

\printf("Redefine diff pol system with random values for  parameters\n\n");
sysrand:= eval(sys,paramval);

\printf("Redefine Diff Ring with random values for  parameters");
R3:=DifferentialRing(blocks = [[y1, u, v]], derivations = [t,x]);
RG3:=RosenfeldGroebner(sysrand,R3);
Equations(RG3);


#dimm:=Dimension(W3drand)[1]:
dimm:=13;

\prinf("Assign random values for parmeters in the substitution procedure as well:");
subsdturand := proc(p) 
local dert_local, pp, pnew, ord, i;
dert_local := eval(v(x, t)*u(x, t)*c1 - diff(u(x, t), x, x)*d1 + u(x, t)^2*b1 - u(x, t)*a1,paramval);
pp := diff(u(x,t),t);
pnew :=p:
ord := difforder(p,t);
#print(ord);
for i from 1 to ord do
    #print("i-th order", i);   
    pnew:= simplify(subs(pp=dert_local,pnew));
    #print("pnew",pnew);
od;
return simplify(pnew);
end proc:

\printf("Next we form a matrix whose entries are the normal form of entries of Wronskian in which derivatives of u wrt t are eliminated. \n\n");
#W3dred := Matrix(dimm):
#for i from 1 to dimm do
# for j from 1 to dimm do
#   print("Computing entry", i,j);
#   #print("Wronskian:", W3drand[i,j]);
#   nf := simplify(NormalForm(W3drand[i,j], RG3)):
#   #print("Normal form of the entries of Wronskian:",nf):
#   W3dred[i,j] := op(subsdturand(nf)): 
#   #print("subst of NF",W3dred[i,j]);
# od;
#od;
\printf("Reusing saved matrix of normal forms instead of recalculating:");
read "NormalFormMatrix":

\printf("Initial Conditions, and with parameters evaluated are:");
ICst:=[u(x,t)=simplify(((a1/b1)*exp(-a*x))/(1+exp(-a*x))), 
v(x,t)=simplify(((a2/b2)*exp(av*x))/(1+exp(-av*x)))]:
ICst1:=eval(ICst,[a=1,av=2]):
ICstrand:=eval(ICst1,paramval);

\printf("Procedure for computing list of numerical values of derivatives of u and v wrt x at a given x=sss:");
xdereval := proc(sss) 
local i,xder;
xder:=list(lhs(ICstrand[1])=evalf(subs(x=sss,rhs(ICstrand[1]))),lhs(ICstrand[2])=evalf(subs(x=sss,rhs(ICstrand[2])))); 
for i from 1 to 30 do
 xder:=[op(xder),diff(lhs(ICstrand[1]),x$i)=evalf(subs(x=sss,diff(rhs(ICstrand[1]),x$i))), diff(lhs(ICstrand[2]),x$i)=evalf(subs(x=sss,diff(rhs(ICstrand[2]),x$i)))]:
od;
return xder;
end proc:
#print("testing at x=3");
#xder3:=xdereval(3);


\printf("Method 1: symbolic computation of assigning IC then evaluation at x=3");
#\printf("Apply IC to the entries of the matrix with x=3");
#W3dredic := Matrix(dimm):
#tstart1:=time():
#for i from 1 to dimm do 
# for j from 1 to dimm do 
#  \print("\n \n \n \n");
#  \print("apply IC to entry", i,j);
#  #print("after Normal Form, before Substitution:",W3dred[i,j]);
#  subsexp:=simplify(subs(ICstrand,W3dred[i,j])):
#  #print(subsexp);
#  print("Evaluating numerivally at x=3:");
#  W3dredic[i,j] := simplify(evalf(subs(x=3,subsexp)));
#  #W3dredic[i,j] := evalf(simplify(subs(x=3,subsexp)));
#  #print("After eval:");
#  #print(i,j);
#  #print(W3dredic[i,j]);
# od;
#od;
#timemethod1:=time()-tstart1;

xder3 := xdereval(3):
xder1 := xdereval(1):
W3drediccopy := Matrix(dimm):
tstart2:=time():
for i from 1 to dimm do 
 for j from 1 to dimm do 
  #\print("\n \n \n \n");
  #\print("apply fast substitution method to entry", i,j);
  W3drediccopy[i,j] := eval(subs(xder1,W3dred[i,j]));
  #print(W3drediccopy[i,j]);
 od;
od;
\printf("Time spent for Method 2 (first evaluate derivatives, then substitute in the matrix):");
timemethod2:=time()-tstart2;

#for i from 1 to dimm do
#  for j from 1 to dimm do
#    print(i,j);
#    print(W3dredic[i,j], W3drediccopy[i,j]);
#  od;
#od;

\printf("Determinant:");
#deter:=Determinant(evalf(W3dredic));
#deter:=Determinant(W3dredic);
deter:=Determinant(W3drediccopy);

\printf("The Case of the parameters in the IC being zero:");
ICst0:=eval(ICst,[a=0,av=0]):
ICstrand0:=eval(ICst0,paramval);
xdereval0 := proc(sss) 
local i,xder;
xder:=list(lhs(ICstrand0[1])=evalf(subs(x=sss,rhs(ICstrand0[1]))),lhs(ICstrand0[2])=evalf(subs(x=sss,rhs(ICstrand0[2])))); 
for i from 1 to 30 do
 xder:=[op(xder),diff(lhs(ICstrand0[1]),x$i)=evalf(subs(x=sss,diff(rhs(ICstrand0[1]),x$i))), diff(lhs(ICstrand0[2]),x$i)=evalf(subs(x=sss,diff(rhs(ICstrand0[2]),x$i)))]:
od;
return xder;
end proc:

xder3 := xdereval0(3):
xder1 := xdereval0(1):
W3drediccopy0 := Matrix(dimm):
for i from 1 to dimm do 
 for j from 1 to dimm do 
  #\print("\n \n \n \n");
  #\print("apply fast substitution method to entry", i,j);
  #W3drediccopy0[i,j] := eval(subs(xder3,W3dred[i,j]));
  W3drediccopy0[i,j] := eval(subs(xder1,W3dred[i,j]));
  #print(W3drediccopy0[i,j]);
 od;
od;

\printf("Determinant:");
deter:=Determinant(W3drediccopy0);
