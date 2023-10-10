#====================Fisher's Equation======

restart:

printf("\n================================================\n");
printf("=============Fisher's Equation======\n\n");
printf("================================================\n\n");

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):
with(VectorCalculus):
with(LinearAlgebra):

#Setting up the PDE:
printf("The PDE system is:\n\n");
eq := simplify(diff(n(x,t),t)-d*diff(n(x,t),x$2) - r*n(x,t)*(1-1/k*n(x,t)));

\printf("Note that if below coefficients are identifiable, then the parameters are identifiable.\n\n");
d, r/k, 1/k;

\printf("Coefficients of the parameters except for 1 are:\n\n");
M1 := simplify(n(x, t)^2);
M2 := simplify(diff(n(x, t), x, x));
M3 := simplify(n(x, t));

#\printf("Setting up differential algebra:\n\n");
R := DifferentialRing(blocks = [n], derivations = [x,t], arbitrary = [d,r,k]): 
RG := RosenfeldGroebner([eq], R): Equations(RG):

#Computing normal form of the determinant of the Wronskian of the monomials
\printf("The Determinant of the Wronskian wrt x is:\n\n");
Wr := simplify(Determinant(Wronskian([M1,M2,M3],x)));
#Normal form of Wronskian:
numer(simplify(NormalForm(Wr, RG))):
#\printf("The coefficients of the normal form of the determiannt of Wronskian are:\n\n");
cf := coeffs(expand(numer(simplify(NormalForm(Wr, RG)))[1]),[k, r, d]):
for i from 1 to nops([cf]) do print(cf[i]); print(\n) od:
\printf("Rosenfeld-Groebner of coefficients is:\n\n");
Equations(RosenfeldGroebner([cf[-1],cf[1..4]],R));

\printf("which results in the following cases:\n\n");
\printf("Case 1.\n\n");
diff(n(x, t), x)=0;
\printf("In this case the PDE become the folling ODE:\n\n");
eqt := simplify(diff(n(t),t))- r*n(t) + (r/k)*(n(t))^2;
\printf("According to the initial condition for t equal to 0, the equation becomes an exponential function of alpha times x. So the equation is independent of x iff alpha is zero, which is not possible according to our assumption.\n\n");

\printf("Case 2.\n\n");
diff(n(x, t), t)=0;
\printf("Then the PDE becomes:\n\n");
eqx := simplify(-d*diff(n(x),x$2) - r*n(x)+(r/k)*(n(x))^2); 

\printf("But according to the IC, n is only a function of x if t is zero, which is the following exponential function:\n\n");
IC:=n(x)=simplify((k*exp(-a*x))/(1+exp(-a*x)));

\printf("which should satisfy our ODE, i.e., plugging in IC to the ODE, we must obtain zero. So below equation must be zero:\n\n");
eqnx:=simplify(eval(eqx,IC));

\printf("the above equation is zero if and only if its numerator is zero, which is the case iff either k is 0 (which is not possible according to our assumptions), or:\n\n");
(a^2*d - r)*exp^(-a*x) - a^2*d - r = 0;
\printf(" This can be zero only if the coeff of the exponential term is zero and the constant term is zero. This leads to the equations:\n\n");
a^2*d - r = 0;
-a^2*d - r = 0;
\printf("which leads to:\n\n"):
2 * a^2 * d = 0;
printf("But this cannot happen as a and d are nonzero. Hence, none of the above cases can happen and therefore, the coefficients are identifiable.");