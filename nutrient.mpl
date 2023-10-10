#====================Nutrient Reduction–Diffusion on a Fixed Domain======

restart:

printf("\n================================================\n");
printf("=============Nutrient Reduction–Diffusion on a Fixed Domain======\n\n");
printf("================================================\n\n");

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):
with(VectorCalculus):
with(LinearAlgebra):

#Setting up the PDE:
printf("The PDE system is:\n\n");
eq := diff(c(x,t),t)-d*diff(c(x,t),x$2)-lambda*c(x,t)/(c_0+c(x,t));

#Clearing out the numerators. This is our input-output equation (called IO-equation below):
printf("The Input-Output equations are:\n \n");
numer(eq);

#Collecting monomials of IO-eq whose coeff is not 1:
\printf("Collecting monomials of IO-eq whose coeff is not 1:\n\n");
M1 := simplify(diff(c(x, t), x, x)*c(x, t));
M2 := simplify(diff(c(x, t), x, x));
M3 := simplify(diff(c(x, t), t));
M4 := c(x,t);

\printf("The Wronskian wrt x is:\n\n");
Wr := simplify(Determinant(Wronskian([M1,M2,M3,M4],x)));

#Setting up differential ring, RosenfelGroebner:
R := DifferentialRing(blocks = [c], derivations = [x,t], arbitrary = [lambda, c_0, d]):
RG := RosenfeldGroebner([eq], R):

#Checking if the Wronskian vanishes: Collect coeff, compute normal form, then RG, check if RG=0:
cf := coeffs(expand(numer(simplify(NormalForm(Wr, RG)))[1]),[lambda, c_0, d]):
#for i from 1 to nops([cf]) do print(cf[i]); print(\n) od;

#Computing Rosenfeld-Groebner of coefficients to check when they vanish (we chose just 10 of those to run the computation faster).
\printf("Computing Rosenfeld-Groebner of coefficients\n\n");
Equations(RosenfeldGroebner([cf[-1],cf[1..10]],R));

#So, the above are out to necessary conditions for the vanishing of the Wronskian. 
#Outside of these conditions, the parameters are identifiable as, by the 
#linear independence of the monomials (due to the non-vanishing of the Wronskian), 
#we can recover the coefficients of the IO-equation, which are d, d * c_0, c_0, and 
#lambda, which shows that we can find d, c_0, and lambda.
\printf("So the Wronskian is nonsingular if the following two cases do not happen:\n\n");
\printf("CASE 1:\n\n"); 
dc/dx=0;
\printf("then:\n\n"); 
d * diff(c(x, t), x, x)*c(x, t)= 0; 
\printf("So c does not depend on x anymore. Therefore equation simplifies to:\n\n");
eq := diff(c(t),t)+lambda*c(t)/(c_0+c(t));
\printf("Using boundary condition c(R,t)=1, we obtain that:\n\n");
diff(c(t),t) = 0; 
\printf("Then eq is zero iff:\n\n");
lambda * c(t) / (c_0 + c(t)) = (lambda)/(c_0+1);
\printf("which can be zero only if is lambda is zero, which is impossible according to our assumption\n\n");

\printf("CASE 2\n\n"); 
diff(c(x,t),t)=0;
\printf("Then c is equal to c(x)  and does not depend on t and is only a function of  x. Plugging this into the PDE we obtain the following ODE:\n\n");
eq := -d*diff(c(x),x$2)+lambda*c(x)/(c_0+c(x));
\printf("Using the boundary condition c(x,0) equal 1 and repeating the argument in Case 1, one can see that this case cannot happen");

#====================

