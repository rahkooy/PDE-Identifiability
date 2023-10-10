#====================Lotka-Volterra Reaction Network======

restart:

printf("\n================================================\n");
printf("=============Lotka-Volterra Reaction Network======\n\n");
printf("================================================\n\n");

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):
with(VectorCalculus):
with(LinearAlgebra):

#Setting up the PDE:
printf("The PDE system is:\n\n");
sys := [
         diff(u(x, t), t) - d1 * diff(u(x, t), x$2) - u(x, t) * (a1 - b1 * u(x, t) - c1 * v(x, t)),
         diff(v(x, t), t) - d2 * diff(v(x, t), x$2) - v(x, t) * (a2 - b2 * u(x, t) - c2 * v(x, t)),
         y1(x, t) - u(x, t),
         y2(x, t) - v(x, t)
       ]; 

#Setting up Differential Ring 
R := DifferentialRing(blocks = [[u, v], [y1, y2]], derivations = [x, t], arbitrary = [a1, b1, c1, d1, a2, b2, c2, d2]):
RG := RosenfeldGroebner(sys, R):

\printf("The output equations are (Note: ignoring input Eq.s):\n\n");
IOeqs := Equations(RG[1])[-2..-1];

#Below we collect coeff(IOeqs[1],a1), ..., coeff(IOeqs[1],d1); and ONLY the linear factors; and NOT coeff(IOeqs[1],1).
Ms1 := map2(coeff, IOeqs[1], [a1, b1, c1, d1], 1):
Ms2 := map2(coeff, IOeqs[2], [a2, b2, c2, d2], 1):

\printf("The Wronskians wrt x are:\n\n");
Wr1 := Wronskian(Ms1, x);
Wr2 := Wronskian(Ms2, x);

#Redefining differential ring removing elimination ranking:
R2 := DifferentialRing(blocks = [[y1, y2, u, v]], derivations = [x, t], arbitrary = [a1, b1, c1, d1, a2, b2, c2, d2]):
RG2 := RosenfeldGroebner(sys, R2):
Equations(RosenfeldGroebner(sys, R2)):

#Computing Normal Forms of the determinats of the Wronskians:
NF1 := NormalForm(simplify(Determinant(Wr1)), RG2):
NF2 := NormalForm(simplify(Determinant(Wr2)), RG2):

#Collecting coeffs of normal forms:
cfs := coeffs(expand(numer(NF1))[1], [a1, b1, c1, d1, a2, b2, c2, d2]):
cfs2 := coeffs(expand(numer(NF2))[1], [a1, b1, c1, d1, a2, b2, c2, d2]):
#for cf in cfs do print(cf); print() od:

\printf("Rosenfeld-Groebner of the coefficints of the determinant of \\
 the wronskian are:\n\n");
Equations(RosenfeldGroebner([cfs[-1],cfs[1..10]],R2));
Equations(RosenfeldGroebner([cfs2[-1],cfs2[1..10]],R2));

\printf("This leads to the following cases:\n\n");
\printf("Case 1. the following three equations:\n\n");
diff(v(x, t), x, t) = 0;
diff(u(x, t), x) = 0;
diff(u(x, t), t) = 0;

\printf("Since derivatives of u wrt both x and t is zero, \\
it means that u is a constant from the boundary conditions, this constant should be:\n\n");
a1/b1 = 0;
\printf("which implies that:\n\n"); 
a1=0; 
\printf("But we assumed that the parameters are nonzero. So this case does not happen.\n\n");

\printf("Case 2:\n\n");
u(x, t) = 0;
\printf("This cannot happen because according to the boundary conditions, \\
limit of u when x goes to infinity is:\n\n");
a1/(b1);

\printf("Case 3:\n\n");
v(x, t) = 0;
\printf("Similar to Case 2\n\n");

\printf("Case 4. the following two qeuations:\n\n");
diff(u(x, t), x) = 0;
diff(v(x, t), x) = 0;
\printf("Similar to Case 1.\n\n");

\printf("Case 5. the following two qeuations:\n\n");
diff(u(x, t), t) = 0;
diff(v(x, t), t) = 0;
\printf("Therefore u and v only depend on x, hence they are given by the initial \\
conditions at t equal to zero, which are exponential functions. Now we plug-in\n\n");
diff(u(x, t), t) = 0;
diff(v(x, t), t) = 0;
\printf("into the PDE system and obtain the following two new equations:\n\n");
equx:= -d1 * diff(u(x), x$2) - u(x) * (a1 - b1 * u(x) - c1 * v(x));
eqvx:= -d2 * diff(v(x), x$2) - v(x) * (a2 - b2 * u(x) - c2 * v(x));

\printf("Initial Conditions are:\n\n");
ICu:=u(x)=simplify(((a1/b1)*exp(-a*x))/(1+exp(-a*x)));
ICv:=v(x)=simplify(((a2/c2)*exp(b*x))/(1+exp(b*x)));

\printf("Plugging in the ICu and ICv into the new PDEs we obtain:\n\n");
S1:=simplify(eval(equx,ICu));
S2:= simplify(eval(S1,ICv));

\printf("The above is zero iff the nominator is zero, which is the case iff either \n\n");
a1*exp(-a*x) = 0;
\printf("or\n\n");
-a2*c1*exp(-x*(2*a - b)) + ((-a^2*d1 + a1)*c2 - 2*a2*c1)*exp(x*(-a + b)) - c2*(a^2*d1 - a1)*exp(-a*x) + (c2*(a^2*d1 + a1) - a2*c1)*exp(b*x) + c2*(a^2*d1 + a1) = 0;

\printf("The first can happen only if a1 is 0, which is not possible; \\
the latter happens if the exponentials are linearly dependent which only can happen \\
if the exponents are the same. Doing case distinction, we should consider all pair of \\
exponents being equal. Consiedring first exponent to be equal to any other exponents\\
 we get 4 cases, each of which are discussed and none can happen. \\
 For*convenience, set: \n\n");
u=-a2 * c1;  
w=-d1 * a^2 + a1;  
t=d1* a^2+a1;
\printf("Then the coefficiens* become\n\n");
[u,   w*c2-2*u,   c2*w,  c2*t-u,  c2*t ];
\printf("Setting:\n\n");
y = e^x;
\printf("leads to the following cases:\n\n");

\printf("Case 5.1\n\n");
-2*a+b=-a+b;
\printf("which implies:");
a=0; 
\printf("which cannot happen.\n\n");


\printf("Case 5.2\n\n");
-2*a+b=-a;
\printf("which implies:");
a=b; 
\printf("The exponents become.\n\n");
y^(-a),y^(0), y^(-a), y^(a), y^(0);
\printf("collecting coeffs of similar exponentials we have \n\n");
u + c2*w = 0; c2*w - 2*u + c2t = 0; c2*t - u = 0;
\printf("which implies\n\n");
2*u = 0;
\printf("which is not possible.\n\n");


\printf("Case 5.3.\n\n");
-2 *a+b = b;
\printf("which implies\n\n");
a=0;
\printf(" which is not possible. \n\n");

\printf("Case 5.4.\n\n");
-2 * a+b=0; 
\printf("which implies \n\n");
b=2* a;
\printf("The exponents become\n\n");
y^0, y^(a), y^(-a), y^(2*a), y^0;
\printf("since the following are linearly independent"\n\n); 
y^(a), y^(2*a) , 
\printf("if their coeffs do not kill each other then we are done. \\
This is the case as the coeffs are \n\n");
c2*w-2*u, c2*w;


\printf("Trying more coefficients, for 22 coefficients, we only obtain four \\
components  out of five obtained using 10 coefficients. For 23 coeff \\
we ran out of memory.\n\n");
#nops([cfs]); nops([cfs2]);
#Equations(RosenfeldGroebner([cfs[-1],cfs[1..15]],R2));
#Equations(RosenfeldGroebner([cfs[-1],cfs[1..20]],R2));
#Equations(RosenfeldGroebner([cfs[-1],cfs[1..22]],R2));
#Equations(RosenfeldGroebner([cfs[-1],cfs[1..23]],R2));