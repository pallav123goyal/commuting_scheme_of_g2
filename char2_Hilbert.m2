t0 = cpuTime();

p = 2;
kk = ZZ/p;

R = kk[x_1..x_14, y_1..y_14];
w = (i,j) -> x_i*y_j - x_j*y_i;

q = {
	w(3,9)+w(5,11)+w(7,13)+w(8,14),
	w(4,10)+w(5,11)+w(6,12)+w(7,13),
	w(2,3)+w(5,10)+w(7,12),
	w(1,4)+w(5,9)+w(8,13),
	w(1,5)+w(2,5)+w(3,4)+w(8,12),
	w(1,6)+w(7,9)+w(8,11),
	w(1,7)+w(2,7)+w(3,6)+w(8,10),
	w(2,8)+w(4,7)+w(5,6),
	w(2,9)+w(4,11)+w(6,13),
	w(1,10)+w(3,11)+w(7,14),
	w(1,11)+w(2,11)+w(6,14)+w(9,10),
	w(1,12)+w(3,13)+w(5,14),
	w(1,13)+w(2,13)+w(4,14)+w(9,12),
	w(2,14)+w(10,13)+w(11,12)};
I = ideal q;

G = groebnerBasis(I, Strategy => "MGB");
forceGB G;
L = ideal G;
print "Dimension of A:";
print dim(L);
print "Hilbert series of A:";
print reduceHilbert hilbertSeries(R/L);

X = matrix{{x_1,x_3,x_5,x_12},{x_9,x_1+x_2,x_4,x_13},{x_11,x_10,x_2,x_14},{x_6,x_7,x_8,0}};
S = R[t];
e2 = sub(coefficient(t^2, det(t*id_(S^4) - sub(X,S))), R);
T = (degreesRing R)_0;
p1 = poincare(R/L);

Ge = gb(I + ideal e2, Hilbert => (1-T^2)*p1);
Le = ideal leadTerm Ge;
print "Hilbert series of A/(e_2):";
print reduceHilbert hilbertSeries(R/Le);

t1 = cpuTime();
print "Total time taken (in seconds)";
print(t1-t0); -- Expected 360 seconds