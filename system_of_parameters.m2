t0 = cpuTime();

p = 32003;
kk = ZZ/p;

sigma = {10,14,6,7,2,13,12,11,3,4,1,9,8,5};

R = kk[x_1..x_14, y_1..y_14];
w = (i,j) -> x_(i)*y_(j)-x_(j)*y_(i);
q = {
	w(3,9)-w(5,11)+2*w(6,12)-w(7,13)+w(8,14),
	w(4,10)-3*w(5,11)+3*w(6,12)-w(7,13)+2*w(8,14),
	2*w(1,3)-w(2,3)+w(5,10)-2*w(6,11)+w(7,12),
	-3*w(1,4)+2*w(2,4)-3*w(5,9)-w(8,13),
	-w(1,5)+w(2,5)+w(3,4)-2*w(6,9)-w(8,12),
	w(1,6)+2*w(3,5)-w(7,9)-w(8,11),
	3*w(1,7)-w(2,7)+3*w(3,6)-w(8,10),
	w(2,8)+w(4,7)-3*w(5,6),
	-2*w(1,9)+w(2,9)-w(4,11)+2*w(5,12)-w(6,13),
	3*w(1,10)-2*w(2,10)+3*w(3,11)+w(7,14),
	w(1,11)-w(2,11)+2*w(3,12)+w(6,14)+w(9,10),
	-w(1,12)+w(3,13)+w(5,14)+2*w(9,11),
	-3*w(1,13)+w(2,13)+w(4,14)+3*w(9,12),
	-w(2,14)+w(10,13)-3*w(11,12)
};
s = i -> x_i - y_(sigma_(i-1));

I = ideal q + ideal apply(1..14, s) + ideal (x_1, x_5);
G = groebnerBasis(I, Strategy => "F4");

forceGB G;
I = ideal G;                    

d = dim I;
print("Dimension of A/(s):");
print d;

H = reduceHilbert hilbertSeries I;
print("Hilbert series:");
print H;

t1 = cpuTime();
print "Total time taken (in seconds)";
print(t1-t0); -- Expected 5 seconds