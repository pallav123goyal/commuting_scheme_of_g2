t0 = cpuTime();

p = 32003;
kk = ZZ/p;

sigma = {10,14,6,7,2,13,12,11,3,4,1,9,8,5};
inverseSigma = apply(toList(1..14),k -> position(sigma,i -> i==k));
c = {13641,22825,27432,6173,21845,18827,20908,8404,5290,9258,14370,10545,3220,26424,1,1};

R = kk[u_1..u_12];
x = {sub(c_(14),R),u_1,u_2,u_3,sub(c_(15),R),u_4,u_5,u_6,u_7,u_8,u_9,u_10,u_11,u_12};
y = apply(inverseSigma,i -> sub(c#i,R)-x_i);

w = (i,j) -> x_(i-1)*y_(j-1)-x_(j-1)*y_(i-1);
delta = 192*w(1,2);
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
I = ideal q;
G = groebnerBasis(I, Strategy => "F4");
forceGB G;
J = ideal G;

d = degree J;
print "dim A_F(c):";
print d;

Gdel = groebnerBasis(J + ideal(delta), Strategy => "F4");
if Gdel == matrix {{1_R}} then print "delta is a unit of A(c)";

t1 = cpuTime();
print "Total time taken (in seconds)";
print(t1-t0); -- Expected 25 seconds