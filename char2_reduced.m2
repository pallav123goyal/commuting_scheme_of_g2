t0 = cpuTime();

p = 2;
kk = ZZ/p;

-- C(gl_2)
R2 = kk[x_1..x_4, y_1..y_4];
X = genericMatrix(R2,x_1,2,2);
Y = genericMatrix(R2,y_1,2,2);

I2 = ideal flatten entries(X*Y - Y*X);
if(radical I2 == I2) then print "C(gl_2) is reduced";

-- Z_(2,2)
R22 = kk[a_1..a_4, b_1..b_4, c_1..c_4, d_1..d_4, l];
X1 = genericMatrix(R22,a_1,2,2);
Y1 = genericMatrix(R22,b_1,2,2);
X2 = genericMatrix(R22,c_1,2,2);
Y2 = genericMatrix(R22,d_1,2,2);

I22 = ideal flatten entries(X1*Y1 - Y1*X1 - l*id_(R22^2)) + ideal flatten entries(X2*Y2 - Y2*X2 - l*id_(R22^2)) + ideal(trace X1 + trace X2, trace Y1 + trace Y2);
if(radical I22 == I22) then print "Z_(2,2) is reduced";

-- C(gl_3)
R3 = kk[x_1..x_9, y_1..y_9];
X = genericMatrix(R3,x_1,3,3);
Y = genericMatrix(R3,y_1,3,3);

I3 = ideal flatten entries(X*Y - Y*X);
if(radical I3 == I3) then print "C(gl_3) is reduced";

t1 = cpuTime();
print "Total time taken (in seconds)";
print(t1-t0); -- Expected 1 second