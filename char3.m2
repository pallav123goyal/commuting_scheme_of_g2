t0 = cpuTime();

p = 3;
kk = ZZ/p;

S = kk[x_1..x_14, y_1..y_14];
w = (i,j) -> x_i*y_j - x_j*y_i;
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

L = matrix(ZZ, {
  { 0, 0, 1, 1, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0,-1, 0, 1, 0, 1, 1,-1, 0, 0,-1},
  { 0, 0, 0, 0, 0,-1,-1, 0, 1,-1, 0, 0,-1, 0,-1, 1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0},
  { 0, 1,-1, 0, 0,-1, 0, 0, 1, 1, 1, 0, 0, 0,-1, 1, 0, 1,-1, 1, 1, 1, 1, 0, 1, 0, 1, 0},
  { 1, 1,-1,-1, 0, 0,-1, 0, 0, 0, 1,-1, 0, 1, 0, 0,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 0, 0},
  { 0, 0,-1, 1, 0, 0, 0, 0,-1, 1,-1, 1, 0, 1, 1, 0, 1, 0,-1, 0, 1, 0, 1, 0,-1, 0, 0, 0},
  { 0, 1, 0, 0, 0,-1,-1,-1, 1, 0,-1, 1, 1, 1,-1,-1,-1,-1, 1, 1,-1, 0, 0, 0,-1, 0, 0,-1},
  { 0,-1,-1, 0, 0,-1, 0, 1,-1, 0, 0, 0, 1, 1, 0,-1, 0, 0, 0, 0, 0,-1,-1,-1,-1, 1, 0, 1},
  { 1,-1, 0, 0,-1, 0, 0, 1, 1,-1,-1, 0,-1, 0, 0, 1, 0,-1, 0, 0,-1, 0, 0, 1, 0, 1, 0, 0},
  { 0, 0, 1, 1, 1,-1, 1, 1, 0, 0,-1, 0, 0,-1, 0, 0,-1, 1, 0,-1, 1,-1,-1, 0,-1, 0, 0,-1},
  { 0,-1, 1,-1,-1,-1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1,-1, 0,-1,-1,-1,-1, 1, 0, 1, 0, 1},
  { 0,-1, 1, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0,-1, 1,-1, 1, 1,-1,-1,-1, 0, 1, 1},
  { 0, 1, 0, 1,-1, 1, 1, 0,-1, 0, 0, 0, 1,-1, 1,-1, 1,-1,-1, 0, 0,-1, 0, 0,-1,-1, 1, 1},
  { 1,-1, 1,-1, 1, 1,-1,-1,-1, 1, 1, 1, 0, 0, 0, 0,-1, 1,-1, 1, 0, 0,-1, 0, 0, 0, 0, 1},
  { 0, 1, 1,-1, 1,-1, 0, 0,-1, 0,-1, 1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0},
  { 1, 0,-1, 0, 0,-1, 0, 1,-1, 0,-1, 0,-1, 0,-1,-1, 0, 0, 0, 1, 0, 0,-1, 0,-1, 0, 0,-1},
  { 0, 1,-1, 1, 0, 1, 0, 0, 1,-1,-1, 0, 0, 0,-1, 0, 1, 1, 0, 1, 0,-1, 0, 1,-1, 0, 0,-1}
 });

ell = flatten entries(sub(L, S) * transpose vars S);
J0 = minimalPresentation(ideal q + ideal ell);
R0 = ring J0;
G0 = groebnerBasis(J0, Strategy => "F4"); 
forceGB G0;
J0 = ideal G0;                                
hs = reduceHilbert hilbertSeries J0;
print(hs);

c = {-1,1,0,-1,-1,1,0,1,0,0,1,1,1,1,-1,-1};
J = minimalPresentation(ideal q + ideal apply(16, i -> ell_i - c_i));
R = ring J;
G = groebnerBasis(J, Strategy => "F4");
Gb = forceGB G;
J = ideal G;   
                                 
A = R/J;
B = basis A;                                    
N = numColumns B;
Bl = lift(B, R);
expo = apply(flatten entries Bl, m -> first exponents m);
pos = hashTable apply(N, j -> expo#j => j);

multMatrix = f -> (
    fB := (f * Bl) % Gb;
    M := mutableMatrix(kk, N, N, Dense => true);
    for j from 0 to N-1 do scan(listForm fB_(0,j), (exp, cf) -> M_(pos#exp, j) = cf);
    M);

Mcube = new MutableList;
for k from 0 to numgens(R)-1 do (
    M := multMatrix R_k;
    Mcube#k = M * M * M;                       
    collectGarbage();)

Frob = mutableMatrix(kk, N, N, Dense => true);
cols = new MutableHashTable;

for j in sort(toList(0..N-1), j -> sum expo#j) do (
    e := expo#j;
    cols#j = if sum e == 0 then (
            u := mutableMatrix(kk, N, 1, Dense => true); u_(j,0) = 1_kk; u)
         else (
            k := position(e, t -> t > 0);
            Mcube#k * cols#(pos#(e - apply(numgens(R), i -> if i == k then 1 else 0))));
    for i from 0 to N-1 do Frob_(i,j) = cols#j_(i,0));

r = rank Frob;
if(r == N) then print("Frobenius is an isomorphism");

t1 = cpuTime();
print "Total time taken (in seconds)";
print(t1-t0); -- Expected 350 seconds