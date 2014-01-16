syms x y 
u1 = (1-x)*(1-y);
u2 = x*(1-y);
u3 = y*(1-x);
u4 = x*y;

gu1 = [diff(u1,x);diff(u1,y)];
gu2 = [diff(u2,x);diff(u2,y)];
gu3 = [diff(u3,x);diff(u3,y)];
gu4 = [diff(u4,x);diff(u4,y)];

b = [1 0]';
bu1 = b'*gu1;
bu2 = b'*gu2;
bu3 = b'*gu3;
bu4 = b'*gu4;

C = zeros(4);
C(1,1) = int(int(bu1*u1,x,0,1),y,0,1);
C(1,2) = int(int(bu2*u1,x,0,1),y,0,1);
C(1,3) = int(int(bu3*u1,x,0,1),y,0,1);
C(1,4) = int(int(bu4*u1,x,0,1),y,0,1);

C(2,1) = int(int(bu1*u2,x,0,1),y,0,1);
C(2,2) = int(int(bu2*u2,x,0,1),y,0,1);
C(2,3) = int(int(bu3*u2,x,0,1),y,0,1);
C(2,4) = int(int(bu4*u2,x,0,1),y,0,1);

C(3,1) = int(int(bu1*u3,x,0,1),y,0,1);
C(3,2) = int(int(bu2*u3,x,0,1),y,0,1);
C(3,3) = int(int(bu3*u3,x,0,1),y,0,1);
C(3,4) = int(int(bu4*u3,x,0,1),y,0,1);

C(4,1) = int(int(bu1*u4,x,0,1),y,0,1);
C(4,2) = int(int(bu2*u4,x,0,1),y,0,1);
C(4,3) = int(int(bu3*u4,x,0,1),y,0,1);
C(4,4) = int(int(bu4*u4,x,0,1),y,0,1);