%Mesh Deformation. 
clear
close all

obj = readObj('man.obj');
FV = obj.f.v;
V = obj.v;
trimesh(FV(:,1:3), V(:,1), V(:,2)); 
w=1000;

[x,y] = ginput(4);

v = zeros(3,1);
for i =1:3
    min_dist= 10000; 
    t=0;
    for j =1:length(V) 
        d = sqrt((x(i)-V(j,1))^2 + (y(i) - V(j,2))^2);
        if (d < min_dist)
            min_dist = d;
            t = j;
        end
    end    
    x(i) = V(t,1);
    y(i) = V(t,2);
    v(i)=t;
end

hold on
plot(x,y,'o');
axis([-1.5 1.5 -1 2])

b1 = zeros(6*length(FV) + 6,1);
b1(6*length(FV)+1) = w*x(4);
b1(6*length(FV)+2) = w*y(4);
b1(6*length(FV)+3) = w*x(2);
b1(6*length(FV)+4) = w*y(2);
b1(6*length(FV)+5) = w*x(3);
b1(6*length(FV)+6) = w*y(3);



A1 = zeros(6*length(FV)+6, 2*length(V));
    Edges = zeros(2,8);
    Edges(1,1)=-1;
    Edges(1,3)=1;
    Edges(2,2)=-1;
    Edges(2,4)=1;
for i=1:length(FV)
    
    v1i=FV(i,1);
    v1j = FV(i,2);
    v1l = FV(i,3);   
    v1r=0;
    e1x = V(FV(i,2),1) - V(FV(i,1),1);
    e1y =  V(FV(i,2),2)-V(FV(i,1),2);
    E1 = [e1x, e1y; e1y -e1x];
    
%     Edges(6*(i-1)+3, 2*(FV(i,2)-1)+1) = 1;
%     Edges(6*(i-1)+4, 2*FV(i,2)) = 1;
%     Edges(6*(i-1)+3, 2*(FV(i,3)-1)+1) = -1;
%     Edges(6*(i-1)+4, 2*FV(i,3)) = -1;
    
    v2i=FV(i,2);
    v2j = FV(i,3);
    v2l = FV(i,1);
    v2r=0;
    e2x = V(FV(i,3),1) - V(FV(i,2),1);
    e2y =  V(FV(i,3),2) - V(FV(i,2),2);
    E2 = [e2x, e2y; e2y -e2x];

%     
%     Edges(6*(i-1)+5, 2*(FV(i,3)-1)+1) = 1;
%     Edges(6*(i-1)+6, 2*FV(i,3)) = 1;
%     Edges(6*(i-1)+5, 2*(FV(i,1)-1)+1) = -1;
%     Edges(6*(i-1)+6, 2*FV(i,1)) = -1;
    
    v3i=FV(i,3);
    v3j = FV(i,1);
    v3l = FV(i,2);
    v3r = 0;
    e3x = V(FV(i,1),1) - V(FV(i,3),1);
    e3y =  V(FV(i,1),2) - V(FV(i,3),2);
    E3 = [e3x, e3y; e3y -e3x];
    

    for j=1:length(FV)
        
        if (((v1i == FV(j,1))||(v1i == FV(j,2))||(v1i == FV(j,3)))&& ((v1j == FV(j,1))||(v1j == FV(j,2))||(v1j == FV(j,3)))&&((v1l ~= FV(j,1))&&(v1l ~= FV(j,2))&&(v1l ~= FV(j,3))))
                     v1r =(v1i ~= FV(j,1))*(v1j ~= FV(j,1))*FV(j,1)+(v1i ~= FV(j,2))*(v1j ~= FV(j,2))*FV(j,2)+(v1i ~= FV(j,3))*(v1j ~= FV(j,3))*FV(j,3);
        end
         if (((v2i == FV(j,1))||(v2i == FV(j,2))||(v2i == FV(j,3)))&& ((v2j == FV(j,1))||(v2j == FV(j,2))||(v2j == FV(j,3)))&&((v2l ~= FV(j,1))&&(v2l ~= FV(j,2))&&(v2l ~= FV(j,3))))
                     v2r =(v2i ~= FV(j,1))*(v2j ~= FV(j,1))*FV(j,1)+(v2i ~= FV(j,2))*(v2j ~= FV(j,2))*FV(j,2)+(v2i ~= FV(j,3))*(v2j ~= FV(j,3))*FV(j,3);
         end
         if (((v3i == FV(j,1))||(v3i == FV(j,2))||(v3i == FV(j,3)))&& ((v3j == FV(j,1))||(v3j == FV(j,2))||(v3j == FV(j,3)))&&((v3l ~= FV(j,1))&&(v3l ~= FV(j,2))&&(v3l ~= FV(j,3))))
                     v3r =(v3i ~= FV(j,1))*(v3j ~= FV(j,1))*FV(j,1)+(v3i ~= FV(j,2))*(v3j ~= FV(j,2))*FV(j,2)+(v3i ~= FV(j,3))*(v3j ~= FV(j,3))*FV(j,3);
        end
    end
  
if (v1r~=0)

%     G1 = [V(v1i,1), V(v1i,2), 1,0 ; V(v1i,2), -V(v1i,1), 0,1; V(v1j,1), V(v1j,2), 1,0 ; V(v1j,2), -V(v1j,1), 0, 1; V(v1l,1), V(v1l,2), 1,0 ; V(v1l,2), -V(v1l,1), 0,1 ; V(v1r,1), V(v1r,2), 1,0 ; V(v1r,2), -V(v1r,1),0,1];
    G1 = [V(v1i,1), V(v1i,2); V(v1i,2), -V(v1i,1); V(v1j,1), V(v1j,2); V(v1j,2), -V(v1j,1); V(v1l,1), V(v1l,2); V(v1l,2), -V(v1l,1); V(v1r,1), V(v1r,2); V(v1r,2), -V(v1r,1)];

G = inv(G1'*G1)*G1'; 
H1 = Edges - E1*G(1:2,:);
    A1(6*(i-1)+1:6*(i-1)+2, 2*(v1r-1)+1) = H1(1:2,7);
    A1(6*(i-1)+1:6*(i-1)+2, 2*(v1r-1)+2) = H1(1:2,8);
else
%     G1 = [V(v1i,1), V(v1i,2), 1,0 ; V(v1i,2), -V(v1i,1), 0, 1; V(v1j,1), V(v1j,2), 1,0 ; V(v1j,2), -V(v1j,1), 0, 1; V(v1l,1), V(v1l,2), 1,0 ; V(v1l,2), -V(v1l,1), 0,1 ];
G1 = [V(v1i,1), V(v1i,2); V(v1i,2), -V(v1i,1); V(v1j,1), V(v1j,2); V(v1j,2), -V(v1j,1); V(v1l,1), V(v1l,2); V(v1l,2), -V(v1l,1) ];

G = inv(G1'*G1)*G1';
    H1 = Edges(:, 1:6) - E1*G(1:2,:);
        
end    
    
if (v2r~=0)
   % G2 = [V(v2i,1), V(v2i,2), 1,0 ; V(v2i,2), -V(v2i,1), 0,1; V(v2j,1), V(v2j,2), 1,0 ; V(v2j,2), -V(v2j,1), 0,1; V(v2l,1), V(v2l,2), 1,0 ; V(v2l,2), -V(v2l,1), 0,1; V(v2r,1), V(v2r,2),  1,0 ; V(v2r,2), -V(v2r,1), 0,1];
        G2 = [V(v2i,1), V(v2i,2) ; V(v2i,2), -V(v2i,1); V(v2j,1), V(v2j,2); V(v2j,2), -V(v2j,1); V(v2l,1), V(v2l,2) ; V(v2l,2), -V(v2l,1); V(v2r,1), V(v2r,2); V(v2r,2), -V(v2r,1)];

    G = inv(G2'*G2)*G2';
    H2 = Edges - E2*G(1:2,:); 
     A1(6*(i-1)+3:6*(i-1)+4, 2*(v2r-1)+1) = H2(1:2,7);
    A1(6*(i-1)+3:6*(i-1)+4, 2*(v2r-1)+2) = H2(1:2,8);
else
%            G2 = [V(v2i,1), V(v2i,2), 1,0 ; V(v2i,2), -V(v2i,1), 0,1; V(v2j,1), V(v2j,2), 1,0 ; V(v2j,2), -V(v2j,1), 0,1; V(v2l,1), V(v2l,2), 1,0 ; V(v2l,2), -V(v2l,1), 0,1];
           G2 = [V(v2i,1), V(v2i,2); V(v2i,2), -V(v2i,1); V(v2j,1), V(v2j,2) ; V(v2j,2), -V(v2j,1); V(v2l,1), V(v2l,2); V(v2l,2), -V(v2l,1)];

G = inv(G2'*G2)*G2';
    H2 = Edges(:,1:6) - E2*G(1:2,:);  
end
if(v3r~=0)
%     G3 = [V(v3i,1), V(v3i,2),  1,0 ; V(v3i,2), -V(v3i,1), 0,1; V(v3j,1), V(v3j,2),  1,0 ; V(v3j,2), -V(v3j,1), 0,1; V(v3l,1), V(v3l,2),  1,0 ; V(v3l,2), -V(v3l,1), 0,1; V(v3r,1), V(v3r,2), 1,0 ; V(v3r,2), -V(v3r,1), 0,1];
       G3 = [V(v3i,1), V(v3i,2); V(v3i,2), -V(v3i,1); V(v3j,1), V(v3j,2); V(v3j,2), -V(v3j,1); V(v3l,1), V(v3l,2); V(v3l,2), -V(v3l,1); V(v3r,1), V(v3r,2); V(v3r,2), -V(v3r,1)];

G = inv(G3'*G3)*G3';
   H3 = Edges - E3*G(1:2,:);
   A1(6*(i-1)+5:6*(i-1)+6, 2*(v3r-1)+1) = H3(1:2,7);
    A1(6*(i-1)+5:6*(i-1)+6, 2*(v3r-1)+2) = H3(1:2,8);
else
%         G3 = [V(v3i,1), V(v3i,2),  1,0 ; V(v3i,2), -V(v3i,1), 0,1; V(v3j,1), V(v3j,2),  1,0 ; V(v3j,2), -V(v3j,1), 0,1; V(v3l,1), V(v3l,2),  1,0 ; V(v3l,2), -V(v3l,1), 0,1];
           G3 = [V(v3i,1), V(v3i,2); V(v3i,2), -V(v3i,1); V(v3j,1), V(v3j,2); V(v3j,2), -V(v3j,1); V(v3l,1), V(v3l,2); V(v3l,2), -V(v3l,1)];

G = inv(G3'*G3)*G3';
    H3 = Edges(:,1:6) - E3*G(1:2,:); 
end   
    A1(6*(i-1)+1:6*(i-1)+2, 2*(v1i-1)+1) = H1(1:2,1);
    A1(6*(i-1)+1:6*(i-1)+2, 2*(v1i-1)+2) = H1(1:2,2);
    A1(6*(i-1)+1:6*(i-1)+2, 2*(v1j-1)+1) = H1(1:2,3);
    A1(6*(i-1)+1:6*(i-1)+2, 2*(v1j-1)+2) = H1(1:2,4);
    A1(6*(i-1)+1:6*(i-1)+2, 2*(v1l-1)+1) = H1(1:2,5);
    A1(6*(i-1)+1:6*(i-1)+2, 2*(v1l-1)+2) = H1(1:2,6);
    
    A1(6*(i-1)+3:6*(i-1)+4, 2*(v2i-1)+1) = H2(1:2,1);
    A1(6*(i-1)+3:6*(i-1)+4, 2*(v2i-1)+2) = H2(1:2,2);
    A1(6*(i-1)+3:6*(i-1)+4, 2*(v2j-1)+1) = H2(1:2,3);
    A1(6*(i-1)+3:6*(i-1)+4, 2*(v2j-1)+2) = H2(1:2,4);
    A1(6*(i-1)+3:6*(i-1)+4, 2*(v2l-1)+1) = H2(1:2,5);
    A1(6*(i-1)+3:6*(i-1)+4, 2*(v2l-1)+2) = H2(1:2,6);
    
    A1(6*(i-1)+5:6*(i-1)+6, 2*(v3i-1)+1) = H3(1:2,1);
    A1(6*(i-1)+5:6*(i-1)+6, 2*(v3i-1)+2) = H3(1:2,2);
    A1(6*(i-1)+5:6*(i-1)+6, 2*(v3j-1)+1) = H3(1:2,3);
    A1(6*(i-1)+5:6*(i-1)+6, 2*(v3j-1)+2) = H3(1:2,4);
    A1(6*(i-1)+5:6*(i-1)+6, 2*(v3l-1)+1) = H3(1:2,5);
    A1(6*(i-1)+5:6*(i-1)+6, 2*(v3l-1)+2) = H3(1:2,6);
  
     
end

for i=1:3
A1(1992+2*(i-1)+1, 2*(v(i)-1)+1) = w;
A1(1992+2*(i-1)+2, 2*v(i)) = w;
end


New = inv(A1'*A1)*A1'*b1;

V=zeros(length(V),2);
for i =1:length(New)
    if (mod(i,2)==0)
        V(i/2, 2) = New(i);
    else
        V((i-1)/2+1,1) = New(i);
    end
    
end

figure
trimesh(FV(:,1:3), V(:,1), V(:,2)); 
hold on
plot(x,y,'o');
axis([-1.5 1.5 -1 2])