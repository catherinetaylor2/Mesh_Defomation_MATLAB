%Mesh Deformation in 3D. 
clear
close all

obj = readObj('cube.obj'); %load mesh info.
FV = obj.f.v;
V = obj.v;
V= V(1:8,:);
FV = FV(1:12,:);
h=figure;
trimesh(FV(:,1:3), V(:,1), V(:,2), V(:,3)); 
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
 w=1;
% % [x,y,z] = ginput(4);
% 
% pts = zeros(4,3);
% datacursormode on
% dcm_obj = datacursormode(h);
% 
% for i=1:4
%     waitforbuttonpress;
%     f = getCursorInfo(dcm_obj);
%     pts(i,:) = f.Position;
% end
% 
% x = pts(:,1);
% y=pts(:,2);
% z=pts(:,3);
%v = zeros(3,1);

v = [1, 4,8];
x = [V(v(1),1), V(v(2),1),V(v(2),1), -1];
y = [V(v(1),2), V(v(2),2),V(v(2),2), -1];
z = [V(v(1),3), V(v(2),3),V(v(2),3), 1];

% for i =1:3
%     min_dist= 10000; 
%     t=0;
%     for j =1:length(V) 
%         d = sqrt((x(i)-V(j,1))^2 + (y(i) - V(j,2))^2)+(z(i) - V(j,2)^2);
%         if (d < min_dist)
%             min_dist = d;
%             t = j;
%         end
%     end    
%     x(i) = V(t,1);
%     y(i) = V(t,2);
%     z(i) = V(t,3);
%     v(i)=t;
% end

hold on
plot3(x(2:3),y(2:3),z(2:3),'o');
plot3(x(1),y(1),z(1),'ro');
plot3(x(4),y(4),z(4),'go');

d = zeros(length(V),1);
D=zeros(length(V));
A = zeros(length(V));
Vertex_neighbours = cell(1, length(V),1);
for i=1:length(V)
    t=0;
    tt=0;
    for j=1:length(FV)
       if ((i == FV(j,1))||(i==FV(j,2))||(i==FV(j,3)))
           t=t+1;
           d(i)=d(i)+1;
           g(3*(t-1)+1:3*(t-1)+3) = FV(j,1:3); 
       end
    end
    D(i,i) = d(i);
    [b,~,~] = unique(g,'first');
    Vertex_neighbours{i} = zeros(length(b)-1,1);
    for k = 1:length(b)
        if (b(k)~=i)
            tt=tt+1;
            Vertex_neighbours{i}(tt) = b(k);
        end
    end
g=[];
end

for i=1:length(V)
    for j=1:length(Vertex_neighbours{i})
        A(i, Vertex_neighbours{i}(j))=1;
    end
end 

L = eye(length(V)) - D\A;
delta = L*V(:, 1:3);

A1 = zeros(3*length(V)+9, 3*length(V));
b1 = zeros(3*length(V)+9,1);
for i=1:length(V)
    A1(3*(i-1)+1, 3*(i-1)+1) = 1;
    A1(3*(i-1)+2, 3*(i-1)+2) = 1;
    A1(3*(i-1)+3, 3*(i-1)+3) = 1;
    b1(3*(i-1)+1) = delta(i,1);
    b1(3*(i-1)+2) = delta(i,2);
    b1(3*(i-1)+3) = delta(i,3);
   for j =1: d(i)
       A1(3*(i-1)+1, 3*(Vertex_neighbours{i}(j)-1)+1) = -1/d(i);
       A1(3*(i-1)+2, 3*(Vertex_neighbours{i}(j)-1)+2) = -1/d(i);
       A1(3*(i-1)+3, 3*(Vertex_neighbours{i}(j)-1)+3) = -1/d(i);
   end
end

   for i=1:3
    A1(3*(i-1)+3*length(V)+1, 3*(v(i)-1)+1) = w;
    A1(3*(i-1)+2+3*length(V), 3*(v(i)-1)+2) = w;
    A1(3*(i-1)+3+3*length(V), 3*(v(i)-1)+3) = w; 
   end
    b1(3*length(V)+1) = w*x(4);
    b1(3*length(V)+2) = w*y(4);
    b1(3*length(V)+3) = w*z(4);
   for i=2:3
    b1(3*length(V)+3*(i-1)+1) = w*x(i);
    b1(3*length(V)+3*(i-1)+2) = w*y(i);
    b1(3*length(V)+3*(i-1)+3) = w*z(i);
   end
   
   V1 = (A1'*A1)\A1'*b1;
   
  Vn=zeros(length(V),3);
for i =1:length(V1)
    if (mod(i,3)==1)
        Vn((i+2)/3, 1) = V1(i);
    elseif (mod(i,3)==2)
        Vn((i+1)/3,2) = V1(i);
    else
            Vn(i/3,3) = V1(i);
    end    
end

figure;
trimesh(FV(:,1:3), Vn(:,1), Vn(:,2), Vn(:,3)); 
hold on
plot3(x(2:3),y(2:3),z(2:3),'o');
plot3(x(1),y(1),z(1),'ro');
plot3(x(4),y(4),z(4),'go');
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])

T = zeros(3*length(V)+9, 3*length(V));
for i =1:length(V)
   
    Ai = zeros(3*length(Vertex_neighbours{i}), 7);
    bi = zeros(3*length(Vertex_neighbours{i}), 1);
    for j =1: length(Vertex_neighbours{i})
        Ai(3*(j-1)+1, :) = [V(Vertex_neighbours{i}(j),1), 0, V(Vertex_neighbours{i}(j),3), -V(Vertex_neighbours{i}(j),2),1,0,0];
        Ai(3*(j-1)+2, :) = [V(Vertex_neighbours{i}(j),2), -V(Vertex_neighbours{i}(j),3), 0 ,V(Vertex_neighbours{i}(j),1),0,1,0];
        Ai(3*(j-1)+3, :) = [V(Vertex_neighbours{i}(j),3), V(Vertex_neighbours{i}(j),2), -V(Vertex_neighbours{i}(j),1),0,0,0,1];
        bi(3*(j-1)+1) =     Vn(Vertex_neighbours{i}(j),1); 
        bi(3*(j-1)+2) =     Vn(Vertex_neighbours{i}(j),2);
        bi(3*(j-1)+3) =     Vn(Vertex_neighbours{i}(j),3);
    end
    Ti = (Ai'*Ai)\Ai';
    Di  = [delta(i,1), 0, delta(i,3), -delta(i,2), 1,0,0; delta(i,2), -delta(i,3), 0 , delta(i,1), 0,1,0; delta(i,3), delta(i,2), - delta(i,1),0,0,0,1];
    TiDi = Di*Ti;
    for j=1:d(i)
        T(3*(i-1)+1, 3*(Vertex_neighbours{i}(j)-1)+1) = TiDi(1,3*(j-1)+1); 
        T(3*(i-1)+1, 3*(Vertex_neighbours{i}(j)-1)+2) = TiDi(1,3*(j-1)+2);
        T(3*(i-1)+1, 3*(Vertex_neighbours{i}(j)-1)+3) = TiDi(1,3*(j-1)+3);
        T(3*(i-1)+2, 3*(Vertex_neighbours{i}(j)-1)+1) = TiDi(2,3*(j-1)+1); 
        T(3*(i-1)+2, 3*(Vertex_neighbours{i}(j)-1)+2) = TiDi(2,3*(j-1)+2);
        T(3*(i-1)+2, 3*(Vertex_neighbours{i}(j)-1)+3) = TiDi(2,3*(j-1)+3);
        T(3*(i-1)+3, 3*(Vertex_neighbours{i}(j)-1)+1) = TiDi(3,3*(j-1)+1); 
        T(3*(i-1)+3, 3*(Vertex_neighbours{i}(j)-1)+2) = TiDi(3,3*(j-1)+2);
        T(3*(i-1)+3, 3*(Vertex_neighbours{i}(j)-1)+3) = TiDi(3,3*(j-1)+3);
    end
    b1(3*(i-1)+1) =0;
    b1(3*(i-1)+2) = 0;
    b1(3*(i-1)+3) = 0;
    
    
end

   V2 = ((A1-T)'*(A1-T))\(A1-T)'*b1;
   
  V2n=zeros(length(V),3);
for i =1:length(V2)
    if (mod(i,3)==1)
        V2n((i+2)/3, 1) = V2(i);
    elseif (mod(i,3)==2)
        V2n((i+1)/3,2) = V2(i);
    else
            V2n(i/3,3) = V2(i);
    end    
end

figure;
trimesh(FV(:,1:3), V2n(:,1), V2n(:,2), V2n(:,3)); 
hold on
plot3(x(2:3),y(2:3),z(2:3),'o');
plot3(x(1),y(1),z(1),'ro');
plot3(x(4),y(4),z(4),'go');
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])