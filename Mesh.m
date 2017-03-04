%Mesh Deformation. 
clear
close all

obj = readObj('man.obj'); %load mesh info.
FV = obj.f.v;
V = obj.v;

figure
subplot(1,2,1)
trimesh(FV(:,1:3), V(:,1), V(:,2)); 
axis([-1.5 1.5 -2 2])
w=1000;

[x,y] = ginput(4); %select points.
hold on
plot(x,y,'o');

v = zeros(3,1); %Find closest vertex.
for i =1:3
    min_dist= 10000;
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

b1 = zeros(6*length(FV) + 6,1);
b1(6*length(FV)+1:6*length(FV)+6) = w*[x(4),y(4), x(2),y(2),x(3),y(3)];
A1 = zeros(6*length(FV)+6, 2*length(V));
Edges = [-1,0,1,0,0,0,0,0;0,-1,0,1,0,0,0,0];
for i=1:length(FV)
    for k =1:3
        vi=FV(i,k); %Build up edge neighbours.
        vj = FV(i,mod(k,3)+1);
        vl = FV(i,(k==1)*3 + (k-1));
        vr=0;
        ex = V(vj,1) - V(vi,1);
        ey =  V(vj,2)-V(vi,2);
        E = [ex, ey; ey -ex];
        for j=1:length(FV) %Find right edge neighbour.
            if (((vi == FV(j,1))||(vi == FV(j,2))||(vi == FV(j,3)))&& ((vj == FV(j,1))||(vj == FV(j,2))||(vj == FV(j,3)))&&((vl ~= FV(j,1))&&(vl ~= FV(j,2))&&(vl ~= FV(j,3))))
                vr =(vi ~= FV(j,1))*(vj ~= FV(j,1))*FV(j,1)+(vi ~= FV(j,2))*(vj ~= FV(j,2))*FV(j,2)+(vi ~= FV(j,3))*(vj ~= FV(j,3))*FV(j,3);
            end
        end
        if (vr~=0)
            G = [V(vi,1), V(vi,2), 1,0 ; V(vi,2), -V(vi,1), 0,1; V(vj,1), V(vj,2), 1,0 ; V(vj,2), -V(vj,1), 0, 1; V(vl,1), V(vl,2), 1,0 ; V(vl,2), -V(vl,1), 0,1 ; V(vr,1), V(vr,2), 1,0 ; V(vr,2), -V(vr,1),0,1];
            G = (G'*G)\G';
            H = Edges - E*G(1:2,:);
            A1(6*(i-1)+2*k-1:6*(i-1)+2*k, 2*(vr-1)+1) = H(1:2,7);
            A1(6*(i-1)+2*k-1:6*(i-1)+2*k, 2*(vr-1)+2) = H(1:2,8);
        else
            G = [V(vi,1), V(vi,2), 1,0 ; V(vi,2), -V(vi,1), 0, 1; V(vj,1), V(vj,2), 1,0 ; V(vj,2), -V(vj,1), 0, 1; V(vl,1), V(vl,2), 1,0 ; V(vl,2), -V(vl,1), 0,1 ];
            G = (G'*G)\G';
            H = Edges(:, 1:6) - E*G(1:2,:);
        end        
        A1(6*(i-1)+2*k-1:6*(i-1)+2*k, 2*(vi-1)+1) = H(1:2,1);
        A1(6*(i-1)+2*k-1:6*(i-1)+2*k, 2*(vi-1)+2) = H(1:2,2);
        A1(6*(i-1)+2*k-1:6*(i-1)+2*k, 2*(vj-1)+1) = H(1:2,3);
        A1(6*(i-1)+2*k-1:6*(i-1)+2*k, 2*(vj-1)+2) = H(1:2,4);
        A1(6*(i-1)+2*k-1:6*(i-1)+2*k, 2*(vl-1)+1) = H(1:2,5);
        A1(6*(i-1)+2*k-1:6*(i-1)+2*k, 2*(vl-1)+2) = H(1:2,6);
    end
end

for i=1:3
    A1(1992+2*(i-1)+1, 2*(v(i)-1)+1) = w;
    A1(1992+2*(i-1)+2, 2*v(i)) = w;
end

V_new = (A1'*A1)\A1'*b1; %find least squares solution.
V=zeros(length(V),2);
for i =1:length(V_new)
    if (mod(i,2)==0)
        V(i/2, 2) = V_new(i);
    else
        V((i-1)/2+1,1) = V_new(i);
    end    
end

subplot(1,2,2)
trimesh(FV(:,1:3), V(:,1), V(:,2)); 
hold on
plot(x,y,'o');
axis([-1.5 1.5 -2 2])


%Not containing 0,1:
%G1 = [V(v1i,1), V(v1i,2); V(v1i,2), -V(v1i,1); V(v1j,1), V(v1j,2); V(v1j,2), -V(v1j,1); V(v1l,1), V(v1l,2); V(v1l,2), -V(v1l,1); V(v1r,1), V(v1r,2); V(v1r,2), -V(v1r,1)];
%G1 = [V(v1i,1), V(v1i,2); V(v1i,2), -V(v1i,1); V(v1j,1), V(v1j,2); V(v1j,2), -V(v1j,1); V(v1l,1), V(v1l,2); V(v1l,2), -V(v1l,1) ];