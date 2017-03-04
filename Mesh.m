%Mesh Deformation. 
clear
close all

obj = readObj('man.obj');
FV = obj.f.v;
V = obj.v;
trimesh(FV(:,1:3), V(:,1), V(:,2)); 
w=1000;

% for i=1:100

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

% hold on
% plot(x,y,'o');
% axis([-1.5 1.5 -1 2])


% bx = zeros(3*length(FV) + 3,1) ; 
% by = zeros(3*length(FV) + 3,1) ; 
% for i =1:length(FV)
%    bx(3*(i-1)+1) = V(FV(i,1),1) - V(FV(i,2),1); 
%    bx(3*(i-1)+2) = V(FV(i,2),1) - V(FV(i,3),1); 
%    bx(3*(i-1)+3) = V(FV(i,3),1) - V(FV(i,1),1); 
%    
%     by(3*(i-1)+1) = V(FV(i,1),2) - V(FV(i,2),2); 
%     by(3*(i-1)+2) = V(FV(i,2),2) - V(FV(i,3),2); 
%     by(3*(i-1)+3) = V(FV(i,3),2) - V(FV(i,1),2); 
% end
% bx(997) = w*x(4);
% bx(998) = w*x(2);
% bx(999) = w*x(3);
% 
% by(997) = w*y(4);
% by(998) = w*y(2);
% by(999) = w*y(3);


b1 = zeros(6*length(FV) + 6);
for i=1:length(FV)
    b(6*(i-1)+1) = V(FV(i,1),1) - V(FV(i,2),1); 
    b(6*(i-1)+2) = V(FV(i,1),2) - V(FV(i,2),2); 
    b1(6*(i-1)+3) = V(FV(i,2),1) - V(FV(i,3),1); 
    b1(6*(i-1)+4) = V(FV(i,2),2) - V(FV(i,3),2); 
    b1(6*(i-1)+5) =  V(FV(i,3),1) - V(FV(i,1),1); 
    b1(6*(i-1)+6) =  V(FV(i,3),2) - V(FV(i,1),2); 
end

% 
% 
% 
% A = zeros(3*length(FV)+3, length(V));
% 
% for i=1:length(FV)
%     A(3*(i-1)+1, FV(i,1)) = 1; 
%     A(3*(i-1)+1, FV(i,2)) = -1; 
%     A(3*(i-1)+2, FV(i,2)) = 1; 
%     A(3*(i-1)+2, FV(i,3)) = -1; 
%     A(3*(i-1)+3, FV(i,3)) = 1; 
%     A(3*(i-1)+3, FV(i,1)) = -1; 
% end
% for i=1:3
% A(996+i, v(i)) = w;
% end
% 
% V(:,1) = (A'*A)\A'*bx;
% V(:,2) = (A'*A)\A'*by;
% 
% hold off
% %figure
% trimesh(FV(:,1:3), V(:,1), V(:,2)); 
% hold on
% plot(x,y,'o');
% axis([-1.5 1.5 -1 2])
% 
% hold off
% % end

