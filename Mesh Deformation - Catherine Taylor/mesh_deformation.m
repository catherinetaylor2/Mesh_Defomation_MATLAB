function V = mesh_deformation( filename )
   
    obj = readObj(filename); %load mesh info.
    FV = obj.f.v;
    V = obj.v;
    
    if (V(1:5,2) == 0)
        a = V(:,1);
        V(:,1) = -V(:,3);
        V(:,2) = -a;  
        V(:,3) = zeros(length(V),1);
    end 
           
    if (V(1:5,3) == 0)          
        
        figure
        subplot(1,2,1)
        trimesh(FV(:,1:3), V(:,1), V(:,2));
        w=1000;

        [x,y] = ginput(4); %select points.
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
        hold on
        plot(x,y,'o');

        %--------------------------------------------------------------------------
        %Algorithm 1:

        b1 = zeros(6*length(FV) + 6,1);
        b1(6*length(FV)+1:6*length(FV)+6) = w*[x(4),y(4), x(2),y(2),x(3),y(3)];
        A1 = zeros(6*length(FV)+6, 2*length(V));
        A2 = zeros(3*length(FV)+3, length(V));
        Edges = [-1,0,1,0,0,0,0,0;0,-1,0,1,0,0,0,0];
        for i=1:length(FV)
            A2(3*(i-1)+1, FV(i,1)) = -1; 
            A2(3*(i-1)+1, FV(i,2)) = 1; 
            A2(3*(i-1)+2, FV(i,2)) = -1; 
            A2(3*(i-1)+2, FV(i,3)) = 1; 
            A2(3*(i-1)+3, FV(i,3)) = -1; 
            A2(3*(i-1)+3, FV(i,1)) = 1; 
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
            A1(6*length(FV)+2*(i-1)+1, 2*(v(i)-1)+1) = w;
            A1(6*length(FV)+2*(i-1)+2, 2*v(i)) = w;
            A2(3*length(FV)+i, v(i)) = w;
        end
        
        V_new = (A1'*A1)\A1'*b1; %find least squares solution.
        V1=zeros(length(V),2);
        for i =1:length(V_new)
            if (mod(i,2)==0)
                V1(i/2, 2) = V_new(i);
            else
                V1((i-1)/2+1,1) = V_new(i);
            end
        end
        
        subplot(1,2,2)
        trimesh(FV(:,1:3), V1(:,1), V1(:,2));
        hold on
        plot(x,y,'o');
        
%--------------------------------------------------------------------------
%Algorithm 2:
        
        b2x = zeros(3*length(FV) + 3,1);
        b2y = zeros(3*length(FV) + 3,1);
        for i=1:length(FV)
            for k=1:3
                vi=FV(i,k); %Build up edge neighbours.
                vj = FV(i,mod(k,3)+1);
                vl = FV(i,(k==1)*3 + (k-1));
                vr=0;
                e = [V(vj,1) - V(vi,1); V(vj,2) - V(vi,2)];
                for j=1:length(FV) %Find right edge neighbour.
                    if (((vi == FV(j,1))||(vi == FV(j,2))||(vi == FV(j,3)))&& ((vj == FV(j,1))||(vj == FV(j,2))||(vj == FV(j,3)))&&((vl ~= FV(j,1))&&(vl ~= FV(j,2))&&(vl ~= FV(j,3))))
                        vr =(vi ~= FV(j,1))*(vj ~= FV(j,1))*FV(j,1)+(vi ~= FV(j,2))*(vj ~= FV(j,2))*FV(j,2)+(vi ~= FV(j,3))*(vj ~= FV(j,3))*FV(j,3);
                    end
                end
                if (vr~=0)
                    G = [V(vi,1), V(vi,2), 1,0 ; V(vi,2), -V(vi,1), 0,1; V(vj,1), V(vj,2), 1,0 ; V(vj,2), -V(vj,1), 0, 1; V(vl,1), V(vl,2), 1,0 ; V(vl,2), -V(vl,1), 0,1 ; V(vr,1), V(vr,2), 1,0 ; V(vr,2), -V(vr,1),0,1];
                    G = (G'*G)\G';
                    t = G(1:2,:)*[V1(vi,1), V1(vi,2), V1(vj,1), V1(vj,2), V1(vl,1), V1(vl,2), V1(vr,1),V1(vr,2)]';
                else
                    G = [V(vi,1), V(vi,2), 1,0 ; V(vi,2), -V(vi,1), 0, 1; V(vj,1), V(vj,2), 1,0 ; V(vj,2), -V(vj,1), 0, 1; V(vl,1), V(vl,2), 1,0 ; V(vl,2), -V(vl,1), 0,1 ];
                    G = (G'*G)\G';
                    t = G(1:2,:)*[V1(vi,1), V1(vi,2), V1(vj,1), V1(vj,2), V1(vl,1), V1(vl,2)]';
                end
                T = 1/(sqrt(t(1)^2+t(2)^2))*[t(1), t(2); -t(2), t(1)];
                b = T*e;
                b2x(3*(i-1)+k) = b(1);
                b2y(3*(i-1)+k) = b(2);
            end
        end
        b2x(3*length(FV)+1) = w*x(4);
        b2y(3*length(FV)+1) = w*y(4);
        for i=2:3
            b2x(3*length(FV)+i) = w*x(i);
            b2y(3*length(FV)+i) = w*y(i);
        end
        
        V2(:,1) = (A2'*A2)\A2'*b2x;
        V2(:,2) = (A2'*A2)\A2'*b2y;
        
        figure
        trimesh(FV(:,1:3), V2(:,1), V2(:,2));
        hold on
        plot(x,y,'o');
        
    else %3D case
        
        figure;
        trimesh(FV(:,1:3), V(:,1), V(:,2), V(:,3)); 
        w=1;
        [x,y,z] = ginput(4);
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
            z(i) = V(t,3);
            v(i)=t;
        end

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

    T = zeros(3*length(V)+9, 3*length(V));
    for i =1:length(V)
        Ai = zeros(3*length(Vertex_neighbours{i}+3), 7);
        Ai(1, :) = [V(i,1), 0, V(i,3), -V(i,2),1,0,0];
        Ai(2, :) = [V(i,2), -V(i,3), 0 ,V(i,1),0,1,0];
        Ai(3, :) = [V(i,3), V(i,2), -V(i,1),0,0,0,1];
        for j =1: length(Vertex_neighbours{i})
            Ai(3*(j-1)+4, :) = [V(Vertex_neighbours{i}(j),1), 0, V(Vertex_neighbours{i}(j),3), -V(Vertex_neighbours{i}(j),2),1,0,0];
            Ai(3*(j-1)+5, :) = [V(Vertex_neighbours{i}(j),2), -V(Vertex_neighbours{i}(j),3), 0 ,V(Vertex_neighbours{i}(j),1),0,1,0];
            Ai(3*(j-1)+6, :) = [V(Vertex_neighbours{i}(j),3), V(Vertex_neighbours{i}(j),2), -V(Vertex_neighbours{i}(j),1),0,0,0,1];
        end

        Ti = (Ai'*Ai)\Ai';
        Di  = [delta(i,1), 0, delta(i,3), -delta(i,2), 1,0,0; delta(i,2), -delta(i,3), 0 , delta(i,1), 0,1,0; delta(i,3), delta(i,2), - delta(i,1),0,0,0,1];

        TiDi = Di*Ti;
        Li =  zeros(3, 3*d(i)+3);
        Li(1,1) =1;
        Li(2,2) =1;
        Li(3,3)=1;
        for k=2:d(i)+1
            Li(1,3*(k-1)+1) =  -1/d(i);
            Li(2, 3*(k-1)+2) = -1/d(i);
            Li(3, 3*(k-1)+3) = -1/d(i);
        end

        TiDi = Li-TiDi;
        T(3*(i-1)+1, 3*(i-1)+ 1) = TiDi(1,1); 
        T(3*(i-1)+1, 3*(i-1)+2) = TiDi(1,2);
        T(3*(i-1)+1, 3*(i-1)+3) = TiDi(1,3);
        T(3*(i-1)+2, 3*(i-1)+1) = TiDi(2,1); 
        T(3*(i-1)+2, 3*(i-1)+2) = TiDi(2,2);
        T(3*(i-1)+2, 3*(i-1)+3) = TiDi(2,3);
        T(3*(i-1)+3, 3*(i-1)+1) = TiDi(3,1); 
        T(3*(i-1)+3, 3*(i-1)+2) = TiDi(3,2);
        T(3*(i-1)+3, 3*(i-1)+3) = TiDi(3,3);
        for j=1:d(i)
            T(3*(i-1)+1, 3*(Vertex_neighbours{i}(j)-1)+1) = TiDi(1,3*(j-1)+4); 
            T(3*(i-1)+1, 3*(Vertex_neighbours{i}(j)-1)+2) = TiDi(1,3*(j-1)+5);
            T(3*(i-1)+1, 3*(Vertex_neighbours{i}(j)-1)+3) = TiDi(1,3*(j-1)+6);
            T(3*(i-1)+2, 3*(Vertex_neighbours{i}(j)-1)+1) = TiDi(2,3*(j-1)+4); 
            T(3*(i-1)+2, 3*(Vertex_neighbours{i}(j)-1)+2) = TiDi(2,3*(j-1)+5);
            T(3*(i-1)+2, 3*(Vertex_neighbours{i}(j)-1)+3) = TiDi(2,3*(j-1)+6);
            T(3*(i-1)+3, 3*(Vertex_neighbours{i}(j)-1)+1) = TiDi(3,3*(j-1)+4); 
            T(3*(i-1)+3, 3*(Vertex_neighbours{i}(j)-1)+2) = TiDi(3,3*(j-1)+5);
            T(3*(i-1)+3, 3*(Vertex_neighbours{i}(j)-1)+3) = TiDi(3,3*(j-1)+6);
        end

        b1(3*(i-1)+1) =0;
        b1(3*(i-1)+2) =0;
        b1(3*(i-1)+3) = 0;

        TiDi =[];
        Ti=[];
        Li=[];
        Ai=[];
    end
    for i=1:3
        T(3*(i-1)+3*length(V)+1, 3*(v(i)-1)+1) = w;
        T(3*(i-1)+2+3*length(V), 3*(v(i)-1)+2) = w;
        T(3*(i-1)+3+3*length(V), 3*(v(i)-1)+3) = w; 
    end
    V2 = (T'*T)\T'*b1;
   
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
    end
    
end