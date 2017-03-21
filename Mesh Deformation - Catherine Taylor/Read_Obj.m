function [V, FV] = Read_Obj( filename )

f=fopen(filename);

V=[];
FV=[];
while 1
    t=fgetl(f);
    if ~ischar(t),   break,   end 
tt=sscanf(t,'%s',1);

if (tt=='v')
 V = [V; sscanf(t(2:end),'%f')'];
end
if (tt=='f')
 FV = [FV; sscanf(t(2:end),'%f')'];
end

end
end

