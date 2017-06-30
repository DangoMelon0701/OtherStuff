function u=u(E,tab,tipo)
Ea=tab{1}';
%ut=tab{2}';
%uf=tab{3}';
switch tipo
    case 1
        ue=tab{2}';
    case 2
        ue=tab{3}';
end
[~,i]=min(abs(Ea-E));
dis=Ea(i)-E;
if dis<0
   u=ue(i)+(ue(i+1)-ue(i))*(E-Ea(i))/(Ea(i+1)-Ea(i));
elseif dis>0
    if i==1
        u=ue(i)+(ue(i+1)-ue(i))*(E-Ea(i))/(Ea(i+1)-Ea(i));
    elseif i>1
        u=ue(i-1)+(ue(i)-ue(i-1))*(E-Ea(i-1))/(Ea(i)-Ea(i-1));
    end
elseif dis==0
   u=ue(i);
end
end