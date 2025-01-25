function [ BSistar] = BSi_star_expo( MUs, BT, K  )

MUs = MUs';

s = length(MUs);
BSistar = zeros(1,s);
if all(MUs == 0)
    BSistar(:) = BT/s;
    return;
end
set = find(MUs>0);

S=sum(1./K(set));

temp = BT./S./K(set) +(1./K(set)).*(log(MUs(set).*K(set))-log(prod((MUs(set).*K(set)).^(1./K(set)/S))));
BSistar(set)=temp;
while ~(all(BSistar>=0))
    set = find(BSistar>0);
    BSistar = zeros(1,s);
    S=sum(1./K(set));
    temp = BT./S./K(set) +(1./K(set)).*(log(MUs(set).*K(set))-log(prod((MUs(set).*K(set)).^(1./K(set)/S))));
    BSistar(set)=temp;
%     pause
end


end
