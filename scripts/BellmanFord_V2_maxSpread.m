function [p,d] = BellmanFord_V2(s,t,w,mx,source,dest)
mx = mx+1;
n=max(t);
distance = Inf(mx,n);
predecessor = zeros(mx,n);
distance(:,source) = 0;


for i = 2:mx
    for j = 1:length(s)
        %         distance(s(j)) + w(j) < distance(t(j))
        if max(distance(i-1,s(j)), w(j)) <= distance(i,t(j))
            distance(i,t(j)) = max(distance(i-1,s(j)) , w(j));
            predecessor(i,t(j)) = s(j);
            %             pause
        end
    end
end

% for j = 1:length(s)
%     %         distance(s(j)) + w(j) < distance(t(j))
%     if distance(s(j)) + w(j) < distance(t(j))
%         error('Graph contains a negative-weight cycle');
%     end
% end
% distance
% predecessor

d = distance(end,dest);
ind = predecessor(end,dest);
p = [dest ind];
i = mx-1;
while(ind ~= source)
    ind = predecessor(i,ind);
    i = i - 1;
    p = [p ind];
end
p = flip(p);

distance;
predecessor;

end