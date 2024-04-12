function [sammlung_dist_sorted,sammlung_ews_sorted]=rearrange(sammlung_dist,sammlung_ews)
%Sortiert die Verteilungen der Hauefigkeit nach aufsteigend
len = size(sammlung_ews,2);
ew_len = zeros(1,len);
for i = 1:len
   ew_len(i)=size(sammlung_ews{i},2);
end
[~,I] = sort(ew_len);
sammlung_ews_sorted = sammlung_ews;
for i = 1:len
    sammlung_ews_sorted{i}=sammlung_ews{I(i)};
end
sammlung_dist_sorted = sammlung_dist;
for i = 1:len
    sammlung_dist_sorted(:,2*i-1:2*i)=sammlung_dist(:,2*I(i)-1:2*I(i));
end
end