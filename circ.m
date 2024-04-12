NX = 100;
sten = stencil(NX,NX,5050,20);
% radius = 50;
% % midpoint = 20100
% x_mid = 100-0.5;
% y_mid = 100-0.5;
% circle_inds = [];
% for row = 1:200
%     for col = 1:200
%         x = col-0.5;
%         y = row-0.5;
%         if (x-x_mid)^2+(y-y_mid)^2<=radius^2
%             circle_inds = [circle_inds, NX*(row-1)+col];
%         end        
%     end
% end
% circle = [] 
% for i = circle_inds
%     circle = [circle, find(sten == i)];
% end
%    
 
mid = 5050;
horizontal_inds = mid+[-7:7];
hori = [];
vertical_inds = mid+[-7:7]*NX;
verti = [];
for i = horizontal_inds
    hori = [hori, find(sten == i)];
end
for i = vertical_inds
    verti = [verti, find(sten == i)];
end