%a=[0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,1,1,0;0,0,1,1,0,0;1,0,1,1,1,1;1,1,1,1,1,1;1,1,1,1,1,1;1,1,1,1,1,1;1,1,1,1,1,1];
function [possible_etch, possible_grow] = pixelanalysis(a)

%%
%a=cat(3,zeros(3),zeros(3),[1,1,1;1,0,1;0,0,1],ones(3),ones(3));

sz=size(a);
above=zeros(sz(1:end-1));
below=ones(sz(1:end-1));
a2=cat(3,above,a,below);

% % X+
% c1=circshift(a2,-1,1);
% % X-
% c2=circshift(a2,1,1);
% % Y+
% c3=circshift(a2,-1,2);
% % Y-
% c4=circshift(a2,1,2);
% Z+
c5=circshift(a2,-1,3);
% Z-
c6=circshift(a2,1,3);

% How many neighbours
%c=c1+c2+c3+c4+c5+c6;

%res=c(:,:,2:end-1);

upnoblock = ~c6(:,:,2:end-1);
downblock = c5(:,:,2:end-1);

% If less than 6 neighbours, is occupied, and has no block above
%possible_etch=(res<6)&(a==1)&(upnoblock==1);

possible_etch=(a==1)&(upnoblock==1);

% If it has a neighbour, is unoccupied, and has a block below
%possible_grow=(res>0)&(a==0)&(downblock==1);

possible_grow=(a==0)&(downblock==1);
%%
end