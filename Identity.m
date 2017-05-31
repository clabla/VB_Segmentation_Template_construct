%generate arrays of X Y Z voxel coordinates--------------------------------

function [Co,d1]=Identity(d,st,varargin)
d  = single(d);
st = single(st);
[e1,e2,e3] = meshgrid(st(2):st(2):d(2),st(1):st(1):d(1),st(3):st(3):d(3));
d1         = size(e1);
if ~isempty(varargin)
    Co      = cat(2,e2(:),e1(:),e3(:));
    Co(:,4) = 1;
else
    Co = cat(4,e2,e1,e3);
end
end