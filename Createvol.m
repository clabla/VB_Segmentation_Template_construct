%create output volumes-----------------------------------------------------

function Out = createvol(name,dir,data,varargin)

nfields = numel(varargin);
if mod(nfields,2)~=0
    error('invalid input number')
end
name = [name '.nii'];
Out  = nifti;
for i=1:0.5*nfields
    Out.(varargin{2*i-1}) = varargin{2*i};
end
Out.dat = file_array(fullfile(dir,name),size(data),'float32');
create(Out);
if numel(size(data))==3
    Out.dat(:,:,:) = data;
elseif numel(size(data))==4
    Out.dat(:,:,:,:) = data;
end
end