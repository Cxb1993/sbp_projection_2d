% SBPWriteVTK2D write finite difference data to a VTK file
%
% inputs:
%   filename: filename
%   N1:       first grid dimension
%   N2:       second grid dimension
%   x1:       first coordinate array: size (N1+1)*(N2+1)
%   x2:       second coordinate array: size (N1+1)*(N2+1)
%   vanme:    cells array containing variable names
%   varargin: rest of arguments are the data fields
function SBPWriteVTK2D(filename, N1, N2, x1, x2, vnames, varargin)
  nfields = nargin - 6;

  Ntotal = (N1+1)*(N2+1);
  assert(Ntotal == length(x1(:)));

  x3 = zeros(size(x1));

  fid = fopen(filename, 'w');
  fprintf(fid, '# vtk DataFile Version 2');
  fprintf(fid, '\nelmers sbp');
  fprintf(fid, '\nASCII');
  fprintf(fid, '\nDATASET STRUCTURED_GRID');
  fprintf(fid, '\nDIMENSIONS %d %d 1', N1+1, N2+1);
  fprintf(fid, '\nPOINTS %d double', (N1+1)*(N2+1));
  fprintf(fid, '\n%25.16e %25.16e %25.16e', [x1(:) x2(:) x3(:)]');
  fprintf(fid, '\nPOINT_DATA %d', Ntotal);
  for n=1:nfields
    fprintf(fid, '\nSCALARS %s double 1', vnames{n});
    fprintf(fid, '\nLOOKUP_TABLE default');
    fprintf(fid, '\n%25.16e', varargin{n});
  end
  fclose(fid);

end
