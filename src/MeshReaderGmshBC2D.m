function [Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D(FileName)

% function [Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D(FileName)
% Purpose  : Read in basic grid information to build grid
% NOTE     : gmsh *.gmsh format is assumed

Globals2D;

vnum = [1 2; 2 3; 3 1];

Fid = fopen(FileName, 'rt');
fline = fgetl(Fid);

while ~strcmp(fline,'$PhysicalNames')
 fline = fgetl(Fid);
end

fline = fgetl(Fid);

npn = sscanf(fline, '%d');

PhyToBCType = zeros(npn, 1);

for i = 1:npn
  fline = fgetl(Fid);
  bck = sscanf(fline, '%d %*d %*s');
  bct = sscanf(fline, '%*d %d %*s');
  bcn = sscanf(fline, '%*s %*s %s');
  bcflag = 0;

  if strfind(bcn, 'In'),        bcflag = In;        end;
  if strfind(bcn, 'Out'),       bcflag = Out;       end;
  if strfind(bcn, 'Wall'),      bcflag = Wall;      end;
  if strfind(bcn, 'Far'),       bcflag = Far;       end;
  if strfind(bcn, 'Cyl'),       bcflag = Cyl;       end;
  if strfind(bcn, 'Dirichlet'), bcflag = Dirichlet; end;
  if strfind(bcn, 'Neuman'),    bcflag = Neuman;    end;
  if strfind(bcn, 'Slip'),      bcflag = Slip;      end;

  PhyToBCType(i) = bcflag;
end

while ~strcmp(fline,'$Nodes')
 fline = fgetl(Fid);
end

fline = fgetl(Fid);

Nv = sscanf(fline, '%d');

VX = (1:Nv);
VY = (1:Nv);
for i = 1:Nv
  fline = fgetl(Fid);
  tmpx = sscanf(fline, '%lf');
  VX(tmpx(1)) = tmpx(2);
  VY(tmpx(1)) = tmpx(3);
end

while (~strcmp(fline,'$Elements'))
 fline = fgetl(Fid);
end
fline = fgetl(Fid);

ka = sscanf(fline, '%d');

EToV = zeros(ka, 3);
BCType = zeros(ka,3);

BCToV = zeros(ka, 2);
BCToPhy = zeros(ka, 1);

bc = 1;
k = 1;
for i = 1:ka
  fline = fgetl(Fid);
  tmpcon = sscanf(fline, '%d');
  elm_number = tmpcon(1);
  elm_type = tmpcon(2);
  number_of_tags = tmpcon(3);
  tags = tmpcon(4:3+number_of_tags);
  node_number_list = tmpcon(4+number_of_tags:length(tmpcon));
  phy = tags(1);

  if (elm_type == 1)
    % bc
    ntags = tmpcon(3);
    BCToV(bc,:) = node_number_list;
    BCToPhy(bc) = phy;
    bc = bc + 1;
  elseif (elm_type == 2)
    % triangle
    EToV(k,:) = node_number_list;
    k = k + 1;
  else
    error('Unkown element');
  end
end

K = k - 1;
EToV = EToV(1:K,:);
BCType = BCType(1:K,:);

Nbc = bc - 1;
BCToV = BCToV(1:Nbc,:);
BCToPhy = BCToPhy(1:Nbc,:);

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);

% Find BCtypes
for i = 1:Nbc
  va = BCToV(i, 1);
  vb = BCToV(i, 2);
  p = BCToPhy(i);

  [rows, cols] = find(EToV == va);

  for j = 1:size(rows)
    f = 0;
    k = rows(j);
    ca = cols(j);
    cb = find(EToV(k,:) == vb);

    if(cb)
      for ftry = 1:3
        if(isempty(setdiff(vnum(ftry,:), [ca cb])))
          f = ftry;
        end
      end

      BCType(k, f) = PhyToBCType(p);
    end
  end
end

%% i = 1;
%% for k = 1:K
%%   for f = 1:3
%%     if(BCType(k, f))
%%       va = EToV(k, vnum(f,:))(1);
%%       vb = EToV(k, vnum(f,:))(2);
%%
%%       fprintf(1, '%d (%g, %g)--(%g, %g) == %d\n', ...
%%         i, VX(va), VY(va), VX(vb), VY(vb), BCType(k, f));
%%       i = i + 1;
%%     end
%%   end
%% end

st = fclose(Fid);
return
