fid=fopen('supercell.def', 'r');
str=fgetl(fid);
fn=str;
ulat=zeros(3);
sc=zeros(3);
for ii=1:3
  str=fgetl(fid);
  ulat(ii,:)=sscanf(str, ' %f %f %f');
end
for ii=1:3
  str=fgetl(fid);
  sc(ii,:)=sscanf(str, ' %d %d %d');
end
fgetl(fid);
str=fgetl(fid);
rrmax=sscanf(str, ' %d %d %d');
str=fgetl(fid);
nsite=sscanf(str, ' %d');
xx=zeros(nsite,3);
for ii=1:nsite
  str=fgetl(fid);
  x(ii,:)=sscanf(str, ' %f %f %f ');
end
xx=x/ulat;
fclose(fid);

for ii=1:nsite
  fprintf('%5d WF center: ( %12.9f, %12.9f, %12.9f ) \n', ii, xx(ii, :));
end

nuc=det(sc);
if (nuc<0) 
  fprintf(' !!! FATAL : supercell definition is not positive definite!');
end

slat=sc*ulat;

xxx=fillsc(xx, sc);

fprintf(' Super cell info\n 1.0\n');
for ii=1:3
  fprintf(' %16.9f%16.9f%16.9f\n', slat(ii, :));
end
fprintf('%d\n WF\n', nuc*nsite);
for ii=1:nuc*nsite
  fprintf('%16.9f %16.9f %16.9f \n', xxx(ii, :));
end
fprintf('********\n');

[rrvec, wwt]=find_ws(slat, rrmax);

nnrpt=size(rrvec, 1);

[ham, rvec, wt]=read_hr(fn);

revlist=revindex(rvec);

nrpt=size(ham, 3);
norb=size(ham, 1);

if (nsite~=norb)
  fprintf(' !!! FATAL ERROR: different # of orbitals!\n');
end

rmax=max(rvec);

rmax_in_sc=ceil(2.0.*max(rvec/sc));

fprintf('Recommended minimum size for supercell : %5d %5d %5d\n', rmax_in_sc);

% start to transform
nnorb=norb*nuc;
hham=zeros(nnorb,nnorb,nnrpt);

for iirpt=1:nnrpt
  for iio=1:nnorb
    for jjo=1:nnorb
      io=mod(iio-1,norb)+1;
      jo=mod(jjo-1,norb)+1;
      [ixx, iRu]=sc2uc(xxx(iio, :), [0 0 0], sc);
      [jxx, jRu]=sc2uc(xxx(jjo, :), rrvec(iirpt, :), sc);
      rr=iRu-jRu;
      if ( abs(rr(1))>rmax(1) || abs(rr(2))>rmax(2) || abs(rr(3))>rmax(3) )
        irpt=0;
      else
        ir=rr+rmax+1;
        irpt=revlist(ir(1), ir(2), ir(3));
      end
      if (irpt==0)
        hham(iio, jjo, iirpt)=0.0;
      else
        hham(iio, jjo, iirpt)=ham(io, jo, irpt);
      end  % if  irpt
    end  % jjo
  end  % iio
end  % iirpt

write_hr(hham, rrvec, wwt); 
