function [ ham, rvec, wt ] = read_hr( fn )
%read_hr Reads real space hamiltonian from _hr.dat file
%   Detailed explanation goes here
  fid=fopen(fn, 'r');
  fgetl(fid);
  norb=fscanf(fid, ' %d\n', 1);
  nrpt=fscanf(fid, ' %d\n', 1);
  
  fprintf('# norb: %d, nrpt: %d\n', norb, nrpt);
  
  rvec=zeros(nrpt, 3);
  ham=zeros(norb, norb, nrpt);
  
  wt=fscanf(fid, ' %d', nrpt);
  fgetl(fid);
    
  for ii=1:nrpt
    for io=1:norb
      for jo=1:norb
        tmp=fscanf(fid, ' %d %d %d %*d %*d %f %f ', 5);
        rvec(ii, :)=tmp(1:3);
        ham(io, jo, ii)=tmp(4)+1i*tmp(5);
      end
    end
  end
  
  fclose(fid);
end

