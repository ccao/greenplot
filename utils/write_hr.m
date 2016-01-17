function write_hr( ham, rvec, wt )
%read_hr Writes real space hamiltonian to _hr.dat file
%   Detailed explanation goes here
  norb=size(ham, 1);
  nrpt=size(ham, 3);
  fid=fopen('sc_hr.dat', 'w');
  fprintf(fid, ' Output from expandwannhr\n%12d\n%12d\n', norb, nrpt);
  
  for ii=1:nrpt
    fprintf(fid, '%5d', wt(ii));
    if (mod(ii,15)==0||ii==nrpt)
      fprintf(fid, '\n');
    end
  end
  
  for ii=1:nrpt
    for io=1:norb
      for jo=1:norb
        fprintf(fid, '%5d%5d%5d%5d%5d%12.6f%12.6f\n', rvec(ii,:), io, jo, real(ham(io, jo, ii)), imag(ham(io, jo, ii)));
      end
    end
  end
    
  fclose(fid);
end

