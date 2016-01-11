function [ xk, ene, dos ] = read_data( fn )
  % This function reads the output of greenplot
  fid=fopen(fn);
  str=fgetl(fid);
  tmp=sscanf(str, ' %d');
  nk=tmp(1);
  nen=tmp(2);
  
  for ik=1:10:nk
    str=fgetl(fid);
    if (ik+9<=nk)
      xk(ik:ik+9)=sscanf(str, '%f');
    else
      xk(ik:nk)=sscanf(str, '%f');
    end
  end
  
  for ie=1:10:nen
    str=fgetl(fid);
    if (ie+9<=nen)
      ene(ie:ie+9)=sscanf(str, '%f');
    else
      ene(ie:nen)=sscanf(str, '%f');
    end
  end
  
  if (length(xk)~=nk) 
    fprintf('Incorrect nk!\n');
  end
  
  if (length(ene)~=nen)
    fprintf('Incorrect nen!\n');
  end
  
  dos=zeros(nen, nk);
  
  for ik=1:nk
    for ie=1:10:nen
      str=fgetl(fid);
      if (ie+9<=nen)
        dos(ie:ie+9, ik)=sscanf(str, '%f');
      else
        dos(ie:nen, ik)=sscanf(str, '%f');
      end
    end
  end
  
  fclose(fid);
end
