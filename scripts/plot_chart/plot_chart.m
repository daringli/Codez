function plot_chart(filename,type="p",titlename=filename,outfilename=strcat(filename,".eps"))

  magics=[8,20,28,50,82,126];
 
  h=figure(1);
  clf
  title(titlename) 
  hold on
 

  [Z, N, E, J, decay, p] = textread(filename, '%f %f %f %f %s %f','Delimiter','\t'); 
  l=length(Z);

  %make plot pretty
  set(gca(),
    'linewidth',2,
    'tickdir','out',
    'ticklength',[0.005,0.005],
    'xtick',[N(1):2:N(l)],
    'xlim',[N(1)-0.5,N(l)+0.5],
    'ytick',[Z(1):2:Z(l)],
    'ylim',[Z(1)-0.5,Z(l)+0.5]
  )   
  axis("equal")


  %find limits of N given Z and Z given N
  minZ=min(Z);
  maxZ=max(Z);
  minN=min(N);
  maxN=max(N);
  limZ=zeros(maxZ,2);
  limN=zeros(maxN,2);
  for z=minN:maxZ
    limN(z,1)=min(N(find(Z==z)));
    limN(z,2)=max(N(find(Z==z)));
  endfor
  for n=minN:maxN
    limZ(n,1)=min(Z(find(N==n)));
    limZ(n,2)=max(Z(find(N==n)));
  endfor

  %truncate list of magic numbers based on this
  magicsZ=magics(find((magics>=minZ) & (magics<=maxZ)));
  magicsN=magics(find((magics>=minN) & (magics<=maxN)));


  %loop over rows in indata
  if((strcmp("p", type) || strcmp("P", type))==1)
    for i=1:l
      plot_elem(Z(i),N(i),decay{i},p(i));
    endfor
  elseif((strcmp("e", type) || strcmp("E", type))==1)
    for i=1:l
      plot_elem(Z(i),N(i),decay{i},E(i));
    endfor
  endif

  %plot the magic lines
  plot_line(magicsZ,"Z",limN);
  plot_line(magicsN,"N",limZ);

  plot_stability();
  Zisn=10:90;
  Nisz=10:90;
  plot(Nisz,Zisn);
  %print to eps with color, a given size and a given font size.
    sizestring=strcat("-S",num2str((maxN-minN)*60),",",num2str((maxZ-minZ)*60));
  print(h,'-deps','-color',outfilename,sizestring,'-F:6')

endfunction
