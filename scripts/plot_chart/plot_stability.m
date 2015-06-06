function plot_stability
  filename="isotopes.dat";
  %read 2,3 and 5th column
  [Z, N, decay] = textread(filename, "%*s %*s %s %s %*s %*s %*s %*s %s %*s %*s %*s %*s %*s","delimiter",'\t',"headerlines",1); 


  Zs=str2double(Z(find(strcmp("STABLE",decay))));
  Ns=str2double(N(find(strcmp("STABLE",decay))));

  plot(Ns,Zs,'ow',"markersize",12,'markerfacecolor', 'white');

endfunction
