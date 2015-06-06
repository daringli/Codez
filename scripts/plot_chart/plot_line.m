function plot_line(i,zn, lims)

  if((strcmp("z", zn) || strcmp("Z", zn))==1)
    x=[lims(i,1)'-0.5;lims(i,2)'+0.5;lims(i,2)'+0.5;lims(i,1)'-0.5];
    y=[i-0.5;i-0.5;i+0.5;i+0.5];
  elseif((strcmp("n", zn) || strcmp("N", zn))==1)
    y=[lims(i,1)'-0.5;lims(i,2)'+0.5;lims(i,2)'+0.5;lims(i,1)'-0.5];
    x=[i-0.5;i-0.5;i+0.5;i+0.5];
  else
    return;
  endif


    
  plot(x,y,'k',"linewidth",3)
endfunction
