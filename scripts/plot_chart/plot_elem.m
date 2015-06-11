function plot_elem(Z,N,decay,data)
  x=[N-0.5,N+0.5,N+0.5,N-0.5];
  y=[Z-0.5,Z-0.5,Z+0.5,Z+0.5];
  if(isnan(data))
    color='w'
  else
    if(strcmp("gamma", decay)==1)
      color='g';
    elseif(strcmp("p", decay)==1)
      color='r';
    elseif(strcmp("n", decay)==1)
      color='b';
    elseif(strcmp("alpha", decay)==1)
      color='c';
    else
      color='w';
    endif
  endif
    
  fill(x,y,color)
  text(N,Z,sprintf("%1.3f", data),'horizontalalignment', 'center','verticalalignment', 'middle')
  %text(N,Z,sprintf("%2.2f", E),'horizontalalignment', 'center','verticalalignment', 'bottom')
endfunction
