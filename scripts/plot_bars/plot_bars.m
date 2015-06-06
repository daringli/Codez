function plot_bars(filename,titlename=filename,outfilename=strcat(filename,".eps"))

fontsize=24;
fontname="/usr/share/fonts/truetype/droid/DroidSans.ttf";

fontstring=strcat("-F",fontname,":",num2str(fontsize));



  h=figure(1);
  clf
  title(titlename) 
  hold on

  [Z, N, E, J, part, prob] = textread(filename, '%f %f %f %f %s %f','Delimiter','\t'); 
  l=length(Z);

  nuclei=unique([Z,N],"rows");
  ln=length(nuclei);
  x=cell(ln,1);
  y=zeros(ln,7);     
  %set the x-label of the nuclei
  for i=1:ln
     %x{i} = ['Z=',num2str(nuclei(i,1)),' N=',num2str(nuclei(i,2))];
     x{i} = [num2str(nuclei(i,1))];

  endfor
  for i=1:l
     %find the index in nuclei corresponding to this row
     nuc_ind=find(nuclei(:,1)==Z(i) & nuclei(:,2)==N(i));
     if strcmp(part{i},"gamma")==1
        y(nuc_ind,1)=prob(i);
     elseif strcmp(part{i},"n")==1
        y(nuc_ind,2)=prob(i);
     elseif strcmp(part{i},"p")==1
        y(nuc_ind,3)=prob(i);
     elseif strcmp(part{i},"alpha")==1
        y(nuc_ind,4)=prob(i);
     elseif strcmp(part{i},"d")==1
        y(nuc_ind,5)=prob(i);
     elseif strcmp(part{i},"t")==1
        y(nuc_ind,6)=prob(i);
     elseif strcmp(part{i},"He3")==1
        y(nuc_ind,7)=prob(i);
     else
       %this is bad and should not happen
     endif
  endfor
 set(gca(),    'xtick',[1:ln],
    'XTickLabel',x,
    'xlim',[0,ln+1],
    'ylim',[0,1])

  b=bar(y,'stacked');
 set (b(1), 'facecolor', 'g');
 set (b(2), 'facecolor', 'b');
 set (b(3), 'facecolor', 'r');
 set (b(4), 'facecolor', 'c');
 set (b(5), 'facecolor', 'k');
 set (b(6), 'facecolor', 'y');
 set (b(7), 'facecolor', 'm');

%sizestring=strcat("-S",num2str((maxN-minN)*60),",",num2str((maxZ-minZ)*60));
%sizestring=strcat("-S",num2str((maxN-minN)*60),",",num2str((maxZ-minZ)*60));
xlabel("Z")
ylabel("p")
  print(h,'-deps','-color',outfilename,fontstring)


endfunction
