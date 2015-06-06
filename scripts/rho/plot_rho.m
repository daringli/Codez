#! /usr/bin/octave -qf
#Usage: reads data to plot from stdin. 
#Flag -J or -j makes it plot files with fixed J and different nuclei
#while the default output makes plots for fixed nuclei and different J.
#Flag -n <name> specifies a name to attach in front of the numberous outfiles

#plot style parameters
ysize=500;
xsize=500;
fontsize=10;
fontname="/usr/share/fonts/truetype/droid/DroidSans.ttf";
%legendlocation="northwest";
legendlocation="southeast";
no_legend=true

xl='E/MeV';
yl='\rho/(MeV^{-1})';

sizestring=strcat("-S",num2str(xsize),",",num2str(ysize));
fontstring=strcat("-F",fontname,":",num2str(fontsize));

%style={'-','--',':','-.','-r','--r',':r','-.r','-g','--g',':g','-.g'};
color=copper(20);

#defaults
useJ=0;
prename="";


arg_list = argv ();
argc=length(arg_list);

for i=1:argc
  flag=arg_list{i}(1:end);
  if(strcmp(flag,"-J") || strcmp(flag,"-j"))
    useJ=1;
  elseif(strcmp(flag,"-N") || strcmp(flag,"-n"))
    i=i+1;
    prename=arg_list{i}(1:end);
  endif
endfor

fid=stdin();

n=1; %index of columns
j=1; %index of rows
E=cell(0);
rho=cell(0);
while (line=fgetl(fid))!=-1
  %process comment lines for metadata
  %assume 2 comment lines in a row, as output from rho gives
  if (line(1)=='#')
    %first comment line says which nuclei it is
    Zindex = strfind(line,'Z=')+2;
    tabindex = strfind(line(Zindex:end),"\t")+Zindex-1;
    Z=str2num(line(Zindex:tabindex));

    Nindex = strfind(line,'N=')+2;
    N=str2num(line(Nindex:end));
    %second comment line explains the data
    line=fgetl(fid);
    Jindices = strfind(line,'=')+1;
    num_J=length(Jindices);
    J=zeros(1,num_J);
    format="%f";
    %read all J
    m=n; %so that we know which index ZNJ starts at for this ZN
    for i=1:num_J
      format=strcat(format,"%f");
      if (i<num_J)
      J(i)=str2double(line(Jindices(i):(Jindices(i+1)-3)));
      else
      J(i)=str2double(line(Jindices(i):end));
      endif
      ZNJ(n,:)=[Z,N,J(i)];
      n=n+1;
    endfor
    line=fgetl(fid);
  endif
  values=sscanf(line,format);
  for i=m:(n-1)
     rho{i}(j)=values(i-m+2);
     E{i}(j)=values(1);
  endfor
  j=j+1;
endwhile

%plot this data in different files according to useJ
if (useJ==1)
  %loop over all unique J to produce different plots
  for j=unique(ZNJ(:,3)')
   h=figure(1);
   clf
   hold on
   xlabel(xl,'interpreter','tex');
   ylabel(yl,'interpreter','tex');
   title(strcat("Level densities for J=",num2str(j)));
   %loop over all nuclei that have been run with this J
   m=1;
   for i=find(ZNJ(:,3)'==j)
      Z=ZNJ(i,1);
      N=ZNJ(i,2);
      p=semilogy(E{i},rho{i});
      set(p,'Color',color(m,:)),; % sets the color of the plot
      legendstring{m}=strcat("(Z,N)=(",num2str(Z),",",num2str(N),")");
      m=m+1;
   endfor
   legend(legendstring,'location',legendlocation);
   if(no_legend==true)
      legend hide
   endif
   %export to eps
   outfilename=strcat(prename,"J",num2str(j),".eps")
   sizestring=strcat("-S",num2str(xsize),",",num2str(ysize));
   print(h,'-deps','-color',outfilename,sizestring,fontstring)
  endfor
else
%loop over all nuclei to produce different plots
  for z=unique(ZNJ(:,1)')
   for n=unique(ZNJ(find(ZNJ(:,1)==z),2)')
     %loop over all J for this nuclei
      h=figure(1);
      clf
      hold on;
      xlabel(xl,'interpreter','tex');
      ylabel(yl,'interpreter','tex');   
      title(strcat("Level densities for (Z,N)=(",num2str(z),",",num2str(n),")"));
      m=1;
      for i=find((ZNJ(:,1)==z).*(ZNJ(:,2)==n))'
        j=ZNJ(i,3);
        p=semilogy(E{i},rho{i}); 
        set(p,'Color',color(m,:)),; % sets the color of the plot
        legendstring{m}=strcat("j=",num2str(j));
        m=m+1;
      endfor
      legend(legendstring,'location',legendlocation);
      if(no_legend==true)
        legend hide
      endif
      copied_legend = findobj(gcf(),"type","axes","Tag","legend");
      set(copied_legend, "FontSize", 5);
      %export to eps
      outfilename=strcat(prename,"Z",num2str(z),"N",num2str(n),".eps")
      print(h,'-deps','-color',outfilename,sizestring,fontstring)
  endfor
 endfor
endif