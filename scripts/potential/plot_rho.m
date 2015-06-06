#! /usr/bin/octave -qf
#Usage: reads data to plot from stdin. 
#Flag -J or -j makes it plot files with fixed J and different nuclei
#while the default output makes plots for fixed nuclei and different J.
#Flag -n <name> specifies a name to attach in front of the numberous outfiles

addpath ("/n/home/bstefan/m/deexcite/git/11-05-2015/scripts/")

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
[metas,datas]=stdin_tables_parse(fid)