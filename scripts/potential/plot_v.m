#! /usr/bin/octave -qf
#Usage: reads data to plot from stdin. 
#Flag -J or -j makes it plot files with fixed J and different nuclei
#while the default output makes plots for fixed nuclei and different J.
#Flag -n <name> specifies a name to attach in front of the numberous outfiles

addpath ("/n/home/bstefan/m/deexcite/git/11-05-2015/scripts/")

#plot style parameters
ysize=500;
ylims=[-15,8];
xsize=500;
fontsize=10;
fontname="/usr/share/fonts/truetype/droid/DroidSans.ttf";
%legendlocation="northeast";
legendlocation="southeast";
no_legend=false;

xl='r/fm';
yl='V/MeV';

sizestring=strcat("-S",num2str(xsize),",",num2str(ysize));
fontstring=strcat("-F",fontname,":",num2str(fontsize));

%style={'-','--',':','-.','-r','--r',':r','-.r','-g','--g',':g','-.g'};
color=prism(20);

#defaults
prename="";


arg_list = argv ();
argc=length(arg_list);

if(argc>0)
  for i=1:argc
    flag=arg_list{i}(1:end);
    if(strcmp(flag,"-N") || strcmp(flag,"-n"))
      i=i+1;
      prename=arg_list{i}(1:end);
    endif
  endfor
endif

fid=stdin();
[metas,datas]=stdin_tables_parse(fid);

%parse meta-data
%store mothers in mothers{j} and corresponding l and evaps in ls{j} evaps{j}
mothers=zeros(1,2);
mother_index=zeros(1);
ls=zeros(1);
evaps=zeros(1);

j=1;
m=1;

for i=1:length(metas)
  off=0;
  metas{i}
  if(length(metas{i})==3)
    if(i~=1)
    legend(legendstring,'location',legendlocation);
      if(no_legend==true)
        legend hide
      endif

      outfilename=strcat(prename,"potZ",num2str(mothers(j,1)),"N",num2str(mothers(j,2)),".eps")
      print(h,'-deps','-color',outfilename,sizestring,fontstring)
      endif

    mother_string=metas{i}{1}(1:end);
    Zindex = strfind(mother_string,'Z=')+2;
    tabindex = strfind(mother_string(Zindex:end),"\t")+Zindex-1;
    mothers(j,1)=str2num(mother_string(Zindex:tabindex));
    Nindex = strfind(mother_string,'N=')+2;
    mothers(j,2)=str2num(mother_string(Nindex:end));
    off=1;

    h=figure(1);
    clf
    hold on;
    xlabel(xl,'interpreter','tex');
    ylabel(yl,'interpreter','tex');   
    title(strcat("Potential to tunnel through under various \nevaporation processes from (Z,N)=(",num2str(mothers(j,1)),",",num2str(mothers(j,2)),")."));  
    ylim(ylims);
  endif

  evap_string=metas{i}{1+off}(1:end);
  Zindex = strfind(evap_string,'Zevap=')+6;
  tabindex = strfind(evap_string(Zindex:end),"\t")+Zindex-1;
  evaps(1)=str2num(evap_string(Zindex:tabindex));
  Nindex = strfind(evap_string,'Nevap=')+6;
  evaps(2)=str2num(evap_string(Nindex:end));
  
  l_string=metas{i}{2+off}(1:end);
  lindex = strfind(l_string,'=')+1;
  num_l=length(lindex);
  for ii=1:num_l
    if(ii<num_l)
      ls(ii)=str2double(l_string(lindex(ii):(lindex(ii+1)-3)));
    else
      ls(ii)=str2double(l_string(lindex(ii):end));
    endif
    i
    p=plot(datas{i}(:,1),datas{i}(:,ii+1));
    set(p,'Color',color(m,:)),; % sets the color of the plot
    legendstring{m}=strcat("(Zev,Nev)=(",num2str(evaps(1)),",", num2str(evaps(2)),") l=",num2str(ls(ii)));
    m=m+1;
  endfor
endfor

    legend(legendstring,'location',legendlocation);
      if(no_legend==true)
        legend hide
      endif
      
      outfilename=strcat(prename,"potZ",num2str(mothers(j,1)),"N",num2str(mothers(j,2)),".eps")
      print(h,'-deps','-color',outfilename,sizestring,fontstring)