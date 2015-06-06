function [metas, datas]=stdin_tables_parse(fid)
  %%parses fid with C-like I/O, assuming it is composed of several tables separated by header lines beginning with #. It can thus parse STDIN
  metas=cell(0);
  datas=cell(0);
  
  data=zeros(0);
  

  was_comment=false;
  metas_index=1; 
  datas_index=0;
				%read from fid.
  while(line=fgetl(fid))~=-1
    %if the line is a comment line
    if (line(1)=='#')
	%add line to meta
	if(!was_comment)
	   %if the previous line was not a comment, 
	   %we end the data cell and start a new metas cell
	   %unless we have no data
	   if(datas_index~=0)
	      datas{datas_index}=data;
           endif
           data=zeros(0);
	   datas_index=datas_index+1;

	   data_row=1;
	   meta_index=0;
        endif
	meta_index=meta_index+1;
	meta{meta_index}=line;
	%mark that this line was comment
	was_comment=true;
	%go to next line
	continue;
     endif
     %we only get here if we did not find a comment
     if(was_comment==true)
       %if a group of comment was recently concluded
       metas{metas_index}=meta;
       meta=cell(0);
       metas_index=metas_index+1;
       was_comment=false;
     endif


     %parse data by looking for tabs
     tabindices = strfind(line,"\t");
     tabindices =[1, tabindices,length(line)];
	for data_col=1:(length(tabindices)-2)
	data(data_row,data_col)=str2num(line(tabindices(data_col):tabindices(data_col+1)));
      endfor
     data_row=data_row+1;
 endwhile
%at end of file, add final data.
datas{datas_index}=data;
endfunction