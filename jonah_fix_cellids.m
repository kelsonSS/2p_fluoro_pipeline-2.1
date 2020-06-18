function jonah_fix_cellids()

% this fixes an error in which the xy columns of the cell definition file
% were swapped
  [fname,pathname] =  uigetfile('//vault3/data/kelson/anlayzed','select CellDef file')
  destPath = fullfile(pathname,fname);
  
  cdef = open(destPath);
  
  xtemp = cdef.ptsIdx(:,2);
 cdef.ptsIdx(:,2) = cdef.ptsIdx(:,3);
 cdef.ptsIdx(:,3) = xtemp;
  
  
  
  
  save(destPath, '-struct','cdef')
