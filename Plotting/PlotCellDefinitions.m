function PlotCellDefinitions(cdef)
% takes a cdef struct and makes a plot of the cell definitions 


 
  for pp = 1:length(cdef.ptsIdx)
      hold on
    xpt = cdef.ptsIdx(pp,2);
    ypt = cdef.ptsIdx(pp,3);
    plot(  xpt + cdef.smRoiBoundaries{pp}(:,3) .* (cos(cdef.smRoiBoundaries{pp}(:,1)))  ,  ypt + cdef.smRoiBoundaries{pp}(:,3) .* (sin(cdef.smRoiBoundaries{pp}(:,1))),'g' ) %plot outer edges of rings
    plot(  xpt + cdef.smRoiBoundaries{pp}(:,2) .* (cos(cdef.smRoiBoundaries{pp}(:,1)))  ,  ypt + cdef.smRoiBoundaries{pp}(:,2) .* (sin(cdef.smRoiBoundaries{pp}(:,1))),'g' ) %plot inner edges of rings
    plot(  xpt + cdef.smNpBoundaries{pp}(:,3) .* (cos(cdef.smNpBoundaries{pp}(:,1)))  ,  ypt + cdef.smNpBoundaries{pp}(:,3) .* (sin(cdef.smNpBoundaries{pp}(:,1))),'r','LineWidth',.5 ) %plot outer edges of rings
    plot(  xpt + cdef.smNpBoundaries{pp}(:,2) .* (cos(cdef.smNpBoundaries{pp}(:,1)))  ,  ypt + cdef.smNpBoundaries{pp}(:,2) .* (sin(cdef.smNpBoundaries{pp}(:,1))),'r','LineWidth',.1 ) %plot inner edges of rings
  end 
end