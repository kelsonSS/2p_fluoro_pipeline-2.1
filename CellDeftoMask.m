function Masks = CellDeftoMask(Cdef)  

if isstring(Cdef)
    Cdef = load(Cdef);
end 



    
for pp = length(Cdef.ptsIdx)
    
   xpt = Cdef.ptsIdx(:,2);
   ypt = Cdef.ptsIdx(:,3);
   ROIxvOut{pp} =  xpt + cdef.smRoiBoundaries{pp}(:,3) .* (cos(cdef.smRoiBoundaries{pp}(:,1))) ;
   ROIyvOut{pp} =  ypt + cdef.smRoiBoundaries{pp}(:,3) .* (sin(cdef.smRoiBoundaries{pp}(:,1))) ;
   ROIxvIn{pp} =  xpt + cdef.smRoiBoundaries{pp}(:,2) .* (cos(cdef.smRoiBoundaries{pp}(:,1))) ;
   ROIyvIn{pp} =  ypt + cdef.smRoiBoundaries{pp}(:,2) .* (sin(cdef.smRoiBoundaries{pp}(:,1))) ;
   roiBWout{pp} = poly2mask( ROIxvOut{pp} , ROIyvOut{pp} , cdef.dimX , cdef.dimX);
   roiBWin{pp} = poly2mask( ROIxvIn{pp} , ROIyvIn{pp} , cdef.dimX , cdef.dimX);
   roiBW2{pp} =  roiBWout{pp} -  roiBWin{pp};
   if sum(roiBW2{pp}(:) < 0) > 0 %accounts for inner diameter extending beyond outer diameter
       roiBW2 {pp} (roiBW2 {pp} < 0 ) = 0;
   end 
end  
   Masks.roiBW2 = roiBW2;
   Masks.roiBWin = roiBWin;
   Masks.roiBWout = roiBWout;