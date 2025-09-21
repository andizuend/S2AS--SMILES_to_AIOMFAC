!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module including subroutines to produce custom 2-D x-y scatter plots, curves,      *
!*   bar graphs, and pie charts using the DISLIN graphics library.                      *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, Switzerland; 2007 - 2009; 2013                                 *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA; 2009 -2012    *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2007                                                            *
!*   -> latest changes: 2025-07-12                                                      *
!*                                                                                      *
!****************************************************************************************
module ModCreatePlot

use Mod_kinds, only : wp
use ModRoundScale, only : roundnicely, scalenicely

implicit none   
!...
integer,public :: NCplot 	        !NCplot = maximum number of distinct plot labels and colors considered; needs to be set during plotting setup
integer,private :: curvenr, curves_l, curves_r, txtl
!...
integer,parameter,public  :: ncprows = 400		        !maximum number of data rows in plot data
integer,parameter,private :: ncurves = int(1E5)         !the set maximum number of different curves in a plot
integer,parameter,private :: nshareas = 20            	!the maximum number of different shaded areas in a plot
integer,parameter,private :: messagetextlength = 200  	!the maximum number of characters per line in the text message
integer,parameter,private :: messagelinesmax = 40    	!the maximum number of lines in the message box
!...
integer,public :: messboxpos_x, messboxpos_y, areapmodeno, barpmodeno, curvepointsno
integer,dimension(6),private :: compbpoints, ncomp
integer,dimension(nshareas),private :: areabpoints
integer,dimension(ncurves),private :: cwaxis, ccolor, cltype, rownr, csymbtype, clstyle
integer,dimension(:,:,:),allocatable,private :: cbarbordcol, cbarshadpat, cbarpatterncol !allocated by call to SetCreateplotArrays
integer,dimension(nshareas),private :: areabordcol, areashadpat, areapatterncol
integer,dimension(:),allocatable,public :: colneutral
integer,dimension(12),public :: colelectrol
integer,dimension(254),private :: customzcolno
!...
real(wp),parameter,public :: rdefault = -8.88888E280_wp !an unusual, small (negative) value...																						  
real(wp),public :: max_y, min_y, lowerylimit, upperylimit, lowerylimit2, upperylimit2, lowerxlimit, upperxlimit											  
!private reals:
real(wp),parameter,private :: deps = epsilon(1.0_wp)
real(wp),private :: zcolbar0, zcolbar254, zcolbarstep													
real(wp),dimension(ncurves),private :: ryraymax, lyraymax, xmaxa, xmina, yminal, yminar, cth
real(wp),dimension(:,:),allocatable,private :: cxvals, cyvals
real(wp),dimension(:,:,:),allocatable,private :: cxyerror
real(wp),dimension(6,ncprows),private :: compbx
real(wp),dimension(nshareas,ncprows),private :: areapx
real(wp),dimension(:,:),allocatable,private :: areapyLow, areapyUp	
real(wp),dimension(:,:,:),allocatable,private :: compby
!...
character(len=3),public :: metaff
character(len=:),allocatable,private :: defaultInpPath, defaultOutPath, defaultColPath
character(len=75),private :: zaxislabel
character(len=150),private :: customcolfile										  
character(len=75),dimension(ncurves),private :: clegtxt
character(len=75),dimension(nshareas),private :: alegtxt					
character(len=75),dimension(6),private :: clylabel,crylabel !bar graph y-axis labels
character(len=85),dimension(6),public :: pietittxt !pie chart title text array
character(len=75),dimension(:),allocatable,public :: pielegtxt !pie chart data legend text array
character(len=85),dimension(ncurves,6),private :: cpietittxt !pie chart title text array
character(len=75),dimension(:,:),allocatable,private :: cpielegtxt !pie chart data legend text array
character(len=messagetextlength),dimension(messagelinesmax),private :: messagetxt
!...
logical,public :: bigDsymbs, biglabelfont, bigsymbs, customCol, equalfontwidth, logxaxis, mirrorAxis, &
    & slimplot, slimplotX, smallsymbs, specialxrange, specialyrange, specialyrange2, squareplot, ternaryLLEoutplot, &
    & ticksinside, bargraphplot, compbars, comptemp, Cstarplot, directory_outputGouri, fixminxax, legplot, logyaxis, &
    & messagebox, newplot, newplot2, newplot3, nogridlines, plotbarborder, plotoverlay, xorg3bar
logical,private :: discreteZcol, zlogscalebar, showcustomcolbar, genericPlotModule												
logical,dimension(12),public :: bargylimits !used to define upper yaxis limits for up to 3 different bar graphs on a plot page
logical,dimension(ncurves),public :: pielabshift
logical,dimension(ncurves),private :: clogyaxis,clegplot,cbigDsymbs
logical,dimension(ncprows),public :: plotbarlabels
logical,dimension(:,:),allocatable,public :: cplotbarlabels

!=================================================================================================================
contains  !the subroutines within this module...
!=================================================================================================================  
        
subroutine SetCreateplotArrays(NCp, genericPlotM, inpPath, outPath, colPath) !utility subroutine to allocate module arrays and set certain plot parameters

implicit none

integer,intent(in) :: NCp                                !the same value as NCplot
logical,intent(in) :: genericPlotM                       !set to .true. if ModCreatePlot is used outside of AIOMFAC program;
character(len=*),intent(in) :: inpPath, outPath, colPath    !the default input/output and colour palette file directory path
!local variables:
integer :: il
!..........................

NCplot = NCp
genericPlotModule = genericplotM

if (allocated(defaultInpPath)) deallocate(defaultInpPath)
il = len_trim(inpPath)
allocate(character(len=il) :: defaultInpPath)
defaultInpPath = trim(inpPath)

if (allocated(defaultOutPath)) deallocate(defaultOutPath)
il = len_trim(outPath)
allocate(character(len=il) :: defaultOutPath)
defaultOutPath = trim(outPath)  !sets the default (relative) output directory path used if no other options are set.  

if (allocated(defaultColPath)) deallocate(defaultColPath)
il = len_trim(colPath)
allocate(character(len=il) :: defaultColPath)
defaultColPath = trim(colPath)

if (.not. genericPlotModule) then !set a default custom color file;
    customcolfile = trim(defaultColPath)//'RGBColTableRealWinXP8BitAZ.dat' 
endif

!allocate arrays that depend on the value of NCp, the number of distinct plot components/colors/curves considered.
allocate( cbarbordcol(6,2*NCp,ncprows), cbarshadpat(6,2*NCp,ncprows), cbarpatterncol(6,2*NCp,ncprows), &
    & colneutral(NCp), compby(6,2*NCp,ncprows), pielegtxt(2*NCp), cpielegtxt(ncurves,2*NCp),  &
    & cplotbarlabels(2*NCp,ncprows), areapyLow(nshareas,ncprows), areapyUp(nshareas,ncprows) )

colneutral = 1
pielegtxt = ""

end subroutine SetCreateplotArrays
!================================================================================================================  

subroutine DeallocCreateplotArrays() !utility subroutine to deallocate module arrays

implicit none
!allocate arrays that depend on the value of NCp, the number of distinct plot components/colors/curves considered.
deallocate( cbarbordcol, cbarshadpat, cbarpatterncol, colneutral, compby, pielegtxt, cpielegtxt, cplotbarlabels, &
    & areapyLow, areapyUp )

end subroutine DeallocCreateplotArrays
!=================================================================================================================  


subroutine DefinePlotData(xvals, yvals, xyerror, waxis, color, ltype, symbtype, lstyle, legtxt, thick)
!meaning of the curve defining parameters:
!DefinePlotData(xvals, yvals, xyerror(the +- error bar size for x and y axis [dimensions 1:4 are x-err neg, x-err pos, y-err neg, y-err pos]), waxis(y-axis 1 or 2?), 
!               color, ltype(1=only line, 2=symb and line, > 2 means symbols and no line), symbtype(-1=no symbol), lstyle(1=solid, 2=dashm, 3=dash, 4=Mydot, 5=dot, 6=dash-dotted),legtxt,cthickn)
    
implicit none
!...
real(wp),dimension(:),intent(in) :: xvals, yvals
real(wp),dimension(1:size(xvals),1:4),intent(in) :: xyerror
integer,intent(in) :: waxis, color, ltype, symbtype, lstyle
real(wp),intent(in) :: thick
!...
real(wp) :: yy
integer :: i, N, npts, nsize, allocstat
real(wp),dimension(5) :: yrayspec
!...
character(len=75),intent(in) :: legtxt
!...
logical :: isdata
!.........................................................................................

if (newplot) then  !this start a new plot page, thus reset the counters and arrays:
    if (.not. allocated(cxvals)) then
        allocate( cxvals(ncurves,ncprows), cyvals(ncurves,ncprows), cxyerror(ncurves,ncprows,4), stat = allocstat)
    endif
    curvenr = 0
    curves_l = 0
    curves_r = 0
    cxvals = 0.0_wp
    cyvals = 0.0_wp
    cxyerror = 0.0_wp
    cwaxis = 0
    ccolor = 0
    cltype = 0
    rownr = 0
    cth = 1.0_wp !curvethickness
    csymbtype = 0
    clstyle = 0
    clegtxt = ""
    xmaxa = 0.0_wp
    xmina = 0.0_wp
    yminal = 0.0_wp
    yminar = 0.0_wp
    lyraymax = 0.0_wp
    ryraymax = 0.0_wp
    clogyaxis = .false.
    clegplot = .true.
    cpietittxt = ""
    cpielegtxt = ""
    cbigDsymbs = .false.
    call defplotcolorlists() !call to define color tables for plots
endif

nsize = min(size(xvals(:), DIM=1), size(yvals(:), DIM=1))

!check data range:
npts = 0
curvepointsno = npts
do i = 1,nsize-1
    if (xvals(i) > (-1.0E8_wp) .and. ((xvals(i) > deps .or. xvals(i) < (-deps)) .or. (xvals(i+1) > deps .or. xvals(i+1) < (-deps)) .or. &
        &  (yvals(i) > deps .or. yvals(i) < (-deps)) .or. (yvals(i+1) > deps .or. yvals(i+1) < (-deps)))) then
    npts = npts+1 !count points worth plotting...
    else
        exit !leave the loop and continue...
    endif
enddo
!check the last array entry:
i = nsize
if (xvals(i) > -1.0E8_wp .and. ((xvals(i) > deps .or. xvals(i) < (-deps)) .or. (yvals(i) > deps .or. yvals(i) < (-deps)))) then
    npts = npts+1 !count points worth plotting...
endif

isdata = .true.
if (npts < 1) then
    !    write(*,*) "Warning: no points of the xvals are different from 0.0_wp!"
    !    write(*,*) "this warning was generated in: subroutine DefinePlotData"
    isdata = .false.
endif

if (isdata) then
    !set the curve attributes in the array:
    curvenr = curvenr+1
    curvepointsno = npts
    cxvals(curvenr,1:npts) = xvals(1:npts)
    cyvals(curvenr,1:npts) = yvals(1:npts)
    do i = 1,npts
        cxyerror(curvenr,i,1:4) = xyerror(i,1:4)
    enddo
    cwaxis(curvenr) = waxis
    ccolor(curvenr) = color  !the color number (0 - 255) (default: in the RAIN colour spectrum of DISLIN)
    cltype(curvenr) = ltype
    rownr(curvenr) = npts
    csymbtype(curvenr) = symbtype
    cbigDsymbs(curvenr) = bigDsymbs
    clstyle(curvenr) = lstyle
    clegtxt(curvenr) = trim(legtxt)
    cpietittxt(curvenr,:) = pietittxt
    cpielegtxt(curvenr,:) = pielegtxt
    cth(curvenr) = thick
    xmaxa(curvenr) = maxval(xvals(1:npts))
    xmina(curvenr) = minval(xvals(1:npts))
    if (waxis == 1) then
        curves_l = curves_l +1
        yminal(curves_l) = minval(yvals(1:npts))
    else
        curves_r = curves_r +1
        yminar(curves_r) = minval(yvals(1:npts))
    endif
    clogyaxis(curvenr) = logyaxis
    clegplot(curvenr) = legplot
    if (legplot .and. len_trim(legtxt) < 1) then
        clegplot(curvenr) = .false.
    endif

    !now find out the appropriate y-axis scaling range for this curve:
    N=npts
    if (ltype == 1) then !ltype is line
        yrayspec(5) = max(min_y, lowerylimit) !0.49_wp
    else
        yrayspec(5) = min_y*0.9_wp
    endif
    yrayspec(3) = maxval(yvals(1:N))
    yrayspec(3) = max(yrayspec(3),yrayspec(5))*1.01_wp
    yrayspec(4) = yrayspec(5)*1.01_wp
    if (yrayspec(3) > yrayspec(4) .and. yrayspec(4) > 1.0E-12_wp .and. ltype == 1 .and. max_y < 10.0_wp) then  !Check for very large values and reduce to a given level...
        if (yrayspec(4) < max_y) then
            yrayspec(3) = yrayspec(4)
        else
            yrayspec(3) = max_y
        endif
    endif
    if (yrayspec(3) > max_y) then
        yrayspec(3) = max_y  !the maximal value on the y-axis...
    else if (yrayspec(3) < min_y) then
        yrayspec(3) = min_y
    endif

    !round up the values to a certain axis step size:
    yy = yrayspec(3)
    if (yy > 1.0_wp .and. yy < 1.001_wp) then
        yy = roundnicely(yy, 2, "down")
    else
        yy = roundnicely(yy, 2, "up")
    endif

    !assign the max value to the corresponding axis:
    if (waxis == 1) then !left axis
        lyraymax(curves_l) = yy !left y-axis
    else !right axis
        ryraymax(curves_r) = yy !right y-axis
    endif
endif !isdata

newplot = .false. !set to read in next curve for the same plot page
logyaxis = .false.
legplot = .true.
bigDsymbs = .false.
pietittxt = ""
pielegtxt = ""

end subroutine DefinePlotData
!=================================================================================================================    


subroutine defbardata(xvals, yvals, lylabel, rylabel, barpatcol, shadingpat, barbordercol)

implicit none
!interface vars:
real(wp),dimension(ncprows),intent(in) :: xvals
real(wp),dimension(:,:),intent(in) :: yvals        !the y-values of the different components
character(len=75),intent(in) :: lylabel, rylabel  !bar graph y-axis labels
integer,dimension(:,:),intent(in) :: barpatcol, shadingpat, barbordercol
!local vars:
integer :: npts, i, nicomp
!.........................................................................................

if (newplot2) then  !this start a new plot page, thus reset the counters and arrays:
    barpmodeno = 0
    compbx = 0.0_wp
    compby = 0.0_wp
    compbpoints = 0.0_wp
    cpietittxt = ""
    cpielegtxt = ""
    cplotbarlabels = .false.
    logxaxis = .false.
    call defplotcolorlists() !call to define color tables for plots
endif

!check data range:
nicomp = 0
do i = 1,2*NCplot
    if (sum(abs(yvals(i,:))) > 1.0E-14_wp) then
        nicomp = i !nicomp+1 !count number of different independent components (nicomp)
    endif
enddo
npts = 0
do i = 1,min( ncprows, size(yvals(1,:)) )
    if (xvals(i) > (-1.0E6_wp) .and. ((xvals(i) > -1.0E-14_wp .or. xvals(i) < (-1.0E-12_wp)) .and. sum(abs(yvals(1:nicomp,i))) > 1.0E-14_wp) &
        & .or. (xvals(i) > 1.0_wp)) then
        npts = i !npts+1 !count number of bars worth plotting...
    endif
enddo

!set the bar attributes in the arrays:
if (npts > 0) then
    barpmodeno = barpmodeno+1 !counts the different number of bar plots on a plot page (the mode number)
    compbx(barpmodeno,1:npts) = xvals(1:npts)
    compby(barpmodeno,1:nicomp,1:npts) = yvals(1:nicomp,1:npts)
    compbpoints(barpmodeno) = npts
    clylabel(barpmodeno) = trim(lylabel)
    crylabel(barpmodeno) = trim(rylabel)
    ncomp(barpmodeno) = nicomp
    cpietittxt(barpmodeno,1:6) = pietittxt(1:6)
    cpielegtxt(barpmodeno,1:2*NCplot) = pielegtxt(1:2*NCplot)
    cplotbarlabels(barpmodeno,1:npts) = plotbarlabels(1:npts)
    cbarpatterncol(barpmodeno,1:nicomp,1:npts) = barpatcol(1:nicomp,1:npts)
    cbarshadpat(barpmodeno,1:nicomp,1:npts) = shadingpat(1:nicomp,1:npts)
    cbarbordcol(barpmodeno,1:nicomp,1:npts) = barbordercol(1:nicomp,1:npts)
endif

curvepointsno = npts
newplot2 = .false. !set to read in next curve for the same plot page
pietittxt = ""
pielegtxt = "" 
plotbarlabels = .false.
legplot = .true.

end subroutine defbardata
!=================================================================================================================  
 

subroutine ShadedArea(xvals, yvalsLow, yvalsUp, areapatcol, shadingpat, areabordercol, legtxt)

implicit none
!interface vars:
real(wp),dimension(:),intent(in) :: xvals, yvalsLow, yvalsUp
integer,intent(in) :: areapatcol, shadingpat, areabordercol !area fill/pattern colour, area fill pattern, area border colour
character(len=75),intent(in) :: legtxt
!local vars:
integer :: npts, i
!.........................................................................................

alegtxt = ""

if (newplot3) then  !this start a new plot page, thus reset the counters and arrays:
    areapmodeno = 0
    areapx = 0.0_wp
    areapyLow = 0.0_wp
    areapyUp = 0.0_wp
    logxaxis = .false.
    call defplotcolorlists() !call to define color tables for plots
endif

!check data range:
npts = 0
do i = 1,min( ncprows, size(yvalsLow(:)) )
    !if (xvals(i) > (-1.0E6_wp) .and. ((xvals(i) > -1.0E-14_wp .or. xvals(i) < (-1.0E-12_wp)) .and. sum(abs(yvalsLow(i))) > 1.0E-14_wp) &
    !    & .or. (xvals(i) > 1.0_wp)) then
    npts = i !npts+1 !count number of points worth plotting...
    !endif
enddo

!set the bar attributes in the arrays:
if (npts > 0) then
    areapmodeno = areapmodeno+1 !counts the different number of bar plots on a plot page (the mode number)
    areapx(areapmodeno,1:npts) = xvals(1:npts)
    areapyLow(areapmodeno,1:npts) = yvalsLow(1:npts)
    areapyUp(areapmodeno,1:npts) = yvalsUp(1:npts)
    areabpoints(areapmodeno) = npts
    areapatterncol(areapmodeno) = areapatcol
    areashadpat(areapmodeno) = shadingpat
    areabordcol(areapmodeno) = areabordercol
    alegtxt(areapmodeno) = trim(legtxt)
endif

curvepointsno = npts
newplot3 = .false. !set to read in next curve for the same plot page
legplot = .true.

end subroutine ShadedArea
!=================================================================================================================


subroutine SetTextlines(text)

implicit none

character(*),intent(in) :: text
!...................................................................

if ((.not. messagebox) .or. (trim(text) == "reset")) then !reset the text strings for the message box in the plotpage (if any)
    messagetxt = "none"
    txtl = 0
endif

if (len_trim(text) > messagetextlength) then
	!$OMP CRITICAL			  
    write(*,*) "The textstring submitted to SetTextlines is too long!"
    write(*,*) "The length is limited to: ",messagetextlength,"characters!"
    write(*,*) ""
	!$OMP end CRITICAL
else if (len_trim(text) > 0 .and. (trim(text) /= "reset")) then
    txtl = txtl+1
    messagetxt(txtl) = trim(text)  !the text string to read in on line txtl for later output on the plotpage in a message box.
    messagebox = .true.
endif
    
end subroutine SetTextlines
!=================================================================================================================


subroutine PlotNow(ptitle, xlabel, lylabel, rylabel, outname_cplot, StandardDev)

use dislin  !the DISLIN plotting library (www.dislin.de)

implicit none
!...
integer :: i, ntop, N, NXL, NYL, NZL, symbsize, NLIN, NMAXLN, ILIN, NXA, &
& NYA, Tsize, Lsize, Namesize, Numbsize, Lcompsize, NTYP, ileg, legstop, &
& bigsymbsize, iclr, ityp, smallpolys, nshadeareas, nxwid, &
& nxpos, nypos, Nheight, ndig
integer,dimension(2) :: NRAY
integer,dimension(4) :: INRAY
integer,dimension(:),allocatable :: barcolset
!...
parameter (NMAXLN = 75)
!...
real(wp) :: StandardDev, xmaxax, xminax, xlab, xstep, yminax, ymaxax, ylab, ystep, logymin, &
    & xminaxin, xmaxaxin, dtiny, rinc, ry, rypos, rypos0
real(wp),dimension(3) :: xrayspec, yrayspec_l, yrayspec_r
real(wp),dimension(:),allocatable :: nxray, nyray												
real(wp),dimension(ncprows) :: xax, yax, E1RAY, E2RAY
!...
character(len=75),intent(in) :: outname_cplot
character(len=75) :: outname_dislin
character(len=(count(clegplot(1:curvenr))*NMAXLN)) :: CBUF
character(len=NMAXLN) :: CSTR
character(len=132),dimension(4) :: ptitle
character(len=75) :: xlabel, lylabel, rylabel
character(len=messagetextlength) :: hline
!...
logical :: ryaxis, logy_l, logy_r, uselegpat
!.........................................................................................

!initialize variables:
xrayspec = 0.0_wp
yrayspec_l = 0.0_wp
yrayspec_r = 0.0_wp
xax = 0.0_wp
yax = 0.0_wp
dtiny = epsilon(1.0_wp)*2.0E2_wp
xminax = 0.0_wp; xmaxax = 10.0_wp; xlab = 0.0_wp; xstep = 1.0_wp;
yminax = 0.0_wp; ymaxax = 10.0_wp; ylab = 0.0_wp; ystep = 1.0_wp;
  
!check whether there is any curve/points to plot at all:
if (curvenr < 1) then !there is nothing to plot
	!$OMP CRITICAL												   
    write(*,*) "There was no curve to plot in this call to PlotNow!"
    write(*,*) "proceeding without producing a plot..."
    write(*,*) ""
	!$OMP end CRITICAL				 
    newplot = .true. !set it back for next plot page
    newplot2 = .true.
    newplot3 = .true.
    messagebox = .false.
    compbars = .false.
    comptemp = .false.
    smallsymbs = .false.
    bigsymbs = .false.
    bigDsymbs = .false.
    specialxrange = .false.
    specialyrange = .false.
    specialyrange2 = .false.
    plotoverlay = .false.
    customCol = .false.
    logyaxis = .false. !set the logarithmic y-axis to false (so use linear).
    logxaxis = .false.
    legplot = .true.
    nogridlines = .false.
    mirrorAxis = .false.
    ticksinside = .false. !if true, tick marks will be plotted inside of axis system
    biglabelfont = .false.
    showcustomcolbar = .false.
    slimplot = .false.
    squareplot = .false.
	equalfontwidth = .false.						
    ptitle = ""
    return !leave the subroutine
endif
  
!paramparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparam
!paramparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparam
!:: adjustable  plot parameters ::  

!squareplot scales the x, y axis length to plot in a square shaped box
    
!Set font sizes etc.: 

if (squareplot .or. slimplot) then
    symbsize = 32 !26 !40 !32 !symbol size
    Tsize = 22   !Title font size
    if (any(len_trim(ptitle(1:4)) > 70)) then
        Tsize = 16
    endif
    Lsize = 19   	!Legends font size
    Lcompsize = 24  !composition plot legend size
    Namesize = 36  	!axis name font
    Numbsize = 32  	!axis numbers size
else
    symbsize = 36 !26 !symbol size
    Tsize = 22   !Title font size
    if (any(len_trim(ptitle(1:4)) > 70)) then
        Tsize = 16
    endif
    Lsize = 19   	!Legends font size
    Lcompsize = 24  !composition plot legend size
    Namesize = 40  	!axis name font
    Numbsize = 32  	!axis numbers size
endif
if (smallsymbs) then
    symbsize = symbsize*2/3 	!symbol size
    bigsymbsize = symbsize*5/4  
else if (bigsymbs) then
    bigsymbsize = 44
    symbsize = bigsymbsize
endif
if (biglabelfont) then
    Namesize = Namesize*4/3							 
endif
if (plotoverlay) then
    Lsize = 18   	!Legends font size
    symbsize = 15 	!small symbols
endif
logymin = 1.0E-3_wp  !the minimum y value when logarithmic y-axis are plottet without special y-range given
    
!paramparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparam
!paramparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparam

    
if (.not. plotoverlay) then !else the DISLIN was already initialized
    !set output file
    call METAFL(metaff) 
    outname_dislin = trim(outname_cplot)//"."//(trim(metaff))
	outname_dislin = trim(defaultOutPath)//trim(outname_dislin) 	!default case
    !special cases (used in AIOMFAC full model):
    if (.not. genericPlotModule) then
        outname_dislin = trim(outname_cplot)//"."//(trim(metaff))
	    if (ternaryLLEoutplot) then
            outname_dislin = "../Output_PhaseSep/"//trim(outname_dislin)
        else if (bargraphplot .or. Cstarplot) then
            outname_dislin = "../Output_Gas2Part/"//trim(outname_dislin)
        else if (directory_outputGouri) then
            outname_dislin = "../Output_Gouri/"//trim(outname_dislin) 
        else
	        outname_dislin = trim(defaultOutPath)//trim(outname_dislin) 	!default case
        endif 
    endif
    call SETFIL(trim(outname_dislin)) !set filename in dislin
    call FILMOD("DELETE") !files will be overwritten if they have the same name..
    call SCRMOD ("NOREV")  !set background color to white and foreground to black
    call SETPAG("DA4L")
    call SCLFAC(1.0_wp)
    call SCLMOD("DOWN")
    call IMGFMT ("RGB")
    call WINSIZ (1600, 1200)
    if (metaff == "png" .or. metaff == "gif") then
        call IMGFMT ("RGB")
        call WINSIZ (1200, 900)
    endif
        
    !Set the overall plot parameters and initialize the plot page:
    call DISINI !initialize dislin
    call SETVLT("RAIN")  !Rainbow colors (256)
    call PAGFLL(255) !background color white
    call SETCLR(0) !black
    if (metaff == "wmf") then
        call WINFNT("Arial")  !FONT in Windows  !Arial  ! Times New Roman
    else if (metaff == "eps" .or. metaff == "ps" .or. metaff == "pdf" .or. metaff == "svg") then
        call PSFONT("Helvetica")  !for ps, eps and pdf output
        call PSMODE('BOTH')  !allow Greek and Italic modes
    else
        call COMPLX
    endif
    call EUSHFT ('GERMAN', '!') !used to print German umlaut characters with Dislin TeX: then !o writes ö
    !date of print and dislin version on the lower right corner
    call HEIGHT(Lsize)
    call PAGHDR("","",2,0) 
    call ERRMOD ("CHECK", "OFF")  !Don't pass warnings of points lying outside of the coordinate systems
    call CHASPC(-0.06_wp)
    call HEIGHT(Namesize) 
    !write some LaTeX code for better axes-text:
    call TEXMOD ("ON")

    !Set Titeltext on line no. 1:
    call TITLIN (trim(ptitle(1)), 1)
    !Set Titeltext on line no. 2:
    call TITLIN (trim(ptitle(2)), 2)
    !Set Titeltext on line no. 3:
    call TITLIN (trim(ptitle(3)), 3)
    !Set Titeltext on line no. 4:
    call TITLIN (trim(ptitle(4)), 4)
    call LINESP(2.6_wp)
    call TITJUS ("LEFT")

    !Set axis labels:
    call GETPAG (NXL, NYL)  !get the page total DIMENSIONs NXL=2970, NYL = 2100
    if (squareplot .and. (.not. slimplot)) then
        call AXSPOS (360, 1400)  !position the axis system (the lower left axis corner)
        call AXSLEN (1000, 1000)  !now set axis length on the page in plot coord.
        if (compbars .or. comptemp) then !for print of composition bars or temperature values
            if (barpmodeno > 1) then !LLE
                call AXSPOS (360, 1400)  !position the axis system (the lower left axis corner)
                call AXSLEN (1000, 1000)  !now set axis length on the page in plot coord.   
            else !only one comp. plot
                call AXSPOS (360, 1400)  !position the axis system (the lower left axis corner)
                call AXSLEN (1000, 1000)  !now set axis length on the page in plot coord.    
            endif
        endif
    else if (slimplot) then
        if (squareplot) then
            call AXSPOS (360, 860)  !position the axis system (the lower left axis corner)
            call AXSLEN (1000, 400)  !now set axis length on the page in plot coord.
        else
            call AXSPOS (360, 1060)  !position the axis system (the lower left axis corner)
            call AXSLEN (1500, 600)  !now set axis length on the page in plot coord.
        endif
    else
        call AXSPOS (360, 1400)  !position the axis system (the lower left axis corner)
        call AXSLEN (1500, 1000)  !now set axis length on the page in plot coord.
        if (compbars .or. comptemp) then !for print of composition bars
            if (barpmodeno > 1) then !LLE
                call AXSPOS (int(NXL*0.10_wp), int(NYL*0.86_wp)-480)  !position the axis system (the lower left axis corner)
                call AXSLEN (int(NXL*0.70_wp), int(NYL*0.64_wp)-480)  !now set axis length on the page in plot coord.   
            else !only one comp. plot
                call AXSPOS (int(NXL*0.10_wp), int(NYL*0.86_wp)-250)  !position the axis system (the lower left axis corner)
                call AXSLEN (int(NXL*0.70_wp), int(NYL*0.64_wp)-250)  !now set axis length on the page in plot coord.    
            endif
        endif
    endif
    call NAMDIS (44, "Y")  !distance between axis text and axis names
    call NAMDIS (44, "X")  
    call HNAME (Namesize)
    call NAME (trim(xlabel), "X")
    call NAME (trim(lylabel), "Y")
endif !plotoverlay
  
!find out if there is a second y-axis...
if (any(cwaxis(1:curvenr) == 2)) then
    ryaxis=.true.
else
    ryaxis=.false.
endif

!1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st
!1st                                   start of 1st y-axis                                    1st
!1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st  

!check for logarithmic y-axis scaling:
logy_l = .false.
logy_r = .false.
do i=1,curvenr
    if (clogyaxis(i) .and. cwaxis(i) == 1) then !log left y-axis
        logy_l = .true.
        call NEGLOG(deps)
    else if (clogyaxis(i) .and. cwaxis(i) == 2) then !log right y-axis
        logy_r = .true.
        call NEGLOG(deps)
    endif
enddo 
 
!the maximum/minimum x-axis/y-axis values and its scaling:
if (.not. specialxrange) then
    xmaxaxin = maxval(xmaxa(1:curvenr))
    xminaxin = minval(xmina(1:curvenr))
    if (fixminxax) then
        xminaxin = min(0.0_wp, xminaxin)
    endif
else
    xminaxin = lowerxlimit
    xmaxaxin = upperxlimit
endif

if (specialxrange) then
    xminax = lowerxlimit
    xmaxax = upperxlimit
    xmaxaxin = upperxlimit
else
    xminax = xminaxin
    xmaxax = xmaxaxin
endif
!automatic scaling of x-axis according to precalculated axis bounds:
xrayspec(1) = xminax
xrayspec(2) = xmaxax
xrayspec(3) = xmaxaxin
if (logxaxis) then
    xrayspec(2) = abs(xrayspec(2))
    if (xrayspec(1) <= 0.0_wp) then
        xrayspec(1) = xrayspec(2)*1.0E-8_wp
    endif
    call AXSSCL("log", "X")  
    call LOGTIC("FULL")    !"AUTO"
else
    call AXSSCL("LIN", "X")
endif
call SETSCL(xrayspec, 3, "X") 

ymaxax = maxval(lyraymax(1:curves_l))
yrayspec_l(3) = ymaxax
yminax = minval(yminal(1:curves_l))
yminax = roundnicely(yminax,2,"down")
if (logy_l .and. logy_r) then
    yminax = min(minval(yminar(1:curves_r)),minval(yminal(1:curves_l)))
    yminax = roundnicely(yminax,2,"down")
    yrayspec_l(3) = max(maxval(ryraymax(1:curves_r)),yrayspec_l(3))
    yrayspec_l(3) = min(10.0_wp**(log10(max(yminax,logymin))+10.0_wp), yrayspec_l(3))  !limit the range on the y-axis to 10 orders of magnitude
else if (logy_l) then
    yrayspec_l(3) = min(10.0_wp**(log10(max(yminax,logymin))+10.0_wp), yrayspec_l(3))  !limit the range on the y-axis to 10 orders of magnitude
endif
if (specialyrange) then
    yminax = lowerylimit
endif

!first find out the y-axis ranges which cover the whole number of curves:
!for the left y-axis:
yrayspec_l(1) = 0.0_wp       !the lowest value on the y-axis
if (yminax < 0.0_wp) then 
    yrayspec_l(1) = yminax
endif
yrayspec_l(2) = min_y
NLIN = 0
call AXSSCL ("LIN", "Y")    !the default choice
if (logy_l) then
    yrayspec_l(2) = max(yrayspec_l(2), logymin)
    call AXSSCL ("log", "Y")
    call LOGTIC ("FULL") !"AUTO"
    if (yminax > -1.0E-12_wp .and. yminax < 1.0E-3_wp .and. yrayspec_l(3) > 1.0E-3_wp) then
        yrayspec_l(1) = logymin
    else
        if (log10(yminax) > 0.0_wp) then !round up
            yrayspec_l(1) = 10.0_wp**(roundnicely(log10(yminax),2,"up"))
        else
            yrayspec_l(1) = 10.0_wp**(roundnicely(log10(yminax),2,"down"))
        endif
    endif
endif
do i=1,curvenr
    if (cwaxis(i) == 1) then
        CSTR = trim(clegtxt(i))
        if (len_trim(CSTR) > 2) then
            NLIN = NLIN+1 !legend line Nr.
        endif
        if (cltype(i) == 1) then !"line"
            csymbtype(i)= -1
        endif
    endif
enddo 

if (specialyrange) then
    yrayspec_l(1) = lowerylimit
    yrayspec_l(2) = upperylimit
    ymaxax = min(ymaxax, upperylimit)
    yminax = min(yminax, lowerylimit)
    if (logy_l) then
        ystep = 1.0_wp
        yrayspec_l(3) = lowerylimit*10.001_wp !cover at least one order of magnitude in log10 scale
        yrayspec_l(2) = roundnicely(upperylimit,1,"up")
        yrayspec_l(2) = max(yrayspec_l(2), 1.0E-15_wp)
    else
        yrayspec_l(3) = upperylimit
    endif
endif
!automatic scaling of axis according to a selected part of the data:
call SETSCL(yrayspec_l, 3, "Y")  
if (mirrorAxis) then !set the exact same values, labels and tick mark properties on x-axis as on y-axis (e.g. for 1:1 plots)
    call SETSCL(yrayspec_l, 3, "X")
endif
    
!plot axis and ticks from (xmin,xmax,first_x_label,increment_x / ymin,ymax,first_y_label,increment_Y):
if (.not. plotoverlay) then
    if (ryaxis) then
        call SETGRF ("NAME", "NAME", "TICKS", "LINE")
    else
        call SETGRF ("NAME", "NAME", "TICKS", "TICKS")
    endif
else
    call SETGRF ("NONE", "NONE", "NONE", "NONE") !set none as the grid and axis will be plotted already before in the calling plot
endif
call HEIGHT(Numbsize)
call CHASPC(-0.07_wp)
call TICKS(4, "XY") !the number of ticks between axis labels (default: 2) 
if (.not. mirrorAxis) then						  
	if (abs(xstep) < 0.01_wp) then
		call LABDIG (4, "X") !number of decimals xaxis
	else if (abs(xstep) < 0.1_wp) then
		call LABDIG (3, "X") !number of decimals xaxis
	else if (abs(xstep) < 0.2_wp) then
		call LABDIG (2, "X") !number of decimals xaxis
	else if (abs(xstep) > 1.0_wp) then
		call LABDIG (1, "X") !number of decimals xaxis
		if (xstep < 10.0_wp) then
			call TICKS(int(xstep), "X")    
		endif
	else
		call LABDIG (1, "X") !number of decimals xaxis
	endif
endif
if ((.not. logy_l) .and. (abs(ymaxax - yminax) < 1.0E-3_wp .or. abs(ymaxax - yminax) > 9.99999E3_wp)) then
    call LABELS("FEXP", "Y") !use powers of 10 "scientific" number format for axis labels
else
    call LABELS("FLOAT", "Y")
endif
!...
!some of the axis use automatic axis scaling and overwrite whatever values are given in the code below:
if (ticksinside) then
    call TICPOS('REVERS', 'XY') !plots tick marks inside plot area
endif					 
call SETVLT("RAIN") !the Dislin rainbow palette
call COLOR("BLACK")
call SOLID
call PENWID(1.0_wp) !set to default value
call LINWID(1) !set to default value
call THKCRV(1) !set to default value
if (mirrorAxis) then !set the exact same values and tick mark properties on x-axis as on y-axis (e.g. for 1:1 plots)
    call TICKS(2, "XY") !the number of ticks between axis labels (default: 2)  
    if (abs(upperylimit - lowerylimit) > 5.0_wp) then
        call INTAX
    endif
    call SETSCL(yrayspec_l, 3, "RESET")
    call GAXPAR (lowerylimit, upperylimit, 'NOEXTEND', 'Y', yminax, ymaxax, ylab, ystep, ndig)
    xminax = yminax
    xmaxax = ymaxax
    xlab = ylab
    xstep = ystep
    call LABDIG (ndig, "XY") !number of decimals xaxis
    call GRAF(xminax, xmaxax, xlab, xstep,  yminax, ymaxax, ylab, ystep) !use parameters defined with GAXPAR
else
    if (logxaxis) then
        !adjust axis scaling parameters with more detail control, first for x-axis:
        call SETSCL(xrayspec, 3, 'RESET')
        call AXSSCL ("log", "X")
        call GAXPAR (log10(lowerxlimit), log10(upperxlimit), 'NOEXTEND', 'X', xminax, xmaxax, xlab, xstep, ndig)
        xlab = ceiling(xminax)
        xstep = min( xstep,  anint((floor(xmaxax) - xlab) / 6.0_wp) )
        call LABDIG (ndig, 'X')     !number of decimals x-axis
        call LOGTIC ("FULL")        !"AUTO"
        call LABELS("log", "X")
        !now for y-axis:
        if (logy_l) then
            call AXSSCL ("log", "Y")
            call GAXPAR (log10(lowerylimit), log10(upperylimit), 'NOEXTEND', 'Y', yminax, ymaxax, ylab, ystep, ndig)
            ylab = ceiling(yminax)
            ystep = min( ystep,  anint((floor(ymaxax) - ylab) / 6.0_wp) )
            call LOGTIC ("FULL")
            call LABELS("log", "Y")
        else
            call SETSCL(yrayspec_l, 3, 'Y')
        endif
        call LABDIG (ndig, 'Y') !number of decimals
		call GRAF(xminax, xmaxax, xlab, xstep,  yminax, ymaxax, ylab, ystep)    !will be scaled automatically using SETSCL    
    else
	    xlab = xminax; ylab = yminax;
	    if (xmaxax > 1.0_wp .and. barpmodeno > 1 .and. xstep > 0.99_wp) then
		    call INTAX
		    yminax = 0.0_wp; ymaxax = 1.0_wp; ylab = 0.0_wp; ystep = 0.1_wp;
		    call GRAF(xminax, xmaxax, xlab, xstep,  yminax, ymaxax, ylab, ystep)
	    else
		    call GRAF(xminax, xmaxax, xlab, xstep,  yminax, ymaxax, ylab, ystep)    !will be scaled automatically using SETSCL
        endif
    endif
endif
    
if (areapmodeno > 0) then
    do nshadeareas=1,areapmodeno
        allocate(nxray(4), nyray(4) )
        smallpolys=0
        do smallpolys = 1,areabpoints(nshadeareas)-1
            ! x-coords for the four corners
            nxray(1) = areapx(nshadeareas,smallpolys)
            nxray(2) = areapx(nshadeareas,smallpolys)
            nxray(3) = areapx(nshadeareas,smallpolys+1)
            nxray(4) = areapx(nshadeareas,smallpolys+1)
            ! y-coords for the four corners
            nyray(1) = areapyLow(nshadeareas,smallpolys)
            nyray(2) = areapyUp(nshadeareas,smallpolys)
            nyray(3) = areapyUp(nshadeareas,smallpolys+1)
            nyray(4) = areapyLow(nshadeareas,smallpolys+1)
            call VLTFIL((trim(defaultColPath)//"RGBDivColTabBlueYellowRed.dat"), "LOAD") 
            call SHDPAT(areashadpat(nshadeareas))
            call SETCLR(areapatterncol(nshadeareas))
            call BARBOR(areabordcol(nshadeareas))
            call RLAREA(nxray, nyray, 4)
        enddo
        deallocate(nxray, nyray)
    enddo
endif
if (.not. plotoverlay) then !else the DISLIN plot title and grid are set elsewhere!
    !Plot now the title:
    call HEIGHT(Tsize)  !Letter height of title
    call CHASPC(-0.04_wp)
	call VLTFIL(trim(defaultColPath)//"RGBColTableRealWinXP8BitAZ.dat", "LOAD") !the Win XP Icon Color Palette plus Grayscale etc. as defined by Andi Zuend
    call SETCLR(0)
    call TITLE
    !plot Gridlines:
    if (.not. nogridlines) then
        call SETVLT("GREY")
        if (logy_l) then
            call SOLID
            call SETCLR(230)
            call PENWID(0.1_wp)                      !the linewidth (allowing values smaller than 1)
            call GRID(0,10)                         !light grey fine gridlines
            call SETCLR(160)
            NRAY(1) = 4                             !pen down
            NRAY(2) = 6                             !pen up
            call MYLINE (NRAY, 2)
            call PENWID(0.2_wp)                      
            call GRID(1,1)
            call LINWID(1)                          !reset the linewidth
        else
            if (logxaxis) then
                call SOLID
                call SETCLR(230)
                call PENWID(0.1_wp)
                call GRID(10,0)                     !light grey fine gridlines on x-axis
                call SETCLR(160)
                NRAY(1) = 4  
                NRAY(2) = 6
                call MYLINE (NRAY, 2)
                call PENWID(0.2_wp)
                call GRID(1,1)
                call LINWID(1)
            else
                call SETCLR(160)
                NRAY(1) = 4
                NRAY(2) = 6
                call MYLINE (NRAY, 2)
                call PENWID(0.2_wp)
                call GRID(1,1)
                call LINWID(1)
            endif
        endif
    endif
    call SOLID
endif !plotoverlay

!calculate and initialize number of lines of 1st legend:
CBUF = ""
CSTR = ""
NLIN = count(clegplot(1:curvenr) .and. cwaxis(1:curvenr) == 1)
if (NLIN > 30) then !dislin currently can only store 30 curve (attributes) for legends
    call LEGINI(CBUF, 30, NMAXLN)
else
    call LEGINI(CBUF, NLIN, NMAXLN)
endif
!check whether a special legend plotting mode is necessary to store curve attributes and change curve color of white curves for visible contrast in legend:
if (any(clegplot(1:curvenr) .and. cwaxis(1:curvenr) == 1 .and. ccolor(1:curvenr) == 255)) then
    uselegpat = .true.
else
    uselegpat = .false.
endif

!plot the curves on the 1st axis ###
ILIN = 0
if (any(.not. clegplot(1:curvenr))) then
    legstop = 2
else
    legstop = 1
endif
if (customCol) then                                     !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
    call VLTFIL(trim(customcolfile), 'LOAD')            !load a customized color palette from the available choices;
else
    call SETVLT('RAIN')                                 !the Dislin rainbow palette
endif
do ileg = 1,legstop                                     !loop two times over the curves, the first time to only plot curves with legend entries.
    do i = 1,curvenr  
        if ((ileg == 1 .and. clegplot(i)) .or. (ileg == 2 .and. (.not. clegplot(i)))) then
            if (cwaxis(i) == 1) then                    !plot it on the left y-axis, the first axis system
                if (ccolor(i) >= 0 .and. ccolor(i) <= 255) then
                    call SETCLR(ccolor(i))              !curve color
                else                                    !special case or error
                    if (ccolor(i) < 0) then             !negative value signals that a WinXP palette color should be used instead of above set color palette
                        call VLTFIL(trim(defaultColPath)//'RGBColTableRealWinXP8BitAZ.dat', 'LOAD')
                        call SETCLR(abs(ccolor(i)))     !set WinXP color
                        if (customCol) then             !reset default color palette choice                          
                            call VLTFIL(trim(customcolfile), 'LOAD')
                        else
                            call SETVLT('RAIN')
                        endif
                    endif
                endif
                call THKCRV(1)                          !curve thickness
                call LINWID(1)                          !the linewidth (will be adjusted by PENWID instead)
                call PENWID(1.0_wp*cth(i))               !the linewidth (allowing also values smaller than 1)
                call INCMRK(0)                          !line
                if (cltype(i) == 2) then                !symbols and line
                    call INCMRK(1)                      !line connecting data points every n=1 datapoints will be plotted
                    call MARKER(csymbtype(i))           !plot symbol type at point coordinates
                    if (smallsymbs) then
                        if (cBigDsymbs(i)) then
                            call HSYMBL(bigsymbsize)    !symbol size
                            call PENWID(1.3_wp*cth(i))   !the linewidth (allowing also values smaller than 1)
                        else
                            call PENWID(0.5_wp*cth(i))   !the linewidth (allowing also values smaller than 1)   
                            call HSYMBL((symbsize*3)/4) !symbol size
                        endif
                    else
                        call PENWID(0.85_wp*cth(i))      !the linewidth (allowing also values smaller than 1)
                        if (cBigDsymbs(i)) then
                            call HSYMBL(bigsymbsize)    !symbol size
                        else   
                            call HSYMBL(symbsize)       !symbol size
                        endif
                    endif
                else if (cltype(i) > 2) then            !symbols without line
                    call INCMRK(-1)                     !no line connecting data points, every n=1 datapoint will be plotted
                    call MARKER(csymbtype(i))           !plot symbol type at point coordinates
                    if (smallsymbs) then 
                        if (cBigDsymbs(i)) then
                            call HSYMBL(bigsymbsize)    !symbol size
                            call PENWID(1.3_wp*cth(i))   !the linewidth (allowing also values smaller than 1)
                        else     
                            call PENWID(0.5_wp*cth(i))   !the linewidth (allowing also values smaller than 1)
                            call HSYMBL(symbsize)       !symbol size
                        endif
                    else
                        if (cBigDsymbs(i)) then
                            call PENWID(1.3_wp*cth(i))   !the linewidth (allowing also values smaller than 1)
                            call HSYMBL(bigsymbsize)    !symbol size
                        else 
                            call PENWID(1.0_wp*cth(i))   !the linewidth (allowing also values smaller than 1)    
                            call HSYMBL(symbsize)       !symbol size
                        endif
                    endif
                endif
                select case (clstyle(i))                !define the line style
                case(1)
                    call LNCAP ("ROUND")                ! ("LONG") ! line caps
                    call SOLID
                    NTYP = 0
                case(2)
                    call LNCAP ("CUT")                  !cut line caps
                    call DASHM
                    NTYP = 5
                case(3)
                    call LNCAP ("ROUND")                !rounded line caps
                    call DASH
                    NTYP = 2
                case(4)                                 !MyDot style
                    NRAY(1) = 1                         !Pen down
                    NRAY(2) = int(1.5_wp*max(cth(i), 4.0_wp)) !NRAY(2) = 2*max(int(cth(i)), 3) !+4  !Pen up
                    call LNCAP ("ROUND")                !rounded line caps
                    call MYLINE (NRAY, 2)               !self defined line style "DOT" with size according to line thickness
                    NTYP = 6
                case(5)
                    call LNCAP ("ROUND")                !rounded line caps
                    call DOT
                    NTYP = 1
                case(6) !. - . - . my dash-dotted line style
                    INRAY(1) = 1                        !Pen down
                    INRAY(2) = 3*max(3, int(cth(i)))    !Pen up
                    INRAY(3) = 4*max(3, int(cth(i)))    !Pen down
                    INRAY(4) = 3*max(3, int(cth(i)))    !Pen up
                    call LNCAP ("ROUND")                !rounded line caps
                    call MYLINE (INRAY(1:4), 4)
                    NTYP = 4
                end select
                call LNJOIN ("SHARP")                   !'SHARP' or 'TRUNC'
                call NOCHEK                             !suppress warning of points lying outside of the plotting area
                !now plot the selected curve:
                N=rownr(i)
                xax = 0.0_wp
                yax = 0.0_wp
                xax(1:N) = cxvals(i,1:N)
                yax(1:N) = cyvals(i,1:N)
                call CURVE(xax(1:N), yax(1:N),N)
                if (ileg == 1) then                     !set the corresponding legend text and color:
                    CSTR = trim(clegtxt(i))						
                    if (len_trim(CSTR) > 1) then
                        if (ILIN < 30) then             !only max. 30 curves attributes can be printed for the legend (current dislin limit)
                            ILIN = ILIN+1
                            !!change color in legends when the line color is white (to have some contrast to the white  background color):
                            if (uselegpat) then
                                if (ccolor(i) == 255) then !white
                                    iclr = 0 !black
                                else 
                                    iclr = ccolor(i)
                                endif 
                                ityp = NTYP             !select a line style from the default styles, closest to the one of the custom types.
                                call LEGPAT(ityp, int(cth(i)), csymbtype(i), iclr, -1, ILIN) !ILIN   
                            endif
                            call LEGVAL (1.0_wp, 'SYMBOL')           !symbol size in legend (default = 0.8_wp)
                            call LEGLIN(CBUF, trim(CSTR), ILIN)     !define the legend line text									
                        endif
                    endif
                endif
                if (any(cxyerror(i,1:N,1:4) > dtiny)) then
                    call LINWID(1) 
                    if (cltype(i) < 2) then             !no symbols, just line
                        call PENWID(min(0.25_wp, max(0.1_wp, cth(i) / 5.0_wp))) 
                        call HSYMBL(max(1,int(cth(i)/3.0_wp))) !set symbol size to tiny size
                        call MARKER(-1)                 !suppress plotting symbol when only error bars are required here... call MARKER(3)
                    else
                        call PENWID(min(0.5_wp, max(0.2_wp, cth(i) / 3.0_wp))) 
                        call MARKER(-1)                 !call MARKER(csymbtype(i))
                    endif
                    !draw vertical (y-axis) error bars:
                    E1RAY(1:N) = cxyerror(i,1:N,3)      !error bars in y-axis
                    E2RAY(1:N) = cxyerror(i,1:N,4)
                    if (any(E1RAY(1:N) > dtiny) .or. any(E2RAY(1:N) > dtiny)) then
                        call BARTYP ('VERT')
                        call ERRBAR(xax(1:N), yax(1:N), E1RAY(1:N), E2RAY(1:N), N)  !vertical error bars
                    endif
                    !draw horizontal (x-axis) error bars:
                    E1RAY(1:N) = cxyerror(i,1:N,1)      !error bars in x-axis
                    E2RAY(1:N) = cxyerror(i,1:N,2)
                    if (any(E1RAY(1:N) > dtiny) .or. any(E2RAY(1:N) > dtiny)) then
                        call BARTYP ('HORI')
                        call ERRBAR(xax(1:N), yax(1:N), E1RAY(1:N), E2RAY(1:N), N)  !horizontal error bars
                    endif
                endif
                call LINWID(1)                          !reset the linewidth default value
                call PENWID(1.0_wp)
            endif !axisc 1 if
        endif !ileg if
    enddo !i, plot curves 1st axis
enddo !ileg
!--.
if (len_trim(CBUF) > 1) then
    call SETVLT("RAIN")
    call SETCLR(0)                                      !black
    call SOLID
    call HEIGHT (Lsize)                                 !Legend font size
    call LEGTIT ("left y-axis:                ")
!    call LEGTIT("")
    call LINESP(2.6_wp)
    if (.not. plotoverlay) then
		if (showcustomcolbar) then
			nxpos = 140
		else
			nxpos = 0
		endif
        if (squareplot) then
            call LEGPOS (1900+nxpos, 40)
        else
			if (curvenr > 1 .and. ryaxis) then
				call LEGPOS (2100+nxpos, 40)
			else
				call LEGPOS (2100+nxpos, 40)
			endif
        endif
    else !plotoverlay
        call LEGPOS (2100, 40)     
    endif
    call FRAME(0)                                       !no box around the legend
    if (customCol) then                                 !read a custom colour table:
        call VLTFIL(trim(customcolfile), "LOAD")
    else
        call SETVLT("RAIN")                             !the Dislin rainbow palette
    endif
    call SETCLR(0) !black
    call LINWID( min( 3, max(1, int(maxval(cth(1:curvenr)))) ) ) !set here curve thickness of curve / points for all legend entries (as no curve-specific choice possible for now).
    call LEGEND (CBUF, 3)
    call LINWID(1) !reset
endif

if (showcustomcolbar) then  !plot a colour bar (z-axis) that has been defined by the used colours
    !plot a custom-made z-axis colour bar based on a collection of stacked coloured rectangles.
    !The colours were determined by the customized colour scaling used in showzcolbar(..) via defzcolval(..).
    call LNCAP ("CUT")
    call GETPOS (NXA, NYA)
    call GETLEN (NXL, NYL, NZL)
    nxpos = NXA + NXL + 100                     !middle of rectangle which is nxwid coordinate points wide in total
    nxwid = 50
    !use 96% of the y-axis length for color bar as default (for continuous scale):
    ry = 0.96_wp*real(NYL, kind=wp)            !scaled preliminary y-axis length
    ntop = min(customzcolno(254), 254)
    N = floor(ry)
    i = (NYL -N)/2
    NYA = NYA - i
    NYL = N
    rinc = real(NYL, kind=wp) / real(ntop, kind=wp) !real-valued increment in plot coordinates
    rypos0 = real(NYA, kind=wp) - 0.5_wp*rinc   !y-coordinate of first point
    rypos = nint(rypos0)
    Nheight = ceiling(rinc) +1              !integer rectangle height
    do i = 1,ntop
        nypos = nint(rypos)
        call POINT(nxpos, nypos, nxwid, Nheight, customzcolno(i))
        rypos = max(rypos - rinc, NYA - NYL + 0.5_wp*rinc)
    enddo
    !plot a secondary y-axis as the actual z-colour axis:
    call SETVLT("RAIN")
    call SETCLR(0) !black
    call SOLID
    call PENWID(1.0_wp)
    call LINWID(1)
    call LNCAP ("LONG")
    call TICKS (2, "Y")
    call HEIGHT(ceiling(Numbsize*0.85))
    call RVYNAM()
    call LABJUS('LEFT', 'Y')
    call CHASPC(-0.07_wp)
    call LABDIS(10, 'Y')
    call NAMDIS(10, 'Y')
    call HNAME (Namesize)
    call LABDIG(1, 'Y')         !number of decimals digits
    if (discreteZcol) then
        call LABDIG(-1, 'Y')    !integer labels
    endif
    call YAXIS(zcolbar0, zcolbar254, zcolbar0, zcolbarstep, NYL, trim(zaxislabel), 0, nxpos+nxwid/2, NYA)
endif

!1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st
!1st                                     end of 1st y-axis                                    1st
!1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st
    

!2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd
!2nd                                    start of 2nd y-axis                                   2nd
!2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd
if (ryaxis .and. (.not. bargraphplot)) then ! 2ndaxis 
  
!set back level to 1 for plotting second axis system
call ENDGRF 
call SETVLT("GREY")
call SETCLR(0)
call SOLID
    
!first find out the y2-axis ranges wich cover the whole number of curves:
ymaxax = maxval(ryraymax(1:curves_r))
yrayspec_r(3) = ymaxax
yminax = minval(yminar(1:curves_r))
yminax = roundnicely(yminax,2,"down")
if (logy_l .and. logy_r) then
    yminax = min(minval(yminar(1:curves_r)),minval(yminal(1:curves_l)))
    yminax = roundnicely(yminax, 2, "down")
    yrayspec_r(3) = max(maxval(ryraymax(1:curves_r)),yrayspec_r(3))
    yrayspec_r(3) = min(10.0_wp**(log10(max(yminax,logymin))+10.0_wp), yrayspec_r(3))  !limit the range on the y-axis to 10 orders of magnitude
else if (logy_r) then
    yrayspec_r(3) = min(10.0_wp**(log10(max(yminax,logymin))+10.0_wp), yrayspec_r(3))  !limit the range on the y-axis to 10 orders of magnitude
endif
  
!for the right y-axis:
yrayspec_r(1) = 0.0_wp !the lowest value on the y-axis
if (yminax < 0.0_wp) then  !xmaxax > 1.1 .and. 
    yrayspec_r(1) = yminax !-10.0
endif
yrayspec_r(2) = 0.01_wp 
call AXSSCL ("LIN", "Y") !the default choice
if (logy_r) then
    yrayspec_r(2) = max(yrayspec_r(2), logymin)
    call AXSSCL ("log", "Y")
    call LOGTIC ("FULL") !"AUTO"
    if (yminax > -1.0E-12_wp .and. yminax < 1.0E-3_wp .and. yrayspec_r(3) > 1.0E-3_wp) then
        yrayspec_r(1) = logymin
    else
        if (log10(yminax) > 0.0_wp) then !round up
            yrayspec_r(1) = 10.0_wp**(roundnicely(log10(yminax),2,"up"))
        else
            yrayspec_r(1) = 10.0_wp**(roundnicely(log10(yminax),2,"down"))
        endif
    endif
endif
do i=1,curvenr
    if (cwaxis(i) == 2) then !right axis
        if (cltype(i) == 1) then !"line"
            yrayspec_r(2) = 1.0_wp
            csymbtype(i)= -1  
        endif
    endif
enddo

if (specialyrange2) then
    yrayspec_r(1) = lowerylimit2
    yrayspec_r(2) = upperylimit2
    ymaxax = min(ymaxax, upperylimit2)
    yminax = min(yminax, lowerylimit2)
    if (logy_r) then
        yrayspec_r(3) = lowerylimit2*10.001_wp !cover at least one order of magnitude in log10 scale
    else
        yrayspec_r(3) = upperylimit2
    endif
endif

if (specialxrange) then
    xminax = lowerxlimit
    xmaxax = upperxlimit
    xrayspec(1) = lowerxlimit
    xrayspec(2) = upperxlimit
    xrayspec(3) = xmaxaxin
endif
if (logxaxis) then
    xrayspec(2) = abs(xrayspec(2))
    if (xrayspec(1) <= 0.0_wp) then
        xrayspec(1) = xrayspec(2)*1.0E-8_wp
    endif
    call AXSSCL ("log", "X")
    call LOGTIC ("FULL") !"AUTO"
else
    call AXSSCL ("LIN", "X")
endif

!automatic scaling of x- and y-axis according to selected data:
call SETSCL(yrayspec_r, 3, "Y")
call SETSCL(xrayspec, 3, "X")
    
!Second axis system for plot:
call NAMDIS (36, "Y")
!Set second y-axis label:
call SETGRF ("LINE", "LINE", "LINE", "NAME")
call TICKS (4, "Y") !the number of ticks between axis labels (default: 2) 

if ((.not. logy_r) .and. (abs(ymaxax-yminax) < 1.0E-3_wp .or. abs(ymaxax-yminax) > 9.99999E3_wp)) then
    call LABELS ("FEXP", "Y") !use powers of 10 "scientific" number format for axis labels
else
    call LABELS ("FLOAT", "Y")
endif
call HNAME (Namesize)
call NAME (trim(rylabel), "Y")
call RVYNAM
call HEIGHT(Numbsize)
call CHASPC(-0.07_wp)
!plot axis and ticks from (xmin,xmax,first_x_label,increment_x / ymin,ymax,first_y_label,increment_Y):
if (ticksinside) then
    call TICPOS('REVERS', 'XY') !plots tick marks inside plot area
endif						
call NEGLOG(deps)
xlab = xminax; ylab = yminax;
if (xmaxax > 1.0_wp .and. barpmodeno > 1 .and. xstep > 0.99_wp) then !LLE relative diff.
    call INTAX
    call GRAF(xminax,xmaxax,xlab,xstep, yminax,ymaxax,ylab,ystep)
else if (xmaxax > 1.0_wp .or. xminax > 0.0_wp .or. (xmaxax-xminax < 0.4_wp)) then !x-axis from xminax to xmaxax
    call GRAF(xminax,xmaxax,xlab,xstep, yminax,ymaxax,ylab,ystep) !will actually be scaled automatically using SETSCL
else
    call GRAF(xminax,xmaxax,xlab,xstep, yminax,ymaxax,ylab,ystep) !will actually be scaled automatically using SETSCL
endif 
   
!plot second Gridlines if the scaling of the right y-axis is different from the left one:
if (abs(maxval(yrayspec_r)-maxval(yrayspec_l)) > 1.0E-12_wp .or. abs(minval(yrayspec_r)-minval(yrayspec_l)) > 1.0E-12_wp) then
    if (.not. nogridlines) then
        call SETVLT("GREY")
        if (logy_r) then
            if (.not. logy_l) then !maybe also then the log-gridlines should be reduced to the main ones to enhance the visualisation!!
                call SETCLR(230)
                NRAY(1) = 12  !Stift unten
                NRAY(2) = 30  !Stift oben
                call MYLINE (NRAY, 2)
                call PENWID(0.1_wp) !the linewidth (allowing values smaller than 1)
                call GRID(0,1)
                call LINWID(1) !reset the linewidth
            endif
        else
            call SETCLR(160)
            NRAY(1) = 12  !pen down
            NRAY(2) = 30  !pen up
            call MYLINE (NRAY, 2)
            call PENWID(0.2_wp) !the linewidth (allowing values smaller than 1)
            call GRID(0,1)
            call LINWID(1) !reset the linewidth
        endif
    endif
    call SOLID
endif

!calculate and initialize number of lines of 2nd legend:
CBUF = ""
CSTR = ""
NLIN = count(clegplot(1:curvenr) .and. cwaxis(1:curvenr) == 2)
if (NLIN > 30) then !dislin currently can only store 30 curve (attributes) for legends
    call LEGINI(CBUF, 30, NMAXLN)
else
    call LEGINI(CBUF, NLIN, NMAXLN)
endif
!check whether a special legend plotting mode is necessary to store curve attributes and change curve color of white curves for visible contrast in legend:
if (any(clegplot(1:curvenr) .and. cwaxis(1:curvenr) == 2 .and. ccolor(1:curvenr) == 255)) then
    uselegpat = .true.
else
    uselegpat = .false.
endif

!plot the curves on the 2nd axis ###
ILIN = 0
if (any(.not. clegplot(1:curvenr))) then
    legstop = 2
else
    legstop = 1
endif
if (customCol) then !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
    call VLTFIL(trim(customcolfile), "LOAD") !load a customized color palette from the available choices;
else
    call SETVLT("RAIN") !the Dislin rainbow palette
endif
do ileg = 1,legstop !loop two times over the curves, the first time to only plot curves with legend entries.
    do i = 1,curvenr  
        if ((ileg == 1 .and. clegplot(i)) .or. (ileg == 2 .and. (.not. clegplot(i)))) then
            if (cwaxis(i) == 2) then !plot it on the right y-axis
                !plot the defined curves:
                call SETCLR(ccolor(i))  !curve colour
                call THKCRV(1)  !curve thickness
                call PENWID(cth(i)) !the linewidth (allowing also values smaller than 1)
                call INCMRK(0) !line
                if (cltype(i) == 2) then !symbols and line
                    call INCMRK(1)  !line connecting data points at every n-th point
                    call MARKER(csymbtype(i))   !plot symbol type at point coordinates
                    if (smallsymbs) then   
                        if (cBigDsymbs(i)) then
                            call PENWID(1.0_wp*cth(i)) !the linewidth (allowing also values smaller than 1)
                            call HSYMBL(bigsymbsize) !symbol size
                        else     
                            call PENWID(0.5_wp*cth(i)) !the linewidth (allowing also values smaller than 1)
                            call HSYMBL((symbsize)/2) !symbol size
                        endif
                    else
                        if (cBigDsymbs(i)) then
                            call PENWID(1.3_wp*cth(i)) !the linewidth (allowing also values smaller than 1)
                            call HSYMBL(bigsymbsize) !symbol size
                        else    
                            call PENWID(1.0_wp*cth(i)) !the linewidth (allowing also values smaller than 1) 
                            call HSYMBL(symbsize) !symbol size
                        endif
                    endif
                else if (cltype(i) > 2) then !symbols without line
                    call INCMRK(-1)  !no line connecting data points
                    call MARKER(csymbtype(i))   !plot symbol type at point coordinates
                    if (smallsymbs) then
                        if (cBigDsymbs(i)) then
                            call HSYMBL(bigsymbsize) !symbol size
                            call PENWID(1.3_wp*cth(i)) !the linewidth (allowing also values smaller than 1)
                        else     
                            call HSYMBL((symbsize)/2) !symbol size
                            call PENWID(1.0_wp*cth(i)) !the linewidth (allowing also values smaller than 1)
                        endif
                    else
                        call PENWID(1.0_wp*cth(i)) !the linewidth (allowing also values smaller than 1)
                        if (cBigDsymbs(i)) then
                            call HSYMBL(bigsymbsize) !symbol size
                        else     
                            call HSYMBL(symbsize) !symbol size
                        endif
                    endif
                endif
                select case (clstyle(i)) !define the line style
                case(1)
                    call LNCAP("ROUND")! ("LONG") !rounded line caps
                    call SOLID
                    NTYP = 0
                case(2)
                    call LNCAP("CUT") !rounded line caps
                    call DASHM
                    NTYP = 5
                case(3)
                    call LNCAP("ROUND") !rounded line caps
                    call DASH
                    NTYP = 2
                case(4) !MyDot style
                    NRAY(1) = 1  !Pen down
                    NRAY(2) = int(1.5_wp*max(cth(i), 4.0_wp)) !NRAY(2) = 2*max(3, int(cth(i))) !+4  !Pen up
                    call LNCAP("ROUND") !rounded line caps
                    call MYLINE (NRAY, 2) !self defined line style "DOT" with size according to line thickness
                    NTYP = 6
                case(5)
                    call LNCAP("ROUND") !rounded line caps
                    call DOT
                    NTYP = 1
                case(6) !. - . - . my dash-dotted line style
                    INRAY(1) = 1  !Pen down
                    INRAY(2) = 3*max(3, int(cth(i)))   !Pen up
                    INRAY(3) = 4*max(3, int(cth(i)))   !Pen down
                    INRAY(4) = 3*max(3, int(cth(i)))   !Pen up
                    call LNCAP("ROUND") !rounded line caps
                    call MYLINE (INRAY(1:4), 4)
                    NTYP = 4
                end select
                call LNJOIN ("SHARP") !'SHARP' or 'TRUNC'
                call NOCHEK  !suppress warning of points lying outside of the plotting area
                !now plot the selected curve:
                N = rownr(i)
                xax(1:N) = cxvals(i,1:N)
                yax(1:N) = cyvals(i,1:N)
                call CURVE(xax(1:N),yax(1:N),N)
                if (ileg == 1) then !set the corresponding legend text and color:
                    CSTR = trim(clegtxt(i))
                    if (len_trim(CSTR) > 1) then
                        if (ILIN < 30) then !only max. 30 curves attributes can be printed for the legend (current dislin limit)
                            ILIN = ILIN+1
                            !!change color in legends when the line color is white (to have some contrast to the white  background color):
                            if (uselegpat) then
                                if (ccolor(i) == 255) then !white
                                    iclr = 0 !black
                                else 
                                    iclr = ccolor(i)
                                endif 
                                ityp = NTYP !select a line style from the default styles, closest to the one of the custom types.
                                call LEGPAT(ityp, int(cth(i)), csymbtype(i), iclr, -1, ILIN) !ILIN   
                            endif
                            call LEGVAL (1.0_wp, 'SYMBOL') !symbol size in legend (default = 0.8_wp)
                            call LEGLIN(CBUF, trim(CSTR), ILIN) !define the legend line text
                        endif
                    endif
                endif
                if (any(cxyerror(i,1:N,1:4) > dtiny)) then
                    call LINWID(1) 
                    if (cltype(i) < 2) then !no symbols, just line
                        call PENWID(min(0.25_wp, max(0.1_wp, cth(i) / 5.0_wp))) 
                        call HSYMBL(max(1,int(cth(i)/3.0_wp))) !set symbol size to tiny size
                        call MARKER(-1) !suppress plotting symbol when only error bars are required here... call MARKER(3)
                    else
                        call PENWID(min(0.5_wp, max(0.2_wp, cth(i) / 3.0_wp))) 
                        call MARKER(-1) !call MARKER(csymbtype(i))
                    endif
                    !draw vertical (y-axis) error bars:
                    E1RAY(1:N) = cxyerror(i,1:N,3) !error bars in y-axis
                    E2RAY(1:N) = cxyerror(i,1:N,4)
                    if (any(E1RAY(1:N) > dtiny) .or. any(E2RAY(1:N) > dtiny)) then
                        call BARTYP ('VERT')
                        call ERRBAR(xax(1:N), yax(1:N), E1RAY(1:N), E2RAY(1:N), N) !vertical error bars
                    endif
                    !draw horizontal (x-axis) error bars:
                    E1RAY(1:N) = cxyerror(i,1:N,1) !error bars in x-axis
                    E2RAY(1:N) = cxyerror(i,1:N,2)
                    if (any(E1RAY(1:N) > dtiny) .or. any(E2RAY(1:N) > dtiny)) then
                        call BARTYP ('HORI')
                        call ERRBAR(xax(1:N), yax(1:N), E1RAY(1:N), E2RAY(1:N), N) !horizontal error bars
                    endif
                endif
                call LINWID(1) !reset the linewidth
                call PENWID(1.0_wp)
            endif !axis if
        endif !legend if
    enddo !plot curves 2nd axis
enddo !ileg
    
!plot legend to the second axis system:
if (len_trim(CBUF) > 1) then
    call SETVLT("RAIN")
    call COLOR("BLACK")
    call SOLID
    call HEIGHT (Lsize)
    call LEGTIT ("right y-axis:             ")
    call LINESP(2.6_wp)
    if (.not. plotoverlay) then
		if (showcustomcolbar) then
            nxpos = 160
        else
            nxpos = 0
        endif			   
        if (squareplot) then
            call LEGPOS (1900, 900)   !position the legend 2 below the other one
        else
            call LEGPOS (2100, 900)   !position the legend on the right of the other one
        endif
    else !plotoverlay
        call LEGPOS (2100, 900)     
    endif
    call FRAME(0)  !no box around the legend
    if (customCol) then !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
        call VLTFIL(trim(customcolfile), "LOAD")
    else
        call SETVLT("RAIN") !the Dislin rainbow palette
    endif
    call SETCLR(0) !black
    call LINWID(min(3, int(maxval(cth(1:curvenr))))) !set here curve thickness of curve / points for all legend entries (as no curve-specific choice possible for now).
    call LEGEND (CBUF, 3)
    call LINWID(1) !reset
endif
  
endif !2ndaxis 
!2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd
!2nd                                       end of 2nd y-axis                                  2nd
!2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd

if(.not. plotoverlay) then
    call SETVLT("RAIN") !the Dislin rainbow palette
    call COLOR("BLACK")
    call SOLID
    call FRAME(1)
    call BOX2D !plots a box around the axis system
endif
    
!Plot a number at coordinates NX,NY
if (trim(ptitle(3)) == "SD of fit:") then !plot standard deviation of fit, SD:
    call GETPOS (NXA, NYA) 
    call GETLEN (NXL, NYL, NZL) 
    call NUMBER (StandardDev, 10, NXA+240, NYA-NYL-126)   !with 10 digits
endif
  
!------------------------------------------------------------------------------------------------
!plot the message box (if there was any line of text defined):
if (messagebox) then
    !set the (initial) plot coordinates and properties of the message box:
    if (messboxpos_x == 0 .and. messboxpos_y == 0) then !the position was not set yet, so use default values
        nxl = 1200
        nyl = 120
    else
        nxl = messboxpos_x
        nyl = messboxpos_y
    endif
	if (equalfontwidth) then
        call PSFONT("Courier")
    else 
        call PSFONT("Helvetica")
    endif						
    call HEIGHT (Lsize)
    call TXTJUS ("LEFT")
    call FRMESS (0) !the thickness of the frame around the box
    !write out message box header horizontal line:
    N = maxval(len_trim(messagetxt(1:txtl)))
	if (equalfontwidth) then
        hline = repeat("-", int(N*2.0_wp)) ! decrease the number of dashes so they fit on the page
    else 
        hline = repeat("-", int(N*5.0_wp))
    endif
    call CHASPC(-0.4_wp)
    call MESSAG (trim(hline), nxl, nyl)
    nyl = nyl + int(Lsize*2.3_wp)
    call CHASPC(-0.01_wp)
    !write out the saved textlines:
    do i=1,txtl 
        call MESSAG ("   "//trim(messagetxt(i)), nxl, nyl)
        nyl = nyl + int(Lsize*2.3_wp)
    enddo
    !write out message box footer horizontal line:
    call CHASPC(-0.4_wp)
    call MESSAG (trim(hline), nxl, nyl)
endif
!------------------------------------------------------------------------------------------------

!end DISLIN plotting:
if (.not. plotoverlay) then   
    call WINMOD("NONE") !for graphic output on XWIN only.:
    call TEXMOD("OFF")
    call DISFIN !end DISLIN
endif

!reset some switches:
newplot = .true.
cxvals = 0.0_wp
cyvals = 0.0_wp
cxyerror = 0.0_wp
call DefinePlotData(cxvals(1,:), cyvals(1,:), cxyerror(1,:,1:4), 1, 0, 1, -1, 3, clegtxt(1), 1.0_wp)    !call here with no actual data to ensure reset of stored curve data.
newplot = .true. !set it back for next plot page
newplot2 = .true.
newplot3 = .true.
messagebox = .false.
compbars = .false.
comptemp = .false.
smallsymbs = .false.
bigsymbs = .false.
specialxrange = .false.
specialyrange = .false.
specialyrange2 = .false.
plotoverlay = .false.
logyaxis = .false. !set the logarithmic y-axis to false (so use linear).
logxaxis = .false.
mirrorAxis = .false.
ticksinside = .false.
biglabelfont = .false.
ptitle = ""
customCol = .false.
bargraphplot = .false.
Cstarplot = .false.
fixminxax = .true.
xorg3bar = .false. 
legplot = .true.
nogridlines = .false.
showcustomcolbar = .false.
slimplot = .false.
squareplot = .false.
equalfontwidth = .false.
if (allocated(barcolset)) deallocate(barcolset)
    
end subroutine PlotNow
!=================================================================================================================


!=================================================================================================================
subroutine PlotPiechart(ptitle,outname_cplot,nneutral,nelectrol)

use dislin  !the DISLIN plotting library (www.dislin.de)

implicit none
!...
integer :: NMAXLN
parameter (NMAXLN = 75)
!...
integer :: i,N,NXL,NYL,symbsize,xpos,ypos,xposstart,yposstart,xlen,ylen,k,lcount
integer :: Tsize,Lsize,Namesize,Numbsize,pielablesize,nneutral,nelectrol
integer,dimension(:),allocatable :: piecolset
!...
real(wp),dimension(ncprows) :: xax,yax
!...
character(len=75),intent(in) :: outname_cplot
character(len=75) :: outname_dislin
character(len=(nneutral+nelectrol+1)*NMAXLN) :: CBUFPIE
character(len=132),dimension(4) :: ptitle
character(len=messagetextlength) :: hline
!...
logical :: ryaxis, plotpielables
!.........................................................................................

!set some variables to zero:
xax = 0.0_wp
yax = 0.0_wp
allocate(piecolset(NCplot))

!check whether there is any curve/points to plot at all:
if (curvenr < 1) then !there is nothing to plot
	!$OMP CRITICAL			  
    write(*,*) "There was no curve to plot in this call to PlotPiechart!"
    write(*,*) "proceeding without producing a plot..."
    write(*,*) ""
	!$OMP end CRITICAL				  
    !set the paramaters back for next plot page
    newplot = .true. 
    newplot2 = .true.
    return !leave the subroutine
endif

!:: adjustable  plot parameters ::  
!Set font sizes etc.: 
symbsize = 32 !40 !26 !symbol size
Tsize = 26   !Title font size
Lsize = 20   !Legends font size
plotpielables = .false. !##############################
pielablesize = 13 !12 !14  !the font size of piechart lables
Namesize = 40  !axis name font
Numbsize = 32  !axis numbers size
call SetCustomColPalette('WinXP') !load custom color palette

!set output file
call METAFL(metaff) 
outname_dislin=trim(outname_cplot)//"."//(trim(metaff))
outname_dislin=trim(defaultOutPath)//trim(outname_dislin)
call SETFIL(trim(outname_dislin)) !set filename in dislin
call FILMOD("DELETE") !files will be overwritten if they have the same name..
call SCRMOD ("NOREV")  !set background color to white and foreground to black
call SETPAG("DA4L")
call SCLFAC(1.0_wp)
call SCLMOD("DOWN")
call IMGFMT ("RGB")
call WINSIZ (1600, 1200)
if (metaff == "png" .or. metaff == "gif") then
    call IMGFMT ("RGB")
    call WINSIZ (1600, 1200)
endif
    
!Set the overall plot parameters and initialize the plot page:
call DISINI !initialize dislin
if (customCol) then !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
    call VLTFIL(trim(customcolfile), "LOAD") !
else
    call SETVLT("RAIN") !the Dislin rainbow palette
endif
call PAGFLL(255) !background color white
call SETCLR(0) !black
call ERRMOD ("CHECK", "OFF")  !Don't pass warnings of points lying outside of the coordinate systems
call TEXMOD ("ON") !allow some LaTeX code for better text

!font style:
if (metaff == "wmf") then
    call WINFNT("Arial")  !FONT in Windows  !Arial  ! Times New Roman
else if (metaff == "eps" .or. metaff == "ps" .or. metaff == "pdf" .or. metaff == "svg") then
    call PSFONT("Helvetica")  !for ps, eps and pdf output
    call PSMODE('BOTH')  !allow Greek and Italic modes
else
    call COMPLX
endif
call EUSHFT ('GERMAN', '!') !used to print German umlaut characters with Dislin TeX: then !o writes ö

!Set Titeltext on line no. 1:
call TITLIN(trim(ptitle(1)), 1)
!Set Titeltext on line no. 2:
call TITLIN(trim(ptitle(2)), 2)
!Set Titeltext on line no. 3:
call TITLIN(trim(ptitle(3)), 3)
!Set Titeltext on line no. 4:
call TITLIN(trim(ptitle(4)), 4)
call LINESP(2.6_wp)
call TITJUS ("LEFT")

call AXSPOS(300, 500)
call AXSLEN(2000, 150)
call SETGRF ("NONE", "NONE", "NONE", "NONE")
call FRAME(0)
call GRAF(0.0_wp, 1.0_wp, 0.0_wp, 0.1_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.2_wp) !pseudo graph for title plot
!plot title above the axis system:
call CHASPC(-0.06_wp)
call SETCLR(0)
call HEIGHT(Tsize)
call TITLE
call ENDGRF

!date of print and dislin version on the lower right corner
call HEIGHT(Lsize)
call PAGHDR("","",2,0) 

!find out if there is a second y-axis...
ryaxis = .false.
do i=1,curvenr
    if (cwaxis(i) == 2) then
        ryaxis = .true.
    endif
enddo

!Set axis labels:
!call GETPAG(NXL, NYL)  !get the page total DIMENSIONs NXL=2970, NYL = 2100
xposstart = 400
yposstart = 1000
ypos = yposstart
xlen = 300
ylen = 300

!change shading pattern cycle:
call PATCYC(1,1)
call PATCYC(2,10)
call PATCYC(3,16)

!define a color set for up to NCplot piechart segments:
piecolset(1:nneutral) = colneutral(1:nneutral)
piecolset(1) = colneutral(1)+2 !lighter blue for water (resp. component 1) color
piecolset(nneutral+1:nneutral+nelectrol) = colelectrol(1:nelectrol)
    
!1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie
!1stPie                               start of 1st set of pie charts                       1stPie
!1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie

!plot the chart on the 1st y position ###
lcount = 0
do i=1,curvenr  
    if (cwaxis(i) == 1) then !plot it on the first "upper" axis system level
        lcount = lcount+1
        !position the chart in the horizontal direction:
        xpos = xposstart + (lcount-1)*(xlen+500)
        call AXSPOS(xpos, ypos)
        call AXSLEN(xlen, ylen)  !now set axis length (pie chart size) on the page in plot coord.
        !set number of segments and the data:
        N = rownr(i)
        yax(1:N) = cyvals(i,1:N)
        !initialize Legend:  
        CBUFPIE = " "
        k = 1 + maxval(len_trim(cpielegtxt(i,1:N)))
        call LEGINI(CBUFPIE,min(N,30), min(NMAXLN, k))
        do k = 1,min(N,30)
            call LEGLIN(CBUFPIE,trim(cpielegtxt(i,k)),k)
        enddo
        !now plot the selected pie chart:
        call SETCLR(0) !black
        call SOLID
        call HEIGHT (pielablesize)  !Legend font size
        call CHASPC(-0.06_wp)
        if (plotpielables) then
            call LABELS("PERCENT","PIE")  !"PERCENT" ; "DATA" ; "BOTH" ; "NONE"
        else
            call LABELS("NONE","PIE") !to skip data in labels
        endif
        call LABDIG(2, "PIE")
        call LABPOS('INTERNAL', 'PIE')  !"EXTERNAL" , 'INTERNAL'
        call FRAME(1) !0 = no frame around segment labels
        !define pie chart properties:
        call CHNPIE("NONE")  !don't cycle piesegment properties
        call SHDPAT(16) !solid fill
        call PIEVEC(0, 'BROKEN')
        k = 1
        call PIECLR (piecolset(k:k+N), piecolset(1:N), N)  !set the piesegment colours: -1 means use current color! 
        !check whether some pie label textboxes have to be shifted in the y-direction due to overlaps:
        pielabshift = .false.
        do k = 1,N-1
            if (sum(yax(k:k+1)) < 0.20_wp) then !shift only one of the labels
                if ((yax(k) < 0.14_wp .and. yax(k) > 0.08_wp .and. yax(k+1) < 0.08_wp) .or. &
                & (yax(k+1) < 0.14_wp .and. yax(k+1) > 0.08_wp .and. yax(k) < 0.08_wp)) then
                    pielabshift(k+1) = .true.
                else if (yax(k) < 0.10_wp .and. yax(k+1) < 0.10_wp) then !shift both labels
                    pielabshift(k:k+1) = .true.
                endif
            endif
        enddo
        do k = 1,N
            call PIECBK(MyPieSettings) !Callbackroutine to shift the pie labels; defined at the end of this file.
        enddo
        call PIEROT(25.0_wp)  !rotate the piechart by 25 degrees counterclockwise
        call LINESP(1.8_wp)
        call PIEVAL(0.5_wp, 'DIST')   !the distance scaling factor for the radial label position; default = 1.0_wp
        call PENWID(0.05_wp) !the linewidth (allowing also values smaller than 1)
        call CIRCSP(1) !smaller arc length to enable plotting of tiny segments
        call LNCAP ("CUT") ! CUT = cut off line caps
        if (plotpielables) then
            call PIEGRF(CBUFPIE,1,yax,N)
        else
            call PIEGRF("",0,yax,N) !for piecharts without labels plotted (when LABELS called with "NONE")
        endif
        !Set and plot pie chart title:
        k = maxval(len_trim(cpietittxt(i,1:4)))
        if (k > 132) then
            write(*,*) 'title text too long'
            read(*,*)
        endif
        call LNCAP("ROUND")
        call PENWID(1.0_wp) !the linewidth (allowing also values smaller than 1)
        call TITLIN (trim(cpietittxt(i,1)),1) !Pie chart title line 1
        call TITLIN (trim(cpietittxt(i,2)),2) !Pie chart title line 2
        call TITLIN (trim(cpietittxt(i,3)),3) !Pie chart title line 3
        call TITLIN (trim(cpietittxt(i,4)),4) !Pie chart title line 4
        call TITJUS ("LEFT")
        call SETCLR(0)
        call HEIGHT(Lsize)
        call LINESP(2.2_wp)
        call TITLE
        call ENDGRF
    endif !axis if
enddo !i, plot curves 1st axis

!1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie
!1stPie                                end of 1st set of pie charts                        1stPie
!1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie	             

!2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie
!2ndPie                               start of 2nd set of pie charts                       2ndPie
!2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie

if (ryaxis) then

!start positions
xposstart = 400
yposstart = 1800
ypos = yposstart

!plot the chart on the 2nd y position ###
lcount = 0
do i=1,curvenr  
    if (cwaxis(i) == 2) then !plot it on the 2nd "lower" y-axis level
        lcount = lcount + 1
        !position the chart in the horizontal direction:
        xpos = xposstart + (lcount-1)*(xlen+500)
        call AXSPOS(xpos, ypos)
        call AXSLEN(xlen, ylen)  !now set axis length (pie chart size) on the page in plot coord.
        !set number of segments and the data:
        N = rownr(i)
        yax(1:N) = cyvals(i,1:N)
        !initialize Legend:  
        CBUFPIE = " "
        k = 1 + maxval(len_trim(cpielegtxt(i,1:N)))
        call LEGINI(CBUFPIE,min(N,30), min(NMAXLN, k))
        do k = 1,min(N,30)
            call LEGLIN(CBUFPIE,trim(cpielegtxt(i,k)),k)
        enddo
        !now plot the selected pie chart:
        call SETCLR(0) !black
        call SOLID
        call HEIGHT (pielablesize)  !segment labels font size
        call CHASPC(-0.06_wp)
        if (plotpielables) then
            call LABELS("PERCENT","PIE")  !"PERCENT" ; "DATA" ; "BOTH" ; "NONE"
        else
            call LABELS("NONE","PIE") !to skip data in labels
        endif
        call LABDIG(2, "PIE")
        call LABPOS('INTERNAL', 'PIE')  !"EXTERNAL" , 'INTERNAL'
        call FRAME(1) !0 = no frame around segment labels
        !define pie chart properties:
        call CHNPIE('NONE')  !don't cycle piesegment properties
        call SHDPAT(16) !solid fill
        call PIEVEC(0, 'BROKEN')  !BROKEN  STRAIGHT
        k = 1
        call PIECLR (piecolset(k:k+N), piecolset(1:N), N)  !set the piesegment colours: -1 means use current color!
        !check whether some pie label textboxes have to be shifted in the y-direction due to overlaps:
        pielabshift = .false.
        do k = 1,N-1
            if (sum(yax(k:k+1)) < 0.20_wp) then !shift only one of the labels
                if ((yax(k) < 0.14_wp .and. yax(k) > 0.08_wp .and. yax(k+1) < 0.08_wp) .or. &
                & (yax(k+1) < 0.14_wp .and. yax(k+1) > 0.08_wp .and. yax(k) < 0.08_wp)) then
                    pielabshift(k+1) = .true.
                else if (yax(k) < 0.10_wp .and. yax(k+1) < 0.10_wp) then !shift both labels
                    pielabshift(k:k+1) = .true.
                endif
            endif
        enddo
        call PIEROT(25.0_wp)
        do k = 1,N
            call PIECBK(MyPieSettings) !Callbackroutine to shift the pie labels; defined at the end of this file.
        enddo
        call PIEVAL(0.5_wp, 'DIST')   !the distance scaling factor for the radial label position; default = 1.0
        call LINESP(1.8_wp)
        call PENWID(0.05_wp) !the linewidth (allowing also values smaller than 1)
        call CIRCSP(1) !smaller arc length to enable plotting of tiny segments
        call LNCAP("CUT")
        if (plotpielables) then
            call PIEGRF(CBUFPIE,1,yax,N)
        else
            call PIEGRF("",0,yax,N) !for piecharts without labels plotted (when LABELS called with "NONE")
        endif
        !Set and plot pie chart title:
        k = maxval(len_trim(cpietittxt(i,1:4)))
        if (k > 132) then
            write(*,*) 'title text too long'
            read(*,*)
        endif
        call PENWID(1.0_wp) !the linewidth (allowing also values smaller than 1)
        call LNCAP("ROUND")
        call TITLIN (trim(cpietittxt(i,1)),1) !Pie chart title line 1
        call TITLIN (trim(cpietittxt(i,2)),2) !Pie chart title line 2
        call TITLIN (trim(cpietittxt(i,3)),3) !Pie chart title line 3
        call TITLIN (trim(cpietittxt(i,4)),4) !Pie chart title line 4
        call TITJUS ("LEFT")
        call SETCLR(0)
        call HEIGHT(Lsize)
        call LINESP(2.2_wp)
        call TITLE
        call ENDGRF
    endif !axis if
enddo !i, plot curves 1st axis

endif !ryaxis

!2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie
!2ndPie                               end of 2nd set of pie charts                         2ndPie
!2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie
  
!------------------------------------------------------------------------------------------------
!plot the message box (if there was any line of text defined):
if (messagebox) then
    call SETCLR(0) !black
    !set the (initial) plot coordinates and properties of the message box:
    if (messboxpos_x == 0 .and. messboxpos_y == 0) then !the position was not set yet, so use default values
        nxl = 1200
        nyl = 120
    else
        nxl = messboxpos_x
        nyl = messboxpos_y
    endif
    call HEIGHT (Lsize)
    call TXTJUS ("LEFT")
    call FRMESS (0) !the thickness of the frame around the box
    !write out message box header horizontal line:
    N = maxval(len_trim(messagetxt(1:txtl)))
    hline = repeat("-", int(N*5.0_wp))
    call CHASPC(-0.4_wp)
    call MESSAG (trim(hline), nxl, nyl)
    nyl = nyl + int(Lsize*2.3_wp)
    call CHASPC(-0.01_wp)
    !write out the saved textlines:
    do i=1,txtl 
        call MESSAG ("   "//trim(messagetxt(i)), nxl, nyl)
        nyl = nyl + int(Lsize*2.3_wp)
    enddo
    !write out message box footer horizontal line:
    call CHASPC(-0.4_wp)
    call MESSAG(trim(hline), nxl, nyl)
endif
!------------------------------------------------------------------------------------------------

!end DISLIN plotting: 
call WINMOD("NONE") !for graphic output on XWIN only.:
call TEXMOD("OFF")
call DISFIN !end DISLIN

!reset some switches:
newplot = .true. !set it back for next plot page
newplot2 = .true.
messagebox = .false.
customCol = .false.
bargraphplot = .false.
Cstarplot = .false.
ptitle = ""
legplot = .true.
slimplot = .false.
deallocate(piecolset)
    
end subroutine PlotPiechart
!=================================================================================================================


!=================================================================================================================
subroutine PlotBarGraphs(ptitle, xlabel, outname_cplot, nneutral, nelectrol)

use dislin
use qsort_c_module

implicit none
!...
integer :: i,N,NXL,NYL,NZL,symbsize,nax,nc,inc
integer :: NMAXLN, ypos, ylen, xlen, ydist
!...
parameter (NMAXLN = 75)
!...
integer :: NXA,NYA,Tsize,Lsize,Namesize,Numbsize,Lcompsize,ifault, &
& nsteps,nneutral,nelectrol,Lheight,allocstat, NXDIG, NYDIG, NZDIG
!...
real(wp) :: xmaxax, xminax, xlab, xstep, logymin, xmaxaxin, xminaxin, &
& yminax, ymaxax, ylab, ystep, davg, sdspace, xmaxc, xminc, yscaler, sc
real(wp),dimension(3) :: yrayspec_l, yrayspec_r
real(wp),dimension(2) :: xrayspec
real(wp),dimension(ncprows) :: xax,yax,xv,y0,xv2
real(wp),dimension(:,:),allocatable :: ysum
real(wp),dimension(:),allocatable :: xray, yray1, yray2
!...
character(len=75),intent(in) :: outname_cplot
character(len=75) :: outname_dislin
character(len=min(max(maxval(ncomp(:))*NMAXLN, (nneutral+nelectrol+1)*NMAXLN), 50*NMAXLN)) :: CBUF
character(len=132),dimension(4) :: ptitle
character(len=75) :: xlabel
character(len=messagetextlength) :: hline
!...
logical :: newgraph

!.........................................................................................
!set some variables to zero:
yrayspec_l = 0.0_wp
yrayspec_r = 0.0_wp
xrayspec = 0.0_wp
xax = 0.0_wp
yax = 0.0_wp
xminax = 0.0_wp; xmaxax = 1.0_wp; xlab = 0.0_wp; xstep = 0.5_wp;
yminax = 0.0_wp; ymaxax = 1.0_wp; ylab = 0.0_wp; ystep = 0.1_wp;
newgraph = .true.
allocate(ysum(2*NCplot,ncprows))
  
!check whether there is any curve/points to plot at all:
if (curvenr < 1 .and. barpmodeno < 1) then !there is nothing to plot
	!$OMP CRITICAL																	  
    write(*,*) "There was no curve to plot in this call to PlotBarGraphs!"
    write(*,*) "proceeding without producing a plot..."
    write(*,*) ""
	!$OMP end CRITICAL												
    newplot = .true. !set it back for next plot page
    newplot2 = .true.
    messagebox = .false.
    compbars = .false.
    comptemp = .false.
    smallsymbs = .false.
    bigsymbs = .false.
    specialxrange = .false.
    specialyrange = .false.
    specialyrange2 = .false.
    plotoverlay = .false.
    logyaxis = .false. !set the logarithmic y-axis to false (so use linear).
    logxaxis = .false.
    ptitle = ""
    customCol = .false.
    bargylimits = .false.
    plotbarborder = .false.
    slimplot = .false.
    slimplotX = .false.
    return !leave the subroutine
endif

!:: adjustable  plot parameters ::  
!squareplot scales the x, y axis length to plot in a square shaped box
    
!Set font sizes etc.: 
symbsize = 38     !symbol size
Tsize = 24        !Title font size
Lsize = 19        !Legends font size
Lcompsize = 24    !composition plot legend size
logymin = 1.0E-3_wp  !the default minimum y value when logarithmic y-axis are plottet
    
!set output file
call METAFL(metaff) 
outname_dislin = trim(outname_cplot)//"."//(trim(metaff))
outname_dislin = trim(defaultOutPath)//trim(outname_dislin)
if (.not. genericPlotModule) then
    outname_dislin = trim(outname_cplot)//"."//(trim(metaff))
    if (ternaryLLEoutplot) then
        outname_dislin = "../Output_PhaseSep/"//trim(outname_dislin)
    else if (bargraphplot .or. Cstarplot) then
        outname_dislin = "../Output_Gas2Part/"//trim(outname_dislin)
    else if (directory_outputGouri) then
        outname_dislin = "../Output_Gouri/"//trim(outname_dislin)   
    else
        outname_dislin = trim(defaultOutPath)//trim(outname_dislin)
    endif 
endif
call SETFIL(trim(outname_dislin)) !set filename in dislin
call FILMOD("DELETE") ! files will be overwritten if they have the same name..
call SCRMOD ("NOREV")  !set background color to white and foreground to black
call SETPAG("DA4L")
call SCLFAC(1.0_wp)
call SCLMOD("DOWN")
call IMGFMT ("RGB")
call WINSIZ (1600, 1200)
if (metaff == "png" .or. metaff == "gif") then
    call IMGFMT ("RGB")
    call WINSIZ (1200, 900)
endif

!axis and label scaling according to the number of bar graphs:
ydist = max(10, 100-(barpmodeno-2)*15)
ylen = (2100-300-300-(barpmodeno-1)*ydist)/barpmodeno !the y-axis length per axis system
if (squareplot) then !use a y:x axis length aspect ratio of 1:3 if possible
    xlen = min(1200, ylen*3) !scale the length of the x-axis charts:
    if (ylen*3 > 1200) then !check limits for desired plot page design
        ylen = 400 != 1200/3
        yscaler = 1.0_wp
    else
        yscaler = real(ylen*3,kind=wp)/1200.0_wp
    endif
    xlen = ylen*3
    Namesize = floor(34*yscaler)
    Numbsize = floor(32*yscaler)
else if (slimplotX) then
    xlen = min(600, nint(ylen*1.5_wp))
    Namesize = max(12, 31-int((barpmodeno-1)*2.0_wp))
    Numbsize = max(11, 28-int((barpmodeno-1)*2.0_wp))
else !use a y:x axis length aspect ratio of 1:4 by default
    if (ylen*4 > 1400) then !check limits for desired plot page design
        ylen = 350 != 1400/4
        yscaler = 1.0_wp
    else
        yscaler = real(ylen*4,kind=wp)/1400.0_wp
    endif
    xlen = ylen*4
    Namesize = floor(36*yscaler)
    Numbsize = floor(34*yscaler)
endif
    
!Set the over all plot parameters and initialize the plot page:
call DISINI
call PAGFLL(255) !background color white
call SETVLT("RAIN")  !Rainbow colors (256)
call SETCLR(0) !black
if (metaff == "wmf") then
    call WINFNT("Arial")  !FONT in Windows  !Arial  ! Times New Roman
else if (metaff == "eps" .or. metaff == "ps" .or. metaff == "pdf" .or. metaff == "svg") then
    call PSFONT("Helvetica")  !for ps, eps and pdf output
    call PSMODE('BOTH')  !allow Greek and Italic modes
else
    call COMPLX
endif
call EUSHFT ('GERMAN', '!') !used to print German umlaut characters with Dislin TeX: then !o writes ö
!date of print, time and dislin version on the lower right corner
call HEIGHT(Lsize)
call PAGHDR("","",2,0) 
call ERRMOD ("CHECK", "OFF")  !Don't pass warnings of points lying outside of the coordinate systems
call CHASPC(-0.06_wp)
call HEIGHT(Namesize) 
!enable some LaTeX code for better axes-text:
call TEXMOD ("ON")

!Set Titeltext on line no. 1:
call TITLIN(trim(ptitle(1)),1)
!Set Titeltext on line no. 2:
call TITLIN(trim(ptitle(2)),2)
!Set Titeltext on line no. 3:							 
call TITLIN(trim(ptitle(3)),3)
!Set Titeltext on line no. 4:							 
call TITLIN(trim(ptitle(4)),4)

!find overall x-axis scaling applicable to all related bar graphs: 
if (specialxrange) then
    xminaxin = lowerxlimit
    xmaxaxin = upperxlimit 
else
    xminaxin = compbx(1,1)
    xmaxaxin = compbx(1,1)
    do nax = 1,barpmodeno
        N = maxval(compbpoints(nax:nax))
        xv = 0.0_wp
        xv(1:N) = compbx(nax,1:N)
        xminaxin = min(xminaxin, minval(xv(1:N), MASK = xv > 0.0_wp))
        xmaxaxin = maxval(xv(1:N), MASK = xv > 0.0_wp)
        xmaxaxin = max(xmaxaxin, min(xminaxin*1.20001_wp, xmaxaxin) )
    enddo   
endif
xminc = xminaxin
xmaxc = xmaxaxin

do nax = 1,barpmodeno !the number of different axis-systems
    !set the axis labels for the current mode
    call NAMDIS (Namesize+10, "Y")  !distance between axis text and axis names
    call NAMDIS (Namesize+10, "X")  
    call HNAME (Namesize)
    call NAME (trim(xlabel), "X")
    call NAME (trim(clylabel(nax)), "Y")
    !position the axis system (the lower left axis corner):
    if (squareplot) then
        ypos = 400 +(nax-1)*ceiling(70.0_wp*yscaler) +ylen*nax
        call AXSPOS(300, ypos)
        call AXSLEN(xlen, ylen)
        newgraph = .true.
    else
        ypos = 400 +(nax-1)*ceiling(70.0_wp*yscaler) +ylen*nax !ypos = 400 +(nax-1)*100 +ylen*nax
        call AXSPOS(300, ypos)
        call AXSLEN(xlen, ylen)
        newgraph = .true.
    endif
    N = maxval(compbpoints(nax:nax))
    xv = 0.0_wp
    !set the number of entries and the x-values:
    xv(1:N) = compbx(nax,1:N)
    call SETCLR(0) !black
    if (newgraph) then !create new axis system 
        xmaxaxin = xmaxc !reset values
        xminaxin = xminc
        !check for appropriate number of steps (scale intervals):
        if (abs(xmaxaxin-xminaxin) < 2.0_wp) then
            if (xmaxaxin < 1.001_wp .and. xmaxaxin > 0.998_wp) then
                xmaxaxin = 1.0_wp
            endif
            nsteps = roundnicely(abs(xmaxaxin-xminaxin)*10.0_wp,2,"up")
            nsteps = max(nsteps, 1)
            if (nsteps > 10) then
                do while (nsteps > 10)
                    nsteps = roundnicely((real(nsteps, kind=wp)/2.0_wp),2,"up")
                enddo
            else if (nsteps < 4) then
                nsteps = roundnicely((real(nsteps, kind=wp)*2.0_wp),2,"up")   
            endif
            call scalenicely(xminaxin, xmaxaxin, nsteps, xminax, xmaxax, xstep, ifault)
        else
            if (abs(xmaxaxin-xminaxin) > 0.08_wp*huge(nsteps)) then
                nsteps = 10
            else
                nsteps = roundnicely(abs(xmaxaxin-xminaxin)*10.0_wp,2,"up")
                nsteps = max(nsteps, 1)
                if (nsteps > 12) then
                    do while (nsteps > 12)
                        nsteps = roundnicely((real(nsteps, kind=wp)/2.0_wp),2,"up")
                    enddo
                else if (nsteps < 4) then
                    nsteps = roundnicely((real(nsteps, kind=wp)*2.0_wp),2,"up")   
                endif
                if (nsteps > int(xmaxaxin)) then
                    nsteps = int(xmaxaxin)
                endif
                do
                    if (mod(int(xmaxaxin-xminaxin),nsteps) == 0) then
                        exit
                    else
                        xmaxaxin = aint(xmaxaxin)+1.0_wp
                    endif
                enddo
            endif
            call scalenicely(xminaxin, xmaxaxin, nsteps, xminax, xmaxax, xstep, ifault)
            if (xmaxaxin-xminaxin+1.0_wp > 2.0_wp) then
                xminax = aint(xminax)
                xmaxax = aint(xmaxax)
                xstep = max(1.0_wp, aint(xstep))
                xstep = min(aint(xstep), aint(xmaxax))
            endif
        endif
        !prepare graph:
        yrayspec_r = 0.0_wp
        yrayspec_r(1) = min(0.0_wp, roundnicely(minval(compby(nax,:,1:N)),2,"down")) !xminax !0.0_wp
        yrayspec_r(2) = min_y
        do i = 1,N
            yrayspec_r(3) = max(yrayspec_r(3), roundnicely(sum(compby(nax,:,i)) - 1.0E-5_wp, 2, "up"))
        enddo
        if (specialyrange) then
            if (bargylimits(nax)) then
                yrayspec_r(1) = lowerylimit
                yrayspec_r(2) = upperylimit
                yrayspec_r(3) = upperylimit
            endif
        endif
        call AXSSCL ("LIN", "Y") !LIN
        ymaxax = maxval(yrayspec_r)
        yminax = minval(yrayspec_r)
        call SETSCL(yrayspec_r, 3, "Y") !scales the axis automatically for GRAF
        !x-axis scaling and mode:
        if (logxaxis) then
            if (specialxrange) then
                xminax = lowerxlimit
                xmaxax = upperxlimit
                xrayspec(1) = xminax*0.999_wp
                xrayspec(2) = xmaxax*1.001_wp 
            else
                xminax = minval(xv(1:N), MASK = xv > 0.0_wp)
                xmaxax = max(xminax*10.0001_wp, maxval(xv(1:N), MASK = xv > 0.0_wp))
                xrayspec(1) = xminax*0.99_wp
                xrayspec(2) = xmaxax*1.01_wp 
            endif
            call AXSSCL ("log", "X") 
            call LOGTIC ("FULL") !"AUTO"
            call SETSCL(xrayspec, 2, "X") !scales the axis automatically for GRAF
        else
            if (specialxrange) then
                xminax = lowerxlimit
                xmaxax = upperxlimit
            endif
            xrayspec(1) = xminax
            xrayspec(2) = xmaxax
            call AXSSCL ("LIN", "X")  
            call SETSCL(xrayspec, 2, "X") !scales the axis automatically for GRAF  
        endif
        
        if (nax == barpmodeno) then !print the x-axis labels only at the last subgraph
            call SETGRF ("NAME", "NAME", "TICKS", "TICKS")
        else
            call SETGRF ("TICKS", "NAME", "TICKS", "TICKS")
        endif
        call LABELS ("FLOAT", "Y")
        call HEIGHT(Numbsize)
        call CHASPC(-0.07_wp)
        call TICKS (2, "Y") !the number of ticks between axis labels (default: 2) 
        call TICKS (4, "X")
        if (abs(xstep) < 0.5_wp) then
            call LABDIG (2, "X") !number of decimals xaxis
        elseif (abs(xstep) > 0.99_wp) then
            call LABDIG (1, "X") !number of decimals xaxis
            if (abs(xmaxax-xminax) > 5.0_wp) then
                call INTAX
            endif
        else
            call LABDIG (1, "X") !number of decimals xaxis
        endif
        if (abs(ymaxax - yminax) < 1.0E-3_wp .or. abs(ymaxax - yminax) > 9.99999E3_wp) then
            call LABELS ("FEXP", "Y") !EXP will use powers of 10 "scientific" number format for axis labels; FEXP will use "1.0E-8" style.
            call LABDIG (1, "Y") !set the number of decimals shown on y-axis with exponential formatting
        else
            call LABELS ("FLOAT", "Y") !FLOAT
            call LABDIG (1, "Y") !set the number of decimals shown on y-axis
        endif
        call GETDIG (NXDIG, NYDIG, NZDIG)
        call LABDIS(Numbsize-8, "XYZ") !default is 24
        call TICLEN (Numbsize-8, ((Numbsize-8)*2)/3) !major and minor tick lengths
        ymaxax = 1.0_wp; yminax = 0.0_wp; ystep = 0.2_wp;
        xlab = xminax; ylab = yminax; !only to initialize when SETSCL is used
        call GRAF(xminax, xmaxax, xlab, xstep,  yminax, ymaxax, ylab, ystep)
        call GETDIG (NXDIG, NYDIG, NZDIG)
        !initialize legend:  
        CBUF = " "
        inc = ncomp(nax)
        call LEGINI(CBUF, min(inc,30), NMAXLN) !dislin can only stor 30 curve attributes currently
        inc = 0
        do nc = 1,ncomp(nax)
            inc = inc+1
            if (nc <= 30) then
                call LEGLIN(CBUF,trim(cpielegtxt(nax,nc)),inc)
            endif
        enddo
        call LEGTIT(trim(cpietittxt(nax,nax)))
        !position legend:
        call GETPOS (NXA, NYA)
        call GETLEN (NXL, NYL, NZL)
        if (squareplot) then
            call LEGPOS(NXA+NXL+50, NYA-NYL)
        else
            call LEGPOS(NXA+NXL+100, NYA-NYL)  !position the legend on the right side of the axis system
        endif
        if (nax == 1) then !plot title
            call SETCLR(0)
            !check the length of the longest title line to adjust font size:
            i = maxval(len_trim(ptitle(1:4)))
            if (i > 100) then
                Lheight = max(int(100.0_wp*real(Tsize, kind=wp) / real(i, kind=wp)), 2)
                call HEIGHT(Lheight)
            else
                call HEIGHT(Tsize)
            endif
            call LINESP(2.6_wp)
            call TITJUS ("LEFT")
            call TITLE
        endif
        newgraph = .false.
    endif
    !prepare the data for bars:
    call BARTYP ('VERT')
    ysum = 0.0_wp !clear
    y0 = 0.0_wp
    do nc = 1,ncomp(nax)
        if (nc == 1) then
            ysum(nc,1:N) = compby(nax,nc,1:N) !x(component nc)
        else
            ysum(nc,1:N) = ysum(nc-1,1:N) + compby(nax,nc,1:N) !x(component nc)
        endif
    enddo !nc
    !calculate average spacing between adjacent bars:
    if (N > 1) then
        !sort xv array:
        xv2 = 0.0_wp
        if (N > 20) then
            xv2(1:N-1) = xv(1:N-1)
            xv2(N) = xv2(N-1)
        else
            xv2(1:N) = xv(1:N)
        endif
        call QsortC(xv2(1:N))
        davg = 0.0_wp
        sdspace = 0.0_wp
        if (logxaxis) then
            do i = 1,N-1
                davg = davg + abs(log10(xv2(i))-log10(xv2(i+1)))
            enddo
        else
            do i = 1,N-1
                davg = davg + abs(xv2(i)-xv2(i+1))
            enddo
        endif
        davg = davg/real(N-1, kind=wp)
        !standared deviation in x-axis spacing between bars:
        if (logxaxis) then
            do i = 1,N-1
                sdspace = sdspace +(abs(log10(xv2(i))-log10(xv2(i+1)))-davg)**2
            enddo
        else
            do i = 1,N-1
                sdspace = sdspace +(abs(xv2(i)-xv2(i+1))-davg)**2
            enddo
        endif
        sdspace = sqrt(sdspace/real(N-1, kind=wp))
    else
        davg = 36.0_wp
    endif
    !-+-+specify bar width:
    if (N > 3) then
        if (sdspace < davg*0.1_wp .or. logxaxis) then !bars are spaced in regular intervals --> use a specified width
            yrayspec_r(1) = davg
            yrayspec_r(2) = -0.95_wp*yrayspec_r(1)*abs(davg*N)/abs(xmaxax-xminax)*(real(xlen, kind=wp)/abs(xmaxax-xminax))
            yrayspec_r(2) = roundnicely(yrayspec_r(2),2,"up")
            if (logxaxis) then 
                yrayspec_r(2) = min(max(yrayspec_r(2)*0.3_wp, -80.0_wp), -24.0_wp)
            endif 
            call BARWTH(yrayspec_r(2))
        else !since bars are spaced at variable x-axis interval distances, use a different approach.
            call BARWTH(0.75_wp)
            call BARMOD ('VARIABLE', 'WIDTH')
        endif
    else
        call BARWTH(-0.5_wp*36.0_wp)
    endif
    !-+-+
    if (customCol) then !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
        call VLTFIL(trim(customcolfile), "LOAD")
    else
        call SETVLT("RAIN") !the Dislin rainbow palette
    endif
    !Plot the Bars:
    allocate(xray(N), yray1(N), yray2(N), STAT=allocstat)
    xray = xv(1:N)
    yray1 = 0.0_wp
    yray2 = 0.0_wp
    inc = ncomp(nax)
    call HEIGHT (Lcompsize/2-1)
    if (logxaxis) then
        i = 0 !for debugging set breakpoint
        where (xray(:) < xminax)
            xray(:) = 0.1_wp*xminax !make sure those values are set to a positive value, but smaller than the lower x-axis limit.
        endwhere
    endif
    do nc = 1,ncomp(nax)
        call THKCRV(1)  !curve thickness
        call LINWID(1)  !the linewidth (will be adjusted by PENWID instead)
        call PENWID(0.3_wp)  !the linewidth (allowing also values smaller than 1)
        call LNCAP ("LONG") !"CUT" "ROUND" or "LONG" : CUT = cut off line caps
        call LNJOIN ("SHARP") !'SHARP' or 'TRUNC'	
        if (nc == 1) then
            yray1 = y0(1:N)
            yray2 = ysum(nc,1:N)
            if (ncomp(nax) == 1 .and. any(cplotbarlabels(nax,1:N))) then !only one component/bar per x-value
                do i = 1,N
                    if (cplotbarlabels(nax,i)) then !plot bar labels: the y-value
                        call LABELS('DELTA', 'BARS')
                        call LABPOS('OUTSIDE', 'BARS')
                    else
                        call LABELS('NONE', 'BARS')
                    endif
                    call SHDPAT(cbarshadpat(nax,nc,i)) !16=solid fill, 8=hatched pattern
                    call SETCLR(cbarpatterncol(nax,nc,i))
                    if (plotbarborder) then
                        call PENWID(0.5_wp) 
                        call LNCAP ("CUT") 
                        call BARBOR(cbarbordcol(nax,nc,i)) !the color value of the bar border line, -1 means use current color
                    endif
                    call BARS(xray(i),yray1(i),yray2(i),1)  !Plot BAR i
                enddo !i
            else
                call SHDPAT(cbarshadpat(nax,nc,1)) !16=solid fill, 8=hatched pattern
                call SETCLR(cbarpatterncol(nax,nc,1))
                if (plotbarborder) then
                    call PENWID(0.5_wp) 
                    call LNCAP ("CUT") 
                    call BARBOR(cbarbordcol(nax,nc,1)) !the color value of the bar border line, -1 means use current color
                endif
                call BARS(xray,yray1,yray2,N)  !Plot BARS
            endif
        else !nc > 1
            yray1 = ysum(nc-1,1:N)
            yray2 = ysum(nc,1:N)
            if (inc == 1 .and. any(cplotbarlabels(nax,1:N))) then !only one component/bar per x-value
                do i = 1,N
                    if (cplotbarlabels(nax,i)) then !plot bar labels: the y-value
                        call LABELS('DELTA', 'BARS')
                        call LABPOS('OUTSIDE', 'BARS')
                    else
                        call LABELS('NONE', 'BARS')
                    endif
                    call SHDPAT(cbarshadpat(nax,nc,i)) !16=solid fill, 8=hatched pattern
                    call SETCLR(cbarpatterncol(nax,nc,i))
                    if (plotbarborder) then
                        call PENWID(0.5_wp) 
                        call LNCAP ("CUT") 
                        call BARBOR(cbarbordcol(nax,nc,i)) !the color value of the bar border line, -1 means use current color
                    endif
                    call BARS(xray(i),yray1(i),yray2(i),1)  !Plot BAR i
                enddo !i
            else !inc > 1
                call SHDPAT(cbarshadpat(nax,nc,1)) !16=solid fill, 8=hatched pattern
                call SETCLR(cbarpatterncol(nax,nc,1))
                if (plotbarborder) then
                    call PENWID(0.5_wp) 
                    call LNCAP ("CUT")
                    call BARBOR(cbarbordcol(nax,nc,1)) !the color value of the bar border line, -1 means use current color
                endif
                call BARS(xray,yray1,yray2,N)  !Plot BARS
            endif
        endif !nc
    enddo !nc
    !Plot Bar-Legends:
    call PENWID(1.0_wp) !the linewidth (allowing also values smaller than 1)
    sc = real(min(maxval(ncomp(:)),30), kind=wp)
    Lheight = ceiling( ylen*0.51_wp / sc ) !estimate the legend text height
    sc = min( 1.0_wp, 4.0_wp/real(barpmodeno, kind=wp) )
    Lheight = min(ceiling(Lcompsize*sc), Lheight)
    call HEIGHT (Lheight)  
    call LINESP (1.8_wp)
    call SETCLR(0)
    call FRAME(0) !no frame
    call LEGEND(CBUF,7)
    call SOLID
    call FRAME(1)
    call BOX2D !plots a box around legend key
    call ENDGRF !set back level to 1
    deallocate(xray, yray1, yray2, STAT=allocstat)
enddo !nax

call SETCLR(0) !black
call SOLID
call FRAME(1)
call BOX2D !plots a box around the axis system
call ENDGRF !set back level to 1 for plotting 3rd axis system!
  
!------------------------------------------------------------------------------------------------
!plot the message box (if there was any line of text defined):
if (messagebox) then
    !set the (initial) plot coordinates and properties of the message box:
    if (messboxpos_x == 0 .and. messboxpos_y == 0) then !the position was not set yet, so use default values
        nxl = 1200
        nyl = 120
    else
        nxl = messboxpos_x
        nyl = messboxpos_y
    endif
    call HEIGHT (int(Lsize*0.75))
    call TXTJUS ("LEFT")
    call FRMESS (0) !the thickness of the frame around the box
    !write out message box header horizontal line:
    N = maxval(len_trim(messagetxt(1:txtl)))
    if (N > 0) then
        i = min(int(N*5.0_wp), 256)
        hline = repeat("-", i)
        call CHASPC(-0.4_wp)
        call MESSAG (trim(hline), nxl, nyl)
        nyl = nyl + int(Lsize*2.3_wp)
        call CHASPC(-0.01_wp)
        !write out the saved textlines:
        do i=1,txtl 
            call MESSAG ("   "//trim(messagetxt(i)), nxl, nyl)
            nyl = nyl + int(Lsize*1.9_wp)
        enddo
        !write out message box footer horizontal line:
        call CHASPC(-0.4_wp)
        call MESSAG (trim(hline), nxl, nyl)
    endif
endif
!------------------------------------------------------------------------------------------------

!end DISLIN plotting:  
call WINMOD("NONE") !for graphic output on XWIN only.:
call TEXMOD("OFF")
call DISFIN !end DISLIN

!reset some switches:
newplot = .true. !set it back for next plot page
newplot2 = .true.
messagebox = .false.
compbars = .false.
comptemp = .false.
smallsymbs = .false.
bigsymbs = .false.
specialxrange = .false.
specialyrange = .false.
specialyrange2 = .false.
plotoverlay = .false.
logyaxis = .false. !set the logarithmic y-axis to false (so use linear).
logxaxis = .false.
ptitle = ""
customCol = .false.
bargraphplot = .false.
Cstarplot = .false.
legplot = .true.
bargylimits = .false.
plotbarborder = .false.
slimplot = .false.
slimplotX = .false.
squareplot = .false.
deallocate(ysum)
    
end subroutine PlotBarGraphs
!=================================================================================================================

!=================================================================================================================
!subroutine to read and initialize a color table for use of up to NCplot plot colors 
!stored in colneutral and colelectrol
subroutine defplotcolorlists()

implicit none

integer :: i, j, n, unitf, istat
character(len=40) :: colfile, dummy
logical :: exists
!................

!check if reading the data from file is necessary at this call (if done previously in this session, it can be skipped):
if (colneutral(1) == 15 .and. colneutral(1) /= colneutral(NCplot) .and. colelectrol(1) == 11) then
    !data already read... so continue  
else !read data
    !define electrolyte color list (for now not read from a file):
    colelectrol(1:12) = [ 11, 8, 162, 21, 241, 2, 220, 194, 25, 4, 30, 253 ] 

    !read neutral component colors from a text-file:
    colneutral = 250 !initialize 
    colfile = trim(defaultColPath)//"plotcolortable1.dat"
    inquire(FILE=trim(colfile), EXIST=exists)
    if (exists) then
        open (NEWUNIT=unitf, FILE=trim(colfile), STATUS="OLD")
        read(unitf,*) dummy !read first line
        read(unitf,*) !read empty line
        read(unitf,*) dummy !read line 3
        do i = 1,min(NCplot, 40)  !read in the data
            read(unitf,*,IOSTAT=istat) n, colneutral(i) !read data using implied do-loop for the values of tLLEwtf_d(1:20,i)
            if (istat /= 0) then
                exit
            endif
        enddo
        close(unitf)
    else
        write(*,*) ""
        write(*,*) "ERROR in 'ModCreatePlot': file listing plot colors not found!"
        write(*,*) ""
        read(*,*) !wait for user action
        colneutral(1:NCplot) = 10 
    endif
    !potentially, add more colors if number of components is large
    if (NCplot > 40) then
        j = 126
        do i = 41,NCplot
            colneutral(i) = j
            j = j +4
            if (j > 250) then
                j = 3
            endif
        enddo
    endif
endif

end subroutine defplotcolorlists
!=================================================================================================================

!=================================================================================================================
!subroutine to retrieve a colour scale number (colno) from the custom set colour bar 
!(i.e. for coloured data points on a virtual z-axis as part of a 2-D plot)
subroutine SetCustomColPalette(colpalette)

implicit none
!interface variables:
character(len=*),intent(in) :: colpalette
!...............................

select case(trim(colpalette))
case('WinXP')
    customcolfile = 'RGBColTableRealWinXP8BitAZ.dat' 
case("Viridis")
    customcolfile = 'RGBColTabViridisPurpleBlueGreenYellow.dat'
case('BlackPurpleRedYellow')
    customcolfile = 'RGBColTabBlackPurpleRedYellow.dat'
case('BluePurpleOrangeYellow')
    customcolfile = 'RGBColTabBluePurpleOrangeYellow.dat'
case('RedBlue')
    customcolfile = 'RedBlueColTableReal_dislin.dat'
case('BalanceBlueGreen')
    customcolfile = 'RGBDivColTabBalanceBlueGreen.dat'
case('BalanceBlueRed')
    customcolfile = 'RGBDivColTabBalanceBlueRed.dat'
case('BlueGreenWhiteOrangePurple')
    customcolfile = 'RGBDivColTabBlueGreenWhiteOrangePurple.dat'
case('BlueYellowRed')
    customcolfile = 'RGBDivColTabBlueYellowRed.dat'
case('CoolWarm')
    customcolfile = 'RGBDivColTabCoolWarm.dat'
case('PurpleOrange')
    customcolfile = 'RGBDivColTabPurpleOrange.dat'
end select

customcolfile = trim(defaultColPath)//trim(customcolfile)
customCol = .true.

end subroutine SetCustomColPalette
!=================================================================================================================

!=================================================================================================================
!subroutine to retrieve a colour scale number (colno) from the custom set colour bar 
!(i.e. for coloured data points on an extra z-axis as part of a 2-D plot)
pure elemental subroutine defzcolval(zval, colno, colbarlowlim, colbaruplim)

implicit none
!interface variables:
real(wp),intent(in) :: zval          !the "z-axis" value for which a colour number should be determined
integer,intent(out) :: colno     !the output colour index on 1 to 254 scale; i.e. the index of the colour number from a (custom) colour table used;
real(wp),intent(in) :: colbarlowlim  !the lowest value mapped to the colour bar, i.e. the z value corresponding to colour 1
real(wp),intent(in) :: colbaruplim   !the highest value mapped to the colour bar, i.e. the z value corresponding to colour 254
!local variables:
!integer :: ncl, ncols
real(wp) :: colslope, intercept
!...............................

!linear colour scale between colors 1 and 254; this range because colours 0 and 255 are reserved for special use; 
colslope = 253.0_wp/(colbaruplim - colbarlowlim)
intercept = 1.0_wp -colslope*colbarlowlim
colno = nint(intercept + colslope *max(min(zval, colbaruplim), colbarlowlim) )     

end subroutine defzcolval
!=================================================================================================================

!=================================================================================================================
!subroutine to define plotting properties of secondary z-axis colour bars (customizable)
subroutine showzcolbar(colbarlowlim, colbaruplim, zlabel, discrete)

implicit none
!interface variables:
real(wp),intent(in) :: colbarlowlim      !the lowest value mapped to the colour bar, i.e. the z value corresponding to colour 1
real(wp),intent(in) :: colbaruplim       !the highest value mapped to the colour bar, i.e. the z value corresponding to colour 254
character(len=75),intent(in) :: zlabel  !z axis label text string
logical,intent(in) :: discrete       !.true. for discrete color scale (e.g. with an integer number of distinct colors; set to .false. for contiuous scale;
!local variables:
integer :: i, idstep, ih, kh, il, isc
integer :: ncols                     !number of distinct colors in case of discrete scale;
integer,dimension(254) :: zcolno
real(wp) :: zinterval, zdstep, rtop, rscale, rdstep, rih
real(wp),dimension(254) :: zv
!...............................

!set private module variables to be used with ZAXIS or ZAXLG in 'Plotnow':
zcolbar0 = colbarlowlim
zcolbar254 = colbaruplim
zinterval = (colbaruplim - colbarlowlim) / 253.0_wp		!continuous 254 color scale
discreteZcol = discrete
if (discrete) then
    ncols = nint(colbaruplim - colbarlowlim)
    zcolbarstep = (colbaruplim - colbarlowlim) / real(ncols, kind=wp)
    zdstep = 254.0_wp / real(ncols, kind=wp)
    if (ncols > 10) then
        i = ncols
        do while (i > 10)
            i = i / 2
            zcolbarstep = zcolbarstep*2.0_wp  
        enddo
    endif
else
    zcolbarstep = (colbaruplim - colbarlowlim) / 10.0_wp       !define the stepping between tick marks
endif
zaxislabel = zlabel
showcustomcolbar = .true. 
!define array of colour values based on custom function defzcolval(...); 
!to be shown as coloured rectangles as a custom-made colour bar:
do i = 1,254
    zv(i) = colbarlowlim + real(i-1, kind=wp)*zinterval
    call defzcolval(zv(i), zcolno(i), colbarlowlim, colbaruplim)
enddo

if (discrete) then    
    idstep = max(254 / ncols, 1)                    !integer division
    rdstep = 254.0_wp / real(ncols, kind=wp)
    rtop = real(idstep*ncols, kind=wp)         !maximum values for ih in this discrete case
    rscale = 254.0_wp / rtop
    ih = 0
    kh = 0
    rih = 0.0_wp
    do i = 1,ncols
        il = ih + 1
        kh = kh + idstep
        rih = rih + rdstep
        ih = max(1, min(nint(rih), 254) )
        isc = max( 1, min(nint(kh*rscale), 254) )
        customzcolno(il:ih) = zcolno(isc)           !discrete color intervals
    enddo
    customzcolno(ih+1:) = zcolno(ih)                !zcolno(ih+1)
else !continuous   
    customzcolno(1:254) = zcolno(1:254)             !as smooth as 254 colors allow...
endif
	 
end subroutine showzcolbar
!=================================================================================================================

!=================================================================================================================
!** user defined subroutine to set some label position options when plotting Pie-Charts **
subroutine MyPieSettings(ISEG, XDAT, XP, NRAD, NOFF, ANG, NVX, NVY, IDRW, IANN)

implicit none
!..
integer,intent(in) :: ISEG
integer,intent(inout) :: NRAD, NOFF, NVX, NVY, IDRW, IANN
!..
real(wp),intent(in) :: XDAT, XP
real(wp),intent(inout) :: ANG
!..........................................................

if (pielabshift(ISEG)) then
    if (ISEG > 1) then
        if (mod(ISEG,2) == 0 .and. (.not. pielabshift(ISEG-1))) then
            NVY = -40  !-40
        else if (pielabshift(ISEG-1) .and. pielabshift(ISEG+1)) then
            NVY = 10
        else
            NVY = 30  !-40
        endif
    else
        NVY = 30  !-40
    endif
endif

end subroutine MyPieSettings
!=================================================================================================================

end module ModCreatePlot