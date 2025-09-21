!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module including subroutines to produce custom 2-D x-y scatter plots, curves,      * 
!*   bar graphs, and pie charts using the DISLIN graphics library.                      *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, Switzerland; 2007 - 2009                                       *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA; 2009 -2012    *
!*                                                                                      *
!*   -> created:        2007                                                            *
!*   -> latest changes: 2019/10/27                                                      *
!*                                                                                      *
!****************************************************************************************

MODULE Mod_createplot

USE dislin

IMPLICIT NONE   
!...
INTEGER(4),PUBLIC :: NTCP,nrows,ncurves
INTEGER(4),PRIVATE :: curvenr,curves_l,curves_r,txtl
INTEGER(4),PRIVATE :: messagetextlength,messagelinesmax
!...
PARAMETER (NTCP = 25)  !NTCP parameter can be changed if more variables / individual components are required
PARAMETER (nrows = 1000)
PARAMETER (ncurves = 40000)  !the maximum number of different curves in a plot
PARAMETER (messagetextlength = 200)   !the maximum number of characters per line in the text message
PARAMETER (messagelinesmax = 40)   !the maximum number of lines in the message box
!...
INTEGER(4),PUBLIC :: messboxpos_x, messboxpos_y, barpmodeno, curvepointsno
INTEGER(4),DIMENSION(6),PRIVATE :: compbpoints, ncomp
INTEGER(4),DIMENSION(ncurves),PRIVATE :: cwaxis, ccolor, cltype, rownr, csymbtype, clstyle
INTEGER(4),DIMENSION(6,2*NTCP,nrows),PRIVATE :: cbarbordcol, cbarshadpat, cbarpatterncol
INTEGER(4),DIMENSION(NTCP),PUBLIC :: colneutral
INTEGER(4),DIMENSION(6),PUBLIC :: colelectrol 
INTEGER(4),DIMENSION(254),PRIVATE :: customzcolno
!...
REAL(8),PARAMETER,PUBLIC :: rdefault = -8.88888D280 !an unusual, small (negative) value...
REAL(8),PUBLIC :: max_y, min_y, lowerylimit, upperylimit, lowerylimit2, upperylimit2, lowerxlimit, upperxlimit
!private reals:
REAL(8),PARAMETER,PRIVATE :: DEPS = EPSILON(1.0D0)
REAL(8),PRIVATE :: zcolbar0, zcolbar254, zcolbarstep
REAL(8),DIMENSION(ncurves),PRIVATE :: ryraymax,lyraymax,xmaxa,xmina,yminal,yminar,xmaxb,xminb,cth
REAL(8),DIMENSION(ncurves,nrows),PRIVATE :: cxvals,cyvals,cyerr
REAL(8),DIMENSION(6,nrows),PRIVATE :: compbx
REAL(8),DIMENSION(6,2*NTCP,nrows),PRIVATE :: compby
!...
CHARACTER(LEN=3),PUBLIC :: metaff
CHARACTER(LEN=75),PRIVATE :: zaxislabel, customcolfile
CHARACTER(LEN=75),DIMENSION(ncurves),PRIVATE :: clegtxt
CHARACTER(LEN=75),DIMENSION(6),PRIVATE :: clylabel,crylabel !bar graph y-axis labels
CHARACTER(LEN=85),DIMENSION(6),PUBLIC :: pietittxt !pie chart title text array
CHARACTER(LEN=32),DIMENSION(2*NTCP),PUBLIC :: pielegtxt !pie chart data legend text array
CHARACTER(LEN=85),DIMENSION(ncurves,6),PRIVATE :: cpietittxt !pie chart title text array
CHARACTER(LEN=32),DIMENSION(ncurves,2*NTCP),PRIVATE :: cpielegtxt !pie chart data legend text array
CHARACTER(LEN=messagetextlength),DIMENSION(messagelinesmax),PRIVATE :: messagetxt
!...
LOGICAL(4),PUBLIC :: newplot,newplot2,squareplot,ternaryLLEoutplot,customCol, &
    & smallsymbs,bigsymbs,bigDsymbs,specialyrange,specialyrange2,specialxrange,logyaxis,messagebox, &
    & bargraphplot,fixminxax,directory_outputGouri,xorg3bar,legplot,logxaxis, &
    & plotbarborder,nogridlines
LOGICAL(4),PRIVATE :: zlogscalebar, showcustomcolbar
LOGICAL(4),DIMENSION(3),PUBLIC :: bargylimits !used to define upper yaxis limits for up to 3 different bar graphs on a plot page
LOGICAL(4),DIMENSION(ncurves),PUBLIC :: pielabshift
LOGICAL(4),DIMENSION(ncurves),PRIVATE :: clogyaxis,clegplot,cbigDsymbs
LOGICAL(4),DIMENSION(nrows),PUBLIC :: plotbarlabels
LOGICAL(4),DIMENSION(2*NTCP,nrows),PUBLIC :: cplotbarlabels

!=================================================================================================================

CONTAINS  !the subroutines within this module...

SUBROUTINE DefinePlotData(xvals,yvals,yerr,waxis,color,ltype,symbtype,lstyle,legtxt,thick)
!meaning of the curve defining parameters:
!DefinePlotData(xvals, yvals, yerr(the +- error bar size), waxis(y-axis 1 or 2?), color, ltype(1=only line, 2=symb and line, > 2 means symbols and no line),
!               symbtype(-1=no symbol), lstyle(1=solid, 2=dashm, 3=dash, 4=Mydot, 5=dot, 6=dash-dotted),legtxt,cthickn)

IMPLICIT NONE
!...
INTEGER(4),INTENT(IN) :: waxis,color,ltype,symbtype,lstyle
REAL(8),INTENT(IN) :: thick
INTEGER(4) :: i,N,npts,nsize
!...
REAL(8) :: yy,roundnicely
REAL(8),DIMENSION(:),INTENT(IN) :: xvals,yvals,yerr
REAL(8),DIMENSION(5) :: yrayspec
!...
CHARACTER(LEN=75),INTENT(IN) :: legtxt
!...
LOGICAL(4) :: isdata
!...
EXTERNAL :: roundnicely
!.........................................................................................

IF (newplot) THEN  !this start a new plot page, thus reset the counters and arrays:
    curvenr = 0
    curves_l = 0
    curves_r = 0
    cxvals = 0.0D0
    cyvals = 0.0D0
    cyerr = 0.0D0
    cwaxis = 0
    ccolor = 0
    cltype = 0
    rownr = 0
    cth = 1.0D0 !curvethickness
    csymbtype = 0
    clstyle = 0
    clegtxt = ""
    xmaxa = 0.0D0
    xmina = 0.0D0
    yminal = 0.0D0
    yminar = 0.0D0
    lyraymax = 0.0D0
    ryraymax = 0.0D0
    clogyaxis = .false.
    clegplot = .true.
    cpietittxt = ""
    cpielegtxt = ""
    cbigDsymbs = .false.
    CALL defplotcolorlists() !call to define color tables for plots
ENDIF

nsize = MIN(SIZE(xvals(:), DIM=1), SIZE(yvals(:), DIM=1))

!check data range:
npts = 0
DO i = 1,nsize-1
    IF (xvals(i) > (-1.0D8) .AND. ((xvals(i) > DEPS .OR. xvals(i) < (-DEPS)) .OR. (xvals(i+1) > DEPS .OR. xvals(i+1) < (-DEPS)) .OR. &
        &  (yvals(i) > DEPS .OR. yvals(i) < (-DEPS)) .OR. (yvals(i+1) > DEPS .OR. yvals(i+1) < (-DEPS)))) THEN
    npts = npts+1 !count points worth plotting...
    ELSE
        EXIT !leave the loop and continue...
    ENDIF
ENDDO
!check the last array entry:
i = nsize
IF (xvals(i) > -1.0D8 .AND. ((xvals(i) > DEPS .OR. xvals(i) < (-DEPS)) .OR. (yvals(i) > DEPS .OR. yvals(i) < (-DEPS)))) THEN
    npts = npts+1 !count points worth plotting...
ENDIF

isdata = .true.
IF (npts < 1) THEN
    !    PRINT*,"Warning: no points of the xvals are different from 0.0D0!"
    !    PRINT*,"this warning was generated in: SUBROUTINE DefinePlotData"
    isdata = .false.
ENDIF

IF (isdata) THEN
    !set the curve attributes in the array:
    curvenr = curvenr+1
    cxvals(curvenr,1:npts) = xvals(1:npts)
    cyvals(curvenr,1:npts) = yvals(1:npts)
    cyerr(curvenr,1:npts) = yerr(1:npts)
    cwaxis(curvenr) = waxis
    ccolor(curvenr) = color  !the color number (0 - 255) (default: in the RAIN colour spectrum of DISLIN)
    cltype(curvenr) = ltype
    rownr(curvenr) = npts
    csymbtype(curvenr) = symbtype
    clstyle(curvenr) = lstyle
    clegtxt(curvenr) = TRIM(legtxt)
    cpietittxt(curvenr,:) = pietittxt
    cpielegtxt(curvenr,:) = pielegtxt
    cth(curvenr) = thick
    xmaxa(curvenr) = MAXVAL(xvals(1:npts))
    xmina(curvenr) = MINVAL(xvals(1:npts))
    IF (waxis == 1) THEN
        curves_l = curves_l +1
        yminal(curves_l) = MINVAL(yvals(1:npts))
    ELSE
        curves_r = curves_r +1
        yminar(curves_r) = MINVAL(yvals(1:npts))
    ENDIF
    clogyaxis(curvenr) = logyaxis
    clegplot(curvenr) = legplot
    IF (legplot .AND. LEN_TRIM(legtxt) < 1) THEN
        clegplot(curvenr) = .false.
    ENDIF

    !now find out the appropriate y-axis scaling range for this curve:
    N=npts
    IF (ltype == 1) THEN !ltype is line
        yrayspec(5) = MAX(min_y, lowerylimit) !0.49D0
    ELSE
        yrayspec(5) = min_y*0.9D0
    ENDIF
    yrayspec(3) = MAXVAL(yvals(1:N))
    yrayspec(3) = MAX(yrayspec(3),yrayspec(5))*1.01D0
    yrayspec(4) = yrayspec(5)*1.01D0
    IF (yrayspec(3) > yrayspec(4) .AND. yrayspec(4) > 1.0D-12 .AND. ltype == 1 .AND. max_y < 10.0D0) THEN  !Check for very large values and reduce to a given level...
        IF (yrayspec(4) < max_y) THEN
            yrayspec(3) = yrayspec(4)
        ELSE
            yrayspec(3) = max_y
        ENDIF
    ENDIF
    IF (yrayspec(3) > max_y) THEN
        yrayspec(3) = max_y  !the maximal value on the y-axis...
    ELSE IF (yrayspec(3) < min_y) THEN
        yrayspec(3) = min_y
    ENDIF

    !round up the values to a certain axis step size:
    yy = yrayspec(3)
    IF (yy > 1.0D0 .AND. yy < 1.1D0) THEN
        yy = roundnicely(yy,5,"down") ! function call rounding on the 5% level
    ELSE
        yy = roundnicely(yy,5,"up") ! function call
    ENDIF

    !assign the max value to the corresponding axis:
    IF (waxis == 1) THEN !left axis
        lyraymax(curves_l) = yy !left y-axis
    ELSE !right axis
        ryraymax(curves_r) = yy !right y-axis
    ENDIF

ENDIF !isdata

newplot = .false. !set to read in next curve for the same plot page
logyaxis = .false.
legplot = .true.
bigDsymbs = .false.
pietittxt = ""
pielegtxt = ""

END SUBROUTINE DefinePlotData
    
!=================================================================================================================    

SUBROUTINE defbardata(xvals,yvals,lylabel,rylabel,barpatcol,shadingpat,barbordercol)

IMPLICIT NONE
!interface vars:
REAL(8),DIMENSION(nrows),INTENT(IN) :: xvals
REAL(8),DIMENSION(2*NTCP,nrows),INTENT(IN) :: yvals  !the y-values of the different components (max. NTCP components)
CHARACTER(LEN=75),INTENT(IN) :: lylabel, rylabel !bar graph y-axis labels
INTEGER(4),DIMENSION(2*NTCP,nrows),INTENT(IN) :: barpatcol, shadingpat, barbordercol
!local vars:
INTEGER(4) :: npts, i, nicomp
!.........................................................................................

IF (newplot2) THEN  !this start a new plot page, thus reset the counters and arrays:
    barpmodeno = 0
    compbx = 0.0D0
    compby = 0.0D0
    compbpoints = 0.0D0
    xmaxb = 0.0D0
    xminb = 0.0D0
    cpietittxt = ""
    cpielegtxt = ""
    cplotbarlabels = .false.
    logxaxis = .false.
    CALL defplotcolorlists() !call to define color tables for plots
ENDIF

!check data range:
nicomp = 0
DO i = 1,2*NTCP
    IF (SUM(ABS(yvals(i,:))) > 1.0D-14) THEN
        nicomp = i !nicomp+1 !count number of different independent components (nicomp)
    ENDIF
ENDDO
npts = 0
DO i = 1,nrows
    IF (xvals(i) > (-1.0D6) .AND. ((xvals(i) > -1.0D-14 .OR. xvals(i) < (-1.0D-12)) .AND. SUM(ABS(yvals(1:nicomp,i))) > 1.0D-14) &
        & .OR. (xvals(i) > 1.0D0)) THEN
        npts = i !npts+1 !count number of bars worth plotting...
    ENDIF
ENDDO

!set the bar attributes in the arrays:
IF (npts > 0) THEN
    barpmodeno = barpmodeno+1 !counts the different number of bar plots on a plot page (the mode number)
    compbx(barpmodeno,1:npts) = xvals(1:npts)
    compby(barpmodeno,1:nicomp,1:npts) = yvals(1:nicomp,1:npts)
    compbpoints(barpmodeno) = npts
    clylabel(barpmodeno) = TRIM(lylabel)
    crylabel(barpmodeno) = TRIM(rylabel)
    ncomp(barpmodeno) = nicomp
    xmaxb(barpmodeno) = MAXVAL(xvals(1:npts))
    xminb(barpmodeno) = MINVAL(xvals(1:npts))
    cpietittxt(barpmodeno,1:6) = pietittxt(1:6)
    cpielegtxt(barpmodeno,1:NTCP) = pielegtxt(1:NTCP)
    cplotbarlabels(barpmodeno,1:npts) = plotbarlabels(1:npts)
    cbarpatterncol(barpmodeno,1:nicomp,1:npts) = barpatcol(1:nicomp,1:npts)
    cbarshadpat(barpmodeno,1:nicomp,1:npts) = shadingpat(1:nicomp,1:npts)
    cbarbordcol(barpmodeno,1:nicomp,1:npts) = barbordercol(1:nicomp,1:npts)
ENDIF

curvepointsno = npts
newplot2 = .false. !set to read in next curve for the same plot page
pietittxt = ""
pielegtxt = "" 
plotbarlabels = .false.
legplot = .true.

END SUBROUTINE defbardata

!=================================================================================================================    

SUBROUTINE SetTextlines(text)

IMPLICIT NONE

CHARACTER(*),INTENT(IN) :: text
!...................................................................

IF (.NOT. messagebox) THEN !reset the text strings for the message box in the plotpage (if any)
    messagetxt = "none"
    txtl = 0
ENDIF

IF (LEN_TRIM(text) > messagetextlength) THEN
    !$OMP CRITICAL
    WRITE(*,*) "The textstring submitted to SetTextlines is too long!"
    WRITE(*,*) "The length is limited to: ",messagetextlength,"characters!"
    WRITE(*,*) ""
    !$OMP END CRITICAL
ELSE
    txtl = txtl+1
    messagetxt(txtl) = TRIM(text)  !the text string to read in on line txtl for later output on the plotpage in a message box.
    messagebox = .true.
ENDIF
    
END SUBROUTINE SetTextlines
    
!=================================================================================================================

SUBROUTINE PlotNow(ptitle,xlabel,lylabel,rylabel,outname_cplot,StandardDev)

USE dislin  !the DISLIN plotting library (www.dislin.de)

IMPLICIT NONE
!...
INTEGER(4) :: i,N,NXL,NYL,NZL,symbsize,NLIN,NMAXLN,ILIN,nsteps, &
    & NXA,NYA,Tsize,Lsize,Namesize,Numbsize,Lcompsize,ifault,NTYP,ileg,legstop, &
    & bigsymbsize,iclr,ityp, zdigits, nxwid, nxpos, nypos, Nheight
INTEGER(4),DIMENSION(2) :: NRAY
INTEGER(4),DIMENSION(4) :: INRAY
!...
PARAMETER (NMAXLN = 75)
!...
REAL(8) :: StandardDev,xmaxax,xminax,xlab,xstep,yminax,ymaxax,ylab,ystep,roundnicely,logymin, &
    & xminaxin,xmaxaxin,dtiny,xte,xte2, rinc, rypos
REAL(8),DIMENSION(3) :: xrayspec, yrayspec_l, yrayspec_r
REAL(8),DIMENSION(nrows) :: xax,yax,E1RAY,E2RAY
!...
CHARACTER(LEN=75),INTENT(IN) :: outname_cplot
CHARACTER(LEN=75) :: outname_dislin
CHARACTER(LEN=(COUNT(clegplot(1:curvenr))*NMAXLN)) :: CBUF
CHARACTER(LEN=NMAXLN) :: CSTR
CHARACTER(LEN=132),DIMENSION(4) :: ptitle
CHARACTER(LEN=75) :: xlabel,lylabel,rylabel
CHARACTER(LEN=messagetextlength) :: hline
!...
LOGICAL(4) :: ryaxis, logy_l, logy_r, uselegpat
!...
EXTERNAL :: roundnicely

!.........................................................................................
!initialize some variables:
xrayspec = 0.0D0
yrayspec_l = 0.0D0
yrayspec_r = 0.0D0
xax = 0.0D0
yax = 0.0D0
dtiny = EPSILON(1.0D0)*2.0D2
xminax = 0.0D0; xmaxax = 1.0D0; xlab = 0.0D0; xstep = 0.5D0;
yminax = 0.0D0; ymaxax = 1.0D0; ylab = 0.0D0; ystep = 0.1D0;

!check whether there is any curve/points to plot at all:
IF (curvenr < 1) THEN !there is nothing to plot
    !$OMP CRITICAL
    WRITE(*,*),"There was no curve to plot in this call to PlotNow!"
    WRITE(*,*),"proceeding without producing a plot..."
    WRITE(*,*),""
    !$OMP END CRITICAL
    newplot = .true. !set it back for next plot page
    newplot2 = .true.
    messagebox = .false.
    smallsymbs = .false.
    bigsymbs = .false.
    bigDsymbs = .false.
    specialxrange = .false.
    specialyrange = .false.
    specialyrange2 = .false.
    customCol = .false.
    logyaxis = .false. !set the logarithmic y-axis to false (so use linear).
    logxaxis = .false.
    legplot = .true.
    nogridlines = .false.
    showcustomcolbar = .false.
    customcolfile = 'RGBColTableRealWinXP8BitAZ.dat' !default value
    ptitle = ""
    RETURN !leave the subroutine
ENDIF

!paramparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparam
!paramparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparam
!:: adjustable  plot parameters ::

!squareplot scales the x, y axis length to plot in a square shaped box

!Set font sizes etc.:

IF (squareplot) THEN
    symbsize = 32 !40 !26 !symbol size
    Tsize = 22   !Title font size
    IF (ANY(LEN_TRIM(ptitle(1:4)) > 70)) THEN
        Tsize = 16
    ENDIF
    Lsize = 19      !Legends font size
    Lcompsize = 24  !composition plot legend size
    Namesize = 36   !axis name font
    Numbsize = 32   !axis numbers size
ELSE
    symbsize = 36   !symbol size
    Tsize = 22      !Title font size
    IF (ANY(LEN_TRIM(ptitle(1:4)) > 70)) THEN
        Tsize = 16
    ENDIF
    Lsize = 19      !Legends font size
    Lcompsize = 24  !composition plot legend size
    Namesize = 40   !axis name font
    Numbsize = 32   !axis numbers size
ENDIF
IF (smallsymbs) THEN
    symbsize = symbsize*3/4 !18  !symbol size
    bigsymbsize = symbsize*5/4  !32
ELSE IF (bigsymbs) THEN
    bigsymbsize = 44
    symbsize = bigsymbsize
ENDIF
logymin = 1.0D-3  !the minimum y value when logarithmic y-axis are plottet without special y-range given

!paramparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparam
!paramparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparamparam

!set output file
CALL METAFL(metaff)
outname_dislin=TRIM(outname_cplot)//"."//(TRIM(metaff))
outname_dislin="../Output_Plots/"//TRIM(outname_dislin)
CALL SETFIL(TRIM(outname_dislin)) !set filename in dislin
CALL FILMOD("DELETE") !files will be overwritten if they have the same name..
CALL SCRMOD ("NOREV")  !set background color to white and foreground to black
CALL SETPAG("DA4L")
CALL SCLFAC(1.0D0)
CALL SCLMOD("DOWN")
CALL IMGFMT ("RGB")
CALL WINSIZ (1600, 1200)
IF (metaff == "png" .OR. metaff == "gif") THEN
    CALL IMGFMT ("RGB")
    CALL WINSIZ (1200, 900)
ENDIF

!Set the overall plot parameters and initialize the plot page:
CALL DISINI !initialize dislin
CALL SETVLT("RAIN")  !Rainbow colors (256)
CALL PAGFLL(255) !background color white
CALL SETCLR(0) !black
IF (metaff == "wmf") THEN
    CALL WINFNT("Arial")  !FONT in Windows  !Arial  ! Times New Roman
ELSE IF (metaff == "eps" .OR. metaff == "ps" .OR. metaff == "pdf" .OR. metaff == "svg") THEN
    CALL PSFONT("Helvetica")  !for ps, eps and pdf output
    CALL PSMODE('BOTH')  !allow Greek and Italic modes
ELSE
    CALL COMPLX
ENDIF
CALL EUSHFT ('GERMAN', '!') !used to print German umlaut characters with Dislin TeX: then !o writes ö
!date of print and dislin version on the lower right corner
CALL HEIGHT(Lsize)
CALL PAGHDR("","",2,0)
CALL ERRMOD ("CHECK", "OFF")  !Don't pass warnings of points lying outside of the coordinate systems
CALL CHASPC(-0.06D0)
CALL HEIGHT(Namesize)
!write some LaTeX code for better axes-text:
CALL TEXMOD ("ON")

!Set Titeltext on line no. 1:
CALL TITLIN (TRIM(ptitle(1)), 1)
!Set Titeltext on line no. 2:
CALL TITLIN (TRIM(ptitle(2)), 2)
!Set Titeltext on line no. 3:
CALL TITLIN (TRIM(ptitle(3)), 3)
!Set Titeltext on line no. 4:
CALL TITLIN (TRIM(ptitle(4)), 4)
CALL LINESP(2.6D0)
CALL TITJUS ("LEFT")

!Set axis labels:
CALL GETPAG (NXL, NYL)          !get the page total DIMENSIONs NXL=2970, NYL = 2100
IF (squareplot) THEN
    CALL AXSPOS (360, 1400)     !position the axis system (the lower left axis corner)
    CALL AXSLEN (1000, 1000)    !now set axis length on the page in plot coord.
ELSE
    CALL AXSPOS (360, 1400)     !position the axis system (the lower left axis corner)
    CALL AXSLEN (1500, 1000)    !now set axis length on the page in plot coord.
    !!        CALL AXSPOS (INT(DFLOAT(NXL)*0.10D0), INT(DFLOAT(NYL)*0.86D0))  !position the axis system (the lower left axis corner)
    !!        CALL AXSLEN (INT(DFLOAT(NXL)*0.70D0), INT(DFLOAT(NYL)*0.64D0))  !now set axis length on the page in plot coord.
ENDIF
CALL NAMDIS (44, "Y")  !distance between axis text and axis names
CALL NAMDIS (44, "X")
CALL HNAME (Namesize)
CALL NAME (TRIM(xlabel), "X")
CALL NAME (TRIM(lylabel), "Y")

!find out if there is a second y-axis...
IF (ANY(cwaxis(1:curvenr) == 2)) THEN
    ryaxis=.true.
ELSE
    ryaxis=.false.
ENDIF

!1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st
!1st                                   start of 1st y-axis                                    1st
!1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st

!check for logarithmic y-axis scaling:
logy_l = .false.
logy_r = .false.
DO i=1,curvenr
    IF (clogyaxis(i) .AND. cwaxis(i) == 1) THEN !log left y-axis
        logy_l = .true.
        CALL NEGLOG(DEPS)
    ELSE IF (clogyaxis(i) .AND. cwaxis(i) == 2) THEN !log right y-axis
        logy_r = .true.
        CALL NEGLOG(DEPS)
    ENDIF
ENDDO

!the maximum/minimum x-axis/y-axis values and its scaling:
IF (.NOT. specialxrange) THEN
    xmaxaxin = MAXVAL(xmaxa(1:curvenr))
    xminaxin = MINVAL(xmina(1:curvenr))
    IF (fixminxax) THEN
        xminaxin = MIN(0.0D0, xminaxin)
    ENDIF
ELSE
    xminaxin = lowerxlimit
    xmaxaxin = upperxlimit
ENDIF

IF (specialxrange) THEN
    xminax = lowerxlimit
    xmaxax = upperxlimit
    xmaxaxin = upperxlimit
ELSE
    xminax = xminaxin
    xmaxax = xmaxaxin
ENDIF
!automatic scaling of x-axis according to precalculated axis bounds:
xrayspec(1) = xminax
xrayspec(2) = xmaxax
xrayspec(3) = xmaxaxin
IF (logxaxis) THEN
    xrayspec(2) = ABS(xrayspec(2))
    IF (xrayspec(1) <= 0.0D0) THEN
        xrayspec(1) = xrayspec(2)*1.0D-8
    ENDIF
    CALL AXSSCL ("LOG", "X")
    CALL LOGTIC ("FULL") !"AUTO"
ELSE
    CALL AXSSCL ("LIN", "X")
ENDIF
CALL SETSCL(xrayspec, 3, "X")

ymaxax = MAXVAL(lyraymax(1:curves_l))
yrayspec_l(3) = ymaxax
yminax = MINVAL(yminal(1:curves_l))
yminax = roundnicely(yminax,5,"down")
IF (logy_l .AND. logy_r) THEN
    yminax = MIN(MINVAL(yminar(1:curves_r)),MINVAL(yminal(1:curves_l)))
    yminax = roundnicely(yminax,5,"down")
    yrayspec_l(3) = MAX(MAXVAL(ryraymax(1:curves_r)),yrayspec_l(3))
    yrayspec_l(3) = MIN(10.0D0**(LOG10(MAX(yminax,logymin))+8.0D0), yrayspec_l(3))  !limit the range on the y-axis to 8 orders of magnitude
ELSE IF (logy_l) THEN
    yrayspec_l(3) = MIN(10.0D0**(LOG10(MAX(yminax,logymin))+8.0D0), yrayspec_l(3))  !limit the range on the y-axis to 8 orders of magnitude
ENDIF
IF (specialyrange) THEN
    yminax = lowerylimit
ENDIF

!first find out the y-axis ranges which cover the whole number of curves:
!for the left y-axis:
yrayspec_l(1) = 0.0D0 !the lowest value on the y-axis
IF (yminax < 0.0D0) THEN
    yrayspec_l(1) = yminax
ENDIF
yrayspec_l(2) = min_y
NLIN = 0
CALL AXSSCL ("LIN", "Y") !the default choice
IF (logy_l) THEN
    yrayspec_l(2) = MAX(yrayspec_l(2), logymin)
    CALL AXSSCL ("LOG", "Y")
    CALL LOGTIC ("FULL") !"AUTO"
    IF (yminax > -1.0D-12 .AND. yminax < 1.0D-3 .AND. yrayspec_l(3) > 1.0D-3) THEN
        yrayspec_l(1) = logymin
    ELSE
        IF (LOG10(yminax) > 0.0D0) THEN !round up
            yrayspec_l(1) = 10.0D0**(roundnicely(LOG10(yminax),5,"up"))
        ELSE
            yrayspec_l(1) = 10.0D0**(roundnicely(LOG10(yminax),5,"down"))
        ENDIF
    ENDIF
ENDIF
DO i=1,curvenr
    IF (cwaxis(i) == 1) THEN
        CSTR = TRIM(clegtxt(i))
        IF (LEN_TRIM(CSTR) > 2) THEN
            NLIN = NLIN+1 !legend line Nr.
        ENDIF
        IF (cltype(i) == 1) THEN !"line"
            csymbtype(i)= -1
        ENDIF
    ENDIF
ENDDO

IF (specialyrange) THEN
    yrayspec_l(1) = lowerylimit
    yrayspec_l(2) = upperylimit
    ymaxax = MIN(ymaxax, upperylimit)
    yminax = MIN(yminax, lowerylimit)
    IF (logy_l) THEN
        ystep = 1.0D0
        yrayspec_l(3) = lowerylimit*10.001D0 !cover at least one order of magnitude in log10 scale
        yrayspec_l(2) = roundnicely(upperylimit,80,"up")
        yrayspec_l(2) = MAX(yrayspec_l(2), 1.0D-15)
    ELSE
        yrayspec_l(3) = upperylimit
    ENDIF
ENDIF
!automatic scaling of axis according to a selected part of the data:
CALL SETSCL(yrayspec_l, 3, "Y")

!plot axis and ticks from (xmin,xmax,first_x_label,increment_x / ymin,ymax,first_y_label,increment_Y):
IF (ryaxis) THEN
    CALL SETGRF ("NAME", "NAME", "TICKS", "LINE")
ELSE
    CALL SETGRF ("NAME", "NAME", "TICKS", "TICKS")
ENDIF
CALL HEIGHT(Numbsize)
CALL CHASPC(-0.07D0)
CALL TICKS(4, "XY") !the number of ticks between axis labels (default: 2)
IF (ABS(xstep) < 0.01D0) THEN
    CALL LABDIG (4, "X") !number of decimals xaxis
ELSE IF (ABS(xstep) < 0.1D0) THEN
    CALL LABDIG (3, "X") !number of decimals xaxis
ELSE IF (ABS(xstep) < 0.2D0) THEN
    CALL LABDIG (2, "X") !number of decimals xaxis
ELSE IF (ABS(xstep) > 1.0D0) THEN
    CALL LABDIG (1, "X") !number of decimals xaxis
    IF (xstep < 10.0D0) THEN
        CALL TICKS(INT(xstep), "X")
    ENDIF
ELSE
    CALL LABDIG (1, "X") !number of decimals xaxis
ENDIF
IF ((.NOT. logy_l) .AND. (ABS(ymaxax-yminax) < 1.0D-3 .OR. ABS(ymaxax-yminax) > 9.99999D3)) THEN
    CALL LABELS ("FEXP", "Y") !EXP will use powers of 10 "scientific" number format for axis labels; FEXP will use "1.0E-8" style.
ELSE
    CALL LABELS ("FLOAT", "Y")
ENDIF
!...
!some of the axis use automatic axis scaling and overwrite whatever values are given in the code below:
CALL SETVLT("RAIN")     !the Dislin rainbow palette
CALL COLOR("BLACK")
CALL SOLID
CALL PENWID(1.0D0)      !set to default value
CALL LINWID(1)          !set to default value
CALL THKCRV(1)          !set to default value
xlab = xminax; ylab = yminax;
IF (xmaxax > 1.0D0 .AND. barpmodeno > 1 .AND. xstep > 0.99D0) THEN !LLE relative diff. and other axis with integer numbering
    CALL INTAX
    yminax = 0.0D0; ymaxax = 1.0D0; ylab = 0.0D0; ystep = 0.1D0;
    CALL GRAF(xminax,xmaxax,xlab,xstep, yminax,ymaxax,ylab,ystep)
ELSE IF (xmaxax > 1.0D0 .OR. xminax > 0.0D0 .OR. (xmaxax-xminax < 0.4D0)) THEN !x-axis from xminax to xmaxax
    CALL GRAF(xminax,xmaxax,xlab,xstep, yminax,ymaxax,ylab,ystep) !will actually be scaled automatically using SETSCL
ELSE !xaxis in the range 0.0D0 to 1.0D0
    CALL GRAF(xminax,xmaxax,xlab,xstep, yminax,ymaxax,ylab,ystep) !will actually be scaled automatically using SETSCL
ENDIF

!Plot now the title:
CALL HEIGHT(Tsize)  !Letter height of title
CALL CHASPC(-0.04D0)
CALL TITLE
!plot Gridlines:
IF (.NOT. nogridlines) THEN
    CALL SETVLT("GREY")
    IF (logy_l) THEN
        CALL SOLID
        CALL SETCLR(230)
        CALL PENWID(0.1D0) !the linewidth (allowing values smaller than 1)
        CALL GRID(0,10) !light grey fine gridlines
        CALL SETCLR(160)
        NRAY(1) = 4  !Stift unten
        NRAY(2) = 6  !Stift oben
        CALL MYLINE (NRAY, 2)
        CALL PENWID(0.2D0) !the linewidth (allowing values smaller than 1)
        CALL GRID(1,1)
        CALL LINWID(1) !reset the linewidth
    ELSE
        IF (logxaxis) THEN
            CALL SOLID
            CALL SETCLR(230)
            CALL PENWID(0.1D0) !the linewidth (allowing values smaller than 1)
            CALL GRID(10,0) !light grey fine gridlines on x-axis
            CALL SETCLR(160)
            NRAY(1) = 4  !Stift unten
            NRAY(2) = 6  !Stift oben
            CALL MYLINE (NRAY, 2)
            CALL PENWID(0.2D0) !the linewidth (allowing values smaller than 1)
            CALL GRID(1,1)
            CALL LINWID(1) !reset the linewidth
        ELSE
            CALL SETCLR(160)
            NRAY(1) = 4  !Stift unten
            NRAY(2) = 6  !Stift oben
            CALL MYLINE (NRAY, 2)
            CALL PENWID(0.2D0) !the linewidth (allowing values smaller than 1)
            CALL GRID(1,1)
            CALL LINWID(1) !reset the linewidth
        ENDIF
    ENDIF
ENDIF
CALL SOLID

!calculate and initialize number of lines of 1st legend:
CBUF = ""
CSTR = ""
NLIN = COUNT(clegplot(1:curvenr) .AND. cwaxis(1:curvenr) == 1)
IF (NLIN > 30) THEN !dislin currently can only store 30 curve (attributes) for legends
    CALL LEGINI(CBUF, 30, NMAXLN)
ELSE
    CALL LEGINI(CBUF, NLIN, NMAXLN)
ENDIF
!check whether a special legend plotting mode is necessary to store curve attributes and change curve color of white curves for visible contrast in legend:
IF (ANY(clegplot(1:curvenr) .AND. cwaxis(1:curvenr) == 1 .AND. ccolor(1:curvenr) == 255)) THEN
    uselegpat = .true.
ELSE
    uselegpat = .false.
ENDIF

!plot the curves on the 1st axis ###
ILIN = 0
IF (ANY(.NOT. clegplot(1:curvenr))) THEN
    legstop = 2
ELSE
    legstop = 1
ENDIF
DO ileg = 1,legstop !loop two times over the curves, the first time to only plot curves with legend entries.
    DO i = 1,curvenr
        IF ((ileg == 1 .AND. clegplot(i)) .OR. (ileg == 2 .AND. (.NOT. clegplot(i)))) THEN
            IF (cwaxis(i) == 1) THEN !plot it on the left y-axis, the first axis system
                !plot the defined curves:
                IF (customCol) THEN !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
                    CALL VLTFIL(TRIM(customcolfile), "LOAD") !the Win XP Icon Color Palette plus Grayscale etc. as defined by Andi Zuend
                ELSE
                    CALL SETVLT("RAIN") !the Dislin rainbow palette
                ENDIF
                CALL SETCLR(ccolor(i))                  !curve colour
                CALL THKCRV(1)                          !curve thickness
                CALL LINWID(1)                          !the linewidth (will be adjusted by PENWID instead)
                CALL PENWID(1.0D0*cth(i))               !the linewidth (allowing also values smaller than 1)
                CALL INCMRK(0)                          !line
                IF (cltype(i) == 2) THEN                !symbols and line
                    CALL INCMRK(1)                      !line connecting data points every n=1 datapoints will be plotted
                    CALL MARKER(csymbtype(i))           !plot symbol type at point coordinates
                    IF (smallsymbs) THEN
                        IF (cBigDsymbs(i)) THEN
                            CALL HSYMBL(bigsymbsize)    !symbol size
                            CALL PENWID(1.3D0*cth(i))   !the linewidth (allowing also values smaller than 1)
                        ELSE
                            CALL PENWID(0.5D0*cth(i))   !the linewidth (allowing also values smaller than 1)
                            CALL HSYMBL((symbsize*3)/4) !symbol size
                        ENDIF
                    ELSE
                        CALL PENWID(1.0D0*cth(i))       !the linewidth (allowing also values smaller than 1)
                        IF (cBigDsymbs(i)) THEN
                            CALL HSYMBL(bigsymbsize)    !symbol size
                        ELSE
                            CALL HSYMBL(symbsize)       !symbol size
                        ENDIF
                    ENDIF
                ELSE IF (cltype(i) > 2) THEN            !symbols without line
                    CALL INCMRK(-1)                     !no line connecting data points, every n=1 datapoint will be plotted
                    CALL MARKER(csymbtype(i))           !plot symbol type at point coordinates
                    IF (smallsymbs) THEN
                        IF (cBigDsymbs(i)) THEN
                            CALL HSYMBL(bigsymbsize)    !symbol size
                            CALL PENWID(1.3D0*cth(i))   !the linewidth (allowing also values smaller than 1)
                        ELSE
                            CALL PENWID(0.5D0*cth(i))   !the linewidth (allowing also values smaller than 1)
                            CALL HSYMBL(symbsize)       !symbol size
                        ENDIF
                    ELSE
                        IF (cBigDsymbs(i)) THEN
                            CALL PENWID(1.3D0*cth(i))   !the linewidth (allowing also values smaller than 1)
                            CALL HSYMBL(bigsymbsize)    !symbol size
                        ELSE
                            CALL PENWID(1.0D0*cth(i))   !the linewidth (allowing also values smaller than 1)
                            CALL HSYMBL(symbsize)       !symbol size
                        ENDIF
                    ENDIF
                ENDIF
                SELECT CASE (clstyle(i)) !define the line style
                CASE(1)
                    CALL LNCAP ("LONG") ! line caps
                    CALL SOLID
                    NTYP = 0
                CASE(2)
                    CALL LNCAP ("CUT") !cut line caps
                    CALL DASHM
                    NTYP = 5
                CASE(3)
                    CALL LNCAP ("ROUND") !rounded line caps
                    CALL DASH
                    NTYP = 2
                CASE(4) !MyDot style
                    NRAY(1) = 1  !Pen down
                    NRAY(2) = INT(1.5D0*MAX(cth(i), 4.0D0)) !2*MAX(INT(cth(i)), 3) !+4  !Pen up
                    CALL LNCAP ("ROUND") !rounded line caps
                    CALL MYLINE (NRAY, 2) !self defined line style "DOT" with size according to line thickness
                    NTYP = 6
                CASE(5)
                    CALL LNCAP ("ROUND") !("ROUND") !rounded line caps
                    CALL DOT
                    NTYP = 1
                CASE(6) !. - . - . my dash-dotted line style
                    INRAY(1) = 1  !Pen down
                    INRAY(2) = 3*MAX(3, INT(cth(i))) !Pen up
                    INRAY(3) = 4*MAX(3, INT(cth(i))) !Pen down
                    INRAY(4) = 3*MAX(3, INT(cth(i))) !Pen up
                    CALL LNCAP ("ROUND") !rounded line caps
                    CALL MYLINE (INRAY(1:4), 4)
                    NTYP = 6
                END SELECT
                CALL LNJOIN ("SHARP") !'SHARP' or 'TRUNC'
                CALL NOCHEK  !suppress warning of points lying outside of the plotting area
                !now plot the selected curve:
                N=rownr(i)
                xax = 0.0D0
                yax = 0.0D0
                xax(1:N) = cxvals(i,1:N)
                yax(1:N) = cyvals(i,1:N)
                CALL CURVE(xax(1:N),yax(1:N),N)
                IF (ileg == 1) THEN !set the corresponding legend text and color:
                    CSTR = TRIM(clegtxt(i))
                    IF (LEN_TRIM(CSTR) > 1) THEN
                        IF (ILIN < 30) THEN !only max. 30 curves attributes can be printed for the legend (current dislin limit)
                            ILIN = ILIN+1
                            !!change color in legends when the line color is white (to have some contrast to the white  background color):
                            IF (uselegpat) THEN
                                IF (ccolor(i) == 255) THEN !white
                                    iclr = 0 !black
                                ELSE
                                    iclr = ccolor(i)
                                ENDIF
                                ityp = NTYP !select a line style from the default styles, closest to the one of the custom types.
                                CALL LEGPAT(ityp, INT(cth(i)), csymbtype(i), iclr, -1, ILIN) !ILIN
                            ENDIF
                            CALL LEGVAL (1.0D0, 'SYMBOL')       !symbol size in legend (default = 0.8D0)
                            CALL LEGLIN(CBUF, TRIM(CSTR), ILIN) !define the legend line text
                        ENDIF
                    ENDIF
                ENDIF
                IF (ANY(cyerr(i,1:N) > dtiny)) THEN
                    E1RAY(1:N) = cyerr(i,1:N) !error bars in y-axis
                    E2RAY(1:N) = cyerr(i,1:N)
                    CALL LINWID(1)
                    IF (cltype(i) < 2) THEN !no symbols, just line
                        CALL PENWID(MIN(0.25D0, MAX(0.1D0, cth(i)/5.0D0)))
                        CALL HSYMBL(MAX(1,INT(cth(i)/3.0D0))) !set symbol size to tiny size
                        CALL MARKER(-1) !suppress plotting symbol when only error bars are required here... CALL MARKER(3)
                    ELSE
                        CALL PENWID(MIN(1.0D0, MAX(0.2D0, cth(i)/3.0D0)))
                        CALL MARKER(-1) !CALL MARKER(csymbtype(i))
                    ENDIF
                ENDIF
                CALL LINWID(1) !reset the linewidth default value
                CALL PENWID(1.0D0)
            ENDIF !axisc 1 IF
        ENDIF !ileg IF
    ENDDO !i, plot curves 1st axis
ENDDO !ileg

!plot legend to the first system (left y-axis):
IF (LEN_TRIM(CBUF) > 1) THEN
    CALL SETVLT("RAIN")
    CALL SETCLR(0) !black
    CALL SOLID
    CALL HEIGHT (Lsize)  !Legend font size
    CALL LEGTIT ("left y-axis:                ")
    !    CALL LEGTIT("")
    CALL LINESP(2.6D0)
    IF (showcustomcolbar) THEN
        nxpos = 160
    ELSE
        nxpos = 0
    ENDIF
    IF (squareplot) THEN
        CALL LEGPOS (1900+nxpos, 40)
    ELSE
        IF (curvenr > 1 .AND. ryaxis) THEN
            CALL LEGPOS (2100+nxpos, 40)
        ELSE
            CALL LEGPOS (2100+nxpos, 40)
        ENDIF
    ENDIF
    CALL FRAME(0)  !no box around the legend
    IF (customCol) THEN !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
        CALL VLTFIL(TRIM(customcolfile), "LOAD") !the Win XP Icon Color Palette plus Grayscale etc. as defined by Andi Zuend
    ELSE
        CALL SETVLT("RAIN") !the Dislin rainbow palette
    ENDIF
    CALL SETCLR(0) !black
    CALL LINWID(MIN(3, MAX(1, INT(MAXVAL(cth(1:curvenr)))))) !set here curve thickness of curve / points for all legend entries (as no curve-specific choice possible for now).
    CALL LEGEND (CBUF, 3)
    CALL LINWID(1) !reset
ENDIF

IF (showcustomcolbar) THEN  !plot a colour bar (z-axis) that has been defined by the used colours
    !plot a custom-made z-axis colour bar based on a collection of stacked coloured rectangles.
    !The colours were determined by the customized colour scaling used in showzcolbar(..) via defzcolval(..).
    CALL GETPOS (NXA, NYA)
    CALL GETLEN (NXL, NYL, NZL)
    nxpos = NXA+NXL+100                         !middle of rectangle which is nxwid coordinate points wide in total
    nxwid = 50
    rinc = REAL(NYL,KIND=8)/253.0D0             !real-valued increment in plot coordinates
    Nheight = CEILING(rinc)                     !rectangle height
    rypos = REAL(NYA, KIND=8) -0.5D0*rinc       !y-coordinate of first point
    DO i = 1,254
        nypos = FLOOR(rypos)
        CALL POINT(nxpos, nypos, nxwid, Nheight, customzcolno(i))
        rypos = MAX(rypos -rinc, NYA-NYL+0.5D0*rinc)
    ENDDO
    !plot a secondary y-axis as the actual z-colour axis:
    CALL SETCLR(0) !black
    CALL SOLID
    CALL PENWID(1.0D0)
    CALL LINWID(1)
    CALL TICKS (2, "Y")
    CALL HEIGHT(CEILING(Numbsize*0.85))
    CALL RVYNAM
    CALL LABJUS('LEFT', 'Y')
    CALL CHASPC(-0.07D0)
    CALL LABDIS(10, 'Y')
    CALL NAMDIS(10, 'Y')
    CALL HNAME (Namesize)
    CALL LABDIG(1, 'Y')    !number of decimals digits
    CALL YAXIS(zcolbar0, zcolbar254, zcolbar0, zcolbarstep, NYL, TRIM(zaxislabel), 0, nxpos+nxwid/2, NYA)
ENDIF

!1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st
!1st                                     end of 1st y-axis                                    1st
!1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st1st


!2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd
!2nd                                    start of 2nd y-axis                                   2nd
!2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd
IF (ryaxis .AND. (.NOT. bargraphplot)) THEN ! 2ndaxis

    !set back level to 1 for plotting second axis system
    CALL ENDGRF
    CALL SETVLT("GREY")
    CALL SETCLR(0)
    CALL SOLID

    !first find out the y2-axis ranges wich cover the whole number of curves:
    ymaxax = MAXVAL(ryraymax(1:curves_r))
    yrayspec_r(3) = ymaxax
    yminax = MINVAL(yminar(1:curves_r))
    yminax = roundnicely(yminax,5,"down")
    IF (logy_l .AND. logy_r) THEN
        yminax = MIN(MINVAL(yminar(1:curves_r)),MINVAL(yminal(1:curves_l)))
        yminax = roundnicely(yminax,5,"down")
        yrayspec_r(3) = MAX(MAXVAL(ryraymax(1:curves_r)),yrayspec_r(3))
        yrayspec_r(3) = MIN(10.0D0**(LOG10(MAX(yminax,logymin))+8.0D0), yrayspec_r(3))  !limit the range on the y-axis to 8 orders of magnitude
    ELSE IF (logy_r) THEN
        yrayspec_r(3) = MIN(10.0D0**(LOG10(MAX(yminax,logymin))+8.0D0), yrayspec_r(3))  !limit the range on the y-axis to 8 orders of magnitude
    ENDIF

    !for the right y-axis:
    yrayspec_r(1) = 0.0D0 !the lowest value on the y-axis
    IF (yminax < 0.0D0) THEN  !xmaxax > 1.1 .AND.
        yrayspec_r(1) = yminax !-10.0
    ENDIF
    yrayspec_r(2) = 0.01D0
    CALL AXSSCL ("LIN", "Y") !the default choice
    IF (logy_r) THEN
        yrayspec_r(2) = MAX(yrayspec_r(2), logymin)
        CALL AXSSCL ("LOG", "Y")
        CALL LOGTIC ("FULL") !"AUTO"
        IF (yminax > -1.0D-12 .AND. yminax < 1.0D-3 .AND. yrayspec_r(3) > 1.0D-3) THEN
            yrayspec_r(1) = logymin
        ELSE
            IF (LOG10(yminax) > 0.0D0) THEN !round up
                yrayspec_r(1) = 10.0D0**(roundnicely(LOG10(yminax),5,"up"))
            ELSE
                yrayspec_r(1) = 10.0D0**(roundnicely(LOG10(yminax),5,"down"))
            ENDIF
        ENDIF
    ENDIF
    DO i=1,curvenr
        IF (cwaxis(i) == 2) THEN !right axis
            IF (cltype(i) == 1) THEN !"line"
                yrayspec_r(2) = 1.0D0
                csymbtype(i)= -1
            ENDIF
        ENDIF
    ENDDO

    IF (specialyrange2) THEN
        yrayspec_r(1) = lowerylimit2
        yrayspec_r(2) = upperylimit2
        ymaxax = MIN(ymaxax, upperylimit2)
        yminax = MIN(yminax, lowerylimit2)
        IF (logy_r) THEN
            yrayspec_r(3) = lowerylimit2*10.001D0 !cover at least one order of magnitude in log10 scale
        ELSE
            yrayspec_r(3) = upperylimit2
        ENDIF
    ENDIF

    IF (specialxrange) THEN
        xminax = lowerxlimit
        xmaxax = upperxlimit
        xrayspec(1) = lowerxlimit
        xrayspec(2) = lowerxlimit
        xrayspec(3) = xmaxaxin
    ENDIF
    IF (logxaxis) THEN
        xrayspec(2) = ABS(xrayspec(2))
        IF (xrayspec(1) <= 0.0D0) THEN
            xrayspec(1) = xrayspec(2)*1.0D-8
        ENDIF
        CALL AXSSCL ("LOG", "X")
        CALL LOGTIC ("FULL") !"AUTO"
    ELSE
        CALL AXSSCL ("LIN", "X")
    ENDIF

    !automatic scaling of x- and y-axis according to selected data:
    CALL SETSCL(yrayspec_r, 3, "Y")
    CALL SETSCL(xrayspec, 3, "X")

    !Second axis system for plot:
    CALL NAMDIS (44, "Y")
    !Set second y-axis label:
    CALL SETGRF ("LINE", "LINE", "LINE", "NAME")
    CALL TICKS (4, "Y") !the number of ticks between axis labels (default: 2)
    IF ((.NOT. logy_r) .AND. (ABS(ymaxax-yminax) < 1.0D-3 .OR. ABS(ymaxax-yminax) > 9.99999D3)) THEN
        CALL LABELS ("FEXP", "Y") !use powers of 10 "scientific" number format for axis labels
    ELSE
        CALL LABELS ("FLOAT", "Y")
    ENDIF
    CALL HNAME(Namesize)
    CALL NAME (TRIM(rylabel), "Y")
    CALL RVYNAM
    CALL HEIGHT(Numbsize)
    CALL CHASPC(-0.07D0)
    !plot axis and ticks from (xmin,xmax,first_x_label,increment_x / ymin,ymax,first_y_label,increment_Y):
    CALL NEGLOG(DEPS)
    xlab = xminax; ylab = yminax;
    IF (xmaxax > 1.0D0 .AND. barpmodeno > 1 .AND. xstep > 0.99D0) THEN !LLE relative diff.
        CALL INTAX
        yminax = 0.0D0; ymaxax = 1.0D0; ystep = 0.1D0;
        CALL GRAF(xminax,xmaxax,xminax,xstep, yminax,ymaxax,yminax,ystep)
    ELSE IF (xmaxax > 1.0D0 .OR. xminax > 0.0D0 .OR. (xmaxax-xminax < 0.4D0)) THEN !x-axis from xminax to xmaxax
        CALL GRAF(xminax,xmaxax,xminax,xstep, yminax,ymaxax,yminax,ystep) !will actually be scaled automatically using SETSCL
    ELSE
        IF (squareplot) THEN
            CALL GRAF(xminax,xmaxax,xminax,xstep, yminax,ymaxax,yminax,ystep) !will actually be scaled automatically using SETSCL
        ELSE
            CALL GRAF(xminax,xmaxax,xminax,xstep, yminax,ymaxax,yminax,ystep) !will actually be scaled automatically using SETSCL
        ENDIF
    ENDIF

    !plot second Gridlines if the scaling of the right y-axis is different from the left one:
    IF (ABS(MAXVAL(yrayspec_r)-MAXVAL(yrayspec_l)) > 1.0D-12 .OR. ABS(MINVAL(yrayspec_r)-MINVAL(yrayspec_l)) > 1.0D-12) THEN
        IF (.NOT. nogridlines) THEN
            CALL SETVLT("GREY")
            IF (logy_r) THEN
                IF (.NOT. logy_l) THEN !maybe also then the log-gridlines should be reduced to the main ones to enhance the visualisation!!
                    CALL SETCLR(230)
                    NRAY(1) = 12  !Stift unten
                    NRAY(2) = 30  !Stift oben
                    CALL MYLINE (NRAY, 2)
                    CALL PENWID(0.1D0) !the linewidth (allowing values smaller than 1)
                    CALL GRID(0,1)
                    CALL LINWID(1) !reset the linewidth
                ENDIF
            ELSE
                CALL SETCLR(160)
                NRAY(1) = 12  !pen down
                NRAY(2) = 30  !pen up
                CALL MYLINE (NRAY, 2)
                CALL PENWID(0.2D0) !the linewidth (allowing values smaller than 1)
                CALL GRID(0,1)
                CALL LINWID(1) !reset the linewidth
            ENDIF
        ENDIF
        CALL SOLID
    ENDIF

    !calculate and initialize number of lines of 2nd legend:
    CBUF = ""
    CSTR = ""
    NLIN = COUNT(clegplot(1:curvenr) .AND. cwaxis(1:curvenr) == 2)
    IF (NLIN > 30) THEN !dislin currently can only store 30 curve (attributes) for legends
        CALL LEGINI(CBUF, 30, NMAXLN)
    ELSE
        CALL LEGINI(CBUF, NLIN, NMAXLN)
    ENDIF
    !check whether a special legend plotting mode is necessary to store curve attributes and change curve color of white curves for visible contrast in legend:
    IF (ANY(clegplot(1:curvenr) .AND. cwaxis(1:curvenr) == 2 .AND. ccolor(1:curvenr) == 255)) THEN
        uselegpat = .true.
    ELSE
        uselegpat = .false.
    ENDIF

    !plot the curves on the 2nd axis ###
    ILIN = 0
    IF (ANY(.NOT. clegplot(1:curvenr))) THEN
        legstop = 2
    ELSE
        legstop = 1
    ENDIF
    DO ileg = 1,legstop !loop two times over the curves, the first time to only plot curves with legend entries.
        DO i = 1,curvenr
            IF ((ileg == 1 .AND. clegplot(i)) .OR. (ileg == 2 .AND. (.NOT. clegplot(i)))) THEN
                IF (cwaxis(i) == 2) THEN !plot it on the right y-axis
                    !plot the defined curves:
                    IF (customCol) THEN !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
                        CALL VLTFIL(TRIM(customcolfile), "LOAD") !the Win XP Icon Color Palette plus Grayscale etc. as defined by Andi Zuend
                    ELSE
                        CALL SETVLT("RAIN") !the Dislin rainbow palette
                    ENDIF
                    CALL SETCLR(ccolor(i))  !curve colour
                    CALL THKCRV(1)  !curve thickness
                    CALL PENWID(cth(i)) !the linewidth (allowing also values smaller than 1)
                    CALL INCMRK(0) !line
                    IF (cltype(i) == 2) THEN !symbols and line
                        CALL INCMRK(1)  !line connecting data points at every n-th point
                        CALL MARKER(csymbtype(i))   !plot symbol type at point coordinates
                        IF (smallsymbs) THEN
                            IF (cBigDsymbs(i)) THEN
                                CALL PENWID(1.0D0*cth(i)) !the linewidth (allowing also values smaller than 1)
                                CALL HSYMBL(bigsymbsize) !symbol size
                            ELSE
                                CALL PENWID(0.5D0*cth(i)) !the linewidth (allowing also values smaller than 1)
                                CALL HSYMBL((symbsize)/2) !symbol size
                            ENDIF
                        ELSE
                            IF (cBigDsymbs(i)) THEN
                                CALL PENWID(1.3D0*cth(i)) !the linewidth (allowing also values smaller than 1)
                                CALL HSYMBL(bigsymbsize) !symbol size
                            ELSE
                                CALL PENWID(1.0D0*cth(i)) !the linewidth (allowing also values smaller than 1)
                                CALL HSYMBL(symbsize) !symbol size
                            ENDIF
                        ENDIF
                    ELSE IF (cltype(i) > 2) THEN !symbols without line
                        CALL INCMRK(-1)  !no line connecting data points
                        CALL MARKER(csymbtype(i))   !plot symbol type at point coordinates
                        IF (smallsymbs) THEN
                            IF (cBigDsymbs(i)) THEN
                                CALL HSYMBL(bigsymbsize) !symbol size
                                CALL PENWID(1.3D0*cth(i)) !the linewidth (allowing also values smaller than 1)
                            ELSE
                                CALL HSYMBL((symbsize)/2) !symbol size
                                CALL PENWID(1.0D0*cth(i)) !the linewidth (allowing also values smaller than 1)
                            ENDIF
                        ELSE
                            CALL PENWID(1.0D0*cth(i)) !the linewidth (allowing also values smaller than 1)
                            IF (cBigDsymbs(i)) THEN
                                CALL HSYMBL(bigsymbsize) !symbol size
                            ELSE
                                CALL HSYMBL(symbsize) !symbol size
                            ENDIF
                        ENDIF
                    ENDIF
                    SELECT CASE (clstyle(i)) !define the line style
                    CASE(1)
                        CALL LNCAP("ROUND")! ("LONG") !rounded line caps
                        CALL SOLID
                        NTYP = 0
                    CASE(2)
                        CALL LNCAP("CUT") !rounded line caps
                        CALL DASHM
                        NTYP = 5
                    CASE(3)
                        CALL LNCAP("ROUND") !rounded line caps
                        CALL DASH
                        NTYP = 2
                    CASE(4) !MyDot style
                        NRAY(1) = 1  !Pen down
                        NRAY(2) =  INT(1.5D0*MAX(cth(i), 4.0D0)) !2*MAX(3, INT(cth(i))) !+4  !Pen up
                        CALL LNCAP("ROUND") !rounded line caps
                        CALL MYLINE (NRAY, 2) !self defined line style "DOT" with size according to line thickness
                        NTYP = 6
                    CASE(5)
                        CALL LNCAP("ROUND") !rounded line caps
                        CALL DOT
                        NTYP = 1
                    CASE(6) !. - . - . my dash-dotted line style
                        INRAY(1) = 1  !Pen down
                        INRAY(2) = 3*MAX(3, INT(cth(i)))   !Pen up
                        INRAY(3) = 4*MAX(3, INT(cth(i)))   !Pen down
                        INRAY(4) = 3*MAX(3, INT(cth(i)))   !Pen up
                        CALL LNCAP("ROUND") !rounded line caps
                        CALL MYLINE (INRAY(1:4), 4)
                        NTYP = 4
                    END SELECT
                    CALL LNJOIN ("SHARP") !'SHARP' or 'TRUNC'
                    CALL NOCHEK  !suppress warning of points lying outside of the plotting area
                    !now plot the selected curve:
                    N = rownr(i)
                    xax(1:N) = cxvals(i,1:N)
                    yax(1:N) = cyvals(i,1:N)
                    CALL CURVE(xax(1:N),yax(1:N),N)
                    IF (ileg == 1) THEN !set the corresponding legend text and color:
                        CSTR = TRIM(clegtxt(i))
                        IF (LEN_TRIM(CSTR) > 1) THEN
                            IF (ILIN < 30) THEN !only max. 30 curves attributes can be printed for the legend (current dislin limit)
                                ILIN = ILIN+1
                                !!change color in legends when the line color is white (to have some contrast to the white  background color):
                                IF (uselegpat) THEN
                                    IF (ccolor(i) == 255) THEN !white
                                        iclr = 0 !black
                                    ELSE
                                        iclr = ccolor(i)
                                    ENDIF
                                    ityp = NTYP !select a line style from the default styles, closest to the one of the custom types.
                                    CALL LEGPAT(ityp, INT(cth(i)), csymbtype(i), iclr, -1, ILIN) !ILIN
                                ENDIF
                                CALL LEGVAL (1.0D0, 'SYMBOL') !symbol size in legend (default = 0.8D0)
                                CALL LEGLIN(CBUF, TRIM(CSTR), ILIN) !define the legend line text
                            ENDIF
                        ENDIF
                    ENDIF
                    IF (ANY(cyerr(i,1:N) > dtiny)) THEN
                        E1RAY(1:N) = cyerr(i,1:N) !error bars in y-axis
                        E2RAY(1:N) = cyerr(i,1:N)
                        CALL LINWID(1)
                        IF (cltype(i) < 2) THEN !no symbols, just line
                            CALL PENWID(MIN(0.3D0, MAX(0.2D0, cth(i)/5.0D0)))
                            CALL HSYMBL(MAX(1,INT(cth(i)/3.0D0)))   !CALL HSYMBL(1) !set symbol size to very tiny
                            CALL MARKER(-1) !CALL MARKER (3)
                        ELSE
                            CALL PENWID(MIN(1.0D0, MAX(0.3D0, cth(i)/3.0D0)))
                            CALL MARKER(-1) !CALL MARKER (csymbtype(i))
                        ENDIF
                        CALL ERRBAR(xax, yax, E1RAY, E2RAY, N)
                    ENDIF
                    CALL LINWID(1) !reset the linewidth
                    CALL PENWID(1.0D0)
                ENDIF !axis IF
            ENDIF !legend IF
        ENDDO !plot curves 2nd axis
    ENDDO !ileg

    !plot legend to the second axis system:
    IF (LEN_TRIM(CBUF) > 1) THEN
        CALL SETVLT("RAIN")
        CALL COLOR("BLACK")
        CALL SOLID
        CALL HEIGHT (Lsize)
        CALL LEGTIT ("right y-axis:             ")
        CALL LINESP(2.6D0)
        IF (showcustomcolbar) THEN
            nxpos = 160
        ELSE
            nxpos = 0
        ENDIF
        IF (squareplot) THEN
            CALL LEGPOS (1900+nxpos, 900)   !position the legend 2 below the other one
        ELSE
            CALL LEGPOS (2100+nxpos, 900)   !position the legend on the right of the other one
        ENDIF
        CALL FRAME(0)  !no box around the legend
        IF (customCol) THEN !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
            CALL VLTFIL(TRIM(customcolfile), "LOAD") !the Win XP Icon Color Palette plus Grayscale etc. as defined by Andi Zuend
        ELSE
            CALL SETVLT("RAIN") !the Dislin rainbow palette
        ENDIF
        CALL SETCLR(0) !black
        CALL LINWID(MIN(3, INT(MAXVAL(cth(1:curvenr))))) !set here curve thickness of curve / points for all legend entries (as no curve-specific choice possible for now).
        CALL LEGEND (CBUF, 3)
        CALL LINWID(1) !reset
    ENDIF

ENDIF !2ndaxis
!2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd
!2nd                                       end of 2nd y-axis                                  2nd
!2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd2nd

CALL SETVLT("RAIN") !the Dislin rainbow palette
CALL COLOR("BLACK")
CALL SOLID
CALL FRAME(1)
CALL BOX2D !plots a box around the axis system

!Plot a number at coordinates NX,NY
IF (TRIM(ptitle(3)) == "SD of fit:") THEN !plot standard deviation of fit, SD:
    CALL GETPOS (NXA, NYA)
    CALL GETLEN (NXL, NYL, NZL)
    CALL NUMBER (StandardDev, 10, NXA+240, NYA-NYL-126)   !with 10 digits
ENDIF

!------------------------------------------------------------------------------------------------
!plot the message box (if there was any line of text defined):
IF (messagebox) THEN
    !set the (initial) plot coordinates and properties of the message box:
    IF (messboxpos_x == 0 .AND. messboxpos_y == 0) THEN !the position was not set yet, so use default values
        nxl = 1200
        nyl = 120
    ELSE
        nxl = messboxpos_x
        nyl = messboxpos_y
    ENDIF
    CALL HEIGHT (Lsize)
    CALL TXTJUS ("LEFT")
    CALL FRMESS (0) !the thickness of the frame around the box
    !write out message box header horizontal line:
    N = MAXVAL(LEN_TRIM(messagetxt(1:txtl)))
    hline = REPEAT("-", INT(N*5.0D0))
    CALL CHASPC(-0.4D0)
    CALL MESSAG (TRIM(hline), nxl, nyl)
    nyl = nyl + INT(Lsize*2.3D0)
    CALL CHASPC(-0.01D0)
    !write out the saved textlines:
    DO i=1,txtl
        CALL MESSAG ("   "//TRIM(messagetxt(i)), nxl, nyl)
        nyl = nyl + INT(Lsize*2.3D0)
    ENDDO
    !write out message box footer horizontal line:
    CALL CHASPC(-0.4D0)
    CALL MESSAG (TRIM(hline), nxl, nyl)
ENDIF
!------------------------------------------------------------------------------------------------

!END DISLIN plotting:
CALL WINMOD("NONE") !for graphic output on XWIN only.:
CALL TEXMOD("OFF")
CALL RESET('ALL') !reset all settings to defaults
CALL DISFIN !end DISLIN

!reset some switches:
newplot = .true.
cxvals = 0.0D0
cyvals = 0.0D0
cyerr = 0.0D0
CALL DefinePlotData(cxvals(1,:),cyvals(1,:),cyerr(1,:),1,0,1,-1,3,clegtxt(1),1.0D0) !call here with no actual data to ensure reset of stored curve data.
newplot = .true. !set it back for next plot page
newplot2 = .true.
messagebox = .false.
smallsymbs = .false.
bigsymbs = .false.
specialxrange = .false.
specialyrange = .false.
specialyrange2 = .false.
logyaxis = .false. !set the logarithmic y-axis to false (so use linear).
logxaxis = .false.
ptitle = ""
customCol = .false.
bargraphplot = .false.
fixminxax = .true.
xorg3bar = .false.
legplot = .true.
nogridlines = .false.
showcustomcolbar = .false.

END SUBROUTINE PlotNow
   
 
!=================================================================================================================
!=================================================================================================================


SUBROUTINE PlotPiechart(ptitle,outname_cplot,nsemivol,nelectrol)

USE dislin  !the DISLIN plotting library (www.dislin.de)

IMPLICIT NONE
!...
INTEGER(4) :: NMAXLN
PARAMETER (NMAXLN = 75)
!...
INTEGER(4) :: i,N,NXL,NYL,symbsize,xpos,ypos,xposstart,yposstart,xlen,ylen,k,lcount
INTEGER(4) :: Tsize,Lsize,Namesize,Numbsize,pielablesize,nsemivol,nelectrol
INTEGER(4),DIMENSION(NTCP) :: piecolset
!...
REAL(8),DIMENSION(nrows) :: xax,yax
!...
CHARACTER(LEN=75),INTENT(IN) :: outname_cplot
CHARACTER(LEN=75) :: outname_dislin
CHARACTER(LEN=(nsemivol+nelectrol+1)*NMAXLN) :: CBUFPIE
CHARACTER(LEN=132),DIMENSION(4) :: ptitle
CHARACTER(LEN=messagetextlength) :: hline
!...
LOGICAL(4) :: ryaxis, plotpielables
!...
EXTERNAL :: MyPieSettings
!.........................................................................................

!set some variables to zero:
xax = 0.0D0
yax = 0.0D0

!check whether there is any curve/points to plot at all:
IF (curvenr < 1) THEN !there is nothing to plot
    !$OMP CRITICAL
    WRITE(*,*),"There was no curve to plot in this call to PlotPiechart!"
    WRITE(*,*),"proceeding without producing a plot..."
    WRITE(*,*),""
    !$OMP END CRITICAL
    !set the paramaters back for next plot page
    newplot = .true. 
    newplot2 = .true.
    RETURN !leave the subroutine
ENDIF

!:: adjustable  plot parameters ::  
!Set font sizes etc.: 
symbsize = 32 !40 !26 !symbol size
Tsize = 26   !Title font size
Lsize = 20   !Legends font size
plotpielables = .false. !##############################
pielablesize = 13 !12 !14  !the font size of piechart lables
Namesize = 40  !axis name font
Numbsize = 32  !axis numbers size
customCol = .true. !load custom color palette

!set output file
CALL METAFL(metaff) 
outname_dislin=TRIM(outname_cplot)//"."//(TRIM(metaff))
outname_dislin="../Output_Plots/"//TRIM(outname_dislin)
CALL SETFIL(TRIM(outname_dislin)) !set filename in dislin
CALL FILMOD("DELETE") !files will be overwritten if they have the same name..
CALL SCRMOD ("NOREV")  !set background color to white and foreground to black
CALL SETPAG("DA4L")
CALL SCLFAC(1.0D0)
CALL SCLMOD("DOWN")
CALL IMGFMT ("RGB")
CALL WINSIZ (1600, 1200)
IF (metaff == "png" .OR. metaff == "gif") THEN
    CALL IMGFMT ("RGB")
    CALL WINSIZ (1600, 1200)
ENDIF
    
!Set the overall plot parameters and initialize the plot page:
CALL DISINI !initialize dislin
IF (customCol) THEN !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
    CALL VLTFIL(TRIM(customcolfile), "LOAD") !the Win XP Icon Color Palette plus Grayscale etc. as defined by Andi Zuend
ELSE
    CALL SETVLT("RAIN") !the Dislin rainbow palette
ENDIF
CALL PAGFLL(255) !background color white
CALL SETCLR(0) !black
CALL ERRMOD ("CHECK", "OFF")  !Don't pass warnings of points lying outside of the coordinate systems
CALL TEXMOD ("ON") !allow some LaTeX code for better text

!font style:
IF (metaff == "wmf") THEN
    CALL WINFNT("Arial")  !FONT in Windows  !Arial  ! Times New Roman
ELSE IF (metaff == "eps" .OR. metaff == "ps" .OR. metaff == "pdf" .OR. metaff == "svg") THEN
    CALL PSFONT("Helvetica")  !for ps, eps and pdf output
    CALL PSMODE('BOTH')  !allow Greek and Italic modes
ELSE
    CALL COMPLX
ENDIF

!Set Titeltext on line no. 1:
CALL TITLIN(TRIM(ptitle(1)), 1)
!Set Titeltext on line no. 2:
CALL TITLIN(TRIM(ptitle(2)), 2)
!Set Titeltext on line no. 3:
CALL TITLIN(TRIM(ptitle(3)), 3)
!Set Titeltext on line no. 4:
CALL TITLIN(TRIM(ptitle(4)), 4)
CALL LINESP(2.6D0)
CALL TITJUS ("LEFT")

CALL AXSPOS(300, 500)
CALL AXSLEN(2000, 150)
CALL SETGRF ("NONE", "NONE", "NONE", "NONE")
CALL FRAME(0)
CALL GRAF(0.0D0,1.0D0,0.0D0,0.1D0,0.0D0,1.0D0,0.0D0,0.2D0) !pseudo graph for title plot
!plot title above the axis system:
CALL CHASPC(-0.06D0)
CALL SETCLR(0)
CALL HEIGHT(Tsize)
CALL TITLE
CALL ENDGRF

!date of print and dislin version on the lower right corner
CALL HEIGHT(Lsize)
CALL PAGHDR("","",2,0) 

!find out if there is a second y-axis...
ryaxis = .false.
DO i=1,curvenr
    IF (cwaxis(i) == 2) THEN
        ryaxis = .true.
    ENDIF
ENDDO

!Set axis labels:
!CALL GETPAG(NXL, NYL)  !get the page total DIMENSIONs NXL=2970, NYL = 2100
xposstart = 400
yposstart = 1000
ypos = yposstart
xlen = 300
ylen = 300

!change shading pattern cycle:
CALL PATCYC(1,1)
CALL PATCYC(2,10)
CALL PATCYC(3,16)

!define a color set for up to NTCP piechart segments:
piecolset(1:nsemivol) = colneutral(1:nsemivol)
piecolset(1) = colneutral(1)+2 !lighter blue for water (resp. component 1) color
piecolset(nsemivol+1:nsemivol+nelectrol) = colelectrol(1:nelectrol)
    
!1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie
!1stPie                               start of 1st set of pie charts                       1stPie
!1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie

!plot the chart on the 1st y position ###
lcount = 0
DO i=1,curvenr  
    IF (cwaxis(i) == 1) THEN !plot it on the first "upper" axis system level
        lcount = lcount+1
        !position the chart in the horizontal direction:
        xpos = xposstart + (lcount-1)*(xlen+500)
        CALL AXSPOS(xpos, ypos)
        CALL AXSLEN(xlen, ylen)  !now set axis length (pie chart size) on the page in plot coord.
        !set number of segments and the data:
        N = rownr(i)
        yax(1:N) = cyvals(i,1:N)
        !initialize Legend:  
        CBUFPIE = " "
        k = 1+MAXVAL(LEN_TRIM(cpielegtxt(i,1:N)))
        CALL LEGINI(CBUFPIE,N,k)
        DO k = 1,N
            CALL LEGLIN(CBUFPIE,TRIM(cpielegtxt(i,k)),k)
        ENDDO
        !now plot the selected pie chart:
        CALL SETCLR(0) !black
        CALL SOLID
        CALL HEIGHT (pielablesize)  !Legend font size
        CALL CHASPC(-0.06D0)
        IF (plotpielables) THEN
            CALL LABELS("PERCENT","PIE")  !"PERCENT" ; "DATA" ; "BOTH" ; "NONE"
        ELSE
            CALL LABELS("NONE","PIE") !to skip data in labels
        ENDIF
        CALL LABDIG(2, "PIE")
        CALL LABPOS('INTERNAL', 'PIE')  !"EXTERNAL" , 'INTERNAL'
        CALL FRAME(1) !0 = no frame around segment labels
        !define pie chart properties:
        CALL CHNPIE("NONE")  !don't cycle piesegment properties
        CALL SHDPAT(16) !solid fill
        CALL PIEVEC(0, 'BROKEN')
        k = 1
        CALL PIECLR (piecolset(k:k+N), piecolset(1:N), N)  !set the piesegment colours: -1 means use current color! 
        !check whether some pie label textboxes have to be shifted in the y-direction due to overlaps:
        pielabshift = .false.
        DO k = 1,N-1
            IF (SUM(yax(k:k+1)) < 0.20D0) THEN !shift only one of the labels
                IF ((yax(k) < 0.14D0 .AND. yax(k) > 0.08D0 .AND. yax(k+1) < 0.08D0) .OR. &
                & (yax(k+1) < 0.14D0 .AND. yax(k+1) > 0.08D0 .AND. yax(k) < 0.08D0)) THEN
                    pielabshift(k+1) = .true.
                ELSE IF (yax(k) < 0.10D0 .AND. yax(k+1) < 0.10D0) THEN !shift both labels
                    pielabshift(k:k+1) = .true.
                ENDIF
            ENDIF
        ENDDO
        DO k = 1,N
            CALL PIECBK(MyPieSettings) !Callbackroutine to shift the pie labels; defined at the end of this file.
        ENDDO
        CALL PIEROT(25.0D0)  !rotate the piechart by 25 degrees counterclockwise
        CALL LINESP(1.8D0)
        CALL PIEVAL(0.5D0, 'DIST')   !the distance scaling factor for the radial label position; default = 1.0D0
        CALL PENWID(0.05D0) !the linewidth (allowing also values smaller than 1)
        CALL CIRCSP(1) !smaller arc length to enable plotting of tiny segments
        CALL LNCAP ("CUT") ! CUT = cut off line caps
        IF (plotpielables) THEN
            CALL PIEGRF(CBUFPIE,1,yax,N)
        ELSE
            CALL PIEGRF("",0,yax,N) !for piecharts without labels plotted (when LABELS called with "NONE")
        ENDIF
        !Set and plot pie chart title:
        k = MAXVAL(LEN_TRIM(cpietittxt(i,1:4)))
        IF (k > 132) THEN
            PAUSE 'title text too long'
        ENDIF
        CALL LNCAP("ROUND")
        CALL PENWID(1.0D0) !the linewidth (allowing also values smaller than 1)
        CALL TITLIN (TRIM(cpietittxt(i,1)),1) !Pie chart title line 1
        CALL TITLIN (TRIM(cpietittxt(i,2)),2) !Pie chart title line 2
        CALL TITLIN (TRIM(cpietittxt(i,3)),3) !Pie chart title line 3
        CALL TITLIN (TRIM(cpietittxt(i,4)),4) !Pie chart title line 4
        CALL TITJUS ("LEFT")
        CALL SETCLR(0)
        CALL HEIGHT(Lsize)
        CALL LINESP(2.2D0)
        CALL TITLE
        CALL ENDGRF
    ENDIF !axis IF
ENDDO !i, plot curves 1st axis

!1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie
!1stPie                                end of 1st set of pie charts                        1stPie
!1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie1stPie	             

!2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie
!2ndPie                               start of 2nd set of pie charts                       2ndPie
!2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie

IF (ryaxis) THEN

!start positions
xposstart = 400
yposstart = 1800
ypos = yposstart

!plot the chart on the 2nd y position ###
lcount = 0
DO i=1,curvenr  
    IF (cwaxis(i) == 2) THEN !plot it on the 2nd "lower" y-axis level
        lcount = lcount+1
        !position the chart in the horizontal direction:
        xpos = xposstart + (lcount-1)*(xlen+500)
        CALL AXSPOS(xpos, ypos)
        CALL AXSLEN(xlen, ylen)  !now set axis length (pie chart size) on the page in plot coord.
        !set number of segments and the data:
        N = rownr(i)
        yax(1:N) = cyvals(i,1:N)
        !initialize Legend:  
        CBUFPIE = " "
        k = 1+MAXVAL(LEN_TRIM(cpielegtxt(i,1:N)))
        CALL LEGINI(CBUFPIE,N,k)
        DO k = 1,N
            CALL LEGLIN(CBUFPIE,TRIM(cpielegtxt(i,k)),k)
        ENDDO
        !now plot the selected pie chart:
        CALL SETCLR(0) !black
        CALL SOLID
        CALL HEIGHT (pielablesize)  !segment labels font size
        CALL CHASPC(-0.06D0)
        IF (plotpielables) THEN
            CALL LABELS("PERCENT","PIE")  !"PERCENT" ; "DATA" ; "BOTH" ; "NONE"
        ELSE
            CALL LABELS("NONE","PIE") !to skip data in labels
        ENDIF
        CALL LABDIG(2, "PIE")
        CALL LABPOS('INTERNAL', 'PIE')  !"EXTERNAL" , 'INTERNAL'
        CALL FRAME(1) !0 = no frame around segment labels
        !define pie chart properties:
        CALL CHNPIE('NONE')  !don't cycle piesegment properties
        CALL SHDPAT(16) !solid fill
        CALL PIEVEC(0, 'BROKEN')  !BROKEN  STRAIGHT
        k = 1
        CALL PIECLR (piecolset(k:k+N), piecolset(1:N), N)  !set the piesegment colours: -1 means use current color!
        !check whether some pie label textboxes have to be shifted in the y-direction due to overlaps:
        pielabshift = .false.
        DO k = 1,N-1
            IF (SUM(yax(k:k+1)) < 0.20D0) THEN !shift only one of the labels
                IF ((yax(k) < 0.14D0 .AND. yax(k) > 0.08D0 .AND. yax(k+1) < 0.08D0) .OR. &
                & (yax(k+1) < 0.14D0 .AND. yax(k+1) > 0.08D0 .AND. yax(k) < 0.08D0)) THEN
                    pielabshift(k+1) = .true.
                ELSE IF (yax(k) < 0.10D0 .AND. yax(k+1) < 0.10D0) THEN !shift both labels
                    pielabshift(k:k+1) = .true.
                ENDIF
            ENDIF
        ENDDO
        CALL PIEROT(25.0D0)
        DO k = 1,N
            CALL PIECBK(MyPieSettings) !Callbackroutine to shift the pie labels; defined at the end of this file.
        ENDDO
        CALL PIEVAL(0.5D0, 'DIST')   !the distance scaling factor for the radial label position; default = 1.0
        CALL LINESP(1.8D0)
        CALL PENWID(0.05D0) !the linewidth (allowing also values smaller than 1)
        CALL CIRCSP(1) !smaller arc length to enable plotting of tiny segments
        CALL LNCAP("CUT")
        IF (plotpielables) THEN
            CALL PIEGRF(CBUFPIE,1,yax,N)
        ELSE
            CALL PIEGRF("",0,yax,N) !for piecharts without labels plotted (when LABELS called with "NONE")
        ENDIF
        !Set and plot pie chart title:
        k = MAXVAL(LEN_TRIM(cpietittxt(i,1:4)))
        IF (k > 132) THEN
            PAUSE 'title text too long'
        ENDIF
        CALL PENWID(1.0D0) !the linewidth (allowing also values smaller than 1)
        CALL LNCAP("ROUND")
        CALL TITLIN (TRIM(cpietittxt(i,1)),1) !Pie chart title line 1
        CALL TITLIN (TRIM(cpietittxt(i,2)),2) !Pie chart title line 2
        CALL TITLIN (TRIM(cpietittxt(i,3)),3) !Pie chart title line 3
        CALL TITLIN (TRIM(cpietittxt(i,4)),4) !Pie chart title line 4
        CALL TITJUS ("LEFT")
        CALL SETCLR(0)
        CALL HEIGHT(Lsize)
        CALL LINESP(2.2D0)
        CALL TITLE
        CALL ENDGRF
    ENDIF !axis IF
ENDDO !i, plot curves 1st axis

ENDIF !ryaxis

!2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie
!2ndPie                               end of 2nd set of pie charts                         2ndPie
!2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie2ndPie
  
!------------------------------------------------------------------------------------------------
!plot the message box (if there was any line of text defined):
IF (messagebox) THEN
    CALL SETCLR(0) !black
    !set the (initial) plot coordinates and properties of the message box:
    IF (messboxpos_x == 0 .AND. messboxpos_y == 0) THEN !the position was not set yet, so use default values
        nxl = 1200
        nyl = 120
    ELSE
        nxl = messboxpos_x
        nyl = messboxpos_y
    ENDIF
    CALL HEIGHT (Lsize)
    CALL TXTJUS ("LEFT")
    CALL FRMESS (0) !the thickness of the frame around the box
    !write out message box header horizontal line:
    N = MAXVAL(LEN_TRIM(messagetxt(1:txtl)))
    hline = REPEAT("-", INT(N*5.0D0))
    CALL CHASPC(-0.4D0)
    CALL MESSAG (TRIM(hline), nxl, nyl)
    nyl = nyl + INT(Lsize*2.3D0)
    CALL CHASPC(-0.01D0)
    !write out the saved textlines:
    DO i=1,txtl 
        CALL MESSAG ("   "//TRIM(messagetxt(i)), nxl, nyl)
        nyl = nyl + INT(Lsize*2.3D0)
    ENDDO
    !write out message box footer horizontal line:
    CALL CHASPC(-0.4D0)
    CALL MESSAG (TRIM(hline), nxl, nyl)
ENDIF
!------------------------------------------------------------------------------------------------

!END DISLIN plotting: 
CALL WINMOD("NONE") !for graphic output on XWIN only.:
CALL TEXMOD("OFF")
CALL DISFIN !end DISLIN

!reset some switches:
newplot = .true. !set it back for next plot page
newplot2 = .true.
messagebox = .false.
customCol = .false.
bargraphplot = .false.
ptitle = ""
legplot = .true.
    
END SUBROUTINE PlotPiechart
   
!=================================================================================================================
!=================================================================================================================

SUBROUTINE PlotBarGraphs(ptitle,xlabel,outname_cplot,nsemivol,nelectrol)

USE dislin  !the DISLIN plotting library (www.dislin.de)
USE qsort_c_module

IMPLICIT NONE
!...
INTEGER(4) :: i,N,NXL,NYL,NZL,symbsize,nax,nc,inc
INTEGER(4) :: NMAXLN, ypos, ylen, iticks, xlen, ydist
!...
PARAMETER (NMAXLN = 75)
!...
INTEGER(4) :: NXA,NYA,Tsize,Lsize,Namesize,Numbsize,Lcompsize,ifault, &
    & nsteps,nsemivol,nelectrol,Lheight,allocstat
!...
REAL(8) :: xmaxax,xminax,roundnicely,logymin,xmaxaxin,xminaxin,xstep, &
    & yminax,ymaxax,ystep,davg
REAL(8),DIMENSION(3) :: yrayspec_l, yrayspec_r
REAL(8),DIMENSION(2) :: xrayspec
REAL(8),DIMENSION(nrows) :: xax,yax,xv,y0,xv2
REAL(8),DIMENSION(2*NTCP,nrows) :: ysum
REAL(8),DIMENSION(:),ALLOCATABLE :: xray, yray1, yray2
!...
CHARACTER(LEN=75),INTENT(IN) :: outname_cplot
CHARACTER(LEN=75) :: outname_dislin
CHARACTER(LEN=MIN(MAX(MAXVAL(ncomp(:))*NMAXLN, (nsemivol+nelectrol+1)*NMAXLN), 30*NMAXLN)) :: CBUF
CHARACTER(LEN=132),DIMENSION(4) :: ptitle
CHARACTER(LEN=75) :: xlabel
CHARACTER(LEN=messagetextlength) :: hline
!...
LOGICAL(4) :: newgraph
!...
EXTERNAL :: roundnicely

!.........................................................................................
!set some variables to zero:
yrayspec_l = 0.0D0
yrayspec_r = 0.0D0
xrayspec = 0.0D0
xax = 0.0D0
yax = 0.0D0
newgraph = .true.
  
!check whether there is any curve/points to plot at all:
IF (curvenr < 1 .AND. barpmodeno < 1) THEN !there is nothing to plot
    !$OMP CRITICAL
    WRITE(*,*),"There was no curve to plot in this call to PlotBarGraphs!"
    WRITE(*,*),"proceeding without producing a plot..."
    WRITE(*,*),""
    !$OMP END CRITICAL
    newplot = .true. !set it back for next plot page
    newplot2 = .true.
    messagebox = .false.
    smallsymbs = .false.
    bigsymbs = .false.
    specialxrange = .false.
    specialyrange = .false.
    specialyrange2 = .false.
    logyaxis = .false. !set the logarithmic y-axis to false (so use linear).
    logxaxis = .false.
    ptitle = ""
    customCol = .false.
    bargylimits = .false.
    plotbarborder = .false.
    RETURN !leave the subroutine
ENDIF

!:: adjustable  plot parameters ::  
!squareplot scales the x, y axis length to plot in a square shaped box
    
!Set font sizes etc.: 
symbsize = 38 !symbol size
Tsize = 24   !Title font size
Lsize = 19   !Legends font size
Lcompsize = 24  !composition plot legend size
IF (.NOT. squareplot) THEN
    IF (barpmodeno < 3) THEN
        Namesize = 48  !axis name font
        Numbsize = 36  !axis numbers size
    ELSE
        Namesize = 40  !axis name font
        Numbsize = 34  !axis numbers size
    ENDIF
ELSE 
    Namesize = 36  !axis name font
    Numbsize = 32  !axis numbers size
ENDIF
logymin = 1.0D-3  !the minimum y value when logarithmic y-axis are plottet
customCol = .true.
    
!set output file
CALL METAFL(metaff) 
outname_dislin = TRIM(outname_cplot)//"."//(TRIM(metaff))
outname_dislin = "../Output_Plots/"//TRIM(outname_dislin)
CALL SETFIL(TRIM(outname_dislin)) !set filename in dislin
CALL FILMOD("DELETE") !files will be overwritten if they have the same name..
CALL SCRMOD ("NOREV")  !set background color to white and foreground to black
CALL SETPAG("DA4L")
CALL SCLFAC(1.0D0)
CALL SCLMOD("DOWN")
CALL IMGFMT ("RGB")
CALL WINSIZ (1600, 1200)
IF (metaff == "png" .OR. metaff == "gif") THEN
    CALL IMGFMT ("RGB")
    CALL WINSIZ (1200, 900)
ENDIF

!axis and label scaling according to the number of bar graphs:
ydist = MAX(10, 100-(barpmodeno-2)*15)
ylen = (2100-400-300-(barpmodeno-1)*ydist)/barpmodeno !the y-axis length per axis system
IF (squareplot) THEN
    xlen = MIN(1200, NINT(ylen*3.0D0)) !scale the length of the x-axis charts:
    Namesize = MAX(12, 34-INT((barpmodeno-1)*2.0D0))
    Numbsize = MAX(11, 32-INT((barpmodeno-1)*2.0D0))
ELSE
    xlen = MIN(1800, NINT(ylen*4.0D0))
    Namesize = MAX(12, 36-INT((barpmodeno-1)*2.1D0))
    Numbsize = MAX(11, 34-INT((barpmodeno-1)*2.1D0))
ENDIF
    
!Set the over all plot parameters and initialize the plot page:
CALL DISINI !initialize dislin
CALL PAGFLL(255) !background color white
CALL SETVLT("RAIN")  !Rainbow colors (256)
CALL SETCLR(0) !black
IF (metaff == "wmf") THEN
    CALL WINFNT("Arial")  !FONT in Windows  !Arial  ! Times New Roman
ELSE IF (metaff == "eps" .OR. metaff == "ps" .OR. metaff == "pdf" .OR. metaff == "svg") THEN
    CALL PSFONT("Helvetica")  !for ps, eps and pdf output
    CALL PSMODE('BOTH')  !allow Greek and Italic modes
ELSE
    CALL COMPLX
ENDIF
!date of print, time and dislin version on the lower right corner
CALL HEIGHT(Lsize)
CALL PAGHDR("","",2,0) 
CALL ERRMOD ("CHECK", "OFF")  !Don't pass warnings of points lying outside of the coordinate systems
CALL CHASPC(-0.06D0)
CALL HEIGHT(Namesize) 
!enable some LaTeX code for better axes-text:
CALL TEXMOD ("ON")

!Set Titeltext on line no. 1:
CALL TITLIN(TRIM(ptitle(1)),1)
!Set Titeltext on line no. 2:
CALL TITLIN(TRIM(ptitle(2)),2)
!Set Titeltext on line no. 3:
CALL TITLIN(TRIM(ptitle(3)),3)
!Set Titeltext on line no. 4:
CALL TITLIN(TRIM(ptitle(4)),4)

!find overall x-axis scaling applicable to all related bar graphs: 
IF (specialxrange) THEN
    xminaxin = lowerxlimit
    xmaxaxin = upperxlimit 
ELSE
    xminaxin = compbx(1,1)
    xmaxaxin = compbx(1,1)
    DO nax = 1,barpmodeno
        N = MAXVAL(compbpoints(nax:nax))
        xv = 0.0D0
        xv(1:N) = compbx(nax,1:N)
        xminaxin = MIN(xminaxin, MINVAL(xv(1:N), MASK = xv > 0.0D0))
        xmaxaxin = MAX(xmaxaxin, MAX(xminaxin*1.20001D0, MAXVAL(xv(1:N), MASK = xv > 0.0D0)))
    ENDDO   
ENDIF

DO nax = 1,barpmodeno !the number of different axis-systems
    !set the axis labels for the current mode
    CALL NAMDIS (Namesize+10, "Y")  !distance between axis text and axis names
    CALL NAMDIS (Namesize+10, "X")  
    CALL HNAME (Namesize)
    CALL NAME (TRIM(xlabel), "X")
    CALL NAME (TRIM(clylabel(nax)), "Y")
    !position the axis system (the lower left axis corner):
    IF (squareplot) THEN
        ypos = 400 +(nax-1)*70 +ylen*nax
        CALL AXSPOS(300, ypos)
        CALL AXSLEN(xlen, ylen)
        newgraph = .true.
    ELSE
        ypos = 400 +(nax-1)*70 +ylen*nax !ypos = 400 +(nax-1)*100 +ylen*nax
        CALL AXSPOS(300, ypos)
        CALL AXSLEN(xlen, ylen)
        newgraph = .true.
    ENDIF
    N = MAXVAL(compbpoints(nax:nax))
    xv = 0.0D0
    !set the number of entries and the x-values:
    xv(1:N) = compbx(nax,1:N)
    CALL SETCLR(0) !black
    IF (newgraph) THEN !create new axis system 
        !check for appropriate number of steps (scale intervals):
        IF (ABS(xmaxaxin-xminaxin) < 2.0D0) THEN
            IF (xmaxaxin < 1.001D0 .AND. xmaxaxin > 0.998D0) THEN
                xmaxaxin = 1.0D0
            ENDIF
            nsteps = roundnicely(ABS(xmaxaxin-xminaxin)*10.0D0,5,"up")
            nsteps = MAX(nsteps, 1)
            IF (nsteps > 10) THEN
                DO WHILE (nsteps > 10)
                    nsteps = roundnicely((DFLOAT(nsteps)/2.0D0),5,"up")
                ENDDO
            ELSE IF (nsteps < 4) THEN
                nsteps = roundnicely((DFLOAT(nsteps)*2.0D0),5,"up")   
            ENDIF
            CALL scalenicely(xminaxin, xmaxaxin, nsteps, xminax, xmaxax, xstep, ifault)
        ELSE
            nsteps = roundnicely(ABS(xmaxaxin-xminaxin)*10.0D0,5,"up")
            nsteps = MAX(nsteps, 1)
            IF (nsteps > 12) THEN
                DO WHILE (nsteps > 12)
                    nsteps = roundnicely((DFLOAT(nsteps)/2.0D0),5,"up")
                ENDDO
            ELSE IF (nsteps < 4) THEN
                nsteps = roundnicely((DFLOAT(nsteps)*2.0D0),5,"up")   
            ENDIF
            IF (nsteps > INT(xmaxaxin)) THEN
                nsteps = INT(xmaxaxin)
            ENDIF
            DO
                IF (MOD(INT(xmaxaxin-xminaxin),nsteps) == 0) THEN
                    EXIT
                ELSE
                    xmaxaxin = AINT(xmaxaxin)+1.0D0
                ENDIF
            ENDDO
            CALL scalenicely(xminaxin, xmaxaxin, nsteps, xminax, xmaxax, xstep, ifault)
            IF (xmaxaxin-xminaxin+1.0D0 > 2.0D0) THEN
                xminax = AINT(xminax)
                xmaxax = AINT(xmaxax)
                xstep = MAX(1.0D0, AINT(xstep))
                xstep = MIN(AINT(xstep), AINT(xmaxax))
            ENDIF
        ENDIF
        !prepare graph:
        yrayspec_r = 0.0D0
        yrayspec_r(1) = MIN(0.0D0, roundnicely(MINVAL(compby(nax,:,1:N)),5,"down")) !xminax !0.0D0
        yrayspec_r(2) = min_y
        DO i = 1,N
            yrayspec_r(3) = MAX(yrayspec_r(3), roundnicely(SUM(compby(nax,:,i))-1.0D-5,5,"up"))
        ENDDO
        IF (specialyrange) THEN
            IF (bargylimits(nax)) THEN
                yrayspec_r(1) = lowerylimit
                yrayspec_r(2) = upperylimit
                yrayspec_r(3) = upperylimit
            ENDIF
        ENDIF
        CALL AXSSCL ("LIN", "Y") !LIN
        ymaxax = MAXVAL(yrayspec_r)
        yminax = MINVAL(yrayspec_r)
        CALL SETSCL(yrayspec_r, 3, "Y") !scales the axis automatically for GRAF
        !x-axis scaling and mode:
        IF (logxaxis) THEN
            IF (specialxrange) THEN
                xminax = lowerxlimit
                xmaxax = upperxlimit
                xrayspec(1) = xminax*0.999D0
                xrayspec(2) = xmaxax*1.001D0 
            ELSE
                xminax = MINVAL(xv(1:N), MASK = xv > 0.0D0)
                xmaxax = MAX(xminax*10.0001D0, MAXVAL(xv(1:N), MASK = xv > 0.0D0))
                xrayspec(1) = xminax*0.99D0
                xrayspec(2) = xmaxax*1.01D0 
            ENDIF
            CALL AXSSCL ("LOG", "X") 
            CALL LOGTIC ("FULL") !"AUTO"
            CALL SETSCL(xrayspec, 2, "X") !scales the axis automatically for GRAF
        ELSE
            xrayspec(1) = xminax
            xrayspec(2) = xmaxax
            CALL AXSSCL ("LIN", "X")  
            !CALL SETSCL(xrayspec, 2, "X") !scales the axis automatically for GRAF  
        ENDIF
        
        IF (nax == barpmodeno) THEN !print the x-axis labels only at the last subgraph
            CALL SETGRF ("NAME", "NAME", "TICKS", "TICKS")
        ELSE
            CALL SETGRF ("TICKS", "NAME", "TICKS", "TICKS")
        ENDIF
        CALL LABELS ("FLOAT", "Y")
        CALL HEIGHT(Numbsize)
        CALL CHASPC(-0.07D0)
        CALL TICKS (2, "Y") !the number of ticks between axis labels (default: 2) 
        CALL TICKS (4, "X")
        IF (ABS(xstep) < 0.1D0) THEN
            CALL LABDIG (3, "X") !number of decimals xaxis
        ELSE IF (ABS(xstep) < 0.2D0) THEN
            CALL LABDIG (2, "X") !number of decimals xaxis
        ELSE IF (ABS(xstep) > 0.99D0) THEN
            CALL LABDIG (1, "X") !number of decimals xaxis
            iticks = INT(xstep)
            IF (iticks < 21) THEN
                CALL TICKS(iticks, "X")
            ELSE
                iticks = 20
                DO
                    IF (MOD(INT(xstep),iticks) == 0) THEN
                        EXIT
                    ELSE
                        iticks = iticks-1
                    ENDIF
                ENDDO
                IF (iticks < 5) THEN
                    iticks = 10
                ENDIF
                CALL TICKS(iticks, "X")
            ENDIF
            IF (ABS(xmaxax-xminax) > 5.0D0) THEN
                CALL INTAX
            ENDIF
        ELSE
            CALL LABDIG (1, "X") !number of decimals xaxis
        ENDIF
        IF (ABS(ymaxax-yminax) < 1.0D-3 .OR. ABS(ymaxax-yminax) > 9.99999D3) THEN
            CALL LABELS ("FEXP", "Y") !EXP will use powers of 10 "scientific" number format for axis labels; FEXP will use "1.0E-8" style.
        ELSE
            CALL LABELS ("FLOAT", "Y")
        ENDIF
        CALL LABDIS(Numbsize-8, "XYZ") !default is 24
        CALL TICLEN (Numbsize-8, ((Numbsize-8)*2)/3) !major and minor tick lengths
        IF (squareplot) THEN
            yminax = 0.0D0; ymaxax = 1.0D0; ystep = 0.2D0;
            CALL GRAF(xminax,xmaxax,xminax,xstep, yminax,ymaxax,yminax,ystep)
        ELSE
            yminax = 0.0D0; ymaxax = 1.0D0; ystep = 0.2D0;
            CALL GRAF(xminax,xmaxax,xminax,xstep, yminax,ymaxax,yminax,ystep)
        ENDIF
        !initialize legend:  
        CBUF = " "
        inc = ncomp(nax)
        CALL LEGINI(CBUF, MIN(inc,30), NMAXLN) !dislin can only stor 30 curve attributes currently
        inc = 0
        DO nc = 1,ncomp(nax)
            inc = inc+1
            IF (nc <= 30) THEN
                CALL LEGLIN(CBUF,TRIM(cpielegtxt(nax,nc)),inc)
            ENDIF
        ENDDO
        CALL LEGTIT(TRIM(cpietittxt(nax,nax)))
        !position legend:
        CALL GETPOS (NXA, NYA)
        CALL GETLEN (NXL, NYL, NZL)
        IF (squareplot) THEN
            CALL LEGPOS(NXA+NXL+50, NYA-NYL)
        ELSE
            CALL LEGPOS(NXA+NXL+100, NYA-NYL)  !position the legend on the right side of the axis system
        ENDIF
        IF (nax == 1) THEN !plot title
            CALL SETCLR(0)
            !check the length of the longest title line to adjust font size:
            i = MAXVAL(LEN_TRIM(ptitle(1:4)))
            IF (i > 100) THEN
                Lheight = MAX(INT(100.0D0*DFLOAT(Tsize)/DFLOAT(i)), 2)
                CALL HEIGHT(Lheight)
            ELSE
                CALL HEIGHT(Tsize)
            ENDIF
            CALL LINESP(2.6D0)
            CALL TITJUS ("LEFT")
            CALL TITLE
        ENDIF
        newgraph = .false.
    ENDIF
    !prepare the data for bars:
    ysum = 0.0D0 !clear
    y0 = 0.0D0
    DO nc = 1,ncomp(nax)
        IF (nc == 1) THEN
            ysum(nc,1:N) = compby(nax,nc,1:N) !x(component nc)
        ELSE
            ysum(nc,1:N) = ysum(nc-1,1:N) + compby(nax,nc,1:N) !x(component nc)
        ENDIF
    ENDDO !nc
    !calculate average spacing between adjacent bars:
    IF (N > 1) THEN
        !sort xv array:
        xv2 = 0.0D0
        xv2(1:N) = xv(1:N)
        CALL QsortC(xv2(1:N))
        davg = 0.0D0
        DO i = 1,N-1
            davg = davg + ABS(xv2(i)-xv2(i+1))
        ENDDO
        davg = davg/DFLOAT(N-1)
    ELSE
        davg = 36.0D0
    ENDIF
    !specify bar width:
    IF (N > 3) THEN
       yrayspec_r(1) = davg 
       ! yrayspec_r(1) = MIN(ABS(xv(2)-xv(3)),ABS(xv(N-2)-xv(N-3)))
       ! yrayspec_r(1) = MAX(yrayspec_r(1),0.5D0)
        IF (squareplot) THEN
            yrayspec_r(2) = -0.77D0*yrayspec_r(1)*(DFLOAT(xlen)/ABS(xmaxax-xminax))
            yrayspec_r(2) = roundnicely(yrayspec_r(2),5,"up")
        ELSE
            yrayspec_r(2) = -0.77D0*yrayspec_r(1)*(DFLOAT(xlen)/ABS(xmaxax-xminax))
            yrayspec_r(2) = roundnicely(yrayspec_r(2),5,"up")   
        ENDIF
        IF (logxaxis) THEN
            yrayspec_r(2) = MIN(yrayspec_r(2)*0.6D0, -24.0D0)
        ENDIF 
        CALL BARWTH(yrayspec_r(2))
    ELSE
        CALL BARWTH(-0.5D0*36.0D0)
    ENDIF
    IF (customCol) THEN !read a custom colour table (generated by ColorTableCreator of Andi Zuend):
        CALL VLTFIL(TRIM(customcolfile), "LOAD") !the Win XP Icon Color Palette plus Grayscale etc. as defined by Andi Zuend
    ELSE
        CALL SETVLT("RAIN") !the Dislin rainbow palette
    ENDIF
    !Plot the Bars:
    ALLOCATE(xray(N), yray1(N), yray2(N), STAT=allocstat)
    xray = xv(1:N)
    yray1 = 0.0D0
    yray2 = 0.0D0
    inc = ncomp(nax)
    CALL HEIGHT (Lcompsize/2-1)
    IF (logxaxis) THEN
        i = 0 !for debugging set breakpoint
        WHERE (xray(:) < xminax) 
            xray(:) = 0.1D0*xminax !make sure those values are set to a positive value, but smaller than the lower x-axis limit.
        ENDWHERE
    ENDIF
    DO nc = 1,ncomp(nax)
        CALL PENWID(0.05D0) !the linewidth (allowing also values smaller than 1)
        CALL LNCAP ("LONG") !"CUT" "ROUND" or "LONG" : CUT = cut off line caps
        CALL LNJOIN ("SHARP") !'SHARP' or 'TRUNC'	
        IF (nc == 1) THEN
            yray1 = y0(1:N)
            yray2 = ysum(nc,1:N)
            IF (ncomp(nax) == 1) THEN !only one component/bar per x-value
                DO i = 1,N
                    IF (cplotbarlabels(nax,i)) THEN !plot bar labels: the y-value
                        CALL LABELS('DELTA', 'BARS')
                        CALL LABPOS('OUTSIDE', 'BARS')
                    ELSE
                        CALL LABELS('NONE', 'BARS')
                    ENDIF
                    CALL SHDPAT(cbarshadpat(nax,nc,i)) !16=solid fill, 8=hatched pattern
                    CALL SETCLR(cbarpatterncol(nax,nc,i))
                    IF (plotbarborder) THEN
                        CALL PENWID(0.5D0) !the linewidth (allowing also values smaller than 1)
                        CALL LNCAP ("CUT") ! CUT = cut off line caps
                        CALL BARBOR(cbarbordcol(nax,nc,i)) !the color value of the bar border line, -1 means use current color
                    ENDIF
                    CALL BARS(xray(i),yray1(i),yray2(i),1)  !Plot BAR i
                ENDDO !i
            ELSE
                CALL SHDPAT(cbarshadpat(nax,nc,1)) !16=solid fill, 8=hatched pattern
                CALL SETCLR(cbarpatterncol(nax,nc,1))
                IF (plotbarborder) THEN
                    CALL PENWID(0.5D0) !the linewidth (allowing also values smaller than 1)
                    CALL LNCAP ("CUT") ! CUT = cut off line caps
                    CALL BARBOR(cbarbordcol(nax,nc,1)) !the color value of the bar border line, -1 means use current color
                ENDIF
                CALL BARS(xray,yray1,yray2,N)  !Plot BARS
            ENDIF
        ELSE !nc > 1
            yray1 = ysum(nc-1,1:N)
            yray2 = ysum(nc,1:N)
            IF (inc == 1) THEN !only one component/bar per x-value
                DO i = 1,N
                    IF (cplotbarlabels(nax,i)) THEN !plot bar labels: the y-value
                        CALL LABELS('DELTA', 'BARS')
                        CALL LABPOS('OUTSIDE', 'BARS')
                    ELSE
                        CALL LABELS('NONE', 'BARS')
                    ENDIF
                    CALL SHDPAT(cbarshadpat(nax,nc,i)) !16=solid fill, 8=hatched pattern
                    CALL SETCLR(cbarpatterncol(nax,nc,i))
                    IF (plotbarborder) THEN
                        CALL PENWID(0.5D0) !the linewidth (allowing also values smaller than 1)
                        CALL LNCAP ("CUT") ! CUT = cut off line caps
                        CALL BARBOR(cbarbordcol(nax,nc,i)) !the color value of the bar border line, -1 means use current color
                    ENDIF
                    CALL BARS(xray(i),yray1(i),yray2(i),1)  !Plot BAR i
                ENDDO !i
            ELSE !inc > 1
                CALL SHDPAT(cbarshadpat(nax,nc,1)) !16=solid fill, 8=hatched pattern
                CALL SETCLR(cbarpatterncol(nax,nc,1))
                IF (plotbarborder) THEN
                    CALL PENWID(0.5D0) !the linewidth (allowing also values smaller than 1)
                    CALL LNCAP ("CUT") ! CUT = cut off line caps
                    CALL BARBOR(cbarbordcol(nax,nc,1)) !the color value of the bar border line, -1 means use current color
                ENDIF
                CALL BARS(xray,yray1,yray2,N)  !Plot BARS
            ENDIF
        ENDIF !nc
    ENDDO !nc
    !Plot Bar-Legends:
    CALL PENWID(1.0D0) !the linewidth (allowing also values smaller than 1)
    IF (ncomp(nax) > 8 .AND. barpmodeno > 1) THEN !scale legend entry sizes
        Lheight = MAX(INT(8.52D0*DFLOAT(Lcompsize)/DFLOAT(ncomp(nax))), 2)
        CALL HEIGHT (Lheight)  
    ELSE
        CALL HEIGHT (Lcompsize)
    ENDIF
    CALL LINESP (1.8D0)
    CALL SETCLR(0)
    CALL FRAME(0) !no frame
    CALL LEGEND(CBUF,7)
    CALL SOLID
    CALL FRAME(1)
    CALL BOX2D !plots a box around legend key
    CALL ENDGRF !set back level to 1
    DEALLOCATE(xray, yray1, yray2, STAT=allocstat)
ENDDO !nax

CALL SETCLR(0) !black
CALL SOLID
CALL FRAME(1)
CALL BOX2D !plots a box around the axis system
CALL ENDGRF !set back level to 1 for plotting 3rd axis system!
  
!------------------------------------------------------------------------------------------------
!plot the message box (if there was any line of text defined):
IF (messagebox) THEN
    !set the (initial) plot coordinates and properties of the message box:
    IF (messboxpos_x == 0 .AND. messboxpos_y == 0) THEN !the position was not set yet, so use default values
        nxl = 1200
        nyl = 120
    ELSE
        nxl = messboxpos_x
        nyl = messboxpos_y
    ENDIF
    CALL HEIGHT (Lsize)
    CALL TXTJUS ("LEFT")
    CALL FRMESS (0) !the thickness of the frame around the box
    !write out message box header horizontal line:
    N = MAXVAL(LEN_TRIM(messagetxt(1:txtl)))
    IF (N > 0) THEN
        i = MIN(INT(N*5.0D0), 256)
        hline = REPEAT("-", i)
        CALL CHASPC(-0.4D0)
        CALL MESSAG (TRIM(hline), nxl, nyl)
        nyl = nyl + INT(Lsize*2.3D0)
        CALL CHASPC(-0.01D0)
        !write out the saved textlines:
        DO i=1,txtl 
            CALL MESSAG ("   "//TRIM(messagetxt(i)), nxl, nyl)
            nyl = nyl + INT(Lsize*2.3D0)
        ENDDO
        !write out message box footer horizontal line:
        CALL CHASPC(-0.4D0)
        CALL MESSAG (TRIM(hline), nxl, nyl)
    ENDIF
ENDIF
!------------------------------------------------------------------------------------------------

!END DISLIN plotting:  
CALL WINMOD("NONE") !for graphic output on XWIN only.:
CALL TEXMOD("OFF")
CALL DISFIN !end DISLIN

!reset some switches:
newplot = .true. !set it back for next plot page
newplot2 = .true.
messagebox = .false.
smallsymbs = .false.
bigsymbs = .false.
specialxrange = .false.
specialyrange = .false.
specialyrange2 = .false.
logyaxis = .false. !set the logarithmic y-axis to false (so use linear).
logxaxis = .false.
ptitle = ""
customCol = .false.
bargraphplot = .false.
legplot = .true.
bargylimits = .false.
plotbarborder = .false.
    
END SUBROUTINE PlotBarGraphs
!=================================================================================================================\


!=================================================================================================================
!SUBROUTINE to read and initialize a color table for use of up to NTCP plot colors 
!stored in colneutral and colelectrol
SUBROUTINE defplotcolorlists()

IMPLICIT NONE

INTEGER(4) :: i, n
CHARACTER(LEN=20) :: colfile, dummy
LOGICAL(4) :: exists, fileend
!................

!check if reading the data from file is necessary at this call (if done previously in this session, it can be skipped):
IF (ALL(colneutral(1:4) > -1) .AND. ALL(colneutral(1:4) < 256) .AND. colneutral(1) /= colneutral(2) .AND. colelectrol(1) == 11) THEN
    !data already read... so continue  
ELSE !read data
    !define electrolyte color list (for now not read from a file):
    colelectrol(1:6) = [ 11, 248, 222, 253, 12, 240 ] 

    !read neutral component colors from a text-file:
    colneutral = 250 !initialize 
    colfile = "plotcolortable1.dat"
    INQUIRE(FILE=TRIM(colfile), EXIST=exists)
    IF (exists) THEN
        OPEN (UNIT=22, FILE=TRIM(colfile), STATUS="OLD")
        READ(22,*) dummy !read first line
        READ(22,*) !read empty line
        READ(22,*) dummy !read line 3
        fileend = .false.
        DO i = 1,NTCP  !read in the data
            fileend = EOF(22)
            IF (.NOT. fileend) THEN
                READ(22,*) n, colneutral(i) !read data using implied do-loop for the values of tLLEwtf_d(1:20,i)
            ELSE
                EXIT
            ENDIF
        ENDDO
        CLOSE(22)
    ELSE
        !$OMP CRITICAL
        WRITE(*,*) ""
        WRITE(*,*) "ERROR in Mod_createplot.f90: file listing plot colors not found!"
        WRITE(*,*) ""
        READ(*,*) !wait for user action
        !$OMP END CRITICAL
        colneutral(1:NTCP) = 10 
    ENDIF
ENDIF

END SUBROUTINE defplotcolorlists
!=================================================================================================================

!=================================================================================================================
!SUBROUTINE to retrieve a colour scale number (colno) from the custom set colour bar 
!(i.e. for coloured data points on a virtual z-axis as part of a 2-D plot)
SUBROUTINE SetCustomColPalette(colpalette)

IMPLICIT NONE
!interface variables:
CHARACTER(LEN=*),INTENT(IN) :: colpalette
!...............................

SELECT CASE(TRIM(colpalette))
CASE('WinXP')
    customcolfile = 'RGBColTableRealWinXP8BitAZ.dat' 
CASE("Viridis")
    customcolfile = 'RGBColTabViridisPurpleBlueGreenYellow.dat'
CASE('BlackPurpleRedYellow')
    customcolfile = 'RGBColTabBlackPurpleRedYellow.dat'
CASE('RedBlue')
    customcolfile = 'RedBlueColTableReal_dislin.dat'
END SELECT
customCol = .true.

END SUBROUTINE SetCustomColPalette
!=================================================================================================================

!=================================================================================================================
!SUBROUTINE to retrieve a colour scale number (colno) from the custom set colour bar 
!(i.e. for coloured data points on a virtual z-axis as part of a 2-D plot)
PURE ELEMENTAL SUBROUTINE defzcolval(zval, colno, colbarlowlim, colbaruplim)

IMPLICIT NONE
!interface variables:
REAL(8),INTENT(IN) :: zval          !the "z-axis" value for which a colour number should be determined
INTEGER(4),INTENT(OUT) :: colno     !the output colour index on 1 to 254 scale; i.e. the index of the colour number from a (custom) colour table used;
REAL(8),INTENT(IN) :: colbarlowlim  !the lowest value mapped to the colour bar, i.e. the z value corresponding to colour 1
REAL(8),INTENT(IN) :: colbaruplim   !the highest value mapped to the colour bar, i.e. the z value corresponding to colour 254
!local variables:
REAL(8) :: colslope, intercept
!...............................

!linear color scale between colors 1 and 254
colslope = 253.0D0/(colbaruplim - colbarlowlim)
intercept = 1.0D0 -colslope*colbarlowlim
colno = INT(intercept + colslope *MAX(MIN(zval, colbaruplim), colbarlowlim) )     

END SUBROUTINE defzcolval
!=================================================================================================================

!=================================================================================================================
!SUBROUTINE to define plotting properties of secondary z-axis colour bars (customizable)
SUBROUTINE showzcolbar(colbarlowlim, colbaruplim, zlabel)

IMPLICIT NONE
!interface variables:
REAL(8),INTENT(IN) :: colbarlowlim      !the lowest value mapped to the colour bar, i.e. the z value corresponding to colour 1
REAL(8),INTENT(IN) :: colbaruplim       !the highest value mapped to the colour bar, i.e. the z value corresponding to colour 254
CHARACTER(LEN=75),INTENT(IN) :: zlabel  !z axis label text string
!local variables:
INTEGER(4) :: i
INTEGER(4),DIMENSION(254) :: zcolno
REAL(8) :: zinterval
REAL(8),DIMENSION(254) :: zv
!...............................

!set private module variables to be used with ZAXIS or ZAXLG in 'Plotnow':
zcolbar0 = colbarlowlim
zcolbar254 = colbaruplim
zcolbarstep = (colbaruplim - colbarlowlim)/10.0D0            !define the stepping between tick marks
zaxislabel = zlabel
showcustomcolbar = .true. 
!define array of colour values based on custom function defzcolval(...); 
!to be shown as coloured rectangles as a custom-made colour bar:
zinterval = (colbaruplim - colbarlowlim)/253.0D0
DO i = 1,254
    zv(i) = colbarlowlim + REAL(i-1, KIND=8)*zinterval
    CALL defzcolval(zv(i), zcolno(i), colbarlowlim, colbaruplim)
ENDDO
customzcolno(1:254) = zcolno(1:254)

END SUBROUTINE showzcolbar
!=================================================================================================================

END MODULE Mod_createplot


!=================================================================================================================
!** user defined SUBROUTINE to set some label position options when plotting Pie-Charts **

SUBROUTINE MyPieSettings(ISEG,XDAT,XP,NRAD,NOFF,ANG,NVX,NVY,IDRW,IANN)

USE Mod_createplot, ONLY : pielabshift

IMPLICIT NONE
!..
INTEGER(4),INTENT(IN) :: ISEG
INTEGER(4),INTENT(INOUT) :: NRAD,NOFF,NVX,NVY,IDRW,IANN
!..
REAL(8),INTENT(IN) :: XDAT,XP
REAL(8),INTENT(INOUT) :: ANG
!..........................................................

IF (pielabshift(ISEG)) THEN
    IF (ISEG > 1) THEN
        IF (MOD(ISEG,2) == 0 .AND. (.NOT. pielabshift(ISEG-1))) THEN
            NVY = -40  !-40
        ELSE IF (pielabshift(ISEG-1) .AND. pielabshift(ISEG+1)) THEN
            NVY = 10
        ELSE
            NVY = 30  !-40
        ENDIF
    ELSE
        NVY = 30  !-40
    ENDIF
ENDIF

END SUBROUTINE MyPieSettings
!=================================================================================================================