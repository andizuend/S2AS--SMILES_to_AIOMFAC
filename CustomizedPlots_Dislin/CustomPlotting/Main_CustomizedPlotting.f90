!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Program to plot custom made 2D x-y scatter plots, curves,                          *
!*   bar graphs, and pie charts using the DISLIN graphics library.                      *
!*   Here specifically for plotting the 2D lumping framework outputs and grid lines.    *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2012-06-08                                                      *
!*   -> latest changes: 2025-07-06                                                      *
!*                                                                                      *
!****************************************************************************************
program CustomizedPlotting_2D_framework

use Mod_kinds, only : wp
use Mod_PlotData, only : ReadDataFiles
use ModCreateplot

implicit none
!...
integer :: i, isymb, k, j, jcent, kset, nsets, icount, inrows, allocstat, info, icol, ny, nx, &
    & clustN, clustN_adj, fileStartNo, fileEndNo, maxClust
integer,parameter :: NCinplot = 30
integer,parameter :: maxinf = 50                                            !set maximum number of different input files (containing data for different curves...)
integer,parameter :: maxrows = int(1E5)                                     !set maximum number of data rows per dataset
integer,parameter :: ndim = 8                                               !number of dimensions (columns) of input file CompPropData array
integer,parameter :: ndl = 6                                                !number of columns of the input file LumpConc
integer,dimension(NCinplot) :: colset1, barcolset
integer,dimension(6) :: colset2
integer,dimension(2*NCinplot,maxrows) :: barpatcol, shadpat, barbordercol
integer,dimension(:),allocatable :: cpv, c2pv, dcol, krows, compCluster, iz, ia
integer,dimension(:),allocatable :: clustSurrogateIndx, clustPop            !dimension covers all clusters
integer,dimension(:,:),allocatable :: key
!...
character(len=5) :: lumpResChar, cchar
character(len=75) :: fileoutname, xlabel, lylabel, rylabel, zlabel, legtxt, ichar, xaxchar, yaxchar, colchar
character(len=200) :: pvchar
character(len=10),dimension(:),allocatable :: filenumber
character(len=30),dimension(:),allocatable :: YaxisChoice
character(len=132),dimension(4) :: ptitle
character(len=75),dimension(:),allocatable :: legendCompPropIn, legendLumpedConcIn
character(len=75),dimension(:),allocatable :: legendCompProp, legendLumpedConc
!...
real(wp) :: axextend, SD, lylim, uylim, uxlim, lxlim, cthickn, colbarlowlim, colbaruplim, sumv
real(wp),dimension(:),allocatable :: xv, yv, zv, wv, xvals, yvals, xgrid, ygrid, x2vals, y2vals
real(wp),dimension(:,:),allocatable :: xyerror0
real(wp),dimension(:,:,:),allocatable :: CompPropDataIn,LumpConcDataIn      !data files, structure: [data row | x1-value | y1-value | y2-value | y3-value |... ]
real(wp),dimension(:,:,:),allocatable :: CompPropData, LumpConcData
real(wp),dimension(:,:),allocatable :: clustCenter 
!...
logical :: dexists, isfullsyst, isKmeans, isSurr, yaxisIsLog, high_volat_surr_present, plot_high_volat_cluster
logical,dimension(:),allocatable :: cPropDataExists                      !reports back whether the input files and read data exist
!..................................


!#### Parameters  ###################################################################################

fileStartNo = 1260                  !ID of first file (must be the FullSystem file) of a sequence to be read and processed 
fileEndNo   = 1264                  !ID of last file of a sequence to be read and processed
lumpResChar = '10x05'               ! '10x05'  '06x03'  '08x04'   !the text string for the lumping grid resolution data files to be selected from Output_lumping

plot_high_volat_cluster = .true.    !set to plot the special high-volatility cluster population and associated surrogate in the k-means cluster ID plot (if present) 

!####################################################################################################


write(*,*) ""
write(*,*) "Customized Plotting:  reading data files... "
write(*,*) ""
!read data from files for the plots to be made; change file names etc. in corresponding ReadDataFiles Subroutine:
k = fileEndNo - fileStartNo + 1
allocate( CompPropDataIn(maxrows,ndim,maxinf), LumpConcDataIn(maxrows,ndl,maxinf), YaxisChoice(k), cPropDataExists(k), filenumber(k), &
    & legendCompPropIn(k), legendLumpedConcIn(k), krows(k), STAT=allocstat )
!-----
call ReadDataFiles(fileStartNo, lumpResChar, nsets, inrows, CompPropDataIn, LumpConcDataIn, xgrid, ygrid, &
                    & legendCompPropIn, legendLumpedConcIn, cPropDataExists, filenumber, YaxisChoice, &
                    & clustN, clustCenter, clustSurrogateIndx, clustPop, compCluster)
!-----
if (nsets > 0 .and. inrows > 0) then
    dexists = .true.
else
    dexists = .false.
endif
!check read in data for existence:
if (.not. dexists) then
    write(*,'(A,/)') "Plot not made due to missing data!"
    write(*,'(A,/)') "Check your input files and their location."
    read(*,*)
    stop !end the program
else
    allocate(cpv(inrows), xv(inrows), yv(inrows), zv(inrows), wv(inrows), xyerror0(inrows,4), xvals(inrows), yvals(inrows), &
        & CompPropData(inrows,ndim,nsets), LumpConcData(inrows,ndim-1,nsets), legendCompProp(k), legendLumpedConc(k), &
        & dcol(inrows), STAT=allocstat)
    CompPropData = 0.0_wp
    LumpConcData = 0.0_wp
    legendCompProp = ""
    legendLumpedConc = ""
    write(*,'(A,/)') "Processing data for plots. This may take a while... "
    write(*,'(A,/)') "(especially in case of thousands of data points)"
endif

if (any(cPropDataExists(1:nsets))) then
    do k = 1,nsets
        do i = 1,ndim
            do j = 1,inrows
                if (CompPropDataIn(j,1,k) <= rdefault) then
                    exit !leave inner loop as no more data for this file
                endif
                CompPropData(j,i,k) = CompPropDataIn(j,i,k)
                if (i <= ndl) then
                    LumpConcData(j,i,k) = LumpConcDataIn(j,i,k)
                endif
            enddo
        enddo
        legendCompProp(k) = legendCompPropIn(k)
        legendLumpedConc(k) = legendLumpedConcIn(k)
        krows(k) = min(inrows, j-1)
    enddo
endif
deallocate( CompPropDataIn, LumpConcDataIn, legendCompPropIn, legendLumpedConcIn )

!initialize:
call SetCreateplotArrays(NCinplot, .true., '../Input_data/', '../Output_Plots/', '../ColourPalettes/')  !initialize ModCreatePlot parameters
xyerror0 = 0.0_wp

!define a color set for up to NCinplot components for use in curves:
colset1(1:NCinplot) = colneutral(1:NCinplot)
colset2(1:6) = colelectrol(1:6)

!define a color set for up to NCinplot bars / bar parts for bargraphs:
barcolset = 35                                      !light grey as default color (e.g. for use of bar data visualisation outside of component specific color no.)
barcolset(1) = colset1(1) + 2                       !use a lighter blue color for bar graphs of component 1
barcolset(2:NCinplot) = colset1(2:NCinplot)
barpatcol = 50                                      !default value
barpatcol(1:NCinplot,1:inrows) = spread(barcolset(1:NCinplot),2,inrows); shadpat = 16; barbordercol = 0;

dcol = 9
i = min(NCinplot, nsets)
dcol(1:i) = colset1(1:i)    !default colors
allocate( key(maxval(krows(:)),nsets) ) 
key = -1  !initialize

!determine whether a special high-volatility surrogate exists:
nx = size(xgrid)
if ( any(CompPropData(1:krows(1),5,1) > 1.001_wp*10.0_wp**xgrid(nx)) ) then
    high_volat_surr_present = .true.
else
    high_volat_surr_present = .false.
endif

!!=========================================================================================================================
!!= Plot of data from read input files
!!=========================================================================================================================
k = 0
do kset = 1,2*nsets  !loop over all files / input data sets
    !!=========================================================================================================================
    !!= y-axis choice (e.g. activity coefficient ratio for each organic with water and hexanediol) vs. pure-component saturation pressure or concentration
    !!=========================================================================================================================
    k = k + 1
    if (kset == nsets + 1) then !reset k
        k = 1
    endif
    if (cPropDataExists(k)) then
        !initial plot values:
        newplot = .true.
        newplot2 = .true.
        messagebox = .false.
        messboxpos_x = 0
        messboxpos_y = 0
        uxlim = -1.0E80_wp
        lxlim = 1.0E80_wp
        uylim = -1.0E80_wp
        lylim = 1.0E80_wp
        specialxrange = .false.
        specialyrange = .false.
        bargylimits = .false.
        logxaxis = .false.
        plotbarborder = .false.
        max_y = 1.0E15_wp
        min_y = 0.0_wp !-1.0E15_wp
        smallsymbs = .false.
        bigsymbs = .false.
        call SetCustomColPalette('Viridis')     !load custom-made colour palette in Mod_createplot
        logyaxis = .false.
        logxaxis = .false.
        legplot = .true.
        nogridlines = .false.
        metaff = "pdf"  !wmf  pdf  eps  png  gif
        ptitle = ""
        cthickn = 3.0_wp !6.0_wp
        !!DefinePlotData(xvals, yvals, xyerror(the +- error bar size for x and y axis [dimensions 1:4 are x-err neg, x-err pos, y-err neg, y-err pos]), waxis(y-axis 1 or 2?), 
        !               color, ltype(1=only line, 2=symb and line, > 2 means symbols and no line), symbtype(-1=no symbol), lstyle(1=solid, 2=dashm, 3=dash, 4=Mydot, 5=dot, 6=dash-dotted),legtxt,cthickn)
        dcol = 0
        xv = 0.0_wp
        yv = 0.0_wp
        wv = 0.0_wp
        zv = 0.0_wp
        cpv = 0
        !determine whether this is the full system or a lumped version:
        isfullsyst = .false.
        isKmeans = .false.
        if (index(legendLumpedConc(k), "_FullSystem") > 0) then     !true if not lumped
            isfullsyst = .true.
        else if (index(legendLumpedConc(k), "_Kmeans") > 0) then
            isKmeans = .true.
        endif
        !set x-y values:
        logxaxis = .true.
        xv(1:krows(k)) = CompPropData(1:krows(k),5,k)           !psat [Pa]
        wv(1:krows(k)) = LumpConcData(1:krows(k),2,k)           !the component mass conc. (after potential lumping)
        cpv(1:krows(k)) = int(CompPropData(1:krows(k),1,k))     !store component number as sorting will be applied later
        select case(YaxisChoice(k))
        case('log10(act.coeff_ratio)')
            logyaxis = .true.
            yv(1:krows(k)) = CompPropData(1:krows(k),6,k)
            yaxchar = "LgActcoeffR"
            lylabel = "$ \gamma_{j,{\rm hex}}/ \gamma_{j,{\rm H_2O}} $"
        case('O:C_ratio')
            logyaxis = .false.
            yv(1:krows(k)) = CompPropData(1:krows(k),2,k)
            yaxchar = "OtoC"
            lylabel = "O:C ratio"
        case('log10(O:C_ratio)')
            logyaxis = .true.
            yv(1:krows(k)) = max(CompPropData(1:krows(k),2,k), 1.0E-2_wp)   !1.0E-2 is used for cases where the O:C ratio is lower to avoid undefined log values
            yaxchar = "LgOtoC"
            lylabel = "$\log_{10}{\rm (O:C)}$"
        case('meanOS_C')
            logyaxis = .false.
            yv(1:krows(k)) = CompPropData(1:krows(k),4,k)
            yaxchar = "MeanOSC"
            lylabel = "$\overline{\rm OS}_{\rm C}$"
        end select
        !select colour axis quantity:
        if (kset <= nsets) then
            call SetCustomColPalette('Viridis')
            zv(1:krows(k)) = CompPropData(1:krows(k),2,k)
            colchar = "OtoC_ratio"
            colbaruplim = max( 2.0_wp, min( 3.0_wp, maxval(zv(1:krows(k))) ) )
            colbarlowlim = 0.0_wp
        else
            call SetCustomColPalette('BlackPurpleRedYellow')    !set a custom colour palette for these plots
            !call SetCustomColPalette('BluePurpleOrangeYellow')
            colchar = "Log10TotMassConc"                        !"TotMassConc"
            sumv = sum(LumpConcData(1:krows(k),2,k))
            colbaruplim = real(floor(log10(0.9_wp*sumv)), kind=wp)
            colbarlowlim = colbaruplim - 6.0_wp                  !for log-scaling, based on upper limit; 0.0_wp
            where(LumpConcData(1:krows(k),2,k) > 0.0_wp)
                zv(1:krows(k)) = log10( LumpConcData(1:krows(k),2,k) )
            elsewhere
                zv(1:krows(k)) = colbarlowlim - 1.0_wp           !make lower than colbarlowlim so it won't be plotted
            endwhere
        endif
        !--- sorting of the plot data arrays such that mass conc. as the sorting condition is increasing with increasing array index:
        if (key(1,k) < 0) then
            key(1:krows(k),k) = [(i, i = 1,krows(k))]                               !the original array index positions before sorting;
            call dlasrt2('I', krows(k), wv(1:krows(k)), key(1:krows(k),k), info)    !use sorting sub with sort order output in key array;
        else
            wv(1:krows(k)) = wv( key(1:krows(k),k) )     
        endif
        !now use the key to sort all the associated arrays:
        xv(1:krows(k)) = xv( key(1:krows(k),k) ) 
        yv(1:krows(k)) = yv( key(1:krows(k),k) )
        zv(1:krows(k)) = zv( key(1:krows(k),k) )
        cpv(1:krows(k)) = cpv( key(1:krows(k),k) ) 
        !--- end of sorting
        yaxisIsLog = logyaxis
        icount = 0
        xvals = 0.0_wp
        yvals = 0.0_wp
        do i = 1,krows(k)     !loop over data rows (component points)
            if (wv(i) > 0.0_wp) then !plot only points of non-zero mass concentration:
                logyaxis = yaxisIsLog
                icount = icount + 1
                xvals(1) = xv(i)
                yvals(1) = yv(i)            
                write(legtxt, '(I0)') cpv(i)    !transfer component number
                write(ichar, '(ES10.3)') zv(i)  !transfer z-axis value
                legtxt = "cp"//trim(legtxt)//", "//trim(colchar)//" = "//trim(ichar)
                if (icount == 28) then          !only max. 30 legend entries are shown
                    legtxt = "plus additional points ..."
                else if (icount > 28) then
                    legtxt = ""
                endif
                !if (krows(k) < 1000) then
                !    bigDsymbs = .true.
                !endif
                !point colors:
                call defzcolval(zv(i), dcol(i), colbarlowlim, colbaruplim)  !map the z value to the colour bar index number dcol
                !--
                call DefinePlotData(xvals, yvals, xyerror0, 1, dcol(i), 3, 21, 5, legtxt, cthickn)
                !--
                lxlim = min( lxlim, xv(i) )
                uxlim = max( uxlim, xv(i) )
                lylim = min( lylim, yv(i) )
                uylim = max( uylim, yv(i) )
            endif
        enddo !i

        !AZ added: plot grid lines of chosen lumping grid:
        if (.not. isKmeans) then                                !only plot grid lines for non-Kmeans cases
            cthickn = 3.0_wp
            icol = -67                                          !set to dark gray (- to use WinXP palette)
            xv = 0.0_wp
            yv = 0.0_wp
            !consider grid line data conversion:
            if (kset == 1) then
                if (logxaxis) then
                    xgrid = 10.0_wp**xgrid
                endif
                if (yaxisIsLog) then
                    ygrid = 10.0_wp**ygrid
                endif
            endif
            nx = size(xgrid)
            ny = size(ygrid)
            !plot x-axis grid lines:
            do j = 1,nx
                if (j == 1) then
                    legtxt = "x-lumping-grid line"
                else
                    legtxt = ""
                endif
                xv(1:2) = xgrid(j)
                yv(1) = ygrid(1)
                yv(2) = ygrid(ny)
                call DefinePlotData(xv(1:2), yv(1:2), xyerror0, 1, icol, 1, -1, 5, legtxt, cthickn)
            enddo
            !plot y-axis grid lines:
            do j = 1,ny
                if (j == 1) then
                    legtxt = "y-lumping-grid line"
                else
                    legtxt = ""
                endif
                xv(1) = xgrid(1)    !start x value
                xv(2) = xgrid(nx)   !end x value
                yv(1:2) = ygrid(j)  !current y-axis grid line value
                call DefinePlotData(xv(1:2), yv(1:2), xyerror0, 1, icol, 1, -1, 5, legtxt, cthickn)
            enddo
            lxlim = min( lxlim, min(xgrid(1), xgrid(nx)) )
            uxlim = max( uxlim, max(xgrid(1), xgrid(nx)) )
            lylim = min( lylim, min(ygrid(1), ygrid(ny)) )
            uylim = max( uylim, max(ygrid(1), ygrid(ny)) )
        endif
        
        ! °°now plot with following parameters °°
        xaxchar = "LgPsat"
        fileoutname = trim(yaxchar)//"_vs_"//trim(xaxchar)//"_"//trim(filenumber(k))//"_res"//trim(lumpResChar)//"_col"//trim(colchar)
        ptitle(1) = trim(fileoutname)
        ptitle(2) = "Activity coefficient ratio vs. pure-component saturation vapor pressure"  !"text describing plot..."
        if (.not. isKmeans) then
            !ptitle(3) = "Input file and method: "//trim(legendCompProp(k))//",  "//"grid resolution: "//trim(lumpResChar)
            if (high_volat_surr_present) then
                write(cchar,'(I0.2)') krows(k) - 1
                ptitle(3) = "Input file and method: "//trim(legendCompProp(k))//",  "//"no. of surrogates: "//trim(cchar)//" (+1 high-volat.)"
            else
                write(cchar,'(I0.2)') krows(k)
                ptitle(3) = "Input file and method: "//trim(legendCompProp(k))//",  "//"no. of surrogates: "//trim(cchar)
            endif
        else
            if (high_volat_surr_present) then
                write(cchar,'(I0.2)') clustN - 1
                ptitle(3) = "Input file and method: "//trim(legendCompProp(k))//",  "//"no. of clusters: "//trim(cchar)//" (+1 high-volat.)"
            else
                write(cchar,'(I0.2)') clustN
                ptitle(3) = "Input file and method: "//trim(legendCompProp(k))//",  "//"no. of clusters: "//trim(cchar)
            endif
        endif
        !lylabel is now set above in the select case block; lylabel = "$ \gamma_{j,{\rm hex}}/ \gamma_{j,{\rm H_2O}} $"
        xlabel = "$p^{\rm o,sat}_j$ / [Pa]"  !or ? "$C^{\rm o}$ / $[\mu {\rm g}/{\rm m}^3$]"
        rylabel = ""
        if (trim(colchar) == "OtoC_ratio") then
            zlabel = "O:C ratio"
        else if (trim(colchar) == "TotMassConc") then
            zlabel = "mass conc. / [$\mu {\rm g}\;{\rm m}^{-3}$]"
        else if (trim(colchar) == "Log10TotMassConc") then
            zlabel = "log$_{\rm 10}$(mass conc. / [$\mu {\rm g}\;{\rm m}^{-3}$])"
        else
            zlabel = " !! undefined !!"
        endif
        SD = 0.0_wp
        squareplot = .false.
        smallsymbs = .true.
        fixminxax = .false.
        messagebox = .false.
        nogridlines = .true.    !to suppress plotting of background gridlines when already the lumping grid is plotted
        call showzcolbar(colbarlowlim, colbaruplim, zlabel, .false.)   !to set z-axis colour bar properties and show colour bar on plot
        specialxrange = .true.  !to plot with a defined x-axis range
        specialyrange = .true.  !to plot with a defined y-axis range
        if (k == 1) then  !set axis limits and use the same for k > 1
            if (logxaxis) then
                lowerxlimit = 10.0_wp**floor(log10(lxlim - 0.05_wp*lxlim))
                if (lxlim - 0.05_wp*lxlim > 4.0_wp*lowerxlimit) then
                    lowerxlimit = lowerxlimit*3.0_wp
                endif
                upperxlimit = 10.0_wp**ceiling(log10(uxlim + 0.05_wp*uxlim))
                if (uxlim + 0.05_wp*uxlim < 0.3_wp*upperxlimit) then
                    upperxlimit = 0.4_wp*upperxlimit  
                endif
            else
                axextend = max(abs(0.05_wp*lxlim), abs(0.05_wp*uxlim))
                lowerxlimit = lxlim -axextend
                upperxlimit = uxlim +axextend
            endif
            if (yaxisIsLog) then
                lowerylimit = 10.0_wp**floor(log10(lylim - 0.05_wp*lylim))
                if (lylim > 4.0_wp*lowerylimit) then
                    lowerylimit = lowerylimit*3.0_wp
                endif
                upperylimit = 10.0_wp**ceiling(log10(uylim + 0.05_wp*uylim))
                if (uylim < 0.3_wp*upperylimit) then
                    upperylimit = 0.4_wp*upperylimit  
                endif
            else
                axextend = max(abs(0.05_wp*lylim), abs(0.05_wp*uylim))
                lowerylimit = lylim -axextend
                upperylimit = uylim +axextend
            endif
        endif
        call PlotNow(ptitle, xlabel, lylabel, rylabel, fileoutname, SD)
    endif    
    !========================================================================================================================
enddo !kset

!========================================================================================================================
do k = 1,nsets
    if (index(legendLumpedConc(k), "_FullSystem") > 0 .and. any(index(legendLumpedConc(:), "_Kmeans") > 0) ) then
        !!=========================================================================================================================
        !!= Plot cluster members and centers for k-means data
        !!=========================================================================================================================
        if (cPropDataExists(k)) then
            !!DefinePlotData(xvals, yvals, xyerror(the +- error bar size for x and y axis [dimensions 1:4 are x-err neg, x-err pos, y-err neg, y-err pos]), waxis(y-axis 1 or 2?),
            !               color, ltype(1=only line, 2=symb and line, > 2 means symbols and no line), symbtype(-1=no symbol), lstyle(1=solid, 2=dashm, 3=dash, 4=Mydot, 5=dot, 6=dash-dotted),legtxt,cthickn)
            allocate( c2pv(clustN), iz(inrows), ia(clustN), x2vals(clustN), y2vals(clustN) )  
            dcol = 0
            x2vals = 0.0_wp
            y2vals = 0.0_wp
            xv = 0.0_wp
            yv = 0.0_wp
            wv = 0.0_wp
            zv = 0.0_wp
            cpv = 0
            c2pv = 0
            iz = 0
            ia = 0
            !set x-y values:
            logxaxis = .true.
            xv(1:krows(k)) = CompPropData(1:krows(k),5,k)
            wv(1:krows(k)) = LumpConcData(1:krows(k),2,k)           !the component mass conc. (after potential lumping)
            cpv(1:krows(k)) = int(CompPropData(1:krows(k),1,k))     !store component number as sorting will be applied later
            select case(YaxisChoice(k))
            case('log10(act.coeff_ratio)')
                logyaxis = .true.
                yv(1:krows(k)) = CompPropData(1:krows(k),6,k)
                yaxchar = "LgActcoeffR"
                lylabel = "$ \gamma_{j,{\rm hex}}/ \gamma_{j,{\rm H_2O}} $"
            case('O:C_ratio')
                logyaxis = .false.
                yv(1:krows(k)) = CompPropData(1:krows(k),2,k)
                yaxchar = "OtoC"
                lylabel = "O:C ratio"
            case('log10(O:C_ratio)')
                logyaxis = .true.
                yv(1:krows(k)) = max(CompPropData(1:krows(k),2,k), 1.0E-2_wp)   !1.0E-2 is used for cases where the O:C ratio is lower to avoid undefined log values
                yaxchar = "LgOtoC"
                lylabel = "$\log_{10}{\rm (O:C)}$"
            case('meanOS_C')
                logyaxis = .false.
                yv(1:krows(k)) = CompPropData(1:krows(k),4,k)
                yaxchar = "MeanOSC"
                lylabel = "$\overline{\rm OS}_{\rm C}$"
            end select
            !select colour axis quantity:
            call SetCustomColPalette('BlueGreenWhiteOrangePurple')  !divergent colour palette	
            !call SetCustomColPalette('BluePurpleOrangeYellow')     !set a custom colour palette for these plots
            colchar = "KmeansClusterID"
            colbaruplim = clustN
            colbarlowlim = 0.0_wp
            iz(1:krows(k)) = compCluster(:)                         !here the cluster no. to which a component belongs;
            maxClust = maxval(compCluster(:))
            !now use the key to sort all the associated arrays:
            xv(1:krows(k)) = xv( key(1:krows(k),k) )
            yv(1:krows(k)) = yv( key(1:krows(k),k) )
            wv(1:krows(k)) = wv( key(1:krows(k),k) )
            iz(1:krows(k)) = iz( key(1:krows(k),k) )
            cpv(1:krows(k)) = cpv( key(1:krows(k),k) )
            !--- end of sorting
            jcent = 0
            icount = 0
            logxaxis = .true.
            do i = 1,krows(k)                               !loop over data rows (component points)
                if (wv(i) > 0.0_wp) then                     !plot only points of non-zero mass concentration:
                    logyaxis = yaxisIsLog
                    xvals(1) = xv(i)
                    yvals(1) = yv(i)
                    write(legtxt, '(I0)') cpv(i)            !transfer component number
                    write(ichar, '(I0.2)') iz(i)            !transfer z-axis value
                    legtxt = "cp"//trim(legtxt)//", "//trim(colchar)//" = "//trim(ichar)
                    if (cpv(i) == clustSurrogateIndx(iz(i)) ) then  !detected surrogate component of a cluster
                        jcent = jcent + 1
                        isSurr = .true.
                        x2vals(jcent) = xv(i)               !store cluster surrogate component properties for plotting later
                        y2vals(jcent) = yv(i)
                        ia(jcent) = iz(i)
                        c2pv(jcent) = cpv(i)
                    else
                        icount = icount + 1
                        isymb = 21
                        isSurr = .false.
                    endif
                    if (icount == 28) then                  !only max. 30 legend entries are shown
                        legtxt = "plus additional points ..."
                    else if (icount > 28) then
                        legtxt = ""
                    endif
                    !if (krows(k) < 1000) then
                    !    bigDsymbs = .true.
                    !endif
                    cthickn = 3.0_wp
                    if (.not. isSurr) then
                        if (high_volat_surr_present .and. (.not. plot_high_volat_cluster) .and. iz(i) == maxClust) then
                            cycle                           !don't plot this point and continue
                        else
                            call defzcolval(real(iz(i),kind=wp), dcol(i), colbarlowlim, colbaruplim)     !map the O:C value to the colour bar index number dcol
                            if (iz(i) == maxClust) dcol(i) = -21    !for the highly-volatile species cluster use a special WinXP color value
                            call DefinePlotData(xvals, yvals, xyerror0, 1, dcol(i), 3, isymb, 1, legtxt, cthickn)
                        endif
                    endif
                endif
            enddo !i
            !now loop over all cluster surrogate components and plot them after the other points, so they are on top and remain visible:
            xvals = 0.0_wp
            yvals = 0.0_wp
            do i = 1,clustN
                logyaxis = yaxisIsLog
                icount = icount+1
                xvals(1) = x2vals(i)
                yvals(1) = y2vals(i)
                write(legtxt, '(I0)') c2pv(i)           !transfer component number
                write(ichar, '(I0.2)') ia(i)            !transfer z-axis value
                legtxt = "cp"//trim(legtxt)//", "//trim(colchar)//" = "//trim(ichar)
                isymb = 19                              !filled diamond symbol
                if (icount > 27) then
                    legtxt = ""
                endif
                bigDsymbs = .true.
                cthickn = 3.0_wp
                if (high_volat_surr_present .and. (.not. plot_high_volat_cluster) .and. ia(i) == maxClust) then
                    cycle                               !don't plot this point and continue
                else
                    call defzcolval(real(ia(i),kind=wp), dcol(i), colbarlowlim, colbaruplim)
                    if (ia(i) == maxClust) dcol(i) = -21    !for the highly-volatile species cluster use a special WinXP color value
                    call DefinePlotData(xvals, yvals, xyerror0, 1, dcol(i), 3, isymb, 1, legtxt, cthickn)
                    !define a second symbol of the same point with the symbol outline in black:
                    isymb = 5                           !open diamond symbol
                    bigDsymbs = .true.
                    call DefinePlotData(xvals, yvals, xyerror0, 1, 0, 3, isymb, 1, legtxt, cthickn)
                endif	
            enddo
            !---
            !add actual k-means cluster center points:
            if (high_volat_surr_present .and. (.not. plot_high_volat_cluster)) then
                clustN_adj = clustN - 1
            else
                clustN_adj = clustN
            endif
            xvals = 0.0_wp
            yvals = 0.0_wp
            icount = icount + 1
            legtxt = "exact cluster centers"
            if (icount > 28) then
                legtxt = ""
            endif
            if (logxaxis) then
                xvals(1:clustN_adj) = 10.0_wp**clustCenter(1,1:clustN_adj)
            else
                xvals(1:clustN_adj) = clustCenter(1,1:clustN_adj)
            endif
            if (yaxisIsLog) then
                yvals(1:clustN_adj) = 10.0_wp**clustCenter(2,1:clustN_adj)
            else
                yvals(1:clustN_adj) = clustCenter(2,1:clustN_adj)
            endif
            cthickn = 3.0_wp
			isymb = 4                   !X symbol
            call DefinePlotData(xvals, yvals, xyerror0, 1, 0, 3, isymb, 1, legtxt, cthickn)
            !---
            
            ! °°now plot with following parameters °°
            xaxchar = "LgPsat"
            fileoutname = trim(yaxchar)//"_vs_"//trim(xaxchar)//"_"//trim(filenumber(k))//"_res"//trim(lumpResChar)//"_col"//trim(colchar)
            ptitle(1) = trim(fileoutname)
            ptitle(2) = "Activity coefficient ratio vs. pure-component saturation vapor pressure"  !"text describing plot..."
            if (high_volat_surr_present) then
                write(cchar,'(I0.2)') clustN - 1
                ptitle(3) = "Method: "//trim(legendCompProp(k))//", "//"# clusters = "//trim(cchar)//" (+1 high-volat.); $\times$ = cluster center; 'diamond' = surrogate comp."
            else
                write(cchar,'(I0.2)') clustN
                ptitle(3) = "Method: "//trim(legendCompProp(k))//", "//"# clusters = "//trim(cchar)//"; $\times$ = cluster center; 'diamond' = surrogate comp."
            endif
            !set above; lylabel = "$ \gamma_{j,{\rm hex}}/ \gamma_{j,{\rm H_2O}} $"
            xlabel = "$p^{\rm o,sat}_j$ / [Pa]"  !or ? "$C^{\rm o}$ / $[\mu {\rm g}/{\rm m}^3$]"
            rylabel = ""
            zlabel = "$k$-means cluster ID"
            SD = 0.0_wp
            squareplot = .false.
            smallsymbs = .true.
            fixminxax = .false.
            messagebox = .false.
            nogridlines = .true.    !to suppress plotting of background gridlines when already the lumping grid is plotted
            call showzcolbar(colbarlowlim, colbaruplim, zlabel, .true.)   !to set z-axis colour bar properties and show colour bar on plot
            specialxrange = .true.
            specialyrange = .true.
            call PlotNow(ptitle, xlabel, lylabel, rylabel, fileoutname, SD)
        endif
    else
        cycle
    endif
    !========================================================================================================================
enddo !kset
deallocate(c2pv, iz, ia, x2vals, y2vals, key) 


do k = 1,1 !nsets
    !!=========================================================================================================================
    !!= Activity coefficient ratio for each organic with water and hexanediol. vs. OSc
    !!=========================================================================================================================
    if (cPropDataExists(k)) then
        uxlim = -1.0E80_wp
        lxlim = 1.0E80_wp
        uylim = -1.0E80_wp
        lylim = 1.0E80_wp
        cthickn = 3.0_wp !6.0_wp
        yaxisIsLog = .false.
        dcol = 0
        xv = 0.0_wp
        yv = 0.0_wp
        zv = 0.0_wp
        wv = 0.0_wp
        wv(1:krows(k)) = LumpConcData(1:krows(k),2,k)                       !the component mass conc. (after potential lumping)
        colbarlowlim = 0.0_wp
        colbaruplim = max( 2.0_wp, min( 3.0_wp, maxval(CompPropData(1:krows(k),2,k)) ) )  !for max. O:C val shown
        logyaxis = .true.
        yaxisIsLog = logyaxis
        do i = 1,krows(k)
            if (wv(i) > 0.0_wp) then !plot only points of non-zero mass concentration:
                logxaxis = .false.
                logyaxis = yaxisIsLog
                xv(1) = CompPropData(i,4,k)                                 !mean OS_C
                yv(1) = CompPropData(i,6,k)
                zv(1) = CompPropData(i,2,k)
                !!point colors:
                !!O:C ratio is zv(1) = CompPropData(i,2,k)
                !if (krows(k) < 1000) then
                !    bigDsymbs = .true.
                !endif
                call defzcolval(zv(1), dcol(i), colbarlowlim, colbaruplim)  !map the O:C value to the colour bar index number dcol(i)
                write(pvchar,'(F10.2)') zv(1)
                legtxt = "O:C ratio is: "//pvchar
                call DefinePlotData(xv, yv, xyerror0, 1, dcol(i), 3, 21, 5, legtxt, cthickn)
                uxlim = max(uxlim, xv(1))
                lxlim = min(lxlim, xv(1))
                uylim = max(uylim, yv(1))
                lylim = min(lylim, yv(1))
            endif
        enddo
        
        ! °°now plot with following parameters °°
        fileoutname = "activitycoeffratio_avgCarbon_OS"//"_res"//trim(lumpResChar)//"_"//trim(filenumber(k))
        ptitle(1) = trim(fileoutname)
        ptitle(2) = "Activity coefficient ratio vs. average Oxidation State of Carbon"     !"text describing plot..."
        ptitle(3) = "Input file and method: "//trim(legendCompProp(k))//",  "//"grid resolution: "//trim(lumpResChar)
        lylabel = "$ \gamma_{j,{\rm hex}}/ \gamma_{j,{\rm H_2O}} $"
        xlabel = "$\overline{\rm OS}_{\rm C}$"
        rylabel = ""
        zlabel = "O:C ratio"
        SD = 0.0_wp
        squareplot = .false.
        smallsymbs = .true.
        fixminxax = .false.
        messagebox = .false.
        nogridlines = .true.    !suppress coordinate gridlines
        call SetCustomColPalette('Viridis')     !load custom-made colour palette in Mod_createplot
        call showzcolbar(colbarlowlim, colbaruplim, zlabel, .false.)  !to set z-axis colour bar properties and show colour bar on plot
        specialxrange = .true.  !to plot with a defined x-axis range
        specialyrange = .true.
        if (k == 1) then    !set the axis limits based on full system and use the same for lumped systems (i.e. k > 1)
            if (logxaxis) then
                lowerxlimit = 10.0_wp**floor(log10(lxlim - 0.05_wp*lxlim))
                if (lxlim - 0.05_wp*lxlim > 4.0_wp*lowerxlimit) then
                    lowerxlimit = lowerxlimit*3.0_wp
                endif
                upperxlimit = 10.0_wp**ceiling(log10(uxlim + 0.05_wp*uxlim))
                if (uxlim + 0.05_wp*uxlim < 0.3_wp*upperxlimit) then
                    upperxlimit = 0.4_wp*upperxlimit  
                endif
            else
                axextend = max(abs(0.05_wp*lxlim), abs(0.05_wp*uxlim))
                lowerxlimit = lxlim -axextend
                upperxlimit = uxlim +axextend
            endif
            if (yaxisIsLog) then
                lowerylimit = 10.0_wp**floor(log10(lylim - 0.05_wp*lylim))
                if (lylim > 4.0_wp*lowerylimit) then
                    lowerylimit = lowerylimit*3.0_wp
                endif
                upperylimit = 10.0_wp**ceiling(log10(uylim + 0.05_wp*uylim))
                if (uylim < 0.3_wp*upperylimit) then
                    upperylimit = 0.4_wp*upperylimit  
                endif
            else
                axextend = max(abs(0.05_wp*lylim), abs(0.05_wp*uylim))
                lowerylimit = lylim -axextend
                upperylimit = uylim +axextend
            endif
        endif
        call PlotNow(ptitle, xlabel, lylabel, rylabel, fileoutname, SD)
    endif
    !========================================================================================================================
enddo !k

do k = 1,1 !nsets    
    !!=========================================================================================================================
    !!= Activity coefficient ratio for each organic with water and hexanediol. vs. O:C
    !!=========================================================================================================================
    if (cPropDataExists(k)) then
        uxlim = -1.0E80_wp
        lxlim = 1.0E80_wp
        uylim = -1.0E80_wp
        lylim = 1.0E80_wp
        specialxrange = .false.
        specialyrange = .false.
        bargylimits = .false.
        logxaxis = .false.
        yaxisIsLog = .false.
        dcol = 0
        xv = 0.0_wp
        yv = 0.0_wp
        wv = 0.0_wp
        wv(1:krows(k)) = LumpConcData(1:krows(k),2,k)  !the component mass conc. (after potential lumping)
        colbarlowlim = 0.0_wp
        colbaruplim = max( 2.0_wp, min( 3.0_wp, maxval(CompPropData(1:krows(k),2,k)) ) )  !for max. O:C val shown
        logyaxis = .true.
        yaxisIsLog = logyaxis
        do i = 1,krows(k)
            if (wv(i) > 0.0_wp) then !plot only points of non-zero mass concentration:
                logxaxis = .false.
                logyaxis = yaxisIsLog
                xv(1) = CompPropData(i,2,k)
                yv(1) = CompPropData(i,6,k)
                !!point colors:
                !!O:C ratio is CompPropData(i,2,k)
                !if (krows(k) < 1000) then
                !    bigDsymbs = .true.
                !endif
                call defzcolval(CompPropData(i,2,k), dcol(i), colbarlowlim, colbaruplim)  !map the O:C value to the colour bar index number dcol(i)
                write(pvchar,'(F10.2)') xv(1)
                legtxt = "O:C ratio is: "//pvchar
                call DefinePlotData(xv, yv, xyerror0, 1, dcol(i), 3, 21, 5, legtxt, cthickn)
                uxlim = max(uxlim, xv(1))
                lxlim = min(lxlim, xv(1))
                uylim = max(uylim, yv(1))
                lylim = min(lylim, yv(1))
            endif
        enddo

        ! °°now plot with following parameters °°
        fileoutname = "activitycoeffratio_OtoC"//"_res"//trim(lumpResChar)//"_"//trim(filenumber(k))
        ptitle(1) = trim(fileoutname)
        ptitle(2) = "Activity coefficient ratio vs. oxygen to carbon ratio of components"!"text describing plot..."
        ptitle(3) = "Input file and method: "//trim(legendCompProp(k))//",  "//"grid resolution: "//trim(lumpResChar)
        lylabel = "$ \gamma_{j,{\rm hex}}/ \gamma_{j,{\rm H_2O}} $"
        xlabel = "O:C ratio"
        rylabel = ""
        zlabel = "O:C ratio"
        SD = 0.0_wp
        squareplot = .false.
        smallsymbs = .true.
        fixminxax = .false.
        messagebox = .false.
        nogridlines = .true.
        call SetCustomColPalette('Viridis')     !load custom-made colour palette in Mod_createplot
        call showzcolbar(colbarlowlim, colbaruplim, zlabel, .false.)   !to set z-axis colour bar properties and show colour bar on plot
        specialxrange = .true.  !to plot with a defined x-axis range
        specialyrange = .true.  !to plot with a defined y-axis range
        if (k == 1) then
            if (logxaxis) then
                lowerxlimit = 10.0_wp**floor(log10(lxlim -0.05_wp*lxlim))
                if (lxlim - 0.05_wp*lxlim > 4.0_wp*lowerxlimit) then
                    lowerxlimit = lowerxlimit*3.0_wp
                endif
                upperxlimit = 10.0_wp**ceiling(log10(uxlim +0.05_wp*uxlim))
                if (uxlim + 0.05_wp*uxlim < 0.3_wp*upperxlimit) then
                    upperxlimit = 0.4_wp*upperxlimit  
                endif
            else
                axextend = max(abs(0.05_wp*lxlim), abs(0.05_wp*uxlim))
                lowerxlimit = lxlim -axextend
                upperxlimit = uxlim +axextend
            endif
            if (yaxisIsLog) then
                lowerylimit = 10.0_wp**floor(log10(lylim -0.05_wp*lylim))
                if (lylim > 4.0_wp*lowerylimit) then
                    lowerylimit = lowerylimit*3.0_wp
                endif
                upperylimit = 10.0_wp**ceiling(log10(uylim +0.05_wp*uylim))
                if (uylim < 0.3_wp*upperylimit) then
                    upperylimit = 0.4_wp*upperylimit  
                endif
            else
                axextend = max(abs(0.05_wp*lylim), abs(0.05_wp*uylim))
                lowerylimit = lylim -axextend
                upperylimit = uylim +axextend
            endif
        endif
        call PlotNow(ptitle,xlabel,lylabel,rylabel,fileoutname,SD)
    endif
    !========================================================================================================================
enddo !k

!deallocate the used dynamic arrays:
call DeallocCreateplotArrays() 
deallocate(CompPropData, LumpConcData, cpv, xv, yv, zv, wv, xvals, yvals, xyerror0, STAT=allocstat)
if (allocated(clustSurrogateIndx)) then
    deallocate( clustCenter, clustSurrogateIndx, clustPop, compCluster )
endif

write(*,*) ""
write(*,'(A,/)') "Program successfully completed"
write(*,'(A,/)') "Press ENTER to leave the program."
read(*,*)

end program CustomizedPlotting_2D_framework