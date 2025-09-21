module Mod_PlotData
    
use Mod_kinds, only : wp

public :: ReadDataFiles

contains
    
    !****************************************************************************************
    !*                                                                                      *
    !*   Subroutine to read data from input text-files for use with custom plots.           *
    !*   Here customized for reading lumping schemes output data as input.                  *
    !*                                                                                      *
    !*   (c) Andi Zuend,                                                                    *
    !*   Dept. Atmospheric & Oceanic Sciences, McGill University                            *
    !*                                                                                      *
    !*   -> created:        2019-10-26                                                      *
    !*   -> latest changes: 2025-07-05                                                      *
    !*                                                                                      *
    !****************************************************************************************

    subroutine ReadDataFiles(fileStartNo, lumpResChar, nsets, inrows, CompPropData, LumpConcData, xgridl, ygridl, &
                            & legendCompProp, legendLumpedConc, cPropDataExists, filenumber, YaxisChoice, &
                            & clustNum, clustCenter, clustSurrogateIndx, clustPop, compCluster)

    use ModCreateplot, ONLY: rdefault

    implicit none

    !interface variables:
    integer,intent(in) :: fileStartNo
    character(len=5),intent(in) :: lumpResChar
    integer,intent(out) :: nsets, inrows
    real(wp),dimension(:,:,:),intent(out) :: CompPropData, LumpConcData
    real(wp),dimension(:),allocatable,intent(out) :: xgridl, ygridl
    character(len=75),dimension(:),intent(out) :: legendCompProp, legendLumpedConc
    logical,dimension(:),intent(out) :: cPropDataExists
    character(len=10),dimension(:),intent(out) :: filenumber
    character(len=30),dimension(:),intent(out) :: YaxisChoice
    integer,intent(out) :: clustNum
    real(wp),dimension(:,:),allocatable,intent(out) :: clustCenter 
    integer,dimension(:),allocatable,intent(out) :: clustSurrogateIndx, clustPop !dimension covers all clusters
    integer,dimension(:),allocatable,intent(out) :: compCluster                  !dimension covers all components
    !...
    !local variables:
    integer :: i, ii, istat, k, kk, krows, maxinf, maxrows, ndim, ndl, norg1, un
    character(len=120) :: inputfolderpath
    character(len=50) :: dummy, namestem1, namestem2
    character(len=50),dimension(5) :: methodchar
    character(len=300) :: namex, fname
    character(len=300),dimension(2*size(legendCompProp)) :: infname
    logical :: fexists
    !..................

    !===== Initialize and set parameters ===========
    inputfolderpath = "../../Output_lumping/"  !relative path from CustomPlots_Dislin project file folder location to input folder "Output_lumping"
    namestem1 = "SystemCompProp_"
    namestem2 = "LumpedConc_"
    methodchar(1) = "_FullSystem"
    methodchar(2) = "_Medoid_"
    methodchar(3) = "_Midpoint_"
    methodchar(4) = "_Weighted_Medoid_"
    methodchar(5) = "_Kmeans_"
    !---
    maxinf = 2*size(legendCompProp)
    maxrows = size(CompPropData, dim=1)
    ndim = size(CompPropData, dim=2)
    ndl = size(LumpConcData, dim=2)
    cPropDataExists = .false.
    infname = ""
    legendCompProp = ""
    legendLumpedConc = ""
    CompPropData = rdefault     !initialize to an unusual value
    LumpConcData = rdefault
    filenumber = ""
    YaxisChoice = ""
    !---
    !determine existing file names to be read:
    ii = 0
    do k = 1,maxinf
        !construct a possible filename and check for its existence
        if (k <= maxinf/2) then
            kk = k
            namex = trim(namestem1)
        else
            kk = k - maxinf/2
            namex = trim(namestem2)
        endif
        write(dummy, '(I4.4)') fileStartNo +kk-1
        do i = 1,2      !check both cases, with resolution indicated or not
            select case(i) 
            case(1)
                fname = trim(namex)//trim(dummy)//trim(methodchar(kk))//trim(lumpResChar)//'.txt'
            case(2)
                fname = trim(namex)//trim(dummy)//trim(methodchar(kk))//'.txt'
            end select
            fname = trim(inputfolderpath)//trim(fname)
            inquire(FILE = trim(fname), EXIST = fexists)
            if (fexists) then
                if (kk <= k) then
                    cPropDataExists(kk) = .true.
                    filenumber(kk) = trim(dummy)
                endif
                infname(k) = trim(fname)       !save file name including relative path
                if (i == 1) then
                    exit  !leave the inner i-loop
                endif
            endif
        enddo !i
    enddo !k

    !===========================================
    !check if input file exists and read its content if true

    !read data from files into an array containing the data for the "curves":
    nsets = 0
    krows = 0
    maxinf = size(infname(:))
    do k = 1,maxinf/2
        kk = maxinf/2 + k   !second half index for infname
        !read CompPropData if it exists:
        fname = trim(infname(k))
        if (len_trim(fname) > 1 .and. index(fname, trim(namestem1)) > 0) then   !a SystemCompProp file...
            i = 1 + index(fname, "/", BACK = .true.)
            ii = index(fname, ".txt")
            nsets = k
            legendCompProp(k) = trim(fname(i:ii-1))
            open (NEWUNIT = un, FILE = trim(fname), STATUS="OLD")
            !---
            read(un,*) dummy    !read first line
            read(un,*) dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy
            read(un,*) dummy, dummy, YaxisChoice(k)  !read the YaxisChoice used for this file
            read(un,*) dummy, dummy, dummy
            read(un,*)          !read empty line
            read(un,*) dummy, dummy
            read(un,*)
            read(un,*) dummy    !read column headers
            do i = 1,maxrows    !read the data rows
                read(un,*,IOSTAT=istat) CompPropData(i,1:ndim,k) !read the x-y-z-...data of each line of the file
                if (istat /= 0) then
                    exit
                endif
            enddo
            inrows = i-1
            krows = max(krows, inrows)
            wait(un)
            close(un)
        endif
        
        !read LumpedConc data if it exists:
        fname = trim(infname(kk))
        if (len_trim(fname) > 1 .and. index(fname, trim(namestem2)) > 0) then !a LumpedConc file...
            i = 1 + index(fname, "/", BACK = .true.)
            ii = index(fname, ".txt")
            nsets = k
            legendLumpedConc(k) = trim(fname(i:ii-1))
            open (NEWUNIT = un, FILE = trim(fname), STATUS="OLD")
            !---
            read(un,*) dummy    !read first line
            read(un,*) dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy
            read(un,*) dummy, dummy, YaxisChoice(k)  !read the YaxisChoice used for this file
            read(un,*) dummy, dummy, dummy
            read(un,*)          !read empty line
            read(un,*) dummy, dummy, dummy
            read(un,*)
            read(un,*) dummy    !read column headers
            do i = 1,maxrows    !read the data rows
                read(un,*,IOSTAT=istat) LumpConcData(i,1:ndl,k) !read the x-y-z-...data of each line of the file
                if (istat /= 0) then
                    exit
                endif
            enddo
            inrows = i-1
            krows = max(krows, inrows)
            wait(un)
            close(un)
        endif
    enddo !k
    
    inrows = krows !for output
    if (inrows+1 > maxrows) then !check for potential issue with too small maxrows value
        write(*,*) "# WARNING # from ReadDataFiles: the value of inrows read it as the limit of set maximum rows." 
        write(*,*) "This likely indicates that parameter 'maxrows' was set at a too low value in CustomPlot."
        write(*,*) "Make sure to increase the value of 'maxrows'."
        read(*,*) !wait for user action
    endif

    !Read grid line coordinates file:
    write(dummy,'(I4.4)') fileStartNo
    fname = 'GridlineCoord_'//trim(dummy)//'_'//trim(lumpResChar)//'.txt'
    fname = trim(inputfolderpath)//trim(fname)
    inquire(FILE = trim(fname), EXIST = fexists)
    if (fexists) then
        open (NEWUNIT = un, FILE = trim(fname), STATUS="OLD")
        read(un,*) dummy    !read first line
        read(un,*)      
        read(un,*) dummy
        read(un,*)  
        read(un,*) dummy, k
        read(un,*) dummy, kk
        allocate( xgridl(k), ygridl(kk) )
        read(un,*)         
        read(un,*) dummy
        read(un,*,IOSTAT=istat) xgridl(1:k)     !read the x-axis grid line coordinates
        read(un,*) dummy
        read(un,*,IOSTAT=istat) ygridl(1:kk) 
        wait(un)
        close(un)
    endif

    !Read additional k-means cluster population data file (when present):
    do k = 2,size(methodchar)
        write(dummy,'(I4.4)') fileStartNo-1 + k
        fname = 'KmeansClusters_'//trim(dummy)//'_'//trim(lumpResChar)//'.txt'
        fname = trim(inputfolderpath)//trim(fname)
        inquire(FILE = trim(fname), EXIST = fexists)
        if (fexists) then
            open (NEWUNIT = un, FILE = trim(fname), STATUS="OLD")
            read(un,*) dummy    !read first line
            read(un,*)      
            read(un,*) dummy
            read(un,*) dummy, dummy, dummy, dummy, clustNum
            read(un,*)      
            read(un,*) dummy
            read(un,*)      
            read(un,*) dummy
            read(un,*) dummy
            !--
            allocate( clustCenter(2,clustNum), clustSurrogateIndx(clustNum), clustPop(clustNum) )
            !read cluster centers and population information:
            do i = 1,clustNum
                read(un,*,IOSTAT=istat) kk, clustCenter(1:2,i), clustSurrogateIndx(i), clustPop(i)
            enddo
            !--
            read(un,*)      
            read(un,*)  
            read(un,*) dummy
            read(un,*) dummy
            read(un,*) dummy
            !read list of components and to which cluster they belong:
            norg1 = LumpConcData(1,1,1)
            allocate( compCluster(norg1:inrows+norg1-1) )
            !allocate(compCluster(1:inrows))
            !kk = 0
            do !until exit
                !kk = kk+1
                read(un,*,IOSTAT=istat) i, compCluster(i) 
                if (istat /= 0) then
                    exit
                endif
            enddo
            !--
            wait(un)
            close(un)
            exit !done with do-loop
        endif
    enddo

    end subroutine ReadDataFiles

end module Mod_PlotData