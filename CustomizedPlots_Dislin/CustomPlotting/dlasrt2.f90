
subroutine DLASRT2( ID, N, D, KEY, INFO )

use Mod_kinds, only : wp
!
!  -- ScaLAPACK routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
character          ID
integer ::            INFO, N
!     ..
!     .. Array Arguments ..
integer ::            KEY(*)
real(wp) ::   D(*)
!     ..
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments
!  =========
!
!  ID      (input) character!1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) integer
!          The length of the array D.
!
!  D       (input/output) real(wp) array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  KEY     (input/output) integer array, dimension (N)
!          On entry, KEY contains a key to each of the entries in D()
!          Typically, KEY(I) = I for all I
!          On exit, KEY is permuted in exactly the same manner as
!          D() was permuted from input to output
!          Therefore, if KEY(I) = I for all I upon input, then
!          D_out(I) = D_in(KEY(I))
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
integer,parameter :: select = 20
!     ..
!     .. Local Scalars ..
integer ::   DIR, ENDD, I, J, START, STKPNT, TMPKEY
real(wp) ::   D1, D2, D3, DMNMX, TMP
!     ..
!     .. Local Arrays ..
integer :: STACK( 2, 32 )
!     ..
!     .. External Functions ..
!logical            LSAME
!EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
!
INFO = 0
DIR = -1
if (trim(ID) == 'D') then
    DIR = 0
else if (trim(ID) == 'I') then
    DIR = 1
end if
if( DIR.EQ.-1 ) then
    INFO = -1
else if( N.LT.0 ) then
    INFO = -2
end if
if( INFO.NE.0 ) then
    !call XERBLA( 'DLASRT2', -INFO )
    return
end if
!
!     Quick return if possible
!
if( N.LE.1 ) return
!
STKPNT = 1
STACK( 1, 1 ) = 1
STACK( 2, 1 ) = N
10  continue
START = STACK( 1, STKPNT )
ENDD = STACK( 2, STKPNT )
STKPNT = STKPNT - 1
if( ENDD - START.GT.0 ) then
    !
    !        Do Insertion sort on D( START:ENDD )
    !
    if( DIR.EQ.0 ) then
        !
        !           Sort into decreasing order
        !
        do 30 I = START + 1, ENDD
            do 20 J = I, START + 1, -1
                if( D( J ).GT.D( J-1 ) ) then
                    DMNMX = D( J )
                    D( J ) = D( J-1 )
                    D( J-1 ) = DMNMX
                    TMPKEY = KEY( J )
                    KEY( J ) = KEY( J-1 )
                    KEY( J-1 ) = TMPKEY
                else
                    go to 30
                end if
20              continue
30          continue
        !
    else
        !
        !           Sort into increasing order
        !
        do 50 I = START + 1, ENDD
            do 40 J = I, START + 1, -1
                if( D( J ).LT.D( J-1 ) ) then
                    DMNMX = D( J )
                    D( J ) = D( J-1 )
                    D( J-1 ) = DMNMX
                    TMPKEY = KEY( J )
                    KEY( J ) = KEY( J-1 )
                    KEY( J-1 ) = TMPKEY
                else
                    go to 50
                end if
40              continue
50          continue
        !
    end if
    !
else if( ENDD - START.GT.select ) then
    !
    !        Partition D( START:ENDD ) and stack parts, largest one first
    !
    !        Choose partition entry as median of 3
    !
    D1 = D( START )
    D2 = D( ENDD )
    I = ( START+ENDD ) / 2
    D3 = D( I )
    if( D1.LT.D2 ) then
        if( D3.LT.D1 ) then
            DMNMX = D1
        else if( D3.LT.D2 ) then
            DMNMX = D3
        else
            DMNMX = D2
        end if
    else
        if( D3.LT.D2 ) then
            DMNMX = D2
        else if( D3.LT.D1 ) then
            DMNMX = D3
        else
            DMNMX = D1
        end if
    end if
    !
    if( DIR.EQ.0 ) then
        !
        !           Sort into decreasing order
        !
        I = START - 1
        J = ENDD + 1
60          continue
70          continue
        J = J - 1
        if( D( J ).LT.DMNMX ) go to 70
80          continue
        I = I + 1
        if( D( I ).GT.DMNMX ) go to 80
        if( I.LT.J ) then
            TMP = D( I )
            D( I ) = D( J )
            D( J ) = TMP
            TMPKEY = KEY( J )
            KEY( J ) = KEY( I )
            KEY( I ) = TMPKEY
            go to 60
        end if
        if( J-START.GT.ENDD - J-1 ) then
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = START
            STACK( 2, STKPNT ) = J
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = J + 1
            STACK( 2, STKPNT ) = ENDD
        else
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = J + 1
            STACK( 2, STKPNT ) = ENDD
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = START
            STACK( 2, STKPNT ) = J
        end if
    else
        !
        !           Sort into increasing order
        !
        I = START - 1
        J = ENDD + 1
90          continue
100         continue
        J = J - 1
        if( D( J ).GT.DMNMX ) go to 100
110         continue
        I = I + 1
        if( D( I ).LT.DMNMX ) go to 110
        if( I.LT.J ) then
            TMP = D( I )
            D( I ) = D( J )
            D( J ) = TMP
            TMPKEY = KEY( J )
            KEY( J ) = KEY( I )
            KEY( I ) = TMPKEY
            go to 90
        end if
        if( J-START.GT.ENDD -J-1 ) then
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = START
            STACK( 2, STKPNT ) = J
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = J + 1
            STACK( 2, STKPNT ) = ENDD
        else
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = J + 1
            STACK( 2, STKPNT ) = ENDD
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = START
            STACK( 2, STKPNT ) = J
        end if
    end if
end if
if( STKPNT.GT.0 ) go to 10
!
return
!
!     End of DLASRT2
!
end