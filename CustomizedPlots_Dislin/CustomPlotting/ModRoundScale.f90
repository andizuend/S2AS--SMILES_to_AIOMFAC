!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module including subroutines for rounding of numbers up/down and for determining   *
!*   scales and ranges of values (e.g. used for plotting axes)                          *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2009                                                            *
!*   -> latest changes: 2020-02-05                                                      *
!*                                                                                      *
!****************************************************************************************
module ModRoundScale

use Mod_kinds, only : wp

implicit none
public

contains
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Function for rounding floating point numbers up/down to a target number of         * 
    !*   significant figures; e.g. useful for axis scaling in plots.                        *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, 2007 - 2009, 2013                                              *
    !*   Dept. Chem. Engineering, California Institute of Technology, 2011 - 2012           *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University, 2013 - present          *
    !*                                                                                      *
    !*   -> created:        2009                                                            *
    !*   -> latest changes: 2020-02-05                                                      *
    !*                                                                                      *
    !****************************************************************************************
    pure elemental function roundnicely(va, nsigfig, updown)  result(val)

    implicit none
    !vars from interface:
    integer,intent(in) :: nsigfig          !number of significant digits to be used for rounding level of numeric value "va"
    real(wp),intent(in) :: va                  !value to be rounded up/down
    character(len=*),intent(in)  :: updown    !'up' or 'down' to indicate rounding direction
    real(wp) :: val                            !the rounded return value (result) 
    !local vars:
    real(wp) :: rscale, magnitude
    logical :: negval, upwards
    !...........................................................................

    if (updown(1:1) == "u" .or. updown(1:1) == "U") then !requesting rounding up to larger value
        upwards = .true.
    else
        upwards = .false.
    endif
    
    if (va < 0.0_wp) then
        negval = .true.
        val = -va
        if (upwards) then !flip the direction for correct rounding of negative values in up/down direction
            upwards = .false.
        else
            upwards = .true.
        endif
    else
        negval = .false.
        val = va
    endif

    !determine lower decadal magnitude level of value:
    if (val >= 1.0_wp) then !leading figure is not a zero
        magnitude = aint( log10(val) )
    else if (val > 0.0_wp) then !leading figure is a zero
        magnitude = aint( log10(val) ) -1.0_wp
    else
        magnitude = aint( log10(tiny(val)) ) +1.0_wp
    endif

    !scale by determined magnitude to isolate significant digits:
    val = val/(10.0_wp**magnitude)

    !determine rounding level:
    if (nsigfig > 3) then
        magnitude = magnitude + 0 !for debugging breakpoint
    endif
    rscale = 10.0_wp**(nsigfig-1)
    if (upwards) then !rounding up to larger value level
        val = anint(val*rscale+0.4999999999_wp) !round scaled value to nearest whole number
    else  !round down
        val = aint(val*rscale)
    endif
    val = (val/rscale)*10.0_wp**magnitude !re-scale
    if (negval) then
        val = -val
    endif

    end function roundnicely
    !=================================================================================================================


    !****************************************************************************
    !     ALGORITHM AS 96  APPL. STATIST. (1976) VOL.25, NO.1
    !     edited and changed by Andi Zuend, IAC, ETH Zurich, 2009
    !
    !     Given extreme values fmn, fmx, and the need for a scale with around n
    !     marks, calculates value for the lowest scale mark (valmin) and
    !     step length (step) and highest scale mark (valmax).
    !****************************************************************************
    pure subroutine scalenicely(fmn, fmx, n, valmin, valmax, step, ifault)

    implicit none
    !interface:
    real(wp),intent(in) :: fmn, fmx
    real(wp),intent(out) :: valmin, step, valmax
    integer,intent(in) :: n
    integer,intent(out) :: ifault
    !local variables:
    real(wp) :: rn, range, x, fmin, fmax
    real(wp),parameter :: deps = 5.0_wp*epsilon(deps), tol = 1.0E-9_wp !the tolerance level

    fmax = fmx
    fmin = fmn

    !Test for valid parameter values:
    ifault = 1
    if (fmax < fmin .or. n  <=  1) return
    ifault = 0
    rn = real(n,kind=wp)
    x = abs(fmax)
    if (x == 0.0_wp) x = 1.0_wp
    if ((fmax - fmin)/x > tol) then
        ! continue
    else  !all values effectively equal
        if (fmax < 0.0_wp) then
            fmax = 0.0_wp
        else if (fmax == 0.0_wp) then
            fmax = 1.0_wp
        else
            fmin = 0.0_wp
        endif
    endif

    !Set valmin:
    if ((fmax - fmin) < 0.1_wp) then
        valmin = roundnicely(fmin, 3, 'down') ! rounding down
    else
        valmin = roundnicely(fmin, 2, 'down')
    endif

    !check whether valmin could be set to zero:
    if (abs(valmin) < (fmax-fmin)*0.05_wp) then !0.05
        valmin = 0.0_wp
    endif

    !Estimate a preliminary valmax:
    if ((fmax - valmin) < 0.01_wp) then
        valmax = roundnicely(fmax, 3, 'up') ! rounding up
    else
        valmax = roundnicely(fmax, 2, 'up') ! rounding up
    endif

    step = (valmax - valmin)/(rn+deps)
    step = roundnicely((step+deps), 1, 'up') ! rounding up
    range = step*rn

    valmax = valmin + range
    !check whether the full range is covered by the found values:
    if (step > 0.0_wp) then
        do  
            if (valmax < fmax) then
                rn = rn+1.0_wp
                range = step*rn
                valmax = valmin + range
            else
                exit
            endif
        enddo
    else
        valmax = 0.0_wp
    endif

    !check whether valmax could be set to zero:
    if (abs(valmax) < (fmax-fmin)*0.05_wp) then
        valmax = 0.0_wp
    endif

    end subroutine scalenicely
    !--------------------------------------------------------------

end module ModRoundScale