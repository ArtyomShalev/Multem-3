!Multem, version 3 is a renewed program for transmission and 
!band-structure calculations of photonic crystals. This source file
!and libmultem2b.f90 are basically modified code of  MULTEM 2: 
!A NEW VERSION OF THE PROGRAM FOR TRANSMISSION AND BAND-STRUCTURE 
!CALCULATIONS OF PHOTONIC CRYSTALS.  
!N. STEFANOU, V. YANNOPAPAS, A. MODINOS. 
!IN COMP. PHYS. COMMUN. 132 (2000) 189
!HERE STARTS THE FORTRAN SOURCE CODE
!=======================================================================

program multem
    use libmultem2b
    implicit none
    !     ------------------------------------------------------------------
    !     A B S T R A C T
    !     this program calculates either the absorbance, reflectivity and
    !     transmittance of light by a finite slab consisting of
    !     homogeneous plates and multilayers of spherical particles
    !     arranged in a two-dimensional bravais lattice, or the  complex
    !     photonic band structure of such an infinite periodic structure.
    !
    !     D E S C R I P T I O N    O F    I N P U T    D A T A
    !     ktype=     1: the direction of an incident em wave is specified
    !                   by the polar angles of incidence "theta" and "fi".
    !                   the program calculates the transmission,reflection
    !                   and  absorption coefficients of a finite  slab
    !                2: the direction of an incident em wave is specified
    !                   by the components of the wavevector parallel to
    !                   the  interfaces of the structure:
    !                   aq(1) and aq(2) (and the frequency). The
    !                   program calculates  the transmission, reflection,
    !                   absorption coefficients of a finite slab
    !                3: the program calculates the photonic  complex band
    !                   structure of such an infinite periodic  structure
    !                   for a wavevector with components parallel to  the
    !                   interfaces of the structure: aq(1) and aq(2)
    !     kscan=     1: scanning over frequencies
    !                2: scanning over wavelengths
    !     kemb        : indicates the presence (=1) or absence (=0) of a
    !                   different embedding medium
    !     lmax        : cutoff in spherical waves expansions
    !     ncomp       : number of different components in the unit slice.
    !                   their type is specified by the integer array
    !                   it(icomp)
    !     it=        1: homogeneous plate of thickness "d"
    !                2: multilayer  of spherical  particles arranged in  a
    !                   2d  bravais lattice. Each layer consists of "nplan"
    !                   non-primitive  planes of spheres with the same 2-d
    !                   periodicity. The number of unit layers is equal to
    !                   2**(nlayer-1).
    !     dl, dr      : position vectors indicating the origin on the left
    !                   and on the right of the  unit, respectively. Both
    !                   are directed from left to right.
    !     al          : primitive  translation  vector  of the  unit slice
    !                   (effective only for band structure calculation).it
    !                   is given in program units.
    !     nunit       : specifies the number of unit slices (2**(nunit-1))
    !                   of the sample
    !     alpha,alphap: length of primitive vectors of the two-dimensional
    !                   lattice. In program units the size of alpha serves
    !                   as the unit length. Thus alpha must be equal to
    !                   1.d0
    !     fab         : angle (in deg) between alpha and alphap
    !     rmax        : upper limit for the length of  reciprocal  lattice
    !                   vectors (in units of 1/alpha) which  must be taken
    !                   into account
    !     zinf,zsup   : minimum  and  maximum  values of frequency (in
    !                   program units: omega*alpha/c), or wavelength (in
    !                   program units: lamda/alpha  ), according to the
    !                   value of kscan. c and lamda refer to vacuum
    !     np          : number of equally spaced points between zinf, zsup
    !     polar       : polarization ('S ' or  'P ') of the incident light
    !     aq(1,2)     : wavevector components parallel  to the  interfaces
    !                   of the structure (xy-plane) in units of 2*pi/alpha
    !     theta,fi    : polar angles of incidence (in deg) of the incident
    !                   light
    !     fein        : angle (in deg) specifying the direction  of  the
    !                   polarization vector for normal  incidence.  not
    !                   effective otherwise
    !     eps*,mu*    : relative dielectric functions and magnetic permea-
    !                   bilities of the various media
    !     ------------------------------------------------------------------
    integer, parameter :: ncompd = 8, npland = 4
    integer       lmax, i, ktype, kscan, ncomp
    integer       np, nunit, icomp, kemb,  ipl
    real(dp)      alpha, rmax
    real(dp)      zinf, zsup, fab, alphap, theta, fi, fein
    complex(dp)   muembl, epsembl, muembr, epsembr, d2, d1
    character(2)  polar
    character(17) text1(2)
    character(5)  dummy

    integer    it(ncompd), nlayer(ncompd), nplan(ncompd)
    real(dp)   dl(3, ncompd, npland), dr(3, ncompd, npland)
    real(dp)   s(ncompd, npland), al(3), d(ncompd), aq(2)
    complex(dp) eps2(ncompd), eps3(ncompd)
    complex(dp) mu1(ncompd), mu2(ncompd), mu3(ncompd), eps1(ncompd)
    complex(dp) musph(ncompd, npland), epssph(ncompd, npland)

    data text1/'homogeneous plate', 'photonic crystal'/

    integer, allocatable ::  multipole_type(:), multipole_order(:), m_projection(:), multipole_combination(:, :)
    integer :: is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected
  
    read(10, 200) ktype, kscan, kemb, lmax, ncomp, nunit
    if(ktype<=0.or.ktype>=4) stop 'illegal input value of ktype'
    if(kscan<=0.or.kscan>=3) stop 'illegal input value of kscan'
    if(kemb<0.or.kemb>=2)   stop 'illegal input value of kemb '
    if(ncomp<=0.or.ncomp>ncompd)&
            stop 'illegal input value of ncomp'
    if(nunit<=0)           stop 'illegal input value of nunit'
    read(10, 202) alpha, alphap, fab, rmax
    fab = fab * pi / 180.0_dp
    read(10, 203) np, zinf, zsup
    if(np<=1)                  stop 'illegal input value of  np '
    if(ktype>=2) then
        read(10, 204) aq(1), aq(2), polar, fein
        fein = fein * pi / 180.d0
        aq(1) = 2.d0 * pi * aq(1)
        aq(2) = 2.d0 * pi * aq(2)
        if(ktype<3) then
            write(6, 222)
        else
            write(6, 223)
        endif
        if(ktype==2) write(6, 207) aq(1), aq(2), polar
        if(ktype==3) write(6, 225) aq(1), aq(2)
    else
        read(10, 204) theta, fi, polar, fein
        write(6, 208) theta, fi, polar
        fein = fein * pi / 180.d0
        theta = theta * pi / 180.d0
        fi = fi * pi / 180.d0
    endif
    do icomp = 1, ncomp
        read(10, 201) it(icomp)
        if(it(icomp)<=0.or.it(icomp)>2)&
                stop 'illegal component type'
        write(6, 209) icomp, text1(it(icomp))
        if(it(icomp)==1) then
            read(10, 204) d(icomp)
            read(10, 205) mu1(icomp), eps1(icomp), mu2(icomp), eps2(icomp), &
                    mu3(icomp), eps3(icomp)
            write(6, 210) mu1(icomp), mu2(icomp), mu3(icomp), eps1(icomp), &
                    eps2(icomp), eps3(icomp)
            read(10, *) dummy, (dl(i, icomp, 1), i = 1, 3)
            read(10, *) dummy, (dr(i, icomp, 1), i = 1, 3)
        else
            read(10, 205) mu1(icomp), eps1(icomp)
            if(dble(mu1(icomp))<=0.d0.or.dble(eps1(icomp))<=0.d0)&
                    then
                write(6, 226)
                stop
            endif
            read(10, 201) nplan(icomp), nlayer(icomp)
            do ipl = 1, nplan(icomp)
                read(10, 206) s(icomp, ipl), musph(icomp, ipl), epssph(icomp, ipl)
                read(10, *) dummy, (dl(i, icomp, ipl), i = 1, 3)
                read(10, *) dummy, (dr(i, icomp, ipl), i = 1, 3)
            end do
            write(6, 211)  mu1(icomp), (musph(icomp, ipl), ipl = 1, nplan(icomp))
            write(6, 220) eps1(icomp), (epssph(icomp, ipl), ipl = 1, nplan(icomp))
            write(6, 224) (s(icomp, ipl), ipl = 1, nplan(icomp))
            write(6, 212) 2**(nlayer(icomp) - 1)
        endif
    end do

    d1 = sqrt(mu1(1) * eps1(1))
    d2 = sqrt(mu1(ncomp) * eps1(ncomp))
    if(it(ncomp)==1) d2 = sqrt(mu3(ncomp) * eps3(ncomp))
    if(aimag(d1)/=0.d0) then
        write(6, 227)
        stop
    endif
    if(aimag(d2)/=0.d0) then
        write(6, 228)
        stop
    endif
    if(ktype/=3) then
        write(6, 221) 2**(nunit - 1)
        if(kemb==1) then
            read(10, 205) muembl, epsembl
            read(10, 205) muembr, epsembr
            d1 = sqrt(muembl * epsembl)
            d2 = sqrt(muembr * epsembr)
            if(aimag(d1)/=0.d0) then
                write(6, 227)
                stop
            endif
            if(aimag(d2)/=0.d0) then
                write(6, 228)
                stop
            endif
        endif
    else
        read(10, *) dummy, (al(i), i = 1, 3)
    endif
      !    multipole decomposition section -----------------------------------
    call cli_parse
    call ini_parse
    is_multipole_type_selected = mrp%is_multipole_type_selected
    is_multipole_order_selected = mrp%is_multipole_order_selected
    is_m_projection_selected = mrp%is_m_projection_selected
    multipole_type = mrp%multipole_type
    multipole_order = mrp%multipole_order
    m_projection = mrp%m_projection

    multipole_combination = get_multipole_combination(lmax, multipole_type, multipole_order, m_projection,&
                            is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected)
    !     ------------------------------------------------------------------
    
    call main_evaluate(ncompd, npland, lmax, i, ktype, kscan, ncomp, np,&
            nunit, icomp, kemb,  ipl, alpha, rmax, zinf, zsup, fab, alphap, theta,&
            fi, fein, d2, d1, polar, &
            it, nlayer, nplan, dl, dr, s, al, d, aq, eps2, eps3, mu1, mu2, mu3,&
            eps1, musph, epssph, multipole_combination)
    stop
    200 format(///, 6(10x, i2))
    201 format(6(10x, i2))
    202 format(4(8x, f12.6))
    203 format(6x, i4, 2(8x, f19.15))
    204 format(2(15x, f19.15), 10x, a2, 10x, f7.2///)
    205 format(2(12x, 2f13.8))
    206 format(10x, f13.8, 2(12x, 2f13.8))
    207 format(3x, 'k_parallel=', 2f12.6, 5x, a2, 'polarization')
    208 format(3x, 'angles of incidence (in rad):  theta=', f7.2, 3x, 'fi=', &
            f7.2, 5x, a2, 'polarization')
    209 format(3x, 'component nr.', i2, 3x, 'type:', 2x, a17)
    210 format(3x, 'mu :', 2f10.5, ' | ', 2f10.5, ' | ', 2f10.5/&
            3x, 'eps:', 2f10.5, ' | ', 2f10.5, ' | ', 2f10.5)
    211 format(3x, 'mu :', 2f10.5, ' | ', 4(3x, 2f10.5))
    212 format(29x, i6, ' unit layers')
    220 format(3x, 'eps:', 2f10.5, ' | ', 4(3x, 2f10.5))
    221 format(3x, 'the sample consists of ', i6, ' unit slices')
    222 format(5x, '****************************************************'/&
            5x, '*** output: transmittance/reflectance/absorbance ***'/&
            5x, '****************************************************')
    223 format(5x, '****************************************************'/&
            5x, '************** output: band structure **************'/&
            5x, '****************************************************')
    224 format(3x, '  s:', 23x, 4(3x, f10.5, 10x))
    225 format(3x, 'k_parallel=', 2f12.6)
    226 format(5x, '----------------------------------'/&
            5x, 'illegal input: spheres embedded in'/&
            5x, 'a  medium  of negative  dielectric'/&
            5x, 'constant. the ewald  summation  in'/&
            5x, 'subroutine xmat does not converge.'/&
            5x, 'direct - space summation is needed'/&
            5x, 'instead.'/&
            5x, '----------------------------------')
    227 format(5x, '----------------------------------'/&
            5x, 'illegal input:semi-infinite medium'/&
            5x, 'of complex refractive index on the'/&
            5x, 'left side of the slab.'/&
            5x, '----------------------------------')
    228 format(5x, '-------------- --------------------'/&
            5x, 'illegal input:semi-infinite medium'/&
            5x, 'of complex refractive index on the'/&
            5x, 'right side of the slab.'/&
            5x, '----------------------------------')
end program