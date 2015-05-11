# coding: utf-8

__author__ = "adrn <adrn@astro.columbia.edu>"

__all__ = ['MiyamotoNagaiPotential', 'HernquistPotential', 'LogarithmicPotential',
           'LeeSutoTriaxialNFWPotential']

def MiyamotoNagaiPotential():
    hblock = """
      REAL*8 mn_a,mn_b,mn_M
      COMMON/miyamotonagai/mn_a,mn_b,mn_M
    """

    readblock = """
      READ(upotpars,'(a)') pcomment
      READ(upotpars,*) mn_a
      READ(upotpars,*) mn_b
      READ(upotpars,*) mn_M
    """

    saveblock = """
      REAL*8 a,b,GM,
     &       sqz2b2,tdr,tdz,phim
      SAVE a,b,GM
    """

    firstc = """
C           Miyamoto-Nagai parameters
            a = mn_a/ru
            b = mn_b/ru
            GM = mn_M/mu
    """

    pot_acc = """
C     Compute potential, acceleration for Miyamoto-Nagai
      sqz2b2 = sqrt(z2 + b*b)
      tdr = GM/(r2 + (a + sqz2b2)**2)**1.5
      tdz = tdr*(a/sqz2b2 + 1.)
      phim = -GM/sqrt(r2+(a+sqz2b2)**2)
      ax(i) = ax(i) - strength*tdr*xx
      ay(i) = ay(i) - strength*tdr*yy
      az(i) = az(i) - strength*tdz*zz
      potext(i) = potext(i) + strength*phim
    """

    SCFPOT = """--- MiyamotoNagaiPotential
6.5                   a [kpc]
0.26                  b [kpc]
100000000000.         mass scale [Msun]
"""

    return dict(hblock=hblock, readblock=readblock, saveblock=saveblock,
                firstc=firstc, pot_acc=pot_acc, SCFPOT=SCFPOT)

def HernquistPotential():
    hblock = """
      REAL*8 hq_c,hq_M
      COMMON/hernquist/hq_c,hq_M
    """

    readblock = """
      READ(upotpars,'(a)') pcomment
      READ(upotpars,*) hq_c
      READ(upotpars,*) hq_M
    """

    saveblock = """
      REAL*8 GMs,hs,
     &       r2,z2,rad,tsrad,phis
      SAVE GMs, hs
    """

    firstc = """
C           Hernquist parameters
            GMs = hq_M/mu
            hs = hq_c/ru
    """

    pot_acc = """
C     Compute potential, acceleration due to spheroid
      r2 = xx*xx + yy*yy
      z2 = zz*zz
      rad = sqrt(r2+z2)
      tsrad = GMs/(rad+hs)**2/rad
      phis = -GMs/(rad+hs)
      ax(i) = ax(i) - strength*tsrad*xx
      ay(i) = ay(i) - strength*tsrad*yy
      az(i) = az(i) - strength*tsrad*zz
      potext(i) = potext(i) + strength*phis
    """

    SCFPOT = """--- HernquistPotential
0.3                   c [kpc]
20000000000.          mass scale [Msun]
"""

    return dict(hblock=hblock, readblock=readblock, saveblock=saveblock,
                firstc=firstc, pot_acc=pot_acc, SCFPOT=SCFPOT)

def LogarithmicPotential():
    hblock = """
      REAL*8 log_rh,log_vc,q1,q2,q3,phi,theta,psi
      COMMON/log/log_rh,log_vc,q1,q2,q3,phi,theta,psi
    """

    readblock = """
      READ(upotpars,'(a)') pcomment
      READ(upotpars,*) log_rh
      READ(upotpars,*) log_vc
      READ(upotpars,*) q1
      READ(upotpars,*) q2
      READ(upotpars,*) q3
      READ(upotpars,*) phi
      READ(upotpars,*) theta
      READ(upotpars,*) psi
    """

    saveblock = """
      REAL*8 xr,yr,zr,axl,ayl,azl,fac,phi0,rh2,log_Ms,
     &       ax_h,ay_h,az_h
      SAVE R11,R12,R13,R21,R22,R23,R31,R32,R33,phi0,rh2
    """

    firstc = """
C           Logarithmic potential parameters
            R11 = -sin(phi)*sin(psi)*cos(theta) + cos(phi)*cos(psi)
            R12 = sin(phi)*cos(psi) + sin(psi)*cos(phi)*cos(theta)
            R13 = sin(psi)*sin(theta)
            R21 = -sin(phi)*cos(psi)*cos(theta) - sin(psi)*cos(phi)
            R22 = -sin(phi)*sin(psi) + cos(phi)*cos(psi)*cos(theta)
            R23 = sin(theta)*cos(psi)
            R31 = sin(phi)*sin(theta)
            R32 = -sin(theta)*cos(phi)
            R33 = cos(theta)

C           Factor to go from (km/s)^2 to correct units
            fac = 232443.89128800036
            log_Ms = fac * log_vc * log_vc * log_rh
            phi0 = log_Ms / mu / (log_rh / ru)
            rh2 = log_rh*log_rh / ru / ru
    """

    pot_acc = """
C     Compute potential, acceleration due to Logarithmic potential
      xr = R11*xx + R12*yy + R13*zz
      yr = R21*xx + R22*yy + R23*zz
      zr = R31*xx + R32*yy + R33*zz

      fac = DLOG(rh2 +
     &           xr*xr/q1/q1 +
     &           yr*yr/q2/q2 +
     &           zr*zr/q3/q3)
      potext(i) = potext(i) + strength * phi0 * fac

      fac = -phi0 / (rh2 + xr*xr/q1/q1 +
     &               yr*yr/q2/q2 + zr*zr/q3/q3)
      ax_h = fac*xr / (q1*q1)
      ay_h = fac*yr / (q2*q2)
      az_h = fac*zr / (q3*q3)

      axl = R11*ax_h + R21*ay_h + R31*az_h
      ayl = R12*ax_h + R22*ay_h + R32*az_h
      azl = R13*ax_h + R23*ay_h + R33*az_h

      ax(i) = ax(i) + strength*axl
      ay(i) = ay(i) + strength*ayl
      az(i) = az(i) + strength*azl

    """

    SCFPOT = """--- LogarithmicPotential
20.                   rh (scale radius) [kpc]
220.                  vc (circular velocity) [km/s]
1.                    q1 (major axis)
0.8                   q2 (intermediate axis)
0.6                   q3 (minor axis)
0                     phi (to misalign halo) [radian]
0                     theta (to misalign halo) [radian]
0                     psi (to misalign halo) [radian]
"""

    return dict(hblock=hblock, readblock=readblock, saveblock=saveblock,
                firstc=firstc, pot_acc=pot_acc, SCFPOT=SCFPOT)

def LeeSutoTriaxialNFWPotential():
    hblock = """
      REAL*8 nfw_rs,nfw_vc,nfw_a,nfw_b,nfw_c,phi,theta,psi
      COMMON/nfw/nfw_rs,nfw_vc,nfw_a,nfw_b,nfw_c,phi,theta,psi
    """

    readblock = """
      READ(upotpars,'(a)') pcomment
      READ(upotpars,*) nfw_rs
      READ(upotpars,*) nfw_vc
      READ(upotpars,*) nfw_a
      READ(upotpars,*) nfw_b
      READ(upotpars,*) nfw_c
      READ(upotpars,*) phi
      READ(upotpars,*) theta
      READ(upotpars,*) psi
    """

    saveblock = """
       REAL*8 xx,yy,zz,xf,yf,zf
       REAL*8 phi,theta,psi
       REAL*8 eb2,ec2,phi0
       REAL*8 phih,p,r,yor2,zor2,axh,ayh,azh,ax_h,ay_h,az_h
       REAL*8 coa,boa,cob
       REAL*8 R11,R12,R13,R21,R22,R23,R31,R32,R33
       REAL*8 x0,x1,x2,x4,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,
      &       x16,x17,x18,x19,x20,x21,x22

       SAVE boa,coa,eb2,ec2,c2,phi0,
      &     R11,R12,R13,R21,R22,R23,R31,R32,R33
    """

    # TODO: phi0 is wrong below, relative to python definition
    firstc = """
C           Logarithmic potential parameters
            boa = nfw_b/nfw_a
            coa = nfw_c/nfw_a

            eb2 = 1.0d0 - boa*boa
            ec2 = 1.0d0 - coa*coa

            R11 = -sin(phi)*sin(psi)*cos(theta) + cos(phi)*cos(psi)
            R12 = sin(phi)*cos(psi) + sin(psi)*cos(phi)*cos(theta)
            R13 = sin(psi)*sin(theta)
            R21 = -sin(phi)*cos(psi)*cos(theta) - sin(psi)*cos(phi)
            R22 = -sin(phi)*sin(psi) + cos(phi)*cos(psi)*cos(theta)
            R23 = sin(theta)*cos(psi)
            R31 = sin(phi)*sin(theta)
            R32 = -sin(theta)*cos(phi)
            R33 = cos(theta)

            phi0 = 0.5 * nfw_vc * nfw_vc
    """

    pot_acc = """
C     Compute potential, acceleration due to Logarithmic potential
      p = r/nfw_rs
      yor2 = yf*yf/rad/rad
      zor2 = zf*zf/rad/rad

      xx = R11*xf + R12*yf + R13*zf
      yy = R21*xf + R22*yf + R23*zf
      zz = R31*xf + R32*yf + R33*zf

      x0 = r + nfw_rs
      x1 = x0**2
      x2 = phi0/(12*r**7*x1)
      x4 = yy*yy
      x6 = zz*zz
      x7 = eb2*x4 + ec2*x6
      x8 = nfw_rs**2
      x9 = 6*x8
      x10 = DLOG(x0/nfw_rs)
      x11 = x0*x10
      x12 = 3*nfw_rs
      x13 = r*x12
      x14 = xx*xx
      x15 = x13 - x14 - x4 - x6
      x16 = x15 + x9
      x17 = 6*nfw_rs*x0*(r*x16 - x11*x9)
      x18 = x1*x10
      x19 = r**2
      x20 = x0*x19
      x21 = 2*r*x0
      x22 = -12*r**5*nfw_rs*x0 + 12*x19*x19*nfw_rs*x18 +
     & x12*x7*(x16*x19 - 18*x18*x8 + x20*(2*r - x12) + x21*
     & (x15 + 9*x8)) - x20*(eb2 + ec2)*(-6*r*nfw_rs*(x19 - x8) +
     & 6*nfw_rs*x11*(x19 - 3*x8) + x20*(-4*r + x12) + x21*
     & (-x13 + 2*x14 + 2*x4 + 2*x6 + x9))

      phih = phi0*((eb2/2 + ec2/2)*((1/p - 1/p**3)*x10 - 1 + (2*x19/x8
     &      - 3*p + 6)/(6*x19/x8)) + (eb2*x4/(2*x19) + ec2*x6/(2*x19))
     &      *((x19/x8 - 3*p - 6)/(2*x19/x8*(p + 1)) + 3*x10/p**3) -
     &      x10/p)

      ax_h = -x2*xx*(x17*x7 + x22)
      ay_h = -x2*yy*(x17*(-x19*eb2 + x7) + x22)
      az_h = -x2*zz*(x17*(-x19*ec2 + x7) + x22)
      axh = R11*ax_h + R21*ay_h + R31*az_h
      ayh = R12*ax_h + R22*ay_h + R32*az_h
      azh = R13*ax_h + R23*ay_h + R33*az_h

      ax(i) = ax(i) + strength*axh
      ay(i) = ay(i) + strength*ayh
      az(i) = az(i) + strength*azh

    """

    SCFPOT = """--- LeeSutoTriaxialNFWPotential
20.                   rh (scale radius) [kpc]
220.                  vc (circular velocity) [km/s]
1.                    a (major axis)
0.8                   b (intermediate axis)
0.6                   c (minor axis)
0                     phi (to misalign halo) [radian]
0                     theta (to misalign halo) [radian]
0                     psi (to misalign halo) [radian]
"""

    return dict(hblock=hblock, readblock=readblock, saveblock=saveblock,
                firstc=firstc, pot_acc=pot_acc, SCFPOT=SCFPOT)
