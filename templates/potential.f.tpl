C***********************************************************************
C
C
                          SUBROUTINE potparams
C
C
C***********************************************************************
C
C
C     Subroutine to read in potential parameters.
C
C     Input parameters:
C
C
C
C=======================================================================

      INCLUDE 'scf.h'
      INCLUDE 'potential.h'

      CHARACTER *1 pcomment
      CHARACTER*8 filepar,filename
C=======================================================================


      filepar = potparsfile
      filename = filepar(1:6)

      OPEN(UNIT=upotpars,FILE=filename,STATUS='OLD')

C     Read parameters, close the file.
C     --------------------------------

{readblock:s}

      CLOSE(UNIT=upotpars)

      RETURN
      END

C***********************************************************************
C
C
                          SUBROUTINE accp_external
C
C
C***********************************************************************
C Triaxial Log halo
C
      INCLUDE 'scf.h'
      INCLUDE 'potential.h'

C   Declaration of local variables.
C   -------------------------------
      INTEGER i
      REAL*8 xx,yy,zz,
     &       gee,msun,cmpkpc,secperyr

      PARAMETER(gee=6.673840d-8,
     &          msun=1.9890999999999997d33,
     &          cmpkpc=3.085677581467192d21,
     &          secperyr=3.15576d7)

      LOGICAL firstc
      DATA firstc/.TRUE./
      SAVE firstc,tu
{saveblock:s}

C Convert parameters to appropriate units....
      IF(firstc)THEN
            firstc = .FALSE.

C           Simulation units
            tu = SQRT((cmpkpc*ru)**3/(msun*mu*gee))
            vu = (ru*cmpkpc*1.e-5)/tu
            WRITE(6,*)ru,tu/secperyr,vu

            strength = 1.
            IF(ntide.GT.0) THEN
                  strength=0.
C
                  xframe = xframe/ru
                  yframe = yframe/ru
                  zframe = zframe/ru
                  vxframe = vxframe/vu
                  vyframe = vyframe/vu
                  vzframe = vzframe/vu
            END IF
C
            DO 5 i=1,nbodies
              vx(i)=vx(i)+vxframe
              vy(i)=vy(i)+vyframe
              vz(i)=vz(i)+vzframe
 5          CONTINUE

 {firstc:s}

      END IF

      DO 10 i=1,nbodies

      xx = (x(i)+xframe)
      yy = (y(i)+yframe)
      zz = (z(i)+zframe)

{pot_acc:s}

 10     CONTINUE

      RETURN
      END
