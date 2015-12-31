      FUNCTION rf(x,y,z)
      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08,TINY=1.5e-38,BIG=3.E37,THIRD=1./3.,
     *    C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
      !Computes Carlson’s elliptic integral of the first kind, RF(x,y,z).
      !Cx, y, and z must be nonnegative, and at most one can be zero.
      !CTINY must be at least 5 times the machine underflow limit, BIG at
      !Cmost one fifth the machine overflow limit.
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0)pause 'invalid arguments in rf1'
      if(min(x+y,x+z,y+z).lt.TINY)pause 'invalid arguments in rf2'
      if(max(x,y,z).gt.BIG)pause 'invalid arguments in rf3'
      xt=x
      yt=y
      zt=z
1     continue
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=.25*(xt+alamb)
         yt=.25*(yt+alamb)
         zt=.25*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END


      FUNCTION rd(x,y,z)
      REAL*8 rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=.05,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,
     * C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)
Computes Carlson’s elliptic integral of the second kind,
CRD(x,y,z). x and y must be nonnegative, and at most one can be
Czero. z must be positive. TINY must be at least twice the
Cnegative 2/3 power of the machine overflow limit. BIG must be at
Cmost 0.1 × ERRTOL times the negative 2/3 power of the machine
Cunderflow limit.
      REAL*8 alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
     sqrtz,sum,xt,yt,zt
      if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.
     * max(x,y,z).gt.BIG)pause 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.
      fac=1.
1     continue
             sqrtx=sqrt(xt)
             sqrty=sqrt(yt)
             sqrtz=sqrt(zt)
             alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
             sum=sum+fac/(sqrtz*(zt+alamb))
             fac=.25*fac
             xt=.25*(xt+alamb)
             yt=.25*(yt+alamb)
             zt=.25*(zt+alamb)
             ave=.2*(xt+yt+3.*zt)
             delx=(ave-xt)/ave
             dely=(ave-yt)/ave
             delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.*eb
      ee=ed+ec+ec
      rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)
     *        + delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
      return
      END


      FUNCTION rj(x,y,z,p)
      REAL*8 rj,p,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (ERRTOL=.05,TINY=2.5e-13,BIG=9.E11,C1=3./14.,C2=1./3.,
     *    C3=3./22.,C4=3./26.,C5=.75*C3,C6=1.5*C4,C7=.5*C2,C8=C3+C3)
C              USES rc,rf
C              Computes Carlson’s elliptic integral of the third kind, RJ (x, y, z, p). x, y, and z must be nonnegative, and at most one can be zero. p must be nonzero. If p < 0, the Cauchy principal value is returned. TINY must be at least twice the cube root of the machine underflow limit, BIG at most one fifth the cube root of the machine overflow limit.
      REAL*8 a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,
     *     ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,
     *     yt,zt,rc,rf
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.
     *     max(x,y,z,abs(p)).gt.BIG)pause 'invalid arguments in rj'
      sum=0.
      fac=1.
      if(p.gt.0.)then
              xt=x
              yt=y
              zt=z
              pt=p
      else
              xt=min(x,y,z)
              zt=max(x,y,z)
              yt=x+y+z-xt-zt
              a=1./(yt-p)
              b=a*(zt-yt)*(yt-xt)
              pt=yt+b
              rho=xt*zt/yt
              tau=p*pt/yt
              rcx=rc(rho,tau)
      endif
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
        beta=pt*(pt+alamb)**2
        sum=sum+fac*rc(alpha,beta)
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        pt=.25*(pt+alamb)
        ave=.2*(xt+yt+zt+pt+pt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        delp=(ave-pt)/ave
      if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.*ec
      ee=eb+2.*delp*(ea-ec)
      rj=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
     *     +delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
      if (p.le.0.) rj=a*(b*rj+3.*(rcx-rf(xt,yt,zt)))
      return
      END

      FUNCTION rc(x,y)
      REAL*8 rc,x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,
     *     C1,C2,C3,C4
      PARAMETER (ERRTOL=.04,TINY=1.69e-38,SQRTNY=1.3e-19,BIG=3.E37,
     * TNBG=TINY*BIG,COMP1=2.236/SQRTNY,COMP2=TNBG*TNBG/25.,
     * THIRD=1./3.,C1=.3,C2=1./7.,C3=.375,C4=9./22.)
!Computes Carlson’s degenerate elliptic integral, RC(x,y). x must be nonnegative and y must be nonzero. If y < 0, the Cauchy principal value is returned. TINY must be at least 5 times the machine underflow limit, BIG at most one fifth the machine maximum overflow limit.
      REAL*8 alamb,ave,s,w,xt,yt
      if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.TINY.or.(x+abs(y)).gt.BIG
     *    .or.(y.lt.-COMP1.and.x.gt.0..and.x.lt.COMP2))
     *    pause 'invalid arguments in rc'
      if(y.gt.0.)then
        xt=x
        yt=y
        w=1.
      else xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
1     continue
        alamb=2.*sqrt(xt)*sqrt(yt)+yt
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        ave=THIRD*(xt+yt+yt)
        s=(yt-ave)/ave
      if(abs(s).gt.ERRTOL)goto 1
      rc=w*(1.+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END
