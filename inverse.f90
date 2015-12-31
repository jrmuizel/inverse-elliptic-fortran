program xigel
        implicit real*8 (a-h,j-z)
        external gel,gel1
        real*8 aigel,bigel,nigel
        integer i,icase,ierr
        common /icase/icase
        do i=1,3
                icase=i
                if(icase.EQ.1) then
                        nc=1.d0;mc=0.5d0
                elseif(icase.EQ.2) then
                        nc=2.d0/3.d0;mc=1.d0/3.d0
                elseif(icase.EQ.3) then
                        nc=1.d0;mc=2.d0/3.d0
                endif
                phi=bigel(nc,mc,1.d-15,1.d-15,gel,B,D,J,sn,cn,dn,f)
                write(*,'(a10,1pe25.15)') "bisection:",phi
                kc=sqrt(mc)
                cB=cel(kc,1.d0,1.d0,0.d0,ierr)
                cD=cel(kc,1.d0,0.d0,1.d0,ierr)
                cJ=cel(kc,nc,0.d0,1.d0,ierr)
                u=aigel(nc,mc,cB,cD,cJ,1.d-15,1.d-15,gel,B,D,J,sn,cn,dn,f)
                phi=atan2(sn,cn)
                write(*,'(a10,1p2e25.15)') "accelerated:",phi,u
                u=nigel(nc,mc,cB,cD,cJ,1.d-15,1.d-15,gel,gel1,B,D,J,sn,cn,dn,f)
                phi=atan2(sn,cn)
                write(*,'(a10,1p2e25.15)') "Newton:",phi,u
        enddo
end program xigel

real*8 function gel(u,nc,mc,B,D,J,sn,cn,dn) 
        real*8 nc,mc,u,B,D,J,sn,cn,dn,n,m
        integer icase
        common /icase/icase
        n=1.d0-nc;m=1.d0-mc
        if(icase.EQ.1) then
                gel=B+mc*D-1.d0
        elseif(icase.EQ.2) then
                gel=u+0.5d0*(u+n*J)-1.d0
        elseif(icase.EQ.3) then
                gel=sn*dn/cn-m*B-1.d0
        endif
        return;end
real*8 function gel1(u,nc,mc,B,D,J,sn,cn,dn,B1,D1,J1,sn1,cn1,dn1)
        real*8 nc,mc,u,B,D,J,sn,cn,dn,B1,D1,J1,sn1,cn1,dn1,n,m
        integer icase
        common /icase/icase
        n=1.d0-nc;m=1.d0-mc
        if(icase.EQ.1) then
                gel1=B1+mc*D1
        elseif(icase.EQ.2) then
                gel1=1.d0+0.5d0*(1.d0+n*J1)
        elseif(icase.EQ.3) then
                gel1=((sn1*dn+sn*dn1)*cn-sn*dn*cn1)/(cn*cn)-m*B1
        endif
        return;end


real*8 function nigel(nc,mc,cB,cD,cJ,uTOL,fTOL,gel,gel1,B,D,J,sn,cn,dn,f)
        implicit real*8 (a-h,j-z)
        integer i;external gel;real*8 aigel,gel,gel1
        m=1.d0-mc;n=1.d0-nc;h=n*nc*(mc-nc)
        u=aigel(nc,mc,cB,cD,cJ,0.2d0,fTOL,gel,B,D,J,sn,cn,dn,f)
        x=cn*cn;y=sn*sn;z=dn*dn
        do i=1,10
                D1=y;B1=x;J1=y/(nc+n*x);cn1=-sn*dn;sn1=cn*dn;dn1= -m*sn*cn
                f1=gel1(u,nc,mc,B,D,J,sn,cn,dn,B1,D1,J1,sn1,cn1,dn1)
                v= -f/f1;u=u+v
                call sersdj(v,n,m,snv,Dv,Jv)
                yv=snv*snv;xv=1.d0-yv;zv=1.d0-m*yv;cnv=sqrt(xv);dnv=sqrt(zv)
                Bv=v-Dv;xi=cn*cnv;eta=sn*snv;zeta=dn*dnv
                nu=1.d0/(1.d0-m*y*yv);cn=(xi-eta*zeta)*nu
                sn=(sn*cnv*dnv+snv*sn1)*nu;dn=(zeta-m*eta*xi)*nu
                x=cn*cn;y=sn*sn;z=dn*dn;W=eta*sn;B=B+Bv-W;D=D+Dv+W
                t=W/(1.d0-n*(y-eta*cn*dn));J=J+Jv+uatan(t,h)
                f=gel(u,nc,mc,B,D,J,sn,cn,dn)
                if(v*v.LT.uTOL.or.abs(f).LT.fTOL) then
                        nigel=u;return
                endif
                enddo
        write(*,*) "(nigel) No convergence"
        return; end

real*8 function aigel(nc,mc,cB,cD,cJ,uTOL,fTOL,gel,B,D,J,sn,cn,dn,f)
        implicit real*8 (a-h,j-z)
        integer i;real*8 gel
        n=1.d0-nc;kc=sqrt(mc);m=1.d0-mc;h=n*nc*(mc-nc)
        uL=0.d0;BL=0.d0;DL=0.d0;JL=0.d0;snL=0.d0;cnL=1.d0;dnL=1.d0
        yL=0.d0;fL=gel(uL,nc,mc,BL,DL,JL,snL,cnL,dnL)
        uU=cB+cD;BU=cB;DU=cD;JU=cJ;snU=1.d0;cnU=0.d0;dnU=kc
        yU=1.d0;fU=gel(uU,nc,mc,BU,DU,JU,snU,cnU,dnU)
        uH=uU;BH=BU;DH=DU;JH=JU;snH=snU;cnH=cnU;dnH=dnU
        yH=yU;fH=fU
        do i=1,60
                uH=0.5d0*uH;y=yH;v=cnH*dnH;p=cnH+dnH;q=1.d0/(1.d0+cnH)
                r=1.d0/(1.d0+dnH);xH=p*r;yH=yH*q*r;zH=p*q;W=yH*snH
                cnH=sqrt(xH);snH=sqrt(yH);dnH=sqrt(zH);BH=0.5d0*(BH+W)
                DH=0.5d0*(DH-W);t=W/(1.d0-n*(y-yH*v))
                JH=0.5d0*(JH-uatan(t,h));u=uL+uH;xi=cnL*cnH;eta=snL*snH
                zeta=dnL*dnH;nu=1.d0/(1.d0-m*yL*yH);cn=(xi-eta*zeta)*nu
                sn=(snL*cnH*dnH+snH*cnL*dnL)*nu;dn=(zeta-m*eta*xi)*nu
                x=cn*cn;y=sn*sn;z=dn*dn;W=eta*sn;B=BL+BH-W;D=DL+DH+W
                t=W/(1.d0-n*(y-eta*cn*dn));J=JL+JH+uatan(t,h)
                f=gel(u,nc,mc,B,D,J,sn,cn,dn)
                if(f.LT.0.d0) then
                        uL=u;BL=B;DL=D;JL=J;snL=sn;cnL=cn;dnL=dn;yL=y;fL=f
                else
                        uU=u;BU=B;DU=D;JU=J;snU=sn;cnU=cn;dnU=dn;yU=y;fU=f
                endif
                if(abs(uU-uL).LT.uTOL.or.abs(f).lt.fTOL) then
                        if(abs(fU).LT.abs(fL)) then
                                u=uU;B=BU;D=DU;J=JU;sn=snU;cn=cnU;dn=dnU;f=fU
                        else
                                u=uL;B=BL;D=DL;J=JL;sn=snL;cn=cnL;dn=dnL;f=fL
                        endif
                        aigel=u;return
                endif
        enddo
        write(*,*) "(aigel) No convergence"
        return;end


real*8 function uatan(t,h)
        real*8 t,h,z,r,y,A3,A5,A7,A9
        parameter (A3=1.d0/3.d0,A5=1.d0/5.d0,A7=1.d0/7.d0,A9=1.d0/9.d0)
        z= -h*t*t
        if(abs(z).lt.1.d-3) then
                uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*A9))))
        elseif(z.lt.0.d0) then
                r=sqrt(h);uatan=atan(r*t)/r
        else
                r=sqrt(-h);y=r*t;uatan=log((1.d0+y)/(1.d0-y))*0.5d0/r
        endif
        return;end

subroutine sersdj(u,n,m,sn,D,J)
        implicit real*8 (a-z)
        parameter (A1=1.d0/6.d0,A2=1.d0/120.d0,A3=1.d0/5040.d0)
        parameter (B1=1.d0/3.d0,B2=1.d0/15.d0,B3=1.d0/315.d0)
        parameter (A4=1.d0/362880.d0,A5=1.d0/39916800.d0)
        parameter (B4=1.d0/2835.d0,B5=1.d0/155925.d0)
        parameter (A6=1.d0/6227020800.d0,B6=1.d0/6080175.d0)
        n2=n*n;n3=n2*n;n4=n2*n2;m2=m*m;m3=m2*m;m4=m2*m2
        m5=m3*m2;m6=m3*m3;mp=1.d0+m;m2p=1.d0+m2;m3p=1.d0+m3
        m4p=1.d0+m4;m5p=1.d0+m5;m6p=1.d0+m6;mm2=m+m2
        mm3=m+m3;mm4=m+m4;mm5=m+m5;m2m3=m2+m3;m2m4=m2+m4
        S1=A1*mp;S2=A2*(m2p+14.d0*m);S3=A3*(m3p+135.d0*mm2)
        S4=A4*(m4p+1228.d0*mm3+5478.d0*m2)
        S5=A5*(m5p+11069.d0*mm4+165826.d0*m2m3)
        S6=A6*(m6p+99642.d0*mm5+4494351.d0*m2m4+13180268.d0*m3)
        D2=B2*mp;D3=B3*(2.d0*m2p+13.d0*m);D4=B4*(m3p+30.d0*mm2)
        D5=B5*(2.d0*m4p+251.d0*mm3+876.d0*m2)
        D6=B6*(2.d0*m5p+1018.d0*mm4+9902.d0*m2m3)
        J2=B2*3.d0;J3=B3*(30.d0*mp-45.d0*n)
        J4=B4*(63.d0*m2p+252.d0*m-315.d0*mp*n+315.d0*n2)
        J5=B5*(510.d0*m3p+5850.d0*mm2-(6615.d0*m2p+21735.d0*m)*n+18900.d0*mp*n2-14175.d0*n3)
        J6=B6*(2046.d0*m4p+59268.d0*mm3+158103.d0*m2-(63360.d0*m3p+497475.d0*mm2)*n+(395010.d0*m2p+ &
                1164240.d0*m)*n2-779625.d0*mp*n3+467775.d0*n4)
        u2=u*u;u3=u2*u;u5=u3*u2
        sn=u*(1.d0-u2*(S1-u2*(S2-u2*(S3-u2*(S4-u2*(S5-u2*S6))))))
        D=u3*(B1-u2*(D2-u2*(D3-u2*(D4-u2*(D5-u2*D6)))))
        J=D+n*u5*(J2-u2*(J3-u2*(J4-u2*(J5-u2*J6))))
        return;end

real*8 function bigel(nc,mc,pTOL,fTOL,gel,B,D,J,sn,cn,dn,f)
        implicit real*8 (a-h,j-z)
        integer i
        real*8 gel
        parameter (PI=3.1415926535897932d0)
        dp=0.5d0*PI;p=dp;p0=0.d0
        B0=0.d0;D0=0.d0;J0=0.d0;sn0=0.d0;cn0=1.d0;dn0=1.d0
        f0=gel(0.d0,nc,mc,B0,D0,J0,sn0,cn0,dn0)
        do i=1,60
                call elbdj(p,nc,mc,B,D,J)
                u=B+D;sn=sin(p);cn=cos(p);dn=sqrt(cn*cn+mc*sn*sn)
                f=gel(u,nc,mc,B,D,J,sn,cn,dn)
                if(abs(p-p0).LT.pTOL.or.abs(f).lt.fTOL) then
                        if(abs(f).GT.abs(f0)) then
                                p=p0;B=B0;D=D0;J=J0;sn=sn0;cn=cn0;dn=dn0;f=f0
                        endif
                        bigel=p;return
                endif
                p0=p;B0=B;D0=D;J0=J;sn0=sn;cn0=cn;dn0=dn;f0=f;dp=dp*0.5d0
                if(f.LT.0.d0) then
                        p=p+dp
                else
                        p=p-dp
                endif
        enddo
        write(*,*) "(bigel) No convergence"
        return;end


subroutine elbdj(phi,nc,mc,B,D,J)
        real*8 phi,nc,mc,B,D,J,sn,sn2,cn2,dn2,sn33,F,rf,rd,rj 
        sn=sin(phi);sn2=sn*sn;cn2=1.d0-sn2;dn2=cn2+mc*sn2;sn33=sn2*sn/3.d0
        !write(*,'(a10,1pe25.15,1pe25.15)') "rf:",cn2,dn2
        F=sn*rf(cn2,dn2,1.d0)
        D=sn33*rd(cn2,dn2,1.d0)
        J=sn33*rj(cn2,dn2,1.d0,cn2+nc*sn2)
        B=F-D
        return;end
