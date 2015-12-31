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


