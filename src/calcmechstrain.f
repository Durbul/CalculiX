!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine calcmechstrain(vkl,vokl,emec,eth,iperturb,nalcon,imat,
     &     xthi,vthj)
!
!     calculates the mechanical strain from the displacement gradients
!     and the thermal stretches (only called if ithermal(1).ne.0
!
      implicit none
!
      integer iperturb(*),nalcon(2,*),imat
!
      real*8 elineng(6),vkl(0:3,3),vokl(3,3),emec(6),eth(6),
     &     wkl(3,3),wokl(3,3),xkl(3,3),xth(3,3),xthi(3,3),vthj
!
      if(iperturb(2).eq.1) then
!
!     nlgeom active:
!     multiplicative decomposition of the total gradient into a
!     mechanical and thermal one
!
!     multiplying with the inverse of the thermal deformation gradient
!     at the end of the increment
!
        if(nalcon(1,imat).le.3) then
!
!         isotropic or orthotropic expansion
!
          wkl(1,1)=(1.d0+vkl(1,1))/(1.d0+eth(1))-1.d0
          wkl(2,2)=(1.d0+vkl(2,2))/(1.d0+eth(2))-1.d0
          wkl(3,3)=(1.d0+vkl(3,3))/(1.d0+eth(3))-1.d0
          wkl(1,2)=(vkl(1,2))/(1.d0+eth(2))
          wkl(1,3)=(vkl(1,3))/(1.d0+eth(3))
          wkl(2,3)=(vkl(2,3))/(1.d0+eth(3))
          wkl(2,1)=(vkl(2,1))/(1.d0+eth(1))
          wkl(3,1)=(vkl(3,1))/(1.d0+eth(1))
          wkl(3,2)=(vkl(3,2))/(1.d0+eth(2))
        else
!
!         anisotropic expansion Fmech=F.(Fth)^{-1}
!
!         Fth
!
          xth(1,1)=eth(1)+1.d0
          xth(2,2)=eth(2)+1.d0
          xth(3,3)=eth(3)+1.d0
          xth(1,2)=eth(4)
          xth(1,3)=eth(5)
          xth(2,3)=eth(6)
!
!         det(Fth)
!
          vthj=xth(1,1)*(xth(2,2)*xth(3,3)-xth(2,3)*xth(2,3))
     &        -xth(1,2)*(xth(1,2)*xth(3,3)-xth(1,3)*xth(2,3))
     &        +xth(1,3)*(xth(1,2)*xth(2,3)-xth(1,3)*xth(2,2))
!
!         Fth^{-1} (inverse of a symmetric matrix is symmetric)
!
          xthi(1,1)=(xth(2,2)*xth(3,3)-xth(2,3)*xth(2,3))/vthj
          xthi(2,2)=(xth(1,1)*xth(3,3)-xth(1,3)*xth(1,3))/vthj
          xthi(3,3)=(xth(1,1)*xth(2,2)-xth(1,2)*xth(1,2))/vthj
          xthi(1,2)=(xth(1,3)*xth(2,3)-xth(1,2)*xth(3,3))/vthj
          xthi(1,3)=(xth(1,2)*xth(2,3)-xth(2,2)*xth(1,3))/vthj
          xthi(2,3)=(xth(1,2)*xth(1,3)-xth(1,1)*xth(2,3))/vthj
          xthi(2,1)=xthi(1,2)
          xthi(3,1)=xthi(1,3)
          xthi(3,2)=xthi(2,3)
!
!         deformation gradient F
!     
          xkl(1,1)=1.d0+wkl(1,1)
          xkl(2,2)=1.d0+wkl(2,2)
          xkl(3,3)=1.d0+wkl(3,3)
          xkl(1,2)=wkl(1,2)
          xkl(1,3)=wkl(1,3)
          xkl(2,3)=wkl(2,2)
          xkl(2,1)=xkl(1,2)
          xkl(3,1)=xkl(1,3)
          xkl(3,2)=xkl(2,3)
!
!         wkl=F.(Fth)^{-1}-1.d0
!
          wkl(1,1)=xkl(1,1)*xthi(1,1)+xkl(1,2)*xthi(1,2)
     &         +xkl(1,3)*xthi(1,3)-1.d0
          wkl(2,2)=xkl(2,1)*xthi(1,2)+xkl(2,2)*xthi(2,2)
     &         +xkl(2,3)*xthi(2,3)-1.d0
          wkl(3,3)=xkl(3,1)*xthi(1,3)+xkl(3,2)*xthi(2,3)
     &         +xkl(3,3)*xthi(3,3)-1.d0
          wkl(1,2)=xkl(1,1)*xthi(1,2)+xkl(1,2)*xthi(2,2)
     &         +xkl(1,3)*xthi(2,3)
          wkl(1,3)=xkl(1,1)*xthi(1,3)+xkl(1,2)*xthi(2,3)
     &         +xkl(1,3)*xthi(3,3)
          wkl(2,3)=xkl(2,1)*xthi(1,3)+xkl(2,2)*xthi(2,3)
     &         +xkl(2,3)*xthi(3,3)
          wkl(2,1)=xkl(2,1)*xthi(1,1)+xkl(2,2)*xthi(1,2)
     &         +xkl(2,3)*xthi(1,3)
          wkl(3,1)=xkl(3,1)*xthi(1,1)+xkl(3,2)*xthi(1,2)
     &         +xkl(3,3)*xthi(1,3)
          wkl(3,2)=xkl(3,1)*xthi(1,2)+xkl(3,2)*xthi(2,2)
     &         +xkl(3,3)*xthi(2,3)
        endif
      else
!
!     nlgeom not active:
!     additive decomposition of the total gradient into a
!     mechanical and thermal one
!
!     subtracting the thermal stretch from the deformation gradients
!     at the end of the increment
!
        wkl(1,1)=vkl(1,1)-eth(1)
        wkl(2,2)=vkl(2,2)-eth(2)
        wkl(3,3)=vkl(3,3)-eth(3)
        wkl(1,2)=vkl(1,2)-eth(4)
        wkl(1,3)=vkl(1,3)-eth(5)
        wkl(2,3)=vkl(2,3)-eth(6)
        wkl(2,1)=vkl(2,1)-eth(4)
        wkl(3,1)=vkl(3,1)-eth(5)
        wkl(3,2)=vkl(3,2)-eth(6)
      endif
!
!     attention! elineng(4),elineng(5) and elineng(6) are engineering strains!
!     
      elineng(1)=wkl(1,1)
      elineng(2)=wkl(2,2)
      elineng(3)=wkl(3,3)
      elineng(4)=wkl(1,2)+wkl(2,1)
      elineng(5)=wkl(1,3)+wkl(3,1)
      elineng(6)=wkl(2,3)+wkl(3,2)
c      write(*,*) 'elineng,wkl ',elineng(1),wkl(1,1)
!     
      if(iperturb(2).eq.1) then
!     
!     Lagrangian strain
!     
         elineng(1)=elineng(1)+
     &        (wkl(1,1)**2+wkl(2,1)**2+wkl(3,1)**2)/2.d0
         elineng(2)=elineng(2)+
     &        (wkl(1,2)**2+wkl(2,2)**2+wkl(3,2)**2)/2.d0
         elineng(3)=elineng(3)+
     &        (wkl(1,3)**2+wkl(2,3)**2+wkl(3,3)**2)/2.d0
         elineng(4)=elineng(4)+wkl(1,1)*wkl(1,2)+wkl(2,1)*wkl(2,2)+
     &        wkl(3,1)*wkl(3,2)
         elineng(5)=elineng(5)+wkl(1,1)*wkl(1,3)+wkl(2,1)*wkl(2,3)+
     &        wkl(3,1)*wkl(3,3)
         elineng(6)=elineng(6)+wkl(1,2)*wkl(1,3)+wkl(2,2)*wkl(2,3)+
     &        wkl(3,2)*wkl(3,3)
!     
!     for frequency analysis or buckling with preload the
!     strains are calculated with respect to the deformed
!     configuration
!     
      elseif(iperturb(1).eq.1) then
!
!     subtracting the thermal stretch from the deformation gradients
!     at the start of the increment
!     
         wokl(1,1)=vokl(1,1)-eth(1)
         wokl(2,2)=vokl(2,2)-eth(2)
         wokl(3,3)=vokl(3,3)-eth(3)
         wokl(1,2)=vokl(1,2)-eth(4)
         wokl(1,3)=vokl(1,3)-eth(5)
         wokl(2,3)=vokl(2,3)-eth(6)
         wokl(2,1)=vokl(2,1)-eth(4)
         wokl(3,1)=vokl(3,1)-eth(5)
         wokl(3,2)=vokl(3,2)-eth(6)
!      
         elineng(1)=elineng(1)+wokl(1,1)*wkl(1,1)+wokl(2,1)*wkl(2,1)+
     &        wokl(3,1)*wkl(3,1)
         elineng(2)=elineng(2)+wokl(1,2)*wkl(1,2)+wokl(2,2)*wkl(2,2)+
     &        wokl(3,2)*wkl(3,2)
         elineng(3)=elineng(3)+wokl(1,3)*wkl(1,3)+wokl(2,3)*wkl(2,3)+
     &        wokl(3,3)*wkl(3,3)
         elineng(4)=elineng(4)+wokl(1,1)*wkl(1,2)+wokl(1,2)*wkl(1,1)+
     &        wokl(2,1)*wkl(2,2)+wokl(2,2)*wkl(2,1)+
     &        wokl(3,1)*wkl(3,2)+wokl(3,2)*wkl(3,1)
         elineng(5)=elineng(5)+wokl(1,1)*wkl(1,3)+wokl(1,3)*wkl(1,1)+
     &        wokl(2,1)*wkl(2,3)+wokl(2,3)*wkl(2,1)+
     &        wokl(3,1)*wkl(3,3)+wokl(3,3)*wkl(3,1)
         elineng(6)=elineng(6)+wokl(1,2)*wkl(1,3)+wokl(1,3)*wkl(1,2)+
     &        wokl(2,2)*wkl(2,3)+wokl(2,3)*wkl(2,2)+
     &        wokl(3,2)*wkl(3,3)+wokl(3,3)*wkl(3,2)
      endif
!     
!     storing the local strains
!     
      if(iperturb(1).ne.-1) then
         emec(1)=elineng(1)
         emec(2)=elineng(2)
         emec(3)=elineng(3)
         emec(4)=elineng(4)/2.d0
         emec(5)=elineng(5)/2.d0
         emec(6)=elineng(6)/2.d0
      else
!     
!        linear iteration within a nonlinear increment:
!
!        subtracting the thermal stretch from the deformation gradients
!        at the start of the increment
!     
         wokl(1,1)=vokl(1,1)-eth(1)
         wokl(2,2)=vokl(2,2)-eth(2)
         wokl(3,3)=vokl(3,3)-eth(3)
         wokl(1,2)=vokl(1,2)-eth(4)
         wokl(1,3)=vokl(1,3)-eth(5)
         wokl(2,3)=vokl(2,3)-eth(6)
         wokl(2,1)=vokl(2,1)-eth(4)
         wokl(3,1)=vokl(3,1)-eth(5)
         wokl(3,2)=vokl(3,2)-eth(6)
!     
         emec(1)=wokl(1,1)+
     &        (wokl(1,1)**2+wokl(2,1)**2+wokl(3,1)**2)/2.d0
         emec(2)=wokl(2,2)+
     &        (wokl(1,2)**2+wokl(2,2)**2+wokl(3,2)**2)/2.d0
         emec(3)=wokl(3,3)+
     &        (wokl(1,3)**2+wokl(2,3)**2+wokl(3,3)**2)/2.d0
         emec(4)=(wokl(1,2)+wokl(2,1)+wokl(1,1)*wokl(1,2)+
     &        wokl(2,1)*wokl(2,2)+wokl(3,1)*wokl(3,2))/2.d0
         emec(5)=(wokl(1,3)+wokl(3,1)+wokl(1,1)*wokl(1,3)+
     &        wokl(2,1)*wokl(2,3)+wokl(3,1)*wokl(3,3))/2.d0
         emec(6)=(wokl(2,3)+wokl(3,2)+wokl(1,2)*wokl(1,3)+
     &        wokl(2,2)*wokl(2,3)+wokl(3,2)*wokl(3,3))/2.d0
      endif
c      write(*,*) 'emec ',emec(1)
!     
      return
      end
