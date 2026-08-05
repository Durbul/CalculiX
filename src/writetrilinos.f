!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine writetrilinos(jq,irow,ad,au,b,neq,nzs,symmetryflag,
     &     inputformat,co,nk,nactdof,jobnamef,mi,idofstart,ndofpernode,
     &     idofoffset,iperm,nfull)
!
      implicit none
!
      character*132 jobnamef,fnlhs,fnrhs,fnrig
!
!     writes the assembled system to <jobnamef>.lhs/.rhs/.rig for the
!     external (Trilinos) solver
!
!     idofstart/ndofpernode select which contiguous slice of
!     nactdof(0:mi(2),node) (0=temperature, 1-3=translations) to
!     expand: (1,3)=mechanical, (0,1)=thermal, (0,mi(2)+1)=coupled.
!     ndofpernode also sets MueLu's amalgamation block size
!
!     idofoffset shifts dof numbers to a local 1..neq range: ccx
!     numbers the mechanical and thermal dofs together as one range
!     split at neq[0], so the thermal-only block needs
!     idofoffset=neq[0] (0 otherwise)
!
!     returns iperm(1:neq), mapping an active dof to its row/column in
!     the expanded nfull x nfull system, and nfull itself, so the
!     caller can un-permute the external solver's solution
!
      integer jq(*),irow(*),neq,nzs,symmetryflag,inputformat,nk,
     &     mi(*),nactdof(0:mi(2),*),idofstart,ndofpernode,idofoffset,
     &     iperm(neq),full2active(ndofpernode*nk),
     &     i,j,k,nnz,node,dir,idof,nfull
!
      real*8 ad(*),au(*),b(*),co(3,*),rigid(ndofpernode*nk,6),cg(3)
!
      fnlhs=jobnamef(1:index(jobnamef,' ')-1)//'.lhs'
      fnrhs=jobnamef(1:index(jobnamef,' ')-1)//'.rhs'
      fnrig=jobnamef(1:index(jobnamef,' ')-1)//'.rig'
!
!     ccx's active-dof numbering has SPC'd/MPC-condensed dofs removed,
!     so nodes don't have a fixed dof count -- incompatible with
!     MueLu's blocked-AMG amalgamation. Re-expand into a node-major
!     numbering (1..nfull) with exactly ndofpernode slots per active
!     node: iperm maps an active dof to its slot, full2active maps
!     back (0 = constrained/condensed out). Constrained slots get an
!     identity row (diagonal 1, rhs 0, no coupling) so they solve to 0
!     without affecting the physical solution
!
      k=0
      full2active(1:ndofpernode*nk)=0
      do i=1,nk
        if(all(nactdof(idofstart:idofstart+ndofpernode-1,i)
     &       .le.idofoffset)) cycle
        do j=idofstart,idofstart+ndofpernode-1
          k=k+1
          if(nactdof(j,i).gt.idofoffset) then
            iperm(nactdof(j,i)-idofoffset)=k
            full2active(k)=nactdof(j,i)-idofoffset
          endif
        enddo
      enddo
      nfull=k
!
!     left hand side
!
      open(50, file=fnlhs, status='unknown')
      if(symmetryflag.eq.0) then
        nnz=nfull+jq(neq+1)-1
        write(50,'(a)')
     &   '%%MatrixMarket matrix coordinate real symmetric'
        write(50,*) nfull,nfull,nnz
        do k=1,nfull
          if(full2active(k).gt.0) then
            write(50,*) k,k,ad(full2active(k))
          else
            write(50,*) k,k,1.d0
          endif
        enddo
        do i=1,neq
          do j=jq(i),jq(i+1)-1
!           irow(j) is always > i (add_sm_st.f stores each coefficient
!           under the smaller of its two dof indices), so au holds the
!           lower triangle; emit as (row=irow(j),col=i) per MM 'symmetric'
            write(50,*) iperm(irow(j)),iperm(i),au(j)
          enddo
        enddo
      elseif(inputformat.eq.1) then
        nnz=nfull+2*(jq(neq+1)-1)
        write(50,'(a)')
     &   '%%MatrixMarket matrix coordinate real general'
        write(50,*) nfull,nfull,nnz
        do k=1,nfull
          if(full2active(k).gt.0) then
            write(50,*) k,k,ad(full2active(k))
          else
            write(50,*) k,k,1.d0
          endif
        enddo
        do i=1,neq
          do j=jq(i),jq(i+1)-1
            write(50,*) iperm(irow(j)),iperm(i),au(j)
            write(50,*) iperm(i),iperm(irow(j)),au(j+nzs)
          enddo
        enddo
      else
        write(*,*) '*ERROR in writetrilinos: input format not known'
        call exit(201)
      endif
      close(50)
!
!     right hand side
!
      open(50,file=fnrhs,status='unknown')
      write(50,'(a)') '%%MatrixMarket matrix array real general'
      write(50,*) nfull,1
      do k=1,nfull
        if(full2active(k).gt.0) then
          write(50,*) b(full2active(k))
        else
          write(50,*) 0.d0
        endif
      enddo
      close(50)
!
!     rigid body modes (nfull x 6 dense, indexed via iperm; constrained
!     slots get an all-zero row). Only meaningful for a pure mechanical
!     block (3 translational dofs/node) -- skip the .rig file otherwise
!     and let the caller omit --nullspace, so MueLu falls back to its
!     default (constant-vector) near-nullspace
!
      if((idofstart.ne.1).or.(ndofpernode.ne.3)) then
        return
      endif
!
      open(50,file=fnrig,status='unknown')
      write(50,'(a)') '%%MatrixMarket matrix array real general'
      write(50,*) nfull,6
!
!     rotation columns (4-6) are taken about the mesh centroid, not
!     the global origin: both span the same nullspace, but centering
!     keeps their magnitude commensurate with the unit-scaled
!     translation columns 1-3 for meshes far from the origin,
!     improving conditioning of MueLu's per-aggregate nullspace QR
!
      cg(1)=sum(co(1,1:nk))/nk
      cg(2)=sum(co(2,1:nk))/nk
      cg(3)=sum(co(3,1:nk))/nk
      rigid(1:nfull,1:6)=0.d0
      do node=1,nk
        do dir=1,3
          idof=nactdof(dir,node)
          if(idof.gt.0) then
            k=iperm(idof)
            if(dir.eq.1) rigid(k,1)=1.d0
            if(dir.eq.2) rigid(k,2)=1.d0
            if(dir.eq.3) rigid(k,3)=1.d0
            if(dir.eq.2) rigid(k,4)=-(co(3,node)-cg(3))
            if(dir.eq.3) rigid(k,4)= (co(2,node)-cg(2))
            if(dir.eq.1) rigid(k,5)= (co(3,node)-cg(3))
            if(dir.eq.3) rigid(k,5)=-(co(1,node)-cg(1))
            if(dir.eq.1) rigid(k,6)=-(co(2,node)-cg(2))
            if(dir.eq.2) rigid(k,6)= (co(1,node)-cg(1))
          endif
        enddo
      enddo
!     rigid is stored column-major (dof varies fastest), matching
!     MatrixMarket dense array order, so the whole block can be
!     written in one statement (one value per output line)
!
      write(50,'(e23.15)') rigid(1:nfull,1:6)
      close(50)
!
      return
      end
