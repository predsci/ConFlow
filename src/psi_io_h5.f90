!#######################################################################
! ****** PSI I/O H5: PSI HDF5 I/O Wrapper code.
!
!     Authors:  PSI
!
!     Predictive Science Inc.
!     www.predsci.com
!     San Diego, California, USA 92121
!#######################################################################
! Copyright 2022 Predictive Science Inc.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!#######################################################################
module rp1d_def
!
!-----------------------------------------------------------------------
! ****** Define a structure to hold a REAL 1D pointer.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      type :: rp1d
        real(r_typ), dimension(:), pointer, contiguous :: f
      end type
!
end module
!#1######################################################################
module ds_def
!
!-----------------------------------------------------------------------
! ****** Definition of the IO data structure.
!-----------------------------------------------------------------------
!
      use number_types
      use rp1d_def
!
      implicit none
!
      integer, parameter, private :: mxdim = 3
!
      type :: ds
        integer :: ndim
        integer, dimension(mxdim) :: dims
        logical :: scale
        logical :: hdf32
        type(rp1d), dimension(mxdim) :: scales
        real(r_typ), dimension(:,:,:), pointer, contiguous :: f
      end type
!
end module
!#######################################################################
subroutine rdh5 (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 1D, 2D, or 3D scientific data set from an HDF5 file.
! ****** The HDF5 file is currently assumed to contain only one
! ****** dataset (1D,2D,or 3D), with or without scales, in group "/",
! ****** and has no other data members.
!
!-----------------------------------------------------------------------
!
! ****** This routine allocates the required memory and returns
! ****** pointers to the data and scale arrays.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    File name to read from.
!
! ****** Output arguments:
!
!          S       : [structure of type DS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was read
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
! ****** Components of structure S:
!
!          NDIM    : [integer]
!                    Number of dimensions found in the data set.
!
!          DIMS    : [integer, dimension(3)]
!                    Number of points in the data set dimensions.
!                    For a 1D data set, DIMS(2)=DIMS(3)=1.
!                    For a 2D data set, DIMS(3)=1.
!
!          SCALE   : [logical]
!                    Flag to indicate the presence of scales (axes)
!                    in the data set.  SCALE=.false. means that scales
!                    were not found; SCALE=.true. means that scales
!                    were found.
!
!          HDF32   : [logical]
!                    Flag to indicate the precision of the data set
!                    read in.  HDF32=.true. means that the data is
!                    32-bit; HDF32=.false. means that the data is
!                    64-bit.
!
!          SCALES  : [structure of type RP1D, dimension(3)]
!                    This array holds the pointers to the scales
!                    when SCALE=.true., and is undefined otherwise.
!
!          F       : [real, pointer to a rank-3 array]
!                    This array holds the data set values.
!
! ****** The storage for the arrays pointed to by F, and the
! ****** scales (if present) in structure SCALES, is allocated by
! ****** this routine.
!
!-----------------------------------------------------------------------
!
      use iso_fortran_env
      use number_types
      use ds_def
      use hdf5
      use h5ds
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(ds) :: s
      character(*) :: fname
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,n_members,obj_type
!
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: dspace_id     ! Dataspace identifier
      integer(HID_T) :: dim_id        ! Dimension identifiers
      integer(HID_T) :: datatype_id   ! Datatype identifiers
!
      integer(SIZE_T) :: prec
!
      integer(HSIZE_T),dimension(:), allocatable :: s_dims,maxpts
      integer(HSIZE_T),dimension(1) :: s_dims_i
!
      real(REAL32), dimension(:,:,:), allocatable :: f4
      real(REAL32), dimension(:),     allocatable :: f4dim
      real(REAL64), dimension(:,:,:), allocatable :: f8
      real(REAL64), dimension(:),     allocatable :: f8dim
!
      character(512) :: obj_name
      character(4), parameter :: cname='RDH5'
!
      logical :: is_scale
!
!-----------------------------------------------------------------------
!
! ****** Initialize dimension count and arrays.
!
      s%ndim = 0
      s%dims(:) = 1
!
! ****** Initialize hdf5 interface.
!
      call h5open_f (ierr)
!
! ****** Open hdf5 file.
!
      call h5Fopen_f (trim(fname),H5F_ACC_RDONLY_F,file_id,ierr)
!
! ****** Get information about the hdf5 file.
!
      call h5Gn_members_f (file_id,"/",n_members,ierr)
!
! ****** Make sure there is (at maximum) one 3D dataset with scales.
!
      if (n_members.eq.0.or.n_members.gt.4) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Input file contains too few/many datasets.'
        write (*,*) 'File name: ',trim(fname)
        return
      endif
!
! ****** Assume the Dataset is in index 0 and get its name.
!
      call h5Gget_obj_info_idx_f (file_id,"/",0,obj_name,obj_type,ierr)
!
! ****** Open Dataset.
!
      call h5Dopen_f (file_id,trim(obj_name),dset_id,ierr)
!
! ****** Make sure the Dataset is not a scale.
!
      call h5DSis_scale_f(dset_id,is_scale,ierr)
      if (is_scale) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Input file Dataset at index 0 is a scale.'
        write (*,*) 'File name: ',trim(fname)
        return
      endif
!
! ****** Get dimensions (need s_dims array for format requirements).
!
      call h5Dget_space_f (dset_id,dspace_id,ierr)
      call h5Sget_simple_extent_ndims_f (dspace_id,s%ndim,ierr)
!
      allocate(s_dims(s%ndim))
!
      allocate(maxpts(s%ndim))
      call h5Sget_simple_extent_dims_f (dspace_id,s_dims,maxpts,ierr)
      deallocate(maxpts)
!
      do j=1,s%ndim
        s%dims(j) = INT(s_dims(j))
      end do
!
! ****** Get the floating-point precision of the data and set flag.
!
      call h5Dget_type_f (dset_id,datatype_id,ierr)
      call h5Tget_precision_f (datatype_id,prec,ierr)
!
      if (prec.eq.32) then
        s%hdf32=.true.
      elseif (prec.eq.64) then
        s%hdf32=.false.
      end if
!
! ****** Allocate the memory for the Dataset array in s.
!
      allocate (s%f(s%dims(1),s%dims(2),s%dims(3)))
!
! ****** Need to read the file in its own datatype, and then convert
! ****** to datatype of s%f.
!
      if (s%hdf32) then
        allocate (f4(s%dims(1),s%dims(2),s%dims(3)))
        call h5Dread_f (dset_id,datatype_id,f4,s_dims,ierr)
        do k=1,s%dims(3)
          do j=1,s%dims(2)
            do i=1,s%dims(1)
              s%f(i,j,k) = REAL(f4(i,j,k),r_typ)
            enddo
          enddo
        enddo
        deallocate (f4)
      else
        allocate (f8(s%dims(1),s%dims(2),s%dims(3)))
        call h5Dread_f (dset_id,datatype_id,f8,s_dims,ierr)
        do k=1,s%dims(3)
          do j=1,s%dims(2)
            do i=1,s%dims(1)
              s%f(i,j,k) = REAL(f8(i,j,k),r_typ)
            enddo
          enddo
        enddo
        deallocate (f8)
      end if
!
      deallocate(s_dims)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDH5:'
        write (*,*) '### Error while reading the dataset.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from H5DREAD_F) = ',ierr,']'
        ierr=4
        return
      end if
!
! ****** Close the hdf5 type descriptor.
!
      call h5Tclose_f (datatype_id,ierr)
!
! ****** Check if there might be scales present, if so, read them.
!
      if (n_members.gt.1) then
!
! ***** First check that the number of scale datasets match the # dim.
!
        if (n_members-1.ne.s%ndim) then
          write (*,*)
          write (*,*) '### ERROR in RDH5:'
          write (*,*) '### # scales does not match # dims.'
          write (*,*) 'File name: ',trim(fname)
          return
        end if
!
        s%scale=.true.
!
! ****** Loop through scales, make sure each is a scale, and read them.
!
        do i=1,n_members-1
!
! ****** Get the name of scale dataset.
!
          call h5Gget_obj_info_idx_f (file_id,"/",i, &
                                     obj_name,obj_type,ierr)
!
! ****** Open scale dataset.
!
          call h5Dopen_f (file_id,trim(obj_name),dim_id,ierr)
!
! ****** Make sure the scale is a scale.
!
          call h5DSis_scale_f (dim_id,is_scale,ierr)
          if (.not.is_scale) then
            write (*,*)
            write (*,*) '### ERROR in RDH5:'
            write (*,*) '### Scale is not a scale.'
            write (*,*) 'File name: ',trim(fname)
            return
          end if
!
! ****** Get dimension of scale.
!
          s_dims_i = s%dims(i)
!
! ****** Allocate scale.
!
          allocate (s%scales(i)%f(s_dims_i(1)))
!
! ****** Get the floating-point precision of the scale.
!
          call h5Dget_type_f (dim_id,datatype_id,ierr)
          call h5Tget_precision_f (datatype_id,prec,ierr)
!
! ****** Read in the scale data.
!
          if (prec.eq.32) then
            allocate (f4dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f4dim,s_dims_i,ierr)
            do j=1,s%dims(i)
              s%scales(i)%f(j) = REAL(f4dim(j),r_typ)
            end do
            deallocate (f4dim)
          elseif (prec.eq.64) then
            allocate (f8dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f8dim,s_dims_i,ierr)
            do j=1,s%dims(i)
              s%scales(i)%f(j) = REAL(f8dim(j),r_typ)
            end do
            deallocate (f8dim)
          end if
!
! ****** Close the scale dataset.
!
          call h5Dclose_f (dim_id,ierr)
!
        enddo
!
! ****** Allocate dummy scales (of length 1) for empty dimensions.
!
        do i=s%ndim+1,3
          allocate (s%scales(i)%f(1))
        enddo
      else
!
! ****** If scales are not present, allocate dummy
! ****** scales (of length 1) so that the pointers to the scales
! ****** are valid.
!
        s%scale = .false.
!
        allocate (s%scales(1)%f(1))
        allocate (s%scales(2)%f(1))
        allocate (s%scales(3)%f(1))
      end if
!
! ****** Close the dataset.
!
      call h5Dclose_f (dset_id,ierr)
!
! ****** Close the file.
!
      call h5Fclose_f (file_id,ierr)
!
! ****** Close FORTRAN interface.
!
      call h5close_f (ierr)
!
end subroutine
!#######################################################################
subroutine wrh5 (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 1D, 2D, or 3D scientific data set to an HDF5 file.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    File name to write to.
!
!          S       : [structure of type DS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
! ****** Output arguments:
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was written
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
!-----------------------------------------------------------------------
!
      use iso_fortran_env
      use number_types
      use ds_def
      use hdf5
      use h5ds
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      type(ds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      character(8) ::   dimname
      integer :: i,j,k
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: dspace_id,dspacedim_id   ! Dataspace identifiers
      integer(HID_T) :: dim_id        ! Dimension identifiers
      integer(HSIZE_T),dimension(3) :: s_dims
      integer(HSIZE_T),dimension(1) :: s_dims_i
!
      real(REAL32), dimension(:,:,:), allocatable :: f4
      real(REAL32), dimension(:),     allocatable :: f4dim
      real(REAL64), dimension(:,:,:), allocatable :: f8
      real(REAL64), dimension(:),     allocatable :: f8dim
!
!-----------------------------------------------------------------------
!
! ****** HDF5 calls are picky about the integer format for the dims
! ****** so the s%dims need to be converted to HSIZE_T integers.
!
! ****** Also, sometimes calls to wrhdf() for 1D and 2D datasets
! ****** do not have the unused dims(i) set to 1 (unset).
! ****** To avoid needing a function call to implicitly reshape
! ****** f(n), set the dims here.
!
      do i=1,3
         if (i.le.s%ndim) then
           s_dims(i) = INT(s%dims(i),HSIZE_T)
         else
           s_dims(i) = 1
         endif
      end do
!
! ****** Initialize hdf5 interface.
!
      call h5open_f (ierr)
!
! ****** Create the file.
!
      call h5Fcreate_f (trim(fname),H5F_ACC_TRUNC_F,file_id,ierr)
!
! ****** Create the dataspace.
!
      call h5Screate_simple_f (s%ndim,s_dims,dspace_id,ierr)
!
! ****** Create and write the dataset (convert s%f to proper type).
!
      if (s%hdf32) then
        allocate (f4(s_dims(1),s_dims(2),s_dims(3)))
        do k=1,s%dims(3)
          do j=1,s%dims(2)
            do i=1,s%dims(1)
              f4(i,j,k) = REAL(s%f(i,j,k),REAL32)
            enddo
          enddo
        enddo
        call h5Dcreate_f (file_id,'Data',H5T_NATIVE_REAL, &
                          dspace_id,dset_id,ierr)
        call h5Dwrite_f (dset_id,H5T_NATIVE_REAL,f4,s_dims,ierr)
        deallocate (f4)
      else
        allocate (f8(s_dims(1),s_dims(2),s_dims(3)))
        do k=1,s%dims(3)
          do j=1,s%dims(2)
            do i=1,s%dims(1)
              f8(i,j,k) = REAL(s%f(i,j,k),REAL64)
            enddo
          enddo
        enddo
        call h5Dcreate_f (file_id,'Data',H5T_NATIVE_DOUBLE, &
                          dspace_id,dset_id,ierr)
        call h5Dwrite_f (dset_id,H5T_NATIVE_DOUBLE,f8,s_dims,ierr)
        deallocate (f8)
      endif
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRH5:'
        write (*,*) '### Could not write the dataset.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from h5Dwrite_f) = ',ierr,']'
        ierr=4
        return
      end if
!
! ****** Check for scales.  If present, add them to the hdf5 dataset.
!
      if (s%scale) then
        do i=1,s%ndim
          if (i.eq.1) then
            dimname='dim1'
          elseif (i.eq.2) then
            dimname='dim2'
          elseif (i.eq.3) then
            dimname='dim3'
          endif
          s_dims_i = s_dims(i)
          call h5Screate_simple_f(1,s_dims_i,dspacedim_id,ierr)
          if (s%hdf32) then
            allocate (f4dim(s_dims_i(1)))
            do j=1,s%dims(i)
              f4dim(j) = REAL(s%scales(i)%f(j),REAL32)
            end do
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_REAL, &
                              dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_REAL, &
                             f4dim,s_dims_i,ierr)
            deallocate (f4dim)
          else
            allocate (f8dim(s_dims_i(1)))
            do j=1,s%dims(i)
              f8dim(j) = REAL(s%scales(i)%f(j),REAL64)
            end do
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_DOUBLE, &
                             dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_DOUBLE, &
                             f8dim,s_dims_i,ierr)
            deallocate (f8dim)
          endif
          call h5DSset_scale_f (dim_id,ierr,dimname)
          call h5DSattach_scale_f (dset_id,dim_id,i,ierr)
          call h5DSset_label_f(dset_id, i, dimname, ierr)
          call h5Dclose_f (dim_id,ierr)
          call h5Sclose_f (dspacedim_id,ierr)
        end do
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in WRH5:'
          write (*,*) '### Could not write the scales.'
          write (*,*) 'File name: ',trim(fname)
          ierr = 5
          return
        endif
      endif
!
! ****** Close the dataset.
!
      call h5Dclose_f (dset_id,ierr)
!
! ****** Close the dataspace.
!
      call h5Sclose_f (dspace_id,ierr)
!
! ****** Close the file.
!
      call h5Fclose_f (file_id,ierr)
!
! ****** Close the hdf5 interface.
!
      call h5close_f (ierr)
!
end subroutine
!#######################################################################
