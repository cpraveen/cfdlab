      subroutine read_input
      use var
      implicit none
      
      open(1,file="input.dat")
      read(1,*)
      read(1,'(A40)') ifilename
      read(1,*)
      read(1,'(A40)') ofilename
      read(1,*)
      read(1,'(F4.2,X,I1,X,I1)') version, file_type, data_size
      read(1,*)
      read(1,'(I1)') ntags
      read(1,*)
      read(1,*) wst, wen
      close(1)

      end subroutine read_input
