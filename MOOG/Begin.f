
      subroutine begin
c***************************************************************************
c     This routine simply starts up MOOG
c     THIS VERSION IS FOR LINUX REDHAT MACHINES
c***************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Pstuff.com'
      character*80 line, systemcall
      integer num


c*****define the number of text screen lines for silent mode;
c     this number is hardwired, since it is not really needed at run time.
c      if (silent .eq. 'y') then
c         maxline = 24
c         go to 10
c      endif


c*****define the number of lines available on the text screen for
c     interactive mode; this number is discovered from the "stty -a"
c     command, for which the output format is unique to the operating
c     system.
c      write (systemcall,*) 'stty -a > tmpsize'
c      call system (systemcall)
c      open (99,file='tmpsize')
c5     read (99,1010,end=15) line
c      do i=1,77
c         if (line(i:i+3) .eq. 'rows') then
c            if     (machine .eq. 'pcl') then
c               read (line(i+4:i+6),1011) maxline
c            elseif (machine .eq. 'mac') then
c               read (line(i-4:i-2),1011) maxline
c            elseif (machine .eq. 'uni') then
c               read (line(i+6:i+8),1011) maxline
c            endif
c            go to 10
c         endif
c      enddo
c      go to 5
c15    array = 'SCREEN ROW COUNT UNKNOWN; USE 24 ([y]/n)? '
c      nchars = 42
c      ikount = 2
c      call getasci (nchars,ikount)
c     choice = chinfo(1:1)
c      if (choice.eq.'y' .or. nchars.le.0) then
c         go to 10
c      else
c         call finish (0)
c      endif
c10    close (99,status='delete')
c      write (systemcall,*) '\\rm -f tmpsize'
c      call system (systemcall)
c      if (maxline .lt. 10) then
c         maxline = 24
c      else
c         maxline = maxline - 2
c      endif


c*****clear the text screen
c      write (systemcall,*) 'clear'
c      call system (systemcall)


c*****open data files carried with the source code: Barklem damping
      nfbarklem = 35
      num = 60
      call getcount (num,moogpath)
      if (moogpath(num:num) .ne. '/') then
         num = num + 1
         moogpath(num:num) = '/'
      endif
      fbarklem(1:num) = moogpath(1:num)
      fbarklem(num+1:num+11) = 'Barklem.dat'
      open (nfbarklem,file=fbarklem)

 
c*****open data files carried with the source code: Barklem UV damping
      nfbarklemUV = 36
      num = 60
      call getcount (num,moogpath)
      if (moogpath(num:num) .ne. '/') then
         num = num + 1
         moogpath(num:num) = '/'
      endif
      fbarklemUV(1:num) = moogpath(1:num)
      fbarklemUV(num+1:num+13) = 'BarklemUV.dat'
      open (nfbarklemUV,file=fbarklemUV)
 

c  write a header and find the appropriate parameter file, and exit normally
      write (array,1001)
      istat = ivwrite (1,1,array,79)
      write (array,1004)
      istat = ivwrite (2,1,array,79)
      array = 'MOOG PARAMETERS? ' 
      nchars = 15
      nfparam = 50     
      lscreen = 4
      if (silent .eq. 'y') then
         fparam = 'batch.par'
      else
         fparam = 'no_filename_given'     
      endif
      call infile ('input  ',nfparam,'formatted  ',0,nchars,
     .             fparam,lscreen)
      read (nfparam,1002) control
      write (array,1003) control
      istat = ivwrite (2,1,array,58)
      write (array,1001)
      istat = ivwrite (3,1,array,79)
      return


c*****format statements
1001  format (79('*'))
1002  format (a7)
1003  format (22x,'MOOG IS CONTROLLED BY DRIVER ',a7)
1004  format (25(' '),'MOOG LTE VERSION (JUN 2014)',26(' '))   
1010  format (a80)
1011  format (i3)


      end
      







