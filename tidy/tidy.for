      PROGRAM TIDY
C
C===================================================================
C                                                                  *
C                    * * *   T I D Y   * * *                       *
C                    Version 7.2, 1999-10-15                       *
C      A FORTRAN PROGRAM TO RENUMBER AND OTHERWISE CLEAN UP        *
C             OLD AND TIRED FORTRAN SOURCE PROGRAMS.               *
C                                                                  *
C                   IN ADDITION TO RENUMBERING,                    *
C             TIDY PROVIDES A LIMITED SET OF FORTRAN               *
C                          DIAGNOSTICS.                            *
C                                                                  *
C                 ANSI FORTRAN  (ANSI X3.9-1978)                   *
C                                                                  *
C===================================================================
C                                                                  *
C  Copyright (C) 1989, The Regents of the University of California *
C                      All Rights Reserved                         *
C                                                                  *
C  THE REGENTS OF THE UNIVERSITY OF CALIFORNIA MAKE NO REPRESENTA- *
C  TION OR WARRANTIES WITH RESPECT TO THE CONTENTS HEREOF AND      *
C  SPECIFICALLY DISCLAIM ANY IMPLIED WARRANTIES OF MERCHANTABILITY *
C  OR FITNESS FOR ANY PARTICULAR PURPOSE.                          *
C                                                                  *
C  Further, the Regents of the University of California reserve the*
C  right to revise this software and/or documentation and to make  *
C  changes from time to time in the content hereof without obliga- *
C  tion of the Regents of the University of California to notify   *
C  any person of such revision or change.                          *
C                                                                  *
C  PERMISSION TO COPY AND DISTRIBUTE THIS PROGRAM, AND TO MAKE     *
C  DERIVATIVE WORKS HEREFROM, IS GRANTED PROVIDED THAT THIS COPY-  *
C  RIGHT NOTICE IS RETAINED IN ALL SOURCE CODE AND USER MANUALS.   *
C                                                                  *
C==================================================================*
C                                                                  *
C                 **************************                       *
C                *         PROGRAM          *                      *
C               *     AND SUBROUTINES BY     *                     *
C              *        HARRY M MURPHY        *                    *
C             *  AIR FORCE WEAPONS LABORATORY  *                   *
C              *   KIRTLAND AIR FORCE BASE    *                    *
C               *         NEW MEXICO         *                     *
C                *         1 9 6 6          *                      *
C                 **************************                       *
C                                                                  *
C     TIDY ACCEPTS ASA FORTRAN WITH 19 CONTINUATION CARDS          *
C     AS WELL AS SOME IBM AND CDC DIALECT FORTRAN STATEMENTS       *
C                                                                  *
C     THIS VERSION MODIFIED FOR USE AT LRL BERKELEY BY             *
C     GERRY TOOL (1967). (STILL CDC/6600)                          *
C                                                                  *
C     THIS PROGRAM HAS BEEN REVISED FOR IBM 360/67 BY ALICE        *
C     V BARLOW, NASA AMES, SUMMER 1972                             *
C                                                                  *
C     ADDITIONS AND REWORKING BY ROGER CHAFFEE, LRL BERKELEY       *
C     AND SLAC COMPUTATIONS RESEARCH GROUP, 1968-1982              *
C                                                                  *
C     CONVERTED TO IBM (RYAN-McFARLAND) PROFESSIONAL FORTRAN       *
C     BY AL STANGENBERGER, DEPT. OF FORESTRY, U.C. BERKELEY        *
C                                                                  *
C     Additions and conversion to command line version             *
C     by Ajit J. Thakkar, Dept. of Chemistry, U. New Brunswick     *
C     E-mail: ajit@unb.ca  WWW  http://www.unb.ca/chem/ajit/       *
C                                                                  *
C==================================================================*
C
C
C  INPUT/OUTPUT
C     FUNCTION          FORTRAN UNIT   CURRENT VALUE
C      CONSOLE OUTPUT     STDERR            6
C      CONSOLE INPUT      STDIN             5  (unused)
C      CONTROL CARD       USRFIL            3
C      INPUT              INFILE            4
C      LIST OUTPUT        OUTFIL            7
C      CARD OUTPUT        PUNFIL            8
C      SCRATCH(NORMAL)    SCFIL1            1
C      SCRATCH(FORMATS)   SCFIL2            2
C
C================================================================
C     I N S T A L L A T I O N   N O T E S
C
C     1.  INCLUDE statements are used to incorporate common block
C         definitions into most subroutines.  Check syntax as these
C         statements are system-dependent.
C
C     2.  CHARACTER SET SPECIFICITY -
C         The code for horizontal tab differs in EBCDIC and ASCII.
C         This value is set (KTAB) in this routine. Fix as needed.
C
C     3.  Non-standard intrinsics.
C         Intrinsic subroutine EXIT(n) is used to abort the program 
C         with a non-zero error code passed to the operating system.
C         Many compilers including Lahey F77L, LF90, LF95 and GNU g77
C         provide an EXIT intrinsic.
C
C         The first and only command line argument is obtained by a
C         call to GetArg. A GetArg subroutine for Lahey F77L, LF90 
C         and LF95 is included in this distribution, and GNU g77 has
C         an intrinsic GetArg.
C
C         Aside from these factors, the rest of the program is
C         standard Fortran-77.
C
C================================================================
C     Programming notes:
C
C     In subroutine holscn, Hollerith characters are changed
C     so they won't be recognized by any other test by
C     changing second character to '@'
C
C     Subroutines holscn and contrl invoke function kupper to convert
C     lower-case alphabetic characters to upper case (except for
C     Hollerith strings).
C
C     The character $ is treated as an alpha in IBM Fortran.
C     the data statement for the special characters, kspk, has
C     been changed so that $ is not recognized as a special
C     character.  this data statement should be changed back
C     on non-IBM systems.
C
C     Subroutine redstr is set up to accommodate an apparent bug
C     in the Ryan-McFarland professional Fortran compiler, that
C     unformatted sequential records seem to be limited to 1024 bytes.
C     Since each record has a 4-byte header and trailer, writes 508
C     character*2 elements, or 254 integer*4 per record.  this may
C     vary for other compilers.
C
C
C  INTERNAL FLAGS (JUST A LIST.  WHERE ELSE TO PUT IT...)
C     MANSI =  0 FLAG ALL NON-ANSI (FORTRAN-77) STATEMENTS
C           =  1 DO NOT FLAG NON-ANSI STATEMENTS
C     MP2   =  1 DO PASS2
C           =  0 NO PASS 2
C     MCOL  = -1 COLLECT FORMAT STATEMENTS AT END
C           =  0 LEAVE THEM IN PLACE
C     MILDO = -1 IF DO-TERMINATOR ALLOWED BUT NON-STANDARD
C           =  0 IF DO-TERMINATOR ALLOWED
C           = +1 IF DO-TERMINATOR FORBIDDEN
C     MCONT =  0 REMOVE CONTINUE CARDS AND DOUBLE BRANCHES
C           =  1 LEAVE THEM
C     MTRAN = -1 CURRENT CARD IS AN UNCONDITIONAL BRANCH
C           =  0 CURRENT CARD NOT NECESSARILY A BRANCH
C     NTRAN =    SAME AS MTRAN, BUT REFERS TO PREVIOUS CARD
C     MLGC  = -1 NORMAL STATEMENT
C           =  0 STATEMENT IS CONTROLLED BY A LOGICAL IF
C     MRIT  =  N LEFT ADJUST TO COLUMN N
C           = -N RIGHT ADJUST TO COLUMN N
C     MDEB  =  0 *NODEBUG
C           =  1 *DEBUG
C     KD15  =    STATEMENT INCREMENT (*STAT=...)
C     KB15  =    STATEMENT BASE (*BASE=...)
C     MPUN  =  0 NO PUNCH OUTPUT
C           =  1 MAKE PUNCH OUTPUT
C     KPUN       SAVES *CARD/*NOCARD (1/0) FOR MPUN VALUE
C     MLIST = -1 (*LIST) LIST PASS 1
C           =  0 (*NOLIST) DONT
C     KPRIN =  1 (*LIST=2) LIST PASS 2
C           =  0 (*NOLIST=2) DONT
C     MPRIN =    KPRIN AT START OF ROUTINE. MAY CHANGE IF ERROR
C                  AT START OF PASS1.
C     KOUNT      COUNTS CARDS IN FOR CURRENT ROUTINE.
C     IQUIT =  0 UNTIL INPUT ENDFILE IS FOUND IN READER.
C           =  1 THEREAFTER
C     MSTOP =  0 NORMALLY
C           = -1 FOR *STOP CARD FOUND--TIME TO FINISH UP
C           =  1 FOR STOP NOW.
C
C================================================================
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
*     LOGICAL DOUSER,SCDISK
      LOGICAL DOUSER
      COMMON /TDYVER/ VERNUM
      CHARACTER*30 VERNUM
      CHARACTER*64 CFILNM
C
      DOUSER=.TRUE.
C
C     SCDISK .TRUE. ALLOWS USER TO SPECIFY DISK TO HOLD SCRATCH FILES.
C          FOR UNIX SYSTEMS, SHOULD SET TO .FALSE.
C     This variable has no meaning in version 7.2
*     SCDISK=.TRUE.
C
C     VALUE FOR TAB AS ASCII
      KTAB=KBL
      KTAB(1:1)=CHAR(9)
C     VALUE FOR TAB AS EBCDIC
C     KTAB(1:1)=CHAR(5)
C
C  Following call is for version 7.2 (two arguments only)
      CALL PCTIDY (DOUSER,CFILNM)
C     CALL PCTIDY (DOUSER,SCDISK,CFILNM)
C
C     INITIALIZE PROGRAM
      CALL INITDY
C     ADJUST ROUTINE NUMBER - PASS1 WILL INCREMENT IT.
      NROUT=NROUT-1
C
C     PROCESS USER CONTROL CARD FILE.
      IF (DOUSER) CALL USRCON (CFILNM)
C
      WRITE (STDERR,30)
      CALL READER
C
 10   CALL PASS1
      IF (MSTOP.NE.0) THEN
           IF (MSTOP.GT.0) GO TO 20
           IF (KOUNT.LE.0) GO TO 20
      END IF
      CALL EDIT
      IF (MP2.EQ.0) GO TO 10
      IF (MREF.NE.0) CALL RDIR
      IF (MDEB.NE.0) WRITE (STDERR,23)
      CALL PASS2
      IF (IQUIT.NE.0) GO TO 20
      IF (MSTOP.EQ.0) GO TO 10
C                            ALL DONE
20    IF (MDEB.NE.0) WRITE (STDERR,25)
      REWIND SCFIL1
      REWIND SCFIL2
      IF (NMSG.GT.0) THEN
           WRITE (OUTFIL,40) NMSG
           WRITE (STDERR,40) NMSG
      ELSE
           WRITE (OUTFIL,50)
           WRITE (STDERR,50)
      END IF
      WRITE (OUTFIL,60) NPUN,VERNUM
C
C     ABNORMAL TERMINATIONS HANDLED BY SUBROUTINE DIAGNO.
      IF (LERR.GT.0) CALL DIAGNO (47)
C
C     GET RID OF SCRATCH FILES UNLESS DEBUGGING
      IF (MDEB.EQ.0) THEN
           CLOSE (SCFIL1,STATUS='DELETE')
           CLOSE (SCFIL2,STATUS='DELETE')
      END IF
C
      STOP
C
30    FORMAT (' Running')
23    format (' Begin Pass2')
25    FORMAT (' Rewinding scratch files - main program')
40    FORMAT (' Warning .',I5,' Diagnostic messages have been generated
     &in this tidy run.')
50    FORMAT (' No diagnostic messages were generated during this tidy r
     1un.')
60    FORMAT (' ',I5,' cards were punched.'/' ',A)
      END
      SUBROUTINE PCTIDY (DOUSER,CFILNM)
C  Version 7.2
C    by Ajit J. Thakkar, Chemistry Dept., U. of New Brunswick
C  Standard F77 file definitions
      INCLUDE 'TIDY.INC'
      INCLUDE 'UNITS.INC'
      COMMON/TDYVER/VERNUM
      CHARACTER*30 VERNUM
      CHARACTER*64 FILNM1, FILNM2, CFILNM
      CHARACTER*64 CLNAME
      INTEGER IDot
      LOGICAL DOUSER, FORFIL

C  Display version number on console
      WRITE (STDERR,'(1x,A)') VERNUM

C  Find and open user control file.
      INQUIRE (FILE='tidy.ini',EXIST=douser)
      IF (DOUSER) THEN
        CFILNM='tidy.ini'
        OPEN(USRFIL,FILE=CFILNM,STATUS='OLD')
        REWIND(USRFIL)
      END IF

C  Find and open source file
C  Get first command line argument using GetArg 
C  A GetArg subroutine (using GetCl) is provided for Lahey compilers
C  whereas g77 has GetArg as a non-standard intrinsic
      CALL GetArG(1,CLNAME)
      FILNM1=CLNAME
      Idot=INDEX(FILNM1,'.')
      IF (Idot.eq.0) THEN
        IDOT=INDEX(FILNM1,' ')
        FILNM1=FILNM1(1:Idot-1)//'.for'
        INQUIRE (FILE=FILNM1,EXIST=forfil)
        IF (.NOT.FORFIL) THEN
          FILNM1=FILNM1(1:Idot-1)//'.f  '
          INQUIRE (FILE=FILNM1,EXIST=forfil)
          IF (.NOT.FORFIL) THEN
            WRITE (STDERR,'(A)') ' Source file does not exist'
            STOP
          END IF
        END IF
      END IF
      OPEN(INFILE,FILE=FILNM1,STATUS='OLD')
      REWIND(INFILE)

C  Open listing and "punched output" files.
C  formatted, sequential, overwrite old files if they exist
      FILNM1=FILNM1(1:Idot)//'lis'
      OPEN(OUTFIL,FILE=FILNM1)
      REWIND(OUTFIL)

      FILNM1=FILNM1(1:Idot)//'tid'
      OPEN(PUNFIL,FILE=FILNM1)
      REWIND(PUNFIL)

C  Open scratch files
C  unformatted, sequential, overwrite old files if they exist
      FILNM1='SCFIL1.TDY'
      OPEN(SCFIL1,FILE=FILNM1,FORM='UNFORMATTED')
      REWIND(SCFIL1)

      FILNM2='SCFIL2.TDY'
      OPEN(SCFIL2,FILE=FILNM2,FORM='UNFORMATTED')
      REWIND(SCFIL2)

      RETURN
      END
      SUBROUTINE QUIT (M)
      INCLUDE 'TIDY.INC'
      INCLUDE 'UNITS.INC'
      IF (MDEB.EQ.0) THEN
        CLOSE(SCFIL1,STATUS='DELETE')
        CLOSE(SCFIL2,STATUS='DELETE')
      END IF
C   All abnormal terminations are handled with Exit which sets the 
C   error condition to its argument.
C   Many compilers including F77L, LF90, LF95 and g77 provide Exit as 
C   a non-standard intrinsic.
      CALL Exit(M)
      RETURN 
      END 
      BLOCK DATA MISDAT
C
C     THIS BLOCK DATA CONTAINS MISCELLANEOUS DATA STATEMENTS FOR TIDY.
C
C     VERSION 6.2 MODIFICATION -----------------------------------------
C     VARIABLES WHICH ARE CONTROLLED BY SUBROUTINE CONTRL ARE SET IN
C     SUBROUTINE INITDY.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
      COMMON /TDYVER/ VERNUM
      CHARACTER*30 VERNUM
C
C     /ALPHA/
      DATA KBL,KDIG/' ','0','1','2','3','4','5','6','7','8','9'/
      DATA KABC/'A','B','C','D','E','F','G','H','I','J','K','L','M','N',
     1'O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA KSPK/'=',',','(','/',')','+','-','*','.','X$','-','''','&','$
     1','!'/
*    1'/
C  $ IN ABOVE STATEMENT REPLACED BY X$, SINCE $ IS NOT SPECIAL
C  CHARACTER IN IBM 360/370 FORTRAN.
      DATA KBL2,KLR2,KLP2,KRP2,KERM/' *','$$','((','))',' $'/
      DATA KAMPR/'& '/,KAT/' @'/,KAPSTR/'''@'/
C
C     /MISCAL/
      DATA KEND/'D','N','E'/
C
C
C     /MISC/
C     LOGICAL UNIT ASSIGNMENTS
      DATA INFILE/4/
      DATA OUTFIL/7/
      DATA PUNFIL/8/
      DATA STDERR/6/
      DATA STDIN/5/
      DATA SCFIL1/1/
      DATA SCFIL2/2/
      DATA USRFIL/3/
C
      DATA IQUIT/0/
      DATA KOUNT/0/
      DATA LERR/0/
      DATA LINE/1/
      DATA MDEB/0/
      DATA MSTOP/0/
      DATA MXREF/256/
      DATA MXRGHT/65/
      DATA NMSG/0/
      DATA NPAGE/0/
      DATA NPUN/0/
C
C     VERSION STRING
      DATA VERNUM/'Tidy 7.2  -  1999-10-15'/
      END
      SUBROUTINE CONTRL
      PARAMETER  (NKTRL=43)
C
C     THIS SUBROUTINE EXECUTES THE TIDY CONTROL STATEMENTS.
C     ALL TIDY CONTROL STATEMENTS MUST HAVE AN * PUNCHED IN COLUMN 1.
C
C     1   BASE   NOBASE   KB15
C     2   IDIN   ======   KD79
C     3   IDST   ======   KD79
C     4   ROUT   ======   NROUT
C     5   STAT   ======   KD15
C     6   CARD   NOCARD   MPUN
C     7   COLL   NOCOLL   MCOL
C     8   COMM   NOCOMM   MCOM
C     9   EXEM   NOEXEM   MEX
C     10  LABE   NOLABE   MLBL
C     11  LAST   ======   MSTOP
C     12  LIST   NOLIST   MLIST
C     13  NEWR   ======   NROUT
C     14  REFE   NOREFE   MREF
C     15  SKIP   ======   MSKP
C     16  STOP   ======   MSTOP
C     17  SERI   NOSERI   MSER  <0 USE KOL73...=0 USE BLANKS >0 SERIAL
C     18  RIGH   ======   MRIT
C     19  LEFT   ======   MRIT
C     20  COLU   NOCOLU   JUST
C     21  INDE   NOINDE   INDENT
C     22  DEBU   NODEBU   MDEB
C     23  CONT   NOCONT   MCONT
C     24  END    ======   SAME AS STOP
C     25  ANSI   NOANSI   MANSI
C     26  FEND   NOFEND   NFEND
C     27  CCHR   ======   KCTCTL
C     28  HTRA   ======   KHTRAN
C     29  DTRA   NODTRA   KDTRAN
C     30  DEL1   ======   KDEL1
C     31  DEL2   ======   KDEL2
C     32  ARET   ======   KALMRK
C     33  ARTR   NOARTR   KALTRN
C     34  BLAN   NOBLAN   KBKCOK (INCLUDE BLANK LINES IN DECK)
C     35  FSPL   NOFSPL   KFSPL  (SPLIT STRINGS IN INDENTED FMTS)
C     36  HLOG   NOHLOG   KHLOG  (LOG TRANSLATED H-FIELDS TO LISTING)
C     37  CASE   NOCASE   MCASE  (TRANSLATE NON-STRINGS TO UPPER CASE)
C     38  UCAS   ======   MCASE  (TRANSLATE NON-STRINGS TO UPPER CASE)
C     39  LCAS   ======   MCASE  (TRANSLATE NON-STRINGS TO LOWER CASE)
C     40  ENDO   NOENDO   MNDOO  (RETAIN END-DO STATEMENTS)
C     41  FMTB   ======   MFMTB  (BASE NUMBER FOR COLLECTED FORMATS)
C     42  MAXC   ======   MAXCNT (MAX. NUMBER OF CONTINUATIONS)
C     43  OLDS   NOOLDS   K72    (OLD DECK HAS SERIALIZATION)
C
      INCLUDE  'tidy.inc'
      INCLUDE  'units.inc'
C
      COMMON /CONTDY/ KTRL(4,NKTRL)
      CHARACTER*2 KTRL
      CHARACTER*2 KUPPER,IT
C
      IF (MDEB.GT.0) WRITE (OUTFIL,'('' CONTRL: '',75A1)') (JINT(I),I=1,
     1JMAX)
      I=14
      ISTAR=-1
      JSW=0
      JL=JMAX-1
C
C     SCAN FOR 'NO' AT START
      DO 10 JB=2,JL
           IT=JINT(JB)
           IF (IT.NE.KBL) THEN
                IT=KUPPER(IT)
                IF (IT.NE.KABC(I)) THEN
                     JC=2
                     GO TO 20
                END IF
C               MATCHED N OR O. INCREMENT CHARACTER.
                I=I+1
C               IF WE'RE TO P, THEN FOUND 'NO'
                IF (I.GT.15) THEN
                     JSW=1
                     JC=JB+1
                     GO TO 20
                END IF
           END IF
 10   CONTINUE
      ISTAR=1
      RETURN
C
 20   DO 40 J=1,NKTRL
           I=1
           DO 30 JCOL=JC,JMAX
                IT=KUPPER(JINT(JCOL))
                IF (IT.EQ.KTRL(I,J)) THEN
                     IF (I.GE.4) GO TO 60
                     I=I+1
                ELSE
                     IF (IT.NE.KBL) GO TO 40
                END IF
 30        CONTINUE
 40   CONTINUE
 50   ISTAR=1
      RETURN
C
C     EXECUTE CONTROL STATEMENT
C
 60   NREC=NREC-1
C                  JSW=1 IF CARD STARTS WITH NO
      IF (JSW.EQ.1) THEN
           GO TO (500,50,50,50,50,110,130,190,290,390,50,530,50,460,50,
     1      50,490,50,50,510,520,230,210,50,90,310,50,370,260,50,250,50,
     2      70,150,340,360,160,170,180,280,320,410,440),J
      ELSE
           GO TO (530,530,530,530,530,100,120,530,530,380,400,530,420,
     1      450,470,400,480,530,530,530,530,220,200,400,80,300,530,530,
     2      240,530,530,530,530,140,330,350,530,530,530,270,530,530,430)
     3      ,J
      END IF
C
C     IF FALL THRU HERE, ABORT - ABOVE LIST NOT COMPLETE.
      CALL DIAGNO (48)
C
C                  NOARTRAN
 70   KALTRN=KBL
      RETURN
C                  ANSI
 80   MANSI=0
      RETURN
C                  NOANSI
 90   MANSI=1
      RETURN
C                  CARD
 100  MPUN=-1
      KPUN=-1
      RETURN
C                  NOCARD
 110  MPUN=0
      KPUN=0
      RETURN
C                  COLL
 120  MCOL=-1
      RETURN
C                  NOCOLL
 130  MCOL=0
      RETURN
C                  BLAN
 140  KBKCOK=1
      RETURN
C                  NOBLAN
 150  KBKCOK=0
      RETURN
C                  NOCASE
 160  MCASE=-1
      RETURN
C
C                  NOUCASE - CHANGE TO *LCASE
 170  J=39
      IDL772=0
      GO TO 570
C
C                  NOLCASE - CHANGE TO *UCASE
 180  J=38
      IDL772=0
      GO TO 570
C                  NOCOMM
 190  MCOM=0
      RETURN
C                  CONT
 200  MCONT=1
      RETURN
C                  NOCONT
 210  MCONT=0
      RETURN
C                  DEBUG
 220  MDEB=1
      RETURN
C                  NODEBUG
 230  MDEB=0
      RETURN
C                  DTRAN
 240  KDTRAN=1
      RETURN
C                  NODEL2 -- IMPLIES *NODTRAN
 250  KDEL2='""'
C                  NODTRAN
 260  KDTRAN=0
      RETURN
C                  ENDO
 270  MNDOO=1
      RETURN
C                  NOENDO
 280  MNDOO=0
      RETURN
C                  NOEXEM
 290  MEX=0
      RETURN
C                  FEND
 300  NFEND=0
      RETURN
C                  NOFEND
 310  NFEND=1
      RETURN
C                  NOFMTB
 320  MFMTB=0
      RETURN
C                  FSPL
 330  KFSPL=0
      RETURN
C                  NOFSPL
 340  KFSPL=1
      RETURN
C                  HLOG
 350  KHLOG=0
      RETURN
C                  NOHLOG
 360  KHLOG=1
      RETURN
C                  NOHTRAN
 370  KHTRAN=0
      RETURN
C                  LABE
 380  MLBL=-1
      RETURN
C                  NOLABE
 390  MLBL=0
      RETURN
C                  LAST/STOP
 400  MSTOP=-1
      RETURN
C                  NOMAXC - DEFAULTS TO NQCNTS
 410  MAXCNT=NQCNTS
      RETURN
C                  NEWR
 420  CALL INITDY
      RETURN
C                  OLDS - OLD DECK HAS SERIAL
 430  K72=72
      RETURN
C                  NOOLDS - OLD DECK HAS SERIAL
 440  K72=80
      RETURN
C                  REFE
 450  MREF=-1
      RETURN
C                  NOREFE
 460  MREF=0
      RETURN
C                  SKIP
 470  MSKP=-1
      RETURN
C                  SERI
 480  MSER=-1
      RETURN
C                  NOSERI
 490  MSER=0
      RETURN
C                  NOBASE
 500  KB15=0
      RETURN
C
C                  NOCOLU
 510  JUST=0
      RETURN
C
C                  NOINDENT
 520  INDENT=0
      RETURN
C
C     GET NUMBER FOLLOWING (=) SIGN.
C
 530  JAVB=JCOL
      DO 540 JCOL=JAVB,JMAX
           IF (JINT(JCOL).EQ.KSPK(1)) GO TO 550
 540  CONTINUE
      L772=1D0
      GO TO 560
 550  JCOL=JCOL+1
      JAVB=JCOL
      CALL RSTAT
C
C     EVERYTHING REQUIRES AN INTEGER - THE DOUBLE PRECISION IS TO
C      HANDLE SOME BIG CONSTANTS IN DATA STATEMENTS WITH RSTAT.
 560  IDL772=int(L772)
      IF (MDEB.GT.0) THEN
           WRITE (OUTFIL,800) IDL772,J
      END IF
 570  GO TO (580,630,630,650,700,50,50,610,620,50,50,770,50,50,50,50,50,
     1750,760,710,740,50,50,50,50,50,780,730,50,780,780,780,780,50,50,
     250,590,590,600,50,720,640,50),J
C
C     IF FALL THRU HERE, ABORT - ABOVE LIST NOT COMPLETE.
      CALL DIAGNO (48)
C
C                  BASE
 580  KB15=IDL772
      RETURN
C
C                  CASE, UCAS
C     0 = ALL, 1 = KEYWORDS ONLY, 2 = NON-KEYWORDS ONLY
 590  MCASE=0
      IF (MOD(IDL772,2).EQ.0) CALL KCTSET (0)
      IF (IDL772.LE.1) CALL KWDWRT (-1)
      RETURN
C                  LCASE
 600  MCASE=0
      IF (MOD(IDL772,2).EQ.0) CALL KCTSET (1)
      IF (IDL772.LE.1) CALL KWDWRT (-2)
      RETURN
C                  COMM
C     KEEP *COMM = *COMM=1   FOR UPWARD COMPATIBILITY
C      *COMM=2 IS ALSO EQUIVALENT...
 610  MCOM=IDL772
      IF (MCOM.LE.2) THEN
           MCOM=-1
      ELSE IF (MCOM.GT.7) THEN
           MCOM=7
      END IF
      RETURN
C                  EXEM
 620  MEX=IDL772
C     KEEP *EXEM = *EXEM=1   FOR UPWARD COMPATIBILITY
      IF (MEX.LE.1) THEN
           MEX=-1
      ELSE IF (MEX.EQ.2) THEN
           MEX=1
      ELSE
           GO TO 50
      END IF
      RETURN
C                  IDIN/IDST
 630  KD79=max(IDL772,1)
      RETURN
C
C                  MAXCNT
 640  MAXCNT=max(IDL772,1)
      IF (MAXCNT.GT.NQCNTS) CALL DIAGNO (50)
      RETURN
C                  ROUT
C     USE TWO LETTERS FOR ROUTINE CODE, CONSTRUCT VALUE OF NROUT.
 650  JCOL=JAVB-1
      NROUT=0
      DO 680 I=1,2
 660       JCOL=JCOL+1
           IT=KUPPER(JINT(JCOL))
           IF (IT.EQ.KBL) GO TO 660
           IF (IT.EQ.KERM) GO TO 690
           DO 670 J=1,26
                IF (IT.NE.KABC(J)) GO TO 670
                NROUT=NROUT*26+J
                GO TO 680
 670       CONTINUE
 680  CONTINUE
C
 690  NROUT=max(NROUT-1,1)
      RETURN
C                  STAT
 700  KD15=max(IDL772,1)
      RETURN
C                  COLU
 710  JUST=max(IDL772,7)
      RETURN
C                  FMTB
 720  MFMTB=max(IDL772,0)
      RETURN
C                  HTRAN
 730  KHTRAN=min(IDL772,3)
      IF (KHTRAN.LT.0) KHTRAN=0
      RETURN
C                            INDENT
 740  INDENT=min(10,IDL772)
      RETURN
C                            RIGHT
 750  MRIT=min(IDL772,5)
      IF (MRIT.EQ.1) MRIT=5
      RETURN
C                            LEFT
 760  MRIT=max(IDL772,1)
      IF (MRIT.GT.5) MRIT=1
      MRIT=-MRIT
      RETURN
C                            LIST/NOLIST
 770  IF (IDL772.EQ.2) THEN
           IF (JSW.EQ.0) THEN
C                            LIST=2.
                KPRIN=1
                MPRIN=1
           ELSE
C                            NOLIST=2.
                MPRIN=0
                KPRIN=0
           END IF
      ELSE
           IF (JSW.EQ.0) THEN
C                            LIST
                MLIST=-1
           ELSE
C                            NOLIST
                MLIST=0
           END IF
      END IF
      RETURN
C
C                  CARDS USING CHARACTER ARGUMENT
 780  JCOL=JAVB-1
 790  JCOL=JCOL+1
      IT=KUPPER(JINT(JCOL))
      IF (IT.EQ.KBL) GO TO 790
      IF (J.EQ.27) THEN
C                            CCHR (CONTINUATION CHAR)
           IF (IT.NE.KERM.AND.IT.NE.KDIG(1)) THEN
                KCTCTL=1
                KCTCHR=JINT(JCOL)
                RETURN
           END IF
C     NO CHARACTER SPECIFIED OR ZERO.
           KCTCTL=0
           KCTCHR=KSPK(10)
           IF (IT.EQ.KDIG(1)) CALL DIAGNO (38)
      ELSE IF (J.EQ.30) THEN
C                            DEL1 (PRIMARY STRING DELIMITER)
           KDEL1=KBL
           KDEL1(1:1)=IT(1:1)
           KAPSTR=KDEL1(1:1)//KAT(2:2)
      ELSE IF (J.EQ.31) THEN
C                            DEL2 (SECONDARY STRING DELIMITER)
           KDEL2=KBL
           KDEL2(1:1)=IT(1:1)
      ELSE IF (J.EQ.32) THEN
C                            ARET (ALT. RETURNS IN CALLS)
           KALMRK=IT
      ELSE IF (J.EQ.33) THEN
C                            ARTR (TRANSLATE KALMRK TO THIS)
           KALTRN=IT
      END IF
      RETURN
C
 800  FORMAT (' CONTRL ASSIGNING VALUE ',I6,' TO INDEX ',I3)
      END
      BLOCK DATA CTLDAT
C
      COMMON /CONTDY/ KTRL1,KTRL2,KTRL3,KTRL4,KTRL5,KTRL6,KTRL7,KTRL8,
     1KTRL9,KTRL10,KTRL11,KTRL12,KTRL13,KTRL14,KTRL15,KTRL16,KTRL17,
     2KTRL18,KTRL19,KTRL20,KTRL21,KTRL22,KTRL23,KTRL24,KTRL25,KTRL26,
     3KTRL27,KTRL28,KTRL29,KTRL30,KTRL31,KTRL32,KTRL33,KTRL34,KTRL35,
     4KTRL36,KTRL37,KTRL38,KTRL39,KTRL40,KTRL41,KTRL42,KTRL43
      CHARACTER*2 KTRL1(4),KTRL2(4),KTRL3(4),KTRL4(4),KTRL5(4),KTRL6(4),
     1KTRL7(4),KTRL8(4),KTRL9(4),KTRL10(4),KTRL11(4),KTRL12(4),KTRL13(4)
     2,KTRL14(4),KTRL15(4),KTRL16(4),KTRL17(4),KTRL18(4),KTRL19(4),
     3KTRL20(4),KTRL21(4),KTRL22(4),KTRL23(4),KTRL24(4),KTRL25(4),
     4KTRL26(4),KTRL27(4),KTRL28(4),KTRL29(4),KTRL30(4),KTRL31(4),
     5KTRL32(4),KTRL33(4),KTRL34(4),KTRL35(4),KTRL36(4),KTRL37(4),
     6KTRL38(4),KTRL39(4),KTRL40(4),KTRL41(4),KTRL42(4),KTRL43(4)
C
C     /CONTDY/
      DATA KTRL1/'B','A','S','E'/
      DATA KTRL2/'I','D','I','N'/
      DATA KTRL3/'I','D','S','T'/
      DATA KTRL4/'R','O','U','T'/
      DATA KTRL5/'S','T','A','T'/
      DATA KTRL6/'C','A','R','D'/
      DATA KTRL7/'C','O','L','L'/
      DATA KTRL8/'C','O','M','M'/
      DATA KTRL9/'E','X','E','M'/
      DATA KTRL10/'L','A','B','E'/
      DATA KTRL11/'L','A','S','T'/
      DATA KTRL12/'L','I','S','T'/
      DATA KTRL13/'N','E','W','R'/
      DATA KTRL14/'R','E','F','E'/
      DATA KTRL15/'S','K','I','P'/
      DATA KTRL16/'S','T','O','P'/
      DATA KTRL17/'S','E','R','I'/
      DATA KTRL18/'R','I','G','H'/
      DATA KTRL19/'L','E','F','T'/
      DATA KTRL20/'C','O','L','U'/
      DATA KTRL21/'I','N','D','E'/
      DATA KTRL22/'D','E','B','U'/
      DATA KTRL23/'C','O','N','T'/
      DATA KTRL24/'E','N','D',' '/
      DATA KTRL25/'A','N','S','I'/
      DATA KTRL26/'F','E','N','D'/
      DATA KTRL27/'C','C','H','R'/
      DATA KTRL28/'H','T','R','A'/
      DATA KTRL29/'D','T','R','A'/
      DATA KTRL30/'D','E','L','1'/
      DATA KTRL31/'D','E','L','2'/
      DATA KTRL32/'A','R','E','T'/
      DATA KTRL33/'A','R','T','R'/
      DATA KTRL34/'B','L','A','N'/
      DATA KTRL35/'F','S','P','L'/
      DATA KTRL36/'H','L','O','G'/
      DATA KTRL37/'C','A','S','E'/
      DATA KTRL38/'U','C','A','S'/
      DATA KTRL39/'L','C','A','S'/
      DATA KTRL40/'E','N','D','O'/
      DATA KTRL41/'F','M','T','B'/
      DATA KTRL42/'M','A','X','C'/
      DATA KTRL43/'O','L','D','S'/
      END
      SUBROUTINE INITDY
C
C     INITIALIZE TIDY -- USED AT START AND WHEN *NEWR EXECUTED.
C
      INCLUDE 'tidy.inc'
C
      INDENT=3
      JUST=7
      K72=72
      KALMRK = '* '
      KALTRN= '  '
      KBKCOK=1
      KBLCMT=' @'
      KB15=0
      KCTCHR=KSPK(13)
      KCTCTL=1
      KD15=10
      KD79=1
      KDEL1 = ''' '
      KDEL2 = '""'
      KDTRAN=0
      KHTRAN=1
      KHLOG=0
      KPRIN=0
      KPUN=-1
      KFSPL=1
      MANSI=1
      MAXCNT=19
      MCASE=0
      MCOL=0
      MCOM=-1
      MCONT=0
      MEX=0
      MFMTB=0
      MLBL=0
      MLIST=0
      MNDOO=1
      MPRIN=0
      MPUN=-1
      MREF=0
      MRIT=5
      MSER=0
      NFEND=0
      NLHTRN=0
      NROUT=1
C     DEFAULT CASE TRANSLATIONS = UPPER
C
C     SUBROUTINE KWDWRT WRITES FORTRAN KEYWORDS IN DESIRED CASE.
C       CHANGE (-1) TO (-2) FOR DEFAULT TRANSLATION TO LOWER-CASE.
      CALL KWDWRT (-1)
C
C     SUBROUTINE KCTSET CONTROLS TRANSLATION OF EXECUTABLE CODE WHICH
C      IS NOT FORTRAN KEYWORDS.
C       CHANGE (1) TO (0) FOR DEFAULT TRANSLATION TO UPPER-CASE
      CALL KCTSET (1)
C
      RETURN
      END
      SUBROUTINE KWSCAN (JT,KSTCR)
      PARAMETER (NKST=86)
C     PARAMETER (NKST=83)
C
C     THIS ROUTINE SCANS FOR FORTRAN KEYWORDS, SETS JT TO CORRECT
C     TYPE IF FOUND, ELSE ZERO.
C
C     INPUT: IF JT = 0, SCANS WHOLE LIST
C               JT > 0, ONLY SCANS THAT WORD.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C
      DIMENSION KSTCR(5)
      COMMON /KSTCOM/ KST(10,NKST)
      CHARACTER*2 KST,WKSTR(10),KUPPER
      COMMON /KSTNUM/ KSTC(6,NKST)
C
      IF (JT.EQ.0) THEN
           NL=1
           NU=NKST
C     ZERO OUT KSTCR FOR NEW SCANS ONLY
           DO 10 I=1,5
                KSTCR(I)=0
10         CONTINUE
      ELSE
           NL=JT
           NU=JT
      END IF
C
C     MAKE UPPER-CASE COPY OF 10 CHARS (MAX STRING LENGTH)
      LAST=JCOL-1
      DO 30 I=1,10
20         LAST=LAST+1
           IF (LAST.GT.JMAX) THEN
                WKSTR(I)=KBL
           ELSE
                IF (JINT(LAST).EQ.KBL) GO TO 20
                WKSTR(I)=KUPPER(JINT(LAST))
           END IF
30    CONTINUE
      IF (MDEB.GT.0) WRITE (OUTFIL,70) WKSTR,JT
C
      DO 60 IT=NL,NU
           NINS=KSTC(1,IT)
C
           DO 40 I=1,NINS
                IF (WKSTR(I).NE.KST(I,IT)) GO TO 60
40         CONTINUE
           JT=KSTC(6,IT)
           DO 50 I=1,5
                KSTCR(I)=KSTC(I,IT)
50         CONTINUE
           IF (MDEB.GT.0) WRITE (OUTFIL,80) KSTCR,JT
           RETURN
C                  LOOP FOR NEXT STATEMENT.
60    CONTINUE
C
C     NO MATCH.
      IF (MDEB.GT.0) WRITE (OUTFIL,90)
      JT=0
C
      RETURN
C
C
70    FORMAT (' KWSCAN checking ',10A1,' mode = ',I2)
80    FORMAT ('   NINS  =',I3,' KLASS  =',I3,' JTYPE =',I3/'   NANSI =',
     1I3,' KSTROK =',I3,' KPOS  =',I3)
90    FORMAT ('  --- no match')
      END
      BLOCK DATA KSTDAT
C
      COMMON /KSTCOM/
     1     KST01 ,KST02 ,KST03 ,KST04 ,KST05
     2    ,KST06 ,KST07 ,KST08 ,KST09 ,KST10
     3    ,KST11 ,KST12 ,KST13 ,KST14 ,KST15
     4    ,KST16 ,KST17 ,KST18 ,KST19 ,KST20
     5    ,KST21 ,KST22 ,KST23 ,KST24 ,KST25
     6    ,KST26 ,KST27 ,KST28 ,KST29 ,KST30
     7    ,KST31 ,KST32 ,KST33 ,KST34 ,KST35
     8    ,KST36 ,KST37 ,KST38 ,KST39 ,KST40
     9    ,KST41 ,KST42 ,KST43 ,KST44 ,KST45
     X    ,KST46 ,KST47 ,KST48 ,KST49 ,KST50
     X    ,KST51 ,KST52 ,KST53 ,KST54 ,KST55
     X    ,KST56 ,KST57 ,KST58 ,KST59 ,KST60
     X    ,KST61 ,KST62 ,KST63 ,KST64 ,KST65
     X    ,KST66 ,KST67 ,KST68 ,KST69 ,KST70
     X    ,KST71 ,KST72 ,KST73 ,KST74 ,KST75
     X    ,KST76 ,KST77 ,KST78 ,KST79 ,KST80
     X    ,KST81 ,KST82 ,KST83 ,KST84 ,KST85
     X    ,KST86
C    X    ,KST81 ,KST82 ,KST83
C
C
      CHARACTER*2 KST01(10),KST02(10),KST03(10),KST04(10),KST05(10)
      CHARACTER*2 KST06(10),KST07(10),KST08(10),KST09(10),KST10(10)
      CHARACTER*2 KST11(10),KST12(10),KST13(10),KST14(10),KST15(10)
      CHARACTER*2 KST16(10),KST17(10),KST18(10),KST19(10),KST20(10)
      CHARACTER*2 KST21(10),KST22(10),KST23(10),KST24(10),KST25(10)
      CHARACTER*2 KST26(10),KST27(10),KST28(10),KST29(10),KST30(10)
      CHARACTER*2 KST31(10),KST32(10),KST33(10),KST34(10),KST35(10)
      CHARACTER*2 KST36(10),KST37(10),KST38(10),KST39(10),KST40(10)
      CHARACTER*2 KST41(10),KST42(10),KST43(10),KST44(10),KST45(10)
      CHARACTER*2 KST46(10),KST47(10),KST48(10),KST49(10),KST50(10)
      CHARACTER*2 KST51(10),KST52(10),KST53(10),KST54(10),KST55(10)
      CHARACTER*2 KST56(10),KST57(10),KST58(10),KST59(10),KST60(10)
      CHARACTER*2 KST61(10),KST62(10),KST63(10),KST64(10),KST65(10)
      CHARACTER*2 KST66(10),KST67(10),KST68(10),KST69(10),KST70(10)
      CHARACTER*2 KST71(10),KST72(10),KST73(10),KST74(10),KST75(10)
      CHARACTER*2 KST76(10),KST77(10),KST78(10),KST79(10),KST80(10)
      CHARACTER*2 KST81(10),KST82(10),KST83(10),KST84(10),KST85(10)
      CHARACTER*2 KST86(10)
C     CHARACTER*2 KST81(10),KST82(10),KST83(10)
C
      COMMON /KSTNUM/
     1     KSTC01 ,KSTC02 ,KSTC03 ,KSTC04 ,KSTC05
     2    ,KSTC06 ,KSTC07 ,KSTC08 ,KSTC09 ,KSTC10
     3    ,KSTC11 ,KSTC12 ,KSTC13 ,KSTC14 ,KSTC15
     4    ,KSTC16 ,KSTC17 ,KSTC18 ,KSTC19 ,KSTC20
     5    ,KSTC21 ,KSTC22 ,KSTC23 ,KSTC24 ,KSTC25
     6    ,KSTC26 ,KSTC27 ,KSTC28 ,KSTC29 ,KSTC30
     7    ,KSTC31 ,KSTC32 ,KSTC33 ,KSTC34 ,KSTC35
     8    ,KSTC36 ,KSTC37 ,KSTC38 ,KSTC39 ,KSTC40
     9    ,KSTC41 ,KSTC42 ,KSTC43 ,KSTC44 ,KSTC45
     X    ,KSTC46 ,KSTC47 ,KSTC48 ,KSTC49 ,KSTC50
     X    ,KSTC51 ,KSTC52 ,KSTC53 ,KSTC54 ,KSTC55
     X    ,KSTC56 ,KSTC57 ,KSTC58 ,KSTC59 ,KSTC60
     X    ,KSTC61 ,KSTC62 ,KSTC63 ,KSTC64 ,KSTC65
     X    ,KSTC66 ,KSTC67 ,KSTC68 ,KSTC69 ,KSTC70
     X    ,KSTC71 ,KSTC72 ,KSTC73 ,KSTC74 ,KSTC75
     X    ,KSTC76 ,KSTC77 ,KSTC78 ,KSTC79 ,KSTC80
     X    ,KSTC81 ,KSTC82 ,KSTC83 ,KSTC84 ,KSTC85
     X    ,KSTC86
C    X    ,KSTC81 ,KSTC82 ,KSTC83
      DIMENSION KSTC01(6),KSTC02(6),KSTC03(6),KSTC04(6),KSTC05(6)
      DIMENSION KSTC06(6),KSTC07(6),KSTC08(6),KSTC09(6),KSTC10(6)
      DIMENSION KSTC11(6),KSTC12(6),KSTC13(6),KSTC14(6),KSTC15(6)
      DIMENSION KSTC16(6),KSTC17(6),KSTC18(6),KSTC19(6),KSTC20(6)
      DIMENSION KSTC21(6),KSTC22(6),KSTC23(6),KSTC24(6),KSTC25(6)
      DIMENSION KSTC26(6),KSTC27(6),KSTC28(6),KSTC29(6),KSTC30(6)
      DIMENSION KSTC31(6),KSTC32(6),KSTC33(6),KSTC34(6),KSTC35(6)
      DIMENSION KSTC36(6),KSTC37(6),KSTC38(6),KSTC39(6),KSTC40(6)
      DIMENSION KSTC41(6),KSTC42(6),KSTC43(6),KSTC44(6),KSTC45(6)
      DIMENSION KSTC46(6),KSTC47(6),KSTC48(6),KSTC49(6),KSTC50(6)
      DIMENSION KSTC51(6),KSTC52(6),KSTC53(6),KSTC54(6),KSTC55(6)
      DIMENSION KSTC56(6),KSTC57(6),KSTC58(6),KSTC59(6),KSTC60(6)
      DIMENSION KSTC61(6),KSTC62(6),KSTC63(6),KSTC64(6),KSTC65(6)
      DIMENSION KSTC66(6),KSTC67(6),KSTC68(6),KSTC69(6),KSTC70(6)
      DIMENSION KSTC71(6),KSTC72(6),KSTC73(6),KSTC74(6),KSTC75(6)
      DIMENSION KSTC76(6),KSTC77(6),KSTC78(6),KSTC79(6),KSTC80(6)
      DIMENSION KSTC81(6),KSTC82(6),KSTC83(6),KSTC84(6),KSTC85(6)
      DIMENSION KSTC86(6)
C     DIMENSION KSTC81(6),KSTC82(6),KSTC83(6)
C
C     /KST/
      DATA KST01/'A','C','C','E','P','T',' ',' ',' ',' '/
      DATA KST02/'A','S','C','E','N','T',' ',' ',' ',' '/
      DATA KST03/'A','S','S','I','G','N',' ',' ',' ',' '/
      DATA KST04/'B','A','C','K','S','P','A','C','E','('/
      DATA KST05/'B','L','O','C','K','D','A','T','A',' '/
      DATA KST06/'B','U','F','F','E','R','I','N','(',' '/
      DATA KST07/'B','U','F','F','E','R','O','U','T','('/
      DATA KST08/'C','A','L','L',' ',' ',' ',' ',' ',' '/
      DATA KST09/'C','H','A','R','A','C','T','E','R',' '/
      DATA KST10/'C','O','M','M','O','N',' ',' ',' ',' '/
      DATA KST11/'C','O','M','P','L','E','X',' ',' ',' '/
      DATA KST12/'C','O','N','T','I','N','U','E',' ',' '/
      DATA KST13/'D','A','T','A',' ',' ',' ',' ',' ',' '/
      DATA KST14/'D','E','C','O','D','E','(',' ',' ',' '/
      DATA KST15/'D','I','M','E','N','S','I','O','N',' '/
      DATA KST16/'D','O','U','B','L','E','P','R','E','C'/
      DATA KST17/'D','O','U','B','L','E',' ',' ',' ',' '/
      DATA KST18/'E','N','C','O','D','E','(',' ',' ',' '/
      DATA KST19/'E','N','D','F','I','L','E','(',' ',' '/
      DATA KST20/'E','N','D','I','F',' ',' ',' ',' ',' '/
      DATA KST21/'E','N','D','F','I','L','E',' ',' ',' '/
      DATA KST22/'E','N','T','R','Y',' ',' ',' ',' ',' '/
      DATA KST23/'E','Q','U','I','V','A','L','E','N','C'/
      DATA KST24/'E','X','T','E','R','N','A','L',' ',' '/
      DATA KST25/'F','I','N','I','S',' ',' ',' ',' ',' '/
      DATA KST26/'F','O','R','M','A','T','(',' ',' ',' '/
      DATA KST27/'F','O','R','T','R','A','N',' ',' ',' '/
      DATA KST28/'I','F','(','U','N','I','T',',',' ',' '/
      DATA KST29/'F','U','N','C','T','I','O','N',' ',' '/
      DATA KST30/'G','O','T','O','(',' ',' ',' ',' ',' '/
      DATA KST31/'G','O','T','O',' ',' ',' ',' ',' ',' '/
      DATA KST32/'I','F','A','C','C','U','M','U','L','A'/
      DATA KST33/'I','F','Q','U','O','T','I','E','N','T'/
      DATA KST34/'I','F','(','D','I','V','I','D','E','C'/
      DATA KST35/'I','F','(','E','N','D','F','I','L','E'/
      DATA KST36/'I','F','(','S','E','N','S','E','L','I'/
      DATA KST37/'I','F','(','S','E','N','S','E','S','W'/
      DATA KST38/'I','F','(',' ',' ',' ',' ',' ',' ',' '/
      DATA KST39/'I','N','T','E','G','E','R',' ',' ',' '/
      DATA KST40/'L','O','G','I','C','A','L',' ',' ',' '/
      DATA KST41/'M','A','C','H','I','N','E',' ',' ',' '/
      DATA KST42/'N','A','M','E','L','I','S','T',' ',' '/
      DATA KST43/'P','A','U','S','E',' ',' ',' ',' ',' '/
      DATA KST44/'P','R','I','N','T',' ',' ',' ',' ',' '/
      DATA KST45/'P','R','O','G','R','A','M',' ',' ',' '/
      DATA KST46/'P','U','N','C','H',' ',' ',' ',' ',' '/
      DATA KST47/'R','E','A','D','I','N','P','U','T','T'/
      DATA KST48/'R','E','A','D','T','A','P','E',' ',' '/
      DATA KST49/'R','E','A','D','(',' ',' ',' ',' ',' '/
      DATA KST50/'R','E','A','D',' ',' ',' ',' ',' ',' '/
      DATA KST51/'R','E','A','L',' ',' ',' ',' ',' ',' '/
      DATA KST52/'R','E','T','U','R','N',' ',' ',' ',' '/
      DATA KST53/'R','E','W','I','N','D','(',' ',' ',' '/
      DATA KST54/'S','E','G','M','E','N','T',' ',' ',' '/
      DATA KST55/'S','E','N','S','E','L','I','G','H','T'/
      DATA KST56/'S','T','O','P',' ',' ',' ',' ',' ',' '/
      DATA KST57/'S','U','B','R','O','U','T','I','N','E'/
      DATA KST58/'T','Y','P','E',' ',' ',' ',' ',' ',' '/
      DATA KST59/'W','R','I','T','E','O','U','T','P','U'/
      DATA KST60/'W','R','I','T','E','T','A','P','E',' '/
      DATA KST61/'W','R','I','T','E','(',' ',' ',' ',' '/
      DATA KST62/'O','V','E','R','L','A','Y',' ',' ',' '/
      DATA KST63/'I','D','E','N','T',' ',' ',' ',' ',' '/
      DATA KST64/'F','R','E','Q','U','E','N','C','Y',' '/
      DATA KST65/'I','M','P','L','I','C','I','T',' ',' '/
      DATA KST66/'L','E','V','E','L',' ',' ',' ',' ',' '/
      DATA KST67/'E','L','S','E','I','F',' ',' ',' ',' '/
      DATA KST68/'E','L','S','E',' ',' ',' ',' ',' ',' '/
      DATA KST69/'T','H','E','N',' ',' ',' ',' ',' ',' '/
      DATA KST70/'C','L','O','S','E','(',' ',' ',' ',' '/
      DATA KST71/'I','N','C','L','U','D','E',' ',' ',' '/
      DATA KST72/'I','N','Q','U','I','R','E','(',' ',' '/
      DATA KST73/'I','N','T','R','I','N','S','I','C',' '/
      DATA KST74/'O','P','E','N','(',' ',' ',' ',' ',' '/
      DATA KST75/'P','A','R','A','M','E','T','E','R',' '/
      DATA KST76/'S','A','V','E',' ',' ',' ',' ',' ',' '/
      DATA KST77/'B','A','C','K','S','P','A','C','E',' '/
      DATA KST78/'E','N','D','D','O',' ',' ',' ',' ',' '/
      DATA KST79/'R','E','W','I','N','D',' ',' ',' ',' '/
      DATA KST80/'C','L','O','S','E',' ',' ',' ',' ',' '/
      DATA KST81/'E','N','D',' ',' ',' ',' ',' ',' ',' '/
      DATA KST82/'D','O','W','H','I','L','E','(',' ',' '/
      DATA KST83/'R','E','P','E','A','T',' ',' ',' ',' '/

      DATA KST84/'C','Y','C','L','E',' ',' ',' ',' ',' '/
      DATA KST85/'E','X','I','T',' ',' ',' ',' ',' ',' '/
      DATA KST86/'N','O','N','E',' ',' ',' ',' ',' ',' '/
C
C     /KSTNUM/
C     ********* NOTE - KPOS IS ADDED TO INSULATE PASS1 FROM ADDITIONS
C     TO ABOVE TABLE.  WHEN ADDING NEW STATEMENTS, SET KPOS TO THE
C     NEW VALUE OF NKST RATHER THAN THE ORDINAL POSITION OF THE NEW
C     ADDITION TO THE TABLE.
C      (NOTE WHEN ADDING - SIMILAR STRINGS MUST BE IN DESCENDING ORDER
C       BY LENGTH, I.E. END MUST FOLLOW ENDIF)
C     WARNING - DO NOT MOVE LINES 69 OR 82 WITHOUT ALTERING PASS1 -
C               THERE ARE EXPLICIT REFERENCES TO THESE LINES.
C
C                KLASS  DESCRIPTION
C                  0.   CONTROL CARD
C                  1.   COMMENT
C                  2.   HEADER
C                  3.   NO STATEMENT NO ALLOWED (NON-EXECTUABLE)
C                  4.   CONTINUE
C                  5.   FORMAT STATEMENT.
C                  6.   STATEMENT NO. ALLOWED, NO REFERENCES
C                  7.   REFERENCES PRESENT, STATEMENT NO. ALLOWED.
C                  8.   END
C                  9.   INTRODUCTORY
C                  10.  DO
C                  11.  ELSE,ENDIF,ELSEIF, UNRECOGNIZED
C                       (TRANSFER CAN GET HERE REGARDLESS OF LABEL)
C
C     KLASS 0.   CONTROL CARD
C             RESERVED FOR FUTURE DEVELOPMENT.
C
C
C                   NINS  KLASS  JTYPE NANSI   KSTROK     KPOS
      DATA KSTC01 /    6,     7,    33,    1,       0,        1/
      DATA KSTC02 /    6,     2,     1,    1,       0,        2/
      DATA KSTC03 /    6,     7,     2,    0,       0,        3/
      DATA KSTC04 /   10,     7,    47,    0,       0,        4/
      DATA KSTC05 /    9,     2,     4,    0,       0,        5/
      DATA KSTC06 /    9,     6,     5,    1,       0,        6/
      DATA KSTC07 /   10,     6,     5,    1,       0,        7/
      DATA KSTC08 /    4,     7,     6,    0,       1,        8/
      DATA KSTC09 /    9,     3,    46,    0,       0,        9/
      DATA KSTC10 /    6,     3,     7,    0,       0,       10/
      DATA KSTC11 /    7,     3,    46,    0,       0,       11/
      DATA KSTC12 /    8,     4,     8,    0,       0,       12/
      DATA KSTC13 /    4,     3,     9,    0,       1,       13/
      DATA KSTC14 /    7,     7,    10,    1,       0,       14/
      DATA KSTC15 /    9,     3,    11,    0,       0,       15/
      DATA KSTC16 /   10,     3,    12,    0,       0,       16/
      DATA KSTC17 /    6,     3,    13,    0,       0,       17/
      DATA KSTC18 /    7,     7,    10,    1,       0,       18/
      DATA KSTC19 /    8,     7,    47,    0,       0,       19/
      DATA KSTC20 /    5,    11,    48,    0,       0,       20/
      DATA KSTC21 /    7,     6,    15,    0,       0,       21/
      DATA KSTC22 /    5,    11,     3,    0,       0,       22/
      DATA KSTC23 /   10,     3,    17,    0,       0,       23/
      DATA KSTC24 /    8,     3,     3,    0,       0,       24/
      DATA KSTC25 /    5,     3,    18,    1,       0,       25/
      DATA KSTC26 /    7,     5,    19,    0,       1,       26/
      DATA KSTC27 /    7,     2,    20,    1,       0,       27/
      DATA KSTC28 /    8,     7,    42,    1,       1,       28/
      DATA KSTC29 /    8,     2,    35,    0,       0,       29/
      DATA KSTC30 /    5,     7,    23,    0,       0,       30/
      DATA KSTC31 /    4,     7,    24,    0,       0,       31/
      DATA KSTC32 /   10,     7,    25,    1,       1,       32/
      DATA KSTC33 /   10,     7,    26,    1,       1,       33/
      DATA KSTC34 /   10,     7,    27,    1,       1,       34/
      DATA KSTC35 /   10,     7,    28,    1,       1,       35/
      DATA KSTC36 /   10,     7,    29,    1,       1,       36/
      DATA KSTC37 /   10,     7,    30,    1,       1,       37/
      DATA KSTC38 /    3,     7,    31,    0,       1,       38/
      DATA KSTC39 /    7,     3,    46,    0,       0,       39/
      DATA KSTC40 /    7,     3,    46,    0,       0,       40/
      DATA KSTC41 /    7,     2,     1,    1,       0,       41/
      DATA KSTC42 /    8,     3,    32,    1,       0,       42/
      DATA KSTC43 /    5,     6,     3,    0,       1,       43/
      DATA KSTC44 /    5,     7,    33,    0,       1,       44/
      DATA KSTC45 /    7,     2,    35,    0,       0,       45/
      DATA KSTC46 /    5,     7,    33,    1,       1,       46/
      DATA KSTC47 /   10,     7,    36,    0,       0,       47/
      DATA KSTC48 /    8,     6,    37,    0,       0,       48/
      DATA KSTC49 /    5,     7,    38,    0,       1,       49/
      DATA KSTC50 /    4,     7,    33,    0,       1,       50/
      DATA KSTC51 /    4,     3,    46,    0,       0,       51/
      DATA KSTC52 /    6,     6,    39,    0,       0,       52/
      DATA KSTC53 /    7,     7,    47,    0,       0,       53/
      DATA KSTC54 /    7,     9,    34,    1,       0,       54/
      DATA KSTC55 /   10,     6,    40,    1,       0,       55/
      DATA KSTC56 /    4,     6,    41,    0,       1,       56/
      DATA KSTC57 /   10,     2,    35,    0,       0,       57/
      DATA KSTC58 /    4,     7,    33,    1,       0,       58/
      DATA KSTC59 /   10,     7,    44,    0,       1,       59/
      DATA KSTC60 /    9,     6,    45,    0,       1,       60/
      DATA KSTC61 /    6,     7,    38,    0,       1,       61/
      DATA KSTC62 /    7,     9,    34,    1,       0,       62/
      DATA KSTC63 /    5,     9,    22,    1,       0,       63/
      DATA KSTC64 /    9,     3,    21,    1,       0,       64/
      DATA KSTC65 /    8,     3,     3,    0,       0,       65/
      DATA KSTC66 /    5,     3,     3,    1,       0,       66/
      DATA KSTC67 /    6,    11,    43,    0,       1,       67/
      DATA KSTC68 /    4,    11,    49,    0,       0,       68/
      DATA KSTC69 /    4,    11,     3,    0,       0,       69/
      DATA KSTC70 /    6,     7,    47,    0,       0,       70/
      DATA KSTC71 /    7,     3,     3,    1,       1,       71/
      DATA KSTC72 /    8,     7,    47,    0,       1,       72/
      DATA KSTC73 /    9,     3,     3,    0,       0,       73/
      DATA KSTC74 /    5,     7,    47,    0,       1,       74/
      DATA KSTC75 /    9,     3,     3,    0,       1,       75/
      DATA KSTC76 /    4,     3,     3,    0,       0,       76/
      DATA KSTC77 /    9,     6,     3,    0,       0,       77/
      DATA KSTC78 /    5,     7,    50,    1,       1,       81/
      DATA KSTC79 /    6,     6,     3,    0,       0,       79/
      DATA KSTC80 /    5,     6,     3,    0,       0,       80/
      DATA KSTC81 /    3,     8,    16,    0,       0,       78/
      DATA KSTC82 /    8,    11,    51,    1,       0,       82/
      DATA KSTC83 /    6,     7,    50,    1,       1,       83/

      DATA KSTC84 /    5,     6,     3,    1,       0,       88/
      DATA KSTC85 /    4,     6,     3,    1,       0,       89/
      DATA KSTC86 /    4,     3,    46,    1,       0,       90/
C                   NINS  KLASS  JTYPE NANSI   KSTROK     KPOS
      END
      SUBROUTINE ADNUM (JERR)
      INCLUDE 'tidy.inc'
C     ADDS STATEMENT NUMBER TO CROSS-REF LIST, AND PLACES SPECIAL
C     MARK IN OUTPUT CARD IMAGE.
C
C     JERR = RETURN CODE 0 = NORMAL TERMINATION
C                        1 = CROSS-REF TABLE OVERFLOW
      ICOL=ICOL+1
      IOUT(ICOL)=KLR2
      IF (NXRF.LE.MXREF) THEN
           JERR=0
           IOUTN(NXRF)=INT(L772)
           NXRF=NXRF+1
           CALL RLIST
      ELSE
           JERR=1
C          TOO MANY CROSS-REFERENCES
           CALL DIAGNO (35)
           MP2=0
      ENDIF
      RETURN
      END
      LOGICAL FUNCTION BAKSCN (C1,C2)
C
C     SCANS A STRING BACKWARD FROM CURRENT POSITION FOR C1 AND C2
      CHARACTER*2 C1, C2, JT, KUPPER, JNT
      INCLUDE 'tidy.inc'
      IP = JCOL
C     FIRST BACK TO LCPY
    5 IF (JINT(IP).NE.LCPY) THEN
           IP = IP-1
           GO TO 5
      END IF
C
C     NOW SCAN FOR C1, C2
      JT = C1
      I = 1
   15 IP = IP-1
      JNT=KUPPER(JINT(IP))
      IF (JNT.EQ.KBL) GO TO 15
      IF (JNT.NE.JT) THEN
           BAKSCN = .FALSE.
           RETURN
      ENDIF
      IF (I.EQ.1) THEN
           JT = C2
           I = 2
           GO TO 15
      ENDIF
      BAKSCN = .TRUE.
      RETURN
      END
      SUBROUTINE COPY (N)
C
C     COPY NON-BLANK CHARACTERS FROM JINT TO IOUT.
C       (UNLESS *EXEM IS SET, THEN COPY BLANKS ALSO)
C
C                        ===   ON ENTRY   ===
C     N .LT. 0 COPIES UNTIL PARENTHESIS COUNT IS ZERO.
C     N .EQ. 0 COPIES ALL REMAINING NON-BLANK DATA FROM JINT TO IOUT.
C     N .GT. 0 COPIES N NON-BLANK DATA FROM JINT TO IOUT.
C     THE FIRST ITEM INSPECTED IS JINT(JCOL).
C     THE FIRST ITEM STORED GOES TO IOUT(ICOL+1).
C
C                        ===   ON EXIT   ===
C     THE LAST ITEM INSPECTED WAS JINT(JCOL-1).
C     THE LAST ITEM STORED WENT TO IOUT(ICOL) AND IS IN LCPY.
C
C     MEOF .LT. 0  FOR NORMAL EXIT.
C     MEOF .EQ. 0  FOR KERM FOUND WHILE COPYING  ALL REMAINING DATA,
C                  OR FOR KERM FOUND BEFORE LEFT PARENTHESIS.
C     MEOF .GT. 0  FOR MISSING RIGHT PARENTHESIS, OR FOR MEOF =0 ON
C                  ENTRY TO COPY.
C
      INCLUDE 'tidy.inc'
      CHARACTER*2 JT
      LOGICAL SAVBLK
C
      IF (MEOF.GE.0.OR.JCOL.GT.JMAX) THEN
          MEOF=1
          LCPY=KERM
          RETURN
      END IF
C
C     SET BLANK STRIP MODE
      SAVBLK=(MEX.GT.0 .OR. (MEX.LT.0.AND.(KLASS.EQ.3.OR.KLASS.EQ.5)))
C
      NT=N
      IF (NT.EQ.0) THEN
C
C     COPY ALL REMAINING NON-BLANK CHARACTERS.
C
10        JT=JINT(JCOL)
          IF (JT.NE.KBL.OR.SAVBLK) THEN
              ICOL=ICOL+1
              IOUT(ICOL)=JT
          END IF
          IF (JT.NE.KERM) THEN
              JCOL=JCOL+1
              GO TO 10
          END IF
          GO TO 70
C
      ELSE IF (NT.GT.0) THEN
C
C     COPY --N-- NON-BLANK CHARACTERS.
C
20        JT=JINT(JCOL)
          IF (JT.NE.KBL) THEN
              ICOL=ICOL+1
              IOUT(ICOL)=JT
              NT=NT-1
              IF (NT.EQ.0) GO TO 80
              IF (JT.EQ.KERM) GO TO 70
          END IF
          JCOL=JCOL+1
          GO TO 20
      ELSE
C
C     COPY TO PARENTHESIS COUNT OF ZERO.
C     LOOK FOR LEFT PARENTHESIS.
C
30        JT=JINT(JCOL)
          IF (JT.NE.KBL) THEN
              ICOL=ICOL+1
              IOUT(ICOL)=JT
              LCPY=JT
              IF (JT.EQ.KSPK(3)) THEN
C        HAVE LEFT PARENTHESIS, COPY UNTIL COUNT OF ZERO.
                  NPAR=1
40                JCOL=JCOL+1
                  JT=JINT(JCOL)
                  IF (JT.NE.KBL) THEN
                      ICOL=ICOL+1
                      IOUT(ICOL)=JT
                      LCPY=JT
                      IF (JT.NE.KSPK(3)) THEN
                          IF (JT.NE.KSPK(5)) THEN
                              IF (JT.NE.KERM) GO TO 40
                              CALL DIAGNO (2)
                              LCPY=KERM
                              GO TO 60
                          END IF
                          NPAR=NPAR-1
C                         IF (NPAR) 50,80,40
                          IF (NPAR.LT.0) THEN
                             GO TO 50
                          ELSE IF (NPAR.EQ.0) THEN
                             GO TO 80
                          ELSE
                             GO TO 40
                          END IF 
                      END IF
                      NPAR=NPAR+1
                  ELSE IF (SAVBLK) THEN
                      ICOL=ICOL+1
                      IOUT(ICOL)=JT
                  END IF
                  GO TO 40
              END IF
              IF (JT.EQ.KSPK(5)) GO TO 50
              IF (JT.EQ.KERM) GO TO 70
          ELSE IF (SAVBLK) THEN
              ICOL=ICOL+1
              IOUT(ICOL)=JT
          END IF
          JCOL=JCOL+1
          GO TO 30
C
50        CALL DIAGNO (3)
60        MEOF=1
          JCOL=JCOL+1
          RETURN
      END IF
C
70    LCPY=KERM
      ICOL=ICOL-1
      MEOF=0
      RETURN
C
80    JCOL=JCOL+1
      LCPY=JT
      RETURN
      END
      SUBROUTINE DIAGNO (N)
      PARAMETER (MXMSG=53)
C
C     THIS ROUTINE WRITES THE GENERAL DIAGNOSTICS FOR TIDY.
C
      DIMENSION LV(MXMSG)
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C     ***                                                            ***
C      1 THE ABOVE STATEMENT IS ILLEGAL AND HAS BEEN DELETED.
C      2 THE ABOVE STATEMENT HAS A MISSING RIGHT PARENTHESIS.
C      3 THE ABOVE STATEMENT HAS AN EXCESS RIGHT PARENTHESIS.
C      4 THE ABOVE STATEMENT INCORRECTLY TERMINATES A DO LOOP.
C      5 THE ABOVE STATEMENT CANNOT BE REACHED BY THE PROGRAM.
C      6 STATEMENT NUMBER TABLE FULL.  RENUMBER PASS DELETED.
C      7 REFERENCE NUMBER TABLE FULL.  RENUMBER PASS DELETED.
C      8 THE ABOVE STATEMENT IS OBSOLETE AND IS DELETED.
C      9 ABOVE STATEMENT HAS AN ILLEGAL FIRST SPECIAL CHARACTER.
C     10 ILLEGAL DATA, FUNCTION, OR SUBROUTINE STATEMENT.
C     11 THE ABOVE COMMON OR DATA STATEMENT IS MISSING A (/).
C     12 THE ABOVE CONTINUE STATEMENT IS REDUNDANT AND IS DELETED.
C     13 THE ABOVE DIMENSION STATEMENT IS NOT COMPLETE.
C     14 W A R N I N G .  STATEMENT SHOULD BE FIRST IN ROUTINE.
C     15 THE ABOVE DO STATEMENT HAS AN INVALID TERMINAL STATEMENT.
C     16 W A R N I N G .  UNSATISFIED DO LOOPS.
C     17 UNNUMBERED OR INVALID FORMAT STATEMENT DELETED.
C     18 WARNING.  ABOVE STATEMENT IS POOR PROGRAMMING PRACTICE.
C     19 ABOVE GO TO STATEMENT IS ILLEGAL.
C     20 ILLEGAL ARITHMETIC IF STATEMENT.   IF (ARITH) 1,2,3
C     21 ABOVE NAMELIST STATEMENT MISSING (/).
C     22 ILLEGAL READ, WRITE , OR PUNCH STATEMENT.
C     23 ILLEGAL READ (12) LIST, OR WRITE (12) LIST, STATEMENT.
C     24 DO LOOP TABLE FULL.  RENUMBER PASS DELETED.
C     25 W A R N I N G .   COMMA FOLLOWING X INSERTED IN ABOVE FORMAT.
C     26 TIDY CANNOT PROCESS THIS CLASS OF PROGRAM.  (COPY EXECUTED.)
C     27 WARNING.  ABOVE DO-LOOP TERMINUS PREVIOUSLY REFERENCED.
C     28 WARNING.  TIDY MAY HAVE CHANGED CARD 2 OF THIS ROUTINE
C     29 W A R N I N G .   END CARD INSERTED.
C     30 THE ABOVE STATEMENT IS TRANSMITTED WITHOUT PROCESSING
C     31 ILLEGAL CLOSE, INQUIRE, OR OPEN STATEMENT
C     32 W A R N I N G .   UNBALANCED ELSE/ELSEIF/ENDIF STATEENT
C     33 W A R N I N G .   UNSATISFIED IF BLOCKS.
C     34 W A R N I N G .   ABOVE STATEMENT NOT ANSI FORTRAN 77
C     35 TOO MANY REFERENCES IN ABOVE. RENUMBER PASS DELETED.
C     36 W A R N I N G .   NON-ANSI (L OR R) HOLLERITH SPEC.
C     37 ABOVE STATEMENT HAS MORE THAN 19 CONTINUATION LINES.
C     38 CCHR CARD IGNORED:   CANNOT USE ZERO.
C     39 >>> HOLLERITH CONSTANT CONVERTED <<<
C     40 W A R N I N G.   *PRECISION ON NUMERIC/LOGICAL VARS NOT ANSI
C     41 W A R N I N G.    VARIABLE NAME LONGER THAN 6 CHARACTERS
C     42 W A R N I N G.    INITIALIZED TYPE DECLARATIONS NOT ANSI
C     43 MORE <END DO> THAN <DO> STATEMENTS
C     44 FATAL ERROR - DO LIST UNDERFLOW
C     45 FATAL ERROR
C     46 FATAL PROBLEM IN DO-LOOP RENUMBERING - SUBROUTINE EDIT
C     47 ABNORMAL TERMINATION
C     48 CONTRL INTERNAL ERROR
C     49 *FMTB CONFLICTS WITH REGULAR STATEMENTS - IGNORED.
C     50 TOO MANY CONTINUATION CARDS - RECOMPILE
C     51 SOMETHING IN COLS 73-80
C     52 BLANKS REMOVED FROM STRING INCLUDING 73-80
C     53 REFERENCE TO MISSING STATEMENT NUMBER
C
      CHARACTER*60 ERMSG (MXMSG)
      DATA (ERMSG(I),I=1,15)/
     1'THE ABOVE STATEMENT IS ILLEGAL AND HAS BEEN DELETED.',
     1'THE ABOVE STATEMENT HAS A MISSING RIGHT PARENTHESIS.',
     1'THE ABOVE STATEMENT HAS AN EXCESS RIGHT PARENTHESIS.',
     1'THE ABOVE STATEMENT INCORRECTLY TERMINATES A DO LOOP.',
     1'THE ABOVE STATEMENT CANNOT BE REACHED BY THE PROGRAM.',
     1'STATEMENT NUMBER TABLE FULL.  RENUMBER PASS DELETED.',
     1'REFERENCE NUMBER TABLE FULL.  RENUMBER PASS DELETED.',
     1'THE ABOVE STATEMENT IS OBSOLETE AND IS DELETED.',
     1'ABOVE STATEMENT HAS AN ILLEGAL FIRST SPECIAL CHARACTER.',
     1'ILLEGAL DATA, FUNCTION, OR SUBROUTINE STATEMENT.',
     1'THE ABOVE COMMON OR DATA STATEMENT IS MISSING A (/).',
     1'THE ABOVE CONTINUE STATEMENT IS REDUNDANT AND IS DELETED.',
     1'THE ABOVE DIMENSION STATEMENT IS NOT COMPLETE.',
     1'W A R N I N G .  STATEMENT SHOULD BE FIRST IN ROUTINE.',
     1'THE ABOVE DO STATEMENT HAS AN INVALID TERMINAL STATEMENT.'/
      DATA (ERMSG(I),I=16,30)/
     1'W A R N I N G .  UNSATISFIED DO LOOPS.',
     1'UNNUMBERED OR INVALID FORMAT STATEMENT DELETED.',
     1'WARNING.  ABOVE STATEMENT IS POOR PROGRAMMING PRACTICE.',
     1'ABOVE GO TO STATEMENT IS ILLEGAL.',
     1'ILLEGAL ARITHMETIC IF STATEMENT.   IF (ARITH) 1,2,3',
     1'ABOVE NAMELIST STATEMENT MISSING (/).',
     1'ILLEGAL READ, WRITE , OR PUNCH STATEMENT.',
     1'ILLEGAL READ (12) LIST, OR WRITE (12) LIST, STATEMENT.',
     1'DO LOOP TABLE FULL.  RENUMBER PASS DELETED.',
     1'W A R N I N G .  COMMA INSERTED FOLLOWING X IN ABOVE FORMAT.',
     1'TIDY CANNOT PROCESS THIS CLASS OF PROGRAM.  (COPY EXECUTED.)',
     1'WARNING.  ABOVE DO-LOOP TERMINUS PREVIOUSLY REFERENCED.',
     1'WARNING.  TIDY MAY HAVE CHANGED CARD 2 OF THIS ROUTINE',
     1'W A R N I N G .  END CARD INSERTED.',
     1'THE ABOVE STATEMENT IS TRANSMITTED WITHOUT PROCESSING.'/
      DATA (ERMSG(I),I=31,45)/
     1'ILLEGAL CLOSE, INQUIRE, OR OPEN STATEMENT',
     1'W A R N I N G .   UNBALANCED ELSE/ELSEIF/ENDIF STATEMENT',
     1'W A R N I N G .   UNSATISFIED IF BLOCKS.',
     1'W A R N I N G .   ABOVE STATEMENT NOT ANSI FORTRAN 77.',
     1'TOO MANY REFERENCES IN ABOVE. RENUMBER PASS DELETED.',
     1'W A R N I N G .   NON-ANSI (L OR R) HOLLERITH SPEC.',
     1'ABOVE STATEMENT HAS TOO MANY CONTINUATION LINES.',
     1'CCHR CARD IGNORED:   CANNOT USE ZERO.',
     1'>>> HOLLERITH CONSTANT CONVERTED <<<',
     1'W A R N I N G. *n PRECISION ON NUMERIC/LOGICAL VARS NOT ANSI',
     1'W A R N I N G.    VARIABLE NAME LONGER THAN 6 CHARACTERS',
     1'W A R N I N G.    INITIALIZED TYPE DECLARATIONS NOT ANSI',
     1'MORE <END DO> THAN <DO> STATEMENTS',
     1'FATAL ERROR - DO LIST UNDERFLOW',
     1'FATAL ERROR'/
      DATA (ERMSG(I),I=46,MXMSG)/
     1'FATAL PROBLEM IN DO-LOOP RENUMBERING - SUBROUTINE EDIT',
     1'ABNORMAL TERMINATION',
     1'CONTRL INTERNAL ERROR - A COMPUTED GO TO LIST IS TOO SHORT',
     1'*FMTB CONFLICTS WITH REGULAR STATEMENTS - IGNORED.',
     1'TOO MANY CONTINUATION CARDS REQUESTED - RECOMPILE',
     1'COLUMNS 73-80 IN ABOVE ADDED TO STATEMENT - CHECK CAREFULLY',
     1'BLANKS REMOVED FROM STRING INCLUDING COLS 73-80',
     1'REFERENCE TO MISSING STATEMENT NUMBER IN NEXT LINE'/
C
C     LV=0 - TIDY USER WARNING - CAUSES NORMAL TERMINATION
C        1 - MINOR FORTRAN ERROR - STOP 1
C        2 - MAJOR FORTRAN ERROR - STOP 2
C        3 - IMMEDIATELY FATAL   - STOP 3
C
C       -1 - TERMINATE WITH PREVIOUS HIGHEST ERROR LEVEL
C
      DATA LV /2,2,2,2,1 ,2,2,2,2,2 ,2,1,2,1,2 ,2,1,1,2,2
     1        ,2,2,2,2,0 ,0,0,1,1,1 ,2,1,2,0,2 ,0,3,0,0,0
     2        ,0,0,2,3,3 ,3,-1,3,0,3 ,2,2,2/
C
      J=N
      IF (J.LE.0.OR.J.GT.MXMSG) J=1
      NMSG=NMSG+1
      IF (LERR.LT.LV(J)) LERR=LV(J)
      IF (MLIST.EQ.-1) THEN
         CALL PAGE (1)
      ELSE
         CALL PAGE ((JMAX-7)/66+4)
         WRITE (OUTFIL,320) (JINT(I),I=1,JMAX)
      END IF
      WRITE (OUTFIL,340) NMSG, ERMSG(J)
C
      IF (MLIST.NE.-1) WRITE (OUTFIL,330) NREC,KBUFF
C
*      IF (LERR.GE.3) STOP 3
      IF (LERR.GE.3) CALL QUIT(3)
      IF (LV(J).LT.0) THEN
*           IF (LERR.EQ.2) STOP 2
*           IF (LERR.EQ.1) STOP 1
           IF (LERR.EQ.2) CALL QUIT(2)
           IF (LERR.EQ.1) CALL QUIT(1)
      END IF
      RETURN
C
C
 320  FORMAT (7X,72A1,19(/12X,'X',66A1))
C330  FORMAT (1X,I4,2X,80A1,/' ')
 330  FORMAT (1X,I4,2X,80A1)
 340  FORMAT (' ******(',I3,') ***',A60,'******',20X,'**********')
      END
      SUBROUTINE DLIST (MERR)
C
C     THIS SUBROUTINE UPDATES THE DEFINED STATEMENT NUMBER LIST, LDEF,
C     BY ADDING THE STATEMENT NUMBER IN L15, IF IT IS UNIQUE.
C              RETURNS MERR = 0 IF LABEL IS OK.
C                             1 IF LABEL TERMINATED A DO-LOOP (OK)
C                            -1 IF ERROR
C                       POSSIBLE ERRORS--
C                            ILLEGAL DO-LOOP NEST
C                            DUPLICATE STATEMENT NUMBER
C                            STATEMENT NUMBER TABLE FULL
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
      MERR=0
      DATA JTYPP/0/
      IF (KLASS.LT.4) THEN
          JTYPP=JTYPE
          RETURN
      END IF
C
C     CHECK FOR FORMAT STATEMENT, WHICH IS LABELED BUT CAN'T HAVE
C      FALL-THRU
      IF (KLASS.EQ.5) THEN
C          PROCESS FORMAT STATEMENT
C           SCAN FOR DUPLICATE STATEMENT NUMBER
          IF (NDEF.GT.0) THEN
              DO 10 I=1,NDEF
                  IF (abs(LDEF(I)).EQ.L15) GO TO 60
 10           CONTINUE
          END IF
C
C          PUT L15 INTO LDEF LIST AFTER LAST NON-NEGATIVE ENTRY
          IF (NDEF.GE.MXPLBL) GO TO 70
          I=NDEF
          NDEF=NDEF+1
 20       IF (I.EQ.0.OR.LDEF(I).GE.0) THEN
              LDEF(I+1)=L15
              LOCDEF(I+1)=NREC
              NNMKLS(I+1)=5
              GO TO 90
          END IF
          LDEF(I+1)=LDEF(I)
          LOCDEF(I+1)=LOCDEF(I)
          NNMKLS(I+1)=NNMKLS(I)
          I=I-1
          GO TO 20
      END IF
C
C     EXECUTABLE STATEMENT (OR END)
      IF (L15.EQ.0) THEN
C          UNLABELLED. IS THERE A FALL-THRU...
          IF (L25.EQ.0) THEN
C
C               UNLABELLED STATEMENT. ERROR IF IT FOLLOWS TRANSFER
C                (EXCEPT COMPUTED GO TO)
              IF (NTRAN.NE.0.AND.JTYPP.NE.23) CALL DIAGNO (5)
          ELSE
C               THERE IS A FALL-THRU LABEL. USE IT.
              L15=L25
              L25=0
              LDEF(NDEF)=abs(LDEF(NDEF))
          END IF
          GO TO 90
      END IF
C               LABELLED. SCRATCH FALL-THRU LABEL
      L25=0
C
C     SCAN FOR DUPLICATE STATEMENT NUMBERS.
C
      IF (NDEF.GT.0) THEN
          DO 30 I=1,NDEF
              IF (abs(LDEF(I)).EQ.L15) GO TO 60
 30       CONTINUE
      END IF
C
      IF (NDEF.GE.MXPLBL) GO TO 70
      NDEF=NDEF+1
      LDEF(NDEF)=L15
      LOCDEF(NDEF)=NREC
      NNMKLS(NDEF)=KLASS
      IF (MDEB.GT.0) WRITE (OUTFIL,100) NDEF,LDEF(NDEF),LOCDEF(NDEF),
     1NNMKLS(NDEF)
C
C     SCAN FOR POSSIBLE DO-LOOP TERMINATIONS.
C
      IF (NDOS.GT.0) THEN
           DO 50 I=1,NDOS
               IF (LDOS(I).EQ.L15) THEN
C                                 ITS IN THE LIST
                   MERR=1
                   IF (I.NE.NDOS) THEN
C                                 ILLEGAL DO-LOOP NEST
                       NMSG=NMSG+1
                       CALL PAGE (1)
                       WRITE (OUTFIL,110) NMSG,I,NDOS
C
C          COMPRESS DO-LOOP TERMINAL LIST AFTER DELETIONS.
C
                       NDOS=NDOS-1
                       DO 40 J=I,NDOS
                           LDOS(J)=LDOS(J+1)
 40                    CONTINUE
                       GO TO 80
                   END IF
C                                 LAST ONE IN LIST. REMOVE IT
                   NDOS=NDOS-1
                   IF (MILDO.NE.0) CALL DIAGNO (4)
                   GO TO 90
               END IF
 50        CONTINUE
      END IF
      GO TO 90
C
C     ERROR DIAGNOSTICS.
C
C                            DUPLICATE STATEMENT NUMBER
 60   NMSG=NMSG+1
      CALL PAGE (1)
      WRITE (OUTFIL,120) NMSG,L15,LOCDEF(I)
      GO TO 80
C                            NUMBER TABLE FULL
 70   CALL DIAGNO (6)
      NDEF=-1
      MP2=0
C                            ERROR EXIT
 80   MPUN=0
      MERR=-1
C                            EXIT
 90   MILDO=0
      NXEQ=NXEQ+1
      JTYPP=JTYPE
      RETURN
C
C
C
 100  FORMAT (' DLIST - NDEF = ',I4,' LDEF()= ',I4,' LOCDEF()=',I4,' NNM
     1KLS()= ',I4)
 110  FORMAT (' ****  (',I3,') *** DO LOOP LEVEL',I2,' TERMINATES WHILE
     1LEVEL',I2,' IS IN EFFECT.     ***')
 120  FORMAT (' ****  (',I3,') *** STATEMENT NUMBER',I6,' DUPLICATES THE
     1 NUMBER AT',I4,'.',8X,'***')
      END
      SUBROUTINE EDIT
C
C     THIS SUBROUTINE EDITS THE DEFINED AND THE REFERENCED STATEMENT
C     NUMBER LIST.
C
C     ON ENTRY, LDEF(I) CONTAINS THE STATEMENT LABELS, IN THE
C     ORDER IN WHICH THEY WERE USED.  THE LABELS OF CONTINUE
C     STATEMENTS WHICH WERE NOT PASSED ON ARE NEGATIVE.
C     LOCDEF(I) CONTAINS THE CARD NUMBER (NREC) OF THE LINE
C     IDENTIFIED BY THAT LABEL.  EXCEPTION FOR DOUBLE BRANCHES--
C     IF LDEF(I)=0, THEN THE STATEMENT WITH THE LABEL LDEF(I-1)
C     WAS A GOTO.  THE TARGET LABEL IS IN LOCDEF(I).
C
C     (1)     DEFINED STATEMENTS THAT ARE NOT REFERENCED ARE DELETED.
C     (2)     THE NEW STATEMENT NUMBERS ARE GENERATED
C     (3)     A STATEMENT NUMBER WHICH IS NEGATIVE IN THE LDEF
C             LIST IS ASSIGNED A NEW STATEMENT NUMBER THE SAME
C             AS THE NEXT POSITIVE LABEL IN THE LDEF LIST
C     (4)     A LABEL FOLLOWED BY A ZERO IN THE LDEF LIST IS
C             ASSIGNED A NEW STATEMENT NUMBER THE SAME AS THE
C             STATEMENT NUMBER ASSIGNED TO THE LABEL GIVEN IN
C             THE LOCREF ARRAY.  (FOR DOUBLE BRANCHES)
C     (5)     PSEUDO-STATEMENT NUMBERS OUTSIDE THE RANGE OF RENUMBERED
C             DEFINED STATEMENT NUMBERS ARE GENERATED FOR EACH
C             REFERENCED STATEMENT WHICH IS NOT DEFINED.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
      IF (NREF.LE.0) NDEF=0
      IF (NDEF.LE.0) RETURN
C
      IF (MDEB.NE.0) THEN
          WRITE (OUTFIL,130) NDEF,NREF
          CALL EDTDMP ('Start EDIT')
      END IF
C
C     SET UP NEWNUM SO THAT IF LDEF(I) NEEDS A NEW NUMBER,
C     NEWNUM(I)=0. IF LDEF(I) WILL REFERENCE LDEF(J), THEN
C     NEWNUM(I)=-LDEF(J).  REMOVE ENTRIES WITH LDEF(I)=0
C
      IT=0
      KFMTS=0
      DO 20 I=1,NDEF
          IF (NNMKLS(I).EQ.5) KFMTS=KFMTS+1
          IF (LDEF(I).GT.0) THEN
C                            POSITIVE IS NORMAL
              IT=IT+1
              NEWNUM(IT)=0
              LDEF(IT)=LDEF(I)
          ELSE IF (LDEF(I).EQ.0) THEN
C                            ZERO MEANS LAST WAS A BRANCH
              NEWNUM(IT)=-LOCDEF(I)
              GO TO 20
          ELSE
C                            NEGATIVE MEANS CONTINUE. LOOK AHEAD
              J=I
 10           J=J+1
              IF (LDEF(J).LT.0.OR.LOCDEF(J).LT.0) GO TO 10
C                            CHECK FOR A FORMAT STATEMENT
              IT=IT+1
              NEWNUM(IT)=-LDEF(J)
              IF (LDEF(J).EQ.0) NEWNUM(IT)=-abs(LDEF(J-1))
              LDEF(IT)=abs(LDEF(I))
          END IF
          LOCDEF(IT)=abs(LOCDEF(I))
 20   CONTINUE
      NDEF=IT
C
      IF (MDEB.NE.0) THEN
          CALL EDTDMP ('Nonzero NEWNUM = refs to be changed')
      END IF
C
C     LDEF NOW CONTAINS DEFINED STATEMENT NUMBERS. LOCDEF(I)
C     HAS LINE NUMBER OF LDEF(I).  NEWNUM(I) HAS ZERO IF LDEF(I)
C     WILL NEED A NEW NUMBER, AND -NNN IF REFERENCES TO LDEF(I)
C     SHOULD BE CHANGED TO REFERENCES TO NNN.
C
C     FOR EACH LREF, SCAN LDEF FOR CHAINS.  BE SURE
C     TARGETS OF GOTOS ARE REFERENCED ALSO.
C
      IT=NREF
      DO 40 I=1,IT
          I1=LREF(I)
C                      GET REFERENCE IN LDEF
          DO 30 J=1,NDEF
              IF (I1.EQ.LDEF(J)) THEN
C                         NEXT LINK IN CHAIN
                  I1=abs(NEWNUM(J))
                  IF (I1.EQ.0) GO TO 40
                  L772=I1
C                      ADD TARGET TO REF LIST
                  CALL RLIST
                  GO TO 40
              END IF
 30       CONTINUE
C                               NOT DEFINED
 40   CONTINUE
C
      IF (MDEB.NE.0) THEN
          CALL EDTDMP ('Chains identified - ready to renumber')
      END IF
C
C     SCAN DEFINED LIST FOR REFERENCES.  DELETE NON-REFERENCED
C     DEFINED STATEMENT NUMBERS.
C
C     FIRST CHECK IF FMTBASE IS LEGAL, IF NOT, IGNORE WITH DIAGNOSTIC.
C       KFMTS COUNTED THE FORMATS, SO CHECK IF HIGHEST NON-FORMAT
C       STATEMENT NUMBER TO BE GENERATED WOULD OVERLAP, OR VICE-VERSA
C       IF FORMATS TO BE NUMBERED BELOW EXECUTABLES (RARE BUT POSSIBLE).
C
      IF (MFMTB.EQ.KB15) THEN
          MFMBAS=0
      ELSE IF (MFMTB.GT.KB15) THEN
          MAXFMT=KB15+(NDEF-KFMTS)*KD15
          IF (MFMTB.LE.MAXFMT) THEN
              MFMBAS=0
              CALL DIAGNO (49)
          ELSE
              MFMBAS=MFMTB
          END IF
      ELSE
          MAXFMT=MFMTB+KFMTS*KD15
          IF (MAXFMT.GT.KB15) THEN
              MFMBAS=0
              CALL DIAGNO (49)
          ELSE
              MFMBAS=MFMTB
          END IF
      END IF
C
      IT=0
      NNUM=0
      NNFM=0
      DO 60 I=1,NDEF
          DO 50 J=1,NREF
              IF (LDEF(I).EQ.LREF(J)) THEN
                  IF (NEWNUM(I).EQ.0) THEN
C                          MAKE NEW NUMBER
                      IF (NNMKLS(I).EQ.5.AND.MFMBAS.GT.0) THEN
                          NNFM=NNFM+1
                          NEWNUM(I)=KD15*NNFM+MFMBAS
                      ELSE
                          NNUM=NNUM+1
                          NEWNUM(I)=KD15*NNUM+KB15
                      END IF
                  END IF
                  IT=IT+1
                  LDEF(IT)=LDEF(I)
                  NEWNUM(IT)=NEWNUM(I)
                  LOCDEF(IT)=LOCDEF(I)
                  NNMKLS(IT)=NNMKLS(I)
                  GO TO 60
              END IF
 50       CONTINUE
C                            NOT REFERENCED
 60   CONTINUE
      NDEF=IT
C
      IF (MDEB.NE.0) THEN
          CALL EDTDMP ('New statement numbers:')
      END IF
C
C     SCAN LDEF FOR INDIRECT REFERENCES AND REPLACE THEM
C
      IT=0
      DO 100 I=1,NDEF
          DO 80 IC=1,10
              IF (NEWNUM(I).GT.0) GO TO 100
              I1=abs(NEWNUM(I))
              DO 70 J=1,NDEF
                  IF (LDEF(J).EQ.I1) THEN
                      NEWNUM(I)=NEWNUM(J)
                      GO TO 80
                  END IF
 70           CONTINUE
              CALL DIAGNO (46)
 80       CONTINUE
C                            LOOP OF GOTO-S. BREAK IT
          IF (IT.NE.0) GO TO 90
          IT=1
          CALL PAGE (-20)
          CALL PAGE (1)
          WRITE (OUTFIL,170)
          WRITE (OUTFIL,160)
 90       NNUM=NNUM+1
          NEWNUM(I)=KD15*NNUM+KB15
          NMSG=NMSG+1
          CALL PAGE (1)
          WRITE (OUTFIL,140) NMSG,I1,NEWNUM(I)
 100  CONTINUE
C
C     SCAN REFERENCED STATEMENT LIST FOR MISSING DEFINITIONS.
C
      IT=0
      DO 120 I=1,NREF
          DO 110 J=1,NDEF
              IF (LREF(I).EQ.LDEF(J)) GO TO 120
 110      CONTINUE
C
C     ADD PSEUDO-STATEMENT NUMBER.
C
          LERR=2
          IF (IT.LE.0) THEN
              IT=1
              CALL PAGE (-20)
              CALL PAGE (4)
              WRITE (OUTFIL,150)
              WRITE (OUTFIL,160)
          END IF
          NDEF=NDEF+1
          IF (NDEF.GT.MXPLBL) THEN
              CALL DIAGNO (6)
              NDEF=-1
              MP2=0
              RETURN
          END IF
          LDEF(NDEF)=LREF(I)
          LOCDEF(NDEF)=0
          NEWNUM(NDEF)=NDEF*KD15+KB15
          NMSG=NMSG+1
          CALL PAGE (1)
          WRITE (OUTFIL,140) NMSG,LREF(I),NEWNUM(NDEF)
 120  CONTINUE
      RETURN
C
C
C
 130  FORMAT (' FOLLOWING *DEBUG OUTPUT FROM SUBR EDIT'/' NDEF = ',I7,'
     1 NREF = ',I7)
 140  FORMAT (7X,'(',I3,') *** STATEMENT NUMBER',I7,' IS ASSIGNED NUMBER
     1',I7,'.',13X,'***')
 150  FORMAT (' ',12X,'*** THE FOLLOWING REFERENCED STATEMENTS ARE NOT D
     1EFINED')
 160  FORMAT (13X,'*** PSEUDO-STATEMENT NUMBERS HAVE BEEN ASSIGNED.'/' '
     1)
 170  FORMAT (' ',12X,'*** THE FOLLOWING STATEMENTS ARE IN ENDLESS CHAIN
     1S OF GOTO''S.')
      END
      SUBROUTINE EDTDMP (MSG)
      CHARACTER*(*) MSG
C
C     DUMP TABLES FOR DEBUG OF SUBROUTINE EDIT.
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C
      WRITE (OUTFIL,10) MSG
      WRITE (OUTFIL,20) (LDEF(I),I=1,NDEF)
      WRITE (OUTFIL,30) (NEWNUM(I),I=1,NDEF)
      WRITE (OUTFIL,40) (LOCDEF(I),I=1,NDEF)
      WRITE (OUTFIL,50) (LREF(I),I=1,NREF)
      WRITE (OUTFIL,60) (NNMKLS(I),I=1,NREF)
      RETURN
C
C
 10   FORMAT (' ',A)
 20   FORMAT (' LDEF  ',9I7)
 30   FORMAT (' NEWNUM',9I7)
 40   FORMAT (' LOCDEF',9I7)
 50   FORMAT (' LREF  ',9I7)
 60   FORMAT (' KLASS ',9I7)
      END
      SUBROUTINE HEADER
C
C                  THIS ROUTINE CENTERS JOB HEADINGS
C
      INCLUDE  'tidy.inc'
      CHARACTER*2 KUPPER
      IF (IPASS.EQ.1) THEN
          DO 10 I=1,72
              JOB(I)=JINT(I)
 10       CONTINUE
      ELSE
C
          DO 20 I=1,80
              JOB(I)=IOUT(I)
 20       CONTINUE
C
          IF (MSER.LT.0) THEN
C
C     SET UP COLUMNS 73-75 BASED ON *LABE OPTION
              IF (MLBL.EQ.0) THEN
C     USE *ROUT VALUE
                  I=(NROUT-1)/26
                  J=NROUT-I*26
                  IF (I.EQ.0) THEN
                      KOL73(3)=KBL
                      KOL73(2)=KABC(J)
                  ELSE
                      KOL73(2)=KABC(I)
                      KOL73(3)=KABC(J)
                  END IF
C
                  KOL73(1)=KBL
              ELSE
C
C     COPY PROGRAM/SUBROUTINE/FUNCTION CARD SERIAL INFORMATION
                  DO 30 I=1,3
                      KOL73(I)=KUPPER(SERIAL(I))
 30               CONTINUE
              END IF
          END IF
      END IF
C
      DO 40 I=73,80
          JOB(I)=KBL
 40   CONTINUE
C
C          COMPRESS STATEMENT BY ELIMINATING MULTIPLE BLANKS
C
      J=1
      K=0
      DO 50 I=1,80
          IF (JOB(I).EQ.KBL) THEN
              IF (K.EQ.1) GO TO 50
              K=1
          ELSE
              K=0
          END IF
          JOB(J)=JOB(I)
          J=J+1
 50   CONTINUE
      DO 60 I=J,80
          JOB(I)=KBL
 60   CONTINUE
C
C                           CENTER HEADING
C
      IB=(80-J)/2
 70   I=J+IB
      JOB(I)=JOB(J)
      J=J-1
      IF (J.GT.0) GO TO 70
C
C                   ELIMINATE REMAINING NON-BLANKS
C
      IB=I-1
      DO 80 I=1,IB
          JOB(I)=KBL
 80   CONTINUE
      RETURN
      END
      SUBROUTINE HOLSCN (LTYPE,LSSCN,LNSTR)
C     THIS SUBROUTINE SCANS ALL FORTRAN CARDS FOR FIELDS OF HOLLERITH-
C     TYPE CONSTANTS.  IN THESE FIELDS,
C     CHARACTERS ARE REPLACED WITH EQUIVALENT CHARACTERS WHICH WILL NOT
C     BE TREATED BY ANALYSIS ROUTINES.
C     THE SEARCH IS MADE BY CHECKING FOR PATTERNS -SNNNL-, WHERE S IS A
C     SPECIAL CHARACTER, NNN IS A DECIMAL NUMBER, AND L IS THE LETTER H,
C     L, OR R.  IN ADDITION, FOR FORMAT STATEMENTS ONLY, IT ACCEPTS THE
C     PATTERN SNNNXNNNL, THE RESULT OF A MISSING -,- AFTER X.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
      CHARACTER*2 IT,KPARAM,KUPPER,KCTRAN
      LOGICAL LHTRN,ISDEL,STRP73
C
      STRP73=.TRUE.
      JCOL=6
      LNSTR=0
      LNTMP=0
      NLHTRN=0
C     IF FORMAT STATEMENT, SKIP FIRST 7 NON-BLANK CHARACTERS
      IF (LTYPE.EQ.26) THEN
           DO 20 N=1,7
10              JCOL=JCOL+1
                IF (JINT(JCOL).EQ.KBL) GO TO 10
                IF (MCASE.EQ.0) JINT(JCOL)=KCTRAN(JINT(JCOL))
20         CONTINUE
           GO TO 130
      END IF
C
C                  *****************************************
C                  *                                       *
C                  *    PROCESS NON-FORMAT STATEMENTS.     *
C                  *                                       *
C                  *****************************************
C
      LFIR=6
      IFIR=14
C                            SET FLAG FOR NON-FORMAT
      IGOOF=-1
C                   LOOK FOR SPECIAL CHARACTERS.
30    I=JCOL
      DO 60 JCOL=I,JMAX
           IT=JINT(JCOL)
           ISDEL=.FALSE.
C          (CHECK FOR SPL CHAR BEFORE DELIMS SINCE NEED J TO SET IFIR.)
C
C     =    ,    (    /    )    +    -    *    .    $    -    '    & NONE
C     1    2    3    4    5    6    7    8    9    10   11   12   13  14
C
           DO 50 J=1,13
                IF (IT.EQ.KSPK(J)) THEN
C                   FOUND ONE.  IS IT THE FIRST...
                     IF (IFIR.EQ.14) THEN
C                   YES
                          IFIR=J
                          LFIR=JCOL
C     QUIT IF THIS STATEMENT TYPE DOESN'T ALLOW STRINGS.  JUST NEEDED
C     IFIR AND LFIR POINTERS.
                          IF (LSSCN.EQ.0.AND.LTYPE.NE.0)
     1                     THEN
                               IF (MCASE.EQ.0) THEN
                                    DO 40 I=JCOL,JMAX
                                         JINT(I)=KCTRAN(JINT(I))
40                                  CONTINUE
                                ENDIF
                               IF (MDEB.GT.0) WRITE (OUTFIL,320) IFIR,
     1                          LFIR
                               RETURN
                          END IF
                     END IF
                     ISDEL=IT.EQ.KDEL1.OR.(IT.EQ.KDEL2.AND.LTYPE.NE.13)
                     IF (ISDEL) GO TO 180
                     GO TO 70
                END IF
50         CONTINUE
C     (DELIMS MAY NOT BE SPECIAL CHARACTER, CHECK TO BE SURE)
           ISDEL=IT.EQ.KDEL1.OR.(IT.EQ.KDEL2.AND.LTYPE.NE.13)
           IF (ISDEL) GO TO 180
           IF (MCASE.EQ.0) JINT(JCOL)=KCTRAN(IT)
60    CONTINUE
      GO TO 310
C                   LOOK FOR FOLLOWING NUMBER.
70    IF (JCOL.EQ.JMAX) GO TO 310
      JCOL=JCOL+1
      CALL RSTAT
C                   REPEAT IF NO NUMBER.
      IF (INT(L772).EQ.0) GO TO 30
C     MAKE IT UPPER CASE
      IF (MCASE.EQ.0) JINT(JCOL)=KCTRAN(JINT(JCOL))
      IT=KUPPER(JINT(JCOL))
C                  IS IT -H-,-L-, OR -R-
      IF (IT.EQ.KABC(8)) THEN
           LHTRN=MOD(KHTRAN,2).EQ.0
      ELSE IF (IT.EQ.KABC(12).OR.IT.EQ.KABC(18)) THEN
           LHTRN=KHTRAN.LT.2
C     COMPLAIN ABOUT L OR R IF ANSI FLAG SET.
           IF (MANSI.EQ.0) CALL DIAGNO (36)
      ELSE
           GO TO 30
      END IF
C                  MARK AS PART OF STRING (FOR INDENTING)
      IF (LHTRN) JINT(JCOL)(2:2)=KAT(2:2)
C
C     ALSO MARK THE NUMBERS.
      KTMP=INT(L772)
      I=JCOL
80    I=I-1
      IF (JINT(I).EQ.KBL) GO TO 80
      IF (LHTRN) JINT(I)(2:2)=KAT(2:2)
      KTMP=KTMP/10
      IF (KTMP.GT.0) GO TO 80
      IP=I
C                  FIND LIMITS OF HOLLERITH FIELD.
      I=JCOL+1
      JCOL=JCOL+INT(L772)
C                   L772 IS THE LENGTH OF THE FIELD, AS FOUND BY RSTAT
C                  CHECK FOR CASE OF HOLLERITH BLANKS SPILLING OFF
C                  END OF CARD. E.G. I=6HXXXXX
      IF (JCOL.GT.JMAX) THEN
C                  REPLACE CURRENT END CARD MARK.
           JINT(JMAX+1)=KBL
C                   AND SET NEW ONE
           JMAX=JCOL
           JINT(JMAX+1)=KERM
      END IF
C                  CHANGE ALL CHARACTERS IN HOLLERITH FIELD.
90    DO 100 J=I,JCOL
           JJ=J
C                  SKIP CHARACTERS IN 73-80 OF ORIGINAL CARD.
95         IF (JINT(JJ)(2:2).EQ.'!') THEN
                JINT(JJ)(2:2)=' '
                JJ=JJ+1
                IF (STRP73) THEN
                     STRP73=.FALSE.
                     CALL DIAGNO(52)
                END IF
                GO TO 95
           ELSE
                JINT(JJ)(2:2)=KAT(2:2)
           END IF
100   CONTINUE
      JCOL=JJ
      IF (.NOT.LHTRN) THEN
C
C     TURN THIS ON IF WANT LOGGING OF H TRANSLATIONS IN FORMATS
           IF (KHLOG.EQ.0) NLHTRN=NLHTRN+1
C
C     IF TRANSLATING H-FIELDS, COPY STRING AND DUPLICATE APOSTROPHES.
           LNTMP=max(int(L772),LNTMP)
           JINT(IP)=KAPSTR
           IP=IP+1
           J=I
C         SKIP ANY EXTRA CHARS FROM 73-80
110        IF(JINT(J)(2:2).EQ.'!') THEN
                JCOL=JCOL+1
                J=J+1
                GO TO 110
           END IF
C
           JINT(IP)=JINT(J)
           IF (JINT(J).EQ.KAPSTR) THEN
                IP=IP+1
                IF (IP.GE.J) CALL MOVSTR (J)
                JINT(IP)=KAPSTR
           END IF
           J=J+1
           IP=IP+1
           IF (J.LE.JCOL) GO TO 110
           JINT(IP)=KAPSTR
120        IP=IP+1
           IF (IP.LE.JCOL) THEN
                JINT(IP)=KBL
                GO TO 120
           END IF
      END IF
      GO TO 30
C
C                  **********************************
C                  *                                *
C                  *   PROCESS FORMAT STATEMENTS.   *
C                  *                                *
C                  **********************************
C
130   IGOOF=0
      IFIR=3
      LFIR=JCOL
      GO TO 170
C
C                  LOOK FOR SPECIAL CHARACTER
140   IF (JCOL.GT.JMAX) GO TO 310
      I=JCOL
      DO 160 JCOL=I,JMAX
           IT=JINT(JCOL)
           ISDEL=IT.EQ.KDEL1.OR.IT.EQ.KDEL2
           IF (ISDEL) GO TO 180
           DO 150 J=1,12
                IF (IT.EQ.KSPK(J)) GO TO 220
150        CONTINUE
           IF (MCASE.EQ.0) JINT(JCOL)=KCTRAN(IT)
160   CONTINUE
      GO TO 310
C
C                  SKIP IF NOT * OR '
170   IF (JINT(JCOL).NE.KDEL1.AND.JINT(JCOL).NE.KDEL2) GO TO 220
C                  CHANGE ALL CHARACTERS BETWEEN *S OR 'S
180   KPARAM=JINT(JCOL)
C                  MARK AS PART OF STRING (FOR INDENTING)
      JINT(JCOL)(2:2)=KAT(2:2)
      IP=JCOL
C
190   IF (JCOL.EQ.JMAX) GO TO 310
      JCOL=JCOL+1
      IT=JINT(JCOL)
C
C     IF CHAR WAS IN COLS 73-80, UNMARK IT SO COPY WILL EAT IT.
      IF (IT(2:2).EQ.'~') THEN
           JINT(JCOL)(2:2)=' '
           IF (STRP73) THEN
                STRP73=.FALSE.
                CALL DIAGNO(52)
           END IF
      ELSE
           JINT(JCOL)(2:2)=KAT(2:2)
      END IF
      IF (IT.EQ.KPARAM) THEN
           IF (JINT(JCOL+1).NE.KPARAM) GO TO 200
C     THIS IS A LITERAL -- NOT TERMINAL DELIMITER
           JCOL=JCOL+1
           JINT(JCOL)(2:2)=KAT(2:2)
      END IF
      GO TO 190
C                            ALL CHANGED, CHANGE DELIMS IF DESIRED.
200   IF (KDTRAN.EQ.1.AND.KPARAM.NE.KDEL1) THEN
           JINT(IP)=KAPSTR
           JINT(JCOL)=KAPSTR
           J=IP
210        J=J+1
           IF (J.LT.JCOL) THEN
                IF (JINT(J).EQ.KAPSTR) THEN
C     DUPLICATE LITERAL VERSION OF DELIMITER
                     CALL MOVSTR (J)
                     JINT(J)=KAPSTR
                END IF
                GO TO 210
           END IF
      END IF
      IF (IGOOF.EQ.-1) GO TO 70
C                  LOOK FOR FOLLOWING NUMBER
220   IF (JCOL.EQ.JMAX) GO TO 310
      JCOL=JCOL+1
      CALL RSTAT
C                  IF NOT A NUMBER, START AGAIN
      IF (INT(L772).EQ.0) GO TO 140
C                  NUMBER FOUND. LOOK AT NEXT CHARACTER.
      IF (MCASE.EQ.0) JINT(JCOL)=KCTRAN(JINT(JCOL))
      IT=KUPPER(JINT(JCOL))
C                  IS IT -H-
      IF (IT.EQ.KABC(8)) THEN
           LHTRN=MOD(KHTRAN,2).EQ.0
           GO TO 250
C                  MAYBE L OR R
      ELSE IF (IT.EQ.KABC(12).OR.IT.EQ.KABC(18)) THEN
           LHTRN=KHTRAN.LT.2
           IF (MANSI.EQ.0) CALL DIAGNO (36)
           GO TO 250
      END IF
C                  IF NOT -X-, START AGAIN.
      IF (IT.NE.KABC(24)) GO TO 140
C                  X FOUND.  LOOK AT NEXT.
230   IF (JCOL.EQ.JMAX) GO TO 310
      JCOL=JCOL+1
      IF (JINT(JCOL).EQ.KBL) GO TO 230
      IF (MCASE.EQ.0) JINT(JCOL)=KCTRAN(JINT(JCOL))
      IT=KUPPER(JINT(JCOL))
C                  IS IT -*-
      IF (IT.EQ.KDEL1.OR.IT.EQ.KDEL2) GO TO 170
C                  IS IT -)- OR -,-
      IF (IT.EQ.KSPK(2)) GO TO 220
      IF (IT.EQ.KSPK(5)) GO TO 220
C
C     INSERT A COMMA
      DO 240 J=JMAX,JCOL,-1
           JINT(J+1)=JINT(J)
240   CONTINUE
      JINT(JCOL)=KSPK(2)
      JMAX=JMAX+1
      JINT(JMAX+1)=KERM
      CALL DIAGNO (25)
      IGOOF=1
      GO TO 220
C
C                  HOLLERITH FOUND.   FIND LIMITS OF FIELD.
250   IF (LHTRN) JINT(JCOL)(2:2)=KAT(2:2)
C
C     ALSO MARK THE NUMBERS.
      J=INT(L772)
      I=JCOL
260   I=I-1
      IF (JINT(I).EQ.KBL) GO TO 260
      IF (LHTRN) JINT(I)(2:2)=KAT(2:2)
      J=J/10
      IF (J.GT.0) GO TO 260
C
      IP=I
      I=JCOL+1
      JCOL=JCOL+INT(L772)
      IF (JCOL.LE.JMAX) GO TO 270
      JINT(JMAX+1)=KBL
      JMAX=JCOL
      JINT(JMAX+1)=KERM
270   DO 280 J=I,JCOL
           JJ=J
275        IF (JINT(JJ)(2:2).EQ.'~') THEN
                JINT(JJ)(2:2)=' '
                JJ=JJ+1
           IF (STRP73) THEN
                STRP73=.FALSE.
                CALL DIAGNO(52)
           END IF
                GO TO 275
           END IF
           JINT(JJ)(2:2)=KAT(2:2)
280   CONTINUE
      IF (.NOT.LHTRN) THEN
C
C     IF TRANSLATING H-FIELDS, COPY STRING AND DUPLICATE APOSTROPHES.
           IF (KHLOG.EQ.0) NLHTRN=NLHTRN+1
           JINT(IP)=KAPSTR
           IP=IP+1
           J=I
C          SKIP ANYTHING IN 73-80 IF MARKED.
290        IF (JINT(J)(2:2).EQ.'~') THEN
                 JCOL=JCOL+1
                 J=J+1
                 GO TO 290
           END IF
           JINT(IP)=JINT(J)
           IF (JINT(J).EQ.KAPSTR) THEN
                IP=IP+1
                IF (IP.GE.J) CALL MOVSTR (J)
                JINT(IP)=KAPSTR
           END IF
           J=J+1
           IP=IP+1
           IF (J.LE.JCOL) GO TO 290
           JINT(IP)=KAPSTR
300        IP=IP+1
           IF (IP.LE.JCOL) THEN
                JINT(IP)=KBL
                GO TO 300
           END IF
      END IF
      GO TO 220
C
310   IF (LNTMP.GT.0) LNSTR=LNTMP
      IF (NLHTRN.GT.0) THEN
           IF (LTYPE.NE.26) CALL DIAGNO (39)
           NLHTRN=0
      END IF
      IF (MDEB.GT.0) WRITE (OUTFIL,320) IFIR,LFIR
      RETURN
 320  FORMAT (' HOLSCN: IFIR = ',I2,' AT COL ',I4)
      END
      SUBROUTINE JTYP19 (JRTCOD,NBCOLD)
C
C                  ***** JTYPE = 19
C     FORMAT (
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C
C     ERROR IF NO STATEMENT NUMBER OR FIRST SPECIAL CHAR NOT (
      IF (L15.EQ.0.OR.JINT(JMAX).NE.KSPK(5)) THEN
           JRTCOD=1
           RETURN
      END IF
C
      IF (MEX.EQ.0) THEN
           IF (MCOL.EQ.-1) THEN
C
C          IF COLLECTING FORMATS, START THEM IN COLUMN 7 (OR JUST).
                ICOL=6
                IF (JUST.GT.0) ICOL=JUST-1
           END IF
C
           CALL KWDWRT(26)
C                            COPY REST OF CARD
           IF (MCOL.EQ.0) THEN
                JRTCOD=3
                RETURN
           END IF
C                            ONTO UNIT 2
           CALL COPY (0)
           IMAX=ICOL
           JTYPE=NREC
           CALL WRBIN (SCFIL2,KILI,SERIAL,IOUT)
           NRT2=NRT2+1
           NBLC=NBCOLD
      ELSE
C
C     EXEMPT FLAG IS ON - TRANSFER TO TAPE1 OR TAPE2 WITHOUT REMOVING
C     ANY BLANKS.
C
C          ITYPE=NREC
           CALL DLIST (MERR)
           NBLC=NBCOLD
           ICOLSV=6
           IF (MERR.GE.0) THEN
                IF (MCOL.NE.0) THEN
                     CALL WRBIN (SCFIL2,KILI,SERIAL,JINT)
                     NRT2=NRT2+1
                ELSE
                     CALL WRBIN (SCFIL1,KILI,SERIAL,JINT)
                     NRT1=NRT1+1
                END IF
           END IF
      END IF
C
      JRTCOD=2
      RETURN
      END
      SUBROUTINE JTYP31(JRTCOD)
C
C                  ***** JTYPE = 31
C     IF (ARITHMETIC) 1,2,3   OR   IF (LOGICAL) STATEMENT.
C
      INCLUDE 'tidy.inc'
      CHARACTER*2 JT
      COMMON /PS1SUB/ KSTC(5), NIFBLK
C
      CALL KWDWRT (38)
C                  COPY UNTIL CLOSED PARENTHESES
      CALL COPY (-1)
      IF (MEOF.GE.0) GO TO 80
      ICOL=ICOL+1
      CALL RSTAT
      IF (INT(L772).NE.0) THEN
C
C     STATEMENT IS    IF (ARITHMETIC) 1,2,3
C
           NCOM=0
           MILDO=-1
           CALL DLIST (MERR)
           IF (MERR.LT.0) GO TO 80
10         CALL ADNUM (JERR)
           IF (JERR.EQ.1) THEN
                JRTCOD=2
                RETURN
           END IF
           CALL COPY (1)
           IF (LCPY.EQ.KSPK(2)) THEN
                NCOM=NCOM+1
                IF (NCOM.GT.3) GO TO 80
                IF (NCOM.EQ.3) CALL DIAGNO (18)
                CALL RSTAT
                IF (INT(L772).EQ.0) GO TO 80
                GO TO 10
           END IF
           IF (LCPY.NE.KERM) GO TO 80
           IF (NCOM.LE.0) GO TO 80
           IF (NCOM.EQ.1) CALL DIAGNO (18)
           MTRAN=MLGC
           JRTCOD=3
           RETURN
      END IF
C
C     STATEMENT IS   IF (LOGICAL) STATEMENT
C
      MLGC=0
C
C        CHECK FOR 'IF () THEN' UNLESS IT IS  ELSEIF () THEN
      IF (JTYPE.EQ.43) GO TO 40
      I=69
      CALL KWSCAN (I,KSTC)
      IF (I.NE.69) GO TO 40
      CALL KWDWRT(I)
C        LOOP TO CHECK REST FOR BLANKS.
      DO 20 I=JCOL,JMAX
           IF (JINT(I).EQ.KERM) GO TO 30
           IF (JINT(I).NE.KBL) GO TO 40
20    CONTINUE
30    NIFBLK=NIFBLK+1
      JRTCOD=4
      RETURN
C
C                   LOOK FOR FIRST SPECIAL CHARACTER.
40    DO 60 LFIR=JCOL,JMAX
           JT=JINT(LFIR)
           DO 50 IFIR=1,11
                IF (JT.EQ.KSPK(IFIR)) GO TO 70
50         CONTINUE
60    CONTINUE
      LFIR=6
      IFIR=14
70    JRTCOD=5
      RETURN
C
80    JRTCOD=1
      RETURN
C
      END
      SUBROUTINE JTYP33 (ITYPE,JRTCOD)
C
C     PROCESS TYPE 33 CARDS - AGS 23 DEC 1993
C
C     JRTCOD IS RETURN CODE - USE COMPUTED GOTO TO BRANCH TO PROPER
C      PLACE IN PASS1.
C
      INCLUDE 'tidy.inc'
C
C                  ***** JTYPE = 33
C     PRINT, TYPE, WRITE, PUNCH, READ, ACCEPT.
C
      CALL KWDWRT (ITYPE)
      CALL RSTAT
      IF (INT(L772).NE.0) GO TO 20
C
C     HAVE WRITE  FMT,LIST
C
C            , AS IN PRINT IFT,XXX
      IF (IFIR.NE.2) THEN
C            *, AS IN PRINT *,XXX
           IF (IFIR.EQ.8.OR.IFIR.EQ.12.OR.IFIR.EQ.14) THEN
                JRTCOD=1
           ELSE
                JRTCOD=2
           END IF
           RETURN
      END IF
C
   10 CALL COPY (1)
      IF (LCPY.EQ.KSPK(2)) THEN
           JRTCOD=3
           RETURN
      END IF
      IF (MEOF.LT.0) GO TO 10
      JRTCOD=2
      RETURN
C
C     HAVE WRITE  12345 LIST
C
   20 CALL ADNUM (JRTCOD)
      IF (JRTCOD.EQ.1) THEN
           JRTCOD=4
           RETURN
      END IF
      IF (IFIR.EQ.2) GO TO 10
      IF (JMAX.GT.JCOL) THEN
           JRTCOD=2
      ELSE
           IMAX=ICOL
           JRTCOD=5
      END IF
      RETURN
      END
      CHARACTER*2 FUNCTION KCTRAN(C)
C
C     CONVERTS ALL LETTERS TO A SINGLE CASE, SELECTED BY USER'S CALL TO
C      SUBROUTINE KCTSET.
C     PORTABLE VERSION - NOT ASCII/EBCDIC DEPENDENT.
C     AGS 12 OCT 93
C
C
      CHARACTER CT
      CHARACTER*2 C
C     COMMON BLOCK FOR CHARACTER TRANSLATION TABLES
      COMMON /CTRAN/ LININ,LINOUT
      CHARACTER*26 LININ,LINOUT
      SAVE
C
C     FIND POSITION OF CHARACTER IN INPUT-CASE ALPHABET
      CT=C(1:1)
      J=INDEX(LININ,CT)
C
C     IF FOUND, RETURN OUTPUT-CASE EQUIVALENT, OTHERWISE ORIGINAL CHAR.
      IF (J.GT.0) THEN
           KCTRAN=LINOUT(J:J)
      ELSE
           KCTRAN=C
      END IF
C
      RETURN
      END
      SUBROUTINE KCTSET (IP)
C
C     SET CHARACTER TRANSLATION TABLE FOR KCTRAN:
C     IP = 0 - LOWER TO UPPER
C     IP = 1 - UPPER TO LOWER
C
C     COMMON BLOCK FOR CHARACTER TRANSLATION TABLES
      COMMON /CTRAN/ LININ,LINOUT
      CHARACTER*26 LININ,LINOUT
      CHARACTER*26 CTBL(0:1)
      SAVE
      DATA CTBL/'abcdefghijklmnopqrstuvwxyz',
     1          'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C
C     ASSIGN INPUT AND OUTPUT ALPHABETS BASED ON VALUE OF IP.
      LININ=CTBL(IP)
      LINOUT=CTBL(1-IP)
C
      RETURN
      END
      CHARACTER*2 FUNCTION KHIDE (C)
      CHARACTER*2 C
C
C     CONVERT CHARACTERS IN HOLSCN STRINGS TO SPECIAL FORM
C      (UNLESS ALREADY SET TO INDICATE EMBEDDED COMMENT STATEMENT)
C      SO THAT BLANKS WILL NOT BE REMOVED FROM STRINGS.
C
      IF (C(2:2).EQ.' ') THEN
           KHIDE=C(1:1)//'@'
      ELSE
           KHIDE=C
      END IF
      RETURN
      END
      SUBROUTINE KIMPAK
C
C     THIS ROUTINE PACKS SUPER-CARD IMAGES FROM IOUT(I) INTO KIM(I,J).
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
      LOGICAL CONIND,SPLSTR,SAVBLK
C
      CONIND=.TRUE.
      SPLSTR=.FALSE.
C
C     SET BLANK STRIP MODE
      SAVBLK=(MEX.GT.0 .OR. (MEX.LT.0.AND.(KLASS.EQ.3.OR.KLASS.EQ.5)))
C
      IF (MDEB.GT.0) WRITE (OUTFIL,130) (IOUT(I),I=1,IMAX)
 10   J=0
C
 20   J=J+1
C
C     KLASS 2 AND BELOW HAVE NO CONTINUATIONS.
      IF (KLASS.LT.2) THEN
           K7=0
           JL=1
           IF (KLASS.EQ.1 .AND. MSER.EQ.0) THEN
                JR=80
           ELSE
                JR=72
           END IF
           GO TO 90
      END IF
C
C     INDENTING COULD MAKE CARD OVERFLOW CONTINUATIONS, IF SO, REPACK.
      IF (J.GT.MAXCNT+1) THEN
C         ISSUE DIAGNOSTIC IF WE HAVE ALREADY TRIED TO REPACK IT.
           IF (.NOT.CONIND) THEN
                CALL DIAGNO (37)
                J=MAXCNT+1
                GO TO 120
           END IF
C         FIRST TIME - SET FLAG AND TRY AGAIN WITHOUT INDENTING.
           CONIND=.FALSE.
           JL=7
           JR=72
           K7=ICOLSV
           GO TO 10
      END IF
C
C     PREPARE COLUMNS 1-6 OF FIRST CARD.
      IF (CONIND) THEN
           IF (J.EQ.1) THEN
                K7=ICOLSV
                DO 30 I=1,6
                     KIM(I,1)=IOUT(I)
 30             CONTINUE
           ELSE
C     BLANK COLUMN 1-5
                DO 40 I=1,5
                     KIM(I,J)=KBL
 40             CONTINUE
C     COLUMN 6 - NUMBER SERIALLY UNLESS CCHR SET OTHERWISE.
                IF (KCTCTL.EQ.0) THEN
                     IF (J.LT.11) THEN
                          KIM(6,J)=KDIG(J)
                     ELSE
                          KIM(6,J)=KSPK(10)
                     END IF
                ELSE
                     KIM(6,J)=KCTCHR
                END IF
           END IF
C
C     SET LEFT EDGE OF TEXT
C      (USE COL 7 IF EXEMPT, NON-INDENTED, OR IF PART OF STRING
           IF (SAVBLK.OR.ICOLSV.EQ.6.OR.(IOUT(K7)(2:2).EQ.KAT(2:2).
     1      AND.IOUT(K7+1)(2:2).EQ.KAT(2:2))) THEN
                JL=7
           ELSE
                JL=ICOLSV
                IF (J.GT.1) JL=JL+1
                DO 50 I=7,JL
                     KIM(I,J)=KBL
 50             CONTINUE
                JL=JL+1
           END IF
C
C     SET RIGHT EDGE OF TEXT
C     FIRST GET RIGHT-MOST POTENTIAL CHAR IN STRING (KRR)
           JR=72
           KRR=K7+JR-JL+1
           IF (KRR.GT.IMAX) THEN
C     IF PAST END OF STATEMENT, STOP AT END.
                JR=JL+IMAX-K7-1
                GO TO 90
           END IF
C
C     NOW CHECK IF WE CAN BREAK IT HERE.
C     BREAK IF PART OF A STRING. KIMPAK PROTECTS DELIMETERS ALSO.
 60        IF (IOUT(KRR)(2:2).EQ.KAT(2:2)) THEN
C
C     FORMAT STATEMENTS - MAY HAVE PROBLEMS WITH QUOTES AT END.
                IF (KLASS.EQ.5) THEN
C          DON'T SPLIT IF TURNED OFF OR AT TOP INDENT LEVEL.
                     IF (KFSPL.EQ.1.OR.ICOLSV.EQ.6) GO TO 90
C          IF NEXT CHAR NOT IN STRING, BREAK IS FINE.
                     IF (IOUT(KRR+1)(2:2).NE.KAT(2:2)) GO TO 90
C
C          COLUMN 72 NOT A QUOTE, CAN SPLIT ON COL 71
                     IF (IOUT(KRR).NE.KAPSTR) THEN
C          INSERT ',' IN STRING
                          JR=JR-1
                          SPLSTR=.TRUE.
                     ELSE
C          COLUMN 72 QUOTE WITHIN A STRING, BACKTRACK.
                          KRR=KRR-1
                          JR=JR-1
                          IF (JR.GT.JL) GO TO 60
                     END IF
C     END FORMAT STRING BREAKER
                END IF
                GO TO 90
           END IF
C
C     BREAK IF IT IS A BLANK (NOT IN STRING)
           IF (IOUT(KRR).EQ.KBL) GO TO 90
C
C     GO BACK IF LEFT PARENTHESIS
 70        IF (IOUT(KRR).EQ.KSPK(3)) THEN
                KRR=KRR-1
                JR=JR-1
                GO TO 70
           END IF
C
C     BREAK FOR SPECIAL CHARACTERS (EXCEPT DECIMAL POINTS)
           DO 80 I=1,14
                IF (IOUT(KRR).EQ.KSPK(I).AND.I.NE.9) GO TO 90
 80        CONTINUE
C
C     OTHERWISE BACK UP ONE, TRY AGAIN.
           KRR=KRR-1
           JR=JR-1
           IF (JR.GT.JL) GO TO 60
C
C     IF GO ALL THE WAY BACK, FORCE IT TO 72
           JR=72
      END IF
C
C     COPY THE TEXT
 90   DO 100 I=JL,JR
           K7=K7+1
           IF (K7.LE.IMAX) THEN
                KIM(I,J)=IOUT(K7)
           ELSE
                KIM(I,J)=KBL
           END IF
 100  CONTINUE
C
C     STRING SPLITTER
      IF (SPLSTR) THEN
           KIM(JR+1,J)=KAPSTR
           IOUT(K7-1)=KSPK(2)
           IOUT(K7)=KAPSTR
           K7=K7-2
           JR=JR+1
           SPLSTR=.FALSE.
      END IF
C
C     SCRUB GARBAGE OFF END IF SHORTER THAN 72
c     IF (JR.LT.72) THEN
           DO 110 I=JR+1,80
                KIM(I,J)=KBL
 110       CONTINUE
c     END IF
C
C     DO ANOTHER CONTINUATION IF NECESSARY.
      IF (K7.LT.IMAX) GO TO 20
C
 120  NCD=J
      IF (MDEB.GT.0) WRITE (OUTFIL,140) ((KIM(I,J),I=1,80),J=1,NCD)
      RETURN
C
 130  FORMAT (' KIMPAK - IOUT IN:'/(1X,75A1))
 140  FORMAT (' KIMPAK - OUTPUT:'/(1X,80A1))
      END
      CHARACTER*2 FUNCTION KUPPER(C)
C
C     CONVERTS LOWER-CASE LETTERS TO UPPER-CASE. PORTABLE VERSION.
C     AGS 23 APR 93
C
      CHARACTER CT
      CHARACTER*2 C
      CHARACTER*26 LC,UC
      SAVE
      DATA LC/'abcdefghijklmnopqrstuvwxyz'/
      DATA UC/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C
C     FIND POSITION OF CHARACTER IN LOWER-CASE ALPHABET
      CT=C(1:1)
      J=INDEX(LC,CT)
C
C     IF FOUND, RETURN UPPER-CASE EQUIVALENT, OTHERWISE ORIGINAL CHAR.
      IF (J.GT.0) THEN
           KUPPER=UC(J:J)
      ELSE
           KUPPER=C
      END IF
C
      RETURN
      END
      SUBROUTINE KWDWRT (KW)
C
C     PRINTS FORTRAN KEYWORDS (KW > 0)
C     SETS CASE OF KEYWORDS (KW < 0) : -1 = UPPER CASE, -2 = LOWER-CASE
C
C     NOTE - SOME KEYWORDS ARE DUPLICATED BECAUSE OF FORTRAN STATEMENT
C      TRANSLATIONS AS FROM FORTRAN II TO FORTRAN-77 FORMS.
C
C     PARAMETER (NKWD=87)
      PARAMETER (NKWD=90)
      CHARACTER CT
      CHARACTER*2 KWD(26,NKWD),JT
      CHARACTER*26 LC(2)
      DIMENSION KWL(NKWD),KSKP(NKWD)
      INTEGER CURCAS
      INCLUDE 'tidy.inc'
C
C     KWL IS LENGTH OF STRING TO COPY TO BUFFER, KSKP IS NUMBER OF
C      ORIGINAL CHARACTERS TO SKIP.
C     MOST KEYS HAVE EXTRA SPACE AT END FOR AUTOMATIC SPACING.
C
      DATA KWL( 1),(KWD(I,1),I=1,  7)/ 7 , 'A','C','C','E','P','T',' '/
      DATA KSKP(  1 ) / 6 /
      DATA KWL( 2),(KWD(I,2),I=1,  7)/ 7 , 'A','S','C','E','N','T',' '/
      DATA KSKP(  2 ) / 6 /
      DATA KWL( 3),(KWD(I,3),I=1,  7)/ 7 , 'A','S','S','I','G','N',' '/
      DATA KSKP(  3 ) / 6 /
      DATA KWL( 4),(KWD(I,4),I=1, 11)/11 , 'B','A','C','K','S','P','A'
     1,'C','E',' ','('/
      DATA KSKP(  4 ) / 10 /
      DATA KWL( 5),(KWD(I,5),I=1, 11)/11 , 'B','L','O','C','K',' ','D'
     1,'A','T','A',' '/
      DATA KSKP(  5 ) / 9 /
      DATA KWL( 6),(KWD(I,6),I=1, 10)/10 , 'B','U','F','F','E','R',' '
     1,'I','N',' '/
      DATA KSKP(  6 ) / 8 /
      DATA KWL( 7),(KWD(I,7),I=1, 11)/11 , 'B','U','F','F','E','R',' '
     1,'O','U','T',' '/
      DATA KSKP(  7 ) /  9 /
      DATA KWL( 8),(KWD(I,8),I=1,  5)/ 5 , 'C','A','L','L',' '/
      DATA KSKP(  8 ) / 4 /
      DATA KWL( 9),(KWD(I, 9),I=1, 10)/10 , 'C','H','A','R','A','C','T'
     1,'E','R',' '/
      DATA KSKP(  9 ) / 9 /
      DATA KWL(10),(KWD(I,10),I=1, 7)/ 7 , 'C','O','M','M','O','N',' '/
      DATA KSKP( 10 ) / 6 /
      DATA KWL(11),(KWD(I,11),I=1,  8)/ 8 , 'C','O','M','P','L','E','X'
     1,' '/
      DATA KSKP( 11 ) / 7 /
      DATA KWL(12),(KWD(I,12),I=1,  9)/ 9 , 'C','O','N','T','I','N','U'
     1,'E',' '/
      DATA KSKP( 12 ) / 8 /
      DATA KWL(13),(KWD(I,13),I=1,  5)/ 5 , 'D','A','T','A',' '/
      DATA KSKP( 13 ) / 4 /
      DATA KWL(14),(KWD(I,14),I=1,  8)/ 8 , 'D','E','C','O','D','E',' '
     1,'('/
      DATA KSKP( 14 ) / 7 /
      DATA KWL(15),(KWD(I,15),I=1, 10)/10 , 'D','I','M','E','N','S','I'
     1,'O','N',' '/
      DATA KSKP( 15 ) / 9 /
      DATA KWL(16),(KWD(I,16),I=1, 17)/17 , 'D','O','U','B','L','E',' '
     1,'P','R','E','C','I','S','I','O','N',' '/
      DATA KSKP( 16 ) / 15 /
      DATA KWL(17),(KWD(I,17),I=1, 17)/17 , 'D','O','U','B','L','E',' '
     1,'P','R','E','C','I','S','I','O','N',' '/
      DATA KSKP( 17 ) / 6 /
      DATA KWL(18),(KWD(I,18),I=1,8)/8 , 'E','N','C','O','D','E',' ',
     1'('/
      DATA KSKP( 18 ) / 7 /
      DATA KWL(19),(KWD(I,19),I=1, 10)/10 , 'E','N','D',' ','F','I','L'
     1,'E',' ','('/
      DATA KSKP( 19 ) / 8 /
      DATA KWL(20),(KWD(I,20),I=1,  6)/ 6 , 'E','N','D',' ','I','F'/
      DATA KSKP( 20 ) / 5 /
      DATA KWL(21),(KWD(I,21),I=1,  9)/ 9 , 'E','N','D',' ','F','I','L'
     1,'E',' '/
      DATA KSKP( 21 ) / 7 /
      DATA KWL(22),(KWD(I,22),I=1,  6)/ 6 , 'E','N','T','R','Y',' '/
      DATA KSKP( 22 ) / 5 /
      DATA KWL(23),(KWD(I,23),I=1, 12)/12 , 'E','Q','U','I','V','A','L'
     1,'E','N','C','E',' '/
      DATA KSKP( 23 ) / 11 /
      DATA KWL(24),(KWD(I,24),I=1,  8)/ 8 , 'E','X','T','E','R','N','A'
     1,'L'/
      DATA KSKP( 24 ) / 8 /
      DATA KWL(25),(KWD(I,25),I=1,  6)/ 6 , 'F','I','N','I','S',' '/
      DATA KSKP( 25 ) / 5 /
      DATA KWL(26),(KWD(I,26),I=1,  8)/ 8 , 'F','O','R','M','A','T',' '
     1,'('/
      DATA KSKP( 26 ) / 7 /
      DATA KWL(27),(KWD(I,27),I=1,  7)/ 7 , 'F','O','R','T','R','A','N'/
      DATA KSKP( 27 ) / 7 /
C
C     ONLY NEED 3 CHARS OF THIS.
      DATA KWL(28),(KWD(I,28),I=1,  8)/ 3 , 'I','F',' ','(','U','N','I'
     1,'T'/
      DATA KSKP( 28 ) / 7 /
      DATA KWL(29),(KWD(I,29),I=1,  9)/ 9 , 'F','U','N','C','T','I','O'
     1,'N',' '/
      DATA KSKP( 29 ) / 8 /
      DATA KWL(30),(KWD(I,30),I=1,  7)/ 7 , 'G','O',' ','T','O',' ','('
     1/
      DATA KSKP( 30 ) / 5 /
      DATA KWL(31),(KWD(I,31),I=1,  6)/ 6 , 'G','O',' ','T','O',' '/
      DATA KSKP( 31 ) / 4 /
      DATA KWL(32),(KWD(I,32),I=1, 26)/26 , 'I','F',' ','(','A','C','C'
     1,'U','M','U','L','A','T','O','R',' ','O','V','E','R','F','L','O',
     1'W',')'     ,' '/
      DATA KSKP( 32 ) / 23 /
      DATA KWL(33),(KWD(I,33),I=1, 23)/23 , 'I','F',' ','(','Q','U','O'
     1,'T','I','E','N','T',' ','O','V','E','R','F','L','O','W',')',' '/
      DATA KSKP( 33 ) / 20 /
      DATA KWL(34),(KWD(I,34),I=1, 18)/18 , 'I','F',' ','(','D','I','V'
     1,'I','D','E',' ','C','H','E','C','K',')',' '/
      DATA KSKP( 34 ) / 15 /
      DATA KWL(35),(KWD(I,35),I=1, 11)/11 , 'I','F',' ','(','E','N','D'
     1,'F','I','L','E'/
      DATA KSKP( 35 ) / 10 /
      DATA KWL(36),(KWD(I,36),I=1, 16)/16 , 'I','F',' ','(','S','E','N'
     1,'S','E',' ','L','I','G','H','T',' '/
      DATA KSKP( 36 ) / 13 /
      DATA KWL(37),(KWD(I,37),I=1, 17)/17 , 'I','F',' ','(','S','E','N'
     1,'S','E',' ','S','W','I','T','C','H',' '/
      DATA KSKP( 37 ) / 14 /
      DATA KWL(38),(KWD(I,38),I=1,  3)/ 3 , 'I','F',' '/
      DATA KSKP( 38 ) / 2 /
      DATA KWL(39),(KWD(I,39),I=1,  8)/ 8 , 'I','N','T','E','G','E','R'
     1,' '/
      DATA KSKP( 39 ) / 7 /
      DATA KWL(40),(KWD(I,40),I=1,  8)/ 8 , 'L','O','G','I','C','A','L'
     1,' '/
      DATA KSKP( 40 ) / 7 /
      DATA KWL(41),(KWD(I,41),I=1,  8)/ 8 , 'M','A','C','H','I','N','E'
     1,' '/
      DATA KSKP( 41 ) / 7 /
      DATA KWL(42),(KWD(I,42),I=1,  9)/ 9 , 'N','A','M','E','L','I','S'
     1,'T',' '/
      DATA KSKP( 42 ) / 8 /
      DATA KWL(43),(KWD(I,43),I=1,  6)/ 6 , 'P','A','U','S','E',' '/
      DATA KSKP( 43 ) / 5 /
      DATA KWL(44),(KWD(I,44),I=1,  6)/ 6 , 'P','R','I','N','T',' '/
      DATA KSKP( 44 ) / 5 /
      DATA KWL(45),(KWD(I,45),I=1,  8)/ 8 , 'P','R','O','G','R','A','M'
     1,' '/
      DATA KSKP( 45 ) / 7 /
      DATA KWL(46),(KWD(I,46),I=1,  6)/ 6 , 'P','U','N','C','H',' '/
      DATA KSKP( 46 ) / 5 /
C     READ INPUT TAPE TRANSLATES TO READ (
      DATA KWL(47),(KWD(I,47),I=1,  6)/ 6 , 'R','E','A','D',' ','('/
C     DATA KWL(47),(KWD(I,47),I=1, 16)/16 , 'R','E','A','D',' ','I','N'
C    1,'P','U','T',' ','T','A','P','E',' '/
      DATA KSKP( 47 ) / 13 /
      DATA KWL(48),(KWD(I,48),I=1,  6)/ 6 , 'R','E','A','D',' ','('/
C     DATA KWL(48),(KWD(I,48),I=1, 10)/10 , 'R','E','A','D',' ','T','A'
C    1,'P','E',' '/
      DATA KSKP( 48 ) / 8 /
      DATA KWL(49),(KWD(I,49),I=1,  6)/ 6 , 'R','E','A','D',' ','('/
      DATA KSKP( 49 ) / 5 /
      DATA KWL(50),(KWD(I,50),I=1,  5)/ 5 , 'R','E','A','D',' '/
      DATA KSKP( 50 ) / 4 /
      DATA KWL(51),(KWD(I,51),I=1,  5)/ 5 , 'R','E','A','L',' '/
      DATA KSKP( 51 ) / 4 /
      DATA KWL(52),(KWD(I,52),I=1,  7)/ 7 , 'R','E','T','U','R','N',' '/
      DATA KSKP( 52 ) / 6 /
      DATA KWL(53),(KWD(I,53),I=1,  8)/ 8 , 'R','E','W','I','N','D',' '
     1,'('/
      DATA KSKP( 53 ) / 7 /
      DATA KWL(54),(KWD(I,54),I=1,  8)/ 8 , 'S','E','G','M','E','N','T'
     1,' '/
      DATA KSKP( 54 ) / 7 /
      DATA KWL(55),(KWD(I,55),I=1, 12)/12 , 'S','E','N','S','E',' ','L'
     1,'I','G','H','T',' '/
      DATA KSKP( 55 ) / 10 /
      DATA KWL(56),(KWD(I,56),I=1,  5)/ 5 , 'S','T','O','P',' '/
      DATA KSKP( 56 ) / 4 /
      DATA KWL(57),(KWD(I,57),I=1, 11)/11 , 'S','U','B','R','O','U','T'
     1,'I','N','E',' '/
      DATA KSKP( 57 ) / 10 /
      DATA KWL(58),(KWD(I,58),I=1,  5)/ 5 , 'T','Y','P','E',' '/
      DATA KSKP( 58 ) / 4 /
      DATA KWL(59),(KWD(I,59),I=1,  7)/ 7 , 'W','R','I','T','E',' ','('/
C     DATA KWL(59),(KWD(I,59),I=1, 18)/18 , 'W','R','I','T','E',' ','O'
C    1,'U','T','P','U','T',' ','T','A','P','E',' '/
      DATA KSKP( 59 ) / 15 /
      DATA KWL(60),(KWD(I,60),I=1,  7)/ 7 , 'W','R','I','T','E',' ','('/
C     DATA KWL(60),(KWD(I,60),I=1, 11)/11 , 'W','R','I','T','E',' ','T'
C    1,'A','P','E',' '/
      DATA KSKP( 60 ) / 9 /
      DATA KWL(61),(KWD(I,61),I=1,  7)/ 7 , 'W','R','I','T','E',' ','('/
      DATA KSKP( 61 ) / 6 /
      DATA KWL(62),(KWD(I,62),I=1,  8)/ 8 , 'O','V','E','R','L','A','Y'
     1,' '/
      DATA KSKP( 62 ) / 7 /
      DATA KWL(63),(KWD(I,63),I=1,  6)/ 6 , 'I','D','E','N','T',' '/
      DATA KSKP( 63 ) / 5 /
      DATA KWL(64),(KWD(I,64),I=1, 10)/10 , 'F','R','E','Q','U','E','N'
     1,'C','Y',' '/
      DATA KSKP( 64 ) / 9 /
      DATA KWL(65),(KWD(I,65),I=1,  9)/ 9 , 'I','M','P','L','I','C','I'
     1,'T',' '/
      DATA KSKP( 65 ) / 8 /
      DATA KWL(66),(KWD(I,66),I=1,  6)/ 6 , 'L','E','V','E','L',' '/
      DATA KSKP( 66 ) / 5 /
      DATA KWL(67),(KWD(I,67),I=1,8)/5 , 'E','L','S','E',' ','I','F'
     1,' '/
      DATA KSKP( 67 ) / 4 /
      DATA KWL(68),(KWD(I,68),I=1,  5)/ 5 , 'E','L','S','E',' '/
      DATA KSKP( 68 ) / 4 /
      DATA KWL(69),(KWD(I,69),I=1,  5)/ 5 , 'T','H','E','N',' '/
      DATA KSKP( 69 ) / 4 /
      DATA KWL(70),(KWD(I,70),I=1, 7)/ 7 , 'C','L','O','S','E',' ','('/
      DATA KSKP( 70 ) / 6 /
      DATA KWL(71),(KWD(I,71),I=1,  8)/ 8 , 'I','N','C','L','U','D','E'
     1,' '/
      DATA KSKP( 71 ) / 7 /
      DATA KWL(72),(KWD(I,72),I=1,  9)/ 9 , 'I','N','Q','U','I','R','E'
     1,' ','('/
      DATA KSKP( 72 ) / 8 /
      DATA KWL(73),(KWD(I,73),I=1, 10)/10 , 'I','N','T','R','I','N','S'
     1,'I','C',' '/
      DATA KSKP( 73 ) / 9 /
      DATA KWL(74),(KWD(I,74),I=1,  6)/ 6 , 'O','P','E','N',' ','('/
      DATA KSKP( 74 ) / 5 /
      DATA KWL(75),(KWD(I,75),I=1, 10)/10 , 'P','A','R','A','M','E','T'
     1,'E','R',' '/
      DATA KSKP( 75 ) / 9 /
      DATA KWL(76),(KWD(I,76),I=1,  5)/ 5 , 'S','A','V','E',' '/
      DATA KSKP( 76 ) / 4 /
      DATA KWL(77),(KWD(I,77),I=1,  9)/ 9 , 'B','A','C','K','S','P','A'
     1,'C','E'/
      DATA KSKP( 77 ) / 9 /
      DATA KWL(78),(KWD(I,78),I=1,  4)/ 4 , 'E','N','D',' '/
      DATA KSKP( 78 ) / 3 /
      DATA KWL(79),(KWD(I,79),I=1,  7)/ 7 , 'R','E','W','I','N','D',' '/
      DATA KSKP( 79 ) / 6 /
      DATA KWL(80),(KWD(I,80),I=1,  6)/ 6 , 'C','L','O','S','E',' '/
      DATA KSKP( 80 ) / 5 /
      DATA KWL(81),(KWD(I,81),I=1,  7)/ 7 , 'E','N','D',' ','D','O',' '
     1/
      DATA KSKP( 81 ) / 5 /
      DATA KWL(82),(KWD(I,82),I=1, 10)/10 , 'D','O',' ','W','H','I','L'
     1,'E',' ','('/
      DATA KSKP( 82 ) / 8 /
      DATA KWL(83),(KWD(I,83),I=1,  7)/ 7 , 'R','E','P','E','A','T',' '/
      DATA KSKP( 83 ) / 6 /
      DATA KWL(84),(KWD(I,84),I=1, 13)/13 , 'S','E','N','S','E',' ','S'
     1,'W','I','T','C','H',' '/
      DATA KSKP( 84 ) / 11 /
      DATA KWL(85),(KWD(I,85),I=1,  3)/ 3 , 'T','O',' '/
      DATA KSKP( 85 ) / 2 /
      DATA KWL(86),(KWD(I,86),I=1,  3)/ 3 , 'D','O',' '/
      DATA KSKP( 86 ) / 2 /
      DATA KWL(87),(KWD(I,87),I=1,  9)/ 9 , 'C','O','N','T','I','N','U'
     1,'E',' '/
      DATA KSKP( 87 ) / 0 /
C
      DATA KWL(88),(KWD(I,88),I=1,  6)/ 6 , 'C','Y','C','L','E',' '/
      DATA KSKP( 88 ) / 5 /
      DATA KWL(89),(KWD(I,89),I=1,  5)/ 5 , 'E','X','I','T',' '/
      DATA KSKP( 89 ) / 4 /
      DATA KWL(90),(KWD(I,90),I=1,  5)/ 5 , 'N','O','N','E',' '/
      DATA KSKP( 90 ) / 4 /
C
C     CURCAS = CURRENT CASE: 1 = UPPER, 2 = LOWER
      DATA CURCAS/1/
      DATA LC/'ABCDEFGHIJKLMNOPQRSTUVWXYZ','abcdefghijklmnopqrstuvwxyz'/
C
C     NEGATIVE VALUES MEAN RESET CASE OF KEYWORDS.
      IF (KW.LE.0) THEN
          NEWCAS=abs(KW)
          IF (NEWCAS.NE.CURCAS.AND.NEWCAS.GT.0.AND.NEWCAS.LE.2) THEN
C          CHANGE CASE OF KEYWORD TABLE.
              DO 20 J=1,NKWD
                  DO 10 I=1,KWL(J)
C
C     FIND POSITION OF CHARACTER IN CURRENT-CASE ALPHABET
                      CT=KWD(I,J)(1:1)
                      K=INDEX(LC(CURCAS),CT)
C
C     IF FOUND, CHANGE TO NEW-CASE EQUIVALENT, OTHERWISE IGNORE.
                      IF (K.GT.0) THEN
                          KWD(I,J)=LC(NEWCAS)(K:K)
                      END IF
 10               CONTINUE
 20           CONTINUE
              CURCAS=NEWCAS
          END IF
          RETURN
      END IF
C
C     WRITE A KEYWORD AND SKIP THE ORIGINAL TEXT
C
C
C     COPY THE KEYWORD TO IOUT
      DO 25 I=1,KWL(KW)
              ICOL=ICOL+1
              IOUT(ICOL)=KWD(I,KW)
 25   CONTINUE
C
C     SKIP NSKIP CHARS OF ORIGINAL TEXT.
      NSKIP=KSKP(KW)
      IF (NSKIP.GT.0) THEN
30        JT=JINT(JCOL)
          IF (JT.NE.KBL) THEN
              NSKIP=NSKIP-1
              IF (NSKIP.EQ.0) THEN
                   JCOL=JCOL+1
                   LCPY=JT
                   RETURN
              ENDIF
              IF (JT.EQ.KERM) THEN
                   LCPY=KERM
                   ICOL=ICOL-1
                   MEOF=0
                   RETURN
               END IF
          END IF
          JCOL=JCOL+1
          GO TO 30
      ENDIF
      RETURN
      END
      SUBROUTINE MOVSTR (J)
      INCLUDE 'tidy.inc'
C
C     ADDS 1 BYTE TO STRING BY SHIFTING UNPROCESSED CHARS RIGHT.
C     USED BY HOLSCN WHEN REPLICATING APOSTROPHES
C
      DO 10 I=JMAX,J,-1
           JINT(I+1)=JINT(I)
 10   CONTINUE
      JMAX=JMAX+1
      JINT(JMAX+1)=KERM
      J=J+1
      JCOL = JCOL+1
      RETURN
      END
      SUBROUTINE NOPRO
C
C     THIS SUBROUTINE EXECUTES A HIGH-SPEED SEARCH FOR AN END STATEMENT.
C     IF MP2 IS ON, CARD IMAGES ARE WRITTEN ON TAPE 1 FOR USE BY PASS2.
C     NO INTERNAL PROCESSING IS DONE ON THE STATEMENTS.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C     SET INITIAL VALUES.
C
      IF (MDEB.NE.0) WRITE (STDERR,5)
 5    FORMAT (' Rewinding scratch files - subroutine Nopro')
      REWIND SCFIL1
      REWIND SCFIL2
      NRT2=0
      NDEF=0
      KLASS=1
*
*     ITYPE=0
      L15=0
      IF (MP2.NE.0) THEN
C
C     WRITE OUT STATEMENT CURRENTLY IN JINT.
C
           IMAX=JMAX
           KLASS=2
           CALL WRBIN (SCFIL1,KILI,SERIAL,JINT)
           NRT1=1
           KLASS=3
           IF (JMAX.GT.72) CALL DIAGNO (28)
      END IF
      GO TO 20
C
C     READ AND COPY CARD IMAGES BY WAY OF KBUFF.
C
 10   CALL READER
 20   NREC=NREC+1
C
C     LOOK FOR LAST NON-BLANK CHARACTER ON CARD.
C
      I=72
 30   IF (KBUFF(I).EQ.KBL) THEN
           I=I-1
           IF (I.GT.7) GO TO 30
      END IF
      IMAX=I
C
C     LOOK FOR END STATEMENT IN INPUT BUFFER KBUFF
C
      J=3
      DO 40 I=7,IMAX
           K=I
           IF (KBUFF(I).NE.KBL) THEN
                IF (KBUFF(I).NE.KEND(J)) GO TO 50
                J=J-1
                IF (J.EQ.0) THEN
C     FOUND AN END CARD IF NEXT CHAR IS BLANK.
                     IF (KBUFF(K+1).EQ.KBL) KLASS=8
                     GO TO 50
                END IF
           END IF
 40   CONTINUE
C
C
C     WRITE OUT CARD IMAGE FOR PASS2.
C
 50   IF (MP2.NE.0) THEN
           CALL WRBIN (SCFIL1,KILI,SERIAL,KBUFF)
           NRT1=NRT1+1
      END IF
C
C     GET NEXT RECORD UNLESS END CARD OR EOF
      IF (IQUIT.NE.1.AND.KLASS.NE.8) GO TO 10
C
C     CLOSE FILE
      IF (MP2.NE.0) REWIND SCFIL1
C
C     LOAD BUFFER, KBUFF, BEFORE EXITING.
C
      IF (IQUIT.EQ.0) CALL READER
      RETURN
      END
      SUBROUTINE PAGE (N)
C
C     THIS SUBROUTINE DOES THE GENERAL PAGE COUNTING FOR TIDY WHILE
C     LIMITING THE OUTPUT TO MAXLIN LINES PER PAGE.
C
C          N>0 -- I WILL WRITE N LINES.  START A NEW PAGE IF NECESSARY.
C          N=0 -- START A NEW PAGE.
C          N<0 -- START A NEW PAGE IF .LT. -N LINES ARE LEFT.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
      DATA MAXLIN/56/
C
      IF (N.LT.0) THEN
C                            CONDITIONAL EJECT (NO LINES WRITTEN)
           IF ((LINE-N).LE.MAXLIN) RETURN
      ELSE IF (N.GT.0) THEN
           LINE=LINE+N
           IF (LINE.LE.MAXLIN) RETURN
      END IF
C                            MAKE NEW PAGE
      IF (LINE.NE.0) THEN
           LINE=0
           IF (N.GT.0) LINE=N
           NPAGE=NPAGE+1
           MPAGE=MPAGE+1
C          WRITE (OUTFIL,10) NROUT,IPASS,MPAGE,NPAGE,JOB
      END IF
      RETURN
C
C
C10   FORMAT (/' ',6X,'* T I D Y *          ROUTINE',I4,4X,'PASS',I2,2X,
C    1'PAGE',I3,21X,'PAGE',I4/7X,80A1/1X)
      END
      SUBROUTINE PASS1
C
C     THIS ROUTINE COLLECTS STATEMENT NUMBERS, MAKES DIAGNOSTIC COMMENTS
C     AND SETS UP THE FORTRAN STATEMENTS IN A FORM SUITABLE FOR PASS2.
C
      INTEGER JTMP(8)
      INCLUDE  'tidy.inc'
      INCLUDE  'units.inc'
*      CHARACTER*2 JNT,JT,ICH,KUPPER,PRVCPY,KCTRAN
      CHARACTER*2 JNT,ICH,KUPPER,PRVCPY
      COMMON /PS1SUB/ KSTC(5),NIFBLK
      DIMENSION KCNDO(MXPLBL)
      LOGICAL BAKSCN
C
C     A    B    C    D    E    F    G    H    I    J    K    L    M
C     1    2    3    4    5    6    7    8    9    10   11   12   13
C
C     N    O    P    Q    R    S    T    U    V    W    X    Y    Z
C     14   15   16   17   18   19   20   21   22   23   24   25   26
C
C     KSPK (SPECIAL CHARACTERS) ---
C     =    ,    (    /    )    +    -    *    .    $    -    '    & NONE
C     1    2    3    4    5    6    7    8    9    10   11   12   13  14
C
C
C     SET UP INITIAL CONDITIONS.
C     REWIND TAPE FILES 1 AND 2.
C
 10   REWIND  SCFIL1
      REWIND  SCFIL2
      DO 20 I=1,MXDO
          LDOS(I)=0
 20   CONTINUE
      IMAX=66*(NQCNTS+1)+6
      IPASS=1
      ICOL=0
      KOUNT=0
      MP2=1
      NBLC=2
      MPUN=KPUN
      MPRIN=KPRIN
      NROUT=NROUT+1
      NRT1=0
      NRT2=0
      MILDO=0
      MLGC=-1
      MSKP=0
      MPAGE=0
      MTRAN=0
      NDEF=0
      NDOS=0
      NFORT=0
      NREC=0
      NREF=0
      L25=0
      NTRAN=0
      NXEQ=0
      NIFBLK=0
      KENDDO=100000
      KCNDP=0
      GO TO 50
C
C                  ILLEGAL FIRST CHARACTER.
 30   JGOOF=9
C                  WRITE DIAGNOSTIC
 40   CALL DIAGNO (JGOOF)
C                  GET NEW CARD.
C     (UNLESS EOF ALREADY)
 50   IF (IQUIT.NE.0) GO TO 850
      CALL SKARD
      NXRF=1
      IF (IMAX.LT.ICOL) IMAX=ICOL
      DO 60 I=1,IMAX
          IOUT(I)=KBL
 60   CONTINUE
      IMAX=0
* new test
*     IF (JINT(1).EQ.KSPK(15)) GO TO 130
C
C     LOOK FOR * IN COLUMN 1
C
      IF (JINT(1).EQ.KSPK(8)) THEN
          CALL CONTRL
          IF (ISTAR.LT.0) THEN
C                  CONTROL CARD FOUND AND EXECUTED.
              IF (MSTOP.NE.0) THEN
C                            *STOP CARD FOUND. QUIT IF FIRST OF ROUTINE
                  IF (NFORT.LE.0) THEN
                      MP2=0
                      RETURN
                  ELSE
C                            OTHERWISE BUILD AN END CARD
                      GO TO 830
                  END IF
              END IF
              IF (MSKP.EQ.0) GO TO 50
              MP2=0
              CALL NOPRO
              GO TO 10
C                  CONTROL CARD FOR DELAYED EXECUTION. SAVE FOR PASS 2.
          ELSE IF (ISTAR.EQ.0) THEN
              KLASS=0
              GO TO 140
          ELSE
C                  * IN COL 1. NOT A CONTROL CARD.  PUT OUT LITERALLY
C                  UNLESS * IN COL 2. ALSO.
              IF (JINT(2).EQ.KSPK(8)) GO TO 50
              IF (MNDOO.eq.1) JINT(1)=KSPK(15)
              GO TO 130
          END IF
      END IF
C
C     *STOP COMMAND EXIT.
C
C     NO * IN COLUMN 1, LOOK FOR C, D, I, F, ., OR $. (UPPER CASE)
C
C
      IF (JINT(1).EQ.KBL) GO TO 160
      JNT=KUPPER(JINT(1))
C
C     COMMENT CARD
C      IF (JNT.EQ.KABC(3)) THEN
* new test
       IF (JNT.EQ.KABC(3).OR.JNT.EQ.KSPK(15)) THEN
          IF (MCOM.EQ.0) GO TO 50
*          IF (MCASE.EQ.0) JINT(1)=KCTRAN(JINT(1))
          IF (MNDOO.eq.1) JINT(1)=KSPK(15)
C
C     SHIFTED COMMENTS OPTION - SHIFT COMMENTS TO BEGIN IN COLUMN MCOM
C      (UNLESS * IN FIELD)
          IF (MCOM.GT.0) THEN
              JCLMX=min(JMAX,MCOM-1)
              DO 90 JCOL=2,JCLMX
                  IF (JINT(JCOL).NE.KBL) THEN
                      IF (JINT(JCOL).EQ.KSPK(8)) GO TO 120
C
C                         SHIFT CARD TO RIGHT
C
                      ICOL=MCOM
                      DO 70 I=JCOL,JMAX
                          ICOL=ICOL+1
                          IOUT(ICOL)=JINT(I)
 70                   CONTINUE
                      IOUT(1)=JINT(1)
                      IF (ICOL.GT.72) ICOL=72
                      IMAX=ICOL
                      DO 80 I=JCOL,JCLMX
                          IOUT(I)=KBL
 80                   CONTINUE
                      KLASS=1
                      JTYPE=0
                      L15=0
                      NBLC=0
C     WRITE THE PROCESSED RECORD.  RAW COMMENTS USE INPUT INSTEAD.
                      GO TO 510
                  END IF
 90           CONTINUE
          END IF
C
C     LOOK FOR BLANK COMMENT
C
          DO 100 JCOL=2,JMAX
              IF (JINT(JCOL).NE.KBL) GO TO 120
 100      CONTINUE
C
C     BLANK COMMENT. TEST IF TWO PREVIOUS CARDS WERE BLANK
C
          NBLC=NBLC+1
          IF (NBLC.GT.2) GO TO 50
          IOUT(1)=JINT(1)
          JMAX=7
          GO TO 130
      END IF
C
C     A BLANK LINE PRESERVED AS A COMMENT WITH NON-PRINTING FIRST CHAR
C      (SET IN SUBROUTINE READER IF *NOSTRIP OPTION TURNED ON)
      IF (JINT(1).EQ.KBLCMT) GO TO 120
C
      IF (JNT.EQ.KABC(4).OR.JNT.EQ.KABC(9).OR.JNT.EQ.KABC(6)) THEN
          CALL DIAGNO (8)
          GO TO 50
      END IF
C
C     LOOK FOR ANY SPECIAL CHARACTER IN COLUMN 1
      DO 110 I=1,14
          IF (JNT.EQ.KSPK(I)) THEN
C
C     SPECIAL CHAR IN COL 1.  GIVE MSG AND TREAT AS COMMENT
C
              CALL DIAGNO (30)
              GO TO 130
          END IF
 110  CONTINUE
      GO TO 160
C
C     NON-BLANK COMMENT.
C
 120  NBLC=0
      IF (JMAX.GT.72.AND.MSER.NE.0) JMAX=72
C
C     COMMENT CARD.  DO WE SAVE THEM...
 130  KLASS=1
 140  JTYPE=0
C
C     WRITE STATEMENT IMAGE ON TAPE 1 FOR PASS 2.
C
 150  L15=0
      IMAX=JMAX
      CALL WRBIN (SCFIL1,KILI,SERIAL,JINT)
      NRT1=NRT1+1
      GO TO 50
C
C               ===============================================
C               *                                             *
C               *      START PROCESSING OF FORTRAN CARDS      *
C               *                                             *
C               ===============================================
C
 160  IF (JMAX.LT.8) GO TO 40
      NFORT=NFORT+1
C     CLASSIFY STATEMENT, THEN CHECK AND CHANGE HOLLERITH FIELDS
C       (DO UNCLASSIFIED (REPLACEMENT, ETC) STATEMENTS, AND ALSO
C       THOSE IN WHICH STRINGS ARE LEGAL PARTS.
      ITYPE=0
      JCOL=6
      CALL KWSCAN (ITYPE,KSTC)
      MPASS1=1
      I=KSTC(5)
      KLASS=KSTC(2)
      NINS=KSTC(1)
      CALL HOLSCN (ITYPE,I,LNGST)
C                  CLEAR FLAGS
      MLGC=-1
      NTRAN=MTRAN
      MTRAN=0
      MEOF=-1
      JGOOF=1
C                  CLEAR STATEMENT AND REFERENCE NUMBERS
      L15=0
      L772=0
C                  CLEAR BLANK COMMENT COUNTER
      NBCOLD=NBLC
      NBLC=0
C                  SET POSITION COUNTERS.
      JCOL=7
      IF (JUST.EQ.0) THEN
C                            NO COLUMN SHIFT
          ICOL=6
 170      IF (JINT(JCOL).NE.KBL) GO TO 180
          JCOL=JCOL+1
          ICOL=ICOL+1
          GO TO 170
      END IF
C                            COLUMN=SOMETHING
      ICOL=JUST-1
C                            ADD INDENT
 180  ICOL=ICOL+INDENT*(NDOS+NIFBLK)
      ICOL=min(ICOL,MXRGHT)
C                            REMEMBER THE STARTING COLUMN
      ICOLSV=ICOL
C                  ANALYSIS OF LOGICAL IF RE-ENTERS HERE.
C
C                  SELECT NEXT COURSE ON BASIS OF FIRST SPECIAL CH.
C             =   ,   (   /  )  +  -  *   .  $  -  '  &  NONE
 190  GO TO (240,350,200,400,30,30,30,400,30,30,30,400,30,400),IFIR
C
C                  FIRST IS (.  LOOK FOR )
 200  NPAR=0
      DO 210 NF=LFIR,JMAX
          IF (JINT(NF).EQ.KSPK(5)) NPAR=NPAR-1
          IF (JINT(NF).EQ.KSPK(3)) NPAR=NPAR+1
          IF (NPAR.EQ.0) GO TO 220
 210  CONTINUE
C                            MISSING )
      JGOOF=2
      GO TO 40
C                  THIS IS THE END OF THE FIRST STACK OF PARENS.
C                  SKIP BLANKS.
C                  FIRST LOOK FOR DO WHILE STATEMENT
 220  IF (KLASS.EQ.3) GO TO 400
      KJ=82
      CALL KWSCAN (KJ,KSTC)
      IF (KJ.EQ.82) GO TO 1390
C
 230  NF=NF+1
      IF (NF.GE.JMAX) GO TO 400
      IF (JINT(NF).EQ.KBL) GO TO 230
C
C                  CHARACTER REPLACEMENT STATEMENTS CAN HAVE 2 SETS OF
C                  PARENS BEFORE =.
      IF (JINT(NF).EQ.KSPK(3)) THEN
          LFIR=NF
          GO TO 200
      END IF
C
      IF (JINT(NF).EQ.KSPK(1)) THEN
C           IF NEXT CHARACTER IS = PROCESS AS ARITHMETIC REPLACEMENT.
          LQUAL=NF
          GO TO 320
      END IF
C           OTHERWISE, PROCESS AS FORTRAN STATEMENT
      GO TO 400
C
C                  FIRST SPECIAL CH. IS =.
 240  LQUAL=LFIR
C                  IS IT A DO STATEMENT.  IF NOT, GO TO ARITHMETIC PROC.
C                  LOOK FOR -D- -O-
      ICH=KABC(4)
      DO 250 J=7,JMAX
          JNT=KUPPER(JINT(J))
          IF (JNT.EQ.KBL) GO TO 250
          IF (JNT.NE.ICH) GO TO 320
          IF (ICH.EQ.KABC(15)) GO TO 260
          ICH=KABC(15)
 250  CONTINUE
      GO TO 320
C                  FOUND -D- -O- NOW LOOK FOR COMMAS.  ALLOW EXACTLY 1
C                  OR 2 COMMAS OUTSIDE OF PARENTHESES, 1 EQUALS.
C                  CERTAIN SPECIAL CHARACTERS NOT ALLOWED.
 260  NCOMA=0
      NLPS=0
      JJ=LQUAL+1
      DO 310 J=JJ,JMAX
          JNT=JINT(J)
          DO 270 I=1,14
              IF (JNT.EQ.KSPK(I)) GO TO (320,300,280,310,290,310,310,
     1         310,310,320,310,320,320,320),I
 270      CONTINUE
          GO TO 310
C
C     COUNT LEFT PARENTHESES
 280      NLPS=NLPS+1
          GO TO 310
C
C     COUNT RIGHT PARENTHESES
 290      NLPS=NLPS-1
          GO TO 310
C
C     A COMMA. DISREGARD IF INSIDE PARENTHESES, ABORT SCAN IF UNBALANCED
 300      IF (NLPS.LT.0) THEN
              GO TO 320
          ELSE IF (NLPS.EQ.0) THEN
              IF (NCOMA.GT.1) GO TO 320
              NCOMA=NCOMA+1
          END IF
 310  CONTINUE
C
      IF (NCOMA.EQ.0) GO TO 320
C                  O.K.  THIS IS A DO STATEMENT.
      KLASS=10
      JTYPE=14
      GO TO 430
C
C              =================================================
C              *                                               *
C              *   START PROCESSING OF ARITHMETIC STATEMENT.   *
C              *                                               *
C              =================================================
 320  KLASS=6
      JTYPE=0
C
C     IF IN ANSI MODE, CHECK LENGTH OF VARIABLE ON LEFT
      IF (MANSI.EQ.0) THEN
          IF (IFIR.EQ.1.OR.IFIR.EQ.3) THEN
              LNGVR=0
              DO 330 J=JCOL,LFIR-1
                  IF (JINT(J).NE.KBL) LNGVR=LNGVR+1
 330          CONTINUE
              IF (LNGVR.GT.6) CALL DIAGNO (41)
          END IF
      END IF
C
 340  CALL COPY (-1)
      IF (MEOF.LT.0) THEN
          GO TO 340
      ELSE IF (MEOF.GT.0.OR.LCPY.EQ.KERM) THEN
          IF (MLGC.NE.0) THEN
              JCOL=1
              CALL RSTAT
              L15=INT(L772)
          END IF
          GO TO 490
      ELSE
          ICOL=ICOL+1
          MEOF=-1
          GO TO 340
      END IF
C
C
C     DO STATEMENTS WITH COMMA BEFORE INDEX VARIABLE

C                  IS IT A DO STATEMENT.  IF NOT, GO TO ARITHMETIC PROC.
C                  LOOK FOR -D- -O-
C                  (UNLESS STATEMENT IS CLASSIFIED)
 350  IF (KLASS.EQ.0) THEN
          ICH=KABC(4)
          DO 360 J=JCOL,JMAX
              JNT=KUPPER(JINT(J))
              IF (JNT.EQ.KBL) GO TO 360
              IF (JNT.NE.ICH) GO TO 400
              IF (ICH.EQ.KABC(15)) THEN
                  JCOLD=JCOL
                  JCOL=J+1
                  GO TO 370
              END IF
              ICH=KABC(15)
 360      CONTINUE
          GO TO 400
C
C          CHECK FOR STATEMENT NUMBER, NEXT NON-BLANK SHOULD BE THE COMM
 370      CALL RSTAT
          IF (INT(L772).NE.0.AND.LFIR.EQ.JCOL) THEN
C          NOW CHECK FOR VARIABLE FOLLOWED BY EQUAL SIGN.  IF FOUND, CHA
C           COMMA TO BLANK AND USE POSITION OF = AS LQUAL, PROCESS AS DO
              JCOL=JCOL+1
              DO 390 J=JCOL,JMAX
                  JNT=JINT(J)
                  DO 380 I=1,13
                      IF (JNT.EQ.KSPK(I)) THEN
                          JCOL=JCOLD
                          IF (I.EQ.1) THEN
                              IFIR=I
                              JINT(LFIR)=KBL
                              LFIR=J
                              LQUAL=LFIR
                              GO TO 260
                          END IF
                          GO TO 400
                      END IF
 380              CONTINUE
 390          CONTINUE
          END IF
      END IF
C
C              ========================================
C              *                                      *
C              *     END OF ARITHMETIC PROCESSING     *
C              *  START FORTRAN STATEMENT PROCESSING  *
C              *                                      *
C              ========================================
C
C                  CHECK EVERY LISTED STATEMENT TYPE.
 400  IF (MPASS1.GT.1) THEN
C     MUST RE-CHECK REST OF IF-STATEMENTS
          ITYPE=0
          CALL KWSCAN (ITYPE,KSTC)
          IF (ITYPE.EQ.0) GO TO 480
      END IF
      NINS=KSTC(1)
      MPASS1=MPASS1+1
C
C                  FOUND IT.
      IF (ITYPE.NE.0) THEN
          KLASS=KSTC(2)
          JTYPE=KSTC(3)
          IF (IFIR.NE.12) THEN
C     COMPLAIN IF NON-ANSI STATEMENT.
              IF (MANSI.EQ.0.AND.KSTC(4).EQ.1) CALL DIAGNO (34)
              IF (MLGC.NE.0) GO TO 410
C                            FOLLOWS LOGICAL IF OR IS FUNCTION DECL.
              IF (KLASS.EQ.3.OR.KLASS.EQ.4.OR.KLASS.EQ.6.OR.KLASS.EQ.7
     1         .OR.KLASS.EQ.11) GO TO 450
              GO TO 40
          ELSE
C        COMPLAIN IF FIRST SPECIAL CHAR ' AND NOT INCLUDE OR PRINT
              IF (ITYPE.NE.71.AND.ITYPE.NE.43.AND.ITYPE.NE.44) GO TO 30
          END IF
      ELSE
C
C                  NOT IN TABLE.  PASS IT WITHOUT PROCESSING.
          CALL DIAGNO (30)
          KLASS=11
          JTYPE=0
      END IF
C
C                  THIS IS A FORTRAN STATEMENT.
C                  SET IMAX IN CASE THIS STATEMENT IS PUT OUT DIRECTLY.
 410  IMAX=JMAX
C                  CHECK FOR EXEMPT STATEMENT.
      IF (KLASS.EQ.3) THEN
          DO 420 J=1,6
              JINT(J)=KBL
 420      CONTINUE
          IF (MEX.EQ.0) GO TO 450
C                  THIS IS A NON-EXECUTABLE (KLASS 3.) FORTRAN STATEMENT
C                  AND THE EXEMPT FLAG IS SET.  SO PUT IT OUT DIRECTLY.
          GO TO 150
      END IF
C
C                  GET STATEMENT NUMBER UNLESS FOLLOWING LOGICAL IF.
      IF (MLGC.EQ.0) GO TO 450
 430  DO 440 I=1,5
          IF (JINT(I).NE.KBL) THEN
              J=INDEX('0123456789',JINT(I)(1:1))
              IF (J.LE.0) GO TO 450
              L15=L15*10+J-1
          END IF
 440  CONTINUE
C
C        IF THIS IS A WEIRD CARD, ALLOW A TRANSFER TO IT
 450  IF (KLASS.EQ.11) NTRAN=0
C
C     GO TO INDIVIDUAL STATEMENT PROCESSING BY JTYPE.
C
      I=JTYPE+1
      GO TO (520,530,560,570,580,590,600,630,660,710,720,740,760,760,
     1770,820,830,890,910,920,930,940,540,950,970,1020,1020,1020,1040,
     21070,1070,1090,1100,1110,1120,1130,1140,1140,1180,1220,1230,1240,
     31250,1080,1140,1140,1270,1350,1360,1370,1380,1390,460),I
C
C     ==================================================================
C     *                                                                *
C     *  AT THIS POINT, COMMENTS AND ARITHMETIC STATEMENTS HAVE BEEN   *
C     *  PROCESSED.  THE STATEMENTS HAVE BEEN CLASSIFIED AS ITYPE AND  *
C     *  KLASS.  THE LAST SYMBOL USED IN SCANNING THE FORTRAN STATE-   *
C     *  MENT IS KST(NINS,ITYPE), AND WAS FOUND AT JINT(LAST).  THE    *
C     *  FIRST SPECIAL CHARACTER, IF ANY, IS KSPK(IFIR), LOCATED AT    *
C     *  JINT(LFIR).  IF A STATEMENT                                   *
C     *  NUMBER IS PERMITTED, IT IS IN L15.  IF NOT, L15=0.            *
C     *  JCOL IS ON THE CURRENT CHARACTER IN THE INPUT STRING (THE     *
C     *  FIRST, UNLESS FOLLOWING A LOGICAL IF).  ICOL IS ON THE MOST   *
C     *  RECENT CHARACTER TO BE PUT INTO THE OUTPUT STRING (E.G. 6.)   *
C     *                                                                *
C     ==================================================================
C
C                  ILLEGAL JTYPE
 460  WRITE (OUTFIL,1420) JTYPE
      CALL DIAGNO (45)
C
C                  COPY REST OF CARD.
 470  ICOL=ICOL+1
 480  CALL COPY (0)
      IF (KLASS.LT.4) GO TO 500
C                  DLIST HANDLES THE STATEMENT NUMBER.
 490  CALL DLIST (MERR)
      IF (MERR.LT.0) GO TO 50
 500  IMAX=ICOL
C                  WRITE STATEMENT IMAGE ON TAPE1 FOR PASS 2.
 510  CALL WRBIN (SCFIL1,KILI,SERIAL,IOUT)
      NRT1=NRT1+1
      GO TO 50
C
C                  ***** JTYPE = 0
C     UNRECOGNIZED FORTRAN CARD
C                  COPY IT, INCLUDING BLANKS
 520  CALL COPY (0)
      GO TO 490
C
C                  ***** JTYPE = 1
C     ASCENT,MACHINE.
 530  I=0
      GO TO 550
C
C                  ***** JTYPE = 22
C     IDENT
C
 540  MP2=1
C            (MUST BE THE FIRST CARD OF THIS PASS.)
 550  IF (NFORT.NE.1) CALL DIAGNO (14)
      CALL DIAGNO (26)
      CALL NOPRO
      CALL HEADER
      RETURN
C
C                  ***** JTYPE = 2
C     ASSIGN
C
 560  CALL KWDWRT (ITYPE)
      CALL RSTAT
      CALL ADNUM (JRET)
      IF (JRET.EQ.1) GO TO 50
      ICOL=ICOL+1
      CALL KWDWRT (85)
      IF (MEOF.LT.0) GO TO 480
      GO TO 40
C
C                  ***** JTYPE = 3
C     BACKSPACE, EXTERNAL, IMPLICIT, PAUSE.
C
 570  CALL KWDWRT (ITYPE)
C     FINISH AN IMPLICIT STATEMENT
      IF (ITYPE.EQ.65) GO TO 400
      GO TO 470
C
C                  ***** JTYPE = 4
C      BLOCK DATA
C
 580  IF (NFORT.NE.1) GO TO 40
      CALL KWDWRT (ITYPE)
      GO TO 480
C
C                  ***** JTYPE = 5
C     BUFFER IN (I,P) (A,B) /// BUFFER OUT (I,P) (A,B)
C
 590  CALL KWDWRT (ITYPE)
      CALL COPY (-1)
      ICOL=ICOL+1
      CALL COPY (-1)
      IF (MEOF.LT.0.AND.JCOL.GT.JMAX) GO TO 490
      GO TO 40
C
C                  ***** JTYPE = 6
C     CALL   (FUNCTION,SUBROUTINE)
C
 600  JGOOF=10
      CALL KWDWRT (ITYPE)
      IF (IFIR.NE.3) GO TO 480
      NLPS=0
 610  CALL COPY (1)
      IF (LCPY.NE.KSPK(3)) THEN
          IF (LCPY.EQ.KSPK(5)) NLPS=NLPS-1
          IF (MEOF.LT.0) GO TO 610
          GO TO 40
      ELSE
          NLPS=NLPS+1
      END IF
      IOUT(ICOL)=KBL2
      JCOL=JCOL-1
 620  PRVCPY=LCPY
      CALL COPY (1)
      IF (MEOF.LT.0) THEN
          IF (LCPY.EQ.KALMRK) THEN
C     ALTERNATE RETURNS MUST BE PRECEDED BY , OR (
              IF (PRVCPY.NE.KSPK(2).AND.PRVCPY.NE.KSPK(3)) GO TO 620
C                            ARGUMENT IS *STATEMENT NUMBER
C     TRANSLATE ALTERNATE RETURN CODE IF DESIRED.
              IF (KALTRN.NE.KBL) IOUT(ICOL)=KALTRN
              CALL RSTAT
C
C     NO NUMBER LEGAL ONLY FOR FUNCTIONS AND SUBROUTINES.
              IF (INT(L772).EQ.0) THEN
                  IF (ITYPE.EQ.29.OR.ITYPE.EQ.57) GO TO 620
                  GO TO 40
              END IF
              CALL ADNUM (JRT)
              IF (JRT.EQ.1) GO TO 50
          END IF
C
C     ADD A SPACE BETWEEN PARAMETERS OR ARGUMENTS.
          IF (LCPY.EQ.KSPK(3)) THEN
              NLPS=NLPS+1
          ELSE IF (LCPY.EQ.KSPK(5)) THEN
              NLPS=NLPS-1
          ELSE IF (LCPY.EQ.KSPK(2).AND.NLPS.LE.2) THEN
              ICOL=ICOL+1
          END IF
C
          GO TO 620
      END IF
C
      IMAX=ICOL
      IF (NPAR.EQ.0) GO TO 490
      GO TO 40
C
C                  ***** JTYPE = 7
C      COMMON
C
 630  CALL KWDWRT (ITYPE)
C          J COUNTS SLASHES
      J=-2
      IF (IFIR.NE.4) GO TO 480
 640  IF (J.EQ.0) GO TO 470
      J=J+1
 650  CALL COPY (1)
      IF (LCPY.EQ.KSPK(4)) GO TO 640
      IF (MEOF.LT.0) GO TO 650
      CALL DIAGNO (11)
      GO TO 510
C
C                  ***** JTYPE = 8
C     CONTINUE
C
 660  JGOOF=12
      IF (L15.EQ.0) GO TO 40
      IF (MLGC.EQ.0) THEN
          DO 670 I=7,ICOL
              IOUT(I)=KBL
 670      CONTINUE
          ICOL=ICOLSV
          MLGC=-1
      END IF
      IF (MCONT.EQ.0) THEN
C                            IS THIS A DO-LOOP TERMINATOR...
          IF (NDOS.GT.0) THEN
              DO 680 I=1,NDOS
                  IF (L15.EQ.LDOS(I)) GO TO 700
 680          CONTINUE
          END IF
C
C                            CHECK IF TRANSFER TO IT BUT NOT DO-TERM.
          DO 690 I=1,NREF
              IF (L15.EQ.LREF(I)) GO TO 700
 690      CONTINUE
C                            COPY THE CARD
          CALL KWDWRT (ITYPE)
C                            PROCESS STATEMENT NUMBER
          CALL DLIST (MERR)
C                            SET A FLAG
          LDEF(NDEF)=-LDEF(NDEF)
          L25=L15
C                            TAKE TRANSFER STATUS OF LAST CARD
          MTRAN=NTRAN
C                            DONT SAVE STATEMENT FOR PASS2
          GO TO 50
      END IF
C                            THIS CONTINUE STATEMENT IS TO BE RETAINED
 700  IF (NDOS.NE.0) THEN
C                            IT TERMINATES THIS DO-LOOP. INDENT
C                            ONE LESS LEVEL
          IF (L15.EQ.LDOS(NDOS).AND.MLGC.NE.0) THEN
              ICOL=ICOL-INDENT
              ICOLSV=ICOL
          END IF
      END IF
      CALL KWDWRT (ITYPE)
      GO TO 490
C
C                  ***** JTYPE = 9
C     DATA
C
 710  CALL KWDWRT (ITYPE)
      IF (IFIR.EQ.4) THEN
          IF (JINT(JMAX).NE.KSPK(4).OR.LFIR.GE.JMAX) CALL DIAGNO (11)
      END IF
      GO TO 480
C
C                  ***** JTYPE = 10
C     DECODE (C,N,V) LIST  ///  ENCODE (C,N,V) LIST
C
 720  JGOOF=23
      CALL KWDWRT (ITYPE)
 730  CALL COPY (1)
      IF (LCPY.NE.KSPK(2)) THEN
          IF (MEOF.LT.0) GO TO 730
          GO TO 40
      END IF
      CALL RSTAT
      IF (INT(L772).GT.0) THEN
          CALL ADNUM (JRT)
          IF (JRT.EQ.1) GO TO 50
      END IF
      GO TO 1200
C
C                  ***** JTYPE = 11
C     DIMENSION
C
 740  JGOOF=13
      CALL KWDWRT (ITYPE)
      NPAR=-1
      DO 750 I=JCOL,JMAX
          CALL COPY (1)
          IF (NPAR.LT.0) THEN
              IF (LCPY.EQ.KSPK(3)) NPAR=NPAR+1
          ELSE IF (NPAR.EQ.0) THEN
              IF (LCPY.EQ.KSPK(5)) NPAR=NPAR+1
          ELSE
              IF (LCPY.NE.KSPK(2)) GO TO 750
              ICOL=ICOL+1
              NPAR=-1
          END IF
 750  CONTINUE
      IF (NPAR.GT.0) GO TO 500
      GO TO 40
C
C                  ***** JTYPE = 12
C     DOUBLE PRECISION
C                  ***** JTYPE = 13
C     DOUBLE, (CONVERT TO DOUBLE PRECISION).
C
 760  CALL KWDWRT (ITYPE)
      GO TO 400
C
C                  ***** JTYPE = 14
C     DO STATEMENT
C
 770  MILDO=1
      CALL KWDWRT (86)
      CALL RSTAT
C
C     IF NO STATEMENT, GIVE IT IMPOSSIBLE (FROM CARDS) NUMBER
C     KCNDO IS STACK OF CURRENTLY-OPEN ENDDO LOOPS
      IF (INT(L772).EQ.0) THEN
C          JUMP IF CONVERSION TO F-77 LOOP NOT DESIRED.
          IF (MNDOO.NE.0) GO TO 1400
          L772=KENDDO
          KCNDP=KCNDP+1
          KCNDO(KCNDP)=KENDDO
          KENDDO=KENDDO+1
      END IF
C
C     BE SURE IT DOESN'T REFERENCE BACKWARD IN PROGRAM.
      IF (NDEF.GT.0) THEN
          DO 780 I=1,NDEF
              IF (abs(LDEF(I)).EQ.INT(L772)) THEN
                  JGOOF=15
                  GO TO 40
              END IF
 780      CONTINUE
      END IF
C
C     ADD STATEMENT NUMBER TO DO-LIST.
C
      IF (NDOS.LT.0) CALL DIAGNO (44)
      IF (NDOS.GT.0) THEN
          IF (LDOS(NDOS).EQ.INT(L772)) GO TO 810
          IF (NDOS.GT.1) THEN
              DO 790 I=2,NDOS
                  IF (LDOS(I-1).EQ.INT(L772)) THEN
                      JGOOF=15
                      GO TO 40
                  END IF
 790          CONTINUE
          END IF
      END IF
C
C     ADD DO-LOOP TO OPEN LIST
      IF (NDOS.EQ.MXDO) THEN
          JGOOF=24
          MPUN=0
          MP2=0
          GO TO 40
      END IF
      NDOS=NDOS+1
      LDOS(NDOS)=INT(L772)
C
      IF (NREF.GT.0) THEN
          DO 800 I=1,NREF
              IF (LREF(I).EQ.INT(L772)) THEN
                  CALL DIAGNO (27)
                  GO TO 810
              END IF
 800      CONTINUE
      END IF
C
 810  CALL ADNUM (JRT)
      IF (JRT.EQ.1) GO TO 50
      GO TO 470
C
C     END DO-LOOP STATEMENT PROCESSING.
C
C
C                  ***** JTYPE = 15
C     END FILE
C
 820  IF (IFIR.NE.14) GO TO 30
      CALL KWDWRT (ITYPE)
      GO TO 480
C
C                  ***** JTYPE = 16
C     END STATEMENT.
C
C                   IS THERE A STATEMENT NUMBER TO USE?
 830  IF (L15.NE.0.OR.L25.NE.0) THEN
C                   YES. MAKE A CONTINUE CARD FOR IT TO FALL TO.
          ICOL=7
          CALL KWDWRT (87)
          MILDO=0
          CALL DLIST (MERR)
          IF (MERR.GE.0) THEN
              JTMP(1)=4
              JTMP(2)=8
              JTMP(3)=L15
              JTMP(4)=14
              JTMP(5)=MTRAN
              JTMP(6)=NXRF
              JTMP(7)=MEX
              JTMP(8)=ICOLSV
              CALL WRBIN (SCFIL1,JTMP,JINT(6),IOUT)
              NRT1=NRT1+1
          END IF
          L15=0
      END IF
      IF (NIFBLK.GT.0) CALL DIAGNO (33)
      IF (NDOS.NE.0) THEN
          CALL DIAGNO (16)
          CALL PAGE (1)
          WRITE (OUTFIL,1410) (LDOS(I),I=1,NDOS)
      END IF
C
C     IF STATEMENT IS NUMBERED AND REFERENCED
      IF (L15.NE.0.AND.NREF.GT.0) THEN
          DO 840 I=1,NREF
              IF (LREF(I).EQ.L15) THEN
                  CALL DIAGNO (18)
C                           GENERATE NEW STOP COMMAND.
                  CALL KWDWRT (56)
                  MILDO=-1
                  CALL DLIST (MERR)
                  IF (MERR.LT.0) GO TO 850
                  JTMP(1)=6
                  JTMP(2)=55
                  JTMP(3)=L15
                  JTMP(4)=10
                  JTMP(5)=MTRAN
                  JTMP(6)=NXRF
                  JTMP(7)=MEX
                  JTMP(8)=ICOLSV
                  CALL WRBIN (SCFIL1,JTMP,JINT(6),IOUT)
                  NRT1=NRT1+1
                  GO TO 850
              END IF
 840      CONTINUE
      END IF
C
C                       PROCESS FORMATS ON TAPE 2
 850  IF (NRT2.GT.0) THEN
          REWIND  SCFIL2
C                                  INSERT BLANK COMMENT CARD.
          IF (NBLC.EQ.0) THEN
              IOUT(1)=KABC(3)
              IF (MNDOO.eq.1) IOUT(1)=KSPK(15)
              DO 860 I=2,7
                  IOUT(I)=KBL
 860          CONTINUE
              KLASS=1
              ITYPE=0
              L15=0
              IMAX=7
              CALL WRBIN (SCFIL1,KILI,SERIAL,IOUT)
              NRT1=NRT1+1
          END IF
C                                TRANSFER FORMAT STATEMENTS
          DO 870 I=1,NRT2
              CALL RDBIN (SCFIL2,KILI,SERIAL,IOUT)
              ICOLSV=6
              NREC=JTYPE
              MILDO=1
              CALL DLIST (MERR)
              IF (MERR.GE.0) THEN
                  CALL WRBIN (SCFIL1,KILI,SERIAL,IOUT)
                  NRT1=NRT1+1
              END IF
 870      CONTINUE
          REWIND  SCFIL2
      END IF
C                                      MAKE END STATEMENT
      IF (NFEND.EQ.0.AND.NFORT.GT.0) THEN
          DO 880 I=1,6
              IOUT(I)=KBL
 880      CONTINUE
          CALL KWDWRT (78)
          KLASS=8
          ITYPE=20
          L15=0
          IMAX=9
          CALL WRBIN (SCFIL1,KILI,SERIAL,IOUT)
          NRT1=NRT1+1
      END IF
      REWIND  SCFIL1
      RETURN
C
C                 ==================================
C                 *   PASS1 NORMALLY EXITS HERE.   *
C                 ==================================
C
C
C                  ***** JTYPE = 17
C     EQUIVALENCE
C
 890  CALL KWDWRT (ITYPE)
 900  CALL COPY (-1)
      IF (MEOF.GE.0) GO TO 500
      GO TO 900
C
C                  ***** JTYPE = 18
C     FINIS.
C
 910  MSTOP=-1
      RETURN
C
C                  ***** JTYPE = 19
C     FORMAT (
C
 920  JGOOF=17
      CALL JTYP19 (JRTCOD,NBCOLD)
      GO TO (40,50,480),JRTCOD
C
C                  ***** JTYPE = 20
C     FORTRAN,ETC
C
 930  CALL KWDWRT (ITYPE)
      IMAX=JMAX
      GO TO 510
C
C                  ***** JTYPE = 21
C     FREQUENCY
C
 940  JGOOF=8
      GO TO 40
C
C                  ***** JTYPE = 23
C     GO TO (***,***),N
C
 950  JGOOF=19
      CALL KWDWRT (ITYPE)
      MILDO=1
      MTRAN=MLGC
C
C     PROCESS --GO TO LIST--.
C
 960  CALL RSTAT
      IF (INT(L772).EQ.0) GO TO 40
      CALL ADNUM (JRT)
      IF (JRT.EQ.1) GO TO 50
      CALL COPY (1)
      IF (LCPY.EQ.KSPK(2)) GO TO 960
      IF (LCPY.NE.KSPK(5)) GO TO 40
      CALL COPY (1)
      IF (LCPY.NE.KSPK(2)) THEN
          IOUT(ICOL+2)=IOUT(ICOL)
          IOUT(ICOL)=KSPK(2)
          ICOL=ICOL+2
      END IF
      GO TO 480
C
C                  ***** JTYPE = 24
C     GO TO ****
C
 970  JGOOF=19
      MILDO=-1
      CALL KWDWRT (ITYPE)
      CALL RSTAT
C
C     TEST REF STATEMENT FOR GO TO N OR GO TO N, (LIST)
C
      IF (INT(L772).EQ.0) GO TO 990
C
C     STATEMENT IS --GO TO 12345--.
C
      IF (L15.EQ.0.AND.L25.EQ.0) GO TO 980
      IF (MLGC.EQ.0) GO TO 980
C     LABELLED GOTO STATEMENT.
      IF (MCONT.EQ.0) THEN
          CALL DLIST (MERR)
          IF (MERR.LT.0) GO TO 40
C          SET UP REFERENCE TRANSLATION
          IF (NDEF.LT.MXPLBL) THEN
              NDEF=NDEF+1
              LDEF(NDEF)=0
              LOCDEF(NDEF)=INT(L772)
              L15=0
C               IF NO WAY TO GET HERE, DELETE IT
              IF (NTRAN.NE.0) GO TO 50
          END IF
      ELSE
          CALL DIAGNO (18)
      END IF
 980  MTRAN=MLGC
      CALL ADNUM (JRT)
      IF (JRT.EQ.1) GO TO 50
      GO TO 490
C
C     GO TO N OR GO TO N,LIST
C
 990  MTRAN=MLGC
      IF (IFIR.NE.2) THEN
C
C          STATEMENT IS --GO TO N--.
C
          IF (IFIR.EQ.14) GO TO 480
          GO TO 40
      END IF
C
C     GO TO N,(LIST)
C
 1000 CALL COPY (1)
      IF (LCPY.NE.KSPK(2)) GO TO 1000
      ICOL=ICOL+1
      CALL COPY (1)
      IF (LCPY.EQ.KSPK(3)) THEN
 1010     CALL RSTAT
          IF (INT(L772).NE.0) THEN
              CALL ADNUM (JRT)
              IF (JRT.EQ.1) GO TO 50
              CALL COPY (1)
              IF (LCPY.EQ.KSPK(2)) GO TO 1010
              IF (LCPY.EQ.KSPK(5)) GO TO 490
          END IF
      END IF
      GO TO 40
C
C                  ***** JTYPE = 25
C     IF ACCUMULATOR OVERFLOW (QUOTIENT, DIVIDE CHECK, END FILE, SENSE)
C                  ***** JTYPE = 26
C     IF QUOTIENT OVERFLOW
C                  ***** JTYPE = 27
C     IF(DIVIDE CHECK)
C
 1020 CALL KWDWRT (ITYPE)
C
C     PROCESS TWO-WAY TRANSFER.
C
 1030 CALL RSTAT
      IF (INT(L772).EQ.0) GO TO 40
      JGOOF=20
      MILDO=-1
      CALL ADNUM (JRT)
      IF (JRT.EQ.1) GO TO 50
      CALL COPY (1)
      IF (LCPY.NE.KSPK(2)) GO TO 40
      CALL RSTAT
      IF (INT(L772).EQ.0) GO TO 40
      GO TO 980
C
C
C                  ***** JTYPE = 28
C     IF(END FILE  I)
C
 1040 CALL KWDWRT (ITYPE)
      DO 1060 I=JCOL,JMAX
          IF (JINT(I).EQ.KSPK(5)) THEN
 1050         CALL COPY (1)
              IF (LCPY.EQ.KSPK(5)) GO TO 1030
              GO TO 1050
          END IF
 1060 CONTINUE
      JGOOF=20
      GO TO 40
C
C                  ***** JTYPE = 29
C     IF(SENSE LIGHT 5) 1,2
C                  ***** JTYPE = 30
C     IF(SENSE SWITCH 5) 1,2
C
 1070 JGOOF=20
      CALL KWDWRT (ITYPE)
      CALL COPY (2)
      ICOL=ICOL+1
      IF (LCPY.EQ.KSPK(5)) GO TO 1030
      GO TO 40
C
C                  ***** JTYPE = 43
C     ELSEIF
C
 1080 IF (NIFBLK.LE.0) THEN
          IOUT(1)=KABC(3)
          CALL DIAGNO (32)
      ELSE
          ICOL=ICOL-INDENT
          ICOLSV=ICOL
      END IF
      CALL KWDWRT (ITYPE)

C          FALL THRU TO IF
C
C                  ***** JTYPE = 31
C     IF (ARITHMETIC) 1,2,3   OR   IF (LOGICAL) STATEMENT.
C
 1090 JGOOF=20
      CALL JTYP31 (JRTCOD)
      GO TO (40,50,500,490,190),JRTCOD
C
C                  ***** JTYPE = 32
C     NAMELIST
C
 1100 JGOOF=21
      CALL KWDWRT (ITYPE)
      J=-1
      IF (IFIR.EQ.4) GO TO 640
      GO TO 40
C
C                  ***** JTYPE = 33
C     PRINT, TYPE, WRITE, PUNCH, READ, ACCEPT.
C
 1110 JGOOF=22
      CALL JTYP33 (ITYPE,JRTCOD)
      GO TO (480,40,470,50,490),JRTCOD
C
C                  ***** JTYPE = 34
C     SEGMENT,OVERLAY
C
 1120 NFORT=NFORT-1
      IF (NFORT.NE.0) CALL DIAGNO (14)
      CALL COPY (NINS)
      CALL HEADER
      IF (IFIR.EQ.3) GO TO 610
      GO TO 40
C                  ***** JTYPE = 35
C     PROGRAM, SUBROUTINE, FUNCTION.
C
 1130 IF (NFORT.NE.1) CALL DIAGNO (14)
      CALL KWDWRT (ITYPE)
      CALL HEADER
      NLPS=0
      IF (IFIR.EQ.3) GO TO 610
      GO TO 480
C
C
C                  ***** JTYPE = 44
C     WRITE OUTPUT TAPE
C                  ***** JTYPE = 36
C     READ INPUT TAPE
C                  ***** JTYPE = 45
C     WRITE TAPE
C                  ***** JTYPE = 37
C     READ TAPE
C
 1140 CALL KWDWRT (ITYPE)
C
C                  CONVERT TO CORRESPONDING READ/WRITE(I,N)LIST
      JGOOF=22
C                  COPY UNTIL COMMA
 1150 CALL COPY (1)
      IF (MEOF.GE.0) GO TO 40
      IF (LCPY.NE.KSPK(2)) GO TO 1150
C                  BINARY READ/WRITE  DO NOT HAVE STATEMENT NUMBERS.
      IF (JTYPE.EQ.36.OR.JTYPE.EQ.44) THEN
C                  PROCESS FORMAT STATEMENT NUMBER
          CALL RSTAT
          IF (INT(L772).NE.0) THEN
              CALL ADNUM (JRT)
              IF (JRT.EQ.1) GO TO 50
              CALL COPY (1)
              IF (LCPY.EQ.KSPK(2)) GO TO 1170
              IF (LCPY.EQ.KERM) THEN
C                      NO COMMA. END WITH )
                  IOUT(ICOL)=KSPK(5)
                  IMAX=ICOL
                  GO TO 490
              END IF
              GO TO 40
          END IF
C                      VARIABLE FORMAT--NO REFERENCE
          KLASS=6
 1160     CALL COPY (1)
C                      LOOK FOR COMMA
          IF (LCPY.EQ.KSPK(2)) GO TO 1170
          IF (MEOF.LT.0) GO TO 1160
      END IF
C                  REPLACE , BY ) AND GO PROCESS LIST
 1170 IOUT(ICOL)=KSPK(5)
      ICOL=ICOL+1
      GO TO 480
C
C                  ***** JTYPE = 38
C     READ ( AND WRITE (
C
 1180 JGOOF=23
 1190 CALL KWDWRT (ITYPE)
      NLPS=0
 1200 CALL COPY (1)
      IF (MEOF.GE.0) GO TO 40
C     LEFT PAREN MEANS START OF AN INTERNAL READ/WRITE SUBSCRIPT
      IF (LCPY.EQ.KSPK(3)) THEN
          NLPS=NLPS+1
          GO TO 1200
      END IF
C     RIGHT PAREN - COPY REST OF CARD UNLESS CLOSING SUBSCRIPT
      IF (LCPY.EQ.KSPK(5)) THEN
          IF (NLPS.LE.0) GO TO 470
          NLPS=NLPS-1
          GO TO 1200
      END IF
C     COMMA - NUMBER WILL FOLLOW UNLESS INTERNAL WRITE SUBSCRIPT
      IF (LCPY.EQ.KSPK(2)) THEN
          IF (NLPS.EQ.0) GO TO 1210
          GO TO 1200
      END IF
C     ACCEPT ANYTHING BUT = SIGN.
      IF (LCPY.NE.KSPK(1)) GO TO 1200
C
C     LAST CHARACTER WAS =.  CHECK KEYWORD FOR NUMBER FOLLOWING.
C      (SKIP FMT AND END FOR TYPE 47)
      IF (JTYPE.NE.47) THEN
C         FMT
          IF (BAKSCN(KABC(20),KABC(13))) GO TO 1210
C         END
          IF (BAKSCN(KABC(4),KABC(14))) GO TO 1210
C         ERR
      END IF
      IF (.NOT.BAKSCN(KABC(18),KABC(18))) GO TO 1200
C
C     GET STATEMENT NUMBER
C
 1210 CALL RSTAT
      IF (INT(L772).GT.0) THEN
          CALL ADNUM (JRT)
          IF (JRT.GT.0) GO TO 50
      END IF
      GO TO 1200
C
C                  ***** JTYPE = 39
C     RETURN
C
 1220 CALL KWDWRT (ITYPE)
      MTRAN=MLGC
      GO TO 480
C
C                  ***** JTYPE = 40
C     SENSE LIGHT
C
 1230 CALL KWDWRT (ITYPE)
      GO TO 470
C
C                  ***** JTYPE = 41
C     STOP
C
 1240 CALL KWDWRT (ITYPE)
      MILDO=-1
      MTRAN=MLGC
      GO TO 480
C
C                  ***** JTYPE = 42
C     IF (UNIT,N) L1,L2,L3,L4
C
 1250 CALL KWDWRT (ITYPE)
      CALL COPY (-1)
      IF (MEOF.LT.0) THEN
          ICOL=ICOL+1
          MILDO=1
          CALL DLIST (MERR)
          IF (MERR.GE.0) THEN
              DO 1260 I=1,4
                  CALL RSTAT
                  IF (INT(L772).EQ.0) GO TO 40
                  CALL ADNUM (JRT)
                  IF (JRT.EQ.1) GO TO 50
                  CALL COPY (1)
                  IF (LCPY.NE.KSPK(2)) THEN
                      IF (I.EQ.4.AND.LCPY.EQ.KERM) GO TO 500
                      GO TO 40
                  END IF
 1260         CONTINUE
          END IF
      END IF
      GO TO 40
C
C                        ***** JTYPE = 46
C     COMPLEX,  INTEGER,  REAL,  LOGICAL,  CHARACTER
C
 1270 CALL KWDWRT (ITYPE)
      KTDCL=0
C
C     CHECK IF HAS BYTE COUNT (REAL*8 ETC)
      IF (IFIR.EQ.8) THEN
          ICOL=ICOL-1
C          STATEMENT IS E.G. REAL*8, I.E. WITH BYTE NUMBER
C          FIRST SWALLOW ANY BLANKS BEFORE IT.
 1280     IF (JCOL.EQ.LFIR) GO TO 1290
          IF (JINT(JCOL).NE.KBL) GO TO 470
          JCOL=JCOL+1
          GO TO 1280
C
C     * WAS NEXT CHARACTER. COPY IT.
 1290     CALL COPY (1)
C
 1300     IF (JINT(JCOL).NE.KBL) THEN
C
C     PROCESS  *(*)
              IF (JINT(JCOL).EQ.KSPK(3)) THEN
                  CALL COPY (3)
                  GO TO 470
              END IF
              GO TO 1320
          END IF
          JCOL=JCOL+1
          GO TO 1300
C
C     GO PAST BYTE COUNT
 1310     CALL COPY (1)
 1320     DO 1330 I=1,10
              IP=INDEX('0123456789',JINT(JCOL)(1:1))
              IF (IP.GT.0) GO TO 1310
 1330     CONTINUE
          ICOL=ICOL+1
C
C     POSSIBLE VIOLATION OF ANSI STANDARD (REAL*8, ETC)
C      (ONLY LEGAL SIZE DECLARATION IS CHARACTER)
          IF (MANSI.EQ.0.AND.ITYPE.NE.9) KTDCL=1
      END IF
C
C     SEE IF IT IS A FUNCTION, IF SO ADD A SPACE AFTER
      I=29
      CALL KWSCAN (I,KSTC)
      IF (I.EQ.29) THEN
          CALL KWDWRT (I)
      ELSE
C
          IF (KTDCL.EQ.1) CALL DIAGNO (40)
C
C     LOOK FOR NON-ANSI INITIALIZED DECLARATIONS.
          IF (MANSI.EQ.0) THEN
              DO 1340 NF=LFIR,JMAX
                  IF (JINT(NF).EQ.KSPK(4)) THEN
                      CALL DIAGNO (42)
                      GO TO 470
                  END IF
 1340         CONTINUE
          END IF
      END IF
      GO TO 480
C
C                        ***** JTYPE = 47
C     OPEN, CLOSE, INQUIRE
 1350 JGOOF=31
      GO TO 1190
C
C                        ***** JTYPE = 48
C     ENDIF
 1360 NIFBLK=NIFBLK-1
      IF (NIFBLK.LT.0) THEN
          NIFBLK=0
          IOUT(1)=KABC(3)
          CALL DIAGNO (32)
      ELSE
          ICOL=ICOL-INDENT
          ICOLSV=ICOL
      END IF
      CALL KWDWRT (ITYPE)
      GO TO 500
C
C                        ***** JTYPE = 49
C     ELSE
 1370 IF (NIFBLK.LE.0) THEN
          IOUT(1)=KABC(3)
          CALL DIAGNO (32)
      ELSE
          ICOL=ICOL-INDENT
          ICOLSV=ICOL
      END IF
      CALL KWDWRT (ITYPE)
      GO TO 500
C
C                        ***** JTYPE = 50
C     ENDDO, REPEAT
C       GET CURRENT END-DO NUMBER
 1380 L15=KCNDO(KCNDP)
      KCNDP=KCNDP-1
      IF (KCNDP.LT.0) CALL DIAGNO (43)
      IF (L15.GT.0) THEN
C     CONVERT TO A CONTINUE STATEMENT
C                            PROCESS STATEMENT NUMBER
          IF (NDOS.NE.0) THEN
C                            IT TERMINATES THIS DO-LOOP. INDENT
C                            ONE LESS LEVEL
              IF (L15.EQ.LDOS(NDOS).AND.MLGC.NE.0) THEN
                  ICOL=ICOL-INDENT
                  ICOLSV=ICOL
              END IF
          END IF
C     CONVERT TO A CONTINUE CARD.
          CALL KWDWRT (87)
          IOUT(ICOL)=KERM
          GO TO 490
      ELSE
C     PASS A DO WHILE LOOP TERMINATOR UNALTERED (BUT PROPERLY INDENTED)
          IF (MLGC.NE.0) THEN
              ICOL=ICOL-INDENT
              ICOLSV=ICOL
          END IF
          NIFBLK=NIFBLK-1
C     END DO
          CALL KWDWRT (ITYPE)
          GO TO 500
      END IF
C
C                        ***** JTYPE = 51
C     DO WHILE
 1390 CALL KWDWRT (ITYPE)
C     TREAT UNNUMBERED DO-LOOP THIS WAY IF DESIRED
 1400 CALL COPY (0)
C     GIVE IT A NEGATIVE PSEUDO-STATEMENT NUMBER IN STACK TO PREVENT
C      CONVERSION TO CONTINUE
      KCNDP=KCNDP+1
      KCNDO(KCNDP)=-KENDDO
      KENDDO=KENDDO+1
      NIFBLK=NIFBLK+1
      GO TO 500
C
C
 1410 FORMAT (13X,'***',10I6,'***')
 1420 FORMAT (' JTYPE =',I3,' IS ILLEGAL.  I AM CONFUSED AND CANNOT GO O
     1N.')
      END
      SUBROUTINE PASS2
C
C     THIS ROUTINE READS THE DATA GENERATED BY PASS1 AND WRITES AND
C      PUNCHES THE RENUMBERED DECK.
C     UNNUMBERED CONTINUE AND FORMAT STATEMENTS ARE DELETED WITHOUT
C      A DIAGNOSTIC.
C     UNREACHABLE STATEMENTS ARE DELETED IF *NO CONTINUES
C      IS IN EFFECT (MCONT=0)
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C     SET UP DIMENSIONED ARRAY FOR EFFICIENT PRINTING
      CHARACTER*2 IOUT72(72),MINUS
      EQUIVALENCE (IOUT72(1),IOUT(1)), (MINUS,KSPK(7))
C        TABLE OF EXECUTABLE(1) OR NON-EXECUTABLE(0) BY KLASS
*      INTEGER IEXFLG(12)
C                                     1 1
C         KLASS   0 1 2 3 4 5 6 7 8 9 0 1
*      DATA IEXFLG/0,0,0,0,1,0,1,1,0,1,1,1/
C
      IF (MP2.EQ.0.OR.NRT1.LE.0) RETURN
C
C     MOVE LIST OF NEW STATEMENT NUMBERS FROM TEMP STORAGE
C
      DO 10 I=1,NDEF
           LOCDEF(I)=NEWNUM(I)
 10   CONTINUE
C
C     SET INITIAL CONSTANTS.
C
      IPASS=2
      MPAGE=0
      NREC=0
      NTRAN=0
      IMAX=66*(MAXCNT+1)+6
      JTYPE=0
C
 20   IF (NRT1.EQ.0) GO TO 200
*      JTYPP=JTYPE
      IOLD=IMAX
      CALL RDBIN (SCFIL1,KILI,SERIAL,IOUT)
C                  BLANK OUT REMAINDER OF PREVIOUS CARD, IF NECESSARY.
      IF (IMAX.LT.IOLD) THEN
           INEW=IMAX+1
           DO 30 I=INEW,IOLD
                IOUT(I)=KBL
 30        CONTINUE
      END IF
C                  LOOK FOR $  (FOR WARNING FLAG)
      IF (KLASS.GT.1) THEN
           DO 40 I=7,IMAX
                IF (IOUT(I).EQ.KSPK(10)) THEN
                     IF (MPRIN.EQ.0) WRITE (OUTFIL,240) IOUT72
                     WRITE (OUTFIL,230)
                     GO TO 50
                END IF
 40        CONTINUE
      END IF
C
 50   NRT1=NRT1-1
      IF (NREC.EQ.0) THEN
           CALL HEADER
           IF (MPRIN.NE.0) CALL PAGE (0)
      END IF
C
      IF (MDEB.NE.0) WRITE (OUTFIL,210) KILI,SERIAL
      I=KLASS+1
C      KLASS=0   1   2   3   4   5   6   7   8   9   10  11
      GO TO (20,130,60,130,100,100,100,70,170,130,70,100),I
C                KLASS  DESCRIPTION
C                  0.   CONTROL CARD
C                  1.   COMMENT
C                  2.   HEADER
C                  3.   NO STATEMENT NO ALLOWED (NON-EXECTUABLE)
C                  4.   CONTINUE
C                  5.   FORMAT STATEMENT.
C                  6.   STATEMENT NO. ALLOWED, NO REFERENCES
C                  7.   REFERENCES PRESENT, STATEMENT NO. ALLOWED.
C                  8.   END
C                  9.   INTRODUCTORY
C                  10.  DO
C                  11.  ELSE,ENDIF,ELSEIF, UNRECOGNIZED
C                       (TRANSFER CAN GET HERE REGARDLESS OF LABEL)
C
C     KLASS 0.   CONTROL CARD
C             RESERVED FOR FUTURE DEVELOPMENT.
C
 60   IF (MPRIN.EQ.0) THEN
           CALL PAGE (2)
C          IF (MPUN.NE.0) THEN
C               WRITE (OUTFIL,280) (KIM(I,1),I=1,72)
C          ELSE
C               WRITE (OUTFIL,290) (KIM(I,1),I=1,72)
C          END IF
      END IF
      GO TO 130
C
C     DO REFERENCES.
C
 70   DO 80 I=7,IMAX
           JINT(I)=IOUT(I)
           IOUT(I)=KBL
 80   CONTINUE
      ICOL=6
      JCOL=7
      JMAX=IMAX
      I=1
C
 90   IF (JINT(JCOL).EQ.KLR2) THEN
C     RENUMBER A REFERENCE
           L772=IOUTN(I)
           JCOL=JCOL+1
           I=I+1
           IF (INT(L772).NE.0) THEN
                CALL RENUM
                IF (INT(L772).EQ.0) CALL DIAGNO (53)
           END IF
      ELSE
C     COPY A CHARACTER
           ICOL=ICOL+1
           IOUT(ICOL)=JINT(JCOL)
           JCOL=JCOL+1
      END IF
      IF (JCOL.LE.JMAX) GO TO 90
      IMAX=ICOL
C
C          DO STATEMENT NUMBER
C
 100  L772=L15
      ICOL=0
      IF (INT(L772).NE.0) THEN
           CALL RENUM
c          IF (INT(L772).EQ.0) CALL DIAGNO (53)
      END IF
C        PRINT ALL LABELLED STATEMENTS, ELSE, ELSEIF, ENDIF
      IF (INT(L772).NE.0.OR.KLASS.EQ.11) THEN
C        REMEMBER THAT THIS STATEMENT HAS A PATH TO IT
          NTRAN=0
          GO TO 130
      END IF
C                 DELETE ALL UNLABELLED CONTINUES AND FORMATS
      IF (KLASS.EQ.4.OR.KLASS.EQ.5) THEN
           IF (MDEB.NE.0) WRITE (OUTFIL,220) KLASS
           GO TO 20
      END IF
C           PUNCH IF THERE IS A PATH TO THIS STATEMENT
c     IF (NTRAN.NE.-1) GO TO 130
c                 *CONTINUE MEANS ALL OTHER KLASSES ARE OK
c     IF (MCONT.NE.0) GO TO 130
C                 PUNCH NON-EXECUTABLE STATEMENTS
c     IF (IEXFLG(KLASS+1).EQ.0) GO TO 130
C     ACCEPT GOTO FOLLOWING A COMPUTED GOTO
c     IF (JTYPE.EQ.24 .AND. JTYPP.EQ.23) GO TO 130
C
C
C
C     WRITE  (PUNCH) NEW STATEMENT.
C
 130  IF (KLASS.EQ.1 .AND.MSER.EQ.0) THEN
           N72R=80
      ELSE
           N72R=72
      END IF
C
      CALL KIMPAK
C
      DO 160 J=1,NCD
           NREC=NREC+KD79
C
C     IF NO SERIAL, DO NOT PRINT TRAILING BLANKS.
           IF (MSER.EQ.0) THEN
                N72=N72R
                DO 140 I=N72R,1,-1
                     IF (KIM(I,J).NE.KBL) THEN
                          N72=I
                          GO TO 150
                     END IF
 140            CONTINUE
           END IF
C
C     PRINT AND/OR PUNCH AN OUTPUT RECORD.
C
C     STRING KOL73 IS SET IN SUBROUTINE HEADER.
C
 150       IF (MPRIN.NE.0) THEN
                CALL PAGE (1)
                IF (MSER.LT.0) THEN
                     WRITE (OUTFIL,240) (KIM(I,J),I=1,72),KOL73,NREC
                ELSE IF (MSER.EQ.0) THEN
                     WRITE (OUTFIL,250) (KIM(I,J),I=1,N72)
                ELSE
                     WRITE (OUTFIL,250) (KIM(I,J),I=1,72),SERIAL
                END IF
           END IF
           IF (MPUN.NE.0) THEN
                NPUN=NPUN+1
                IF (MSER.LT.0) THEN
                     WRITE (PUNFIL,260) (KIM(I,J),I=1,72),KOL73,NREC
                ELSE IF (MSER.EQ.0) THEN
                     WRITE (PUNFIL,270) (KIM(I,J),I=1,N72)
                ELSE
                     WRITE (PUNFIL,270) (KIM(I,J),I=1,72),SERIAL
                END IF
           END IF
C
 160  CONTINUE
C           REMEMBER IF THIS IS AN UNCONDITIONAL TRANSFER
      IF (MTRAN.EQ.-1) NTRAN=-1
      GO TO 20
C
C     END STATEMENT.
C
 170  NREC=NREC+KD79
C
C     IF NO SERIAL, DO NOT PRINT TRAILING BLANKS.
      IF (MSER.EQ.0) THEN
           DO 180 I=72,1,-1
                IF (IOUT72(I).NE.KBL) THEN
                     N72=I
                     GO TO 190
                END IF
 180       CONTINUE
      END IF
 190  IF (MPRIN.NE.0) THEN
           CALL PAGE (1)
           IF (MSER.LT.0) THEN
                WRITE (OUTFIL,240) IOUT72,KOL73,NREC,MINUS
           ELSE IF (MSER.EQ.0) THEN
                WRITE (OUTFIL,250) (IOUT72(I),I=1,N72)
           ELSE
                WRITE (OUTFIL,250) IOUT72,SERIAL
           END IF
      END IF
      IF (MPUN.NE.0) THEN
           NPUN=NPUN+1
           IF (MSER.LT.0) THEN
                WRITE (PUNFIL,260) IOUT72,KOL73,NREC,MINUS
           ELSE IF (MSER.EQ.0) THEN
                WRITE (PUNFIL,270) (IOUT72(I),I=1,N72)
           ELSE
                WRITE (PUNFIL,270) IOUT72,SERIAL
           END IF
      END IF
 200  RETURN
C
C
 210  FORMAT (' KLASS',I3,' JTYPE',I3,' L15',I7,' IMAX',I4,' TRAN',I2,'
     1NXRF: ',I4/'  MEX=',I4,' ICOLSV = ',I3,' SERIAL:',8A2)
 220  FORMAT (' DELETING A KLASS=',I3,' STATEMENT')
 230  FORMAT (' ',70X,'$ $ $ $ $')
 240  FORMAT (7X,75A1,I4,A1)
 250  FORMAT (7X,80A1)
 260  FORMAT (75A1,I4,A1)
 270  FORMAT (80A1)
C280  FORMAT (' ',15X,72A1,5X,'--PUNCHED')
C290  FORMAT (' ',15X,72A1,5X,'--NOT PUNCHED')
      END
      SUBROUTINE RDBIN (NUNIT,KV,SER,LIST)
C
C     READ UNFORMATTED INTERMEDIATE FILES.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
      INTEGER KV(8)
      CHARACTER*2 SER(8),LIST(1)
C
C
      READ (NUNIT) KV,SER
      IF (MDEB.NE.0) WRITE (STDERR,10) NUNIT,KV
      CALL REDSTR (NUNIT,LIST,KV(4),IOUTN,KV(6),2)
C                            NORMAL EXIT
      RETURN
C
 10   FORMAT (' Unit ',i2,'  read: ',8I8)
      END
      SUBROUTINE RDIR
C
C     THIS SUBROUTINE GENERATES A REFERENCE DIRECTORY OF STATEMENT
C     NUMBERS SHOWING THE OLD STATEMENT NUMBER, ITS LOCATION IN THE
C     ROUTINE, AND THE NEW STATEMENT NUMBER GENERATED BY TIDY.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C
      DIMENSION INDEX(MXPLBL)
C
      IF (NDEF.LE.0) RETURN
      CALL PAGE (-(8+NDEF))
      CALL PAGE (4)
      WRITE (OUTFIL,60)
C
      DO 10 I=1,NDEF
           INDEX(I)=I
 10   CONTINUE
C
C     ADDRESS-SORT STATEMENT NUMBERS
C
      IF (NDEF.EQ.1) GO TO 40
      M=NDEF+1
 20   NR=0
      M=M-1
      DO 30 I=2,M
           J=INDEX(I-1)
           K=INDEX(I)
           IF (LDEF(J).EQ.LDEF(K)) THEN
                INDEX(I-1)=K
                INDEX(I)=J
                NR=1
           END IF
 30   CONTINUE
      IF (NR.NE.0) GO TO 20
C
C     WRITE  DIRECTORY
C
 40   DO 50 I=1,NDEF
           NW1=NEWNUM(I)
           NO1=LDEF(I)
           LO1=LOCDEF(I)
           J=INDEX(I)
           NW2=NEWNUM(J)
           NO2=LDEF(J)
           LO2=LOCDEF(J)
           CALL PAGE (1)
           WRITE (OUTFIL,70) NW1,NO1,LO1,NO2,LO2,NW2
 50   CONTINUE
C
      CALL PAGE (3)
      WRITE (OUTFIL,80)
C
      RETURN
C
 60   FORMAT (' ',32X,'STATEMENT NUMBER DIRECTORY'/' ',22X,'NEW    OLD
     1 LOC',13X,'OLD   LOC      NEW'/1X)
 70   FORMAT (21X,I5,' = ',I6,',(',I4,').',8X,I6,',(',I4,') = ',I5,'.')
 80   FORMAT (' ',20X,'OLD STATEMENT NUMBERS NOT APPEARING IN THIS DIREC
     1TORY'/21X,'WERE NOT REFERENCED AND HENCE ARE DELETED.')
      END
      SUBROUTINE READER
C     THIS ROUTINE READS CARDS ONE BY ONE, UNTIL IT FINDS A
C     NON-BLANK ONE, THEN RETURNS.   IF IT FINDS AN END-OF-FILE, OR IF
C     IQUIT IS NON-ZERO, IT GENERATES A *STOP CARD.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C
      IF (IQUIT.EQ.0) THEN
10         READ (INFILE,60,END=30) KBUFF
C
C          QUICK CHECK IF THERE IS SOMETHING THERE...
           IF (KBUFF(7).NE.KBL) RETURN
C
C          LOOK FOR A TOTALLY BLANK CARD.
           DO 20 I=1,72
                IF (KBUFF(I).NE.KBL) RETURN
20         CONTINUE
C
C          BLANK CARD. IF INCLUDE FLAG IS SET, MAKE FIRST CHARACTER
C           SPECIAL CODE SO CAN BE RECOGNIZED AS A BLANK COMMENT.
C           OTHERWISE ISSUE MESSAGE AND GET NEXT CARD.
           IF (KBKCOK.EQ.1) THEN
                KBUFF(1)=KBLCMT
                KBUFF(2)=KERM
                RETURN
           ELSE
                CALL PAGE (1)
                WRITE (OUTFIL,70)
                GO TO 10
           END IF
      END IF
C                            NO MORE INPUT
30    IQUIT=1
C
C     MAKE A *STOP CONTROL CARD.
      KBUFF(1)=KSPK(8)
      KBUFF(2)=KABC(19)
      KBUFF(3)=KABC(20)
      KBUFF(4)=KABC(15)
      KBUFF(5)=KABC(16)
C
C     BLANK OUT REST OF LINE.
      DO 40 I=6,72
           KBUFF(I)=KBL
40    CONTINUE
      L15=0
      L25=0
      RETURN
C
60    FORMAT (80A1)
70    FORMAT (35X,'( B L A N K   C A R D )')
      END
      SUBROUTINE REDSTR (LU,LIST,NCHR,IRF,NR,IOP)
      CHARACTER*2 LIST(NCHR)
      DIMENSION IRF(NR)
C
C     WRITE OUT STRING AS SERIES OF 508-(CHAR*2) RECS
C      (APPARENTLY 1024 BYTES IS MAGIC NUMBER FOR PROFORT, AND EACH REC
C       HAS 4-BYTE HEADER AND TRAILER)
C
      DATA MXCHR/508/,MXINT/254/
      NL=1
      MU=MXCHR
 10   NU=min(NCHR,MU)
      NB=NU-NL+1
      CALL IOSTR (LU,LIST(NL),NB,IOP)
      IF (NCHR.GT.NU) THEN
           MU=MU+MXCHR
           NL=NU+1
           GO TO 10
      END IF
C
C     NOW DO THE CROSS-REFERENCE TABLE (253 REFS?~)
      NL=1
      MU=MXINT
 20   NU=min(NR,MU)
      NB=NU-NL+1
      CALL IONUM (LU,IRF(NL),NB,IOP)
      IF (NR.GT.NU) THEN
           MU=MU+MXINT
           NL=NU+1
           GO TO 20
      END IF
C
      RETURN
      END
      SUBROUTINE IOSTR (LU,LIST,NB,IOP)
C
C     READ OR WRITE A STRING
C
      CHARACTER*2 LIST(NB)
      IF (IOP.EQ.1) THEN
           WRITE (LU) LIST
      ELSE
           READ (LU) LIST
      END IF
      RETURN
      END
      SUBROUTINE IONUM (LU,IRF,NR,IOP)
C
C     READ OR WRITE AN INTEGER ARRAY.
C
      DIMENSION IRF(NR)
      IF (IOP.EQ.1) THEN
           WRITE (LU) IRF
      ELSE
           READ (LU) IRF
      END IF
      RETURN
      END
      SUBROUTINE RENUM
C
C     THIS SUBROUTINE INSPECTS THE OLD STATEMENT NUMBER IN L772 AND
C     INSERTS THE NEW NUMBER CORRESPONDING TO L772 IN IOUT STARTING AT
C     ICOL+1.  ON EXIT, L772 CONTAINS THE NEW STATEMENT NUMBER.
C
      INCLUDE 'tidy.inc'
C
C     SEARCH DEFINED STATEMENT TABLE FOR L772.
C
      IF (NDEF.NE.0) THEN
           DO 50 II=1,NDEF
                IF (LDEF(II).EQ.INT(L772)) THEN
C
C     ASSEMBLE NEW STATEMENT NUMBER CHARS IN REVERSE ORDER.
C
                     I=NEWNUM(II)
                     L772=I
                     DO 10 L=1,5
                          IT=I/10
                          K=I-IT*10
                          J=L
                          NTEMP(J)=KDIG(K+1)
                          I=IT
                          IF (I.EQ.0) GO TO 20
10                   CONTINUE
                     J=5
C
C     INSERT STATEMENT NUMBER DIGITS.
C
20                   IF (ICOL.EQ.0) THEN
C                            COLUMNS 1-5
                          DO 30 IK=1,5
                               IOUT(IK)=KBL
30                        CONTINUE
                          IF (MRIT.GE.0) THEN
C                            RIGHT ADJUST TO COLUMN -MRIT
                               ICOL=dim(MRIT,J)
                          ELSE
C                            LEFT ADJUST TO COLUMN MRIT
                               ICOL=min(-MRIT,6-J)
                               ICOL=dim(ICOL,1)
                          END IF
                     END IF
40                   ICOL=ICOL+1
                     IOUT(ICOL)=NTEMP(J)
                     J=J-1
                     IF (J.NE.0) GO TO 40
                     RETURN
                END IF
50         CONTINUE
      END IF
C
C     NOT IN STATEMENT NUMBER LIST. DELETE NUMBER.
C
      L772=0
      RETURN
      END
      SUBROUTINE RLIST
C
C     THIS SUBROUTINE UPDATES THE REFERENCED STATEMENT NUMBER LIST.
C     L772 CONTAINS THE REFERENCED STATEMENT NUMBER.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C
      IF (INT(L772).EQ.0) RETURN
C                  POOR PROGRAMMING PRACTICE.
      IF (INT(L772).EQ.L15) CALL DIAGNO (18)
      IF (NREF.LT.0) RETURN
C
C     CHECK IF IT IS ALREADY REFERENCED.
      IF (NREF.GT.0) THEN
           DO 10 I=1,NREF
                IF (LREF(I).EQ.INT(L772)) RETURN
 10        CONTINUE
      END IF
C
C     ADD REFERENCED STATEMENT TO TABLE.
C
      NREF=NREF+1
      IF (NREF.LE.MXPLBL) THEN
           LREF(NREF)=INT(L772)
      ELSE
C                  TABLE FULL
           CALL DIAGNO (7)
           NREF=-1
           MP2=0
      END IF
C
      IF (MDEB.GT.0) WRITE (OUTFIL,20) int(L772),NREF
C
      RETURN
C
 20   FORMAT (' RLIST added #',I6,' NREF = ',I5)
      END
      SUBROUTINE RSTAT
C
C     THIS SUBROUTINE GETS THE STATEMENT NUMBER REFERENCED AT LOCATION
C     JCOL AND PUTS IT IN L772.  JCOL IS LEFT SET AT THE LOCATION OF THE
C     NEXT SYMBOL ON JINT.
C
      INCLUDE 'tidy.inc'
      L772=0
      IF (JCOL.GT.JMAX) THEN
           JCOL=JMAX
      ELSE
C
           I=JCOL
           DO 20 JCOL=I,JMAX
C              SKIP BLANKS
c>                IF (JINT(JCOL).NE.KBL) THEN
                IF (JINT(JCOL)(1:1).NE.' ') THEN
C                   CHECK IF A DIGIT
                     IP=INDEX('0123456789',JINT(JCOL)(1:1))-1
                     IF (IP.LT.0) RETURN
C                   ADD DIGIT TO NUMBER
                     L772=L772*10+IP
                END IF
 20        CONTINUE
           JCOL=JMAX
           LCPY=KERM
           MEOF=0
      END IF
      RETURN
      END
      SUBROUTINE SKARD
C
C     SUPER-CARD INPUT ROUTINE.
C     THIS ROUTINE READS FORTRAN STATEMENTS WITH UP TO 19 CONTINUATION
C     CARDS AND PACKS THE STATEMENT INTO THE SUPER-CARD --JINT--.
C
      INCLUDE  'tidy.inc'
      INCLUDE  'units.inc'
      LOGICAL RSHFT,LNG72
      CHARACTER*2 KB1,KB6,KZERO,KC,KSTAR,KDOL,KPER,KUPPER,KB1CR1
      EQUIVALENCE (KB1,KBUFF(1)),(KB6,KBUFF(6))
      EQUIVALENCE (KZERO,KDIG(1)),(KC,KABC(3)),(KSTAR,KSPK(8))
      EQUIVALENCE (KDOL,KSPK(14)),(KPER,KSPK(9))
C
      RSHFT=.TRUE.
      KREND=K72
      LNG72=.FALSE.
C
C     TEST FOR A CONTINUATION CARD - SHOULD NOT BE HERE
C      (ANSI F77 ALLOWS EMBEDDED COMMENTS IN CONTINUED STATEMENTS, SO
C       THIS PATCH SHOULD BE REMOVED IF A WAY TO DO THEM IS FOUND)
      IF (KBUFF(1).EQ.KAMPR.OR.(KBUFF(1).EQ.KBL.AND.(KBUFF(6)
     1.NE.KBL.AND.KBUFF(6).NE.KZERO))) THEN
           WRITE (OUTFIL,140)
           CALL DIAGNO (45)
      END IF
C
C     SAVE FIRST CHARACTER OF CARD
      KB1CR1=KUPPER(KBUFF(1))
C
      JMAX=1
      DO 30 I=1,KREND
           IF (I.GT.72.AND.KBUFF(I).NE.KBL) LNG72=.TRUE.
           IF (KBUFF(I).EQ.KTAB) THEN
                IF (I.LT.7.AND.RSHFT) THEN
C                  BLANK REST OF NUMBER FIELD
                     DO 10 L=JMAX,6
                          JINT(L)=KBL
 10                  CONTINUE
                     JMAX=7
                     RSHFT=.FALSE.
C     BLANK THE SERIAL FIELD
                     DO 20 L=1,8
                          SERIAL(L)=KBL
 20                  CONTINUE
C     SET LINE LENGTH TO 80
                     KREND=80
                     GO TO 30
                ELSE
C     TABS PAST COLUMN 6 TRANSLATE TO SPACES WITH F77
                     KBUFF(I)=KBL
                END IF
           END IF
C     IDENTIFY CHARS IN COLS 73-80
           IF (I.GT.72) THEN
                JINT(JMAX)=KBUFF(I)(1:1)//'~'
           ELSE
                JINT(JMAX)=KBUFF(I)
           END IF
           JMAX=JMAX+1
 30   CONTINUE
C
C     TRIM OFF TRAILING BLANKS PAST 72.
 40   IF (JINT(JMAX-1).EQ.' ~') THEN
           JMAX=JMAX-1
           GO TO 40
      END IF
C
C     GRAB EXISTING SERIAL NUMBER IF NEEDED.
      IF (MSER.NE.0.AND.RSHFT) THEN
           DO 50 I=1,8
                SERIAL(I)=KBUFF(I+72)
 50        CONTINUE
      END IF
C
C     SKIP PAGE HEADER IF NOT BEGINNING.
      IF (KOUNT.LE.0) THEN
           CALL HEADER
           IF (MLIST.NE.0) CALL PAGE (0)
      END IF
C
      MEOF=-1
      KOUNT=KOUNT+1
      NREC=NREC+1
      IF (MLIST.NE.0) THEN
           CALL PAGE (1)
           WRITE (OUTFIL,150) NREC,KBUFF
      END IF
C
      NXRF=2
      J=1
C
C     LOOK FOR CONTINUATION CARDS AND TRANSFER THEM TO IOUT VIA KBUFF.
C
      IF (IQUIT.NE.1) THEN
C     IF FIRST CARD WAS A COMMENT, DO NOT TRY TO CONTINUE IT...
           IF (KB1CR1.EQ.KC.OR.KB1CR1.EQ.KBLCMT.OR.KB1CR1.EQ.KSTAR.OR.
     1      KB1CR1.EQ.KDOL.OR.KB1CR1.EQ.KPER) THEN
                LNG72=.FALSE.
                CALL READER
                GO TO 110
           END IF
C
C     NOT COMMENT, CONTINUATIONS ARE LEGAL.
           DO 100 J=2,NQCNTS
                CALL READER
                IF (IQUIT.EQ.1) GO TO 110
C     AMPERSAND MEANS CONTINUATION.
                IF (KB1.EQ.KAMPR) THEN
                     K7=2
                     KREND=80
                     GO TO 70
                ELSE
                     K7=7
                     KREND=K72
                END IF
C     CHECK FOR A TAB IN NUMBER FIELD. IF SO, NOT A CONTINUATION
                DO 60 I=1,6
                     IF (KBUFF(I).EQ.KTAB) GO TO 110
 60             CONTINUE
C     CHECK FOR CONTINUATION OR COMMENTS
C      (DO NOT ALTER COLUMN 1 OF THE CARD.)
                KB1CR1=KUPPER(KB1)
                IF (KB1CR1.EQ.KC) GO TO 110
* new test
                IF (KB1CR1.EQ.'!') GO TO 110
                IF (KB1CR1.EQ.KBLCMT) GO TO 110
                IF (KB1CR1.EQ.KSTAR) GO TO 110
                IF (KB1CR1.EQ.KDOL) GO TO 110
                IF (KB1CR1.EQ.KPER) GO TO 110
                IF (KB6.EQ.KBL) GO TO 110
                IF (KB6.EQ.KZERO) GO TO 110
C
 70             DO 80 I=K7,KREND
                     IF (I.GT.72.AND.KBUFF(I).NE.KBL) LNG72=.TRUE.
                     IF (KBUFF(I).EQ.KTAB) KBUFF(I)=KBL
                     IF (I.GT.72) THEN
                          JINT(JMAX)=KBUFF(I)(1:1)//'~'
                     ELSE
                          JINT(JMAX)=KBUFF(I)
                     END IF
                     JMAX=JMAX+1
 80             CONTINUE
C
C               TRIM OFF TRAILING BLANKS IN 73-80
 90             IF (JINT(JMAX-1).EQ.' ~') THEN
                     JMAX=JMAX-1
                     GO TO 90
                END IF
C
                IF (MLIST.EQ.0) GO TO 100
                CALL PAGE (1)
                WRITE (OUTFIL,160) KBUFF
 100       CONTINUE
C
C     MAYBE TOO MANY CONTINUATION CARDS - WILL TRAP THIS IN KIMPAK
C      AFTER EXHAUSTING ALL REMEDIES...
C
           J=NQCNTS
           CALL READER
      END IF
C
C     LOCATE LAST NON-BLANK COLUMN IN CARD AND EXIT.
C
 110  NCD=J-1
      JMAX=JMAX-1
      DO 120 I=JMAX,1,-1
           IF (JINT(I).NE.KBL.AND.JINT(I).NE.' ~') THEN
                JMAX=I
                GO TO 130
           END IF
 120  CONTINUE
      JMAX=1
 130  JINT(JMAX+1)=KERM
      IF (LNG72) CALL DIAGNO (51)
      RETURN
C
C
 140  FORMAT (' FATAL ERROR - STATEMENT BEGINS WITH CONTINUATION LINE.'/
     1'  POSSIBLY COMMENT WITHIN CONTINUED STATEMENT.'/'  TIDY CANNOT PR
     2OCESS THESE ALTHOUGH THEY ARE LEGAL IN FORTRAN-77.')
 150  FORMAT (1X,I4,2X,80A1)
 160  FORMAT (7X,80A1)
      END
      SUBROUTINE USRCON (CFILNM)
C
C     READS A SEPARATE FILE OF TIDY CONTROL CARDS SO USER DOES NOT
C     HAVE TO EDIT THEM INTO SOURCE FILE.
C
      CHARACTER*64 CFILNM
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
C
      WRITE (OUTFIL,30) CFILNM
C
 10   READ (USRFIL,40,END=20) (JINT(I),I=1,75)
      WRITE (OUTFIL,50) (JINT(I),I=1,75)
      IF (JINT(1).NE.KSPK(8)) THEN
           WRITE (OUTFIL,60)
      ELSE
           JMAX=75
           CALL CONTRL
      END IF
      GO TO 10
C
 20   CLOSE (USRFIL,STATUS='KEEP')
      RETURN
C
C
 30   FORMAT ('       ** T I D Y **  SPECIAL CONTROL CARD FILE: ',A)
 40   FORMAT (75A1)
 50   FORMAT (' ',75A1)
 60   FORMAT (' CONTROL CARDS MUST HAVE * IN COLUMN 1.')
      END
      SUBROUTINE WRBIN (NUNIT,KV,SER,LIST)
C
C     WRITE UNFORMATTED INTERMEDIATE FILES.
C
      INCLUDE 'tidy.inc'
      INCLUDE 'units.inc'
      INTEGER KV(8)
      CHARACTER*2 SER(8),LIST(1)
C
C
      WRITE (NUNIT) KV,SER
      IF (MDEB.NE.0) WRITE (STDERR,10) NUNIT,KV
      CALL REDSTR (NUNIT,LIST,KV(4),IOUTN,KV(6),1)
      RETURN
C
 10   FORMAT (' Unit ',i2,' write: ',8I8)
      END
