module f_wrapper

use iso_c_binding, only: c_double, c_int, c_char, c_null_char

implicit none

contains

function c_to_f_string(s) result(str)
  use iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  integer i, nchars
  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str)
  str = transfer(s(1:nchars), str)
end function c_to_f_string

subroutine broaden(Ai, X1, Bi, X2, N, resol, wli, fluxi, fluxo)bind(c,name='broaden')
	use iso_c_binding, only: c_double, c_int, c_char, c_null_char
	integer, parameter :: dp=kind(0.d0)
	character(kind=c_char,len=1), intent(in) :: Ai(*)
  character(kind=c_char,len=1), intent(in) :: Bi(*)
  character(len=:), allocatable :: str

	real(c_double), intent(in), value :: X1
  real(c_double), intent(in), value :: X2
	integer(c_int), intent(in), value :: N
	real(c_double), intent(in), value :: resol
	real(c_double), intent(in) :: wli(N)
	real(c_double), intent(in) :: fluxi(N)
	real(c_double), intent(out) :: fluxo(N)

	DIMENSION H(4100000)
	DIMENSION RED(40000),BLUE(40000)
	DIMENSION RED1(40000),BLUE1(40000),RED2(40000),BLUE2(40000)
	EQUIVALENCE (RED(1),RED1(1)),(BLUE(1),BLUE1(1))

  INTEGER :: IWL, I, IWL999, IWL1001, IWLNMU, NH, NPROF
  REAL :: RED,BLUE,H,RED1,RED2,BLUE1,BLUE2, WT1, WT2
  REAL :: FWHM, FWHM1, FWHM2
  REAL :: ratio, Wend, Wbegin, Wcen, vstep
  REAL :: SUM, VMAC
  character(len=10) :: A,B

  A = c_to_f_string(Ai)
  B = c_to_f_string(Bi)

	Wbegin = wli(1)

	ratio=1._dp+1._dp/resol
	Wend=Wbegin*ratio**(N-1)
	Wcen=(Wbegin+Wend)*.5_dp
	vstep=2.99792458e5_dp/resol

	FWHM=-1._dp
	IF(B.EQ.'PM        ')FWHM=X1/WCEN/1000._dp*299792.458_dp
	IF(B.EQ.'KM        ')FWHM=X1
	IF(B.EQ.'RESOLUTION')FWHM=299792.458_dp/X1
	IF(B.EQ.'CM-1      ')FWHM=X1/(1.e7_dp/Wcen)*299792.458_dp
	IF(FWHM.LT.0.)THEN
    print *, B
		print *, 'BAD B INPUT'
		CALL EXIT
	ENDIF
	FWHM1=FWHM
	FWHM2=FWHM

  IF(X2.GT.0.)THEN
  IF(B.EQ.'PM        ')FWHM2=X2/WCEN/1000.*299792.458D0
  IF(B.EQ.'KM        ')FWHM2=X2
  IF(B.EQ.'RESOLUTION')FWHM2=299792.458D0/X2
  IF(B.EQ.'CM-1      ')FWHM2=X2/(1.D7/WCEN)*299792.458D0
  ENDIF

	IF(A.EQ.'MACRO     ')GO TO 10
	IF(A.EQ.'GAUSSIAN  ')GO TO 20
	! IF(A.EQ.'SINX/X    ')GO TO 60
	! IF(A.EQ.'RECT      ')GO TO 30
	! IF(A.EQ.'PROFILE   ')GO TO 40

  print *, A
	print *, 'BAD A INPUT'
	CALL EXIT

! MACROTURBULENT VELOCITY IN KM
!   10 print *, 'Calc VMAC Kernel w/ VMAC = ', X1, 'KM/S'
	 10   VMAC=X1
      DO 11 I=1,40000
      	RED(I)=EXP(-(FLOAT(I-1)*VSTEP/VMAC)**2)
      	IF(RED(I).LT.1.E-5)GO TO 12
   11 CONTINUE
   12 NPROF=I
      RED(1)=RED(1)/2.
      SUM=0.
      DO 13 I=1,NPROF
   13 SUM=SUM+RED(I)
      SUM=SUM*2.
      DO 14 I=1,NPROF
	      RED(I)=RED(I)/SUM
   14 BLUE(I)=RED(I)
      GO TO 50

! GAUSSIAN INSTRUMENTAL PROFILE HALF WIDTH IN KM  FWHM
   ! 20 print *, 'Calc GAUSSIAN Kernel w/ FWHM = ', FWHM, 'KM/S', X1, X2
   20 DO 21 I=1,40000
        RED(I)=EXP(-(FLOAT(I-1)*VSTEP/FWHM*.8325546_dp*2.)**2)
        IF(RED(I).LT.1.D-5)GO TO 22
   21 CONTINUE
   22 NPROF=I
      RED(1)=RED(1)/2.
      SUM=0.
      DO 23 I=1,NPROF
   23 SUM=SUM+RED(I)
      SUM=SUM*2.
      DO 24 I=1,NPROF
        RED(I)=RED(I)/SUM
   24 BLUE(I)=RED(I)
      IF(X2.EQ.0.)GO TO 50
      DO 25 I=1,40000
        RED2(I)=EXP(-(FLOAT(I-1)*VSTEP/FWHM2*.8325546_dp*2.)**2)
        IF(RED2(I).LT.1.D-5)GO TO 26
      25 CONTINUE
   26 NPROF=I
      RED2(1)=RED2(1)/2.
      SUM=0.
      DO 27 I=1,NPROF
      27 SUM=SUM+RED2(I)
      SUM=SUM*2.
      DO 28 I=1,NPROF
      RED2(I)=RED2(I)/SUM
   28 BLUE2(I)=RED2(I)
      GO TO 50

   ! 50 print *, 'Doing the broadening'

   50 NH=(N+39999+39999)
      DO 52 I=1,NH
   52 H(I)=0.

  150 DO 157 IWL=1,N
        IWL1001=IWL+40001
        IWL999=IWL+39999
        DO 153 I=1,NPROF
            H(IWL1001-I)=H(IWL1001-I)+BLUE(I)*fluxi(IWL)
        153 H(IWL999+I)=H(IWL999+I)+RED(I)*fluxi(IWL)
  157 CONTINUE

      IF(X2.EQ.0.)GO TO 160
      ! print *, 'Doing left side'
      DO 358 IWL=1,N
          WT2=FLOAT(IWL-1)/FLOAT(N-1)
          WT1=1.-WT2
      358 H(IWL+40000)=H(IWL+40000)*WT1
      ! print *, 'Doing right side'
      DO 357 IWL=1,N
          WT2=FLOAT(IWL-1)/FLOAT(N-1)
          IWL1001=IWL+40001
          IWL999=IWL+39999
          DO 354 I=1,NPROF
          353 H(IWL1001-I)=H(IWL1001-I)+BLUE2(I)*fluxi(IWL)*WT2
          354 H(IWL999+I)=H(IWL999+I)+RED2(I)*fluxi(IWL)*WT2
      357 CONTINUE

  160 DO 170 IWL=1,N
	      IWLNMU=(IWL+39999)
	      fluxo(IWL)=H(IWLNMU)
  170 CONTINUE



end subroutine broaden

end module