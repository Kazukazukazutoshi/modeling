 program simulator
	implicit none
	integer,parameter :: yr = 3	! MAXIMUM AGE
	integer,parameter :: po = 3 ! NUMBER OF PORTS
	integer,parameter :: sp = 2 ! NUMBER OF SPECIES
	integer,parameter :: mx = 76 !�i�q�̐�(i=j)
	integer,parameter :: nn = 1000 !�̐�
	double precision	:: ground(0:mx+1,0:mx+1)=0.,catchable(po,mx,mx)=0.,pvalue(po,mx,mx)=0.,sumini(sp)=0.
	double precision	:: inipop(sp)=0.,iniplace(sp,mx,mx)=0.,pop(sp,yr,mx,mx)=0.,speed(sp,mx,mx)=0.
	double precision	:: pref(sp,yr,0:mx+1,0:mx+1)=0.,pref2(2,sp,yr,0:mx+1,0:mx+1)=0.,npop(sp,yr,0:mx+1,0:mx+1)=0.,qq(po,sp,yr)=0.
    double precision    :: qq2(po,sp,yr,mx,mx)=0.,quantity(po,sp,yr+1)=0.,psum=0., ptsum=0.
	double precision	:: blen(sp,yr)=0.,blparam(sp,po)=0.,sele(sp,po,yr)=0.,sparam(sp,2)=0.,sparam2(po)=0.,bwei(sp,yr)=0.
    double precision    :: bwparam(sp,2)=0.,price(po,sp,yr)=0.
	double precision	:: mrate(sp,yr)=0.,deathcause(sp,yr,5)=0.,shipl=0.,gr(3,po,12)=0.,gr2(sp,mx,mx)=0.,op(sp)=0.
    double precision	:: redeath(sp,12)=0.,saleable(sp,yr)=0.
	double precision 	:: para(4,mx,mx)=0.,temp(5)=0.,temp2(2)=0.,p(mx,mx)=0.,p1(sp,4)=0.,p2(sp,2)=0.,tship=0.,neighbor=0.
	double precision 	:: effort(po,12)=0.,haulmap(po,mx,mx)=0.,maxsp(2)=0.,border(mx,mx)=0.,cpue(po,sp,yr,mx,mx)=0.
	double precision	:: earnmoney(po,sp,yr)=0.,totalcatch(po,sp,yr)=0.,ground2(mx,mx)=0.

	double precision	:: xfish(sp,yr+1,nn)=0.,yfish(sp,yr+1,nn)=0.,xxfish(sp,yr+1,nn)=0.,yyfish(sp,yr+1,nn)=0.
	double precision	:: spop(sp,yr,mx,mx)=0.,rpop(sp,yr,mx,mx)=0.
	double precision	:: pmax=0.,alpha=0.,kr=0.				!kr�͌덷�W��

	double precision	:: r=0.,s=0.
	double precision	:: z=0.
	integer :: seedsize
	integer,allocatable	:: seed(:)

	integer :: d,i,j,k,t,y,yy,pp,mo,tt,n,fi,fj,ffi,ffj,nm !yy�͌v�Z�����N�ڂ��Cpp�͋���̔ԍ�,mo�͌���,n�͌̔ԍ�,fi,fj�͋��̂���}�X��
	integer	:: sumday(12)=0.,portloc(po,2)=0.,interval(po)=0.,shipgo(po,2)=0.,shakoprice(12)
	character(len=32) :: fname

!�v���O�����S�̗̂�����L�q���܂�

	call random_seed(size=seedsize)
	allocate(seed(seedsize))
	do i = 1, seedsize
		call system_clock(count=seed(i)) 
	enddo
	call random_seed(put=seed)

	call initiate     !���������̃T�u���[�`��
	write(*,*) 'Hello'
	do yy= 1, yr 	  !���̍ő�N���Q+1�N���v�Z���܂�
	do d = 1, 365     !yy�͉��N�ڂ��Ad��4��1�����N�_�Ƃ��ĉ����ڂ�
	call spawn(d)     !���̉����̃T�u���[�`��
	do t = 1, 48      !�^�C���X�e�b�v��0.5���Ԃł�
	call prefer(d,t)  !�I�D���x�v�Z�̃T�u���[�`��
	call migrate(d,t,yy)   !�ړ��̃T�u���[�`��
!	call naturaldeath !���R���S�̃T�u���[�`��		!�ЂƂ܂�
!	call fishing(d)   !���l�̃T�u���[�`��			!�ЂƂ܂�
	enddo
	call results(d)   !���ʏo�͂̃T�u���[�`��
	enddo
	call aging        !�N���Q�J��グ�̃T�u���[�`��
	enddo
	call terminate    !�I�������̃T�u���[�`��
	
	contains
	subroutine initiate
	
	open(11,file='input/parameter.dat') !�p�����[�^�̓��̓t�@�C���ł�
	do i = 1, 3
	read(11,*) (effort(i,j),j=1,12)     !�e���`��4���`��3���̓w�͗�(��/���Ɠ�)
	enddo
	read(11,*) (sumday(i),i=1,12)       !�e�����܂ł�4��1����1���ڂƂ����Ƃ��̗݌v����
	do i = 1, 3
	read(11,*) (portloc(i,j),j=1,2)     !�e���`��(i,j)
	enddo
	read(11,*) (interval(i),i=1,3)      !�����Ԋu�ő��Ƃ��邩(����؂����ɏo������)
	do i = 1, 2
	read(11,*) (mrate(i,y),y=1,yr)      !���킲�ƂɁA�N���Q���ƂɁA���R���S��
	enddo
	do pp = 1, 2
	read(11,*) (blparam(pp,j),j=1,3)    !���킲�ƂɁA�̒��Ɋւ���p�����[�^
	enddo
	do pp = 1, 2
	read(11,*) (bwparam(pp,j),j=1,2)    !���킲�ƂɁA�̏d�Ɋւ���p�����[�^
	enddo		
	do i = 1, 2
	read(11,*) (sparam(i,j),j=1,2)      !�ԖڑI�𗦂Ɋւ���p�����[�^
	enddo
	
	read(11,*) (sparam2(i),i=1,3)       !���`���Ƃ̖Ԗڂ̃T�C�Y(cm)
	
	read(11,*) inipop           !���킲�ƂɁA�������̌̐��̏����l(���) !read(11,*) (inipop(pp)) �ŃG���[���o�����ߏ���������
	
	
	do i = 1, 2
	do j = 1, 3
	read(11,*) (gr(i,j,k),k=1,12)       !����(i)���ƂɁA���`(j)���ƂɁA��(k)���Ƃ̋���\��
	enddo
	enddo	
	read(11,*) (op(j),j=1,2)            !���킲�Ƃ̕␳�W��(�p�S�̂ł́A�v�Z�Ώۋ��`�̉��{�̋��l������Ă��邩)
	do i = 1, 2
	read(11,*) (redeath(i,j),j=1,12)    !����(i)���ƁA��(j)���Ƃ̍ĕ������̎��S��
	enddo
	do pp = 1, sp
	read(11,*) (p1(pp,i),i=1,4)         !����(pp)���ƂɁA��������(��)�A�����œK�I�D����(��)�A�㑤�œK�I�D����(��)�A�������(��)
	enddo
	do pp = 1, sp
	read(11,*) (p2(pp,i),i=1,2)         !����(pp)���ƂɁA�����n���_�f�Z�x(mg/L)�A�\���n���_�f�Z�x(mg/L)
	enddo	
	read(11,*) (maxsp(i),i=1,2)         !����(i)���ƂɁA1�^�C���X�e�b�v�ňړ��\�ȋ��̊���
	read(11,*) (shakoprice(i),i=1,12)   !���ʂ̃V���R�̒l�i(�~/kg)
	close(11)	
	
!	READ BENTHOS DATA
	
	open(22,file='input/benthos1.dat')
	do i = 1, mx
	read(22,*) (para(3,i,j), j=1,mx) !�A�i�S�p�̒ꎿ�I�D�f�[�^
	enddo
	close(22)

	open(23,file='input/benthos2.dat')
	do i = 1, mx
	read(23,*) (para(4,i,j), j=1,mx) !�V���R�p�̒ꎿ�I�D�f�[�^
	enddo
	close(23)
	
	open(24,file='input/iniplace1.dat')
	do i = 1, mx
	read(24,*) (iniplace(1,i,j), j=1,mx) !�A�i�S�̏������z��
	enddo
	close(24)	

	open(25,file='input/iniplace2.dat')
	do i = 1, mx
	read(25,*) (iniplace(2,i,j), j=1,mx) !�V���R�̏������z��
	enddo
	close(25)			
	
	open(10,file='input/chikei1.dat')
	do i = 1, mx
	read(10,*) (ground(i,j), j=1,mx)     !�n�`�f�[�^(1���v�Z�Ώۂ̊C,0����)
	enddo
	close(10)

	open(12,file='input/chikei2.dat')
	do i = 1, mx
	read(12,*) (ground2(i,j), j=1,mx)    !�n�`�f�[�^(1���C,0����)
	enddo
	close(12)
	
	open(80,file='input/farea1.dat')
	open(81,file='input/farea2.dat')	
	open(82,file='input/farea3.dat')	
	do i = 1, mx
	read(80,*) (catchable(1,i,j), j=1,mx) !�鎭�̑��Ƌ��(1�����Ɖ\)
	read(81,*) (catchable(2,i,j), j=1,mx) !�L��̑��Ƌ��(1�����Ɖ\)
	read(82,*) (catchable(3,i,j), j=1,mx) !�L�l�̑��Ƌ��(1�����Ɖ\)
	enddo
	close(82)
	close(81)
	close(80)
	
	open(83,file='input/border.dat') 
	do i = 1, mx
	read(83,*) (border(i,j),j=1,mx) !�A�i�S���p�O�ڏo���鋫�E
	enddo
	close(83)
	
	open(51,file='output/fishing.csv')      !�ȉ��o�͗p�t�@�C���A���l��
	open(52,file='output/biology.csv')	    !������
	open(53,file='output/deathcause1.csv')  !�A�i�S����
	open(54,file='output/deathcause2.csv')	!�V���R����
	open(56,file='output/earnmoney.csv')	!���l���z
	open(74,file='output/migration.csv')	!�̈ړ�

	end subroutine

	subroutine spawn(d) !�����̃T�u���[�`��
	integer d
	
	do n = 1, nn
		if(mod(n,10) == 1) then					!(O58)
			if(d == 1) xfish(1,1,n) = 13050.
			if(d == 1) yfish(1,1,n) = 51750.
			if(d == sumday(6)+1) xfish(2,1,n) = 13050.
			if(d == sumday(6)+1) yfish(2,1,n) = 51750.
		endif
		if(mod(n,10) == 2) then					!(C29)
			if(d == 1) xfish(1,1,n) = 2250.
			if(d == 1) yfish(1,1,n) = 25650.
			if(d == sumday(6)+1) xfish(2,1,n) = 2250.
			if(d == sumday(6)+1) yfish(2,1,n) = 25650.
		endif
		if(mod(n,10) == 3) then					!(X12)
			if(d == 1) xfish(1,1,n) = 21150.
			if(d == 1) yfish(1,1,n) = 10350.
			if(d == sumday(6)+1) xfish(2,1,n) = 21150.
			if(d == sumday(6)+1) yfish(2,1,n) = 10350.
		endif
		if(mod(n,10) == 4) then					!(AO1)
			if(d == 1) xfish(1,1,n) = 36450.
			if(d == 1) yfish(1,1,n) = 450.
			if(d == sumday(6)+1) xfish(2,1,n) = 36450.
			if(d == sumday(6)+1) yfish(2,1,n) = 450.
		endif
		if(mod(n,10) == 5) then					!(BE9)
			if(d == 1) xfish(1,1,n) = 50850.
			if(d == 1) yfish(1,1,n) = 7650.
			if(d == sumday(6)+1) xfish(2,1,n) = 50850.
			if(d == sumday(6)+1) yfish(2,1,n) = 7650.
		endif
		if(mod(n,10) == 6) then					!(AS24)
			if(d == 1) xfish(1,1,n) = 40050.
			if(d == 1) yfish(1,1,n) = 21150.
			if(d == sumday(6)+1) xfish(2,1,n) = 40050.
			if(d == sumday(6)+1) yfish(2,1,n) = 21150.
		endif
		if(mod(n,10) == 7) then					!(AM47)
			if(d == 1) xfish(1,1,n) = 34650.
			if(d == 1) yfish(1,1,n) = 41850.
			if(d == sumday(6)+1) xfish(2,1,n) = 34650.
			if(d == sumday(6)+1) yfish(2,1,n) = 41850.
		endif
		if(mod(n,10) == 8) then					!(AY42)
			if(d == 1) xfish(1,1,n) = 45450.
			if(d == 1) yfish(1,1,n) = 37350.
			if(d == sumday(6)+1) xfish(2,1,n) = 45450.
			if(d == sumday(6)+1) yfish(2,1,n) = 37350.
		endif
		if(mod(n,10) == 9) then					!(BT67)
			if(d == 1) xfish(1,1,n) = 64350.
			if(d == 1) yfish(1,1,n) = 59850.
			if(d == sumday(6)+1) xfish(2,1,n) = 64350.
			if(d == sumday(6)+1) yfish(2,1,n) = 59850.
		endif
		if(mod(n,10) == 0) then					!(AL62)
			if(d == 1) xfish(1,1,n) = 33750.
			if(d == 1) yfish(1,1,n) = 55350.
			if(d == sumday(6)+1) xfish(2,1,n) = 33750.
			if(d == sumday(6)+1) yfish(2,1,n) = 55350.
		endif
	enddo

	
!	if(d == 1) then !�A�i�S��4��1���ɉ���
!		sumini(1) = sum(iniplace(1,:,:))
!		do i = 1, mx
!		do j = 1, mx
!		pop(1,1,i,j) = inipop(1) * iniplace(1,i,j) / sumini(1) * 1000.
!		enddo
!		enddo
!		
!	else if(d == sumday(6)+1) then !�V���R��10��1���ɉ���
!		sumini(2) = sum(iniplace(2,:,:))	
!		do i = 1, mx
!		do j = 1, mx
!		pop(2,1,i,j) = inipop(2) * iniplace(2,i,j) / sumini(2) * 1000.
!		enddo
!		enddo		
!	endif
	
	end subroutine

	
	subroutine prefer(d,t) !�I�D���x�̃T�u���[�`��
	integer d
	integer t
		
! READ TEMP,SAL,DO DATA EVERY DAY

	write(fname,60),(d-1)*48+t
	open(60,file=fname)
	do i = 1, 2839 !�ǂݍ��݃f�[�^�̍s������
	read(60,*) (temp(j), j=1,5)
	para(1,mx-(temp(2)-1),temp(1)) = temp(3) !�����f�[�^�̓ǂݍ���
	para(2,mx-(temp(2)-1),temp(1)) = temp(5) !�n���_�f�f�[�^�̓ǂݍ���
	enddo
	close(60)
60	format('input/field00/bwq-',i5.5,'.dat')	
		
!   CALUCULATE THE PREFERENCE STRENGTH	
	
!   �����I�D���x���v�Z
	do pp= 1, 2
	do y = 1, yr
	do i = 1, mx
	do j = 1, mx	
	
	if(pp == 1) then
		if(para(1,i,j) < p1(pp,1)) then
			pref2(1,pp,y,i,j) = 0.0
		else if (para(1,i,j) < p1(pp,2)) then
			pref2(1,pp,y,i,j) = 0.0 + 1.0 * (para(1,i,j)-p1(pp,1))/(p1(pp,2)-p1(pp,1))
		else if (para(1,i,j) < p1(pp,3)) then
			pref2(1,pp,y,i,j) = 1.0
		else if (para(1,i,j) < p1(pp,4)) then
			pref2(1,pp,y,i,j) = 1.0 - 1.0 * (para(1,i,j)-p1(pp,3))/(p1(pp,4)-p1(pp,3))
		else 
			pref2(1,pp,y,i,j) = 0.0
		endif
		
	else
	pref2(1,pp,y,i,j) = 1.0
	endif

!   �n���_�f�Z�x�I�D���x���v�Z	
		if(para(2,i,j) < p2(pp,1)) then
			pref2(2,pp,y,i,j) = 0.
		else if (para(2,i,j) < p2(pp,2)) then
			pref2(2,pp,y,i,j) = 0.0 + 1.0 * (para(2,i,j)-p2(pp,1))/(p2(pp,2)-p2(pp,1))
		else
			pref2(2,pp,y,i,j) = 1.0
		endif
		
!   �����I�D���x���v�Z

	if(pp == 1) then
		pref(pp,y,i,j) = ground(i,j) * pref2(1,pp,y,i,j) * para(2,i,j)			!�A�i�S�́A�C���ǂ����~�����I�D�~�n���_�f���l
	else
		pref(pp,y,i,j) = ground(i,j) * pref2(1,pp,y,i,j) * para(2,i,j) * para(pp+2,i,j)	!�V���R�́A�C���ǂ����~�����I�D�~�n���_�f���l�~�ꎿ�I�D
	endif

!	if( pp == 1 ) then
!	pref(pp,y,i,j) = ground(i,j) * pref2(1,pp,y,i,j) * pref2(2,pp,y,i,j) !�A�i�S�́A�C���ǂ����~�����I�D�~�n���_�f�I�D 
!	else
!	pref(pp,y,i,j) = ground(i,j) * pref2(1,pp,y,i,j) * pref2(2,pp,y,i,j) * (para(pp+2,i,j)) !�V���R�́A�C���ǂ����~�����I�D�~�n���_�f�I�D�~�ꎿ�I�D
!	endif


	enddo ! j-loop
	enddo ! i-loop  	
	enddo ! y-loop
	enddo ! pp-loop
	
	
	end subroutine

	subroutine migrate(d,t,yy) !�I�D�ړ��̃T�u���[�`��
	integer d
	integer t
	integer yy

	z = 3.14159265358979 / 180.

	do n = 1, nn
	do pp = 1, 1
	do y = 1,yy
   		fi = int( yfish(pp,y,n) / 900 ) + 1
		fj = int( xfish(pp,y,n) / 900 ) + 1


   20 		if(n <= 500) then
			pmax = 0.
			if(pref(pp,y,fi-1,fj) > pmax) pmax = pref(pp,y,fi-1,fj)
			if(pref(pp,y,fi,fj-1) > pmax) pmax = pref(pp,y,fi,fj-1)
			if(pref(pp,y,fi+1,fj) > pmax) pmax = pref(pp,y,fi+1,fj)
			if(pref(pp,y,fi,fj+1) > pmax) pmax = pref(pp,y,fi,fj+1)
			if(pref(pp,y,fi,fj) > pmax) pmax = pref(pp,y,fi,fj)
			
			if(pmax /= pref(pp,y,fi,fj)) then
			
				if(pmax == pref(pp,y,fi+1,fj)) alpha = 0.
				if(pmax == pref(pp,y,fi,fj+1)) alpha = 90.
				if(pmax == pref(pp,y,fi-1,fj)) alpha = 180.
				if(pmax == pref(pp,y,fi,fj-1)) alpha = 270.
				
				call random_number(r)
				call random_number(s)
				xxfish(pp,y,n) = xfish(pp,y,n) + 900. * ( r + 0.5 ) * DSIN( z * ( alpha + ( s - 0.5 ) * 90. ))
				yyfish(pp,y,n) = yfish(pp,y,n) + 900. * ( r + 0.5 ) * DCOS( z * ( alpha + ( s - 0.5 ) * 90. ))
			else
				call random_number(r)
				call random_number(s)
				xxfish(pp,y,n) = xfish(pp,y,n) + r * 450. * DSIN( z * s * 360. )
				yyfish(pp,y,n) = yfish(pp,y,n) + r * 450. * DCOS( z * s * 360. )
			endif

		else
			kr = 900. - ( ( 450. / p2(pp,2) ) * pref(pp,y,fi,fj) )		!�덷�̏���B����𗐐��ɂ�����
			if(kr <= 450.) kr = 450.

			if(kr == 450.) then
				call random_number(r)
				call random_number(s)
				xxfish(pp,y,n) = xfish(pp,y,n) + r * kr * DSIN( z * s * 360. )
				yyfish(pp,y,n) = yfish(pp,y,n) + r * kr * DCOS( z * s * 360. )
			else
				call random_number(r)
				call random_number(s)
				xxfish(pp,y,n) = xfish(pp,y,n) + ( r + 0.5 ) * kr * DSIN( z * s * 360. )
				yyfish(pp,y,n) = yfish(pp,y,n) + ( r + 0.5 ) * kr * DCOS( z * s * 360. )
			endif


			!�덷�͈͂́A�I�l���x��0����ő�l�܂ŁB�ő�l��1�~�\���n���_�f�Z�x
			
			
		endif


		ffi = int( yyfish(pp,y,n) / 900 ) + 1
		ffj = int( xxfish(pp,y,n) / 900 ) + 1
		
!		write(*,*) t,n,alpha,fi,fj
		
		if(ground(ffi,ffj) == 0) goto 20

		xfish(pp,y,n) = xxfish(pp,y,n)
		yfish(pp,y,n) = yyfish(pp,y,n)
		


	enddo
	enddo
	enddo

	
!	do i = 1, mx
!	do j = 1, mx
!	neighbor = ground(i-1,j) + ground(i,j-1) + ground(i,j+1) + ground(i+1,j) !����Z���̏㉺���E�ɊC�̊i�q���������邩
!	do pp = 1, 2
!	speed(pp,i,j) = maxsp(pp) !���킲�Ƃ̋��̃X�s�[�h
!	do y  = 1, yr
!	
!	psum = pref(pp,y,i-1,j) + pref(pp,y,i,j-1) + (5. - neighbor) * pref(pp,y,i,j) +  pref(pp,y,i,j+1) + pref(pp,y,i+1,j) !����Z���Ƃ��̑O�㍶�E�̑I�D���x�̍��v
!	if(psum == 0.) then
!		npop(pp,y,i,j)   = npop(pp,y,i,j) + pop(pp,y,i,j) !�ړ��\�ȏꏊ�̑I�D���x���S��0�Ȃ�΁A���̏ꏊ�ɂƂǂ܂�
!	else
!		npop(pp,y,i-1,j) = npop(pp,y,i-1,j) + pop(pp,y,i,j) * (speed(pp,i,j) * pref(pp,y,i-1,j)) / psum !��Ɉړ����鋛�̐�
!		npop(pp,y,i,j-1) = npop(pp,y,i,j-1) + pop(pp,y,i,j) * (speed(pp,i,j) * pref(pp,y,i,j-1)) / psum !���Ɉړ����鋛�̐�
!		npop(pp,y,i,j)   = npop(pp,y,i,j)	+ pop(pp,y,i,j) * (1 - speed(pp,i,j) * (psum - pref(pp,y,i,j) * (5. - neighbor)) / psum ) !�Ƃǂ܂鋛�̐�
!		npop(pp,y,i,j+1) = npop(pp,y,i,j+1) + pop(pp,y,i,j) * (speed(pp,i,j) * pref(pp,y,i,j+1)) / psum !�E�Ɉړ����鋛�̐�
!		npop(pp,y,i+1,j) = npop(pp,y,i+1,j) + pop(pp,y,i,j) * (speed(pp,i,j) * pref(pp,y,i+1,j)) / psum !���Ɉړ����鋛�̐�
!	endif
!	enddo ! y-loop
!	enddo ! pp-loop
!	enddo ! j-loop
!	enddo ! i-loop
!	
!	do i  = 1, mx
!	do j  = 1, mx
!	do pp = 1, 2
!	do y  = 1, yr
!	pop(pp,y,i,j) = npop(pp,y,i,j)
!	npop(pp,y,i,j) = 0.
!	enddo
!	enddo
!	enddo
!	enddo
!	
!	do y = 1, yr
!	if(y == 3 .or. (y == 2 .and. month(d) >= 6 ) ) then !�A�i�S�̘p�O�ړ��A2�N�ڋ��Q��9��1���ȍ~�܂���3�N�ڋ��Q���Ώ�
!		do i = 1, mx
!		do j = 1, mx
!			if(border(i,j) == 1) then
!				deathcause(1,y,3) = deathcause(1,y,3) + pop(1,y,i,j) !�ɐ��p�̋��E��ɂ���A�i�S��deathcause��p�O�ڏo(������3�Ԗ�)�ɃJ�E���g
!				pop(1,y,i,j) =  0.                                   !�ɐ��p�̋��E��ɂ���A�i�S��0�ɂ���
!			endif
!		enddo
!		enddo
!	endif
!	enddo
!		
11	format(76(f10.2))			
	
	end subroutine

!	subroutine naturaldeath !���R���S�̃T�u���[�`��
!	do i  = 1, mx
!	do j  = 1, mx
!	do pp = 1, 2
!	do y  = 1, yr
!	deathcause(pp,y,1) = deathcause(pp,y,1) + (1-exp(-1.*mrate(pp,y)/(365.*48.))) * pop(pp,y,i,j) !���R���S��deathcause(������1�Ԗ�)�ɃJ�E���g
!	pop(pp,y,i,j) = exp(-1.*mrate(pp,y)/(365.*48.)) * pop(pp,y,i,j)	!���R���S�ɂ�錸���𔽉f
!	
!	if(para(2,i,j) <= p2(pp,1)) then !�n�_�f����̉e���ɂ�鎀�S�A�n���_�f�Z�x�̉�����������Ă���ꍇ
!	deathcause(pp,y,5) = deathcause(pp,y,5) + pop(pp,y,i,j) !�n���_�f�ɂ�鎀�S��deathcause(������5�Ԗ�)�ɃJ�E���g
!	pop(pp,y,i,j) = 0. !���̏ꍇ�A���Y�Z���̌̐���0�ɂȂ�
!	endif
!	
!	enddo
!	enddo
!	enddo
!	enddo
!	
!
!	
!	end subroutine
	
!	subroutine fishing(d) !���l�Ɋւ���T�u���[�`��
!	integer d
!	
!	if( t <= 10 ) then !1���̂����^�C���X�e�b�v��1�`10�̂Ƃ��ɑ��Ƃ���(30�����Ƃ�10��)
!	do y = 1, yr
!	do pp= 1, 2	
!	if(pp == 1) then
!	blen(pp,y) = MAX(1.,blparam(pp,1) * ( 1- exp(blparam(pp,2) * ( (y-1) + d/365. + blparam(pp,3) )))) !�A�i�S�̑̒��̌v�Z
!	else
!	blen(pp,y) = MAX(1., 14.81 * ( 1- exp( -0.739 * ( (y-1) + d/365. - 0.214 ) + 0.032 * sin(2.*3.141592*(d/365. - 0.141 ))))) !�V���R�̑̒��̌v�Z
!	endif
!	bwei(pp,y) = bwparam(pp,1) * ( blen(pp,y) ) ** bwparam(pp,2) !�A�i�S����уV���R�̑̏d�̌v�Z
!	if( pp == 1) then
!		saleable(pp,y) = ndist((blen(pp,y)-25.)/2.39) !ndist�͐��K���z�̗ݐϊ����A���ϑ̒�blen�̂Ƃ��A�����T�C�Y�ȏ�(saleable)�Ȍ̂̊������v�Z
!		price(1,pp,y) = 0.3586 * blen(pp,y) ** 2.341 !1�N�ڋ��Q�A�i�S�̉��i
!		price(2,pp,y) = 0.1628 * blen(pp,y) ** 2.341 !2�N�ڋ��Q�A�i�S�̉��i
!		price(3,pp,y) = 0.1507 * blen(pp,y) ** 2.341 !3�N�ڋ��Q�A�i�S�̉��i
!	else
!		saleable(pp,y) = ndist((blen(pp,y)-10.)/0.99) !ndist�͐��K���z�̗ݐϊ����A���ϑ̒�blen�̂Ƃ��A�����T�C�Y�ȏ�(saleable)�Ȍ̂̊������v�Z
!		price(1,pp,y) = 0. !1�N�ڋ��Q�V���R�̉��i
!		price(2,pp,y) = 0. !2�N�ڋ��Q�V���R�̉��i
!		price(3,pp,y) = shakoprice(month(d)) !3�N�ڋ��Q�V���R�̉��i
!	endif
!	enddo
!	enddo
!
!	do k = 1, 3 !�ԖڑI�𗦂̌v�Z
!	do y = 1, yr
!	if( y == 1 .and. month(d) <= 6) then !�A�i�S�́A1�N�ڋ��Q��9���܂ł͋��l����Ȃ��Ɖ���
!		sele(1,k,y) = 0.
!		else
!		sele(1,k,y) = MAX(0., 1 /( 1 + exp(-1 * sparam(1,1) * (blen(1,y) / sparam2(k) - sparam(1,2) )))) !�V���R�ԖڑI�𗦂̌v�Z
!	endif
!
!	if( y == 1 .or. ( y == 2 .and. month(d) <= 2)) then !�V���R�́A2�N�ڋ��Q��5���܂ł͋��l����Ȃ��Ɖ���
!		sele(2,k,y) = 0.
!		else
!		sele(2,k,y) = MAX(0., 1 - exp(sparam(2,1) - sparam(2,2) * (blen(2,y) / sparam2(k)))) !�V���R�ԖڑI�𗦂̌v�Z
!	endif
!	
!	enddo
!	enddo
!	
!	do i = 1, mx
!	do j = 1, mx
!	gr2(1,i,j) = min(1.,MAX(0., (para(1,i,j) - 10.4)/4.82 )) !���l�\��(�̋G�ߕϓ���)�A�A�i�S�̋���\���͐����ŕω�����Ɖ���
!	gr2(1,i,j) = 1.00
!	if(d>=221 .and. d<=345) then
!	gr2(2,i,j) = 0.50 !���l�\��(�̋G�ߕϓ���)�A�V���R��4��1�����琔����221���ڂ���345���ڂ܂ŋ���\������������Ɖ���
!	else
!	gr2(2,i,j) = 1.00
!	endif
!	enddo !j-loop
!	enddo !i-loop
!	
!	tship = 0.
!	
!	do k = 1, 3 !�������狙�l
!	if( MOD(d,interval(k)) == 0 ) then	!d(4��1������̓���)���Ainterval(���ƊԊu)�Ŋ���؂����ɑ��Ƃ���
!	shipl = effort(k,month(d)) !shipl  !shipl�͑��Ƃ��s���c��ǐ�
!	tship = tship + effort(k,month(d)) !���̓��̑S���`�ł̑��Ɛǐ��𐔂���
!	do while(shipl > 0.) !shipl��0�ɂȂ�܂ő��Ƃ��J��Ԃ�
!	
!	do i = 1, mx
!	do j = 1, mx
!	if(catchable(k,i,j) == 0.) then !catchable�ł͂Ȃ��Ƃ���̑��Ɖ��l�̓}�C�i�X����
!	pvalue(k,i,j) = -999.
!	else !catchable�ȏꏊ�ł�
!	do pp = 1,2
!	do y = 1, yr
!	qq2(k,pp,y,i,j)  =  pop(pp,y,i,j) * ( 0.08 * sele(pp,k,y) * gr(pp,k,month(d)) * gr2(pp,i,j)) * saleable(pp,y)  * bwei(pp,y) * price(k,pp,y) !���ҋ��l���z���v�Z,0.08��1�g�Ԃň����ʐς�1�i�q�ʐςɑ΂��銄��
!	!gr�͌��ʋ���\���Cgr2�͋���\���̕ϓ���
!	if(yy == yr + 1 .and. t == 1 .and. shipl == effort(k,month(d))) then
!	cpue(k,pp,y,i,j) =  pop(pp,y,i,j) * ( 0.08 * sele(pp,k,y) * gr(pp,k,month(d)) * gr2(pp,i,j)) * saleable(pp,y)  * bwei(pp,y) / 1000. !����CPUE(kg/haul)���v�Z
!	endif
!	
!	enddo
!	enddo
!	
!	pvalue(k,i,j) = MAX(0.,sum(qq2(k,:,:,i,j)) * ( 1. - 0.02 *  ((i - portloc(k,1)) ** 2. + (j - portloc(k,2)) ** 2.) ** ( 1. / 2. ))) !���ꉿ�l���v�Z
!
!	endif
!	enddo
!	enddo
!	shipgo(k,:) = maxloc(pvalue(k,:,:)) !�D�̍s���ꏊ�́Apvalue���ő�ƂȂ�ꏊ�ɂ���
!	haulmap(k,shipgo(k,1),shipgo(k,2)) = haulmap(k,shipgo(k,1),shipgo(k,2)) + 1. !�D�̍s�����ꏊ���J�E���g����
!
!	do y = 1 , yr
!	do pp= 1 , 2
!	qq(k,pp,y) =  pop(pp,y,shipgo(k,1),shipgo(k,2)) * ( 0.08 * sele(pp,k,y) * gr(pp,k,month(d)) * gr2(pp,shipgo(k,1),shipgo(k,2)) * MIN(1.,shipl) * op(pp) ) !
!	if(yy == yr + 1) then
!	totalcatch(k,pp,y) = totalcatch(k,pp,y) +  qq(k,pp,y) / op(pp)  * saleable(pp,y) * bwei(pp,y) / 1000. !���l�ʂ��v�Z
!	earnmoney(k,pp,y) = earnmoney(k,pp,y) + qq(k,pp,y) / op(pp)  * saleable(pp,y) * bwei(pp,y) / 1000. * price(k,pp,y) !���l���z���v�Z
!	endif
!	quantity(k,pp,y) = quantity(k,pp,y) +  qq(k,pp,y) / op(pp)  * saleable(pp,y) * bwei(pp,y) / 1000. !���l�ʂ��v�Z
!	quantity(k,pp,yr+1) = quantity(k,pp,yr+1) +  qq(k,pp,y) / op(pp) * (1.-saleable(pp,y)) * bwei(pp,y) / 1000. !�����ʂ͋��l�ʂ�yr+1�N�ڂƂ��ĕۑ�
!	deathcause(pp,y,2) = deathcause(pp,y,2) +  qq(k,pp,y) * saleable(pp,y)  !����(���l�ɂ�鎀�S)���J�E���g
!	deathcause(pp,y,3) = deathcause(pp,y,3) +  qq(k,pp,y) * (1.-saleable(pp,y)) * redeath(pp,month(d)) !����(�ĕ����ɂ�鎀�S)���J�E���g
!	pop(pp,y,shipgo(k,1),shipgo(k,2)) = pop(pp,y,shipgo(k,1),shipgo(k,2)) - qq(k,pp,y) * (saleable(pp,y) + (1.-saleable(pp,y)) * redeath(pp,month(d)))	!���l�ƍĕ����ɂ�鎀�S�̐��������ʂɔ��f�����܂�
!	enddo	
!	enddo		
!	shipl = MAX(shipl - 1,0.)	!�c��̑D�̐���1���炵�܂�
!	enddo	
!	endif
!	enddo
!	endif
!	end subroutine
	


	
	subroutine results(d) !���ʕ\���̃T�u���[�`��
	integer d
	
	if( yy == yr ) then !�ŏI�N�x�̂݌��ʂ��o��
!	RESOURCE DENSITY	!��ԕ��z���o��

	write(fname,75),d
	open(75,file=fname)      !���W�t�@�C��

	write(fname,40),d
	open(40,file=fname)	


		spop(:,:,:,:) = 0
		rpop(:,:,:,:) = 0


	do pp = 1, sp
	do y = 1, yr
	do n = 1, nn
		if(n <= 500) then
			fi = int( yfish(pp,y,n) / 900 ) + 1
			fj = int( xfish(pp,y,n) / 900 ) + 1
			write(75,70) pp,y,n,xfish(pp,y,n),yfish(pp,y,n),fi,fj
			spop(pp,y,fi,fj) = spop(pp,y,fi,fj) + 1
		else
			fi = int( yfish(pp,y,n) / 900 ) + 1
			fj = int( xfish(pp,y,n) / 900 ) + 1
			write(75,70) pp,y,n,xfish(pp,y,n),yfish(pp,y,n),fi,fj
			rpop(pp,y,fi,fj) = rpop(pp,y,fi,fj) + 1
		endif
	enddo
	enddo
	enddo


	do i = 1, mx
	do j = 1, mx
	if(ground(i,j) > 0.) then !�v�Z�̈�̊C�ł�
		write(40,50) i,j,para(1,i,j),para(2,i,j),pref(1,3,i,j),pref(2,3,i,j),spop(1,3,i,j),&
            rpop(1,3,i,j),haulmap(1,i,j),haulmap(2,i,j),haulmap(3,i,j)
!                    i�Aj�A�����A�n���_�f�A�A�i�S�I�D�A�V���R�I�D�A3�N�ڑI�l�A�i�S�̐��A3�N�ڃ����_���A�i�S�̐��A�鎭�D���z�A�L��D���z�A�L�l�D���z

	else if(ground2(i,j) > 0.) then !�v�Z�̈�ł͂Ȃ����ǊC�̏ꏊ�ł�-888
		write(40,50) i,j,-888.,-888.,-888.,-888.,-888.,-888.,-888.,-888.,-888.,-888.	
	else !���̏ꏊ�ł�-999
		write(40,50) i,j,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.
	endif
	enddo
	enddo
	close(40)

	close(75)
	
!	CPUE TRANSITION
!	write(51,51) d,(quantity(1,1,1)+quantity(1,1,2)+quantity(1,1,3)),(quantity(2,1,1)+quantity(2,1,2)+quantity(2,1,3)),(quantity(3,1,1)+quantity(3,1,2)+quantity(3,1,3)), &
!				 (quantity(3,2,1)+quantity(3,2,2)),(quantity(3,2,3)),quantity(3,2,4),&
!				 sum(quantity(:,1,1)),sum(quantity(:,1,2)),sum(quantity(:,1,3))	

!	RESOURCE TRANSITION		 
!	write(52,52) d,sum(pop(1,1,:,:))/1000.,sum(pop(1,2,:,:))/1000.,sum(pop(1,3,:,:))/1000.,0.,sum(pop(2,1,:,:))/1000.,sum(pop(2,2,:,:))/1000.,sum(pop(2,3,:,:))/1000.,0.

!	DEATHCAUSE OF CONGER EEL	
!	write(53,53) d,sum(pop(1,1,:,:))/1000.,deathcause(1,1,1)/1000.,deathcause(1,1,2)/1000.,deathcause(1,1,3)/1000. ,deathcause(1,1,5)/1000., &
!				   sum(pop(1,2,:,:))/1000.,deathcause(1,2,1)/1000.,deathcause(1,2,2)/1000.,deathcause(1,2,3)/1000. ,deathcause(1,2,5)/1000., &
!				   sum(pop(1,3,:,:))/1000.,deathcause(1,3,1)/1000.,deathcause(1,3,2)/1000.,deathcause(1,3,3)/1000. ,deathcause(1,3,5)/1000.

!	DEATHCAUSE OF MANTIS SHRIMP				   
!	write(54,54) d,sum(pop(2,1,:,:))/1000.,deathcause(2,1,1)/1000.,deathcause(2,1,2)/1000.,deathcause(2,1,3)/1000., deathcause(2,1,5)/1000., &
!				   sum(pop(2,2,:,:))/1000.,deathcause(2,2,1)/1000.,deathcause(2,2,2)/1000.,deathcause(2,2,3)/1000., deathcause(2,2,5)/1000., &
!				   sum(pop(2,3,:,:))/1000.,deathcause(2,3,1)/1000.,deathcause(2,3,2)/1000.,deathcause(2,3,3)/1000., deathcause(2,3,5)/1000.

!	EARN MONEY BY FISHING
!	write(56,56) d,sum(earnmoney(1,1,:)),sum(earnmoney(2,1,:)),sum(earnmoney(3,1,:)),sum(earnmoney(3,2,:)), &
!				   sum(totalcatch(1,1,:)),sum(totalcatch(2,1,:)),sum(totalcatch(3,1,:)),sum(totalcatch(3,2,:))


!	MIGRATION

	nm = 10

	fi = int( yfish(1,3,nm) / 900 ) + 1
	fj = int( xfish(1,3,nm) / 900 ) + 1

	write(74,74) nm,d,xfish(1,3,nm),yfish(1,3,nm),fi,fj

				   
!	quantity(:,:,:) = 0.
!	if(d == sumday(month(d))) then
!	haulmap(:,:,:) = 0.

!	endif


	endif
	if( mod(d,10) == 0 ) then
		write(*,60) yy,"[years]",d,"[days]" !10�����Ƃɉ�ʂɐi���󋵂�\��
	endif	

40	format('output/result/',i3.3,'.csv')
50	format(2(I3,','),7(f12.3,','),3(f11.0,','))
51	format(I3,',',8(f11.2,','),f11.2)
52	format(I3,',',7(f11.2,','),f11.2)
53	format(I3,',',14(f10.0,','),f10.0)
54	format(I3,',',14(f10.0,','),f10.0)
56	format(I3,',',50(f10.0,','))
60	format(2(I4,A))
70	format(3(I3,','),2(f11.2,','),2(I3,','))
74	format(2(I3,','),2(f12.3,','),2(I3,','))
75	format('output/coordinate/cdn',i3.3,'.csv')

	end subroutine

	subroutine aging !�N���Q��1�N�Y�����T�u���[�`��
	
	do pp = 1, sp
	do n = 1, nn
		xfish(pp,yr,n) = xfish(pp,yr-1,n)
		yfish(pp,yr,n) = yfish(pp,yr-1,n)
		xfish(pp,yr-1,n) = xfish(pp,yr-2,n)
		yfish(pp,yr-1,n) = yfish(pp,yr-2,n)
		xfish(pp,1,n) = 0.
		yfish(pp,1,n) = 0.
	enddo
	enddo


!	do pp= 1, 2
!	do i = 1,mx
!	do j = 1,mx
!		pop(pp,yr,i,j) = pop(pp,yr-1,i,j)	
!		if(yr >= 3) then
!			do y = 1,yr-2
!			pop(pp,yr-y,i,j) = pop(pp,yr-y-1,i,j)
!			enddo
!			pop(pp,1,i,j) = 0.
!		endif
!	enddo
!	enddo
!	enddo
	
!	quantity(:,:,:) = 0.
!	deathcause(:,:,:) = 0.
!	haulmap(:,:,:) = 0.
	end	subroutine

	subroutine terminate !�I���������s���T�u���[�`��
	close(51)
	close(52)
	close(53)
	close(54)
	close(74)
	write(*,*) 'Complete!'
	end subroutine
	
	
	double precision function ndist(a) !a�ɕ΍����W���΍�����͂���ƁA���̒l�܂ł̐��K���z�̗ݐϊ�����Ԃ��֐�
	
	double precision :: a,x,dt
	double precision, parameter :: pi = 3.141592
	integer :: i,sig
	ndist = 0.5
	
	if (a > 0.) then
		sig = 1
	else 
		sig = -1
	endif

	dt = abs(a)/100.

	do  i = 1, 100

	x = dt * (i - 0.5) 

	ndist = ndist + sig * dt / sqrt(2.*pi) * exp(-0.5 * x ** 2. ) 

	enddo
	end function ndist
	
	real function month(d) !d(4��1������̓���)����͂���ƁA����(1��4���A12��3��)��Ԃ��֐�
	integer d
	do mo = 1, 12
	if( d <= sumday(mo) ) then
	month = mo
	exit
	endif
	enddo
	end function month

end program simulator