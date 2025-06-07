 program simulator
	implicit none
	integer,parameter :: yr = 3	! MAXIMUM AGE
	integer,parameter :: po = 3 ! NUMBER OF PORTS
	integer,parameter :: sp = 2 ! NUMBER OF SPECIES
	integer,parameter :: mx = 76 !格子の数(i=j)
	integer,parameter :: nn = 1000 !個体数
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
	double precision	:: pmax=0.,alpha=0.,kr=0.				!krは誤差係数

	double precision	:: r=0.,s=0.
	double precision	:: z=0.
	integer :: seedsize
	integer,allocatable	:: seed(:)

	integer :: d,i,j,k,t,y,yy,pp,mo,tt,n,fi,fj,ffi,ffj,nm !yyは計算が何年目か，ppは魚種の番号,moは月数,nは個体番号,fi,fjは魚のいるマス目
	integer	:: sumday(12)=0.,portloc(po,2)=0.,interval(po)=0.,shipgo(po,2)=0.,shakoprice(12)
	character(len=32) :: fname

!プログラム全体の流れを記述します

	call random_seed(size=seedsize)
	allocate(seed(seedsize))
	do i = 1, seedsize
		call system_clock(count=seed(i)) 
	enddo
	call random_seed(put=seed)

	call initiate     !初期処理のサブルーチン
	write(*,*) 'Hello'
	do yy= 1, yr 	  !魚の最大年級群+1年分計算します
	do d = 1, 365     !yyは何年目か、dは4月1日を起点として何日目か
	call spawn(d)     !魚の加入のサブルーチン
	do t = 1, 48      !タイムステップは0.5時間です
	call prefer(d,t)  !選好強度計算のサブルーチン
	call migrate(d,t,yy)   !移動のサブルーチン
!	call naturaldeath !自然死亡のサブルーチン		!ひとまず
!	call fishing(d)   !漁獲のサブルーチン			!ひとまず
	enddo
	call results(d)   !結果出力のサブルーチン
	enddo
	call aging        !年級群繰り上げのサブルーチン
	enddo
	call terminate    !終了処理のサブルーチン
	
	contains
	subroutine initiate
	
	open(11,file='input/parameter.dat') !パラメータの入力ファイルです
	do i = 1, 3
	read(11,*) (effort(i,j),j=1,12)     !各漁港の4月～翌3月の努力量(隻/操業日)
	enddo
	read(11,*) (sumday(i),i=1,12)       !各月末までの4月1日を1日目としたときの累計日数
	do i = 1, 3
	read(11,*) (portloc(i,j),j=1,2)     !各漁港の(i,j)
	enddo
	read(11,*) (interval(i),i=1,3)      !何日間隔で操業するか(割り切れる日に出漁する)
	do i = 1, 2
	read(11,*) (mrate(i,y),y=1,yr)      !魚種ごとに、年級群ごとに、自然死亡率
	enddo
	do pp = 1, 2
	read(11,*) (blparam(pp,j),j=1,3)    !魚種ごとに、体長に関するパラメータ
	enddo
	do pp = 1, 2
	read(11,*) (bwparam(pp,j),j=1,2)    !魚種ごとに、体重に関するパラメータ
	enddo		
	do i = 1, 2
	read(11,*) (sparam(i,j),j=1,2)      !網目選択率に関するパラメータ
	enddo
	
	read(11,*) (sparam2(i),i=1,3)       !漁港ごとの網目のサイズ(cm)
	
	read(11,*) inipop           !魚種ごとに、加入時の個体数の初期値(千尾) !read(11,*) (inipop(pp)) でエラーが出たため書き換えた
	
	
	do i = 1, 2
	do j = 1, 3
	read(11,*) (gr(i,j,k),k=1,12)       !魚種(i)ごとに、漁港(j)ごとに、月(k)ごとの漁具能率
	enddo
	enddo	
	read(11,*) (op(j),j=1,2)            !魚種ごとの補正係数(湾全体では、計算対象漁港の何倍の漁獲がされているか)
	do i = 1, 2
	read(11,*) (redeath(i,j),j=1,12)    !魚種(i)ごと、月(j)ごとの再放流時の死亡率
	enddo
	do pp = 1, sp
	read(11,*) (p1(pp,i),i=1,4)         !魚種(pp)ごとに、下限水温(℃)、下側最適選好水温(℃)、上側最適選好水温(℃)、上限水温(℃)
	enddo
	do pp = 1, sp
	read(11,*) (p2(pp,i),i=1,2)         !魚種(pp)ごとに、下限溶存酸素濃度(mg/L)、十分溶存酸素濃度(mg/L)
	enddo	
	read(11,*) (maxsp(i),i=1,2)         !魚種(i)ごとに、1タイムステップで移動可能な魚の割合
	read(11,*) (shakoprice(i),i=1,12)   !月別のシャコの値段(円/kg)
	close(11)	
	
!	READ BENTHOS DATA
	
	open(22,file='input/benthos1.dat')
	do i = 1, mx
	read(22,*) (para(3,i,j), j=1,mx) !アナゴ用の底質選好データ
	enddo
	close(22)

	open(23,file='input/benthos2.dat')
	do i = 1, mx
	read(23,*) (para(4,i,j), j=1,mx) !シャコ用の底質選好データ
	enddo
	close(23)
	
	open(24,file='input/iniplace1.dat')
	do i = 1, mx
	read(24,*) (iniplace(1,i,j), j=1,mx) !アナゴの初期分布域
	enddo
	close(24)	

	open(25,file='input/iniplace2.dat')
	do i = 1, mx
	read(25,*) (iniplace(2,i,j), j=1,mx) !シャコの初期分布域
	enddo
	close(25)			
	
	open(10,file='input/chikei1.dat')
	do i = 1, mx
	read(10,*) (ground(i,j), j=1,mx)     !地形データ(1が計算対象の海,0が陸)
	enddo
	close(10)

	open(12,file='input/chikei2.dat')
	do i = 1, mx
	read(12,*) (ground2(i,j), j=1,mx)    !地形データ(1が海,0が陸)
	enddo
	close(12)
	
	open(80,file='input/farea1.dat')
	open(81,file='input/farea2.dat')	
	open(82,file='input/farea3.dat')	
	do i = 1, mx
	read(80,*) (catchable(1,i,j), j=1,mx) !鈴鹿の操業区域(1が操業可能)
	read(81,*) (catchable(2,i,j), j=1,mx) !有滝の操業区域(1が操業可能)
	read(82,*) (catchable(3,i,j), j=1,mx) !豊浜の操業区域(1が操業可能)
	enddo
	close(82)
	close(81)
	close(80)
	
	open(83,file='input/border.dat') 
	do i = 1, mx
	read(83,*) (border(i,j),j=1,mx) !アナゴが湾外移出する境界
	enddo
	close(83)
	
	open(51,file='output/fishing.csv')      !以下出力用ファイル、漁獲量
	open(52,file='output/biology.csv')	    !資源量
	open(53,file='output/deathcause1.csv')  !アナゴ死因
	open(54,file='output/deathcause2.csv')	!シャコ死因
	open(56,file='output/earnmoney.csv')	!漁獲金額
	open(74,file='output/migration.csv')	!個体移動

	end subroutine

	subroutine spawn(d) !加入のサブルーチン
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

	
!	if(d == 1) then !アナゴは4月1日に加入
!		sumini(1) = sum(iniplace(1,:,:))
!		do i = 1, mx
!		do j = 1, mx
!		pop(1,1,i,j) = inipop(1) * iniplace(1,i,j) / sumini(1) * 1000.
!		enddo
!		enddo
!		
!	else if(d == sumday(6)+1) then !シャコは10月1日に加入
!		sumini(2) = sum(iniplace(2,:,:))	
!		do i = 1, mx
!		do j = 1, mx
!		pop(2,1,i,j) = inipop(2) * iniplace(2,i,j) / sumini(2) * 1000.
!		enddo
!		enddo		
!	endif
	
	end subroutine

	
	subroutine prefer(d,t) !選好強度のサブルーチン
	integer d
	integer t
		
! READ TEMP,SAL,DO DATA EVERY DAY

	write(fname,60),(d-1)*48+t
	open(60,file=fname)
	do i = 1, 2839 !読み込みデータの行数総数
	read(60,*) (temp(j), j=1,5)
	para(1,mx-(temp(2)-1),temp(1)) = temp(3) !水温データの読み込み
	para(2,mx-(temp(2)-1),temp(1)) = temp(5) !溶存酸素データの読み込み
	enddo
	close(60)
60	format('input/field00/bwq-',i5.5,'.dat')	
		
!   CALUCULATE THE PREFERENCE STRENGTH	
	
!   水温選好強度を計算
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

!   溶存酸素濃度選好強度を計算	
		if(para(2,i,j) < p2(pp,1)) then
			pref2(2,pp,y,i,j) = 0.
		else if (para(2,i,j) < p2(pp,2)) then
			pref2(2,pp,y,i,j) = 0.0 + 1.0 * (para(2,i,j)-p2(pp,1))/(p2(pp,2)-p2(pp,1))
		else
			pref2(2,pp,y,i,j) = 1.0
		endif
		
!   総合選好強度を計算

	if(pp == 1) then
		pref(pp,y,i,j) = ground(i,j) * pref2(1,pp,y,i,j) * para(2,i,j)			!アナゴは、海がどうか×水温選好×溶存酸素数値
	else
		pref(pp,y,i,j) = ground(i,j) * pref2(1,pp,y,i,j) * para(2,i,j) * para(pp+2,i,j)	!シャコは、海がどうか×水温選好×溶存酸素数値×底質選好
	endif

!	if( pp == 1 ) then
!	pref(pp,y,i,j) = ground(i,j) * pref2(1,pp,y,i,j) * pref2(2,pp,y,i,j) !アナゴは、海がどうか×水温選好×溶存酸素選好 
!	else
!	pref(pp,y,i,j) = ground(i,j) * pref2(1,pp,y,i,j) * pref2(2,pp,y,i,j) * (para(pp+2,i,j)) !シャコは、海がどうか×水温選好×溶存酸素選好×底質選好
!	endif


	enddo ! j-loop
	enddo ! i-loop  	
	enddo ! y-loop
	enddo ! pp-loop
	
	
	end subroutine

	subroutine migrate(d,t,yy) !選好移動のサブルーチン
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
			kr = 900. - ( ( 450. / p2(pp,2) ) * pref(pp,y,fi,fj) )		!誤差の上限。これを乱数にかける
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


			!誤差範囲は、選考強度が0から最大値まで。最大値は1×十分溶存酸素濃度
			
			
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
!	neighbor = ground(i-1,j) + ground(i,j-1) + ground(i,j+1) + ground(i+1,j) !あるセルの上下左右に海の格子がいくつあるか
!	do pp = 1, 2
!	speed(pp,i,j) = maxsp(pp) !魚種ごとの魚のスピード
!	do y  = 1, yr
!	
!	psum = pref(pp,y,i-1,j) + pref(pp,y,i,j-1) + (5. - neighbor) * pref(pp,y,i,j) +  pref(pp,y,i,j+1) + pref(pp,y,i+1,j) !あるセルとその前後左右の選好強度の合計
!	if(psum == 0.) then
!		npop(pp,y,i,j)   = npop(pp,y,i,j) + pop(pp,y,i,j) !移動可能な場所の選好強度が全部0ならば、元の場所にとどまる
!	else
!		npop(pp,y,i-1,j) = npop(pp,y,i-1,j) + pop(pp,y,i,j) * (speed(pp,i,j) * pref(pp,y,i-1,j)) / psum !上に移動する魚の数
!		npop(pp,y,i,j-1) = npop(pp,y,i,j-1) + pop(pp,y,i,j) * (speed(pp,i,j) * pref(pp,y,i,j-1)) / psum !左に移動する魚の数
!		npop(pp,y,i,j)   = npop(pp,y,i,j)	+ pop(pp,y,i,j) * (1 - speed(pp,i,j) * (psum - pref(pp,y,i,j) * (5. - neighbor)) / psum ) !とどまる魚の数
!		npop(pp,y,i,j+1) = npop(pp,y,i,j+1) + pop(pp,y,i,j) * (speed(pp,i,j) * pref(pp,y,i,j+1)) / psum !右に移動する魚の数
!		npop(pp,y,i+1,j) = npop(pp,y,i+1,j) + pop(pp,y,i,j) * (speed(pp,i,j) * pref(pp,y,i+1,j)) / psum !下に移動する魚の数
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
!	if(y == 3 .or. (y == 2 .and. month(d) >= 6 ) ) then !アナゴの湾外移動、2年目級群の9月1日以降または3年目級群が対象
!		do i = 1, mx
!		do j = 1, mx
!			if(border(i,j) == 1) then
!				deathcause(1,y,3) = deathcause(1,y,3) + pop(1,y,i,j) !伊勢湾の境界上にいるアナゴのdeathcauseを湾外移出(死因の3番目)にカウント
!				pop(1,y,i,j) =  0.                                   !伊勢湾の境界上にいるアナゴを0にする
!			endif
!		enddo
!		enddo
!	endif
!	enddo
!		
11	format(76(f10.2))			
	
	end subroutine

!	subroutine naturaldeath !自然死亡のサブルーチン
!	do i  = 1, mx
!	do j  = 1, mx
!	do pp = 1, 2
!	do y  = 1, yr
!	deathcause(pp,y,1) = deathcause(pp,y,1) + (1-exp(-1.*mrate(pp,y)/(365.*48.))) * pop(pp,y,i,j) !自然死亡をdeathcause(死因の1番目)にカウント
!	pop(pp,y,i,j) = exp(-1.*mrate(pp,y)/(365.*48.)) * pop(pp,y,i,j)	!自然死亡による減少を反映
!	
!	if(para(2,i,j) <= p2(pp,1)) then !貧酸素水塊の影響による死亡、溶存酸素濃度の下限を下回っている場合
!	deathcause(pp,y,5) = deathcause(pp,y,5) + pop(pp,y,i,j) !溶存酸素による死亡をdeathcause(死因の5番目)にカウント
!	pop(pp,y,i,j) = 0. !この場合、当該セルの個体数は0になる
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
	
!	subroutine fishing(d) !漁獲に関するサブルーチン
!	integer d
!	
!	if( t <= 10 ) then !1日のうちタイムステップが1～10のときに操業する(30分ごとに10回)
!	do y = 1, yr
!	do pp= 1, 2	
!	if(pp == 1) then
!	blen(pp,y) = MAX(1.,blparam(pp,1) * ( 1- exp(blparam(pp,2) * ( (y-1) + d/365. + blparam(pp,3) )))) !アナゴの体長の計算
!	else
!	blen(pp,y) = MAX(1., 14.81 * ( 1- exp( -0.739 * ( (y-1) + d/365. - 0.214 ) + 0.032 * sin(2.*3.141592*(d/365. - 0.141 ))))) !シャコの体長の計算
!	endif
!	bwei(pp,y) = bwparam(pp,1) * ( blen(pp,y) ) ** bwparam(pp,2) !アナゴおよびシャコの体重の計算
!	if( pp == 1) then
!		saleable(pp,y) = ndist((blen(pp,y)-25.)/2.39) !ndistは正規分布の累積割合、平均体長blenのとき、放流サイズ以上(saleable)な個体の割合を計算
!		price(1,pp,y) = 0.3586 * blen(pp,y) ** 2.341 !1年目級群アナゴの価格
!		price(2,pp,y) = 0.1628 * blen(pp,y) ** 2.341 !2年目級群アナゴの価格
!		price(3,pp,y) = 0.1507 * blen(pp,y) ** 2.341 !3年目級群アナゴの価格
!	else
!		saleable(pp,y) = ndist((blen(pp,y)-10.)/0.99) !ndistは正規分布の累積割合、平均体長blenのとき、放流サイズ以上(saleable)な個体の割合を計算
!		price(1,pp,y) = 0. !1年目級群シャコの価格
!		price(2,pp,y) = 0. !2年目級群シャコの価格
!		price(3,pp,y) = shakoprice(month(d)) !3年目級群シャコの価格
!	endif
!	enddo
!	enddo
!
!	do k = 1, 3 !網目選択率の計算
!	do y = 1, yr
!	if( y == 1 .and. month(d) <= 6) then !アナゴは、1年目級群の9月までは漁獲されないと仮定
!		sele(1,k,y) = 0.
!		else
!		sele(1,k,y) = MAX(0., 1 /( 1 + exp(-1 * sparam(1,1) * (blen(1,y) / sparam2(k) - sparam(1,2) )))) !シャコ網目選択率の計算
!	endif
!
!	if( y == 1 .or. ( y == 2 .and. month(d) <= 2)) then !シャコは、2年目級群の5月までは漁獲されないと仮定
!		sele(2,k,y) = 0.
!		else
!		sele(2,k,y) = MAX(0., 1 - exp(sparam(2,1) - sparam(2,2) * (blen(2,y) / sparam2(k)))) !シャコ網目選択率の計算
!	endif
!	
!	enddo
!	enddo
!	
!	do i = 1, mx
!	do j = 1, mx
!	gr2(1,i,j) = min(1.,MAX(0., (para(1,i,j) - 10.4)/4.82 )) !漁獲能率(の季節変動項)、アナゴの漁具能率は水温で変化すると仮定
!	gr2(1,i,j) = 1.00
!	if(d>=221 .and. d<=345) then
!	gr2(2,i,j) = 0.50 !漁獲能率(の季節変動項)、シャコは4月1日から数えて221日目から345日目まで漁具能率が半減すると仮定
!	else
!	gr2(2,i,j) = 1.00
!	endif
!	enddo !j-loop
!	enddo !i-loop
!	
!	tship = 0.
!	
!	do k = 1, 3 !ここから漁獲
!	if( MOD(d,interval(k)) == 0 ) then	!d(4月1日からの日数)が、interval(操業間隔)で割り切れる日に操業する
!	shipl = effort(k,month(d)) !shipl  !shiplは操業を行う残り隻数
!	tship = tship + effort(k,month(d)) !その日の全漁港での操業隻数を数える
!	do while(shipl > 0.) !shiplが0になるまで操業を繰り返す
!	
!	do i = 1, mx
!	do j = 1, mx
!	if(catchable(k,i,j) == 0.) then !catchableではないところの操業価値はマイナス無限
!	pvalue(k,i,j) = -999.
!	else !catchableな場所では
!	do pp = 1,2
!	do y = 1, yr
!	qq2(k,pp,y,i,j)  =  pop(pp,y,i,j) * ( 0.08 * sele(pp,k,y) * gr(pp,k,month(d)) * gr2(pp,i,j)) * saleable(pp,y)  * bwei(pp,y) * price(k,pp,y) !期待漁獲金額を計算,0.08は1曳網で引く面積の1格子面積に対する割合
!	!grは月別漁具能率，gr2は漁具能率の変動項
!	if(yy == yr + 1 .and. t == 1 .and. shipl == effort(k,month(d))) then
!	cpue(k,pp,y,i,j) =  pop(pp,y,i,j) * ( 0.08 * sele(pp,k,y) * gr(pp,k,month(d)) * gr2(pp,i,j)) * saleable(pp,y)  * bwei(pp,y) / 1000. !期待CPUE(kg/haul)を計算
!	endif
!	
!	enddo
!	enddo
!	
!	pvalue(k,i,j) = MAX(0.,sum(qq2(k,:,:,i,j)) * ( 1. - 0.02 *  ((i - portloc(k,1)) ** 2. + (j - portloc(k,2)) ** 2.) ** ( 1. / 2. ))) !漁場価値を計算
!
!	endif
!	enddo
!	enddo
!	shipgo(k,:) = maxloc(pvalue(k,:,:)) !船の行く場所は、pvalueが最大となる場所にする
!	haulmap(k,shipgo(k,1),shipgo(k,2)) = haulmap(k,shipgo(k,1),shipgo(k,2)) + 1. !船の行った場所をカウントする
!
!	do y = 1 , yr
!	do pp= 1 , 2
!	qq(k,pp,y) =  pop(pp,y,shipgo(k,1),shipgo(k,2)) * ( 0.08 * sele(pp,k,y) * gr(pp,k,month(d)) * gr2(pp,shipgo(k,1),shipgo(k,2)) * MIN(1.,shipl) * op(pp) ) !
!	if(yy == yr + 1) then
!	totalcatch(k,pp,y) = totalcatch(k,pp,y) +  qq(k,pp,y) / op(pp)  * saleable(pp,y) * bwei(pp,y) / 1000. !漁獲量を計算
!	earnmoney(k,pp,y) = earnmoney(k,pp,y) + qq(k,pp,y) / op(pp)  * saleable(pp,y) * bwei(pp,y) / 1000. * price(k,pp,y) !漁獲金額を計算
!	endif
!	quantity(k,pp,y) = quantity(k,pp,y) +  qq(k,pp,y) / op(pp)  * saleable(pp,y) * bwei(pp,y) / 1000. !漁獲量を計算
!	quantity(k,pp,yr+1) = quantity(k,pp,yr+1) +  qq(k,pp,y) / op(pp) * (1.-saleable(pp,y)) * bwei(pp,y) / 1000. !放流量は漁獲量のyr+1年目として保存
!	deathcause(pp,y,2) = deathcause(pp,y,2) +  qq(k,pp,y) * saleable(pp,y)  !死因(漁獲による死亡)をカウント
!	deathcause(pp,y,3) = deathcause(pp,y,3) +  qq(k,pp,y) * (1.-saleable(pp,y)) * redeath(pp,month(d)) !死因(再放流による死亡)をカウント
!	pop(pp,y,shipgo(k,1),shipgo(k,2)) = pop(pp,y,shipgo(k,1),shipgo(k,2)) - qq(k,pp,y) * (saleable(pp,y) + (1.-saleable(pp,y)) * redeath(pp,month(d)))	!漁獲と再放流による死亡個体数を資源量に反映させます
!	enddo	
!	enddo		
!	shipl = MAX(shipl - 1,0.)	!残りの船の数を1減らします
!	enddo	
!	endif
!	enddo
!	endif
!	end subroutine
	


	
	subroutine results(d) !結果表示のサブルーチン
	integer d
	
	if( yy == yr ) then !最終年度のみ結果を出力
!	RESOURCE DENSITY	!空間分布を出力

	write(fname,75),d
	open(75,file=fname)      !座標ファイル

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
	if(ground(i,j) > 0.) then !計算領域の海では
		write(40,50) i,j,para(1,i,j),para(2,i,j),pref(1,3,i,j),pref(2,3,i,j),spop(1,3,i,j),&
            rpop(1,3,i,j),haulmap(1,i,j),haulmap(2,i,j),haulmap(3,i,j)
!                    i、j、水温、溶存酸素、アナゴ選好、シャコ選好、3年目選考アナゴ個体数、3年目ランダムアナゴ個体数、鈴鹿船分布、有滝船分布、豊浜船分布

	else if(ground2(i,j) > 0.) then !計算領域ではないけど海の場所では-888
		write(40,50) i,j,-888.,-888.,-888.,-888.,-888.,-888.,-888.,-888.,-888.,-888.	
	else !陸の場所では-999
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
		write(*,60) yy,"[years]",d,"[days]" !10日ごとに画面に進捗状況を表示
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

	subroutine aging !年級群を1年ズラすサブルーチン
	
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

	subroutine terminate !終了処理を行うサブルーチン
	close(51)
	close(52)
	close(53)
	close(54)
	close(74)
	write(*,*) 'Complete!'
	end subroutine
	
	
	double precision function ndist(a) !aに偏差÷標準偏差を入力すると、その値までの正規分布の累積割合を返す関数
	
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
	
	real function month(d) !d(4月1日からの日数)を入力すると、月数(1＝4月、12＝3月)を返す関数
	integer d
	do mo = 1, 12
	if( d <= sumday(mo) ) then
	month = mo
	exit
	endif
	enddo
	end function month

end program simulator