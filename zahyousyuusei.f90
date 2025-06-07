program simulator
	implicit none
    integer :: ios
    
	! ========================================
	! ファイルパス設定（汎用性のため最初に定義）
	! ========================================
	! === 入力ファイル ===
	character(len=100), parameter :: TERRAIN_FILE = 'input/0629chikei.txt'
	character(len=100), parameter :: TEMPERATURE_FILE = 'input/bTMP1-00012.dat'
	character(len=100), parameter :: PARAMETER_FILE = 'input/parameter_basic.dat'
	
	! === 出力ディレクトリ・ファイル ===
	character(len=100), parameter :: OUTPUT_DIR = 'zahyousyuuseiresult/'
	character(len=100), parameter :: COORD_FILE_PREFIX = 'min'
	character(len=100), parameter :: FISH_NUMBER_FILE = 'fish_number.csv'
	character(len=100), parameter :: PARAMETER_OUTPUT_FILE = 'parameter.csv'
	
	! ========================================
	! シミュレーションパラメータ
	! ========================================
	integer,parameter :: yr = 1	! 最大年数
	integer,parameter :: po = 3 ! NUMBER OF PORTS
	integer,parameter :: sp = 1 ! 魚種の数
    integer,parameter :: nx = 73, ny = 42 ! グリッドの数 (x方向=73, y方向=42)
	integer,parameter :: nn = 5000 ! 個体数
	integer,parameter :: co = 1 ! 構造物の数
	real(8) :: random1, random2, random3, random4
	
	! 地形データ：ground(x方向, y方向) = ground(nx, ny)
	double precision	:: ground(0:nx+1,0:ny+1)=0.
	
	! 魚の座標：x_fish(x方向), y_fish(y方向)
	double precision	:: x_fish(sp,yr+1,nn)=0., y_fish(sp,yr+1,nn)=0.
	
	! 環境データ：pref(魚種, 年, x方向, y方向)
	double precision    :: pref(sp,yr,0:nx+1,0:ny+1)=0., feed(co,nn)=0., preffeed(co,nn)=0.
	
	! パラメータ配列：para(パラメータ種類, x方向, y方向)
	double precision 	:: para(4,nx,ny)=0., temp(5)=0., p1(sp,3)=0., tp(nn)=0.
	
	double precision	:: random_dir=0., angle1=0., angle2=0., angle3(nn)=0., dx(nn)=0., dy(nn)=0.
	double precision	:: dx_temp(nn)=0., dy_temp(nn)=0., dx_feed(nn)=0., dy_feed(nn)=0., dec(nn)=0.
	double precision	:: fx1(nn)=0., fy1(nn)=0.
	double precision	:: Rk, kijun, r(nn)=0., c=0.
	
	double precision	:: distance(nn,nn)=0.		! 距離
	! 結合・分離・整列
	double precision	:: x_co(nn)=0., y_co(nn)=0., num_co(nn)=0., dx_co(nn)=0., dy_co(nn)=0.
	double precision	:: x_sep(nn)=0., y_sep(nn)=0., num_sep(nn)=0., dx_sep(nn)=0., dy_sep(nn)=0.
	double precision	:: x_ali(nn)=0., y_ali(nn)=0., num_ali(nn)=0., dx_ali(nn)=0., dy_ali(nn)=0.
	double precision	:: random_x(nn)=0., random_y(nn)=0.
	double precision	:: dxx(nn)=0., dyy(nn)=0., dis(nn)=0.		! 整列以外の移動距離/軌跡と構造物の距離
	double precision	:: pattern=0., num_con=0., num_temp=0., num_near(nn)=0., near=0., school=0.
	double precision	:: num1=0., num2=0., num3=0., num4=0., num5=0., num6=0., num7=0., num8=0.
	
	double precision	:: q(co)=0., xcon(co)=0., ycon(co)=0., vel(co)=0., ang(co)=0., W(co)=0., kk(co)=0.
	double precision	:: bunsi(nn,nn)=0., bunbo(nn,nn)=0., sig(nn,nn)=0.
	double precision	:: temp_start=0., temp_max=0., temp_last=0.
	double precision	:: para_co=0., para_sep=0., para_ali=0.	! 結合・分離・整列の移動距離のパラメータ
	double precision	:: z=3.14159265358979/180., pi=3.14159265358979
	double precision	:: Ka=0., deltat=60., size=2000., time=100., bl=0.3	! 拡散係数/タイムステップ/グリッドサイズ/計算時間(h)/マアジの体長
	
	! グリッドインデックス：gi(x方向), gj(y方向)
	integer :: lasttime, v, d, i, j, k, t, n, m, gi(nn), gj(nn)
	integer :: pp=1, y=1, l=1		! 魚は1種類、構造物は1つ
	character(len=32) :: fname
 
 
 
	
	write(*,*) '========================================================'
	write(*,*) '          Fish School Simulation Started'
	write(*,*) '========================================================'
	write(*,'(A,I0,A)') ' Time Step          : ', int(deltat), ' minutes'
	write(*,'(A,I0,A)') ' Simulation Time    : ', int(time), ' hours'
	write(*,'(A,I0,A,I0)') ' Grid Size          : nx=', nx, ', ny=', ny
	write(*,'(A,F6.1,A,F6.1,A)') ' Physical Domain    : x=', nx*size/1000., 'km, y=', ny*size/1000., 'km'
	write(*,'(A,I0)') ' Number of Fish     : ', nn
	write(*,*) ''
	write(*,*) '=== Input Files ==='
	write(*,'(A,A)') ' Terrain File       : ', trim(TERRAIN_FILE)
	write(*,'(A,A)') ' Temperature File   : ', trim(TEMPERATURE_FILE)
	write(*,'(A,A)') ' Parameter File     : ', trim(PARAMETER_FILE)
	write(*,*) ''
	write(*,*) '=== Output Files ==='
	write(*,'(A,A)') ' Output Directory   : ', trim(OUTPUT_DIR)
	write(*,'(A,A,A,A)') ' Coordinate Files   : ', trim(OUTPUT_DIR), trim(COORD_FILE_PREFIX), 'XXXXX.csv'
	write(*,'(A,A,A)') ' Fish Number File   : ', trim(OUTPUT_DIR), trim(FISH_NUMBER_FILE)
	write(*,'(A,A,A)') ' Parameter File     : ', trim(OUTPUT_DIR), trim(PARAMETER_OUTPUT_FILE)
	write(*,*) '========================================================'
	write(*,*) ''
	
! 地形読み込み
	write(*,*) ' > Loading terrain data...'
	open(10,file=TERRAIN_FILE, iostat=ios)
	if (ios == 0) then
		do j = 1, ny  ! y方向（行）
			read(10,*) (ground(i,j), i=1,nx)  ! x方向（列）
		enddo
		close(10)
		write(*,'(A,A)') '   OK: Terrain data loaded: ', trim(TERRAIN_FILE)
	else
		write(*,'(A,A)') '   WARNING: Cannot open terrain file: ', trim(TERRAIN_FILE)
		write(*,*) '   -> Using default values (all ocean)'
		! デフォルト値として全て海域に設定
		do i = 1, nx
			do j = 1, ny
				ground(i,j) = 1.0
			enddo
		enddo
	endif
	
! パラメータ読み込み
	write(*,*) ' > Checking parameter file...'
	open(20,file=PARAMETER_FILE,status="old", action="read", iostat=ios)
	if (ios == 0) then
		close(20)
		write(*,'(A,A)') '   OK: Parameter file found: ', trim(PARAMETER_FILE)
	else
		write(*,'(A,A)') '   WARNING: Parameter file not found: ', trim(PARAMETER_FILE)
		write(*,*) '   -> Using default values'
	endif
	
! 初期位置を計算
	do n = 1, nn
10		continue
		call random_number(random1)
		call random_number(random2)
		x_fish(1,1,n) = size * (nx-1) * random1  ! x方向
		y_fish(1,1,n) = size * (ny-1) * random2  ! y方向
		
		! 物理座標からグリッドインデックスへ変換
		gi(n) = int(x_fish(1,1,n)/size) + 1  ! x方向のグリッドインデックス
		gj(n) = int(y_fish(1,1,n)/size) + 1  ! y方向のグリッドインデックス
		
		if (ground(gi(n), gj(n)) /= 1) goto 10  ! 陸地ならやり直し
	enddo

! 水温選好強度を計算
! 水温データファイルを読み込む
	write(*,*) ' > Loading temperature data...'
	open(60, file=TEMPERATURE_FILE, status="old", action="read", iostat=ios)
	if (ios == 0) then
		do
			read(60,*, iostat=ios) i, j, temp(3), temp(4), temp(5)
			if (ios /= 0) exit
			! 注意：ファイルのi,jがどの軸に対応するか確認が必要
			! ここではi=x方向, j=y方向と仮定
			if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
				para(1, i, j) = temp(3)
			endif
		end do
		close(60)
		write(*,'(A,A)') '   OK: Temperature data loaded: ', trim(TEMPERATURE_FILE)
	else
		write(*,'(A,A)') '   WARNING: Cannot open temperature file: ', trim(TEMPERATURE_FILE)
		write(*,*) '   -> Using default values (20 degrees C)'
		! デフォルト値として一様な水温分布を設定
		do i = 1, nx
			do j = 1, ny
				para(1, i, j) = 20.0  ! デフォルト水温20℃
			enddo
		enddo
	endif

	! 水温読み込み直後のチェック
	write(*,*) ''
	write(*,*) ' === Data Verification ==='
	write(*,'(A)') '   Temperature samples:'
	write(*,'(A,I2,A,I2,A,F6.2,A)') '     - Point (', 10, ',', 10, ') = ', para(1,10,10), ' °C'
	write(*,'(A,I2,A,I2,A,F6.2,A)') '     - Point (', 31, ',', 5,  ') = ', para(1,31,5),  ' °C'
	write(*,'(A,I2,A,I2,A,F6.2,A)') '     - Point (', 50, ',', 20, ') = ', para(1,50,20), ' °C'
	write(*,'(A)') '   Terrain samples:'
	write(*,'(A,I2,A,I2,A,F3.0,A)') '     - Point (', 10, ',', 10, ') = ', ground(10,10), ' (1=sea, 0=land)'
	write(*,'(A,I2,A,I2,A,F3.0,A)') '     - Point (', 31, ',', 5,  ') = ', ground(31,5),  ' (1=sea, 0=land)'
	write(*,'(A,I2,A,I2,A,F3.0,A)') '     - Point (', 50, ',', 20, ') = ', ground(50,20), ' (1=sea, 0=land)'
	write(*,*) ''

! 各グリッドに当てはめ、マダイの生息水温を適用
	temp_start = 12.
	temp_max = 25.3
	temp_last = 30.
    
	do i = 1, nx
		do j = 1, ny
			if(para(1,i,j) < temp_start) then
				pref(pp,y,i,j) = 0.0
			else if (para(1,i,j) < temp_max) then
				pref(pp,y,i,j) = 0.0 + 1.0 * (para(1,i,j)-temp_start)/(temp_max-temp_start)
			else if (para(1,i,j) < temp_last) then
				pref(pp,y,i,j) = 1.0 - 1.0 * (para(1,i,j)-temp_max)/(temp_last-temp_max)
			else 
				pref(pp,y,i,j) = 0.0
			endif
		enddo
	enddo
	
! 繰り返し開始(10分ごと)
	lasttime = nint(time * 60. / deltat )		! 繰り返し回数
	do t = 1, lasttime     ! 1タイムステップで繰り返し
	
! 水温選好強度に則って忌避型で移動(方向：ランダム・距離：選好強度が小さいほど大きくなる)
		Ka=(bl*3*60)**2/2*deltat
		do n = 1, nn
			if(x_fish(1,1,n) >= size*0. .and. x_fish(1,1,n) <= size*DBLE(nx) .and.&
			   y_fish(1,1,n) >= size*0. .and. y_fish(1,1,n) <= size*DBLE(ny) ) then

				call random_number(random_dir)
				angle1 = random_dir*360.
				tp(n) = 1. - pref(pp,y,gi(n),gj(n))
				dx_temp(n) = tp(n) * DCOS( z * angle1 ) * sqrt(2*Ka*deltat)
				dy_temp(n) = tp(n) * DSIN( z * angle1 ) * sqrt(2*Ka*deltat)
			endif
		enddo
	
! 誘引型の移動
! エサの選好強度を計算・移動
		do n = 1, nn
			v = bl*10*60/(deltat**0.1)-bl*3*10/sqrt(2.)	!m/min
			
			dx_feed(n) = 0.
			dy_feed(n) = 0.
			r(n) = 0.
			dec(n) = 0.
			
			do l = 1, co
				q(l) = 1000.		! 濃度の基準値
				W(l) = 50.		! 構造物の幅
				kk(l) = 20. * ((W(l)/1000.)**(1.333))		! 拡散係数
				ang(l) = 0.
				vel(l) = 0.		! 流速
				
				! 構造物の位置（領域の中央）
				xcon(1) = size * nx / 2.  ! x方向の中央
				ycon(1) = size * ny / 2.  ! y方向の中央
				
				if(sqrt((ycon(1)- y_fish(1,1,n))**2 + (xcon(1) - x_fish(1,1,n))**2) > (W(l)/2.)) then
					Rk = sqrt((x_fish(1,1,n)-xcon(1))**2 + (y_fish(1,1,n)-ycon(1))**2)		! 構造物からの距離
					
					feed(l,n) = q(l) * exp(- vel(l) * (Rk - abs(y_fish(1,1,n)-ycon(1))) / (2 * kk(l))) / (4 * pi * kk(l) * Rk)
					kijun =  q(l) * exp(- vel(l) * (W(l)/2) / (2 * kk(l))) / (4 * pi * kk(l) * (W(l)/2) )		! 構造物表面でのエサの濃度
					
					preffeed(l,n) = feed(l,n) / kijun
					angle2 = atan2(ycon(1)- y_fish(1,1,n), xcon(1) - x_fish(1,1,n))
					r(n) = sqrt((ycon(1)- y_fish(1,1,n))**2 + (xcon(1) - x_fish(1,1,n))**2)
					dec(n) = r(n) / (c + r(n))	! 減衰項：卒論では未使用
		
					dx_feed(n) = dx_feed(n) + dec(n) * v * preffeed(l,n) * deltat * DCOS(angle2)
					dy_feed(n) = dy_feed(n) + dec(n) * v * preffeed(l,n) * deltat * DSIN(angle2)
				endif
			enddo
		enddo
	
! 距離行列生成/重心計算/移動
		do n = 1, nn
			dx(n) = dx_temp(n) + dx_feed(n)
			dy(n) = dy_temp(n) + dy_feed(n)
			
! 距離計算
			do m = 1, nn
				if (n /= m) then  ! 自分自身を除く
					distance(n,m) = sqrt((x_fish(1,1,m) - x_fish(1,1,n))**2 + (y_fish(1,1,m) - y_fish(1,1,n))**2)
					if (distance(n,m) > 0) then  ! ゼロ除算を避ける
						bunsi(n,m) = dx(n)*(x_fish(1,1,m) - x_fish(1,1,n))+dy(n)*(y_fish(1,1,m) - y_fish(1,1,n))
						bunbo(n,m) = sqrt(dx(n)**2+dy(n)**2)*distance(n,m)
						if (bunbo(n,m) > 0) then
							sig(n,m) = acos(max(-1.0, min(1.0, bunsi(n,m)/bunbo(n,m))))
						else
							sig(n,m) = 0.0
						endif
					else
						sig(n,m) = 0.0
					endif
				else
					distance(n,m) = 0.0
					sig(n,m) = 0.0
				endif
			enddo
			
! 結合
			num_co(n)=0.
			x_co(n)=0.
			y_co(n)=0.
			dx_co(n)=0.
			dy_co(n)=0.
			
			do m = 1, nn
				if(abs(sig(n,m))<2.79) then
					if(distance(n,m)>=bl .and. distance(n,m)<=(bl*5./tan(1./12.))) then
						x_co(n) = x_co(n) + x_fish(1,1,m)
						y_co(n) = y_co(n) + y_fish(1,1,m)
						num_co(n) = num_co(n) + 1.
					endif
				endif
			enddo
			if(num_co(n) < 1. ) then
				dx_co(n) = 0.
				dy_co(n) = 0.
			else
				dx_co(n) = x_co(n) / num_co(n) - x_fish(1,1,n)
				dy_co(n) = y_co(n) / num_co(n) - y_fish(1,1,n)
			endif
			
! 分離
			num_sep(n)=0.
			x_sep(n)=0.
			y_sep(n)=0.
			dx_sep(n)=0.
			dy_sep(n)=0.
			do m = 1, nn
				if(abs(sig(n,m))<2.79) then
					if(distance(n,m)>=0. .and. distance(n,m)<bl) then
						x_sep(n) = x_sep(n) + x_fish(1,1,m)
						y_sep(n) = y_sep(n) + y_fish(1,1,m)
						num_sep(n) = num_sep(n) + 1.
					endif
				endif
			enddo
			if(num_sep(n) < 1. ) then
				dx_sep(n) = 0.
				dy_sep(n) = 0.
			else
				dx_sep(n) = x_sep(n) / num_sep(n) - x_fish(1,1,n)
				dy_sep(n) = y_sep(n) / num_sep(n) - y_fish(1,1,n)
			endif
			
! 整列
			num_ali(n) = 0.
			x_ali(n) = 0.
			y_ali(n) = 0.
			dx_ali(n)=0.
			dy_ali(n)=0.
			dxx(n)=0.
			dyy(n)=0.
			
			do m = 1, nn
				if(abs(sig(n,m))<2.79) then
					if(distance(n,m)>=0. .and. distance(n,m)<=(bl*5./tan(1./12.))) then
						dxx(m) = dx_temp(m) + dx_feed(m) + dx_co(m) * para_co + dx_sep(m) * para_sep
						dyy(m) = dy_temp(m) + dy_feed(m) + dy_co(m) * para_co + dy_sep(m) * para_sep
						x_ali(n) = x_ali(n) + dxx(m)
						y_ali(n) = y_ali(n) + dyy(m)
						num_ali(n) = num_ali(n) + 1.
					endif
				endif
			enddo
			if(num_ali(n) < 1. ) then
				dx_ali(n)=0.
				dy_ali(n)=0.
			else
				dx_ali(n)=x_ali(n)/num_ali(n)-dxx(n)
				dy_ali(n)=y_ali(n)/num_ali(n)-dyy(n)
			endif
		enddo
	
	! 移動パラメータ初期化
	para_co=1.		! 結合の移動パラメータ
	para_sep=-1.	! 分離の移動パラメータ
	para_ali=1.
	pattern=0.
	W(1) = 50.		! 構造物の幅
	xcon(1) = size * nx / 2.  ! x方向の中央
	ycon(1) = size * ny / 2.  ! y方向の中央
	
	do n = 1, nn
		random_x(n) = 0.
		random_y(n) = 0.
			
! 餌・温度
		dx(n) = dx_temp(n) + dx_feed(n)
		dy(n) = dy_temp(n) + dy_feed(n)
		pattern = 4.

		! 陸地チェック用の一時保存
		fx1(n) = x_fish(1,1,n)
		fy1(n) = y_fish(1,1,n)
		
		x_fish(1,1,n) = x_fish(1,1,n) + dx(n)
		y_fish(1,1,n) = y_fish(1,1,n) + dy(n)
		
! 構造物内部に魚が入らず表面に留める・構造物を通過しない
		if(sqrt((x_fish(1,1,n)-xcon(1))**2+(y_fish(1,1,n)-ycon(1))**2)<W(1)/2) then
			fx1(n) = xcon(1) + DCOS(atan2(y_fish(1,1,n)-ycon(1), x_fish(1,1,n)-xcon(1))) * W(1)/2
			fy1(n) = ycon(1) + DSIN(atan2(y_fish(1,1,n)-ycon(1), x_fish(1,1,n)-xcon(1))) * W(1)/2
			x_fish(1,1,n) = fx1(n)
			y_fish(1,1,n) = fy1(n)
		endif
		
! 周期境界
		if(x_fish(1,1,n) >= size*dble(nx)) then
			x_fish(1,1,n) = x_fish(1,1,n) - size*dble(nx)
		endif
		if(x_fish(1,1,n) <= 0.) then
			x_fish(1,1,n) = x_fish(1,1,n) + size*dble(nx)
		endif
		if(y_fish(1,1,n) >= size*dble(ny)) then
			y_fish(1,1,n) = y_fish(1,1,n) - size*dble(ny)
		endif
		if(y_fish(1,1,n) <= 0.) then
			y_fish(1,1,n) = y_fish(1,1,n) + size*dble(ny)
		endif

		! グリッドインデックスの更新
		gi(n) = int( x_fish(1,1,n) / size ) + 1
		gj(n) = int( y_fish(1,1,n) / size ) + 1
		
		! 境界チェック
		if (gi(n) < 1) gi(n) = 1
		if (gi(n) > nx) gi(n) = nx
		if (gj(n) < 1) gj(n) = 1
		if (gj(n) > ny) gj(n) = ny
		
		! 陸地チェック
		if (ground(gi(n), gj(n)) /= 1) then
			x_fish(1,1,n) = fx1(n)
			y_fish(1,1,n) = fy1(n)
			! インデックスも元に戻す
			gi(n) = int( fx1(n) / size ) + 1
			gj(n) = int( fy1(n) / size ) + 1
		endif
	enddo
		
! 近い距離にある個体数
	school = 0.
	near = 0.
	do n = 1, nn
		num_near(n) = 0.
		do m = 1, nn
			distance(n,m) = sqrt((x_fish(1,1,m) - x_fish(1,1,n))**2 + (y_fish(1,1,m) - y_fish(1,1,n))**2)
			if(distance(n,m)<=(5.*bl/tan(1./12.))) then		! 整列で感知できる範囲内
				num_near(n) = num_near(n) + 1.
			endif
		enddo
		near = near + num_near(n)
	enddo
	school = near/(nn*(nn-1))
		
! 構造物周辺・最適温度帯の個体数
	num_con=0.
	num_temp=0.
	num1=0.
	num2=0.
	num3=0.
	num4=0.
	num5=0.
	num6=0.
	num7=0.
	num8=0.
	do n = 1, nn
		if(sqrt((ycon(1)- y_fish(1,1,n))**2 + (xcon(1) - x_fish(1,1,n))**2) <= 500.) then
			num_con = num_con + 1.
		endif
		if(abs(y_fish(1,1,n)-size* DBLE(ny) *(temp_max-12.)/(30.-12.))<=500.) then
			num_temp = num_temp + 1.
		endif
	enddo
		
	if(t==lasttime) then
		do n = 1, nn
			r(n) = sqrt((ycon(1)- y_fish(1,1,n))**2 + (xcon(1) - x_fish(1,1,n))**2)
			if(r(n)<100.) then
				num1 = num1 + 1.
			endif
			if(r(n)>=100. .and. r(n)<200.) then
				num2 = num2 + 1.
			endif
			if(r(n)>=200. .and. r(n)<300.) then
				num3 = num3 + 1.
			endif
			if(r(n)>=300. .and. r(n)<400.) then
				num4 = num4 + 1.
			endif
			if(r(n)>=400. .and. r(n)<500.) then
				num5 = num5 + 1.
			endif
			if(r(n)>=500. .and. r(n)<600.) then
				num6 = num6 + 1.
			endif
			if(r(n)>=600. .and. r(n)<700.) then
				num7 = num7 + 1.
			endif
			if(r(n)>=700. .and. r(n)<800.) then
				num8 = num8 + 1.
			endif
		enddo
	endif

70	format(1(I3,','),2(f11.2,','),2(I3,','),3(f11.2,','))
74	format(2(I3,','),2(f12.3,','),2(I3,','))
77	format(3(I6),1(f10.8))
			   
! 結果表示
	write(fname,'(A,A,I5.5,A)') trim(OUTPUT_DIR), trim(COORD_FILE_PREFIX), int(t*deltat), '.csv'
	! 出力ディレクトリが存在しない場合は現在のディレクトリに出力
	open(75,file=fname, iostat=ios)
	if (ios /= 0) then
		write(fname,'(A,I5.5,A)') trim(COORD_FILE_PREFIX), int(t*deltat), '.csv'
		open(75,file=fname)
	endif
	do n = 1, nn
		write(75,70) n, x_fish(1,1,n), y_fish(1,1,n), gi(n), gj(n)
	enddo
	close(75)
		
! 魚個体数ファイル出力
	open(unit=76, file=trim(OUTPUT_DIR)//trim(FISH_NUMBER_FILE), status='unknown', position='append', action='write', iostat=ios)
	if (ios /= 0) then
		open(unit=76, file=trim(FISH_NUMBER_FILE), status='unknown', position='append', action='write')
	endif
	write(76,77) int(t*deltat), int(num_con), int(num_temp), school
	close(76)
		
	enddo
	
	! パラメータファイル出力
	open(unit=78, file=trim(OUTPUT_DIR)//trim(PARAMETER_OUTPUT_FILE), status='replace', action='write', iostat=ios)
	if (ios /= 0) then
		open(unit=78, file=trim(PARAMETER_OUTPUT_FILE), status='replace', action='write')
	endif
		if(pattern .EQ. 1.) then
			write(78, '(A)') 'all'
			write(*,*) 'all'
		else if(pattern .EQ. 2.) then
			write(78, '(A)') 'feed + boids'
			write(*,*) 'feed + boids'
		else if(pattern .EQ. 3.) then
			write(78, '(A)') 'temp + boids'
			write(*,*) 'temp + boids'
		else if(pattern .EQ. 4.) then
			write(78, '(A)') 'feed + temp'
			write(*,*) 'feed + temp'
		else if(pattern .EQ. 5.) then
			write(78, '(A)') 'feed'
			write(*,*) 'feed'
		else if(pattern .EQ. 6.) then
			write(78, '(A)') 'temp'
			write(*,*) 'temp'
		endif
		write(78, '(I10)') nn
		write(78, '(F10.3)') size
		write(78, '(F10.3)') deltat
		write(78, '(F10.3)') time
		write(78, '(F10.3)') temp_start
		write(78, '(F10.3)') temp_max
		write(78, '(F10.3)') temp_last
		write(78, '(F10.3)') q(1)
		write(78, '(F10.3)') W(1)
		write(78, '(F10.3)') kk(1)
		write(78, '(F10.3)') xcon(1)
		write(78, '(F10.3)') ycon(1)
		write(78, '(F10.3)') Ka
		write(78, '(I5)') v
		write(78, '(F10.3)') para_co
		write(78, '(F10.3)') para_sep
		write(78, '(F10.3)') para_ali
		write(78, '(F10.3)') bl
		write(78, '(F10.3)') bl
		write(78, '(F10.3)') bl/2./tan(1./12.)
		write(78, '(I3)') 100
		write(78, '(F3.1)') c
		write(78, '(F6.1)') num1
		write(78, '(F6.1)') num2
		write(78, '(F6.1)') num3
		write(78, '(F6.1)') num4
		write(78, '(F6.1)') num5
		write(78, '(F6.1)') num6
		write(78, '(F6.1)') num7
		write(78, '(F6.1)') num8
	
	close(78)

	close(74)
	write(*,*) ''
	write(*,*) '========================================================'
	write(*,*) '           Simulation Completed!'
	write(*,*) '========================================================'
	write(*,*) ' Result file locations:'
	write(*,'(A,A,A,A)') '  - Coordinate data  : ', trim(OUTPUT_DIR), trim(COORD_FILE_PREFIX), 'XXXXX.csv'
	write(*,'(A,A,A)') '  - Fish numbers     : ', trim(OUTPUT_DIR), trim(FISH_NUMBER_FILE)
	write(*,'(A,A,A)') '  - Parameters       : ', trim(OUTPUT_DIR), trim(PARAMETER_OUTPUT_FILE)
	write(*,'(A,F6.1,A,F6.1,A)') ' Final domain       : x=[0,', nx*size/1000., ']km, y=[0,', ny*size/1000., ']km'
	write(*,*) '========================================================'
	
end program simulator