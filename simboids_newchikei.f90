program simulator
	implicit none
	integer,parameter :: yr = 1	! 最大年数
	integer,parameter :: po = 3 ! NUMBER OF PORTS
	integer,parameter :: sp = 1 ! 魚種の数
    integer,parameter :: nx = 42, ny = 73 !格子の数
	!integer,parameter :: mx = 50 !格子の数(i=j)
	integer,parameter :: nn = 5000 !個体数
	integer,parameter :: co = 1 !構造物の数
	real(8) :: random1, random2, random3, random4
	double precision	:: ground(0:nx+1,0:ny+1)=0.
    !double precision	:: ground(0:mx+1,0:mx+1)=0.
	double precision	:: x1fish(sp,yr+1,nn)=0.,y1fish(sp,yr+1,nn)=0.
	double precision    :: pref(sp,yr,0:nx+1,0:ny+1)=0.,feed(co,nn)=0.,preffeed(co,nn)=0.
    !double precision    :: pref(sp,yr,0:mx+1,0:mx+1)=0.,feed(co,nn)=0.,preffeed(co,nn)=0.
	double precision 	:: para(4,nx,ny)=0.,temp(5)=0., p1(sp,3)=0., tp(nn)=0.
	!double precision 	:: para(4,mx,mx)=0.,temp(5)=0., p1(sp,3)=0., tp(nn)=0.
	double precision	:: random_dir=0., angle1=0., angle2=0., angle3(nn)=0., dx(nn)=0., dy(nn)=0.
	double precision	:: dx_temp(nn)=0., dy_temp(nn)=0., dx_feed(nn)=0., dy_feed(nn)=0., dec(nn)=0.
	double precision	:: fx1(nn)=0., fy1(nn)=0.
	double precision	:: Rk, kijun, r(nn)=0., c=0.
	
	double precision	:: distance(nn,nn)=0.		!距離
	!結合・分離・整列
	double precision	:: x_co(nn)=0., y_co(nn)=0., num_co(nn)=0., dx_co(nn)=0., dy_co(nn)=0.
	double precision	:: x_sep(nn)=0., y_sep(nn)=0., num_sep(nn)=0., dx_sep(nn)=0., dy_sep(nn)=0.
	double precision	:: x_ali(nn)=0., y_ali(nn)=0., num_ali(nn)=0., dx_ali(nn)=0., dy_ali(nn)=0.
	double precision	:: random_x(nn)=0., random_y(nn)=0.
	double precision	:: dxx(nn)=0., dyy(nn)=0., dis(nn)=0.		!整列以外の移動距離/軌跡と構造物の距離
	double precision	:: pattern=0., num_con=0., num_temp=0., num_near(nn)=0., near=0., school=0.
	double precision	:: num1=0., num2=0., num3=0., num4=0., num5=0., num6=0., num7=0., num8=0.
	
	double precision	:: q(co)=0.,xcon(co)=0.,ycon(co)=0.,vel(co)=0.,ang(co)=0.,W(co)=0.,kk(co)=0.
	double precision	:: bunsi(nn,nn)=0., bunbo(nn,nn)=0., sig(nn,nn)=0.
	double precision	:: temp_start=0., temp_max=0., temp_last=0.
	double precision	:: para_co=0., para_sep=0., para_ali=0.	!結合・分離・整列の移動距離のパラメータ
	double precision	:: z=3.14159265358979/180., pi=3.14159265358979
	double precision	:: Ka=0., deltat=5., size=100., time=70., bl=0.3	!拡散係数/タイムステップ/グリッドサイズ/計算時間(h)/マアジの体長
	integer :: lasttime, v, d, i, j, k, t, n, m, fi(nn), fj(nn)
	integer :: pp=1, y=1, l=1		!魚は1種類、構造物は1つ
	character(len=32) :: fname
 
 
 
	write(*,*) 'Hello', int(deltat), "分ごと", int(time), "時間"
	
!	地形読み込み
	open(10,file='input/0629chikei.txt')
	do i = 1, nx
    !do i = 1, mx
	read(10,*) (ground(i,j), j=1,ny)     !地形データ(1が計算対象の海,0が陸)
	!read(10,*) (ground(i,j), j=1,mx)     !地形データ(1が計算対象の海,0が陸)
    enddo
	close(10)
	
!	パラメータ読み込み
	open(20,file='input/parameter basic.dat',status="old", action="read")
	
	
!	初期位置を計算
	do n = 1, nn
	    call random_number(random1)
	    call random_number(random2)
	    x1fish(1,1,n) = size + size * (ny-2) * random1	
		y1fish(1,1,n) = size + size * (nx-2) * random2
        !x1fish(1,1,n) = size + size * 48. * random1	
		!y1fish(1,1,n) = size + size * 48. * random2
!	    x1fish(1,1,n) = size*15. + size*2. * random1
!	    y1fish(1,1,n) = size*15. + size*2. * random2
		fi(n) = int( x1fish(1,1,n) / size ) + 1
		fj(n) = int( y1fish(1,1,n) / size ) + 1
!	write(*,*) random1, random2, x1fish(1,1,n), y1fish(1,1,n),fi(n),fj(n)
	enddo


!   水温選好強度を計算
!	水温データの読み込み
    open(60,file='input/grid_temperature.csv')
    do i = 1, 3066	!読み込みデータの行数総数
    read(60,*) (temp(j), j=1,3)
!	write(*,*) i, temp(1)
    para(1,int(temp(1)+1.),int(temp(2)+1.)) = temp(3)	!水温データ
    enddo
!	write(*,*) temp(3)
    close(60)

!	各グリッドに当てはめ
	temp_start = 18.
	temp_max = 21.
	temp_last = 24.
	do i = 1, nx
	do j = 1, ny
    !do i = 1, mx
	!do j = 1, mx
		if(para(1,i,j) <temp_start) then
			pref(pp,y,i,j) = 0.0
		else if (para(1,i,j) < temp_max) then
			pref(pp,y,i,j) = 0.0 + 1.0 * (para(1,i,j)-temp_start)/(temp_max-temp_start)
		else if (para(1,i,j) < temp_last) then
			pref(pp,y,i,j) = 1.0 - 1.0 * (para(1,i,j)-temp_max)/(temp_last-temp_max)
		else 
			pref(pp,y,i,j) = 0.0
		endif
		
!		write(*,*) i,j,para(1,i,j),pref(pp,y,i,j),ground(i,j)
	enddo
	enddo
	
	
	
!	繰り返し開始(10分ごと)
	lasttime = nint(time * 60. / deltat )		!繰り返し回数 70h
	do t = 1, lasttime     !1タイムステップで繰り返し
!	write(*,*) int(t*deltat)
	
	
!	水温選好強度に則って忌避型で移動(方向：ランダム・距離：選好強度が小さいほど大きくなる)
	Ka=(bl*3*60)**2/2*deltat
	do n = 1, nn
	
		if(x1fish(1,1,n)>=size*0. .and. x1fish(1,1,n)<=size*DBLE(ny) .and.&
            y1fish(1,1,n)>=size*0. .and. y1fish(1,1,n)<=size*DBLE(nx) ) then
        !if(x1fish(1,1,n)>=size*0. .and. x1fish(1,1,n)<=size*50. .and. y1fish(1,1,n)>=size*0. .and. y1fish(1,1,n)<=size*50. ) then

			call random_number(random_dir)
			angle1 = random_dir*360.
			tp(n) = 1. - pref(pp,y,fi(n),fj(n))
			dx_temp(n) = tp(n) * DCOS( z * angle1 ) * sqrt(2*Ka*deltat)
			dy_temp(n) = tp(n) * DSIN( z * angle1 ) * sqrt(2*Ka*deltat)
			
!			write(*,*) "temp", n, x1fish(1,1,n), y1fish(1,1,n), dx_temp(n), dy_temp(n), sqrt(2*Ka*deltat)
			
		endif
	enddo
	
	
!   誘引型の移動
!	エサの選好強度を計算・移動
	do n = 1, nn
	v = bl*10*60/(deltat**0.1)-bl*3*10/sqrt(2.)	!m/min
!	v=10.
	
	dx_feed(n) = 0.
	dy_feed(n) = 0.
	r(n) = 0.
	dec(n) = 0.
	
	do l = 1, co
		q(l) = 1000.		!濃度の基準値
		W(l) = 50.		!構造物の幅
		kk(l) = 20. * ((W(l)/1000.)**(1.333))		!拡散係数
		ang(l) = 0.
		vel(l) = 0.		!流速
        xcon(1)= size * ny / 2.
		ycon(1)= size * nx / 2.
		!xcon(1)= size * 25.
		!ycon(1)= size * 25.
		
!		write(*,*) "feed",n, x1fish(1,1,n), y1fish(1,1,n)
		
			if(sqrt((ycon(1)- y1fish(1,1,n))**2 + (xcon(1) - x1fish(1,1,n))**2) > (W(l)/2.)) then
				Rk = sqrt((x1fish(1,1,n)-xcon(1))**2 + (y1fish(1,1,n)-ycon(1))**2)		!構造物からの距離
				
				feed(l,n) = q(l) * exp(- vel(l) * (Rk - abs(y1fish(1,1,n)-ycon(1))) / (2 * kk(l))) / (4 * pi * kk(l) * Rk)
				kijun =  q(l) * exp(- vel(l) * (W(l)/2) / (2 * kk(l))) / (4 * pi * kk(l) * (W(l)/2) )		!構造物表面でのエサの濃度
				
				preffeed(l,n) = feed(l,n) / kijun
				angle2 = atan2(ycon(1)- y1fish(1,1,n), xcon(1) - x1fish(1,1,n))
				r(n) = sqrt((ycon(1)- y1fish(1,1,n))**2 + (xcon(1) - x1fish(1,1,n))**2)
				dec(n) = r(n) / (c + r(n))	!減衰項：卒論では未使用
	
				dx_feed(n) = dx_feed(n) + dec(n) * v * preffeed(l,n) * deltat * DCOS(angle2)
				dy_feed(n) = dy_feed(n) + dec(n) * v * preffeed(l,n) * deltat * DSIN(angle2)
				
!				if(dx_feed(n)>100. .or. dy_feed(n)>100.) then
!				write(*,*) t, n, x1fish(1,1,n), y1fish(1,1,n), v, preffeed(l,n), dx_feed(n), dy_feed(n)
!				endif
			endif
		
	enddo
	enddo
	
	
	
!	距離行列生成/重心計算/移動
	do n = 1, nn
	dx(n) = dx_temp(n) + dx_feed(n)
	dy(n) = dy_temp(n) + dy_feed(n)
	
!		距離計算
		do m = 1, nn
			distance(n,m) = sqrt((x1fish(1,1,m) - x1fish(1,1,n))**2 + (y1fish(1,1,m) - y1fish(1,1,n))**2)
			bunsi(n,m) = dx(n)*(x1fish(1,1,m) - x1fish(1,1,n))+dy(n)*(y1fish(1,1,m) - y1fish(1,1,n))
			bunbo(n,m) = sqrt(dx(n)**2+dy(n)**2)*sqrt((x1fish(1,1,m) - x1fish(1,1,n))**2 + (y1fish(1,1,m) - y1fish(1,1,n))**2)
			sig(n,m) = acos(bunsi(n,m)/bunbo(n,m))
			
!			write(*,*) "distance", n, x1fish(1,1,n), y1fish(1,1,n), m, x1fish(1,1,m), y1fish(1,1,m), distance(n,m)
		enddo
		
		
!		結合
		num_co(n)=0.
		x_co(n)=0.
		y_co(n)=0.
		dx_co(n)=0.
		dy_co(n)=0.
		
		do m = 1, nn
		if(abs(sig(n,m))<2.79) then
			if(distance(n,m)>=bl .and. distance(n,m)<=(bl*5./tan(1./12.))) then
				x_co(n) = x_co(n) + x1fish(1,1,m)
				y_co(n) = y_co(n) + y1fish(1,1,m)
				num_co(n) = num_co(n) + 1.
			endif
		endif
		enddo
		if(num_co(n) < 1. ) then
		  dx_co(n) = 0.
		  dy_co(n) = 0.
!		  write(*,*) n, "num_co=",num_co(n)
		 else
		  dx_co(n) = x_co(n) / num_co(n) - x1fish(1,1,n)
		  dy_co(n) = y_co(n) / num_co(n) - y1fish(1,1,n)
		endif
!		write(*,*) "cen", n,  dx_co(n), dy_co(n), x1fish(1,1,n), y1fish(1,1,n)
		
		
!		分離
		num_sep(n)=0.
		x_sep(n)=0.
		y_sep(n)=0.
		dx_sep(n)=0.
		dy_sep(n)=0.
		do m = 1, nn
		if(abs(sig(n,m))<2.79) then
			if(distance(n,m)>=0. .and. distance(n,m)<bl) then
				x_sep(n) = x_sep(n) + x1fish(1,1,m)
				y_sep(n) = y_sep(n) + y1fish(1,1,m)
				num_sep(n) = num_sep(n) + 1.
			endif
		endif
		enddo
		if(num_sep(n) < 1. ) then
		  dx_sep(n) = 0.
		  dy_sep(n) = 0.
!		  write(*,*) n, "num_sep=",num_sep(n)
		 else
		  dx_sep(n) = x_sep(n) / num_sep(n) - x1fish(1,1,n)
		  dy_sep(n) = y_sep(n) / num_sep(n) - y1fish(1,1,n)
		endif
!		write(*,*) "sep", n,  dx_sep(n), dy_sep(n), x1fish(1,1,n), y1fish(1,1,n)
		
		
!		整列
        num_ali(n) = 0.
		x_ali(n) = 0.
		y_ali(n) = 0.
		dx_ali(nn)=0.
		dy_ali(nn)=0.
		dxx(nn)=0.
		dyy(nn)=0.
		
		do m = 1, nn
		if(abs(sig(n,m))<2.79) then
			if(distance(n,m)>=0. .and. distance(n,m)<=(bl*5./tan(1./12.))) then
			
				dxx(m) = dx_temp(m) + dx_feed(m) + dx_co(m) * para_co + dx_sep(m) * para_sep
				dyy(m) = dy_temp(m) + dy_feed(m) + dy_co(m) * para_co + dy_sep(m) * para_sep
				x_ali(n) = x_ali(n) + dxx(m)
				y_ali(n) = y_ali(n) + dyy(m)
				num_ali(n) = num_ali(n) + 1.
				
!				write(*,*) t, n, x1fish(1,1,m), y1fish(1,1,m)
				
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
!		write(*,*) "ali", n, ang_ali(n), num_ali(n)
	enddo
	
!	移動
	
	para_co=1.		!結合の移動パラメータ
	para_sep=-1.	!分離の移動パラメータ
	para_ali=1.
	pattern=0.
	random_x(n) = 0.
	random_y(n) = 0.
	W(1) = 50.		!構造物の幅
	xcon(1)= size * 25.
	ycon(1)= size * 25.
	
	do n = 1, nn
		
!		write(*,*) t, n, x1fish(1,1,n), y1fish(1,1,n), ang_ali(n), num_ali(n)
		
!		全部
!		dx(n) = dx_temp(n) + dx_feed(n) + dx_co(n) * para_co + dx_sep(n) * para_sep + dx_ali(n) * para_ali
!		dy(n) = dy_temp(n) + dy_feed(n) + dy_co(n) * para_co + dy_sep(n) * para_sep + dy_ali(n) * para_ali
!		pattern = 1.
		
!		餌・ボイド
!		dx(n) = dx_feed(n) + dx_co(n) * para_co + dx_sep(n) * para_sep + dx_ali(n) * para_ali
!		dy(n) = dy_feed(n) + dy_co(n) * para_co + dy_sep(n) * para_sep + dy_ali(n) * para_ali
!		pattern = 2.
		
!		温度・ボイド
!		dx(n) = dx_temp(n) + dx_co(n) * para_co + dx_sep(n) * para_sep + dx_ali(n) * para_ali
!		dy(n) = dy_temp(n) + dy_co(n) * para_co + dy_sep(n) * para_sep + dy_ali(n) * para_ali
!		pattern = 3.

!		餌・温度
		dx(n) = dx_temp(n) + dx_feed(n)
		dy(n) = dy_temp(n) + dy_feed(n)
		pattern = 4.

!		餌
!		dx(n) = dx_feed(n)
!		dy(n) = dy_feed(n)
!		pattern = 5.
		
!		温度
!		dx(n) = dx_temp(n)
!		dy(n) = dy_temp(n)
!		pattern = 6.

		x1fish(1,1,n) = x1fish(1,1,n) + dx(n)
		y1fish(1,1,n) = y1fish(1,1,n) + dy(n)
		
!		write(*,*) x1fish(1,1,n), y1fish(1,1,n), dx(n), dy(n)
		
!		構造物内部に魚が入らず表面に留める・構造物を通過しない

		if(sqrt(x1fish(1,1,n)**2+y1fish(1,1,n)**2)<w(1)/2) then
		
!		if(abs(dx(n))>0. .and. abs(dy(n))>0.) then
!		dis(n) = abs(dy(n)*xcon(1)-dx(n)*ycon(1)-dy(n)*(x1fish(1,1,n)-dx(n))+dx(n)*(y1fish(1,1,n)-dy(n)))/sqrt(dx(n)**2+dy(n)**2)
!		write(*,*) t, n, dis(n), dx(n), dy(n)
		
		
!		if(dis(n) <= (W(1)/2)) then
		
			fx1(n) = xcon(1) + DCOS(atan2(ycon(1)- (y1fish(1,1,n)+dy(n)), xcon(1) - (x1fish(1,1,n)+dx(n)))) * W(1)/2
			fy1(n) = ycon(1) + DSIN(atan2(ycon(1)- (y1fish(1,1,n)+dy(n)), xcon(1) - (x1fish(1,1,n)+dx(n)))) * W(1)/2
			x1fish(1,1,n) = fx1(n)
			y1fish(1,1,n) = fy1(n)
			
!			write(*,*) t, n, xcon(1), ycon(1), x1fish(1,1,n), y1fish(1,1,n)
			
		endif
!		endif
		
!		周期境界
		if(x1fish(1,1,n)>=size*dble(ny)) then
			x1fish(1,1,n) = x1fish(1,1,n) - size*dble(ny)
		endif
		if(x1fish(1,1,n)<=0.) then
			x1fish(1,1,n) = x1fish(1,1,n) + size*dble(ny)
		endif
		if(y1fish(1,1,n)>=size*dble(nx)) then
			y1fish(1,1,n) = y1fish(1,1,n) - size*dble(nx)
		endif
		if(y1fish(1,1,n)<=0.) then
			y1fish(1,1,n) = y1fish(1,1,n) + size*dble(nx)
        endif
		!
  !      if(x1fish(1,1,n)>=size*dble(mx)) then
		!	x1fish(1,1,n) = x1fish(1,1,n) - size*dble(mx)
		!endif
		!if(x1fish(1,1,n)<=0.) then
		!	x1fish(1,1,n) = x1fish(1,1,n) + size*dble(mx)
		!endif
		!if(y1fish(1,1,n)>=size*dble(mx)) then
		!	y1fish(1,1,n) = y1fish(1,1,n) - size*dble(mx)
		!endif
		!if(y1fish(1,1,n)<=0.) then
		!	y1fish(1,1,n) = y1fish(1,1,n) + size*dble(mx)
		!endif
!		境界線上に留まる
!		x1fish(1,1,n) = min (max ( x1fish(1,1,n) + dx(n), 0.), size*dble(mx) )
!		y1fish(1,1,n) = min (max ( y1fish(1,1,n) + dy(n), 0.), size*dble(mx) )

		fi(n) = int( x1fish(1,1,n) / size ) + 1
		fj(n) = int( y1fish(1,1,n) / size ) + 1
		

!		ボイドモデルのみ
!		if(t == 1) then
	!		dx(n) = dx_co(n) * para_co + dx_sep(n) * para_sep + dx_ali(n) * para_ali + 50.
	!		dy(n) = dy_co(n) * para_co + dy_sep(n) * para_sep + dy_ali(n) * para_ali + 50.
			
	!		x1fish(1,1,n) = min (max ( x1fish(1,1,n) + dx(n), 0.), size*dble(mx) )
	!		y1fish(1,1,n) = min (max ( y1fish(1,1,n) + dy(n), 0.), size*dble(mx) )
			
	!		fi(n) = int( x1fish(1,1,n) / size ) + 1
	!		fj(n) = int( y1fish(1,1,n) / size ) + 1
			
	!		write(*,*) t, n, dx(n), dy(n), x1fish(1,1,n), y1fish(1,1,n)
			
!		else
	!		call random_number(random3)
	!		dx(n) = dx_co(n) * para_co + dx_sep(n) * para_sep + dx_ali(n) * para_ali + DCOS( z * random3 * 360. ) * 10.
	!		dy(n) = dy_co(n) * para_co + dy_sep(n) * para_sep + dy_ali(n) * para_ali + DSIN( z * random3 * 360. ) * 10.
			
	!		x1fish(1,1,n) = min (max ( x1fish(1,1,n) + dx(n), 0.), size*dble(mx) )
	!		y1fish(1,1,n) = min (max ( y1fish(1,1,n) + dy(n), 0.), size*dble(mx) )
			
	!		fi(n) = int( x1fish(1,1,n) / size ) + 1
	!		fj(n) = int( y1fish(1,1,n) / size ) + 1
			
	!		write(*,*) t, n, DCOS( z * random3 * 360. ) * 10., DSIN( z * random3 * 360. ) * 10., x1fish(1,1,n), y1fish(1,1,n)
			
!		endif
	enddo
	
	
!	近い距離にある個体数
	school = 0.
	near = 0.
	do n = 1, nn
		num_near(n) = 0.
		do m = 1, nn
			distance(n,m) = sqrt((x1fish(1,1,m) - x1fish(1,1,n))**2 + (y1fish(1,1,m) - y1fish(1,1,n))**2)
			if(distance(n,m)<=(5.*bl/tan(1./12.))) then		!整列で感知できる範囲内
				num_near(n) = num_near(n) + 1.
			endif
!			write(*,*) "distance", n, x1fish(1,1,n), y1fish(1,1,n), m, x1fish(1,1,m), y1fish(1,1,m), distance(n,m)
		enddo
		near = near + num_near(n)
	enddo
	school = near/(nn*(nn-1))
!	write(*,*) near, nn*(nn-1), school
	
!	構造物周辺・最適温度帯の個体数
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
		if(sqrt((ycon(1)- y1fish(1,1,n))**2 + (xcon(1) - x1fish(1,1,n))**2) <= 500.) then
			num_con = num_con + 1.
		endif
        if(abs(y1fish(1,1,n)-size* DBLE(nx) *(temp_max-12.)/(30.-12.))<=500.) then
		!if(abs(y1fish(1,1,n)-size*50.*(temp_max-12.)/(30.-12.))<=500.) then
			num_temp = num_temp + 1.
		endif
	enddo
	
	if(t==lasttime) then
		do n = 1, nn
		r(n) = sqrt((ycon(1)- y1fish(1,1,n))**2 + (xcon(1) - x1fish(1,1,n))**2)
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
75	format('harimalagrangeresult/min', I4.4, '.csv')
76  format('harimalagrangeresult/para_num/fish_number.csv')
77  format(3(I6),1(f10.8))
78  format('harimalagrangeresult/para_num/parameter.csv')
		   
!	結果表示
	    
    write(fname,75)  int(t*deltat)
	open(75,file=fname)      !座標ファイル
    do n = 1, nn
		write(75,70) n,x1fish(1,1,n),y1fish(1,1,n),fi(n),fj(n)
	enddo
	close(75)
	
!	
	open(unit=76, file='harimalagrangeresult/para_num/fish_number.csv', status='old', position='append', action='write')
		write(76,77) int(t*deltat), int(num_con), int(num_temp), school
	close(76)
	
	enddo
	
	
	open(unit=78, file='harimalagrangeresult/para_num/parameter.csv', status='replace', action='write', encoding='UTF-8')
!		個体数,グリッドサイズ,タイムステップ,実行時間(h),温度選好下限,温度選好最適,温度選好上限
!		エサの濃度,構造物の幅(m),拡散係数,構造物の座標(x),構造物の座標(y),魚の速度(m/min),結合,分離,整列
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
	write(*,*) 'Complete!'
	
end program simulator