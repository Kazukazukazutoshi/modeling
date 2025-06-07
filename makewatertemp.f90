program grid_temperature
  implicit none
  integer, parameter :: nx = 42, ny = 73
  double precision :: min_temp, max_temp, temp_step
  double precision :: temps(ny)
  integer :: i, j
  character(len=100) :: filename
  integer :: unit

  ! 温度範囲設定
  min_temp = 15.0d0
  max_temp = 30.0d0
  temp_step = (max_temp - min_temp) / dble(ny - 1)

  ! 温度を列方向に線形変化させる
  do j = 0, ny - 1
    temps(j+1) = min_temp + temp_step * dble(j)
  end do

  ! CSVファイル出力
  filename = 'grid_temperature.csv'
  open(unit=10, file=filename, status='replace')

  !write(10,*) 'row,col,temperature'
  do i = 0, nx - 1
    do j = 0, ny - 1
      write(10,'(I3,",",I3,",",F8.2)') i, j, temps(j+1)
    end do
  end do

  close(10)
  print *, 'CSVファイルに出力しました：', trim(filename)

end program grid_temperature