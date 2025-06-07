program grid_temperature
  implicit none
  integer, parameter :: nx = 42, ny = 73
  double precision :: min_temp, max_temp, temp_step
  double precision :: temps(ny)
  integer :: i, j
  character(len=100) :: filename
  integer :: unit

  ! ���x�͈͐ݒ�
  min_temp = 15.0d0
  max_temp = 30.0d0
  temp_step = (max_temp - min_temp) / dble(ny - 1)

  ! ���x�������ɐ��`�ω�������
  do j = 0, ny - 1
    temps(j+1) = min_temp + temp_step * dble(j)
  end do

  ! CSV�t�@�C���o��
  filename = 'grid_temperature.csv'
  open(unit=10, file=filename, status='replace')

  !write(10,*) 'row,col,temperature'
  do i = 0, nx - 1
    do j = 0, ny - 1
      write(10,'(I3,",",I3,",",F8.2)') i, j, temps(j+1)
    end do
  end do

  close(10)
  print *, 'CSV�t�@�C���ɏo�͂��܂����F', trim(filename)

end program grid_temperature