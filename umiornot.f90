program convert_dat
  implicit none
  integer, parameter :: nx = 42, ny = 73
  ! nx�͏c�̃O���b�h���Any�͉��̃O���b�h��
  double precision :: A(nx, ny)
  integer :: i, j
  character(len=100) :: infile, outfile

  infile = 'new0629depth.dat'
  outfile = '0629chikei.dat'

  ! ���̓t�@�C���ǂݍ���
  open(unit=10, file=infile, status='old')
  do i = 1, nx
    read(10,*) (A(i,j), j = 1, ny)
  enddo
  close(10)

  ! �ϊ��F0�ȊO��1��
  A = merge(1.0d0, 0.0d0, A /= 0.0d0)

  ! �o�̓t�@�C���i�����Łj
  open(unit=11, file=outfile, status='replace')
  do i = 1, nx
  write(11, '(73(I1,1X))') (int(A(i,j)), j = 1, ny)
  enddo
  close(11)

  print *, '�ϊ������F', trim(outfile)

end program convert_dat
