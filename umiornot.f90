program convert_dat
  implicit none
  integer, parameter :: nx = 42, ny = 73
  ! nxは縦のグリッド数、nyは横のグリッド数
  double precision :: A(nx, ny)
  integer :: i, j
  character(len=100) :: infile, outfile

  infile = 'new0629depth.dat'
  outfile = '0629chikei.dat'

  ! 入力ファイル読み込み
  open(unit=10, file=infile, status='old')
  do i = 1, nx
    read(10,*) (A(i,j), j = 1, ny)
  enddo
  close(10)

  ! 変換：0以外を1に
  A = merge(1.0d0, 0.0d0, A /= 0.0d0)

  ! 出力ファイル（整数で）
  open(unit=11, file=outfile, status='replace')
  do i = 1, nx
  write(11, '(73(I1,1X))') (int(A(i,j)), j = 1, ny)
  enddo
  close(11)

  print *, '変換完了：', trim(outfile)

end program convert_dat
