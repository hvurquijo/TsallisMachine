!#######################_PROGRAM: PLOTING STUFFS_#######################################
program Tsallis_graphs
   use OMP_LIB
  implicit none
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_VARIABLES DEFINITION_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double precision :: beta
  double precision, allocatable:: beta_mat(:)

  real :: q_inf, q_sup, tmp, q
  real, allocatable :: q_mat(:)
  integer :: k, i, ii, n, m, cycl, iii, fi, mbeta, threadnum              !n is the number of spins, m is the size o the interval for variating q
                                                   ! cycl is the number of times the program will aproximate the mean values
  integer :: flag, flag2
  character*30 :: file_name, spec1
  character*8 :: date
  character*6 :: time_char
  character*5 :: q_char, beta_char
  character*2 :: n_char, threadnumchar
  character :: CR
  logical :: exist

  double precision, allocatable :: si_mean_exp(:,:,:),si_mean_TS(:,:,:),h_TS(:,:,:)
  double precision, allocatable :: sij_mean_exp(:,:,:,:),sij_mean_TS(:,:,:,:), J_TS(:,:,:,:)
  double precision, allocatable :: sijk_mean_exp(:,:,:,:,:),sijk_mean_TS(:,:,:,:,:)
  double precision, allocatable :: Cij_exp(:,:,:,:), Cij_TS(:,:,:,:)


  double precision, allocatable :: J(:,:), h(:)
  double precision, allocatable :: si_mean(:,:,:),sij_mean(:,:,:)
  double precision, allocatable :: sijk_mean(:,:,:), E_J_mean(:,:,:), E_h_mean(:,:,:)
  real :: T_start, T_finish
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_INTERFACES_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  interface
     subroutine matwrite(MATRIX)
       double precision, intent(in) :: MATRIX(:,:)
     end subroutine matwrite
 
     subroutine Jhgen(N,TYPE,INTERACTION,FIELD)
       integer, intent(in) :: N
       character*2, intent(in) :: TYPE !'++' for J,h in [0,1] ; '--' for J,h in [-1,1]
       double precision, dimension(N,N), intent(out) :: INTERACTION
       double precision, dimension(N), intent(out) :: FIELD
     end subroutine Jhgen

     subroutine simple_save(data, name_des,spec, date, time)
       CHARACTER*30, INTENT(IN) :: NAME_DES,spec
       DOUBLE PRECISION, INTENT(IN) :: DATA(:) !vector
       CHARACTER*8, INTENT(IN) :: date
       CHARACTER*6, INTENT(IN) :: time
     end subroutine simple_save

     subroutine simple_save2(data, name_des,spec, date, time)
       CHARACTER*30, INTENT(IN) :: NAME_DES,spec
       DOUBLE PRECISION, INTENT(IN) :: DATA(:,:) !vector
       CHARACTER*8, INTENT(IN) :: date
       CHARACTER*6, INTENT(IN) :: time
     end subroutine simple_save2

     subroutine save3D_as_list(n,data,name_des,spec, date, time)!save all non zero values in a vector of square matrix
         integer, intent(in) :: n
         CHARACTER*30, INTENT(IN) :: NAME_DES,spec
         DOUBLE PRECISION, INTENT(IN) :: DATA(:,:,:) !vector
         CHARACTER*8, INTENT(IN) :: date
         CHARACTER*6, INTENT(IN) :: time
      end subroutine save3D_as_list

      subroutine save2D_as_list(n,data,name_des,spec, date, time)!save all non zero values in a vector of square matrix
         integer, intent(in) :: n
         CHARACTER*30, INTENT(IN) :: NAME_DES,spec
         DOUBLE PRECISION, INTENT(IN) :: DATA(:,:) !vector
         CHARACTER*8, INTENT(IN) :: date
         CHARACTER*6, INTENT(IN) :: time
      end subroutine save2D_as_list

     subroutine TS_mean_calculator(h,J,q,q_exp,beta,n,cycl,si_mean,sij_mean,sijk_mean,err_J,err_h,&
      mean_si_exp,mean_si_TS,mean_sij_exp,mean_sij_TS,mean_sijk_exp, mean_sijk_TS,Cij_exp,Cij_TS,J_TS,h_TS)
      double precision, intent(in) :: beta, J(:,:), h(:)
      real, intent(in) :: q, q_exp
      integer, intent(in) :: n, cycl
      double precision, intent(out), dimension(cycl) :: si_mean,sij_mean,sijk_mean,err_J,err_h
      double precision, intent(out), dimension(n)     :: mean_si_exp,   mean_si_TS, h_TS
      double precision, intent(out), dimension(n,n)   :: mean_sij_exp,  mean_sij_TS, J_TS
      double precision, intent(out), dimension(n,n)   :: Cij_exp,       Cij_TS
      double precision, intent(out), dimension(n,n,n) :: mean_sijk_exp, mean_sijk_TS
     end subroutine TS_mean_calculator   

     function dec2bin(ARG1,ARG2)
      integer, intent(in) :: ARG1
      integer, intent(in) :: ARG2
      integer, dimension(ARG2) :: DEC2BIN
    end function dec2bin
  end interface
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%_END INTERFACES_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_PROGRAM BODY_ENERGY DISTRIBUITION_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM(), OMP_GET_MAX_THREADS()
  threadnum=OMP_GET_MAX_THREADS()
  write(threadnumchar, "(i2)") threadnum-1
  call system("export OMP_NUM_THREADS="//trim(threadnumchar))
  call system("export OMP_NESTED=TRUE")

  !print *, "Write the interval of q to search (Ex. 0.5):  q = "
  !read *, q_inf, q_sup !reads a value for q<1, and the upper and lower limits for the q matrix
  q_inf = 0.1
  q_sup = 2
  print*,"Write q"
  read*,q
  !print*,q==1
    
  !print*, "Write the resolution of q to search"
  !read*, m  !number of messures
  !m = 19

  print*, "Write n"
  read*, n  !number of messures
  
  print *, "Write a number of cycles the program will run cycl>0 (type 0 for default cycl = 200):  cycl = "
  read *, cycl
  if (cycl < 1) then
     cycl = 200
     print *, "Setting cycl = 200"
  endif


 !allocate(q_mat(m+1))
 allocate(J(n,n),h(n))
 call Jhgen(n,'++',J,h) !This soubroutine generates random values of J,h
 call matwrite(J)
 print*,h


!------------------------CORREÇÃO_FLUTUAÇÃO_DE_ENERGIA-------------------------
!  if(q==0.1)then
!      J = 2.349*J !Obtido para 10 spins
!      h = 2.349*h !Obtido para 10 spins
!  elseif(q==0.2)then
!      J = 1.849*J !Obtido para 10 spins
!      h = 1.849*h !Obtido para 10 spins
!      !J = 86*J !Obtido para 10 spins
!      !h = 86*h !Obtido para 10 spins
!  elseif(q==0.3)then
!      J = 1.464*J !Obtido para 10 spins
!      h = 1.464*h !Obtido para 10 spins
!   elseif(q==0.4)then
!      J = 1.165*J !Obtido para 10 spins
!      h = 1.165*h !Obtido para 10 spins
!   elseif(q==0.5)then
!      J = 0.932*J !Obtido para 10 spins
!      h = 0.932*h !Obtido para 10 spins
!   elseif(q==0.6)then
!      J = 0.750*J !Obtido para 10 spins
!      h = 0.750*h !Obtido para 10 spins
!   elseif(q==0.7)then
!      J = 0.607*J !Obtido para 10 spins
!      h = 0.607*h !Obtido para 10 spins
!   elseif(q==0.8)then
!      J = 0.495*J !Obtido para 10 spins
!      h = 0.495*h !Obtido para 10 spins
!   elseif(q==0.9)then
!      J = 0.405*J !Obtido para 10 spins
!      h = 0.405*h !Obtido para 10 spins
!   elseif(q==1)then
!      J = 0.335*J !Obtido para 10 spins
!      h = 0.335*h !Obtido para 10 spins
!   elseif(q==1.1)then
!      J = 0.717*J !Obtido para 10 spins
!      h = 0.717*h !Obtido para 10 spins
!   elseif(q==1.2)then
!      J = 2.002*J !Obtido para 10 spins
!      h = 2.002*h !Obtido para 10 spins
!   elseif(q==1.3)then
!   elseif(q==1.4)then
!   elseif(q==1.5)then
!   elseif(q==1.6)then
!   elseif(q==1.7)then
!   elseif(q==1.8)then
!   elseif(q==1.9)then
!   elseif(q==2)then
!   else
!      J = 1*J 
!      h = 1*h 
!   endif
!-------------------------FIN_CORREÇÃO_ENERGIA---------------------------------------


 
  tmp =0
  flag = 0
  flag2 = 0
  m=1
  allocate(q_mat(m))
  !q_mat = (/0.1,0.5,0.7,0.8,0.9,1.0,1.2/)
  q_mat = (/0.8/)
  mbeta=8
  allocate(beta_mat(mbeta))
  beta_mat = (/0.5,0.38,0.30,0.21,0.16,0.12,0.09,0.07/)
  print*,q_mat

  !m=m+1

  !print*,"antes de alocar"
  allocate(si_mean(mbeta,m,cycl),sij_mean(mbeta,m,cycl),sijk_mean(mbeta,m,cycl), E_J_mean(mbeta,m,cycl), E_h_mean(mbeta,m,cycl))
  allocate(si_mean_exp(mbeta,m,n),sij_mean_exp(mbeta,m,n,n),sijk_mean_exp(mbeta,m,n,n,n),Cij_exp(mbeta,m,n,n))
  allocate(si_mean_TS(mbeta,m,n),sij_mean_TS(mbeta,m,n,n),sijk_mean_TS(mbeta,m,n,n,n), Cij_TS(mbeta,m,n,n))
  allocate(h_TS(mbeta,m,n),J_TS(mbeta,m,n,n))
  !print*,"depois de alocar"
  call cpu_time(T_start)
  do ii = 1,size(beta_mat)
   BETA =beta_mat(ii)
      !$OMP PARALLEL
      !$OMP DO
         do i=1,m
            call TS_mean_calculator(H,J,q_mat(i),q,BETA,n,CYCL,si_mean(ii,i,:),sij_mean(ii,i,:),&
                  sijk_mean(ii,i,:),E_J_mean(ii,i,:),E_h_mean(ii,i,:),si_mean_exp(ii,i,:),si_mean_TS(ii,i,:),sij_mean_exp(ii,i,:,:)&
                  ,sij_mean_TS(ii,i,:,:),sijk_mean_exp(ii,i,:,:,:),sijk_mean_TS(ii,i,:,:,:),Cij_exp(ii,i,:,:),Cij_TS(ii,i,:,:),&
                  J_TS(ii,i,:,:),h_TS(ii,i,:))
                  
            !write(6,'(A,F5.1,A)', advance='no') '<----------------------------------------',REAL(i)/(m_menor+m_maior)*100.0,&
             !     "% Completed---------------------------------------------->"//CR
            !call SLEEP(1)
         enddo
      !$OMP END DO
      !$OMP END PARALLEL
  print *, "<----------------------------------------100% Completed---------------------------------------------->"
  !description = "RESULTADOS COM REGRA MODIFICADA: REGRA*Q/(Q-1)"

  
   if (n<10)then
      write(n_char,"(i1)") n
      n_char="0"//n_char
   else
      write(n_char, "(i2)") n !converts the integer value of n to a char
   end if

   if (beta>=10 .and. beta<=100)then
      write(beta_char,"(F4.1)") beta
   elseif(beta<10 .and. beta >= 0.01)then
      write(beta_char,"(F4.2)") beta
   endif

  if(q<1)then
      write(q_char, "(F3.2)") q
      spec1 = "Qme_0"//trim(q_char)//"Beta"//trim(beta_char)//"_N"//trim(n_char)
   endif
   if(q>=1)then
      write(q_char, "(F4.2)") q
      spec1 = "Qme_"//trim(q_char)//"Beta"//trim(beta_char)//"_N"//trim(n_char)
   endif
  
  
  print *, "Folder ",spec1," created"
  

  call date_and_time(date, time_char) !calls date and time and saves them at the character type variables date and time
  
  
  file_name = "<Si>_vs_Cycl"
  call simple_save2(si_mean(ii,:,:), file_name,spec1, date, time_char) !file_name, spec, date and time are used to construct the proper file name 
  file_name = "<Sij_vs_Cycl>"
  call simple_save2(sij_mean(ii,:,:),file_name,spec1, date, time_char)
  file_name = "<Sijk_vs_Cycl>"
  call simple_save2(sijk_mean(ii,:,:),file_name,spec1, date, time_char)
  file_name = "<E_J>_vs_Cycl"
  call simple_save2(E_J_mean(ii,:,:),file_name,spec1, date, time_char)
  file_name = "<E_h>_vs_Cycl"
  call simple_save2(E_h_mean(ii,:,:),file_name,spec1, date, time_char)
  file_name = "q_matrix"
  call simple_save(DBLE(q_mat),file_name,spec1, date, time_char)
  file_name = "si_mean_exp_vs_q"
  call simple_save2(si_mean_exp(ii,:,:),file_name,spec1, date, time_char)
  file_name = "si_mean_TS_vs_q"
  call simple_save2(si_mean_TS(ii,:,:),file_name,spec1, date, time_char)
   !print*,i
  do i=1,m
     write(file_name,"(f4.2)") q_mat(i)
     file_name="Cij_exp_"//trim(file_name)
     !print*, "inside loop"
     call save2D_as_list(n,Cij_exp(ii,i,:,:),file_name,spec1, date, time_char)
  enddo
  !print*,"fez o do1"
  do i=1,m
     write(file_name,"(f4.2)") q_mat(i)
     file_name="Cij_TS_"//trim(file_name)
     call save2D_as_list(n,Cij_TS(ii,i,:,:),file_name,spec1, date, time_char)
  enddo
  !print*,"fez o do
  do i=1,m
   write(file_name,"(f4.2)") q_mat(i)
   file_name="si_mean_exp_"//trim(file_name)
   call simple_save(si_mean_exp(ii,i,:),file_name,spec1, date, time_char)
enddo
!print*,"fez o do3"
   do i=1,m
      write(file_name,"(f4.2)") q_mat(i)
      file_name="si_mean_TS_"//trim(file_name)
      call simple_save(si_mean_TS(ii,i,:),file_name,spec1, date, time_char)
   enddo
  do i=1,m
     write(file_name,"(f4.2)") q_mat(i)
     file_name="sij_mean_exp_"//trim(file_name)
     call save2D_as_list(n,sij_mean_exp(ii,i,:,:),file_name,spec1, date, time_char)
  enddo
  !print*,"fez o do3"
  do i=1,m
     write(file_name,"(f4.2)") q_mat(i)
     file_name="sij_mean_TS_"//trim(file_name)
     call save2D_as_list(n,sij_mean_TS(ii,i,:,:),file_name,spec1, date, time_char)
  enddo
  do i=1,m
     write(file_name,"(f4.2)") q_mat(i)
     file_name="Tijk_mean_exp_"//trim(file_name)
     call save3D_as_list(n,sijk_mean_exp(ii,i,:,:,:),file_name,spec1, date, time_char)
  enddo
  do i=1,m
     write(file_name,"(f4.2)") q_mat(i)
     file_name="Tijk_mean_TS_"//trim(file_name)
     call save3D_as_list(n,sijk_mean_TS(ii,i,:,:,:),file_name,spec1, date, time_char)
  enddo
   do i=1,m
   write(file_name,"(f4.2)") q_mat(i)
   file_name="J_TS_"//trim(file_name)
   !print*, "inside loop"
   call save2D_as_list(n,J_TS(ii,i,:,:),file_name,spec1, date, time_char)
   enddo
   do i=1,m
      write(file_name,"(f4.2)") q_mat(i)
      file_name="h_TS_"//trim(file_name)
      call simple_save(h_TS(ii,i,:),file_name,spec1, date, time_char)
   enddo


   file_name="J"
   !print*, "inside loop"
   call save2D_as_list(n,J(:,:),file_name,spec1, date, time_char)

   file_name="h"
   call simple_save(h(:),file_name,spec1, date, time_char)






  call cpu_time(T_finish)
  print*,(T_finish-T_start)
  INQUIRE( FILE = "win31.mp3", EXIST = exist )
  if (exist)then
   call system("mplayer win31.mp3 > /dev/null 2>&1")
  else
   print*, char(7)
  endif

  INQUIRE(FILE = "send_email.py", EXIST=exist)
  if (exist)then
   open(fi, FILE = "tmpf95.txt", status = 'new')
   write(fi,*) "The Tsallis machine with especifications:\ "//spec1//",\ has finished"
   close(fi)
   call system("python3 send_email.py")
   call system("rm tmpf95.txt")
  endif

enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_END PROGRAM BODY_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end program Tsallis_Graphs
!###########################_END PROGRAM_###############################################




















!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&                                                                                                                                               &
!&                  HERE STARTS THE DEFINITIONS OF THE FUNCTIONS AND SUBROUTINES REFERED AT THE INTERFACE MAIN PROGRAM                           &
!&                                                                                                                                               &
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gsign(NUMBER)
   implicit none

   double precision, intent(in) :: NUMBER
   double precision :: gsign

   if(NUMBER == 0.d0)then
      gsign = 0
   elseif(NUMBER<0.d0)then
      gsign = -1
   else
      gsign = 1
   endif

   return
end function   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%_CRIATES DAIGONAL MATRIX FROM VECTOR_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diag(VECTOR)
  implicit none

  double precision, intent(in) :: VECTOR(:)
  integer :: i,j,k
  double precision, dimension(size(VECTOR), size(VECTOR)) :: diag
  k = size(VECTOR)
  diag=0.d0
  forall (i=1:k, j=1:k, i==j)
     diag(i,j)=VECTOR(i)
  end forall
  return
end function diag
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_GENERATES A TRIANGULAR MATRIX_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#################################################################################################################################################
!#                                                                                                                                               #
!# INPUT: ARG1 is the matrix to triangularize                                                                                                    #
!#        ARG2 is an option (could be "u" or "U", "d" or "D") that declares if output matrix would be up("u") triangular or down("d") triangular #
!#        ARG3 indicates the size of the non zero indices, that is: [SIZE of non-zero elements = numer of COLUMNS - ARG3].                       #
!# This implies that ARG3 must be lower than the number of COLUMNS                                                                               #
!#                                                                                                                                               #
!# OUTPUT: A triangular (up or down) matrix with same dimensions as ARG1 and size of the non zero elemnts [ = numer of COLUMNS - ARG3]           #
!#                                                                                                                                               #
!#################################################################################################################################################
function tri(MATRIX, OPTION, DISPLACEMENT)
  implicit none

  double precision, intent(in) :: MATRIX(:,:)
  character, intent(in) :: OPTION
  integer, intent(in) :: DISPLACEMENT
  integer :: i,j, size_(2)
  double precision, dimension(ubound(MATRIX,1),ubound(MATRIX,2)) :: tri

  size_=shape(MATRIX)
  !allocate(double precision :: tri(size_(1),size_(2)))
  tri = 0.d0

  if (DISPLACEMENT<size_(2))then
     if (option=="u" .OR. option=="U") then
        forall (i=1:ubound(MATRIX,1), j=1:ubound(MATRIX,2), j>i+DISPLACEMENT-1)
           tri(i,j)=MATRIX(i,j)
        end forall
     elseif (option=="d" .OR. option=="D") then
        forall (i=1:ubound(MATRIX,1), j=1:ubound(MATRIX,2), j<i-DISPLACEMENT+1)
           tri(i,j)=MATRIX(i,j)
        end forall
     end if
  else
     PRINT *, "####################"
     print *, "#   WARNING!!!!!   #"
     PRINT *, "####################"
     print *, " "
     print *, "The third parameter must be an integer smaller than the number of columns!!!"
     print *, " "
  end if
  return
  !deallocate(C)
end function tri
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tri3(MATRIX, OPTION, DISPLACEMENT)
  implicit none

  double precision, intent(in) :: MATRIX(:,:,:)
  character, intent(in) :: OPTION
  integer, intent(in) :: DISPLACEMENT
  integer :: i,j, k, size_(3)
  double precision, dimension(ubound(MATRIX,1),ubound(MATRIX,2),ubound(MATRIX,3)) :: tri3

  size_=shape(MATRIX)
  !allocate(double precision :: tri(size_(1),size_(2)))
  tri3 = 0.d0

  if (DISPLACEMENT<size_(2))then
     if (option=="u" .OR. option=="U") then
        do i=1,ubound(MATRIX,1)
           do j=1,ubound(MATRIX,2)
              do k=1,ubound(MATRIX,3)
                 if(j>i+DISPLACEMENT-1 .AND. k>j+DISPLACEMENT-1)then
                    tri3(i,j,k)=MATRIX(i,j,k)
                 end if
              end do
           end do
        end do
     elseif (option=="d" .OR. option=="D") then
        do i=1,ubound(MATRIX,1)
           do j=1,ubound(MATRIX,2)
              do k=1,ubound(MATRIX,3)
                 if (k<j-DISPLACEMENT+1 .AND. j<i-DISPLACEMENT+1) then
                    tri3(i,j,k)=MATRIX(i,j,k)
                 end if
              end do
           end do
        end do
     endif
  else
     PRINT *, "####################"
     print *, "#   WARNING!!!!!   #"
     PRINT *, "####################"
     print *, " "
     print *, "The third parameter must be an integer smaller than the number of columns!!!"
     print *, " "
  end if
  return
  !deallocate(C)
end function tri3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%_ISING HAMILTONIAN_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#################################################################################################################################################
!#                                                                                                                                               #
!# INPUT: ARG1 is the matrix of J, wich is an upper triangular matrix of dimension NxN                                                           #
!#        ARG2 is the vector conteining the values of h with dimension n                                                                         #
!#        ARG3 is the vector conteining the values of h with dimension n                                                                         #
!#        ARG4 are the options                                                                                                                   #
!#                                                                                                                                               #
!# OPTIONS: "++" para hamiltonianno positivo, "--" para hamiltoniano negativo, "s+" para sigma no intervalo [0,1]                                #
!#                                                                                                                                               #
!# OUTPUT: A double precision scalar                                                                                                             #
!#                                                                                                                                               #
!#########################0########################################################################################################################

function ising(J, h, sigma, q)
implicit none

double precision, intent(in) :: J(:,:)
double precision, intent(in) :: h(:)
double precision, intent(in) :: sigma(:)

double precision, intent(in) :: q
integer :: n, i, ii
double precision :: ising

interface
    function gsign(NUMBER)
       double precision, intent(in) :: NUMBER
       double precision :: gsign
    end function   
end interface

n = size(sigma)

ising = 0.d0
do i=1,n
 ising = ising - h(i)*(sigma(i)-gsign((1-q)*h(i)))
 do ii = 1,n
    if(ii>i)then
       ising = ising - J(i,ii)*(sigma(i)*sigma(ii)-gsign((1-q)*J(i,ii)))
    endif
 enddo
enddo

return
end function ising



!%%%%%%%%%%%%%%%%%%%%%%%%_J,H GENERATORS_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Jhgen(n,type,J,h)
  implicit none

  integer, intent(in) :: n
  character*2, intent(in) :: type
  double precision, dimension(n,n), intent(out) :: J
  double precision, dimension(n), intent(out) :: h

  interface
     function tri(A, option1, option2)
       double precision, intent(in) :: A(:,:)
       character, intent(in) :: option1
       integer, intent(in) :: option2
       double precision, DIMENSION(UBOUND(A,1),UBOUND(A,2)) :: TRI
     end function tri
  end interface

  
  if (type == '++') then
     call random_number(J)
     J = tri(J,"u",1)
     call random_number(h)
     h = 2*h-1
  elseif (type == '-+') then
     call random_number(J)
     J = 2*J-1
     
     J = tri(J,"u",1)
     call random_number(h)
     h = 2*h-1
     print*,h
  elseif (type == '--') then
     call random_number(J)
     J =J-1 
     J = tri(J,"u",1)
     call random_number(h)
     !h = h-1
     h = 2*h-1
  end if

  return
  !deallocate(J,h)
  
end subroutine Jhgen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%_SIGMA MATRIX GENERATOR_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigmamat(sigma)
  implicit none

  double precision, intent(in) :: sigma(:)
  double precision, dimension(size(sigma),size(sigma)) :: sigmamat
  integer :: n, i, j

  n=size(sigma)
  sigmamat=0
  do i=1,n
     do j=1,n
        if(j>i)then
           sigmamat(i,j) = sigma(i)*sigma(j)
        endif
     end do
  end do
  return
  
end function sigmamat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%_SIGMA MATRIX3 GENERATOR_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigmamat3(sigma)
  implicit none

  double precision, intent(in) :: sigma(:)
  double precision, dimension(size(sigma), size(sigma), size(sigma)) :: sigmamat3
  integer :: n, i, j, k

  n=size(sigma)
  sigmamat3=0

  
  do i=1,n
     do j=1,n
        do k=1,n
           if(k>j .AND. j>i)then
              sigmamat3(i,j,k)=sigma(i)*sigma(j)*sigma(k)
           end if
        end do
     end do
  end do
 
  return
  
end function sigmamat3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_TRANSFORMS FROM DECIAML TO BINARY _%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dec2bin(NUMBER,REPRESENTATION)
  implicit none

  integer, intent(in) :: NUMBER
  integer, intent(in) :: REPRESENTATION !NUMBER OF BITS IN THE REPRESENTATION
  integer, dimension(REPRESENTATION) :: dec2bin
  integer, dimension(REPRESENTATION) :: tmp_bin_vect
  integer :: i,j, tmp_num
  
  tmp_num = NUMBER
  
  do i=1,REPRESENTATION
     dec2bin(i)=mod(tmp_num,2)
     tmp_num = ishft(tmp_num,-1)
  end do
  
  do j=1,REPRESENTATION
     tmp_bin_vect(j)=dec2bin(1+REPRESENTATION-j)
  end do
  dec2bin=tmp_bin_vect
  return
end function dec2bin
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%_FOR DEBUGGING. WRITES A GIVEN MxN MATRIX AT THE SCREEN_%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine matwrite(MATRIX)
  implicit none
  
  double precision, intent(in) :: MATRIX(:,:)
  !double precision :: mat(2,2)
  integer :: i
   
    do i=lbound(MATRIX,1), ubound(MATRIX,1)
        print *, MATRIX(i,:)
    end do
  !deallocate(mat)
  return
end subroutine matwrite
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_SAVES DATA TO FILE_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine simple_save(data, name_des,spec, date, time)
  implicit none

  CHARACTER*30, INTENT(IN) :: NAME_DES,spec
  DOUBLE PRECISION, INTENT(IN) :: DATA(:) !vector
  CHARACTER*8, INTENT(IN) :: date
  CHARACTER*6, INTENT(IN) :: time
  

  character*255 :: folder_name, file_name, cwd, des
  !character*8 :: date
  !character*6 :: time
  integer :: stat, fi, ST
  logical :: exist
  ST = 0

  des = ' '

  !call date_and_time(date, time)
  call getcwd(cwd)
  
  folder_name = date//trim(spec)
  !print *, time
  !print *, trim(name_des)
  file_name = trim(cwd)//'/'//trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//'/'//trim(name_des)//'.dat'
  !print *, file_name
  !print *, folder_name
  !print *, trim(cwd)
  fi=1
  
  INQUIRE( FILE = trim(folder_name)//'/.', EXIST = exist )
  
  if (.not. exist) then
     call system('mkdir '//trim(folder_name),stat)
     print *, "Status = ",stat, " (0 means no problems found while saving)"
     if (stat == 0) then
        print *, "Folder", trim(folder_name)//" created"
     else
        print *, "Error ocurred during creation of folder "//folder_name
        ST = -1
     end if
  end if
  
  INQUIRE( FILE = trim(folder_name)//'/'//trim(spec)//'/.', EXIST = exist )
  
  if (.not. exist) then
     call system('mkdir '//trim(folder_name)//'/'//trim(spec),stat)
     print *, "Status = ",stat, "(0 means no problems found while saving)"
     if (stat == 0) then
        print *, trim(folder_name)//'/'//trim(spec)//" folder created"
     else
        print *, "Error ocurred during creation of folder "//trim(folder_name)//'/'//trim(spec)
        ST = -1
     end if
  end if

  INQUIRE( FILE = trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//'/.', EXIST = exist )
  
  if (.not. exist) then
     call system('mkdir '//trim(folder_name)//'/'//trim(spec)//'/'//trim(time),stat)
     print *, "Status = ",stat, " (0 means no problems found while saving)"
     if (stat == 0) then
        print *, trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//" folder created"
     else
        print *, "Error ocurred during creation of folder "//trim(folder_name)//'/'//trim(spec)//'/'//trim(time)
        ST = -1
     end if
  end if
  
  open(fi, FILE = file_name, status = 'new')
  
  write(1,*) data
   
  close(fi)
  
  
  
end subroutine simple_save

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_SAVES DATA TO FILE_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine simple_save2(data, name_des,spec, date, time)
  implicit none

  CHARACTER*30, INTENT(IN) :: NAME_DES,spec
  DOUBLE PRECISION, INTENT(IN) :: DATA(:,:) !vector
  CHARACTER*8, INTENT(IN) :: date
  CHARACTER*6, INTENT(IN) :: time
  

  character*255 :: folder_name, file_name, cwd, des
  !character*8 :: date
  !character*6 :: time
  integer :: stat, fi, ST, k, size(2)
  logical :: exist
  ST = 0

  des = ' '

  !call date_and_time(date, time)
  call getcwd(cwd)
  
  folder_name = date//trim(spec)
  !print *, time
  !print *, trim(name_des)
  file_name = trim(cwd)//'/'//trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//'/'//trim(name_des)//'.dat'
  !print *, file_name
  !print *, folder_name
  !print *, trim(cwd)
  fi=1
  
  INQUIRE( FILE = trim(folder_name)//'/.', EXIST = exist )
  
  if (.not. exist) then
     call system('mkdir '//trim(folder_name),stat)
     print *, "Status = ",stat, " (0 means no problems found while saving)"
     if (stat == 0) then
        print *, "Folder ",trim(folder_name)//" created"
     else
        print *, "Error ocurred during creation of folder "//folder_name
        ST = -1
     end if
  end if
  
  INQUIRE( FILE = trim(folder_name)//'/'//trim(spec)//'/.', EXIST = exist )
  
  if (.not. exist) then
     call system('mkdir '//trim(folder_name)//'/'//trim(spec),stat)
     print *, "Status = ",stat, " (0 means no problems found while saving)"
     if (stat == 0) then
        print *, trim(folder_name)//'/'//trim(spec)//" folder created"
     else
        print *, "Error ocurred during creation of folder "//trim(folder_name)//'/'//trim(spec)
        ST = -1
     end if
  end if

  INQUIRE( FILE = trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//'/.', EXIST = exist )
  
  if (.not. exist) then
     call system('mkdir '//trim(folder_name)//'/'//trim(spec)//'/'//trim(time),stat)
     print *, "Status = ",stat, " (0 means no problems found while saving)"
     if (stat == 0) then
        print *, trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//" folder created"
     else
        print *, "Error ocurred during creation of folder "//trim(folder_name)//'/'//trim(spec)//'/'//trim(time)
        ST = -1
     end if
  end if
  
  open(fi, FILE = file_name, status = 'new')
  size = shape(data)
  do k = 1,size(1)
     write(1,*) data(k,:)
  end do
   
  close(fi)
  
  
  
end subroutine simple_save2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%_SAVE MATRIX DATA TO FILE AS A LIST_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine save2D_as_list(n,data,name_des,spec, date, time)!save all non zero values in a vector of square matrix
   implicit none
   integer, intent(in) :: n

   CHARACTER*30, INTENT(IN) :: NAME_DES,spec
   DOUBLE PRECISION, INTENT(IN) :: DATA(:,:) !vector
   CHARACTER*8, INTENT(IN) :: date
   CHARACTER*6, INTENT(IN) :: time
  

  character*255 :: folder_name, file_name, cwd, des
  integer :: stat, fi, ST, k, size(2), i, ii
  logical :: exist
  double precision, allocatable :: array(:)

  allocate(array(0))

   do i=1,n
      do ii=1,n
         if(data(i,ii)==0.d0)then
            array = [array,data(i,ii)]
         endif
      enddo
   enddo

   ST = 0

   des = ' '
 
   !call date_and_time(date, time)
   call getcwd(cwd)
   
   folder_name = date//trim(spec)
   !print *, time
   !print *, trim(name_des)
   file_name = trim(cwd)//'/'//trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//'/'//trim(name_des)//'.dat'
   !print *, file_name
   !print *, folder_name
   !print *, trim(cwd)
   fi=1
   
   INQUIRE( FILE = trim(folder_name)//'/.', EXIST = exist )
   
   if (.not. exist) then
      call system('mkdir '//trim(folder_name),stat)
      print *, "Status = ",stat, " (0 means no problems found while saving)"
      if (stat == 0) then
         print *, "Folder", trim(folder_name)//" created"
      else
         print *, "Error ocurred during creation of folder "//folder_name
         ST = -1
      end if
   end if
   
   INQUIRE( FILE = trim(folder_name)//'/'//trim(spec)//'/.', EXIST = exist )
   
   if (.not. exist) then
      call system('mkdir '//trim(folder_name)//'/'//trim(spec),stat)
      print *, "Status = ",stat, "(0 means no problems found while saving)"
      if (stat == 0) then
         print *, trim(folder_name)//'/'//trim(spec)//" folder created"
      else
         print *, "Error ocurred during creation of folder "//trim(folder_name)//'/'//trim(spec)
         ST = -1
      end if
   end if
 
   INQUIRE( FILE = trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//'/.', EXIST = exist )
   
   if (.not. exist) then
      call system('mkdir '//trim(folder_name)//'/'//trim(spec)//'/'//trim(time),stat)
      print *, "Status = ",stat, " (0 means no problems found while saving)"
      if (stat == 0) then
         print *, trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//" folder created"
      else
         print *, "Error ocurred during creation of folder "//trim(folder_name)//'/'//trim(spec)//'/'//trim(time)
         ST = -1
      end if
   end if
   
   open(fi, FILE = file_name, status = 'new')
   
   write(1,*) data
    
   close(fi)
   
   
end subroutine save2D_as_list
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%_SAVE MATRIX DATA TO FILE AS A LIST_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine save3D_as_list(n,data,name_des,spec, date, time)!save all non zero values in a vector of square matrix
   implicit none

   integer, intent(in) :: n

   CHARACTER*30, INTENT(IN) :: NAME_DES,spec
   DOUBLE PRECISION, INTENT(IN) :: DATA(:,:,:) !vector
   CHARACTER*8, INTENT(IN) :: date
   CHARACTER*6, INTENT(IN) :: time
  

  character*255 :: folder_name, file_name, cwd, des
  integer :: stat, fi, ST, k, size(2), i, ii, iii
  logical :: exist
  double precision, allocatable :: array(:)

  allocate(array(0))

   do i=1,n
      do ii=1,n
         do iii = 1,n
            if(data(i,ii,iii)==0.d0)then
               array = [array, data(i,ii,iii)]
            endif
         enddo
      enddo
   enddo

   ST = 0

   des = ' '
 
   !call date_and_time(date, time)
   call getcwd(cwd)
   
   folder_name = date//trim(spec)
   !print *, time
   !print *, trim(name_des)
   file_name = trim(cwd)//'/'//trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//'/'//trim(name_des)//'.dat'
   !print *, file_name
   !print *, folder_name
   !print *, trim(cwd)
   fi=1
   
   INQUIRE( FILE = trim(folder_name)//'/.', EXIST = exist )
   
   if (.not. exist) then
      call system('mkdir '//trim(folder_name),stat)
      print *, "Status = ",stat, " (0 means no problems found while saving)"
      if (stat == 0) then
         print *, "Folder", trim(folder_name)//" created"
      else
         print *, "Error ocurred during creation of folder "//folder_name
         ST = -1
      end if
   end if
   
   INQUIRE( FILE = trim(folder_name)//'/'//trim(spec)//'/.', EXIST = exist )
   
   if (.not. exist) then
      call system('mkdir '//trim(folder_name)//'/'//trim(spec),stat)
      print *, "Status = ",stat, "(0 means no problems found while saving)"
      if (stat == 0) then
         print *, trim(folder_name)//'/'//trim(spec)//" folder created"
      else
         print *, "Error ocurred during creation of folder "//trim(folder_name)//'/'//trim(spec)
         ST = -1
      end if
   end if
 
   INQUIRE( FILE = trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//'/.', EXIST = exist )
   
   if (.not. exist) then
      call system('mkdir '//trim(folder_name)//'/'//trim(spec)//'/'//trim(time),stat)
      print *, "Status = ",stat, " (0 means no problems found while saving)"
      if (stat == 0) then
         print *, trim(folder_name)//'/'//trim(spec)//'/'//trim(time)//" folder created"
      else
         print *, "Error ocurred during creation of folder "//trim(folder_name)//'/'//trim(spec)//'/'//trim(time)
         ST = -1
      end if
   end if
   
   open(fi, FILE = file_name, status = 'new')
   
   write(1,*) data
    
   close(fi)
   
   
end subroutine save3D_as_list
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#######################_SUBROUTINE: MAQUINA DE TSALLIS_#######################################
subroutine TS_mean_calculator(h,J,q,q_exp,beta,n,cycl,si_mean,sij_mean,sijk_mean,err_J,err_h,&
   mean_si_exp,mean_si_TS,mean_sij_exp,mean_sij_TS,mean_sijk_exp,mean_sijk_TS,Cij_exp,Cij_TS, J_REACHED,H_REACHED)
 implicit none
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_VARIABLES DEFINITION_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 double precision, intent(in) :: beta, J(:,:), h(:)
 real, intent(in) :: q,q_exp
 integer, intent(in) :: n, cycl

 double precision, intent(out), dimension(cycl) :: si_mean,sij_mean,sijk_mean,err_J,err_h !<|Delta<si>|>, <|Delta<sij>|>, <|Delta<Tij>|>,<|Delta EJ|>, <|Delta Eh|>

  double precision, intent(out), dimension(n)     :: mean_si_exp,   mean_si_TS,h_reached
  double precision, intent(out), dimension(n,n)   :: mean_sij_exp,  mean_sij_TS,J_reached
  double precision, intent(out), dimension(n,n)   :: Cij_exp,       Cij_TS
  double precision, intent(out), dimension(n,n,n) :: mean_sijk_exp, mean_sijk_TS
 
 integer :: ns
 !double precision, dimension(n,n,n) :: mean_sijk_exp, mean_sijk_TS !<Tijk>exp, <Tijk>TS
 double precision, dimension(n,n) :: delta_J!, mean_sij_exp, mean_sij_TS !J_TS(aproximação), <Sij>exp, <Sij>TS, erro de J
 double precision, dimension(n) :: delta_h!, mean_si_exp, mean_si_TS !h_TS(aproximação), <Si>exp, <Si>TS, erro de h
 double precision :: tmp, tmp1, Zq_exp, Zq, prob !temporal variables and partition functions
 integer :: i,ii,k !for "do" cycles


double precision, allocatable :: Zq_exp_mat(:), states(:,:), Zq_mat(:) !matrices of partition fucntions and all the states, good for parallelization
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_INTERFACES_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 interface
    subroutine matwrite(MATRIX) !for debugging only
      double precision, intent(in) :: MATRIX(:,:)
    end subroutine matwrite

    function ising(INTERACTION, FIELD, SIGMA,q)
      double precision, intent(in) :: INTERACTION(:,:)
      double precision, intent(in) :: FIELD(:), SIGMA(:)
      double precision :: ising
      double precision, intent(in)::q
    end function ising

    function sigmamat(SIGMA)
      double precision, intent(in) :: SIGMA(:)
      double precision, dimension(size(sigma),size(sigma)) :: sigmamat
    end function sigmamat

    function sigmamat3(SIGMA)
      double precision, intent(in) :: SIGMA(:)
      double precision, dimension(size(sigma),size(sigma),size(sigma)) :: sigmamat3
    end function sigmamat3

    function qZpart(INTERACTION, FIELD, BETA, Q, DIST_TYPE, HAM_TYPE, states)
      double precision, intent(in) :: INTERACTION(:,:)
      double precision, intent(in) :: FIELD(:)
      double precision, intent(in) :: BETA, states(:,:)
      integer, intent(in) :: DIST_TYPE
      character*2, intent(in) :: HAM_TYPE !'s+' for calculates ising with (sigma+1)/2
      real, intent(in) :: Q
      double precision  :: qZpart
    end function qZpart
    
    function dec2bin(NUMBER,REPRESENTATION)
      integer, intent(in) :: NUMBER
      integer, intent(in) :: REPRESENTATION
      integer, dimension(REPRESENTATION) :: dec2bin
    end function dec2bin
 end interface
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%_END INTERFACES_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_PROGRAM BODY_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ns = 2**n
!Messures = size(data,1)
allocate(Zq_exp_mat(ns), states(n,ns), Zq_mat(ns)) !allocating big temporal matrices, good for parallelization 
 !-----------Defines starting values of the J, h aproximations----------------------------

 J_reached = 0.d0           !Value of the interaction calculated by the model
 h_reached = 0.d0           !Value of the field calculated by the  model
 prob = 0.d0                !Value of experimental probability, for debugging
 mean_si_exp = 0.d0         !Value of the experimental mean of sigma_i
 mean_sij_exp = 0.d0        !Value of the experimental mean of (sigma_i*sigma_j)
 mean_sijk_exp = 0.d0       !Value of the experimental mean of triplets
 Zq = 0.d0                  !Experimental partition function
 Zq_exp_mat = 0.d0
 Zq_mat = 0.d0
 !---------------------------------------------------------------------
 
 !------------------------CALCULATES EXPERIMENTAL PROBABILITY AND MEANS-------------------------
do i=0,ns-1                         !Calculates the q-partition function 
   states(:,i+1) = 2*DEC2BIN(i,N)-1
   if (Q_exp==1.0)then
      Zq_exp_mat(i+1) = exp(-beta*ising(J,H,states(:,i+1),DBLE(Q_exp)))
   else
      Zq_exp_mat(i+1) = (1-beta*(Q_exp-1)*ising(J,H,states(:,i+1),DBLE(Q_exp)))
      if(Zq_exp_mat(i+1)<=0)then
         Zq_exp_mat(i+1)=0
      else
         Zq_exp_mat(i+1) = Zq_exp_mat(i+1)**(1/(Q_exp-1))
      endif
   endif
enddo
Zq_exp = sum(Zq_exp_mat)

 do i=0,ns-1     
   if(Q_exp==1.0)then                               
      tmp = exp(-beta*ising(J,H,states(:,i+1),DBLE(Q_exp)))/Zq_exp
   else
      tmp = (1-beta*(Q_exp-1)*ising(J,H,states(:,i+1),DBLE(Q_exp)))
      if(tmp<=0)then
         tmp=0
      else
         tmp = (1-beta*(Q_exp-1)*ising(J,H,states(:,i+1),DBLE(Q_exp)))**(1/(Q_exp-1))/Zq_exp      !Experimental probability
      endif
   endif
    
    prob = prob+tmp                                              !For calculation of the sum of the experimental probability, for debugging

       mean_si_exp = mean_si_exp + states(:,i+1)*tmp                 !Calculates <Si>exp
       mean_sij_exp = mean_sij_exp + sigmamat(states(:,i+1))*tmp     !Calculates <Sij>exp
       mean_sijk_exp = mean_sijk_exp + sigmamat3(states(:,i+1))*tmp  !Calculates <Tijk>exp
 end do
 !-----------------------------------------------------------------------------------------


!##################################_HERE STARTS THE MACHINE_###########################################
 label1:do k=1,cycl  !this "do" controls the number of aproximations the time will make
    
    mean_si_TS = 0.d0      !Value of the sigma_i's mean calculated by second model
    mean_sij_TS = 0.d0     !Value of the mean of (sigma_i*sigma_j) calculated by second model
    mean_sijk_TS = 0.d0
  
    !----------------Calculates the (sigma_i), (sigma_i*sigma_j) and triplets by the model-------------------
   !Zq_mat=0.d0
   Zq=0.d0

   do i=0,ns-1                         !Calculates the q-partition function 
      if(q==1.0)then
         Zq_mat(i+1)=exp(-beta*ising(J_REACHED,H_REACHED,states(:,i+1),DBLE(Q)))
         !print*, Zq_mat(i+1)
      else
         Zq_mat(i+1) = (1-beta*(q-1)*ising(J_REACHED,H_REACHED,states(:,i+1),DBLE(Q)))
         if(Zq_mat(i+1)<=0)then
            Zq_mat(i+1)=0
         else
            Zq_mat(i+1) = Zq_mat(i+1)**(1/(q-1))
         endif
         
         !if(q>1 .AND. q<1.1 .AND. Zq_mat(i+1)>1)then
         !   print*, Zq_mat(i+1)
         !   print*,(q-1)
         !   print*,ising(J_REACHED,H_REACHED,states(:,i+1),DBLE(Q))   
         !endif
      endif
   enddo
   Zq=sum(Zq_mat)
   !if(q>1 .AND. q<1.1)then
   !   print*,Zq
   !endif
   
   !---------------------------------------------------------------------------------------------------------
   !---------------------------------CALCULATES TSALLIS PROBABILITY AND MEANS-------------------------
    do i=0,ns-1
      if(Q==1.0)then
         tmp1=exp(-beta*ising(J_REACHED,H_REACHED,states(:,i+1),DBLE(Q)))/Zq
         !print*, tmp1
      else
         tmp1 = (1-beta*(Q-1)*ising(J_REACHED,H_REACHED,states(:,i+1),DBLE(Q)))
         if(tmp1<=0)then
            tmp1=0
         else
            tmp1 = tmp1**(1/(Q-1))/Zq
         endif

         !if(q>1 .AND. q<1.1)then
         !   print*, Zq_mat(i+1)
         !endif
      endif
       !print*, tmp1
        mean_si_TS = mean_si_TS + states(:,i+1)*tmp1                 !Calculates <Si>TS
        mean_sij_TS = mean_sij_TS + sigmamat(states(:,i+1))*tmp1     !Calculates <Sij>TS
        mean_sijk_TS = mean_sijk_TS + sigmamat3(states(:,i+1))*tmp1  !Calculates <Tijk>TS
    end do
   !-----------------------------------------------------------------------------------------------------

    if(k==cycl)then !for debuggin
       print *, " "
       print *, "q experimental: ",q_exp
       print *, "q aproximação: ",q
       print *, "A diferença das probabilidades é ", log(abs(tmp-tmp1))
       print *, " "
    end if
    
    !-------------------------------------------------------------------------------------------

    !----------Update rules--------------------------------------------------------------
    delta_h = -(mean_si_exp-mean_si_TS)    !Final part of the second rule's calculation
    delta_J = -(mean_sij_exp-mean_sij_TS)  !Final part of the second rule's calculation
    !-------------------------------------------------------------
    
    !------------------Updating--------------------------------
    H_REACHED = H_REACHED-DELTA_H 
    J_REACHED = J_REACHED-DELTA_J

    !------------------------------------------------------------------------------------

    sijk_mean(k) = sum(abs(mean_sijk_exp-mean_sijk_TS))/((n-1)*(n-2)*(n)/6) !tripletos
    si_mean(k) = sum(abs(mean_si_exp-mean_si_TS))/n
    sij_mean(k) = sum(abs(mean_sij_exp-mean_sij_TS))/((n-1)*(n)/2)
    err_J(k) = sum(abs(J_REACHED-J))/((n-1)*(n)/2)
    err_h(k) = sum(abs(H_REACHED-H))/n
    
    
  enddo label1
!###########################_HERE ENDS THE MACHINE_########################################

do i=1,n
   do ii=1,n
      Cij_exp(i,ii) = mean_sij_exp(i,ii)-mean_si_exp(i)*mean_si_exp(ii)
   enddo
enddo

do i=1,n
   do ii=1,n
      Cij_TS(i,ii) = mean_sij_TS(i,ii)-mean_si_TS(i)*mean_si_TS(ii)
   enddo
enddo

 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_END PROGRAM BODY_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end subroutine TS_Mean_Calculator
!###########################_END PROGRAM_###################################################