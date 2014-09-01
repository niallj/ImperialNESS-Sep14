      SUBROUTINE INVERT(B,A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Inverts a real 3 by 3 matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(15,99)
      REAL(kind = double), DIMENSION(3,3) :: A, B 
      REAL(kind = double) :: COM 
      INTEGER :: I,J,K

      DO  I = 1,3              
        DO  J = 1,3         
          A(I,J) = B(I,J)                                        
        END DO
      END DO

      DO K = 1,3                                           
        COM = A(K,K)                                           
        A(K,K) = 1.0D0                                           
        DO J = 1,3                                                    
          A(K,J) = A(K,J)/COM
        END DO
        DO 5 I = 1,3                                    
        IF((I-K) == 0) THEN                                
          GOTO 5
        END IF
        COM = A(I,K)                                    
        A(I,K) = 0.0
        DO J = 1,3                                    
          A(I,J) = A(I,J) - COM*A(K,J)                      
        END DO
    5   CONTINUE                                      
      END DO

      END SUBROUTINE   
