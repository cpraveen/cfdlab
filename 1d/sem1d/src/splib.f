      SUBROUTINE VALEPO(N,X,Y,DY,D2Y)                                   
**************************************************************          
*   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N           
*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
*   N  = DEGREE OF THE POLYNOMIAL                                       
*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
**************************************************************          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
                                                                        
         Y   = 1.D0                                                     
         DY  = 0.D0                                                     
         D2Y = 0.D0                                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
         Y   = X                                                        
         DY  = 1.D0                                                     
         D2Y = 0.D0                                                     
      IF(N .EQ. 1) RETURN                                               
                                                                        
         YP   = 1.D0                                                    
         DYP  = 0.D0                                                    
         D2YP = 0.D0                                                    
      DO 1 I=2,N                                                        
         C1 = DFLOAT(I)                                                 
         C2 = 2.D0*C1-1.D0                                              
         C4 = C1-1.D0                                                   
         YM = Y                                                         
         Y  = (C2*X*Y-C4*YP)/C1                                         
         YP = YM                                                        
         DYM  = DY                                                      
         DY   = (C2*X*DY-C4*DYP+C2*YP)/C1                               
         DYP  = DYM                                                     
         D2YM = D2Y                                                     
         D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1                       
         D2YP = D2YM                                                    
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZELEGL(N,ET,VN)                                        
*********************************************************************   
*   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA   
*   N  = ORDER OF THE FORMULA                                           
*   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*)                                        
      IF (N .EQ. 0) RETURN                                              
                                                                        
         N2 = (N-1)/2                                                   
         SN = DFLOAT(2*N-4*N2-3)                                        
         ET(0) = -1.D0                                                  
         ET(N) = 1.D0                                                   
         VN(0) = SN                                                     
         VN(N) = 1.D0                                                   
      IF (N .EQ. 1) RETURN                                              
                                                                        
         ET(N2+1) = 0.D0                                                
         X = 0.D0                                                       
      CALL VALEPO(N,X,Y,DY,D2Y)                                         
         VN(N2+1) = Y                                                   
      IF(N .EQ. 2) RETURN                                               
                                                                        
         PI = 3.14159265358979323846D0                                  
         C  = PI/DFLOAT(N)                                              
      DO 1 I=1,N2                                                       
         ETX = DCOS(C*DFLOAT(I))                                        
      DO 2 IT=1,8                                                       
      CALL VALEPO(N,ETX,Y,DY,D2Y)                                       
         ETX = ETX-DY/D2Y                                               
2     CONTINUE                                                          
         ET(I) = -ETX                                                   
         ET(N-I) = ETX                                                  
         VN(I) = Y*SN                                                   
         VN(N-I) = Y                                                    
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DMLEGL(N,NM,ET,VN,DMA)                                 
************************************************************************
*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE        
*  LEGENDRE GAUSS-LOBATTO NODES                                         
*  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX              
*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
*  ET  = VECTOR OF THE NODES, ET(I), I=0,N                              
*  VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N                      
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), DMA(0:NM,0:*)                         
          DMA(0,0) = 0.D0                                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
      DO 1 I=0,N                                                        
          VI = VN(I)                                                    
          EI = ET(I)                                                    
      DO 2 J=0,N                                                        
      IF (I .NE. J) THEN                                                
          VJ = VN(J)                                                    
          EJ = ET(J)                                                    
          DMA(I,J) = VI/(VJ*(EI-EJ))                                    
      ELSE                                                              
          DMA(I,I) = 0.D0                                               
      ENDIF                                                             
2     CONTINUE                                                          
1     CONTINUE                                                          
                                                                        
          DN = DFLOAT(N)                                                
          C  = .25D0*DN*(DN+1.D0)                                       
          DMA(0,0) = -C                                                 
          DMA(N,N) = C                                                  
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
      SUBROUTINE WELEGL(N,ET,VN,WT)                                     
*********************************************************************** 
*   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA 
*   N  = ORDER OF THE FORMULA                                           
*   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N                       
*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                            
*********************************************************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), WT(0:*)                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          N2 = (N-1)/2                                                  
          DN = DFLOAT(N)                                                
          C  = 2.D0/(DN*(DN+1.D0))                                      
      DO 1 I=0,N2                                                       
          X = ET(I)                                                     
          Y = VN(I)                                                     
          WTX = C/(Y*Y)                                                 
          WT(I) = WTX                                                   
          WT(N-I) = WTX                                                 
1     CONTINUE                                                          
                                                                        
      IF(N-1 .EQ. 2*N2) RETURN                                          
          X = 0.D0                                                      
          Y = VN(N2+1)                                                  
          WT(N2+1) = C/(Y*Y)                                            
                                                                        
      RETURN                                                            
      END                                                               
