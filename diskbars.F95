PROGRAM diskbars

! ESTE CODIGO FORTRAN RESOLVE A CINEMATICA DE UM MECANISMO BIELA-MANIVELA.
! TRES OPÇOES DE CONFIGURAÇÃO SÃO POSSIVEIS.
! Referencia: SHAMES, Irving H. Engineering Mechanics.
! Programação inicial por Luciano FLEISCHFRESSER em 23 JULHO DE 2013.
! Última atualização em 21 OUTUBRO DE 2014.

INTEGER :: I, J, K
INTEGER, PARAMETER :: T = 20
REAL, PARAMETER :: PI = 3.141592654
REAL :: AB, OA, BC, OC, PSS_NGLR, MG_PNT
REAL, DIMENSION(0:21) :: L, M, N, O, P, Q, LF, BT, GM_NM, GM_DN, GM, &
MG, FTR, V_A, V_B, MG_AB, MG_BC, &
BTRM1, BTRM2, BTRM3, BTRM4, DNAB, MG_AB_PNT, &
BCTRM1, BCTRM2, BCTRM3, BCTRM4, DNBC, MG_BC_PNT

20 WRITE(*,*) "Digite:"
WRITE(*,*) "1: se deseja simular o mecanismo original"
WRITE(*,*) "2: se deseja simular o mecanismo que permite giro completo"
WRITE(*,*) "3: se deseja simular o mecanismo com velocidade variavel"
WRITE(*,*) "4: se deseja sair do programa"
READ(*,*) J
SELECT CASE (J)
CASE (1)

OPEN (unit = 10, file = "mecanismo1.dat", status = "old",  &
         action = "write", form = "formatted", position = "rewind")

OA = 0.3000000000000
BC = 0.3000000000000
OC = 1.3000000000000
AB = SQRT((OC - BC)**2 + OA**2)
PSS_NGLR = 4.5 * ( PI / 180.0 )
LF(0) = 0.0
DO I = 0, T
MG(I) = 2.0
V_A(I) = MG(I) * OA
L(I) = SIN (LF(I))
M(I) = COS (LF(I))
BT(I) = ACOS ( 1.0 - L(I))
N(I) = SIN (BT(I))
O(I) = COS (BT(I))
GM_NM(I) = BC * N(I) + OA * M(I)
GM_DN(I) = AB
GM(I) = ASIN (GM_NM(I) / GM_DN(I))
P(I) = SIN(GM(I))
Q(I) = COS(GM(I))
FTR(I) = (M(I)*Q(I) + L(I)*P(I)) / (N(I)*Q(I) + O(I)*P(I))
V_B(I) = V_A(I) * FTR(I)

MG_AB(I) = (V_A(I) / AB) * (M(I) * O(I) - L(I) * N(I)) / (O(I) * P(I) + N(I) * Q(I))
MG_BC(I) = -V_B(I) / BC

BTRM1(I) = OA * MG(I)**2 * (L(I) * O(I) + M(I) * N(I))

BTRM2(I) = BC * MG_BC(I)**2 * ( O(I)**2 + N(I)**2)

BTRM3(I) = AB * MG_AB(I)**2 * ( O(I) * Q(I) - N(I) * P(I))

DNAB(I) = AB * (O(I) * P(I) + N(I) * Q(I))

MG_AB_PNT(I) = - (BTRM1(I) + BTRM2(I) + BTRM3(I)) / DNAB(I)

BCTRM1(I) = OA * MG(I)**2 * (M(I) * P(I) - L(I) * Q(I))
BCTRM2(I) = MG_AB(I)**2 * AB * (P(I)**2 + Q(I)**2)
BCTRM3(I) = MG_BC(I)**2 * BC * (N(I) * P(I) - O(I) * Q(I))

DNBC(I) = BC * (O(I)*P(I) + N(I)*Q(I))
MG_BC_PNT(I) = (BCTRM1(I) - BCTRM2(I) + BCTRM3(I)) / DNBC(I)
LF(I+1) = LF(I) + PSS_NGLR
END DO
PRINT *
PRINT *, "Mecanismo original"
PRINT *, "Compare os valores para ALFA = 0 com o exemplo resolvido no livro..."
PRINT *
PRINT *, " Alfa ", "      Vel. Ang. BC ", " Acel. Ang. BC "
PRINT *, "(graus)", "      (rad/s) ", "      (rad/s²) "
PRINT *, "_________________________________________________"
PRINT *

WRITE (unit = 10, fmt = "(/, a, /, /, a, a, a, a, /, /, a, a, a, a, /)") "    RESULTADOS P/ GIRO CONST. ANTI-HORARIO", &
              "         ANGULOS (graus)   "," VELOC. LINEARES (m/s) ","   VELOC. ANGULARES (rad/s)  ", &
              " ACEL. ANGULARES (rad/s²)  ",   "    alfa    beta    gama ", "       V_a       V_b    ", &
              "    W        W_ab     W_bc   ", "    dW_ab     dW_bc "

DO I = 0, T
PRINT "(f5.1, f15.3, f15.3)", LF(I) * 180.0 / PI, MG_BC(I), MG_BC_PNT(I)

WRITE (unit = 10, fmt = "(3F8.1, 2x, 2F10.2, 2x, 3f9.2, 2x, F10.2, 2x, F10.2)") & 
       LF(I)*180.0/PI, BT(I)*180.0/PI, GM(I)*180.0/PI, V_A(I), V_B(I), MG(I), MG_AB(I), &
       MG_BC(I), MG_AB_PNT(I), MG_BC_PNT(I)

END DO
PRINT *

CASE (2)

OPEN (unit = 10, file = "mecanismo2.dat", status = "old",  &
         action = "write", form = "formatted", position = "rewind")

OA = 0.3000000000000
BC = 0.3000000000000
AB = 1.044030651
OC = SQRT(AB**2 - (OA+BC)**2)
PSS_NGLR = 4.5 * ( PI / 180.0 )
LF(0) = 0.0
DO I = 0, T-1
MG(I) = -2.0
V_A(I) = MG(I)*OA
L(I) = SIN(LF(I))
M(I) = COS(LF(I))
BT(I) = ACOS (L(I))
N(I) = SIN(BT(I))
O(I) = COS(BT(I))
GM_NM(I) = BC * N(I) + OA * M(I)
GM(I) = ASIN (GM_NM(I) / AB)
If (GM(I) .le. 1.0e-5) Then
GM(I) = 0.0
ENDIf
P(I) = SIN(GM(I))
Q(I) = COS(GM(I))

FTR(I) = (M(I) * Q(I) + L(I) * P(I)) / (O(I) * P(I) - N(I) * Q(I))
V_B(I) = V_A(I) * FTR(I)

MG_AB(I) = -(V_A(I)/AB) * ((M(I) * O(I) + L(I) * N(I))/(N(I) * Q(I) - O(I) * P(I)))
MG_BC(I) = V_B(I) / BC
BTRM1(I) = MG(I)**2 * OA * (L(I)*O(I) - M(I)*N(I))
BTRM2(I) = BC * MG_BC(I)**2 * (O(I)**2 + N(I)**2)
BTRM3(I) = MG_AB(I)**2 * AB * (O(I)*Q(I) - N(I)*P(I))
DNAB(I) = AB * (O(I)*P(I) + N(I)*Q(I))
MG_AB_PNT(I) = (1.0 / DNAB(I)) * (BTRM1(I) - BTRM2(I) - BTRM3(I))
BCTRM1(I) = -OA * MG(I)**2 * (L(I) * Q(I) + M(I) * P(I))
BCTRM2(I) = MG_AB(I)**2 * AB * (Q(I)**2 + P(I)**2)
BCTRM3(I) = MG_BC(I)**2 * BC * (O(I) * Q(I) - N(I) * P(I))
DNBC(I) = BC * (N(I) * Q(I) + O(I) * P(I))
MG_BC_PNT(I) = (BCTRM1(I) + BCTRM2(I) + BCTRM3(I)) / DNBC(I)
LF(I+1) = LF(I) + PSS_NGLR
END DO
PRINT *
PRINT *, "Mecanismo modificado para o projeto computacional com giro completo"
PRINT *,
PRINT *, " Alfa ", "      Vel. Ang. BC ", " Acel. Ang. BC "
PRINT *, "(graus)", "      (rad/s) ", "      (rad/s²) "
PRINT *, "_________________________________________________"
PRINT *

WRITE (unit = 10, fmt = "(/, a, /, /, a, a, a, a, /, /, a, a, a, a, /)") "    RESULTADOS P/ GIRO COMPLETO HORARIO", &
              "         ANGULOS (graus)   "," VELOC. LINEARES (m/s) ","   VELOC. ANGULARES (rad/s)  ", &
              " ACEL. ANGULARES (rad/s²)  ",   "    alfa    beta    gama ", "       V_a       V_b    ", &
              "    W        W_ab     W_bc   ", "    dW_ab     dW_bc "

DO I = 0, T-1
PRINT "(f5.1, f15.3, f15.3)", LF(I) * 180.0 / PI, MG_BC(I), MG_BC_PNT(I)

WRITE (unit = 10, fmt = "(3F8.1, 2x, 2F10.2, 2x, 3f9.2, 2x, F10.2, 2x, F10.2)") & 
       LF(I)*180.0/PI, BT(I)*180.0/PI, GM(I)*180.0/PI, V_A(I), V_B(I), MG(I), MG_AB(I), &
       MG_BC(I), MG_AB_PNT(I), MG_BC_PNT(I)

END DO
PRINT *
WRITE(*,*)"NESTE CASO, A BARRA DE CONEXÃO FICA HORIZONTAL PARA ALFA = 90 GRAUS."
WRITE(*,*)"ENTÃO O PROGRAMA É TERMINADO UM PASSO ANGULAR ANTES PARA EVITAR DIVISÃO POR ZERO"
PRINT *
CASE (3)

OPEN (unit = 10, file = "mecanismo3.dat", status = "old",  &
         action = "write", form = "formatted", position = "rewind")

OA = 0.3000000000000
BC = 0.3000000000000
OC = 1.3000000000000
AB = SQRT((OC - BC)**2 + OA**2)

MG_PNT = 4.37
PSS_NGLR = 4.5 * ( PI / 180.0 )
MG(0) = 0.0
LF(0) = 0.0
DO I = 0, T
L(I) = SIN (LF(I))
M(I) = COS (LF(I))
V_A(I) = MG(I) * OA
BT(I) = ACOS ( 1.0 - L(I))
N(I) = SIN (BT(I))
O(I) = COS (BT(I))
GM_NM(I) = BC * N(I) + OA * M(I)
GM_DN(I) = OC - BC * O(I) - OA * L(I)
GM(I) = ATAN (GM_NM(I) / GM_DN(I))
P(I) = SIN(GM(I))
Q(I) = COS(GM(I))

FTR(I) = (M(I)*Q(I) + L(I)*P(I)) / (N(I)*Q(I) + O(I)*P(I))
MG_AB(I) = (V_A(I) / (AB * Q(I))) * (O(I) * FTR(I) - L(I))
V_B(I) = V_A(I) * FTR(I)
MG_BC(I) = -V_B(I) / BC

BTRM1(I) = MG(I)**2 * OA * (L(I)*O(I) + M(I)*N(I))

BTRM2(I) = MG_PNT * OA * (M(I) * O(I) - L(I) * N(I))

BTRM3(I) = MG_BC(I)**2 * BC * (O(I)**2 + N(I)**2)

BTRM4(I) = MG_AB(I)**2 * AB * (O(I) * Q(I) - N(I) * P(I))

DNAB(I) = AB * (O(I) * P(I) + N(I) * Q(I))

MG_AB_PNT(I) = 1.0/DNAB(I) * (-BTRM1(I)+BTRM2(I)-BTRM3(I)-BTRM4(I))

BCTRM1(I) = MG(I)**2 * OA * (L(I)*Q(I) - M(I)*P(I))

BCTRM2(I) = MG_PNT * OA * (M(I) * Q(I) + L(I) * P(I))

BCTRM3(I) = MG_BC(I)**2 * BC * (N(I)*P(I) - O(I)*Q(I))

BCTRM4(I) = MG_AB(I)**2 * AB * (P(I)**2 + Q(I)**2)

DNBC(I) = BC * (O(I) * P(I) + N(I) * Q(I))

MG_BC_PNT(I) = 1.0/DNBC(I) * (-BCTRM1(I)+BCTRM2(I)+BCTRM3(I)-BCTRM4(I))

LF(I+1) = LF(I) + PSS_NGLR
MG(I+1) = 0.10 + MG_PNT * LF(I+1)
END DO
PRINT *
PRINT *, "Mecanismo modificado para o projeto computacional com velocidade variavel"
PRINT *
PRINT *, " Alfa ", "      Vel. Ang. BC ", " Acel. Ang. BC "
PRINT *, "(graus)", "       (rad/s) ", "      (rad/s²) "
PRINT *, "__________________________________________"
PRINT *

WRITE (unit = 10, fmt = "(/, a, /, /, a, a, a, a, /, /, a, a, a, a, /)") "    RESULTADOS P/ GIRO VAR. ANTI-HORARIO", &
              "         ANGULOS (graus)   "," VELOC. LINEARES (m/s) ","   VELOC. ANGULARES (rad/s)  ", &
              " ACEL. ANGULARES (rad/s²)  ",   "    alfa    beta    gama ", "       V_a       V_b    ", &
              "    W        W_ab     W_bc   ", "    dW_ab     dW_bc "

DO I = 0, T
PRINT "(f5.1, f15.3, f15.3)", LF(I) * 180.0 / PI, MG_BC(I), MG_BC_PNT(I)

WRITE (unit = 10, fmt = "(3F8.1, 2x, 2F10.2, 2x, 3f9.2, 2x, F10.2, 2x, F10.2)") & 
       LF(I)*180.0/PI, BT(I)*180.0/PI, GM(I)*180.0/PI, V_A(I), V_B(I), MG(I), MG_AB(I), &
       MG_BC(I), MG_AB_PNT(I), MG_BC_PNT(I)

END DO
PRINT *
CASE (4)
WRITE(*,*)
WRITE(*,*)"ATÉ A PRÓXIMA!"
WRITE(*,*)
STOP
CASE DEFAULT
WRITE(*,*)
WRITE(*,*)"OPÇÃO INVÁLIDA: ESCOLHA NOVAMENTE"
WRITE(*,*)
END SELECT

GO TO 20
END PROGRAM diskbars
