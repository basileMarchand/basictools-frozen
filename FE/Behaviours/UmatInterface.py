# -*- coding: utf-8 -*-
matstring = """
***behavior gen_evp
 **elasticity isotropic
   young 2.1e5
   poisson 0.3
# **potential gen_evp ep
#  *criterion mises
#  *flow plasticity
#  *isotropic nonlinear
#   R0 180.
#   Q  360.
#   b   30.
***return
"""

DMD0479_C002B = """
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loi de comportement EVP sinh anisotrope de l'AM1 T1R1R2
% DMD0479-C002B
% Applicabilité : pièces coulés
% Traitement thermique : T1R1R2
% Version : B
% Validité : K
% Catégorie :
% Niveau Stat (moyen ou mini) : moyen
% Date : 05/2014
% Expert validation :
% Note : FCT_2014_00079
% Eprouvettes : cylindriques
% Prélèvements : barreaux massifs
% Bornes du domaine d'identification : 20°C-1150°C
% Bornes du domaine extrapolé : 20°C-1250°C
% Remarques :
% Historique :
%     * DMD0479-C002A : loi identifiée par l'Onera dans le cadre du PRC chaud
%     * DMD0479-C002B : correction d'erreurs d'arrondis sur les constantes élastiques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
***material
   %*integration runge_kutta 1.e-3 dtmin 1.e-12
   *integration theta_method_a 0.5 1.e-6 100
***return
***behavior gen_evp
  **coefficient
                masvol 8.9e-9
  **elasticity      cubic
                y1111    temperature
                295918.8125                      0.
                296350.4063                      20.
                297291.0938                      100.
                297837.0938                      150.
                297831.                                               200.
                297768.5                             250.
                297313.6875                      300.
                296454.                                               350.
                295536.4063                      400.
                294050.1875                      450.
                292509.0938                      500.
                290746.5                             550.
                288255.9063                      600.
                285726.6875                      650.
                282665.8125                      700.
                279097.0938                      750.
                275352.8125                      800.
                270996.0938                      850.
                266356.9063                      900.
                261174.2969                      950.
                255498.9063                      1000.
                249389.2031                      1050.
                242644.7031                      1100.
                235340.7969                      1150.
                227185.                                               1200.
                217569.7031                      1250.

                y1122    temperature
                203905.5                             0.
                204683.                                               20.
                206973.2031                      100.
                208479.5938                      150.
                209507.4063                      200.
                210569.4063                      250.
                211314.4063                      300.
                211744.2031                      350.
                212203.                                               400.
                212170.2969                      450.
                212168.4063                      500.
                212024.9063                      550.
                211238.0938                      600.
                210493.5                             650.
                209298.2031                      700.
                207673.7031                      750.
                205966.2031                      800.
                203715.                                               850.
                201269.7031                      900.
                198360.2031                      950.
                195039.4063                      1000.
                191361.                                               1050.
                187129.2031                      1100.
                182319.5938                      1150.
                176936.2969                      1200.
                170388.7969                       1250.

                y1212    temperature
                124875.1016                      0.
                124300.7969                      20.
                122100.1016                      100.
                120670.8984                      150.
                119189.5                             200.
                117688.6016                      250.
                116130.5                             300.
                114547.5                             350.
                112790.3984                      400.
                111247.1016                      450.
                109541.                                               500.
                107805.1016                      550.
                106022.1016                      600.
                104199.2031                      650.
                102343.7031                      700.
                100411.7031                      750.
                98512.5                                               800.
                96543.70313                      850.
                94535.79688                      900.
                92489.79688                      950.
                90407.70313                      1000.
                88292.39844                      1050.
                86132.60156                      1100.
                83941.89844                      1150.
                81719.39844                      1200.
                79264.39844                      1250.

   **thermal_strain isotropic
                alpha     temperature
                1.11E-05              -60.00
                1.11E-05              0.00
                1.12E-05              50.00
                1.16E-05              108.08
                1.18E-05              166.16
                1.21E-05              224.24
                1.23E-05              293.94
                1.25E-05              352.02
                1.26E-05              410.10
                1.27E-05              468.18
                1.29E-05              537.88
                1.31E-05              595.96
                1.32E-05              654.04
                1.34E-05              712.12
                1.36E-05              781.82
                1.39E-05              839.90
                1.42E-05              897.98
                1.45E-05              956.06
                1.50E-05              1025.76
                1.55E-05              1083.84
                1.60E-05              1141.92
                1.66E-05              1200.00
                1.71E-05              1250.00
                ref_temperature            20.0
  **potential octahedral ev
   *flow hyperbolic
                K             temperature
                8.00E+01             -60.0
                8.00E+01             20.0
                8.00E+01             400.0
                1.10E+02             650.0
                1.95E+02             700.0
                3.10E+02             750.0
                3.60E+02             800.0
                3.60E+02             850.0
                3.20E+02             900.0
                2.25E+02             950.0
                1.55E+02             1000.0
                1.20E+02             1050.0
                1.00E+02             1100.0
                8.50E+01             1150.0
                5.50E+01             1250.0
                n             temperature
                2.00E+01             -60.0
                2.00E+01             20.0
                1.80E+01             400.0
                1.50E+01             650.0
                1.20E+01             700.0
                8.00E+00             750.0
                7.00E+00             800.0
                7.00E+00             850.0
                7.00E+00             900.0
                7.00E+00             950.0
                7.00E+00             1000.0
                7.00E+00             1050.0
                7.00E+00             1100.0
                7.00E+00             1150.0
                7.00E+00             1250.0
                eps0      temperature
                5.00E-02              -60.0
                5.00E-02              20.0
                5.00E-02              400.0
                5.00E-02              650.0
                5.00E-02              850.0
                5.00E-02              900.0
                5.00E-02              950.0
                5.00E-02              1000.0
                5.00E-02              1050.0
                5.00E-02              1100.0
                5.00E-02              1150.0
                5.00E-02              1250.0

   Kinf 1.e-4

   *isotropic        constant
                R0           temperature
                2.70E+02             -60.0
                2.70E+02             20.0
                3.10E+02             400.0
                3.10E+02             650.0
                1.95E+02             750.0
                1.60E+02             800.0
                1.00E+02             850.0
                7.50E+01             900.0
                5.00E+01             950.0
                4.50E+01             1000.0
                3.50E+01             1050.0
                2.50E+01             1100.0
                1.50E+01             1150.0
                5.00E+00             1250.0
   *kinematic      nonlinear
                C             temperature
                5.00E+05             -60.0
                5.00E+05             20.0
                8.50E+05             400.0
                8.50E+05             650.0
                7.20E+05             750.0
                6.00E+05             800.0
                4.40E+05             850.0
                3.50E+05             900.0
                2.40E+05             950.0
                1.95E+05             1000.0
                1.60E+05             1050.0
                1.35E+05             1100.0
                7.50E+04             1150.0
                3.85E+04             1250.0
                D             temperature
                2.50E+03             -60.0
                2.50E+03             20.0
                1.80E+04             400.0
                1.20E+04             650.0
                5.00E+03             750.0
                3.20E+03             800.0
                2.40E+03             850.0
                2.05E+03             900.0
                1.60E+03             950.0
                1.65E+03             1000.0
                1.75E+03             1050.0
                1.80E+03             1100.0
                2.50E+03             1150.0
                3.90E+03             1250.0
                m            temperature
                2.19E+00             -60.0
                2.19E+00             20.0
                2.19E+00             400.0
                2.19E+00             650.0
                2.19E+00             850.0
                2.19E+00             900.0
                2.19E+00             950.0
                2.19E+00             1000.0
                2.19E+00             1050.0
                2.19E+00             1100.0
                2.19E+00             1150.0
                2.19E+00             1250.0
                M           temperature
                9.00E+06             -60.0
                9.00E+06             20.0
                2.00E+06             400.0
                1.00E+06             650.0
                5.10E+05             750.0
                3.00E+05             850.0
                2.40E+05             900.0
                1.60E+05             950.0
                7.70E+04             1000.0
                5.00E+04             1050.0
                6.00E+04             1100.0
                7.00E+04             1150.0
                1.00E+05             1250.0
   *interaction    slip
                h1           1.0
                h2           0.0
                h3           0.0
                h4           0.0
                h5           0.0
                h6           0.0
  **potential      cubic cu
   *flow hyperbolic
                K             temperature
                1.00E+02             -60.0
                1.00E+02             20.0
                1.00E+02             650.0
                2.40E+02             750.0
                3.20E+02             800.0
                3.90E+02             850.0
                3.90E+02             900.0
                3.60E+02             950.0
                2.20E+02             1000.0
                1.70E+02             1050.0
                1.60E+02             1100.0
                1.50E+02             1150.0
                1.30E+02             1250.0
                n             temperature
                2.00E+01             -60.0
                2.00E+01             20.0
                1.80E+01             400.0
                1.50E+01             650.0
                1.20E+01             700.0
                8.00E+00             750.0
                6.50E+00             800.0
                5.00E+00             850.0
                5.00E+00             900.0
                5.00E+00             950.0
                5.00E+00             1000.0
                5.00E+00             1050.0
                5.00E+00             1100.0
                5.00E+00             1150.0
                5.00E+00             1250.0
                eps0      temperature
                5.00E-02              -60.0
                5.00E-02              20.0
                5.00E-02              400.0
                5.00E-02              650.0
                5.00E-02              850.0
                5.00E-02              900.0
                5.00E-02              950.0
                5.00E-02              1000.0
                5.00E-02              1050.0
                5.00E-02              1100.0
                5.00E-02              1150.0
                5.00E-02              1250.0

   Kinf 1.e-4

   *isotropic        constant
                R0           temperature
                2.30E+02             -60.0
                2.30E+02             20.0
                2.00E+02             650.0
                1.80E+02             750.0
                1.50E+02             800.0
                1.10E+02             850.0
                8.00E+01             900.0
                6.50E+01             950.0
                5.00E+01             1000.0
                4.50E+01             1050.0
                3.50E+01             1100.0
                2.00E+01             1150.0
                1.00E+01             1250.0
   *kinematic      nonlinear
                C             temperature
                4.15E+05             -60.0
                4.15E+05             20.0
                4.00E+05             650.0
                2.40E+05             750.0
                1.80E+05             850.0
                1.20E+05             900.0
                1.10E+05             950.0
                9.00E+04             1000.0
                7.00E+04             1050.0
                5.00E+04             1100.0
                3.00E+04             1150.0
                1.00E+00             1250.0
                D             temperature
                9.00E+02             -60.0
                9.00E+02             20.0
                2.00E+03             650.0
                2.00E+03             750.0
                2.00E+03             850.0
                1.80E+03             900.0
                1.50E+03             950.0
                1.30E+03             1000.0
                1.20E+03             1050.0
                8.00E+02             1100.0
                8.00E+02             1150.0
                8.00E+02             1250.0
                m            temperature
                2.19E+00             -60.0
                2.19E+00             20.0
                2.19E+00             400.0
                2.19E+00             650.0
                2.19E+00             850.0
                2.19E+00             900.0
                2.19E+00             950.0
                2.19E+00             1000.0
                2.19E+00             1050.0
                2.19E+00             1100.0
                2.19E+00             1150.0
                2.19E+00             1250.0
                M           temperature
                7.00E+06             -60.0
                7.00E+06             20.0
                2.00E+06             650.0
                4.00E+05             750.0
                2.50E+05             850.0
                1.50E+05             900.0
                1.30E+05             950.0
                8.00E+04             1000.0
                5.00E+04             1050.0
                2.00E+04             1100.0
                2.00E+04             1150.0
                1.00E+04             1250.0
   *interaction slip
                h1           1.0
                h2           0.0
                h3           0.0
%            **interaction    iso          ev           cu
%            h             0.0
****return
"""