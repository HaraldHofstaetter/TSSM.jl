#some schemes
#Emb 4/3 AK p
embedded_scheme_43 = EmbeddedScheme(
          ( 0.125962888700250514,  0.333588446797901933,
            0.751193431379145450, -0.338296598434303506,
            0.127551831557005609,  0.127551831557005609,
           -0.338296598434303506,  0.751193431379145450,
            0.333588446797901933,  0.125962888700250514 ),

          ( 0.125962888700250514,  0.333588446797901933,
            0.751193431379145450, -0.338296598434303506,
            0.0                 ,  0.261153550449697153,
           -0.242703571757396124,  0.596114052266110425,
            0.365547251678000160,  0.147440548920593995 ),

            4 )
#PP 3/4 A
palindromic_scheme_34 = PalindromicScheme( 
          ( 0.268330095781759925,  0.919661523017399857, 
           -0.187991618799159782, -0.187991618799159782, 
            0.919661523017399857,  0.268330095781759925 ),
            3 )
#PP 5/6 A
palindromic_scheme_56 = PalindromicScheme(
          ( 0.201651044312324230,   0.578800656272664932, 
            0.562615975356569200,   0.273128836056524479, 
            0.253874038247554845,  -0.102733803148432142, 
           -0.835351693190370636,   0.068014946093165092, 
            0.068014946093165092,  -0.835351693190370636,
           -0.102733803148432142,   0.253874038247554845, 
            0.273128836056524479,   0.562615975356569200, 
            0.578800656272664932,   0.201651044312324230 ),            
            5 )
#defectbased
embedded_scheme_43_D = DefectBasedScheme(
            embedded_scheme_43.scheme2,
            embedded_scheme_43.order-1 )
palindromic_scheme_34_D = DefectBasedScheme( 
            palindromic_scheme_34.scheme,
            palindromic_scheme_34.order )
palindromic_scheme_56_D = DefectBasedScheme(
            palindromic_scheme_56.scheme,            
            palindromic_scheme_56.order )
schemes=[embedded_scheme_43  ,palindromic_scheme_34  ,palindromic_scheme_56,
         embedded_scheme_43_D,palindromic_scheme_34_D,palindromic_scheme_56_D];