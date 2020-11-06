****************
**************** Adaptation of Extention of the Mantel test
**************** presented by Sang-Il Lee in his paper
**************** "Lee, S. I. (2004). A generalized significance testing method for global 
****************  measures of spatial association: an extension of the Mantel test. 
****************  Environment and Planning A, 36(9), 1687-1703."
** @Author: Thor69

*========================================================
*==== First Bivariate spatial association  ==============
*========================================================
** First: Matrix weight V, accoording notation of Sang-Il Lee 
** A vector of Z, a z-standarized of values
**

**** Import matrix to MATA
mata: X = st_matrix("Z_x")
mata: Y = st_matrix("Z_y")
mata: V   = st_matrix("V")
mata: N = rows(V)
mata: v1 = J(N,1,1)

**** Standardize varibles
mata: mean_X = mean(X)
mata: mean_Y = mean(Y)
mata: std_X	 = ((1/(N-1))*trace((X - v1*mean_X)*(X - v1*mean_X)'))^(1/2)
mata: std_Y  = ((1/(N-1))*trace((Y - v1*mean_Y)*(Y - v1*mean_Y)'))^(1/2)
mata: Z_x    = (X - v1*mean_X)*(1/std_X)
mata: Z_y    = (Y - v1*mean_Y)*(1/std_Y)

****
**** Change diagonal of V matrix: zero diagonal to 1's.
mata:
V_1 = V
for(i=1; i<=N; i++){
V_1[i,i] = 1
}
end

mata: V = V_1

*** Statistict of Bivariate L:
***
mata: BB_num = Z_x'*(V'*V)*Z_y
mata: BB_den = v1'*(V'*V)*v1
mata: BB = BB_num/(BB_den)

** Definition of 
mata: Q = Z_x*Z_y'
mata: P = (V'*V)/(v1'*(V'*V)*v1)

****Table 1 of Paper:
mata: F_off_0 = v1'*P*v1 - trace(P)
mata: F_on_0  = trace(P)
mata: F_off_1 = trace(P'*P) - ( diagonal(P)'*diagonal(P))
mata: F_on_1  = diagonal(P)'*diagonal(P)
mata: F_off_2 = (P*v1 - diagonal(P))'*(P*v1 - diagonal(P))
mata: F_all_2 = (P*v1)'*(P*v1) 
mata: G_off_0 = v1'*Q*v1 - trace(Q)
mata: G_on_0  = trace(Q)
*** Replaca trace(Q'*Q) for other algorithm method 
mata: Q_1 = Z_x
mata
for(i=1; i<=N; i++){
Q_1[i] = Q[,i]'*Q[,i]
}
end
mata: Q_trace = v1'*Q_1
**
mata: G_off_1 =  Q_trace - ( diagonal(Q)'*diagonal(Q))
mata: G_on_1  = diagonal(Q)'*diagonal(Q)
mata: G_off_2 = (Q*v1 - diagonal(Q))'*(Q*v1 - diagonal(Q))
mata: G_all_2	= (Q*v1)'*(Q*v1) 


*** Tabla 2
****Expected 
mata: E_gamma_off = F_off_0*G_off_0/(N*(N-1))
mata: E_gamma_on  = F_on_0*G_on_0/N
****Variance
mata: var_gamma_off = 2*F_off_1*G_off_1/(N*(N-1)) + 4*( F_off_2 - F_off_1 )*(G_off_2 - G_off_1)/(N*(N-1)*(N-2)) + (F_off_0^2 + 2*F_off_1 - 4*F_off_2)*(G_off_0^2 + 2*G_off_1 - 4*G_off_2)/(N*(N-1)*(N-2)*(N-3)) - E_gamma_off^2
mata: var_gamma_on  = F_on_1*G_on_1/N  + (F_on_0^2 - F_on_1)*(G_on_0^2 - G_on_1)/(N*(N-1)) - E_gamma_on^2
mata: cov_var       = (F_all_2 - F_on_1 - F_off_2)*(G_all_2 - G_on_1 - G_off_2)/(2*N*(N-1)) + ( F_on_0*F_off_0 - (F_all_2 - F_on_1 - F_off_2) )*( G_on_0*G_off_0 - (G_all_2 - G_on_1 - G_off_2) )/(N*(N-1)*(N-2)) - (E_gamma_off)*(E_gamma_on)

*** Estadisticos de Test:
***

mata: E_gamma   = E_gamma_off + E_gamma_on
mata: Var_gamma = var_gamma_off + var_gamma_on + 2*cov_var
mata: std_gamma = Var_gamma^(1/2)
mata: Z_gamma   = (BB - E_gamma)/std_gamma

*** Export to Stata:
mata: st_matrix("E_gamma", E_gamma)
mata: st_matrix("std_gamma", std_gamma)
mata: st_matrix("stat", BB)
mata: st_matrix("Z_gamma", Z_gamma)

*** Adicional Statistics:
***SSS:
mata: SSS_x = (Z_x'*(V'*V)*Z_x)/(v1'*(V'*V)*v1)
mata: SSS_y = (Z_y'*(V'*V)*Z_y)/(v1'*(V'*V)*v1)
***Pearson
mata: pearson_xy = (1/N)*trace(Q)

***Pearson hat
mata:
V_n = V*v1
for(i=1; i<=N; i++){
V_n[i] = V_n[i]^(-1)
}
V_n=editmissing(V_n, 0)
V_std = V
for(i=1; i<=N; i++){
for(j=1; j<=N; j++){
V_std[i,j] = V_n[i]*V_std[i,j]
}
}
end
***
mata: Y_hat 		 = (V_std*Y)
mata: X_hat 		 = (V_std*X)
mata: mean_X_hat 	 = mean(X_hat)
mata: mean_Y_hat 	 = mean(Y_hat)
mata: std_X_hat 	 = ((1/(N-1))*trace((X_hat - v1*mean_X_hat)*(X_hat - v1*mean_X_hat)'))^(1/2)
mata: std_Y_hat 	 = ((1/(N-1))*trace((Y_hat - v1*mean_Y_hat)*(Y_hat - v1*mean_Y_hat)'))^(1/2)
mata: Z_x_hat 		 = (X_hat - v1*mean_X_hat)*(1/std_X_hat)
mata: Z_y_hat 		 = (Y_hat - v1*mean_Y_hat)*(1/std_Y_hat)
mata: Q_hat 		 = Z_x_hat*Z_y_hat'
mata: pearson_xy_hat = (1/N)*trace(Q_hat)
***


****SSS_x
mata: A1 =(((X_hat - v1*mean_X_hat)'*(X_hat - v1*mean_X_hat))^(1/2))/(((X_hat - v1*mean_X)'*(X_hat - v1*mean_X))^(1/2))
mata: A2 =(((Y_hat - v1*mean_Y_hat)'*(Y_hat - v1*mean_Y_hat))^(1/2))/(((Y_hat - v1*mean_Y)'*(Y_hat - v1*mean_Y))^(1/2))
mata: A = A1*A2
mata: B1 = (mean_X_hat - mean_X)*(v1'*(Y_hat - v1*mean_Y))
mata: B2 = ( ( (X - v1*mean_X)'*(X - v1*mean_X) )*( (Y - v1*mean_Y)'*(Y - v1*mean_Y) ) )^(1/2) 
mata: B = B1/B2
mata: L_aprox = (SSS_x^(1/2))*(SSS_y^(1/2))*pearson_xy_hat 
mata: L_aprox_large = (SSS_x^(1/2))*(SSS_y^(1/2))*A*pearson_xy_hat + B

****Moran_An
mata: Moran_An = (Z_x'*V_std*Z_y)/(Z_x'*Z_x)
mata: Wartenberg = (Z_x'*V*Z_y)/(v1'*V*v1)
****Export to 
mata: st_matrix("SSS_x"   	, SSS_x)
mata: st_matrix("SSS_y"		, SSS_y)
mata: st_matrix("pearson_xy", pearson_xy)
mata: st_matrix("pearson_xy_hat",	pearson_xy_hat)
mata: st_matrix("L_aprox"	, L_aprox)
mata: st_matrix("Moran_An"	, Moran_An)
mata: st_matrix("L_aprox_large"	, L_aprox_large)
mata: st_matrix("Wartenberg"	, Wartenberg)


mat Results_M = [ stat, E_gamma, std_gamma, Z_gamma, SSS_x, SSS_y, pearson_xy, pearson_xy_hat, L_aprox, L_aprox_large, Moran_An, Wartenberg ] 
mat colnames Results_M = stat E std z SSS_x SSS_y pearson_xy peason_xy_hat L_aprox L_aprox_large Moran_An Wartenberg


*** Vectores cuyo punto producto son represendos por el estadistico
mata: z_x_hat1 		 = (X_hat - v1*mean_X)*(1/std_X)
mata: z_y_hat1 		 = (Y_hat - v1*mean_Y)*(1/std_Y)
mata: st_matrix("z_x_hat1"	, z_x_hat1)
mata: st_matrix("z_y_hat1"	, z_y_hat1)
mata: st_matrix("z_x"	, Z_x)
mata: st_matrix("z_y"	, Z_y)

/*
preserve
clear
mat ZZ = [z_x_hat1 , z_y_hat1 , z_x, z_y]
mat colnames ZZ = z_x_hat  z_y_hat z_x z_y
svmat ZZ, names(col)
scatter z_x_hat z_y_hat, yline(0) xline(0) title("Scatter of elements of Lee's Stat.") name(graphs_vector1, replace) 
scatter z_x z_y, yline(0) xline(0) title("Scatter of elements of Pearson Corr.") name(graphs_vector2, replace) 
restore
*/

