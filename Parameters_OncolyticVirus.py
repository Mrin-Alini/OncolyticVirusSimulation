
# Parameters of Ad1d24.P19

# lam == Untreated tumor growth per day
lam_Ad1 = 0.0996
# beta == rate of viral infection (of the tumor) per day
beta_Ad1 = 14.51
# k == the transition rate from eclipse to infectious cells (per day)
# # i.e. Infected cells become actively infectious after 1/k time
k_Ad1 = 0.07406
# ro (delta in paper) == rate of infected tumor cell lysis (per day)
ro_Ad1 = 0.07406
# p == rate of virus production (by infected tumor cells)
p_Ad1 = 1.243 * 10**(-4)
# e == virus' sensitivity to interferon (immune response molecule)
e_Ad1 = 47.58
# c == rate of virus clearance (per day)
c_Ad1 = 0.06857
# alpha == decay rate of interferon (immune response molecule) (per day)
alpha_Ad1 = 4.854 * 10**(-13)

params_Ad1 = [lam_Ad1, beta_Ad1, k_Ad1, ro_Ad1, p_Ad1, e_Ad1, c_Ad1, alpha_Ad1]
y0_Ad1 =[41, 0,0, 8.793 * 10**(-3), 0.0]


# Parameters of Ad2d24.P19

lam_Ad2 = 0.0996
beta_Ad2 = 7.297
k_Ad2 = 0.2346
ro_Ad2 = 0.2333
p_Ad2 = 0.02650
e_Ad2 = 58910
c_Ad2 = 0.07022
alpha_Ad2 = 5.804 ** (-8)

params_Ad2 = [lam_Ad2, beta_Ad2, k_Ad2, ro_Ad2, p_Ad2, e_Ad2, c_Ad2, alpha_Ad2]
y0_Ad2 = [41, 0, 0, 6.354 * 10 ** (-4), 0.0]

# Parameters of Ad5.P19

lam_Ad5 = 0.0996
beta_Ad5 = 7.297
k_Ad5 = 0.3198
ro_Ad5 = 0.3198
p_Ad5 = 4.659*10**(-4)
e_Ad5 = 22.87
c_Ad5 = 0.2127
alpha_Ad5 = 5.053*10**(-12)

params_Ad5 = [lam_Ad5, beta_Ad5, k_Ad5, ro_Ad5, p_Ad5, e_Ad5, c_Ad5, alpha_Ad5]
y0_Ad5 = [41, 0, 0, 8.440*10**(-3), 0.0]

# Parameter of Ad6d24.P19

lam_Ad6 = 0.0996
beta_Ad6 = 45.35
k_Ad6 = 0.3189
ro_Ad6 = 0.3162
p_Ad6 = 2.073*10**(-3)
e_Ad6 = 6295
c_Ad6 = 0.04982
alpha_Ad6 = 2.421*10**(-7)

params_Ad6 = [lam_Ad6, beta_Ad6, k_Ad6, ro_Ad6, p_Ad6, e_Ad6, c_Ad6, alpha_Ad6]
y0_Ad6 = [41, 0, 0, 7.078*10**(-4), 0.0]



