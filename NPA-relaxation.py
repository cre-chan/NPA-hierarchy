import ncpol2sdpa as ncr
import sympy

# n_vars = 2 
# level = 2
# x = ncr.generate_variables('x',n_vars)

# obj = x[0]*x[1]+x[1]*x[0]
# inequalities = [-x[1]**2 + x[1] + 1/2 >= 0]
# substitutions = {x[0]**2:x[0]}

# sdp = ncr.SdpRelaxation(x)
# sdp.get_relaxation(level, 
#                    objective=obj, 
#                    inequalities=inequalities, 
#                    substitutions=substitutions)

# sdp.write_to_file("example.dat-s")

# try to solve the MSG
level = 1
F = ncr.generate_operators('F',3, hermitian=True)
G = ncr.generate_operators('G',3,hermitian=True)

Poly = 9/2 + 1/2 * F[0] *G[0] + 1/2 * F[1] *G[0] + 1/2 * F[2] * G[0] + 1/2 * F[0] * G[1] + 1/2 * F[1] * G[1] + 1/2 * F[2] * G[1] + 1/2 * F[0] * G[2] + 1/2 * F[1] * G[2] + 1/2 * F[2] * G[2]

Substitutions = {
    op**2: 1 for op in F+G
} | {
   f*g:g*f for f in F for g in G 
}

sdp = ncr.SdpRelaxation(F+G)
sdp.get_relaxation(level, 
                   objective=Poly,
                   substitutions=Substitutions)

sdp.write_to_file("MSG-1.csv")
sdp.write_to_file("MSG-1.dat-s")

# # 変数を定義します。
# F0, F1 = ncr.generate_operators('F', 2,hermitian=True)
# G0, G1 = ncr.generate_operators('G', 2,hermitian=True)
# [I] = ncr.generate_variables('I',1)


# # ゲーム多項式を定義します。
# poly = 0.5 * F0 * G0 + 0.5 * F1 * G0 + 0.5 * F0 * G1 - 0.5 * F1 * G1 + 2*I

# # 制約を定義します。
# constraints = {
#     F0**2 :1,
#     F1**2 : 1,
#     G0**2 : 1,
#     G1**2 : 1,
#     F0 * G0 : G0 * F0,
#     F0 * G1 : G1 * F0,
#     F1 * G0 : G0 * F1,
#     F1 * G1 : G1 * F1
# }

# equalities = [
#     sympy.Eq(I , 1)
# ]

# # SDP問題を生成し、実行します。
# sdp = ncr.SdpRelaxation([F0,F1,G0,G1,I])
# sdp.get_relaxation(1, 
#                    objective=-poly, 
#                    equalities=equalities,
#                    substitutions=constraints)

# sdp.write_to_file("CHSH-1.dat-s")