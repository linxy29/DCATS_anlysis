## This file is used to check the situation of mu
library(DCATS)

ctrl_corrected = Haber2017$count_ctrl
for (i in seq_len(nrow(Haber2017$count_ctrl))) {
  ctrl_corrected[i, ] <- sum(Haber2017$count_ctrl[i, ]) *
    multinom_EM(Haber2017$count_ctrl[i, ], Haber2017$svm_mat, verbose = FALSE)$mu
}
print(ctrl_corrected)

Hpoly3_corrected = Haber2017$count_Hpoly3
for (i in seq_len(nrow(Haber2017$count_Hpoly3))) {
  Hpoly3_corrected[i, ] <- sum(Haber2017$count_Hpoly3[i, ]) *
    multinom_EM(Haber2017$count_Hpoly3[i, ], Haber2017$svm_mat, verbose = FALSE)$mu
}
print(Hpoly3_corrected)

Hpoly10_corrected = Haber2017$count_Hpoly10
for (i in seq_len(nrow(Haber2017$count_Hpoly10))) {
  Hpoly10_corrected[i, ] <- sum(Haber2017$count_Hpoly10[i, ]) *
    multinom_EM(Haber2017$count_Hpoly10[i, ], Haber2017$svm_mat, verbose = FALSE)$mu
}
print(Hpoly10_corrected)

Salma_corrected = Haber2017$count_Salma
for (i in seq_len(nrow(Haber2017$count_Salma))) {
  Salma_corrected[i, ] <- sum(Haber2017$count_Salma[i, ]) *
    multinom_EM(Haber2017$count_Salma[i, ], Haber2017$svm_mat, verbose = FALSE)$mu
}
print(Salma_corrected)

