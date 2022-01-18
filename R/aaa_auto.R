# This file has been created automatically by 'rfGenerateConstants'


 ## from  src/AutoRandomFieldsUtils.h
 auto_rfutils_h 	<- as.integer(1)



 MAXUNITS 	<- as.integer(4)
 MAXCHAR 	<- as.integer(18)
 RFOPTIONS 	<- "RFoptions"
 CLASS_TRYERROR 	<- "try-error"

 WARN_UNKNOWN_OPTION_ALL 	<- as.integer(4)
 WARN_UNKNOWN_OPTION_SINGLE 	<- as.integer(3)
 WARN_UNKNOWN_OPTION_CAPITAL 	<- as.integer(2)
 WARN_UNKNOWN_OPTION_NONE1 	<- as.integer(1)
 WARN_UNKNOWN_OPTION_NONE 	<- as.integer(0)

 CONTACT 	<- " Please contact the maintainer martin.schlather@math.uni-mannheim.de.\n"




 ## from  src/AutoRandomFieldsUtilsLocal.h




 LA_INTERN 	<- as.integer(0)
 LA_R 	<- as.integer(1)
 LA_AUTO 	<- as.integer(2)
 LA_GPU 	<- as.integer(3)
 LA_QUERY 	<- as.integer(4)

 LA_LAST 	<- as.integer(LA_QUERY)

 PIVOT_NONE 	<- as.integer(0)
 PIVOT_DO 	<- as.integer(1)
 PIVOT_AUTO 	<- as.integer(2)
 PIVOT_IDX 	<- as.integer(3)
 PIVOT_UNDEFINED 	<- as.integer(4)

 PIVOT_LAST 	<- as.integer(PIVOT_UNDEFINED)

 PIVOTSPARSE_MMD 	<- as.integer(1)
 PIVOTSPARSE_RCM 	<- as.integer(2)

 Inone 	<- as.integer(0)
 Iinstall 	<- as.integer(1)
 Iask 	<- as.integer(2)
 Isse 	<- as.integer(3)
 Isse2 	<- as.integer(4)
 Isse3 	<- as.integer(5)
 Issse3 	<- as.integer(6)
 Iavx 	<- as.integer(7)
 Iavx2 	<- as.integer(8)
 Iavx512f 	<- as.integer(9)
 Igpu 	<- as.integer(10)

 INSTALL_LAST 	<- as.integer(Igpu)


 Cholesky 	<- as.integer(0)
 SVD 	<- as.integer(1)
 Eigen 	<- as.integer(2)
 Sparse 	<- as.integer(3)
 NoInversionMethod 	<- as.integer(4)
 QR 	<- as.integer(5)
 LU 	<- as.integer(6)
 NoFurtherInversionMethod 	<- as.integer(7)
 GPUcholesky 	<- as.integer(8)
 Rcholesky 	<- as.integer(9)
 direct_formula 	<- as.integer(10)
 Diagonal 	<- as.integer(11)



 nr_InversionMethods 	<- as.integer((Diagonal+1))
 nr_user_InversionMethods 	<- as.integer((NoFurtherInversionMethod+1))

 LAST_R_TYPE_NAME 	<- as.integer(32)



