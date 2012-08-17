
SET(PERL "${PERL_EXECUTABLE}")

MACRO(ADD_PERL_TEST name command )
 ADD_TEST( ${name}
    ${CMAKE_CURRENT_BINARY_DIR}/../${command} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )
 
 IF(MINC_TEST_ENVIRONMENT)
  set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")
 ENDIF(MINC_TEST_ENVIRONMENT)
ENDMACRO(ADD_PERL_TEST)

MACRO(ADD_CUSTOM_PERL_TEST name command )

 ADD_TEST( ${name}
    ${CMAKE_CURRENT_SOURCE_DIR}/${command} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )
 
 IF(MINC_TEST_ENVIRONMENT)
  set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")
 ENDIF(MINC_TEST_ENVIRONMENT)
 
ENDMACRO(ADD_CUSTOM_PERL_TEST)


ADD_CUSTOM_PERL_TEST(classify_bayes compare_classify_results.pl ${CMAKE_CURRENT_SOURCE_DIR}/bar.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/bar_cls_bayes.mnc.gz --method bayes)
ADD_CUSTOM_PERL_TEST(classify_ann compare_classify_results.pl ${CMAKE_CURRENT_SOURCE_DIR}/bar.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/bar_cls_ann.mnc.gz --method ann)
ADD_CUSTOM_PERL_TEST(classify_bayes compare_classify_results.pl ${CMAKE_CURRENT_SOURCE_DIR}/bar.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/bar_cls_min.mnc.gz --method min)
# knn seem to be broken!
# ADD_CUSTOM_PERL_TEST(classify_knn compare_classify_results.pl ${CMAKE_CURRENT_SOURCE_DIR}/bar.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/bar_cls_knn.mnc.gz --method knn)
