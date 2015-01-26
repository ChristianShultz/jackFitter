
#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
								@CXXFLAGS@  @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@ 

#utils
nobase_include_HEADERS =   \
													 jackFitter/jackFitter.h \
													 jackFitter/ensem_data.h    \
													 jackFitter/ensem_ensem_chisq_fitter.h    \
													 jackFitter/fit_correlators.h    \
													 jackFitter/fit_forms.h    \
													 jackFitter/jackknife_fitter.h    \
													 jackFitter/linear_algebra.h    \
													 jackFitter/multi_princorr_fit.h    \
													 jackFitter/plot.h    \
													 jackFitter/three_point_fit_bias_functions.h    \
													 jackFitter/three_point_fit_compatators.h    \
													 jackFitter/three_point_fit_forms.h    \
													 jackFitter/three_point_fit_functions.h    \
													 jackFitter/three_point_fit.h    \
													 jackFitter/three_point_xml.h   

# the lib
lib_LIBRARIES = libjackFitter.a


libjackFitter_a_SOURCES =	\
													jackFitter/ensem_data.cc    \
													jackFitter/ensem_ensem_chisq_fitter.cc    \
													jackFitter/fit_correlators.cc    \
													jackFitter/fit_forms.cc    \
													jackFitter/jackknife_fitter.cc    \
													jackFitter/linear_algebra.cc    \
													jackFitter/multi_princorr_fit.cc    \
													jackFitter/plot.cc    \
													jackFitter/three_point_fit_bias_functions.cc    \
													jackFitter/three_point_fit.cc    \
													jackFitter/three_point_fit_compatators.cc    \
													jackFitter/three_point_xml.cc    
