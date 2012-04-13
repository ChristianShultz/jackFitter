#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@


#jackFitter
nobase_include_HEADERS =	jackFitter/ensem_data.h		\
				jackFitter/fit_correlators.h	\
				jackFitter/fit_forms.h		\
				jackFitter/jackknife_fitter.h	\
				jackFitter/linear_algebra.h	\
				jackFitter/multi_princorr_fit.h	\
				jackFitter/plot.h

#
#Library
#

lib_LIBRARIES = libjackFitter.a

#jackFitter 
libjackFitter_a_SOURCES = 	jackFitter/ensem_data.cc		\
				jackFitter/fit_correlators.cc		\
				jackFitter/fit_forms.cc			\
				jackFitter/jackknife_fitter.cc		\
				jackFitter/linear_algebra.cc		\
				jackFitter/multi_princorr_fit.cc	\
				jackFitter/plot.cc
