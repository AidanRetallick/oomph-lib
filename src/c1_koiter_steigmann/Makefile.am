# A few file definitions
#-----------------------

# Define the sources
sourcesC1 =large_displacement_plate_equations.cc large_displacement_plate_curved_bell_elements.cc

# Define the headers
headersC1 = large_displacement_plate_models.h large_displacement_plate_elements.h large_displacement_plate_curved_bell_elements.h \
koiter_steigmann_plate_elements.h foeppl_von_karman_correction_elements.h

# Define name of library
libnameC1 = C1_large_displacement_plate_models

# Combine headers and sources
headersC1_and_sourcesC1 = $(headersC1) $(sourcesC1)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libC1_large_displacement_plate_models.la

# Sources that the library depends on: 
#-------------------------------------
libC1_large_displacement_plate_models_la_SOURCES = $(headersC1_and_sourcesC1)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
C1_large_displacement_plate_modelsinclude_HEADERS = $(libnameC1).h 
C1_large_displacement_plate_modelsincludedir = $(includedir)

#Here's the subdirectory where the actual header files are placed
library_C1_large_displacement_plate_modelsincludedir=$(includedir)/C1_large_displacement_plate_models

#These are the header files that are to be placed in subdirectory
library_C1_large_displacement_plate_modelsinclude_HEADERS=$(headersC1)

# Required libraries -- [assuming that we want to link against stuff in generic
#-------------------     add other oomph-lib libraries if you need them....]
# Include path for library headers -- need to refer to include files
# in their respective source directories as they will not have been 
# installed yet! 
AM_CPPFLAGS += -I$(top_builddir)/src/generic


# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libnameC1).h: dummy_$(libnameC1).h
# $(libnameFvK).h: dummy_$(libnameFvK).h

dummy_$(libnameC1).h:  $(headersC1)
	echo  $(libnameC1) $(headersC1) > all_$(libnameC1).aux
	$(AWK) -f $(top_builddir)/bin/headers.awk < \
	       all_$(libnameC1).aux > $(libnameC1).h
	rm all_$(libnameC1).aux

# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook: 
	(cd $(library_C1_large_displacement_plate_modelsincludedir) && rm -f $(headersC1) )
	(echo "$(headersC1)" > include_files.list )
	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_C1_large_displacement_plate_modelsincludedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


# Tidy up
#--------
clean-local:
	rm -f $(libnameC1).h
