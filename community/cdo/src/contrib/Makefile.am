EXTRA_DIST = makecompl.rb cdoCompletion.bash cdoCompletion.tcsh cdoCompletion.zsh
CLEANFILES = 

completions:
	if (hash ruby &> /dev/null) ; then \
	  ruby -KN $(srcdir)/makecompl.rb -o cdoCompletion; \
	fi

test:
	if (hash ruby &> /dev/null) ; then \
	    (cd ruby;ruby test/test_cdo.rb) \
	fi
	if (hash python &> /dev/null) ; then \
	    (cd python; python test/test_cdo.py) \
	fi

if MAINTAINER_MODE
all-local: completions

CLEANFILES  += `ls cdoCompletion.*`
endif
