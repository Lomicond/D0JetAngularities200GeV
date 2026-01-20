#!/bin/csh -f

echo "Are you sure that you want to delete all logs, reports and production root files? (y/n): "
set yn = $<
while ("$yn" != "y"&&"$yn" != "n")
echo "Hmm?"
set yn = $<
end

if ( "$yn" == "y" ) then
	rm ./csh/*
	rm ./err/*
	rm ./list/*
	rm ./log/*
	rm ./out/*
	rm ./report/*
	rm ./production/*
	rm ./LocalLibraries.zip	
	rm sched*.dataset
	rm sched*.session.xml
	rm *dataset.tmp
	rm D0JetAnalysis*.dataset
	rm D0JetAnalysis*.session.xml
	rm -r D0JetAnalysis*.package
	rm -r LocalLibraries.package
	echo "It is cleaned up."
else
    echo "Ok, I'll leave it alone" 
endif
endif
