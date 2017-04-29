project(Utils)

set(HEADER
	agtkCommandLineArgumentParser.h
	agtkCommandLineArgumentParser.hxx
)

add_custom_target(Utils ALL SOURCES ${HEADER})