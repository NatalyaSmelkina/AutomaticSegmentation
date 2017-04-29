project(Core)

set(HEADER
	HessianAnalysisFilter.h
	HessianAnalysisFilter.hxx
	AirwaysMeasureFilter.h
	AirwaysMeasureFilter.hxx
)

add_custom_target(Core ALL SOURCES ${HEADER})