SET (src_dir ${CMAKE_CURRENT_SOURCE_DIR})

SET (tpc_bottom_point_extractor_sources
  ${src_dir}/TPCBottomPointExtractorAB.cc
  ${src_dir}/TPCBottomPointExtractorAB.h
  PARENT_SCOPE
)

SET (tpc_bottom_point_extractor_data
  ${src_dir}/TPCBottomPointExtractorAB.xml
  ${src_dir}/TPCBottomPointExtractorAB.xsd
  PARENT_SCOPE
)
